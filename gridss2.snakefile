'''
ml SAMtools
ml BWA
ml GATK/4.1.4.1-GCCcore-8.3.0-Java-11
ml Python/3.7.4-foss-2019b-fh1
ml Pysam/0.15.4-GCC-8.3.0-Python-3.7.4
ml PyYAML/5.1.2-GCCcore-8.3.0-Python-3.7.4
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml R/3.6.2-foss-2019b-fh1
ml Java/11.0.2

#To run this snakemake on the cluster (remove -np): 
snakemake -s gridss2.snakefile --latency-wait 60 --restart-times 3 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} --constraint=gizmok" -j 20 -np
'''

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

preprocess_endings = (
    ".cigar_metrics",
    ".coverage.blacklist.bed",
    ".idsv_metrics",
    ".insert_size_histogram.pdf",
    ".insert_size_metrics",
    ".mapq_metrics",
    ".sv.bam",
    ".sv.bam.bai",
    ".tag_metrics",
    )

assembly_endings = (
    ".cigar_metrics",
    ".coverage.blacklist.bed",
    ".downsampled_0.bed",
    ".excluded_0.bed",
    ".idsv_metrics",
    ".mapq_metrics",
    ".quality_distribution.pdf",
    ".quality_distribution_metrics",
    ".subsetCalled_0.bed",
    ".sv.bam",
    ".sv.bam.bai",
    ".tag_metrics",
    )


reference_index_endings = (".amb",".ann", ".bwt", ".pac", ".sa", ".gridsscache", ".img")

groups = list(config["groups"].keys())
#tumors = [val.split(',')[1] for val in config["groups"].values()]
#normals = [val.split(',')[0] for val in config["groups"].values()]

def get_group_samples_bams(wildcards):
    sample_list = config["groups"][wildcards.group].split(',')
    return [config["samples"][sample] for sample in sample_list]

def get_group_samples(wildcards):
    return config["groups"][wildcards.group].split(',')


print(f"unique groups: {groups}")

rule all:
    input: 
        expand("results/{sample}.bam.gridss.working/{sample}.bam{ending}", sample=config["samples"], ending=preprocess_endings),
        expand("results/assembly/group.{group}.bam", group=groups),
        expand("results/group.{group}.bam.gridss.working/group.{group}.bam{ending}", group=groups, ending=assembly_endings),
        expand("results/vcf/group.{group}.vcf", group=groups),
        expand("results/vcf/{group}.gripss.filtered.vcf.gz", group=groups),
        expand("results/vcf/{group}.gripss.filtered.vcf", group=groups),
        expand("results/bedpe_with_simple_annotation/{group}.gripss.filtered.bedpe", group=groups),


#Perfrom preprocess, all samples can be run in parallel
rule gridss_preprocess:
    input:
        bam=lambda wildcards: config["samples"][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.sample]+".bai",
        reference=config["reference_fasta"],
        dictionary=config["reference_dict"],
        refindex=multiext(config["reference_fasta"], *reference_index_endings)
    output: #per sample output
        multiext("results/{sample}.bam.gridss.working/{sample}.bam", *preprocess_endings)
    params:
        gridss_jar=config["gridss_jar"],
        path_to_gridss=config["path_to_gridss"]
    resources:
        mem=8
    log:
        "logs/preprocess/{sample}.preprocess.log"
    shell:
        "{params.path_to_gridss} -s preprocess {input.bam} -w results -j {params.gridss_jar} -r {input.reference} > {log} 2> {log}"

#Perform GRIDSS breakend assembly- pefrom this for each group (N-T paired)
rule gridss_assemble:
    input:
        bams=lambda wildcards: get_group_samples_bams(wildcards),
        reference=config["reference_fasta"],
        #refindex=multiext(config["reference_fasta"], *reference_index_endings),
        preprocess=lambda wc: expand("results/{sample}.bam.gridss.working/{sample}.bam{ending}", sample=get_group_samples(wc), ending=preprocess_endings),
    output:
        assembly="results/assembly/group.{group}.bam",
        assembly_others=multiext("results/group.{group}.bam.gridss.working/group.{group}.bam", *assembly_endings)
    params:
        gridss_jar=config["gridss_jar"],
        path_to_gridss=config["path_to_gridss"],
        input_args=lambda wildcards: " ".join(get_group_samples_bams(wildcards)),
    threads:
        10
    log:
        "logs/assemble/group.{group}.log"
    shell: 
        "{params.path_to_gridss} -s assemble -a {output.assembly} -w results --threads {threads} {params.input_args} -j {params.gridss_jar} -r {input.reference} > {log} 2> {log}"

rule gridss_call:
    input:
        bams=lambda wildcards: get_group_samples_bams(wildcards),
        reference=config["reference_fasta"],
        assembly="results/assembly/group.{group}.bam",
    output:
        vcf="results/vcf/group.{group}.vcf",
        idx="results/vcf/group.{group}.vcf.idx",
        #tmpidx=temp("results/group.{group}.vcf.gridss.working/group.{group}.vcf.allocated.vcf.idx") # be aware the group occurs multiple times here
    params:
        gridss_jar=config["gridss_jar"],
        path_to_gridss=config["path_to_gridss"],
        input_args=lambda wildcards: " ".join(get_group_samples_bams(wildcards)),
    log:
        "logs/call/group.{group}.log"
    threads:
        10
    shell:
        "{params.path_to_gridss} -s call -a {input.assembly} -w results --threads {threads} {params.input_args} -j {params.gridss_jar} -r {input.reference}  -o {output.vcf} > {log} 2> {log}"

#Perform normal/tumor somatic variant calling
#https://github.com/hartwigmedical/hmftools/tree/master/gripss
rule gridss_somatic_filter:
    input:
        vcf="results/vcf/group.{group}.vcf",
        idx="results/vcf/group.{group}.vcf.idx",
        pon_sgl_file=config["pon_dir"]+"gridss_pon_single_breakend.bed",
        pon_sv_file=config["pon_dir"]+"gridss_pon_breakpoint.bedpe",
        reference=config["reference_fasta"],
    params:
        output_dir="results/vcf",
        path_to_gridss_somatic_filter=config["path_to_gripss"],
        java=config["java"],   
        normal_sample=lambda wildcards: get_group_samples(wildcards)[0],
        tumor_sample=lambda wildcards: get_group_samples(wildcards)[1],
    output:
        output="results/vcf/{group}.gripss.filtered.vcf.gz", #this doesn't work when I change group to tumor
        output_all="results/vcf/{group}.gripss.vcf.gz" #this doesn't work when I change group to tumor
    log:
        "logs/call/somatic_filter.{group}.log"
    threads:
        10
    shell:     
        """
        {params.java} -jar {params.path_to_gridss_somatic_filter} -vcf {input.vcf} -sample {params.tumor_sample} -reference {params.normal_sample} -ref_genome {input.reference} -pon_sgl_file {input.pon_sgl_file} -pon_sv_file {input.pon_sv_file} -output_dir {params.output_dir} > {log} 2> {log}
        
        # Use this only when group name is different from tumor sample name (uncomment this line in that case)
        # mv results/vcf/{params.tumor_sample}.gripss.filtered.vcf.gz {output.output}
        # mv results/vcf/{params.tumor_sample}.gripss.filtered.vcf.gz.tbi {output.output}.tbi
        # mv results/vcf/{params.tumor_sample}.gripss.vcf.gz {output.output_all}
        # mv results/vcf/{params.tumor_sample}.gripss.vcf.gz.tbi {output.output_all}.tbi

        """      

rule unzip_somatic_vcf:
    input:
        "results/vcf/{group}.gripss.filtered.vcf.gz"
    output:
        "results/vcf/{group}.gripss.filtered.vcf"
    shell:
        """
        gunzip -c {input} > {output}
        """
        
rule create_bedpe_with_simpleSVannotation:
    input:
        "results/vcf/{group}.gripss.filtered.vcf"
    output:
        "results/bedpe_with_simple_annotation/{group}.gripss.filtered.bedpe"
    params:
        path_to_bedpe_src=config["path_to_bedpe_src"],
        bs_genome=config["bs_genome"],
        output_dir="results/bedpe_with_simple_annotation/"   
    shell:
        "Rscript {params.path_to_bedpe_src} {input} {params.output_dir} {params.bs_genome}"

