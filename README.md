# gridss2_somatic_sv_call
* This is a snakemake version of gridss2 (https://github.com/PapenfussLab/gridss).

* GRIDSS execution procedures are well described at 
https://github.com/PapenfussLab/gridss/blob/master/QuickStart.md#optimise-gridss-execution

* After assembly, it performs Normal/Tumor somatic variant calling using `gripss`.
https://github.com/hartwigmedical/hmftools/tree/master/gripss

* GRIDSS being an SV breakpoint caller initially reports all variants as BND but can be post-annotated as Trans, Inv, Ins, Del and Dup with a R script (simple-event-annotation.R) which is used in the last rule `rule create_bedpe_with_simpleSVannotation` in snakemake.

* Why all calls BND?
 https://github.com/PapenfussLab/gridss/blob/master/Readme.md#why-are-all-calls-bnd
  
* Please specify the samples to be analyzed in config/samples.yaml, following the format explained therein.

