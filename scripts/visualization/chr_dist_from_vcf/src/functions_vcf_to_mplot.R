#library
##################################################
library(qqman)
library(vcfR)
library(tidyverse)
##################################################
.str_vaf = function(y) {
    vcf_gt_TC[y,2] %>%
        str_split(pattern = ":")  %>%
            .[[1]]  %>%
                tail(.,1) -> return_object
    return(return_object)
}
.prepare_manh = function(x) {
##### load filter#####
# load
read_vcf <- read.vcfR(x)
# add ID
read_vcf@fix |>
    as.data.frame() |>
    mutate(ID = row_number()) -> vcf_fix
# get T-C ID
vcf_fix %>%
    dplyr::filter(REF==args_param_REF) %>%
    dplyr::filter(ALT==args_param_ALT) %>%
    .$ID -> id_T_C

#### Filter REF T ALT C####
vcf_fix[id_T_C,] -> vcf_fix_TC
read_vcf@gt %>%
    as.data.frame() %>%
    .[id_T_C,]-> vcf_gt_TC

#### transform df####
seq(1:length(vcf_gt_TC[,2])) %>%
    purrr::map(., .str_vaf) %>%
    unlist() -> vaf_value

df_trim_TC <- data.frame(CHR =vcf_fix_TC$CHROM,
                         POS = vcf_fix_TC$POS,
                         VAF_value = vaf_value,
                         stringsAsFactors = FALSE,
                         row.names = NULL
)
# remove chrY row
df_trim_TC[!(df_trim_TC$CHR=="chrY"), ] -> df_trim_Y

# chrX -> chr23
str_replace_all(df_trim_Y$CHR,
                pattern = c("chrX" = "chr23")
) -> df_trim_Y$CHR
# tirm chr
str_replace_all(df_trim_Y$CHR,
                pattern = c("chr" = "")
) -> df_trim_Y$CHR
df_trim_Y$CHR <- as.integer(df_trim_Y$CHR)

df <- data.frame(CHR =df_trim_Y$CHR,
                 BP = as.numeric(df_trim_Y$POS),
                 P = as.numeric(df_trim_Y$VAF_value)
)
return(df)
}