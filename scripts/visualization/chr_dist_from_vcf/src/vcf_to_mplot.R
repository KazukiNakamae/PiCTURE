source("src/functions_vcf_to_mplot.R")
#### args####
args <- commandArgs(trailingOnly = TRUE)
args_input <- args[1]
args_input_2 <- args[2]
args_output <- args[3]
args_output_2 <- args[4]
args_param_REF <- args[5]
args_param_ALT <- args[6]
args_param_col <- args[7]
#### test args####
# args_input <-c("data/vcfdata2/positive/SRR11561273.hg38.identified.snp.fltr.vaf.headerfixed.0.0_0.8vaf.vcf")
# args_input_2 <-c("data/vcfdata2/positive/SRR11561273.hg38.identified.snp.fltr.vaf.headerfixed.0.8_1.0vaf.vcf")
# args_output <-c("output/vcfdata2/vis_manhattan/positive/SRR11561273_all.png")
# args_output_2 <-c("output/vcfdata2/vis_manhattan/positive/SRR11561273_all.csv")
# args_param_REF <- c("T")
# args_param_ALT <- c("C")
# args_param_col <- c("orange2")
#### title####
args_output %>%
    str_replace(., "output/vcfdata2/vis_manhattan/", "") %>%
    str_replace(., "/", "_") -> args_title
#### df_1#####
# load
read_vcf <- read.vcfR(args_input)
# add ID
read_vcf@fix |>
    as.data.frame() |>
        mutate(ID = row_number()) -> vcf_fix
# get T-C ID
vcf_fix %>%
    dplyr::filter(REF==args_param_REF) %>%
        dplyr::filter(ALT==args_param_ALT) %>%
            .$ID -> id_T_C
# Filter REF T ALT C
vcf_fix[id_T_C,] -> vcf_fix_TC
read_vcf@gt %>%
    as.data.frame() %>%
        .[id_T_C,]-> vcf_gt_TC

# transform df
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

df_1 <- data.frame(CHR =df_trim_Y$CHR,
                 BP = as.numeric(df_trim_Y$POS),
                 P = as.numeric(df_trim_Y$VAF_value)
                 )
#### df_2#####
# load
read_vcf <- read.vcfR(args_input_2)
# add ID
read_vcf@fix |>
    as.data.frame() |>
    mutate(ID = row_number()) -> vcf_fix
# get T-C ID
vcf_fix %>%
    dplyr::filter(REF==args_param_REF) %>%
    dplyr::filter(ALT==args_param_ALT) %>%
    .$ID -> id_T_C
# Filter REF T ALT C
vcf_fix[id_T_C,] -> vcf_fix_TC
read_vcf@gt %>%
    as.data.frame() %>%
    .[id_T_C,]-> vcf_gt_TC

# transform df
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

df_2 <- data.frame(CHR =df_trim_Y$CHR,
                   BP = as.numeric(df_trim_Y$POS),
                   P = as.numeric(df_trim_Y$VAF_value)
                   )

#### df####
bind_rows(df_1, df_2)->df
#### manhattan plot####
mutate(df,
       SNP = paste(CHR, BP)
) ->df_SNP
str_replace_all(df_SNP$SNP,
                pattern = c(" " = "_")
) ->df_SNP$SNP

png(filename=args_output,
    width=800,
    height=600
)

manhattan(df_SNP,
          main =args_title,
          logp = FALSE,
          col = args_param_col,
          ylab =c("VAF")
)

dev.off()

df %>% group_by(CHR) %>% count() -> count_chr
# CSVファイルとして保存
write_csv(count_chr, args_output_2)
