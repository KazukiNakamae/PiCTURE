source("src/functions_vcf_to_mggplot.R")
#### args####
args <- commandArgs(trailingOnly = TRUE)
args_input <- args[1]
args_output <- args[2]
args_param_REF <- args[3]
args_param_ALT <- args[4]

#### test args####
# args_input <-c("data/vcfdata2/positive/SRR11561273.hg38.identified.snp.fltr.vaf.headerfixed.0.0_0.8vaf.vcf")
# args_input_2 <-c("data/vcfdata2/positive/SRR11561273.hg38.identified.snp.fltr.vaf.headerfixed.0.8_1.0vaf.vcf")
# args_output <-c("output/vcfdata2/vis_manhattan/positive/SRR11561273.png")
# args_param_REF <- c("T")
# args_param_ALT <- c("C")

#### title####
args_output %>%
    str_replace(., "output/vcfdata2/vis_manhattan/", "") %>%
    str_replace(., ".hg38.identified.snp.fltr.vaf.headerfixed.", "") %>%
    str_replace(., "vaf.vcf", "") %>%
    str_replace(., "/", "_") -> plot_title

#### load filter#####
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

df <- data.frame(chromosome =df_trim_Y$CHR,
                 position = as.numeric(df_trim_Y$POS),
                 pvalue = as.numeric(df_trim_Y$VAF_value)
                 )
#### manhattan ggplot####
df <- df %>% 
    arrange(chromosome, position) %>%
    mutate(chromosome = factor(chromosome))
.ggmanhatten(df) -> ggm_1

#### load filter#####
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

df_2 <- data.frame(chromosome =df_trim_Y$CHR,
                 position = as.numeric(df_trim_Y$POS),
                 pvalue = as.numeric(df_trim_Y$VAF_value)
)
#### manhattan ggplot####
df_2 <- df_2 %>% 
    arrange(chromosome, position) %>%
    mutate(chromosome = factor(chromosome))
.ggmanhatten(df_2) -> ggm_2

#### patchwork####
gg <- ggm_1 +
    ggm_2 +
    plot_layout(ncol = 1) +
    plot_annotation(title = plot_title,
                    caption = 'made with patchwork',
                    theme = theme(plot.title = element_text(size = 50, hjust = 0.5))
    )
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 10.0, 
       height = 10.0,
       limitsize = FALSE
)