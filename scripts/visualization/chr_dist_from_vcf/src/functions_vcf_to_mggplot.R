#library
##################################################
library(qqman)
library(vcfR)
library(tidyverse)
library(ggplot2)
library(patchwork)
##################################################
.str_vaf = function(x) {
    vcf_gt_TC[x,2] %>%
        str_split(pattern = ":")  %>%
            .[[1]]  %>%
                tail(.,1) -> return_object
    return(return_object)
}
.ggmanhatten = function(x) {
    ggplot(x, 
           aes(x = chromosome, 
               y = pvalue, 
               color = chromosome)
    ) +
        geom_jitter(width = 0.3, 
                    size = 1, 
                    alpha = 0.6) +
        scale_color_discrete(name  = "Chromosome") +
        labs(x = "Chromosome", y = "VAF") +
        theme_minimal() +
        theme(legend.position = "none") +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) -> return_object
    return(return_object)
}