# vcf_to_mplot
###################################################
pos_sample = ["SRR11561273","SRR11561298","SRR11561324"]
neg_sample = ["SRR11561289","SRR11561314","SRR11561326"]
additional = ["SRR11561288","SRR11561291","SRR11561293"]

rule all:
    input:
         expand('output/vcfdata2/vis_manhattan/positive/{pn_s}_all.png',
            pn_s=pos_sample
            ),
         expand('output/vcfdata2/vis_manhattan/negative/{pn_s}_all.png',
            pn_s=neg_sample
            ),
         expand('output/vcfdata2/vis_manhattan/additional/{pn_s}_all.png',
            pn_s=additional
            ),
         expand('output/vcfdata2/vis_manhattan/positive/{pn_s}_chr_count.csv',
            pn_s=pos_sample
            ),
         expand('output/vcfdata2/vis_manhattan/negative/{pn_s}_chr_count.csv',
            pn_s=neg_sample
            ),
         expand('output/vcfdata2/vis_manhattan/additional/{pn_s}_chr_count.csv',
            pn_s=additional
            )


rule vcf_to_mplot_pos:
    input:
        'data/vcfdata2/positive/{pn_s}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_0.8vaf.vcf',
        'data/vcfdata2/positive/{pn_s}.hg38.identified.snp.fltr.vaf.headerfixed.0.8_1.0vaf.vcf'
    output:
        'output/vcfdata2/vis_manhattan/positive/{pn_s}_all.png',
        'output/vcfdata2/vis_manhattan/positive/{pn_s}_chr_count.csv'
    params:
        args_param_REF = 'C',
        args_param_ALT = 'T',
        args_param_col = 'orange2'
    benchmark:
        'benchmarks/output/vcfdata2/vis_manhattan/positive/{pn_s}_all.txt'
    container:
        'docker://yamaken37/vcf_to_mplot:20230627'
    resources:
        mem_gb=200
    log:
        'logs/output/vcfdata2/vis_manhattan/positive/{pn_s}_all.log'
    shell:
        'src/vcf_to_mplot.sh {input} {output} {params} >& {log}'

rule vcf_to_mplot_neg:
    input:
        'data/vcfdata2/negative/{pn_s}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_0.8vaf.vcf',
        'data/vcfdata2/negative/{pn_s}.hg38.identified.snp.fltr.vaf.headerfixed.0.8_1.0vaf.vcf'
    output:
        'output/vcfdata2/vis_manhattan/negative/{pn_s}_all.png',
        'output/vcfdata2/vis_manhattan/negative/{pn_s}_chr_count.csv'
    params:
        args_param_REF = 'C',
        args_param_ALT = 'T',
        args_param_col = 'orange2'
    benchmark:
        'benchmarks/output/vcfdata2/vis_manhattan/negative/{pn_s}_all.txt'
    container:
        'docker://yamaken37/vcf_to_mplot:20230627'
    resources:
        mem_gb=200
    log:
        'logs/output/vcfdata2/vis_manhattan/negative/{pn_s}_all.log'
    shell:
        'src/vcf_to_mplot.sh {input} {output} {params} >& {log}'

rule vcf_to_mplot_add:
    input:
        'data/vcfdata2/additional/{pn_s}.hg38.identified.snp.fltr.vaf.headerfixed.0.0_0.8vaf.vcf',
        'data/vcfdata2/additional/{pn_s}.hg38.identified.snp.fltr.vaf.headerfixed.0.8_1.0vaf.vcf'
    output:
        'output/vcfdata2/vis_manhattan/additional/{pn_s}_all.png',
        'output/vcfdata2/vis_manhattan/additional/{pn_s}_chr_count.csv'
    params:
        args_param_REF = 'C',
        args_param_ALT = 'T',
        args_param_col = 'orange2'
    benchmark:
        'benchmarks/output/vcfdata2/vis_manhattan/additional/{pn_s}_all.txt'
    container:
        'docker://yamaken37/vcf_to_mplot:20230627'
    resources:
        mem_gb=200
    log:
        'logs/output/vcfdata2/vis_manhattan/additional/{pn_s}_all.log'
    shell:
        'src/vcf_to_mplot.sh {input} {output} {params} >& {log}'