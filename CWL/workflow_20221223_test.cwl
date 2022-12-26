class: Workflow
cwlVersion: v1.0
id: workflow_20221223
label: workflow_20221223
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: fastq2
    type: File?
    'sbg:x': -596.597900390625
    'sbg:y': -127
  - id: fastq1
    type: File
    'sbg:x': -568.597900390625
    'sbg:y': 11
  - id: reference
    type: File
    'sbg:x': 88.84902954101562
    'sbg:y': 279.6493835449219
  - id: reference_1
    type: File
    'sbg:x': 156.00564575195312
    'sbg:y': -235.32728576660156
outputs: []
steps:
  - id: trim_galore
    in:
      - id: fastq1
        source: fastq1
      - id: fastq2
        source: fastq2
    out:
      - id: fastq1_trimmed
      - id: fastq1_trimmed_unpaired
      - id: fastq2_trimmed
      - id: fastq2_trimmed_unpaired
      - id: trim_galore_log
      - id: trimmed_fastqc_html
      - id: trimmed_fastqc_zip
    run: ./trim_galore.cwl
    'sbg:x': -441.3398742675781
    'sbg:y': -64.79248809814453
  - id: star_mapping
    in:
      - id: file1
        source: trim_galore/fastq1_trimmed
      - id: file2
        source: trim_galore/fastq2_trimmed
    out:
      - id: output
      - id: unmapped_reads
    run: ./star_mapping.cwl
    'sbg:x': -477.4384765625
    'sbg:y': 145.74627685546875
  - id: picard__add_or_replace_read_groups
    in:
      - id: INPUT
        source: star_mapping/output
    out:
      - id: sequences_with_new_read_group
    run: ./picard_AddOrReplaceReadGroups.cwl
    'sbg:x': -297.0946350097656
    'sbg:y': 74.38215637207031
  - id: _g_a_t_k__mark_duplicates
    in:
      - id: InputFile
        source: picard__add_or_replace_read_groups/sequences_with_new_read_group
    out:
      - id: alignment
      - id: index
      - id: metrics
      - id: vcf
    run: ./GATK-MarkDuplicates.cwl
    'sbg:x': -206.45425415039062
    'sbg:y': -92.01351928710938
  - id: gatk_split
    in:
      - id: input
        source: _g_a_t_k__mark_duplicates/alignment
    out:
      - id: output
    run: ./gatk_split.cwl
    'sbg:x': -15.105905532836914
    'sbg:y': -101.7102279663086
  - id: gatk_variant
    in:
      - id: input
        source: gatk_split/output
    out:
      - id: assembly-region-out
      - id: bam-output
      - id: graph-output
      - id: output
    run: ./gatk_variant.cwl
    'sbg:x': -55
    'sbg:y': 129.9819793701172
  - id: gatk_dbimport
    in:
      - id: vcffile
        source: gatk_variant/output
    out:
      - id: output
    run: ./gatk_dbimport.cwl
    'sbg:x': 163.4736328125
    'sbg:y': -92.41595458984375
  - id: gatk__genotype_g_v_c_fs_v4_0_11_0
    in:
      - id: reference
        source: reference
      - id: variant
        source: gatk_dbimport/output
    out:
      - id: vcf
    run: ./gatk_genotypegvcfs.cwl
    label: gatk-GenotypeGVCFs-v4.0.11.0
    'sbg:x': 192
    'sbg:y': 95.31230163574219
  - id: gatk__variant_filtration_v4_0_11_0
    in:
      - id: reference
        source: reference_1
      - id: variant
        source: gatk__genotype_g_v_c_fs_v4_0_11_0/vcf
    out:
      - id: vcf
    run: ./gatk_variantfilter.cwl
    label: gatk-VariantFiltration-v4.0.11.0
    'sbg:x': 339.3415832519531
    'sbg:y': -76.66967010498047
  - id: bcftools_vaf
    in:
      - id: variantfile
        source: gatk__variant_filtration_v4_0_11_0/vcf
    out:
      - id: output
    run: ./bcftools_vaf.cwl
    'sbg:x': 373.9954833984375
    'sbg:y': 121.9842300415039
  - id: picard_fix_vcf_header
    in:
      - id: variantfile
        source: bcftools_vaf/output
    out:
      - id: output
    run: ./picard_fix_vcf_header.cwl
    'sbg:x': 490.3258361816406
    'sbg:y': 109.99324035644531
  - id: bcftools_view
    in:
      - id: variantfile
        source: picard_fix_vcf_header/output
    out:
      - id: output
    run: ./bcftools_view.cwl
    'sbg:x': 609.0045166015625
    'sbg:y': 9.013520240783691
  - id: bcftools_view_1
    in:
      - id: variantfile
        source: picard_fix_vcf_header/output
    out:
      - id: output
    run: ./bcftools_view.cwl
    'sbg:x': 600.3077392578125
    'sbg:y': 181.31231689453125
requirements: []
