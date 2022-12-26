#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

id: gatk-VariantFiltration-v4.0.11.0
label: gatk-VariantFiltration-v4.0.11.0

requirements:
    InlineJavascriptRequirement: {}
    DockerRequirement:
        # dockerPull: broadinstitute/gatk:4.0.11.0
        dockerPull: quay.io/biocontainers/gatk4:4.1.6.0--py38_0


# baseCommand: [ gatk, --java-options, -Xmx4G, VariantFiltration ]
baseCommand: [ gatk, VariantFiltration ]

inputs:
  variant:
    type: File
    doc: Input VCF file
    inputBinding:
      prefix: --variant
    # secondaryFiles:
    #   - .tbi
  reference:
    type: File
    doc: Reference FASTA file
    inputBinding:
      prefix: --reference
    # secondaryFiles:
    #   - .fai
    #   - ^.dict
  output:
    type: string
    doc: Output VCF file name
    # default: variants.filter.genotype.vcf.gz
    default: variants.filter.genotype.vcf
    inputBinding:
      prefix: --output
  
  filter-expression:
    type: string
    doc: Filter condition
    default: "QD < 2.0 || QUAL < 30.0 || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
    inputBinding:
      prefix: --filter-expression


arguments:
  - id: filter-name
    prefix: --filter-name
    valueFrom: FILTER

outputs:
  vcf:
    type: File
    outputBinding:
      glob: $(inputs.output)
    # secondaryFiles:
    #   - .tbi