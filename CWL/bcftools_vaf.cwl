cwlVersion: v1.0
class: CommandLineTool
baseCommand: [bcftools]
arguments: [+fill-tags, $(inputs.variantfile), -Ov, -o, $(inputs.output), --, -t, VAF]
inputs:
  - id: variantfile
    type: File
  - id: output
    type: string


outputs:
  output:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*"