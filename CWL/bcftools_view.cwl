cwlVersion: v1.0
class: CommandLineTool
baseCommand: [bcftools]
arguments: [view, -O, z, -o, $(inputs.output), -e, $(inputs.vafsetting), $(inputs.variantfile)]
inputs:
  - id: variantfile
    type: File
  - id: output
    type: string
  - id: vafsetting
    type: string


outputs:
  output:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*"