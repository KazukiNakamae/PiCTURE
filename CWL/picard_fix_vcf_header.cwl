cwlVersion: v1.0
class: CommandLineTool
baseCommand: [picard]
arguments: [FixVcfHeader, -INPUT, $(inputs.variantfile), -OUTPUT, $(inputs.output)]
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