cwlVersion: v1.0
class: CommandLineTool
baseCommand: [gatk]
arguments: [SplitNCigarReads, -R, /home/bonohulab2/Nakamae/RNA_offtarget_analysis/Homo_sapiens.GRCh38.dna.primary_assembly.fa, -I, $(inputs.input), -O, $(inputs.output)]
inputs:
  - id: input
    type: File
  - id: output
    type: string


outputs:
  # output:
  #   type:
  #     type: array
  #     items: File
  #   outputBinding:
  #     glob: "*"
  - id: output
    doc: Output file from corresponding to the input argument output-filename
    type: File
    outputBinding:
      glob: $(inputs.output)