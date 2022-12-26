cwlVersion: v1.0
class: CommandLineTool
baseCommand: [STAR]
arguments: [--genomeDir, /home/bonohulab2/Nakamae/RNA_offtarget_analysis/hg38_index/, --readFilesIn, $(inputs.file1), $(inputs.file2), --runThreadN $(inputs.runThreadN), --outSAMtype, BAM, SortedByCoordinate, --outFileNamePrefix, $(inputs.outFileNamePrefix)]
inputs:
  - id: file1
    type: File
  - id: file2
    type: File
  - id: runThreadN
    type: int
  - id: outFileNamePrefix
    type: string

outputs:
  output:
    type:
     - File
     - File[]
    outputBinding:
      glob: "*.bam"

  unmapped_reads:
    type: ["null", File]
    outputBinding:
      glob: "Unmapped.out*"