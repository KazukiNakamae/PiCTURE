cwlVersion: v1.0
class: CommandLineTool
baseCommand: [gatk]
arguments: [GenomicsDBImport, -R, ../../Nakamae/RNA_offtarget_analysis/Homo_sapiens.GRCh38.dna.primary_assembly.fa, -V, $(inputs.vcffile), -L, intervals.list, --genomicsdb-workspace-path, $(runtime.outdir)]
inputs:
  - id: vcffile
    type: File
  - id: path
    type: string

outputs:
  output:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*"