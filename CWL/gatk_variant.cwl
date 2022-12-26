cwlVersion: v1.0
class: CommandLineTool
baseCommand: [gatk]
arguments: [--java-options, "-Djava.io.tmpdir=/lscratch/$SLURM_JOBID -Xms200G -Xmx200G -XX:ParallelGCThreads=14", HaplotypeCaller, -R, /home/bonohulab2/Nakamae/RNA_offtarget_analysis/Homo_sapiens.GRCh38.dna.primary_assembly.fa, -I, $(inputs.input), -O, $(inputs.output), -ERC, GVCF, --tmp-dir, ./tmp, --sample-name, $(inputs.samplename)]
inputs:
  - id: input
    type: File
  - id: output
    type: string
  - id: samplename
    type: string

outputs:
  vcf:
    type: File
    outputBinding:
      glob: $(inputs.output)
    secondaryFiles:
      - .tbi  


