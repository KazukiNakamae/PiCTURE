class: CommandLineTool
cwlVersion: v1.0

baseCommand:
  - STAR
inputs:
  - default: zcat
    id: readFilesCommand
    type: string
  - default:
      - BAM
      - SortedByCoordinate
    id: outSAMtype
    type: "string[]"
    inputBinding:
      position: 1
      prefix: "--outSAMtype"
      shellQuote: false
    label: type of SAM/BAM output
    doc: >-
      1st word: BAM: output BAM without sorting, SAM: output SAM without
      sorting, None: no SAM/BAM output, 2nd, 3rd: Unsorted: standard unsorted,
      SortedByCoordinate: sorted by coordinate. This option will allocate extra
      memory for sorting which can be specified by -limitBAMsortRAM

  # - id: readFilesIn
  #   type:
  #    - File
  #    - File[]
  - id: fq1
    type: File
    inputBinding:
      position: 50
      prefix: "--readFilesIn"
    label: path to file that contain input read1
    doc: path to file that contain input read1
  - id: fq2
    type: File
    inputBinding:
      position: 51
    label: path to file that contain input read2
    doc: path to file that contain input read2
  - id: outFileNamePrefix
    type: string?  
  - id: runThreadN
    type: int
  - id: genomeDir
    type: Directory

outputs:
  - id: aligned
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.outSAMtype == "BAM Unsorted")
            return p+"Aligned.out.bam";
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Aligned.sortedByCoord.out.bam";
        }
    # format: edam:format_1930  # FASTQ
  - id: bamRemDups
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.bamRemoveDuplicatesType != "UniqueIdentical")
            return null;
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Processed.out.bam";
        }
  - id: mappingstats
    type: File?
    outputBinding:
      loadContents: true
      glob: |
        ${
          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Log.final.out";
        }
    # format: edam:format_1930  # FASTQ
  - id: readspergene
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"ReadsPerGene.out.tab";
        }
  - id: transcriptomesam
    type: 'File[]'
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Aligned.toTranscriptome.out.bam";
        }
  - id: Log.out
    type: File?
    outputBinding:
      glob: >-
        ${   var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return p+"Log.out"; }
  - id: Log.progress.out
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Log.progress.out";
        }
  - id: SJ.out.tab
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"SJ.out.tab";
        }
doc: >-
  STAR: Spliced Transcripts Alignment to a Reference.
  https://github.com/alexdobin/STAR/blob/main/doc/STARmanual.pdf
label: "STAR mapping: running mapping jobs."
requirements:
  - class: ShellCommandRequirement
  # - class: InitialWorkDirRequirement
  #   listing:
  #     - entry: $(inputs.genomeDir)
  #       writable: true
  #     - entry: $(inputs.readFilesIn)
  #       writable: true
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: "kazukinakamae/star_for_human_gatk:latest"