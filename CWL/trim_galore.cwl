class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
  edam: http://edamontology.org/
$schemas:
  - https://edamontology.org/EDAM_1.18.owl
  
baseCommand:
  - trim_galore
inputs:
  - id: adapter1
    type: string?
    doc: |
      Adapter sequence for first reads.
      if not specified, trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3' Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the string to "illumina", "nextera", or "small_rna".
  - id: adapter2
    type: string?
    doc: |
      Adapter sequence for second reads - only for paired end data.
      if not specified, trim_galore will try to autodetect whether ...
      - Illumina universal adapter (AGATCGGAAGAGC)
      - Nextera adapter (CTGTCTCTTATA)
      - Illumina Small RNA 3' Adapter (TGGAATTCTCGG)
      ... was used.
      You can directly choose one of the above configurations
      by setting the adapter1 string to "illumina", "nextera", or "small_rna".
  - id: fastq1
    type: File
    inputBinding:
      position: 10
    doc: |
      raw reads in fastq format; can be gzipped;
      if paired end, the file contains the first reads;
      if single end, the file contains all reads
  - id: fastq2
    type: File?
    inputBinding:
      position: 11
    doc: |
      (optional) raw reads in fastq format; can be gzipped;
      if paired end, the file contains the second reads;
      if single end, the file does not exist
  - default: 1
    id: min_adapter_overlap
    type: int
    inputBinding:
      position: 1
      prefix: '--stringency'
    doc: minimum overlap with adapter seq in bp needed to trim
  - default: 20
    id: min_read_length
    type: int
    inputBinding:
      position: 1
      prefix: '--length'
    doc: discard reads that get shorter than this value
  - default: 35
    id: min_unpaired_read_rescue_length
    type: int
    doc: |
      if only one read of a pair passes the qc and adapter trimming,
      it needs at least this length to be rescued
  - default: 20
    id: qual_trim_cutoff
    type: int
    inputBinding:
      position: 1
      prefix: '--quality'
    doc: trim all base with a phred score lower than this valueFrom
outputs:
  - id: fastq1_trimmed
    type: File
    outputBinding:
      glob: |
        ${
            if ( inputs.fastq2 == null  ){ return "*trimmed.fq*" }
            else { return "*val_1.fq*" }
        }
    format: edam:format_1930  # FASTQ
  - id: fastq1_trimmed_unpaired
    type: File?
    outputBinding:
      glob: '*unpaired_1.fq*'
  - id: fastq2_trimmed
    type: File?
    outputBinding:
      glob: '*val_2.fq*'
    format: edam:format_1930  # FASTQ
  - id: fastq2_trimmed_unpaired
    type: File?
    outputBinding:
      glob: '*unpaired_2.fq*'
  - id: trim_galore_log
    type: 'File[]'
    outputBinding:
      glob: '*trimming_report.txt'
  - id: trimmed_fastqc_html
    doc: html report of post-trimming fastqc
    type: 'File[]'
    outputBinding:
      glob: '*fastqc.html'
  - id: trimmed_fastqc_zip
    doc: all data of post-trimming fastqc e.g. figures
    type: 'File[]'
    outputBinding:
      glob: '*fastqc.zip'
doc: |
  Adaptor trimming of reads (single or paired end) in fastq format.
arguments:
  - position: 1
    prefix: '--fastqc_args'
    valueFrom: '"--noextract"'
  - position: 1
    prefix: '--gzip'
  - position: 1
    valueFrom: |
      ${
        if ( inputs.adapter1 == "illumina" ){ return "--illumina" }
        else if ( inputs.adapter1 == "nextera" ){ return "--nextera" }
        else if ( inputs.adapter1 == "small_rna" ){ return "--small_rna" }
        else { return null }
      }
  - position: 1
    prefix: '--adapter'
    valueFrom: |
      ${
        if ( inputs.adapter1 != null && inputs.adapter1 != "illumina" && inputs.adapter1 != "nextera" && inputs.adapter1 != "small_rna" ){
          return inputs.adapter1
        } else {
          return null
        }
      }
  - position: 1
    prefix: '--adapter2'
    valueFrom: |
      ${
        if ( inputs.fastq2 != null && inputs.adapter2 != null && inputs.adapter1 != "illumina" && inputs.adapter1 != "nextera" && inputs.adapter1 != "small_rna" ){
          return inputs.adapter2
        } else {
          return null
        }
      }
  - position: 1
    valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return "--paired" }
      }
  - position: 1
    valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return "--retain_unpaired" }
      }
  - position: 1
    prefix: '--length_1'
    valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return inputs.min_unpaired_read_rescue_length }
      }
  - position: 1
    prefix: '--length_2'
    valueFrom: |
      ${
        if ( inputs.fastq2 == null ){ return null }
        else { return inputs.min_unpaired_read_rescue_length }
      }
requirements:
  - class: ResourceRequirement
    ramMin: 7000
    coresMin: 6
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7'
