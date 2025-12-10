import pysam


# open vcf file
vcf = pysam.VariantFile("GSE121668_vcf_snp_fltr/SRR8096262.hg38.identified.snp.fltr.vcf")
# open fasta file
genome = pysam.FastaFile("Homo_sapiens.GRCh38.dna.primary_assembly.fa")
# define by how many bases the variant should be flanked
flank = 50

# iterate over each variant
for record in vcf:
  # extract sequence
  #
  # The start position is calculated by subtract the number of bases
  # given by 'flank' from the variant position. The position in the vcf file
  # is 1-based. pysam's fetch() expected 0-base coordinate. That's why we
  # need to subtract on more base.
  #
  # The end position is calculated by adding the number of bases
  # given by 'flank' to the variant position. We also need to add the length
  # of the REF value and subtract again 1 due to the 0-based/1-based thing.
  #
  # Now we have the complete sequence like this:
  # [number of bases given by flank]+REF+[number of bases given by flank]
  start_pos = record.pos-1-flank
  if(start_pos < 0):
    continue
  seq = genome.fetch(record.chrom, start_pos, record.pos-1+len(record.ref)+flank)

  # print out tab seperated columns:
  # CRHOM, POS, REF, ALT, flanking sequencing with variant given in the format '[REF/ALT]'
  with open('variant_region_seq/SRR8096262.hg38.identified.snp.fltr.all.fa', 'a') as f:
    f.write('>'+str(record.chrom)+'_'+str(record.pos)+'_'+str(record.ref)+'_'+str(record.alts[0]) + '\n' + '{}{}{}'.format(seq[:flank], record.ref, seq[flank+len(record.ref):]) + '\n')
  with open('variant_region_seq/SRR8096262.hg38.identified.snp.fltr.' + str(record.ref) + 'to' + str(record.alts[0]) + '.fa', 'a') as f:
    f.write('>'+str(record.chrom)+'_'+str(record.pos)+'_'+str(record.ref)+'_'+str(record.alts[0]) + '\n' + '{}{}{}'.format(seq[:flank], record.ref, seq[flank+len(record.ref):]) + '\n')
