#!/usr/bin/python
import os
import sys
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import defaultdict, Counter

def extract_motif(alnfile, positions, length_filter):
  aa = set(['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','-'])
  motifs = []
  count_seq = 0
  for record in SeqIO.parse(alnfile,"fasta"):
    ID    = str(record.id)
    seq   = str(record.seq)
    if not set(seq).issubset(aa): continue
    if length_filter != 'NA':
      if len(seq) != length_filter: continue
    count_seq += 1
    motif = ''
    for pos in positions:
      motif += seq[pos]
    motifs.append(motif)
  print (count_seq)
  return motifs

def make_sequence_logo(sequence_list, figname):
  CDRH3_len = len(sequence_list[0])
  logo_width = (CDRH3_len-0)*1
  logo_width = 11 if logo_width > 11 else logo_width
  fig, ax = plt.subplots(1,1,figsize=[logo_width,2])
  seqlogo_matrix = logomaker.alignment_to_matrix(sequence_list)
  seqlogo = logomaker.Logo(seqlogo_matrix, font_name="Arial", color_scheme="weblogo_protein", width=1, ax=ax)
  seqlogo.style_spines(visible=False)
  seqlogo.ax.set_xticks([])
  seqlogo.ax.set_yticks([])
  seqlogo.fig.tight_layout()
  plt.savefig(figname, dpi=150)
  plt.close()
  print('Written %s' % figname, file = sys.stdout)

def main():
  alnfile_N2       = 'Fasta/Human_H3N2_NA_2020.aln'
  alnfile_avian_N2 = 'Fasta/Avian_N2_NA.aln'
  alnfile_subtypes = 'Fasta/NA_subtypes.aln'
  N2_motifs  = extract_motif(alnfile_N2, [282, 287, 303, 354, 382], 469)
  avian_N2_motifs  = extract_motif(alnfile_avian_N2, [302, 308, 325, 377, 405], 'NA')
  all_motifs = extract_motif(alnfile_subtypes, [299,304,320,374,402], 'NA')
  make_sequence_logo(N2_motifs, 'graph/seqlogo_cluster2_N2.png')
  make_sequence_logo(avian_N2_motifs, 'graph/seqlogo_cluster2_avian_N2.png')
  make_sequence_logo(all_motifs, 'graph/seqlogo_cluster2_all.png')

if __name__ == "__main__":
  main()
