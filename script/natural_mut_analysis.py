#!/usr/bin/python
import sys
import numpy as np
from Bio import SeqIO
from collections import defaultdict, Counter

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def Read_aln(alnfile):
  aa = set(['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W'])
  records   = [record for record in SeqIO.parse(alnfile,"fasta")]
  alndict = defaultdict(list)
  for record in records:
    header = str(record.id)
    if header.count('|')!=5: continue
    ID    = header.rsplit('|')[0]
    PSG   = header.rsplit('|')[1]
    year  = header.rsplit('|')[-1][0:4]
    seq   = str(record.seq)
    if not set(seq).issubset(aa): continue
    assert(isInt(year))
    alndict[year].append(seq)
  return alndict

def reading_fit_file(filename):
  infile   = open(filename, 'r')
  fit_dict = {}
  for line in infile.readlines():
    if 'Mutation' in line: continue
    line = line.rsplit(",")
    mut  = line[1]
    fit  = float(line[-2])
    if 'WT' in mut or '-' in mut or 'Amp' in mut: continue
    fit_dict[mut] = fit
  return fit_dict

def seqs2consensus(seqlist):
  consensus = ''
  for n in range(len(seqlist[0])):
    resi = []
    for seq in seqlist:
      resi.append(seq[n])
    most_common,num_most_common = Counter(resi).most_common(1)[0]
    consensus+=most_common
  return consensus

def call_muts(con_seq, ref_seq):
  assert(len(con_seq)==len(ref_seq))
  muts = [ref_aa+str(pos+1)+con_aa for con_aa, ref_aa, pos in zip(con_seq, ref_seq, range(len(ref_seq))) if con_aa != ref_aa]
  return (muts)

def Consensus_mut(alndict, outfile, ref_seq):
  print("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write('year'+"\t"+'mut'+"\n")
  for year in sorted(map(int, alndict.keys())):
    seqs     = alndict[str(year)]
    con_seq  = seqs2consensus(seqs)
    muts     = call_muts(con_seq, ref_seq)
    outfile.write(str(year)+"\t"+str(",".join(muts))+"\n")
  outfile.close()

def Compute_freq(alndict, fit_dict, ref_seq, outfile):
  print ("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['mut', 'year', 'freq', 'fit'])+"\n")
  mut_dict_all = defaultdict(float)
  for year in sorted(alndict.keys()):
    total_seq = len(alndict[year])
    mut_dict_year = defaultdict(int)
    print ('processing year: %s (%i seqs)' % (year, total_seq))
    for n in range(len(ref_seq)):
      pos  = str(n + 1)
      WT_aa = ref_seq[n]
      for nat_seq in alndict[year]:
        if WT_aa != nat_seq[n]:
          mut_dict_year[WT_aa+pos+nat_seq[n]] += 1
    for mut in mut_dict_year.keys():
      if mut not in fit_dict.keys(): continue
      fit  = fit_dict[mut]
      freq = float(mut_dict_year[mut])/float(total_seq)
      outfile.write("\t".join(map(str,[mut, year, freq, fit]))+"\n")
  outfile.close()

def fitness_inference(alndict, fit_dict, ref_seq, outfile_con, outfile_all):
  print ("writing: %s" % outfile_con)
  print ("writing: %s" % outfile_all)
  outfile_con = open(outfile_con, 'w')
  outfile_all = open(outfile_all, 'w')
  outfile_con.write("\t".join(['year','fit'])+"\n")
  outfile_all.write("\t".join(['year','fit'])+"\n")
  for year in sorted(alndict.keys()):
    seqs     = alndict[str(year)]
    con_seq  = seqs2consensus(seqs)
    fit = 0
    fits = []
    muts = call_muts(con_seq, ref_seq)
    for mut in muts:
      if mut in fit_dict.keys():
        fit += float(fit_dict[mut])
        fits.append(float(fit_dict[mut]))
    if len(fits) == 0: fits.append(0)
    mean_fit = np.mean(fits)
    outfile_con.write("\t".join(map(str,[year, np.mean(fits)]))+"\n")
    #outfile_con.write("\t".join(map(str,[year, fit]))+"\n")

    for seq in seqs:
      fit = 0
      fits = []
      muts = call_muts(seq, ref_seq)
      for mut in muts:
        if mut in fit_dict.keys():
          fit += float(fit_dict[mut])
          fits.append(float(fit_dict[mut]))
      if len(fits) == 0: fits.append(0)
      mean_fit = np.mean(fits)
      outfile_all.write("\t".join(map(str,[year, np.mean(fits)]))+"\n")
      #outfile_all.write("\t".join(map(str,[year, fit]))+"\n")
  outfile_con.close()
  outfile_all.close()

def main():
  ref_seq    = open('Fasta/Mos99_NA.pep','r').readlines()[1].strip()
  alnfile    = 'Fasta/Human_H3N2_NA_2020.aln'
  fitfile    = 'result/Mos99_fit.csv'
  #outfile1   = 'result/N2_mutation_year_from_Mos99.tsv'
  outfile2   = 'result/N2_mutation_freq.tsv'
  #outfile3   = 'result/N2_natural_seq_con_fit.tsv'
  #outfile4   = 'result/N2_natural_seq_all_fit.tsv'
  alndict    = Read_aln(alnfile)
  fit_dict   = reading_fit_file(fitfile)
  #Consensus_mut(alndict, outfile1, ref_seq)
  Compute_freq(alndict,fit_dict,ref_seq,outfile2)
  #fitness_inference(alndict, fit_dict, ref_seq, outfile3, outfile4)

if __name__ == "__main__":
    main()
