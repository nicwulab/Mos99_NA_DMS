#!/usr/bin/python
#Required dssp_parser: http://openwetware.org/wiki/Wilke:ParseDSSP
import os
import sys
import glob
from collections import defaultdict
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

def read_fitfile(file_file):
  infile = open(file_file, 'r')
  fit_dict = {}
  for line in infile.readlines():
    if 'pos' in line: continue
    pos, fit = line.rstrip().rsplit()
    fit_dict[pos] = float(fit)
  infile.close()
  return (fit_dict)

def read_sitefile(site_file):
  infile = open(site_file, 'r')
  site_dict = {}
  for line in infile.readlines():
    if 'pos' in line: continue
    site, positions = line.rstrip().rsplit("\t")
    site_dict[site] = list(map(int, positions.rsplit(',')))
  return site_dict

def reading_ASA(file_asa):
  infile = open(file_asa,'r')
  dict_asa = {}
  for line in infile.readlines():
    if 'aa' in line: continue
    line = line.rstrip().rsplit("\t")
    aa  = line[0]
    asa = line[3]
    dict_asa[aa] = asa
  return dict_asa

def parse_dssp(dssp, dict_asa, chainID, positions):
  RSA_dict = {}
  for pos in positions:
    dssp_info = list(dssp[chainID,(' ', pos, ' ')])
    aa = dssp_info[0]
    SA = dssp_info[2]
    RSA = float(SA)/float(dict_asa[aa])
    RSA_dict[pos] = RSA
  return RSA_dict

def classify_pos_type(pos, RSA_monomer, delta_RSA, site_dict):
  pos_type = []
  for site in site_dict.keys():
    if int(pos) in site_dict[site]:
      pos_type.append(site)
  if len(pos_type) == 0:
    if RSA_monomer < 0.05:
      pos_type.append('buried')
    elif delta_RSA > 0.5*RSA_monomer:
      pos_type.append('interface')
  if len(pos_type) > 1:
    print ('position %s belongs to more than 1 classification (%s)' % (pos, ', '.join(pos_type)))
  if len(pos_type) == 0:
    pos_type.append('exposed')
  pos_type = '-'.join(pos_type)
  return pos_type

def compile_out(outfile, RSA_dict_monomer, RSA_dict_tetramer, fit_dict, positions, site_dict):
  print ('writing: %s' % outfile)
  outfile = open(outfile,'w')
  outfile.write("\t".join(['pos', 'RSA_monomer', 'RSA_tetramer', 'delta_RSA', 'type', 'fit'])+"\n")
  for pos in positions:
    RSA_monomer  = RSA_dict_monomer[pos]
    RSA_tetramer = RSA_dict_tetramer[pos]
    delta_RSA    = RSA_monomer-RSA_tetramer
    fit          = fit_dict[str(pos)]
    pos_type     = classify_pos_type(pos, RSA_monomer, delta_RSA, site_dict)
    outfile.write("\t".join(map(str, [pos, RSA_monomer, RSA_tetramer, delta_RSA, pos_type, fit]))+"\n")
  outfile.close()

def main():
  outfile   = 'result/position_type_vs_fit.tsv'
  dssp_monomer  = dssp_dict_from_pdb_file('PDB/N2_NA_monomer.pdb', DSSP='mkdssp')[0]
  dssp_tetramer = dssp_dict_from_pdb_file('PDB/N2_NA_tetramer.pdb', DSSP='mkdssp')[0]
  file_asa  = 'data/ASA.table'
  dict_asa  = reading_ASA(file_asa)
  fit_file  = 'result/Mos99_mean_mut_fit.tsv'
  fit_dict  = read_fitfile(fit_file)
  site_file = 'data/sites_info.tsv'
  site_dict = read_sitefile(site_file)
  chainID   = 'A'
  positions = list(range(82, 466))
  RSA_dict_monomer  = parse_dssp(dssp_monomer, dict_asa, chainID, positions)
  RSA_dict_tetramer = parse_dssp(dssp_tetramer, dict_asa, chainID, positions)
  compile_out(outfile, RSA_dict_monomer, RSA_dict_tetramer, fit_dict, positions, site_dict)

if __name__ == "__main__":
  main()
