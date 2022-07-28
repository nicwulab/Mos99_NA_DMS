#!/usr/bin/python
import sys
from prody import *
from scipy.spatial import distance 

def extract_resi_centroid(p):
  centroid_dict = {}
  for chain in sorted(list(set(p.getChids()))):
    p_chain = p.select('chain '+chain)
    for resnum in sorted(list(set(p_chain.getResnums()))):
      res = p_chain.select('resnum '+str(resnum))
      centroid = (calcCenter(res))
      centroid_dict[chain+'-'+str(resnum)] = centroid
  return (centroid_dict)

def min_distance(A_centroid_dict, B_centroid_dict):
  A_min_dist_to_B = {}
  for A_res in A_centroid_dict.keys():
    A_res_centroid = A_centroid_dict[A_res]
    min_dist = 9999
    for B_res in B_centroid_dict.keys():
      B_res_centroid = B_centroid_dict[B_res]
      dist = distance.euclidean(A_res_centroid, B_res_centroid)
      if dist < min_dist: min_dist = dist
    A_min_dist_to_B[A_res] = min_dist
  return (A_min_dist_to_B)

def writing_dist_file(outfile, NA_min_dist_to_SIA_dict):
  print ("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['pos', 'dist_to_SIA'])+"\n")
  NA_residues = list(NA_min_dist_to_SIA_dict.keys())
  NA_residues = sorted(list(map(lambda x:int(x.rsplit('-')[1]), NA_residues)))
  for NA_resi in NA_residues:
    if 82 <= NA_resi <= 465:
      dist_to_SIA = NA_min_dist_to_SIA_dict['A-'+str(NA_resi)]
      outfile.write("\t".join(map(str,[NA_resi, dist_to_SIA]))+"\n")
  outfile.close()

def main():
  outfile = 'result/Dist_to_active_site.tsv'
  p = parsePDB('PDB/Mos99_WT_sialic_acid_final.pdb')
  SIA = p.select('chain A resnum 803')
  NA  = p.select('chain A resnum 82 to 469')
  SIA_centroid_dict = extract_resi_centroid(SIA)
  NA_centroid_dict = extract_resi_centroid(NA)
  NA_min_dist_to_SIA_dict = min_distance(NA_centroid_dict, SIA_centroid_dict)
  writing_dist_file(outfile, NA_min_dist_to_SIA_dict)

if __name__ == "__main__":
  main()
