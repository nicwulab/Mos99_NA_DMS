#!/usr/bin/python
import sys
import pandas as pd

def main():
  outfile = 'result/Mos99_mean_mut_fit.tsv'
  df = pd.read_csv('result/Mos99_fit.csv')
  df = df.rename(columns={'0': 'pos'})
  df = df[~df['Mutation'].str.contains('_')]
  df = df[~df['Mutation'].str.contains('silent')]
  df = df[['pos','Fitness']]
  df = df.groupby(['pos']).mean()
  print ('writing: %s' % outfile)
  df.to_csv(outfile, sep="\t")

if __name__ == "__main__":
  main()
