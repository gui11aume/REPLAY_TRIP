#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

from math import sqrt

csize = {
   '2L': 23011544,
   '2R': 21146708,
   '3L': 24543557,
   '3R': 27905053,
   'X':  22422827,
   '4':  1351857,
}

wsz = 2000

def main():
   # Create a dictionary of vectors.
   Vterm = dict()
   Vprom = dict()
   for chrom in csize:
      Vterm[chrom] = [0] * (1 + csize[chrom] / wsz)
      Vprom[chrom] = [0] * (1 + csize[chrom] / wsz)
   # Read coordinates of active genes.
   with open('act_genes_r5.57.txt') as f:
      for line in f:
         if line[0] == '#' or line.startswith('seqname'):
            continue
         items = line.split()
         # Find promoters / terminators depending on the orientation.
         (p,t) = (int(items[1]), int(items[2])) if items[4] == '+' else \
               (int(items[2]), int(items[1]))
         Vprom[items[0]][p / wsz] += 1
         Vterm[items[0]][t / wsz] += 1

   sys.stderr.write('parsing allprom.txt\n')
   CHROM = dict()
   NEXP = dict()
   # Replicates.
   VECT1 = dict()
   VECT2 = dict()
   with open('../allprom.txt') as f:
      discard = next(f)
      for line in f:
         bcd,chrom,strand,pos,nexp = line.split()[:5]
         CHROM[bcd] = chrom
         NEXP[bcd] = nexp
         VECT1[bcd] = [0] * (1 + csize[chrom] / wsz)
         VECT2[bcd] = [0] * (1 + csize[chrom] / wsz)

   # Convenience function to fill the dictionaries of vectors.
   def update_dict(fname, VDICT):
      sys.stderr.write('parsing %s\n' % fname)
      with open(fname) as f:
         for line in f:
            bcd,nreads,chrom,strand,pos = line.split()
            try:
               VDICT[bcd][int(pos) / wsz] += int(nreads)
            except (KeyError, IndexError):
               continue

   update_dict('/data/91_4C_12679_ACTTGA.fastq.gz_final.tsv', VECT1)
   update_dict('/data/91_4C_rep_15204_ATCACG.fastq.gz_final.tsv', VECT2)
   
   update_dict('/data/92_4C_12680_GATCAG.fastq.gz_final.tsv', VECT1)
   update_dict('/data/92_4C_rep_15205_CGATGT.fastq.gz_final.tsv', VECT2)
   
   update_dict('/data/93_4C_12681_TAGCTT.fastq.gz_final.tsv', VECT1)
   update_dict('/data/93_4C_rep_15206_TTAGGC.fastq.gz_final.tsv', VECT2)
   
   update_dict('/data/94_4C_12682_GGCTAC.fastq.gz_final.tsv', VECT1)
   update_dict('/data/94_4C_rep_15207_TGACCA.fastq.gz_final.tsv', VECT2)
   
   update_dict('/data/95_4C_12683_CTTGTA.fastq.gz_final.tsv', VECT1)
   update_dict('/data/95_4C_rep_15208_ACAGTG.fastq.gz_final.tsv', VECT2)

   def normalize(V1, V2):
      # Copy first vector.
      prob = V1[:]
      for i in range(len(prob)):
         prob[i] = sqrt(prob[i] * V2[i])
      S = sum(prob)
      for i in range(len(prob)):
         prob[i] /= S
      return prob

   sys.stdout.write('bcd\tnexp\tprom\tterm\n')
   for bcd in VECT1:
      if bcd not in VECT2: continue
      V1 = VECT1[bcd]
      V2 = VECT2[bcd]
      # Use only the barcodes with a minimum number of contacts.
      n_contacts_1 = sum([1 for a in V1 if a > 0])
      n_contacts_2 = sum([1 for a in V2 if a > 0])
      if n_contacts_1 < 20 or n_contacts_2 < 20: continue
      pp = normalize(V1, V2)

      pt = pp[:]
      for i in range(len(pp)):
         pp[i] *= Vprom[CHROM[bcd]][i]
         pt[i] *= Vterm[CHROM[bcd]][i]
      Tp = sum(pp)
      Tt = sum(pt)
      sys.stdout.write('%s\t%s\t%.3f\t%.3f\n' % (bcd, NEXP[bcd], Tp, Tt))

if __name__ == '__main__':
   main()
