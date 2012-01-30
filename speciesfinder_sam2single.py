#!/panvol1/simon/bin/python

import argparse
import sys

def revcomp(seq):
   '''Reverse complement sequence'''
   
   comp = []
   for s in seq:
      if s == 'A': comp.append('T')
      if s == 'T': comp.append('A')
      if s == 'C': comp.append('G')
      if s == 'G': comp.append('C')
      if s == 'N': comp.append('N')
      
   return comp[::-1]

def extract_single(i, prefix):
   '''Extract pairs from sam file (unfiltered)'''
   
   if i == '-':
      fh = sys.stdin
   else:
      fh = open(i, 'r')
   
   fh_out3 = open(prefix+'.single.fq', 'w')
   
   # parse
   for line in fh:
      if line.startswith('@'): continue
      
      line = line.rstrip()
      fields = line.split('\t')
      if fields[1] == '':
         end = '1'
      else:
         end = fields[1][-1]
      
      # reverse read
      if fields[1].find('r') > -1:
         seq = revcomp(fields[9])
         q = fields[10][::-1]
         fh_out3.write('@%s/%s\n%s\n+\n%s\n' % (fields[0], end, ''.join(seq), q))
      else:
         fh_out3.write('@%s/%s\n%s\n+\n%s\n' % (fields[0], end, fields[9], fields[10]))
   
   fh_out3.close()

if __name__ == '__main__':
   
   # create the parser
   parser = argparse.ArgumentParser(prog='cluster_sam2single.py', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=50, width=130), usage='%(prog)s [options]', description='''Extract reads from sam (stdin) and write .single.fq. Sam must be viewed with -X''')
   
   parser.add_argument('--i', help='input [-]', default='-')
   parser.add_argument('--prefix', help='output prefix [extracted]', default='extracted')
   
   args = parser.parse_args()
   #args = parser.parse_args('--i Campy-05-0121.sam --prefix Campy-05-0121.bam'.split())
   
   extract_single(args.i, args.prefix)
