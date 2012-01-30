#!/panvol1/simon/bin/python2.7

import subprocess
import argparse
import logging
import pipelinemod

def bwa_index(fa):
   '''Index fasta file for bwa alignment'''
   
   import subprocess
   import pipelinemod
   paths = pipelinemod.setSystem()
   
   cmd = '/panvol1/simon/bin/bwa-0.5.9/bwa'
   #cmd = paths['bwa_home'] + 'bwa'
   arg = ' index %s' % fa
   call = cmd+arg
   logger.info(call)
   subprocess.check_call(call, shell=True)

def create_genomefile(fa, genomefile):
   '''Create genome file from fasta input'''
   
   from Bio import SeqIO
   
   outhandle = open(genomefile, 'w')
   for record in SeqIO.parse(fa, 'fasta'):
      outhandle.write('%s\t%i\n' % (record.id, len(record)))
   
   outhandle.close()


parser = argparse.ArgumentParser(description=
   '''Prepare files for FMAP-MLST. Does bwa index on the fasta files and
   creates genome files for coverage predictions
   '''
)

parser.add_argument('--fa', help='allele fasta files', nargs='+', action='append')
parser.add_argument('--log', help='log level [info]', default='info')

args = parser.parse_args()

# set logger
logger = logging.getLogger('mlst_fmap_prepare.py')
hdlr = logging.FileHandler('mlst_fmap_prepare.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
if args.log == 'info':
   logger.setLevel(logging.INFO)

# run
for fa in args.fa[0]:
   bwa_index(fa)
   create_genomefile(fa, fa+'.genome')

