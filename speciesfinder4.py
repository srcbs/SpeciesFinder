#!/tools/opt/python/python2.7.2/bin/python

import argparse
import subprocess
import sys
import os

paths = {}
paths['bwa_home'] = '/panfs1/cge-servers/SpeciesFinder/scripts/'
paths['samtools_home'] = '/panfs1/cge-servers/SpeciesFinder/scripts/'
paths['trinity_home'] = '/panfs1/cge-servers/SpeciesFinder/scripts/trinity/'
paths['sam2single'] = '/panfs1/cge-servers/SpeciesFinder/scripts/speciesfinder_sam2single.py'
paths['blastall'] = '/panfs1/cge-servers/SpeciesFinder/scripts/blastall'
paths['nc_tax_name'] = '/panfs1/cge-servers/SpeciesFinder/scripts/db/nc_tax_name.tab'
paths['taxonomy_ncbi'] = '/panfs1/cge-servers/SpeciesFinder/scripts/db/taxonomy/'
paths['R'] = '/tools/bin/R-2.12'
paths['speciesfinder_home'] = '/panfs1/cge-servers/SpeciesFinder/scripts/'

def set_abspath():
   '''Returns absolute path of file to argparse'''
   class SetAbspath(argparse.Action):
      def __call__(self, parser, args, filenames, option_string=None):
         import os
         if filenames == '-':
            setattr(args, self.dest, '-')
         elif type(filenames) == str:
            f_abs = os.path.abspath(filenames)
            setattr(args, self.dest, f_abs)
         elif type(filenames) == list:
            new_list = []
            for f in filenames:
               new_list.append(os.path.abspath(f))
            setattr(args, self.dest, new_list)
         else:
            setattr(args, self.dest, filenames)
   return SetAbspath

def check_norecords(args):
   '''Check the number of records in input fasta'''
   
   cmd = 'grep -c ">" %s' % (args.i)
   #cmd = 'grep -c ">" %s' % (i)
   p_grep = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   stdout,stderr = p_grep.communicate()
   records = int(stdout.rstrip())
   return records

def uniqueList(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def average(values):
    '''Computes the arithmetic mean of a list of numbers.'''
    return sum(values, 0.0) / len(values)

def get_readlengths(i):
   '''Count readlengths from first 10 reads'''
   
   read_lengths=[]
   c = 0
   fh = open(i, 'r')
   for n, line in enumerate(fh):
      pos = n%4
      line = line.rstrip()
      if pos==1: 
         read_lengths.append(len(line))
         c = c +1
      if c > 100: break
   avg = average(read_lengths)
   return int(avg)
   

def assembly_trinity(args):
   '''Performs bwasw alignment of input reads vs. pre-indexed fasta and then trinity assembly on resulting bam-file'''
   
   if args.i == '-': bam = '%s.bam' % args.sample
   else: bam = os.path.split(args.i)[1] + '.bam'
   fq = bam + '.single.fq'
   
   # alignment #
   sys.stderr.write('Running alignment\n')
   cmd = paths['bwa_home'] + 'bwa'
   arg = ' bwasw -t %s %s %s' % (args.n, args.db, args.i)
   bwa_call = cmd+arg
   
   cmd = paths['samtools_home'] + 'samtools'
   arg = ' view -Sb -F 4 - > %s' % bam
   sam_call = cmd+arg
   
   if args.Mbases == 0:
      # use all reads
      call = '%s | %s' % (bwa_call, sam_call)
   else:
      # check that the actual no. of lines is there in the file #
      # calc read length and determine no. of reads to get X Mbases
      avg = get_readlengths(args.i)
      
      # calculate the number of reads used to get X Mbases
      bases = args.Mbases*1e6
      reads_no = int(float(bases)/avg)
      
      # only take first x reads
      bwa_call = '%sbwa bwasw -t %s %s - '% (paths['bwa_home'], args.n, args.db)
      call = 'head -n %i %s | %s | %s' % (reads_no*4, args.i, bwa_call, sam_call)
   
   print call
   #p_bwasw = subprocess.Popen(call, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   p_bwasw = subprocess.call(call, shell=True)
   
   # extract fq
   cmd = paths['samtools_home'] + 'samtools'
   arg = ' view -X -F 4 %s | %s --prefix %s ' % (bam, paths['sam2single'], bam)
   #arg =  ''' view %s | perl -ane 'print "@", $F[0], "\n", $F[9], "\n+\n", $F[10], "\n"' > %s ''' % (bam, fq)
   print cmd+arg
   p_fq = subprocess.call(cmd+arg, shell=True)
   
   # run Trinity
   cmd = paths['trinity_home'] + 'Trinity.pl'
   arg = ' --seqType fq --single %s --run_butterfly --CPU %s --output trinity' % (fq, args.n)
   print cmd+arg
   p_trinity = subprocess.Popen(cmd+arg, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   stdout,stderr = p_trinity.communicate()
   #print stdout
   print stderr
   
   return 'trinity/Trinity.fasta'


def assembly_map(args):
   '''Performs bwasw alignment of input contigs'''
   
   if args.i == '-': bam = '%s.bam' % args.sample
   else: bam = os.path.split(args.i)[1] + '.bam'
   fq = bam + '.single.fq'
   
   # alignment #
   sys.stderr.write('Running alignment\n')
   cmd = paths['bwa_home'] + 'bwa'
   arg = ' bwasw -t %s %s %s' % (args.n, args.db, args.i)
   bwa_call = cmd+arg
   
   cmd = paths['samtools_home'] + 'samtools'
   arg = ' view -Sb -F 4 - > %s' % bam
   sam_call = cmd+arg
   
   call = '%s | %s' % (bwa_call, sam_call)
   print call
   #p_bwasw = subprocess.Popen(call, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   p_bwasw = subprocess.call(call, shell=True)
   
   # extract the hit #
   hits = 'ssu.fa'
   cmd = paths['samtools_home'] + 'samtools'
   arg = ''' view %s | perl -ane ' print ">$F[0]\n$F[9]\n"' > %s ''' % (bam, hits)
   p_extract = subprocess.call(cmd+arg, shell=True)
   
   return hits


def blast(hits_fa, fa, n):
   '''Run blast and parse output'''
   
   sys.stderr.write('Running blast\n')
   
   # run blast
   call = '%s -p blastn -d %s -a %s -m8 -i %s -o ssu.blast' % (paths['blastall'], fa, n, hits_fa)
   p = subprocess.Popen(call, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   stdout, stderr = p.communicate()
   
   # open nuc mapping
   nuc2genome = {}
   nuc2tax = {}
   fh = open(paths['nc_tax_name'], 'r')
   for line in fh:
      line = line.rstrip()
      fields = line.split('\t')
      nuc2genome[fields[0]] = fields[2]
      nuc2tax[fields[0]] = fields[1]
   
   fh.close()
   
   # parse blast output
   # query_id, subject_id, identity, aln_length, nmismatches, ngapopenings, query_start, query_end, subject_start, subject_end, evalue, bitscore
   
   # get ten best hit per contig
   fh = open('ssu.blast', 'r')
   c = 0
   hits = {}
   prev_genome = ''
   for line in fh:
      line = line.rstrip()
      fields = line.split()
      query_id, subject_id, identity, aln_length, nmismatches, ngapopenings, query_start, query_end, subject_start, subject_end, evalue, bitscore = fields
      
      if hits.has_key(query_id):
         if fields[1] == prev_genome:
            pass
         else:
            prev_genome = fields[1]
            c = c + 1
            if c < 50:
               hits[query_id].append(fields[1:])
      else:
         c = 0
         hits[query_id] = [fields[1:]]
         prev_genome = fields[1]
   
   if len(hits) > 1: sys.stderr.write('More than 1 transcript were created (%i), perhaps the ssu is divided into subunits\n' % len(hits))
   
   # look up contigs in nc_tax_name
   lines = []
   c = 0
   for key in hits.keys():
      for v in hits[key]:
         c = c + 1
         status=""
         nc = v[0].split('|')[1]
         rrna_int=v[0].split('|')[2][1:]
         rrna_length = int(rrna_int.split('-')[1]) - int(rrna_int.split('-')[0])
         if rrna_length < 0: rrna_length = -1*rna_length
         coverage = float(v[2])/(rrna_length+1)
         try: 
            name = nuc2genome[nc]
         except:
            name = "NA"
            sys.stderr.write('NC was not in nuc2genome: %s\n' % nc)
         try:
            tax = nuc2tax[nc]
         except:
            tax = "NA"
            sys.stderr.write('NC was not in nuc2tax: %s\n' % nc)
         
         #if c == 1: besttax = tax
         if args.contigs:
            if (coverage > 0.995 and float(v[1]) > 99.5): status="PASS"
            else: status="FAIL"
         else:
            if (coverage > 0.98 and float(v[1]) > 98.0): status="PASS"
            else: status="FAIL"
         
         lines.append([str(c), status, tax, name, nc, key, str(coverage*100), v[1], v[2], v[3], v[4], v[9], v[10]])
   
   # sort lines based sum of coverage and identity # (actually sorting is not need when R-script is used) #
   # make sum
   s = []
   for l in lines:
      s.append(float(l[6])+float(l[7]))
   
   # zip and sort
   lines_zip = zip(s, lines)
   lines_zip.sort(reverse=True)
   lines_sorted = zip(*lines_zip)[1]
   
   # write out
   fh_out = open('ssu.out', 'w')
   fh_out.write('Hit#\tStatus\tTaxid\tName\tNC\tTranscript\tCoverage\tIdentity\tAln_length\tnmismatches\tngapopenings\tevalue\tbitscore\n')
   for line in lines_sorted:
      fh_out.write('%s\n' % '\t'.join(line))
   fh_out.close()
   
   # sort instead using R
   cmd = '%s --vanilla ssu.out < %sspeciesfinder_parse.R' % (paths['R'], paths['speciesfinder_home'])
   exitcode = subprocess.call(cmd, shell=True)
   
   # read in best tax from ssu.out.tab
   fh = open('ssu.out.tab', 'r')
   header = fh.readline()
   line = fh.readline()
   fields = line.split('\t')
   besttax = fields[2]
   
   return besttax

def parse_taxonomy(taxid):
   '''Traverse taxonomy'''
   
   if taxid == '' or taxid == 'NA': return None
   
   # create lineage
   lineage = {}
   tax2family_handle = open(paths['taxonomy_ncbi'] + 'nodes.dmp', 'r')
   for line in tax2family_handle:
      s = line.split('\t|\t')
      lineage[s[0]] = (s[1], s[2])
   
   # parse through to find lineage
   taxonomy = [taxid]
   full_tax = []
   tax = taxid
   while 1:
      tax_tuple = lineage[tax]
      tax = tax_tuple[0]
      if tax == '1': break
      taxonomy.append(tax)
      full_tax.append(tax_tuple)
      
   # get names of the taxids
   tax_dict = {}
   for t in taxonomy:
      tax_dict[t] = 1
   
   fh = open(paths['taxonomy_ncbi'] + 'names.dmp', 'r')
   for line in fh:
      fields = line.split('\t|\t')
      if fields[3][:-3] == 'scientific name':
         if tax_dict.has_key(fields[0]):
            tax_dict[fields[0]] = fields[1]
   
   # order and print
   from collections import OrderedDict
   
   od = OrderedDict( (k,tax_dict[k]) for k in taxonomy)
   fh_out = open('ssu.lineage', 'w')
   fh_out.write('%s\n' % (' --> '.join(od.values()[::-1][1:])))
   fh_out.write('%s\n' % (' --> '.join(od.keys()[::-1][1:])))
   fh_out.close()


def main(args):
   '''Main'''
   
   sys.stderr.write('Processing file: %s\n' % args.i)
   
   # move to directory
   if not os.path.exists(args.sample):
      os.makedirs(args.sample)
   
   os.chdir(args.sample)
   
   # assembly
   if args.contigs == False:
      hits_fa = assembly_trinity(args)
   else:
      records = check_norecords(args)
      if records > 10: 
         #sys.stderr.write('Skipping BWA\n')
         #hits_fa = assembly_map(args)
         hits_fa = args.i
      else: hits_fa = args.i
   
   # 16S rrna (SSU)
   sys.stderr.write('Running SSU\n')
   besttax = blast(hits_fa, args.db, args.n)
   
   if args.tax: parse_taxonomy(besttax)
   
   sys.stderr.write('Done\n')


if __name__ == '__main__':
   parser = argparse.ArgumentParser(prog='speciesfinder.py', formatter_class=lambda prog: argparse.HelpFormatter(prog,max_help_position=35, width=110),
      description='''Perform SW alignment of reads (fastq) or contigs vs database (alleles). Performs local assembly of the mapped reads or extract mapped contigs. Then runs blast vs input database''',
      usage='%(prog)s in.fastq|contigs.fa [options]')
   
   parser.add_argument('i', help='input fastq file, contigs or - for stdin [-]', default='-', action=set_abspath())
   parser.add_argument('--db', help='input reference alleles database [db/ssu.all.filt.fa]', action=set_abspath(), default='/panfs1/cge-servers/SpeciesFinder/scripts/db/ssu.all.filt.fa')
   parser.add_argument('--contigs', help='input is not reads but contigs [False]', default=False, action='store_true')
   parser.add_argument('--sample', help='sample name and output dir [fmap/]', default='fmap/')
   parser.add_argument('--Mbases', help='Mbases to use for assembly [0=all]', default=0, type=float)
   parser.add_argument('--n', help='number of threads for bwasw [4]', default=4, type=int)
   parser.add_argument('--tax', help='parse taxonomy [False]', default=False, action='store_true')
   
   args = parser.parse_args()
   #args = parser.parse_args('Campy-05-0121_7_1_sequence.txt ../16SRNA/wgs_ssu.fa --sample Campy-05-0121 --ssu --reads_no 10000'.split())
   #args = parser.parse_args('Campy-05-0121/assembly/contigs.fa ../16SRNA/wgs_ssu.fa --sample Campy-05-0121 --ssu --contigs'.split())
   #args = parser.parse_args('/panfs1/cge/people/salvocos/MLST/complete_genome_results/Abaumannii/Abaumannii_FSA/CP000521.fsa --sample Abaumannii --contigs'.split())
   
   main(args)
