from shutil import move
import sys, string
import numpy as np
from utilities import tag_file_name

SESSION_KEY = ''.join(np.random.choice(tuple(string.ascii_letters), 8))
TEMP_EXT = '_temp_%s' % (SESSION_KEY)
MIN_READ_LEN = 20
RE_SITES = {'MboI':'^GATC_', 'AluI':'AG^CT', 'BglII':'A^GATC_T', 'MseI': 'T^TA_A', 'HaeIII': 'GG^CC'}

def get_ligation_junction(re_site):
  """
  AG^CT gives AG + CT = AG + CT 
              TC   GA            
  
  ^GATC_ gives ---- + GATC = GATC + GATC
               CTAG   ----
  
  A^GATC_T gives A---- + GATCT = A + GATC + GATC + T
                 TCTAG   ----A
                  
  A_GATC^T gives AGATC + ----T = A + GATC + GATC + T
                 T----   CTAGA
  """
  
  if '_' in re_site:
    i, j = sorted([re_site.index('^'), re_site.index('_')])
    seq = re_site[:i] + re_site[i+1:j] + re_site[i+1:j] + re_site[j+1:]
    
  else: # Blunt
    seq = re_site.replace('^', '')
  
  return seq  

def clip_reads(fastq_file, file_root, junct_seq, replaced_seq, is_second=False, min_len=MIN_READ_LEN):
  """
  Clips reads at ligation junctions, removes anbigous ends and discards very short reads
  """
  tag = 'reads2_clipped' if is_second else 'reads1_clipped'
  clipped_file = tag_file_name(file_root, tag, '.fastq')
  clipped_file_temp = clipped_file + TEMP_EXT
    
  in_file_obj = open(fastq_file, 'r')
  n_rep = len(replaced_seq)
  n_junc = len(junct_seq)
  
  n_reads = 0
  n_clip = 0
  n_short = 0
  mean_len = 0
  
  out_file_obj = open(clipped_file_temp, 'w', 2**16)
  write = out_file_obj.write
  readline = in_file_obj.readline
  
  line1 = readline()
  while line1[0] != '@':
    line1 = readline()
  
  while line1:
    n_reads += 1
    line2 = readline()[:-1]
    line3 = readline()
    line4 = readline()[:-1]

    if junct_seq in line2:
      n_clip += 1
      i = line2.index(junct_seq)
      line2 = line2[:i] + replaced_seq
      line4 = line4[:i+n_rep]

    if not any(base in line2 for base in ('A', 'C', 'G', 'T')):
      line1 = readline()
      continue

    while line2[-1] == 'N':
      line2 = line2[:-1]
      line4 = line4[:-1]
    
    if n_junc < n_rep:
      while len(line4) < len(line2):
        line4 += line4[-1]
    
    n = len(line2)
    
    if n < min_len:
      n_short += 1
      
    mean_len += n
    write('%s%s\n%s%s\n' % (line1, line2, line3, line4))
 
    line1 = readline()
  
  if n_reads:
    mean_len /= n_reads

  # to do later: stat file
  stats = [('input_reads',n_reads),
           ('clipped',(n_clip, n_reads)),
           ('too_short',(n_short, n_reads)),
           ('mean_length',mean_len)]
  
  stat_key = 'clip_2' if is_second else 'clip_1'
  
  move(clipped_file_temp, clipped_file)
  return clipped_file

if __name__ == '__main__':
  re_site = RE_SITES[sys.argv[3]]
  lig_junc = get_ligation_junction(re_site)
  re_seq = re_site.replace('^', '').replace('_','')

  outfile = clip_reads(sys.argv[1], sys.argv[2], junct_seq=lig_junc, replaced_seq=re_seq, is_second=(sys.argv[4]=="True"))


