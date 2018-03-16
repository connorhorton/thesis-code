from shutil import move
import sys, string
import numpy as np
from utilities import tag_file_name
from collections import defaultdict

SESSION_KEY = ''.join(np.random.choice(tuple(string.ascii_letters), 8))
TEMP_EXT = '_temp_%s' % (SESSION_KEY)

def remove_promiscuous(pairs_file, resolve_limit=1e3):

  # resolve_limit : Allow two suitably close promiscous ends if long range cis or trans
 
  clean_pairs_file = tag_file_name(pairs_file, 'clean')
  clean_pairs_file_temp = clean_pairs_file + TEMP_EXT

  promiscous_pairs_file = tag_file_name(pairs_file, 'promisc')  
  promisc_file_obj = open(promiscous_pairs_file, 'w')
  write_promisc = promisc_file_obj.write 
  
  in_file_obj = open(pairs_file, 'r')
  clean_file_obj = open(clean_pairs_file_temp, 'w')
  write_clean = clean_file_obj.write
  
  frag_counts = defaultdict(set)
  n_promiscuous = 0
  n_resolved = 0
  n_clean = 0
  
  for line in in_file_obj:
    line_data = line.split()
    if '#' not in line_data[0]:
      chr_a = line_data[1]
      chr_b = line_data[3]
      re_start_a = int(line_data[9])
      re_start_b = int(line_data[13])
      frag_counts[(chr_a, re_start_a)].add((chr_b, re_start_b))
      frag_counts[(chr_b, re_start_b)].add((chr_a, re_start_a))

  remove = set()
  for re_start, re_ends in frag_counts.items():
    if len(re_ends) > 1: # we only care if the same fragment is ligating two diff places 
        
      if len(re_ends) == 2: 
        chr_a, re_start_a = re_start
        chr_b1, re_start_b1 = re_ends.pop()
        chr_b2, re_start_b2 = re_ends.pop()
        
        if (chr_b1 == chr_b2) and abs(re_start_b1-re_start_b2) < resolve_limit: # Other ends are very close to each other
            n_resolved += 1
        
        else: # Other ends too far apart
          remove.add(re_start)
       
      else:
        remove.add(re_start)    
  
  n_trans = 0
  n_cis_near = 0
  n_cis_far = 0
      
  in_file_obj.seek(0)
  for line in in_file_obj:
    line_data = line.split()
    if '#' in line_data[0]:
      write_clean(line)
      write_promisc(line)

    else:
      chr_a = line_data[1]
      chr_b = line_data[3]
      re_start_a = int(line_data[9])
      re_start_b = int(line_data[13])

      if (chr_a, re_start_a) in remove or (chr_b, re_start_b) in remove:
        n_promiscuous += 1    
        write_promisc(line)
      else:
        n_clean += 1
        write_clean(line)      

  in_file_obj.close()
  promisc_file_obj.close()
  
  n = n_promiscuous + n_clean
  n_clean -= n_resolved
  
  # stats = [('input_pairs', n),
  #          ('clean',(n_clean, n)),
  #          ('promiscuous',(n_promiscuous, n)),
  #          ('resolved',(n_resolved, n)),
  #          ('accepted',(n_clean+n_resolved, n)),
  #          ]

  # stat_key = 'promsic'

  move(clean_pairs_file_temp, clean_pairs_file)
  
  return clean_pairs_file

if __name__ == '__main__':
  remove_promiscuous(sys.argv[1])