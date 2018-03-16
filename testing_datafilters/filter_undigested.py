import sys

def filter_undigested(infile, accepted, discarded):
	pairs = open(infile, 'r')
	accepted_file = open(accepted, 'w')
	discarded_file = open(discarded, 'w')

	for line in pairs:
		if "#" in line:
			accepted_file.write(line)
			discarded_file.write(line)
		else:
			readID, chr_a, pos_a, chr_b, pos_b, strand_a, strand_b, pair_type, re1_a_id, re1_a_start, re1_a_end, re1_b_id, re1_b_start, re1_b_end = line.split()
			if chr_a == chr_b and re1_a_id == re1_b_id: # Same Re1 fragment
				discarded_file.write(line)
			else:
				accepted_file.write(line)
			

	pairs.close()
	accepted_file.close()
	discarded_file.close()

if __name__ == '__main__':
	filter_undigested(sys.argv[1], sys.argv[2], sys.argv[3]) # arg 1: input pairs file; arg 2: accepted pairs file; arg 3: discarded pairs file