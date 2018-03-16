import sys
NCC_FORMAT = '%s %d %d %d %d %s %s %d %d %d %d %s %d %d %d\n'

def main():
	pairs = open(sys.argv[1], 'r')
	sam = open(sys.argv[2], 'r')
	ncc = open(sys.argv[3], 'w')

	readID = []
	chr_a = {} 
	start_a = {}
	end_a = {}
	re1_a_start = {}
	re1_a_end = {}
	strand_a = {}
	chr_b = {}
	start_b = {}
	end_b = {}
	re1_b_start = {}
	re1_b_end = {}
	strand_b = {}

	for line in pairs:
		if "#" not in line:
			split = line.split('\t')
			ID = split[0].split('.')[1]
			readID.append(ID)
			chr_a[ID] = split[1]
			start_a[ID] = split[2]
			chr_b[ID] = split[3]
			start_b[ID] = split[4]
			strand_a[ID] = split[5]
			strand_b[ID] = split[6]
			re1_a_start[ID] = split[9]
			re1_a_end[ID] = split[10]
			re1_b_start[ID] = split[12]
			re1_b_end[ID] = split[13]

	for line in sam:
		split = line.split('\t')
		if '@' not in split[0]:
			ID = split[0].split('.')[1]
			seq = split[9]
			if strand_a == "+":
				end_a[ID] = int(start_a[ID]) + len(seq)
			else:
				end_a[ID] = start_a[ID]
				start_a[ID] = int(start_a[ID]) - len(seq)
			if strand_b == "+":
				end_b[ID] = int(start_b[ID]) + len(seq)
			else:
				end_b[ID] = start_b[ID]
				start_b[ID] = int(start_b[ID]) - len(seq)

	for ID in readID:
		ncc.write(NCC_FORMAT % (chr_a[ID], int(start_a[ID]), int(end_a[ID]), int(re1_a_start[ID]), int(re1_a_end[ID]), strand_a[ID], \
			chr_b[ID], int(start_b[ID]), int(end_b[ID]), int(re1_b_start[ID]), int(re1_b_end[ID]), strand_b[ID], int(ID), int(ID), 0))

if __name__ == '__main__':
	main()