def collect_data(files, folders, name):
	reads_arr = {}
	ratio_arr = {}
	gini_arr = {}
	cells = []

	failures = defaultdict(list)
	incomplete = defaultdict(list)
	count = 1
	for folder in folders:
		reads_arr[folder] = []
		ratio_arr[folder] = []
		gini_arr[folder] = []
		for file in files:
			if "cis_to_trans.out" in file:
				if count == 1:
					cells.append(file.split('.')[0])
				try:
					file_obj = open(folder+'/'+file)
				except IOError:
					incomplete[file].append(folder)
				try:
					reads, short_cis, long_cis, trans, ratio, percent_long, gini = file_obj.readlines()
				except ValueError:
					failures[file].append(folder)
				reads_arr[folder].append(reads.strip('\n').split('\t')[1].replace(",", ""))
				ratio_arr[folder].append(ratio.strip('\n').split('\t')[1])
				gini_arr[folder].append(gini.strip('\n').split('\t')[1])
				file_obj.close()
		count += 1
	print(incomplete)

def main():

if __name__ == '__main__':
	main()