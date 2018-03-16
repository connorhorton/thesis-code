import os, sys
from natsort import natsorted, ns
from scipy import stats
import numpy as np
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
from decimal import Decimal
from collections import defaultdict

# labels = {"1xreplicatesonly": ("Dedup\nAll reads", "Discard singletons,\nthen dedup"), \
# 	"discard1x": ("Discard\nSingletons","Keep\nSingletons"), \
# 	"keepallreads": ("Discard\nDuplicates", "Keep\nAll reads"), \
# 	"promiscuous": ("Keep\nPromiscuous reads", "Discard\nPromiscuous reads"), \
# 	"tooshort": ("Discard\nReads <20 bp", "Keep\nReads <20 bp"), \
# 	"toobig": ("Discard Reads\n>2000 bp", "Keep Reads\n>2000 bp"), \
# 	"undigested": ("Keep reads\nw/ restriction site", "Discard reads\nw/ restriction site"), \
# 	"UUonly": ("Keep\n\"Rescued reads\"", "Discard\n\"Rescued reads\"") }

labels = {"1xreplicatesonly": ("discard", "keep singletons"), \
	"discard1x": ("discard", "keep singletons"), \
	"keepallreads": ("keep", "discard duplicates"), \
	"promiscuous": ("discard", "keep promiscuous reads"), 
	"multi": ("keep", "discard multimapping reads"), \
	"tooshort": ("discard", "keep reads <20 bp"), \
	"toobig": ("discard", "keep reads >2000 bp"), \
	"undigested": ("discard", "keep undigested reads"), \
	"UUonly": ("discard", "keep rescued reads") }

def make_subplot(folder, axes, control, exp, xlabels):
	axes.boxplot([control, exp], positions=(1,2), showmeans=True)
	t_stat, p_val = stats.ttest_rel(control, exp)
	plt.xticks((1,2), xlabels, rotation=30, fontsize=13)
	plt.yticks(fontsize=12)
	ymin, ymax = axes.get_ylim()
	axes.plot([1, 1, 2, 2], [ymax, ymax*1.05, ymax*1.05, ymax], lw=1, c='k')
	if p_val < 1e-5:
		axes.text(1.5, ymax*1.05, "p<1e-5", ha='center', va='bottom', color='k')
	elif p_val < 1e-3:
		axes.text(1.5, ymax*1.05, "p=%.2e" % Decimal(p_val), ha='center', va='bottom', color='k')
	elif p_val < 0.05:
		axes.text(1.5, ymax*1.05, "p=%.3f" % p_val, ha='center', va='bottom', color='k')
	else:
		axes.text(1.5, ymax*1.05, "n.s.", ha='center', va='bottom', color='k')
	axes.set_ylim((ymin,ymax*1.12))
	return t_stat, p_val

def make_plot(name, figure_array, folder, gini_arr, reads_arr, ratio_arr, gini_control, reads_control, ratio_control):
	print("Making figures for %s!" % folder)
	figure_array[folder] = plt.figure(figsize=(10,3.65)) #figsize(width,height)

	xlabels = labels[folder]

	gini_exp = np.array(list(gini_arr[folder].values()),dtype=np.float64)
	gini_ax = figure_array[folder].add_subplot(1,3,1)
	gini_ax.set_ylabel("Gini coefficient", fontsize=13)
	gini_t, gini_p = make_subplot(folder, gini_ax, gini_control, gini_exp, xlabels)
	print("gini t-statistic: %f, p-val: %.2e for %s" % (gini_t, Decimal(gini_p), folder))

	reads_exp = np.array(list(reads_arr[folder].values()),dtype=np.float64)
	reads_ax = figure_array[folder].add_subplot(1,3,2)
	reads_ax.set_ylabel("Number of contacts", fontsize=13)
	reads_t, reads_p = make_subplot(folder, reads_ax, reads_control, reads_exp, xlabels)
	print("reads t-statistic: %f, p-val: %.2e for %s" % (reads_t, Decimal(reads_p), folder))

	ratio_exp = np.array(list(ratio_arr[folder].values()),dtype=np.float64)
	ratio_ax = figure_array[folder].add_subplot(1,3,3)
	ratio_ax.set_ylabel("Percentage of contacts in cis", fontsize=13)
	ratio_t, ratio_p = make_subplot(folder, ratio_ax, ratio_control, ratio_exp, xlabels)
	print("ratio t-statistic: %f, p-val: %.2e for %s" % (ratio_t, Decimal(ratio_p), folder))

	figure_array[folder].tight_layout()
	# figure_array[folder].subplots_adjust(top=0.85)
	# figure_array[folder].suptitle("Effect of %s on data from %s et al." % (folder, name), fontsize=16)
	figure_array[folder].savefig("%s_%s.png" % (name,folder))
	plt.close(figure_array[folder])

def collect_data(dirname, files, folders, name):
	reads_arr = {}
	ratio_arr = {}
	gini_arr = {}

	failures = defaultdict(list)
	incomplete = defaultdict(list)
	for folder in folders:
		reads_arr[folder] = {}
		ratio_arr[folder] = {}
		gini_arr[folder] = {}
		for file in files:
			if "cis_to_trans.out" in file:
				try:
					file_obj = open(dirname+folder+"/"+file)
				except IOError:
					incomplete[file].append(folder)
					continue
				try:
					reads, short_cis, long_cis, trans, ratio, percent_long, gini = file_obj.readlines()
					reads_arr[folder][file] = reads.strip('\n').split('\t')[1].replace(",", "")
					ratio_arr[folder][file] = ratio.strip('\n').split('\t')[1]
					gini_arr[folder][file] = gini.strip('\n').split('\t')[1]
				except ValueError:
					failures[file].append(folder)
				file_obj.close()

	# print(len(failures))
	for file in failures:
		# print("removing...", file, failures[file])
		for folder in folders:
			if file in reads_arr[folder]:
				del reads_arr[folder][file] 
			if file in ratio_arr[folder]:
				del ratio_arr[folder][file] 
			if file in gini_arr[folder]:
				del gini_arr[folder][file] 

	for file in incomplete:
		print("deleting incomplete outfiles")
		for folder in folders:
			if file in reads_arr[folder]:
				del reads_arr[folder][file]
			if file in ratio_arr[folder]:
				del ratio_arr[folder][file] 
			if file in gini_arr[folder]:
				del gini_arr[folder][file]

	return reads_arr, gini_arr, ratio_arr

def make_statfile(reads, gini, ratio, folder, outfile):
	# print gini index, then reads, then ratio for all cells
	outfile.write(",".join(gini[folder].values()))
	outfile.write("\n")
	outfile.write(",".join(reads[folder].values()))
	outfile.write("\n")	
	outfile.write(",".join(ratio[folder].values()))


def main():
	dirs = []
	dataset_names = {}
	for i in range(1,len(sys.argv)):
		if i > 0:
			if i % 2 == 1:
				dirs.append(sys.argv[i])
			else:
				dataset_names[dirs[-1]] = sys.argv[i]

	folders = defaultdict(list)
	files = {}
	reads, gini, ratio = {}, {}, {}
	gini_control, reads_control, ratio_control = {}, {}, {}
	xlabels = []
	for dir_ in dirs:
		for (dirpath, dirnames, filenames) in os.walk(dir_):
			folders[dir_].extend(dirnames)
			files[dir_] = filenames
		folders[dir_].sort()
		files[dir_] = natsorted(files[dir_], alg=ns.IGNORECASE)
		reads[dir_], gini[dir_], ratio[dir_] = collect_data(dir_, files[dir_], folders[dir_], dataset_names[dir_])
		gini_control[dir_] = np.array([gini[dir_]["control"][file] for file in gini[dir_]["control"]],dtype=np.float64)
		reads_control[dir_] = np.array([reads[dir_]["control"][file] for file in reads[dir_]["control"]],dtype=np.float64)
		ratio_control[dir_] = np.array([ratio[dir_]["control"][file] for file in ratio[dir_]["control"]],dtype=np.float64)
		
		# figure_array = {}
		# print("making figures for %s" % dataset_names[dir_])
		# for folder in folders[dir_]:
		# 	if folder != "control":
		# 		make_plot(dataset_names[dir_], figure_array, folder, \
		# 			gini[dir_], reads[dir_], ratio[dir_], \
		# 			gini_control[dir_], reads_control[dir_], ratio_control[dir_])

		statfile_name = dataset_names[dir_] + "_stats.txt"
		outfile = open(statfile_name, 'w')
		make_statfile(reads[dir_], gini[dir_], ratio[dir_], "control", outfile)
		outfile.close()

		if "_" in dataset_names[dir_]:
			split = dataset_names[dir_].split("_")
			label = split[0]+(" ${et}$ ${al.}$\n(%s cells)" % split[1])
			xlabels.append(label)
		else:
			xlabels.append(dataset_names[dir_]+" ${et}$ ${al.}$")

	common_folders = set(folders[dirs[0]])
	for folder in folders:
		common_folders = common_folders.intersection(set(folders[folder]))


	for folder in common_folders:
		gini_diff, reads_diff, ratio_diff = {}, {}, {}
		gini_t, gini_p, reads_t, reads_p, ratio_t, ratio_p = {}, {}, {}, {}, {}, {}
		if folder != "control":
			if folder in ["toobig", "tooshort"]:
				for dir_ in dirs:
					gini_exp = np.array([gini[dir_][folder][file] for file in gini[dir_][folder]],dtype=np.float64)
					gini_diff[dir_] = np.subtract(gini_control[dir_], gini_exp) # elementwise subtraction (similar to paired t-test)
					gini_t[dir_], gini_p[dir_] = stats.ttest_rel(gini_control[dir_], gini_exp) # paired t-test

					reads_exp = np.array([reads[dir_][folder][file] for file in reads[dir_][folder]],dtype=np.float64)
					reads_diff[dir_] = np.subtract(reads_control[dir_], reads_exp)
					reads_t[dir_], reads_p[dir_] = stats.ttest_rel(reads_control[dir_], reads_exp)

					ratio_exp = np.array([ratio[dir_][folder][file] for file in ratio[dir_][folder]],dtype=np.float64)
					ratio_diff[dir_] = np.subtract(ratio_control[dir_], ratio_exp)
					ratio_t[dir_], ratio_p[dir_] = stats.ttest_rel(ratio_control[dir_], ratio_exp)
			else:
				for dir_ in dirs:
					gini_exp = np.array([gini[dir_][folder][file] for file in gini[dir_][folder]],dtype=np.float64)
					gini_diff[dir_] = np.subtract(gini_exp, gini_control[dir_]) # elementwise subtraction (similar to paired t-test)
					gini_t[dir_], gini_p[dir_] = stats.ttest_rel(gini_exp, gini_control[dir_]) # paired t-test

					reads_exp = np.array([reads[dir_][folder][file] for file in reads[dir_][folder]],dtype=np.float64)
					reads_diff[dir_] = np.subtract(reads_exp, reads_control[dir_])
					reads_t[dir_], reads_p[dir_] = stats.ttest_rel(reads_exp, reads_control[dir_])

					ratio_exp = np.array([ratio[dir_][folder][file] for file in ratio[dir_][folder]],dtype=np.float64)
					ratio_diff[dir_] = np.subtract(ratio_exp, ratio_control[dir_])
					ratio_t[dir_], ratio_p[dir_] = stats.ttest_rel(ratio_exp, ratio_control[dir_])

			gini_fig = plt.figure(figsize=(8,5))
			plt.ylabel("Difference in GiniQC\n(%s - %s)" % labels[folder], fontsize=13)
			plt.boxplot([gini_diff[dir_] for dir_ in dirs])
			for dir_ in dirs:
				print("gini stats %s and %s" % (dataset_names[dir_], folder))
				print(stats.describe(gini_diff[dir_]))
			plt.plot(plt.xlim(), (0,0), '--')
			plt.xticks((1,2,3,4), xlabels, fontsize=13)
			plt.yticks(fontsize=12)
			ymin, ymax = plt.ylim()
			plt.ylim((ymin,ymax*1.1))
			count = 1
			for dir_ in dirs:
				p_val = gini_p[dir_]
				if p_val < 1e-5:
					plt.text(count, ymax*1.03, "p<1e-5", ha='center', va='top', color='k', fontsize=12)
				elif p_val < 1e-3:
					plt.text(count, ymax*1.03, "p=%.2e" % Decimal(p_val), ha='center', va='top', color='k', fontsize=12)
				elif p_val < 0.05:
					plt.text(count, ymax*1.03, "p=%.3f" % p_val, ha='center', va='top', color='k', fontsize=12)
				else:
					plt.text(count, ymax*1.03, "n.s.", ha='center', va='top', color='k', fontsize=12)
				count += 1
			gini_fig.tight_layout()
			plt.savefig("ginifig_%s.png" % folder)
			plt.close(gini_fig)

			reads_fig = plt.figure(figsize=(8,5))
			plt.ylabel("Difference in total number of contacts\n(%s - %s)" % labels[folder], fontsize=13)
			# plt.yscale('symlog')
			plt.boxplot([reads_diff[dir_] for dir_ in dirs])
			for dir_ in dirs:
				print("reads stats %s and %s" % (dataset_names[dir_], folder))
				print(stats.describe(reads_diff[dir_]))
			plt.plot(plt.xlim(), (0,0), '--')
			plt.xticks((1,2,3,4), xlabels, fontsize=13)
			plt.yticks(fontsize=12)
			ymin, ymax = plt.ylim()
			plt.ylim((ymin,ymax*1.1))
			count = 1
			for dir_ in dirs:
				p_val = reads_p[dir_]
				if p_val < 1e-5:
					plt.text(count, ymax*1.03, "p<1e-5", ha='center', va='top', color='k', fontsize=12)
				elif p_val < 1e-3:
					plt.text(count, ymax*1.03, "p=%.2e" % Decimal(p_val), ha='center', va='top', color='k', fontsize=12)
				elif p_val < 0.05:
					plt.text(count, ymax*1.03, "p=%.3f" % p_val, ha='center', va='top', color='k', fontsize=12)
				else:
					plt.text(count, ymax*1.03, "n.s.", ha='center', va='top', color='k', fontsize=12)
				count += 1
			reads_fig.tight_layout()
			plt.savefig("readsfig_%s.png" % folder)
			plt.close(reads_fig)

			ratio_fig = plt.figure(figsize=(8,5))
			plt.ylabel("Difference in percentage of contacts in ${cis}$\n(%s - %s)" % labels[folder], fontsize=13)
			plt.boxplot([ratio_diff[dir_] for dir_ in dirs])
			for dir_ in dirs:
				print("cis/trans stats for %s and %s" % (dataset_names[dir_], folder))
				print(stats.describe(ratio_diff[dir_]))
			plt.plot(plt.xlim(), (0,0), '--')
			plt.xticks((1,2,3,4), xlabels, fontsize=13)
			plt.yticks(fontsize=12)
			ymin, ymax = plt.ylim()
			plt.ylim((ymin,ymax*1.1))
			count = 1
			for dir_ in dirs:
				p_val = ratio_p[dir_]
				if p_val < 1e-5:
					plt.text(count, ymax*1.03, "p<1e-5", ha='center', va='top', color='k', fontsize=12)
				elif p_val < 1e-3:
					plt.text(count, ymax*1.03, "p=%.2e" % Decimal(p_val), ha='center', va='top', color='k', fontsize=12)
				elif p_val < 0.05:
					plt.text(count, ymax*1.03, "p=%.3f" % p_val, ha='center', va='top', color='k', fontsize=12)
				else:
					plt.text(count, ymax*1.03, "n.s.", ha='center', va='top', color='k', fontsize=12)
				count += 1
			ratio_fig.tight_layout()
			plt.savefig("ratiofig_%s.png" % folder)
			plt.close(ratio_fig)

			print("\n")


if __name__ == '__main__':
	main()