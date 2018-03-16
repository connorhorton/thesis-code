import sys
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
import cooler
import numpy as np
from gini import normalize_matrix

def main():
	matrix = cooler.Cooler(sys.argv[1])
	normalized = normalize_matrix(matrix)
	normalized.sort()

	cumulative = np.array(normalized,dtype=np.float64)
	for i in range(1,len(cumulative)):
		cumulative[i] += cumulative[i-1]

	fig = plt.figure(figsize=(5.5,4.5))
	lorenz = plt.plot(cumulative, 'b')
	equal = plt.plot(range(0,len(cumulative)), np.linspace(min(cumulative),max(cumulative),num=len(cumulative)), 'r')

	plt.ylabel("Cumulative number of contacts", fontsize=16)
	plt.xlabel("Cumulative number of genomic bins\n (sorted by number of contacts)", fontsize=16)
	plt.legend((lorenz[0], equal[0]), ("Lorenz Curve", "Line of Equality"), fontsize=13)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	name = sys.argv[1].split(".")[0]

	fig.tight_layout()
	fig.savefig("cumulativetrans_%s.png" % name, dpi=500)


if __name__ == '__main__':
	main()