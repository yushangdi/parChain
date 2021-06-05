import sys

import scipy.cluster.hierarchy as sch
import numpy as np
import pandas as pd
import time


ADDR = "./datasets/"

def main(dataset, method, metric='euclidean'):
	data_addr = ADDR + dataset + ".pbbs"
	data = pd.read_csv(data_addr,  sep = " ", header=None)
	# print(data.head())
	start_time = time.time()
	Z = sch.linkage(data, method, metric=metric)
	end_time = time.time()
	print(dataset, method + "-" + metric, np.sum(Z[:,2]), end_time-start_time)
	


if __name__ == "__main__":
	dataset = sys.argv[1]
	method = sys.argv[2]
	if method not in ["avg", "avgsq", 'complete', "average", "single", "ward", "all"]:
		print("invalid method ", method)
	else:
		if method == "all":
			main(dataset, 'complete')
			main(dataset, "ward")
			main(dataset, "average")
			main(dataset, "average", "sqeuclidean")
		elif method == "avgsq":
			main(dataset, "average", "sqeuclidean")
		elif method == "avg":
			main(dataset, "average")
		else:
			main(dataset, method)
                        
