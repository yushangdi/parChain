import sys

from sklearn.cluster import AgglomerativeClustering
import scipy
import numpy as np
import pandas as pd
import time


ADDR = "~/datasets/"

def main(dataset, method):
	data_addr = ADDR + dataset + ".pbbs"
	data = pd.read_csv(data_addr,  sep = " ", header=None)
	start_time = time.time()
	model = AgglomerativeClustering(distance_threshold=0, n_clusters=None, linkage=method, affinity = "euclidean", compute_full_tree = True)
	model = model.fit(data)
	end_time = time.time()
	print(dataset, method + "-euclidean", np.sum(model.distances_), end_time-start_time)

def main_sqeuc(dataset):
	data_addr = ADDR + dataset + ".pbbs"
	data = pd.read_csv(data_addr,  sep = " ", header=None)
	start_time = time.time()
	matrix = scipy.spatial.distance.cdist(data,data, metric = "sqeuclidean")
	print("compute matrix", time.time() - start_time)
	model = AgglomerativeClustering(distance_threshold=0, n_clusters=None, linkage="average", affinity = "precomputed", compute_full_tree = True)
	model = model.fit(matrix)
	end_time = time.time()
	print(dataset, "average-sqeuclidean", np.sum(model.distances_), end_time-start_time)
	


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
			if dataset != "10D_UCI4_100K": # oom
				main_sqeuc(dataset)
			else:
				print(dataset, "will oom")
		elif method == "avgsq":
			if dataset != "10D_UCI4_100K": # oom
				main_sqeuc(dataset)
			else:
				print(dataset, "will oom")
		elif method == "avg":
			main(dataset, "average")
		else:
			main(dataset, method)
