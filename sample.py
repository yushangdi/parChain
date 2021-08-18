import pandas as pd
ADDR = "./datasets/"
# data_addr = ADDR + "2D_UniformFill_10M" + ".pbbs"
# data = pd.read_csv(data_addr, sep = " ", header=None).drop(2, axis = 1) #skiprows =1, 
# print(data.head())
# print()
# data_sample = data.sample(int(1e5))
# print(data_sample .head())
# data_sample.to_csv( ADDR + "2D_UniformFill_100K" + ".pbbs", sep = " ", header = None, index = False)

# data_sample = data.sample(int(1e4))
# print(data_sample .head())
# data_sample.to_csv( ADDR + "2D_UniformFill_10K" + ".pbbs", sep = " ", header = None, index = False)


data_addr = ADDR + "3D_GeoLife_24M" + ".pbbs"
data = pd.read_csv(data_addr, sep = " ", header=None).drop(3, axis = 1) #skiprows =1, 
print(data.head())
print()
# data_sample = data.sample(int(1e5))
# print(data_sample .head())
# data_sample.to_csv( ADDR + "3D_GeoLife_100K" + ".pbbs", sep = " ", header = None, index = False)

# data_sample = data.sample(int(1e4))
# print(data_sample .head())
# data_sample.to_csv( ADDR + "3D_GeoLife_10K" + ".pbbs", sep = " ", header = None, index = False)

data_sample = data.sample(int(1e6))
print(data_sample .head())
data_sample.to_csv( ADDR + "3D_GeoLife_1M" + ".pbbs", sep = " ", header = None, index = False)

data_sample = data.sample(int(1e7))
print(data_sample .head())
data_sample.to_csv( ADDR + "3D_GeoLife_10M" + ".pbbs", sep = " ", header = None, index = False)

# data_sample = data.sample(int(1e3))
# print(data_sample .head())
# data_sample.to_csv( ADDR + "3D_GeoLife_1K" + ".pbbs", sep = " ", header = None, index = False)