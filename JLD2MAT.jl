using MAT
using JLD

data_name = readline()
data = load(data_name)
file = matopen(data_name[1:end-3] * "mat", "w")
write(file, "data", data)
close(file)
