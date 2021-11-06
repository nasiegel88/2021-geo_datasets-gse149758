source("R/packages.r")

if (length(ls(pattern = "_h5")) == 0) {source("R/download-data.r")}

if (length(list.files(path = "./output/rdata")) == 0) {source("R/combine-objects.r")} 
