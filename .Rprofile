# Set library path
.libPaths(c("/srv/conda/envs/notebook/lib/R/library" , .libPaths() ) )

# Install TinyTex
if (!tinytex::is_tinytex()) tinytex::install_tinytex()

dir.create('bin', showWarnings = FALSE)

# Causes issues if loaded before xfun
#library('utils')
#remove.packages("xfun")