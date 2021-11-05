## Download data from GEO database

temp.control_LD09 <- paste(tempfile(),".h5",sep = "")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4511886&format=file&file=GSM4511886%5Fcontrol%5FLD09%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5", 
              temp.control_LD09, mode = "wb")
control_31810_h5 <- Read10X_h5(temp.control_LD09)

temp.control_31810 <- paste(tempfile(),".h5",sep = "")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4511885&format=file&file=GSM4511885%5Fcontrol%5F31810%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5", 
              temp.control_31810, mode = "wb")
control_31810_h5 <- Read10X_h5(temp.control_31810)

temp.control_LE99 <- paste(tempfile(),".h5",sep = "")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM4511887&format=file&file=GSM4511887%5Fcontrol%5FLE99%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5", 
              temp.control_LE99, mode = "wb")
control_LE99_h5 <- Read10X_h5(temp.control_LE99)
