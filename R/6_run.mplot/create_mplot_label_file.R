list.files("data-raw/sam.files")
filenames <- list.files("data-raw/sam.files")
write.csv(filenames,file="data-raw/mplot_labels/Epac_labels_all.csv",row.names = FALSE, col.names = FALSE)
