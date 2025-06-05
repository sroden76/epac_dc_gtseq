library(vcfR)
library(microhaplot)

run.label <- "all_Epac"

sam.path <- paste0("data-raw/sam.files")
label.path <- file.path("data-raw/mplot_labels/EPac_labels.txt")
vcf.path <- "vcf/Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"
out.path <- "results-R/microhaplot"
app.path <- "C:/Users/suzanne.roden/Documents/R/Microhaplot/Dc_mex_fem/microhaplot/shiny/microhaplot"

vcf <- read.vcfR(vcf.path)
#locus.ignore.file <- read.csv(paste0("microhaplot/",run.label, ".locus_annotation.csv"))

# I've prepped the data, so can just jump straight to running the Shiny app
haplo.read.tbl <- prepHaplotFiles(run.label = run.label,
                                  sam.path = sam.path,
                                  out.path = out.path,
                                  label.path = label.path,
                                  vcf.path = vcf.path,
                                  app.path = app.path,
                                  n.jobs = 4) # use all the cores!
#})
runShinyHaplot(app.path)
