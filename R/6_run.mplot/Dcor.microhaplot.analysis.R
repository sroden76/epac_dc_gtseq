library(vcfR)
library(microhaplot)

run.label <- "all_Epac"

sam.path <- paste0("data-raw\\sam.files")
label.path <- file.path("data-raw\\mplot_labels\\EPac_labels_all.txt")
vcf.path <- "vcf\\Dcor_DcPanel_205_maf.targetSNPs_012224.recode.vcf"
#out.path <- "results-R\\microhaplot"
app.path <- "C:\\Users\\suzanne.roden\\Documents\\Shiny\\microhaplot"

vcf <- read.vcfR(vcf.path)


# I've prepped the data, so can just jump straight to running the Shiny app
haplo.read.tbl <- prepHaplotFiles(run.label = run.label,
                                  sam.path = sam.path,
                                  #out.path = out.path,
                                  label.path = label.path,
                                  vcf.path = vcf.path,
                                  app.path = app.path,
                                  n.jobs = 4) # use all the cores!
#})
runShinyHaplot(app.path)
