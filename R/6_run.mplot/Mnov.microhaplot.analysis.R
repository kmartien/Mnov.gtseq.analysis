library(vcfR)
library(microhaplot)

run.label <- "RunMS58"

# for your dataset: customize the following paths
# sam.path <- paste0("data-raw/sam.files/", run.label, ".all")
# label.path <- "data-raw/mplot_labels/RunMS51.label.txt"
 vcf.path <- "vcf/RunMS58.filtered.recode.vcf"
# app.path <- "~/Shiny/microhaplot"

sam.path <- paste0("data-raw/sam.files/", run.label)
label.path <- file.path("data-raw/mplot_labels", paste0(run.label, ".label.txt"))
#vcf.path <- paste0("vcf/", run.label, ".vcf")
out.path <- "results-R/microhaplot"
app.path <- "/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot"

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

runShinyHaplot(app.path)
