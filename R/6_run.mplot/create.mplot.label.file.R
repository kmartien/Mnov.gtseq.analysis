library(tidyverse)

sam.folder <- 'final.sams'

labels <- tibble(fname = list.files(path = paste0("data-raw/sam.files/", sam.folder), pattern = ".sam")) |> 
  mutate(LABID = substr(fname, start = 1, stop = 8),
         stratum = 'NA')

write.table(labels, file = paste0("data-raw/mplot_labels/", sam.folder, ".label.txt"), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

