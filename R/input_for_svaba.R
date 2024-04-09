library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--breakpoints_bed"), type="character", default=NULL, help="output_of_breakpoints_to_bed.py"),
  make_option(c("--aliquot_barcode"), type="character", default=NULL, help="aliquot_barcode"),
  make_option(c("--output_path"), type="character", default=NULL, help="output_path")
)

opt = parse_args(OptionParser(option_list=option_list))

breakpoints <- read_tsv(opt$breakpoints_bed, col_names = F)

breakpoints_1kb <- rbind(breakpoints %>% 
        mutate(chr = gsub('chr', '', X1),
               start = X2 %>% as.numeric() - 500,
               end = X2 %>% as.numeric() + 500) %>% 
        select(chr, start, end),
      breakpoints %>% 
        mutate(chr = gsub('chr', '', X3),
               start = X4 %>% as.numeric() - 500,
               end = X4 %>% as.numeric() + 500) %>% 
        select(chr, start, end))

write_tsv(breakpoints_1kb, 
          paste0(opt$output_path, opt$aliquot_barcode, "_aa_1kb_bp_svaba_input.bed"),
          col_names = FALSE)
