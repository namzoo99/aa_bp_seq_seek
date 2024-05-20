library(tidyverse)
library(optparse)

option_list <- list(
  make_option(c("--aa_output_path"), type="character", default=NULL, help="aa_output_path"),
  make_option(c("--aliquot_barcode"), type="character", default=NULL, help="aliquot_barcode"),
  make_option(c("--output_path"), type="character", default=NULL, help="output_path")
)

opt = parse_args(OptionParser(option_list=option_list))


# this is for optoin '-i' in 'breakpoints_to_bed.py'
cycle_path <- list.files(opt$aa_output_path,
           pattern = '_cycles.txt$', recursive = TRUE, full.names = TRUE) %>% 
  .[str_detect(., "_AA_results")]

graph_path <- list.files(opt$aa_output_path,
           pattern = '_graph.txt$', recursive = TRUE, full.names = TRUE) %>% 
  .[str_detect(., "_AA_results")]

write_tsv(data.frame(rep(opt$aliquot_barcode, length(cycle_path)),
                     cycle_path,
                     graph_path),
          paste0(opt$output_path, opt$aliquot_barcode, "_bp.input.txt"),
          col_names = FALSE)


# this is for option '-r' in 'breakpoints_to_bed.py'
aa_summary_file <- list.files(opt$aa_output_path,
                            pattern = '_summary.txt$', recursive = TRUE, full.names = TRUE) %>% 
  .[str_detect(., "_AA_results")] %>% 
  lapply(., function(x){read_tsv(x, col_names=F)})


reshape_sumFile <- function(x){
  
  x %>%
    .[!str_detect(.$X1, "----"),] %>% .[-1,] %>%  # remove rows that --- row and indicates number of amplicons per sample
    separate(col=X1, into=c('amplicon_number', 'info', "sep", 'value'), sep='\\s+') %>%  # separate rows by 'space(\\s+)'
    .[c(1,2,4)] %>% spread(data=., key="info", value="value")
  
}

reshpaed_aa_summary <- lapply(aa_summary_file, reshape_sumFile)

aa_summary_intervals <- sapply(reshpaed_aa_summary, function(x){ separate_longer_delim(x %>% select(Intervals), Intervals, delim = ",") }) %>% 
  sapply(function(x){paste0("chr", x)})

write_tsv(aa_summary_intervals %>% t %>% data.frame,
          paste0(opt$output_path, opt$aliquot_barcode, "_aa_summary_intervals.txt"),
          col_names=FALSE)
