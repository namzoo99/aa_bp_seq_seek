library(tidyverse)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("--breakpoints_bed"), type="character", default=NULL, help="output_of_breakpoints_to_bed.py"),
  make_option(c("--aa_output_path"), type="character", default=NULL, help="aa_output_path"),  
  make_option(c("--svaba_sv_vcf"), type="character", default=NULL, help="filtered_svaba_sv_vcf_file"), 
  make_option(c("--svaba_unfiltered_sv_vcf"), type="character", default=NULL, help="unfiltered_svaba_sv_vcf_file"), 
  make_option(c("--aliquot_barcode"), type="character", default=NULL, help="aliquot_barcode"),
  make_option(c("--output_path"), type="character", default=NULL, help="output_path")
)

opt = parse_args(OptionParser(option_list=option_list))

################################################################################
#
# Define
#       input
# 
################################################################################
# output of 'process test3'
bp.input_breakpoints <- read_tsv(opt$breakpoints_bed,
                                 col_names = FALSE)


# directory 'cycle / annotated cycle' file
all_cycle_file_list <- list.files(opt$aa_output_path,
                                  pattern = "_cycles.txt", recursive = TRUE, full.names = TRUE) %>%
  .[str_detect(., paste(unique(bp.input_breakpoints$X6), collapse = "|"))] %>% 
  .[!str_detect(., "/files/")]


# svaba output
svaba.sv.vcf <- read_tsv(opt$svaba_sv_vcf, 
                         skip = 144)

svaba.unfiltered.sv.vcf <- read_tsv(opt$svaba_unfiltered_sv_vcf, 
                                    skip = 144)


################################################################################
#
# Add location info.('Start', 'End') to each segments
# 'Start' and 'End' columns reflect direction info. of segments(- or +)
# 
################################################################################
# Split 'annotated_cycle' file and 'cycle' file
annotated_cycle <- all_cycle_file_list[str_detect(all_cycle_file_list, "_annotated_cycles.txt")] %>%
  map_df(~read_tsv(., id = "file_name", col_names=FALSE))

cycle <- all_cycle_file_list[!str_detect(all_cycle_file_list, "_annotated_cycles.txt")] %>%
  map_df(~read_tsv(., id = "file_name", col_names=FALSE, col_types = "cccccc"))


# edit annotated cycle file 
edited_annotated_cycle <- annotated_cycle %>% .[str_detect(.$X1, "Cycle="),] %>% 
  separate(., col=X1,
           into = c("Cycle", "Copy_count", "Length", "IsCyclicPath", "CycleClass", "Segments"),
           sep = ";") %>% 
  #filter(CycleClass %in% "CycleClass=ecDNA-like") %>% 
  mutate(Segments = sapply(str_split(.$Segments, "="), '[[', 2)) %>% 
  separate_rows(., Segments, sep= ",")

edited_annotated_cycle$abs_Segments <- substr(edited_annotated_cycle$Segments, 1, nchar(edited_annotated_cycle$Segments)-1)

annotated_cycle_SegInfo <- annotated_cycle %>% .[str_detect(.$X1, "Segment\t"),] %>% 
  separate(., col=X1,
           into = c("Segments", "abs_Segments", "Chr", "Start", "End"),
           sep = "\t") 

annotated_cycle_bp <- left_join(edited_annotated_cycle, annotated_cycle_SegInfo %>% select(-Segments), by=c("abs_Segments", "file_name"))


# edit cycle file 
edited_cycle <- cycle %>% .[str_detect(.$X1, "Cycle="),] %>% 
  separate(., col=X1,
           into = c("Cycle", "Copy_count", "Segments"),
           sep = ";") %>% 
  .[1:4]%>% 
  #filter(CycleClass %in% "CycleClass=ecDNA-like") %>% 
  mutate(Segments = sapply(str_split(.$Segments, "="), '[[', 2)) %>% 
  separate_rows(., Segments, sep= ",")

edited_cycle$abs_Segments <- substr(edited_cycle$Segments, 1, nchar(edited_cycle$Segments)-1)

cycle_SegInfo <- cycle %>% .[str_detect(.$X1, "Segment"),] %>%
  mutate(Segments=X1,
         abs_Segments=X2 %>% as.character(),
         Chr=X3,
         Start=X4,
         End=X5) %>% 
  select(-X1, -X2, -X3, -X4, -X5)

cycle_bp <- left_join(edited_cycle, cycle_SegInfo %>% select(-Segments), by=c("abs_Segments", "file_name"))



################################################################################
#
# Create columns(bp_front, bp_back) about bp facing each other in cycle
# then +1 to these columns (svaba -> 1-based coord, aa -> 0-based coord)
#
################################################################################
matching_aa_svaba <- function(x, y){
  
  aa_bp <- x
  svaba_output <- y
  
  aa_bp$bp_front <- ifelse(str_detect(aa_bp$Segments, "-$") == TRUE, aa_bp$End, aa_bp$Start)
  aa_bp$bp_back <- ifelse(str_detect(aa_bp$Segments, "-$") == TRUE, aa_bp$Start, aa_bp$End)
  
  aa_bp$bp_front <- as.numeric(aa_bp$bp_front)
  aa_bp$bp_back <- as.numeric(aa_bp$bp_back)
  
  
  facing_aabp <- lapply(split(aa_bp, list(aa_bp$Cycle, aa_bp$file_name)) %>% .[which(lapply(., nrow) != 0)], function(x){
    
    if(nrow(x)==1){
      
      x2 <- x
      x2[nrow(x),"Chr_bp_front"] <- x[1,"Chr"]
      print(x2)
      
    }else{
      
      x2 <- x
      
      x2[nrow(x),"bp_front"] <- x[1,"bp_front"]
      x2[nrow(x),"Chr_bp_front"] <- x[1,"Chr"]
      
      x2[1:(nrow(x)-1), "bp_front"] <- x[2:nrow(x),"bp_front"]
      x2[1:(nrow(x)-1),"Chr_bp_front"] <- x[2:nrow(x),"Chr"]
      
      print(x2)
      
    } }) %>% do.call(rbind, .)
  
  rownames(facing_aabp) <- NULL
  
  facing_aabp <- facing_aabp %>% mutate(bp_front_1coords = bp_front+1,
                                        bp_back_1coords = bp_back+1,
                                        bp = paste0(Chr_bp_front, ":", (bp_front+1)),
                                        bp_mate = paste0(Chr, ":", (bp_back+1)))
  
  # Extract contig info. from 'svaba.sv.vcf' 
  contig.svaba <- sapply(str_split(svaba_output$INFO, ";"), function(x){ x[grep("SCTG", x)] }) %>% 
    cbind(svaba_output %>% select(`#CHROM`, POS, ALT), SCTG = .)
  
  contig.svaba <- contig.svaba %>% 
    mutate(Svaba_bp_mate = gsub("\\[|\\]", "", contig.svaba$ALT) %>% 
             gsub('[ATCG]', "", .))
  contig.svaba$`#CHROM` <- as.character(contig.svaba$`#CHROM`)
  
  contig.svaba <- contig.svaba %>% 
    mutate(Svaba_bp_mate_chrom = sapply(str_split(.$Svaba_bp_mate, ":"), '[', 1),
           Svaba_bp_mate_pos = sapply(str_split(.$Svaba_bp_mate, ":"), '[', 2)) %>% 
    relocate("#CHROM", POS, Svaba_bp_mate_chrom, Svaba_bp_mate_pos)
  
  contig.svaba$Svaba_bp_mate_pos <- as.numeric(contig.svaba$Svaba_bp_mate_pos)
  
  # Check for overlapping SVs in aa bp extended by 1, 10, 100, 500bp
  match_sctg <- function(facing_row, contig_data, bp_offset) {
    matched <- contig_data %>%
      filter(`#CHROM` %in% facing_row$Chr_bp_front,
             between(POS, facing_row$bp_front_1coords - bp_offset, facing_row$bp_front_1coords + bp_offset),
             Svaba_bp_mate_chrom %in% facing_row$Chr,
             between(Svaba_bp_mate_pos, facing_row$bp_back_1coords - bp_offset, facing_row$bp_back_1coords + bp_offset)) %>%
      select(SCTG) %>% ifelse(nrow(.) == 0, data.frame(SCTG = NA), .) %>%
      data.frame() %>% 
      apply(., 1, function(x){gsub("SCTG=", "", x)}) %>% 
      paste0(., collapse = ",")
    
    return(matched)
  }
  
  Matched_AA_Svaba <- data.frame()
  for (i in 1:nrow(facing_aabp)) {
    mathedDF_1bp <-  bind_cols(facing_aabp[i,], match_sctg(facing_aabp[i,], contig.svaba, 1))
    mathedDF_10bp <- bind_cols(facing_aabp[i,], match_sctg(facing_aabp[i,], contig.svaba, 10))
    mathedDF_100bp <- bind_cols(facing_aabp[i,], match_sctg(facing_aabp[i,], contig.svaba, 100))
    mathedDF_500bp <- bind_cols(facing_aabp[i,], match_sctg(facing_aabp[i,], contig.svaba, 500))
    
    merged_DF <- cbind(mathedDF_1bp,
                       mathedDF_10bp[grep("\\...", colnames(mathedDF_10bp))],
                       mathedDF_100bp[grep("\\...", colnames(mathedDF_100bp))],
                       mathedDF_500bp[grep("\\...", colnames(mathedDF_500bp))])
    
    Matched_AA_Svaba <- rbind(Matched_AA_Svaba, merged_DF)
  }
  
  colnames(Matched_AA_Svaba)[grep("\\...", colnames(Matched_AA_Svaba))] <- c("SCTG_1bp", "SCTG_10bp", "SCTG_100bp", "SCTG_500bp")
  
  return(Matched_AA_Svaba)
  
}

Matched_annotated_cycle_Svaba <- matching_aa_svaba(annotated_cycle_bp, svaba.sv.vcf)
Matched_cycle_Svaba <- matching_aa_svaba(cycle_bp, svaba.sv.vcf)

Matched_annotated_cycle_Svaba_unfiltered <-  matching_aa_svaba(cycle_bp, svaba.unfiltered.sv.vcf)
Matched_cycle_Svaba_unfiltered <- matching_aa_svaba(cycle_bp, svaba.unfiltered.sv.vcf)


################################################################################
#
# Save output file
# 
# 
################################################################################
write_tsv(Matched_annotated_cycle_Svaba, paste0(opt$output_path, opt$aliquot_barcode, "_matched_svaba_annotated_aa.tsv"), na = "")
write_tsv(Matched_cycle_Svaba, paste0(opt$output_path, opt$aliquot_barcode, "_matched_svaba_aa.tsv"), na = "")

write_tsv(Matched_annotated_cycle_Svaba_unfiltered, paste0(opt$output_path, opt$aliquot_barcode, "_matched_unfiltered_svaba_annotated_aa.tsv"), na = "")
write_tsv(Matched_cycle_Svaba_unfiltered, paste0(opt$output_path, opt$aliquot_barcode, "_matched_unfiltered_svaba_aa.tsv"), na = "")
