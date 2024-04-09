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
# then +1 to these columns to match with svaba result
#
################################################################################
matching_aa_svaba <- function(x){
  
  aa_bp <- x
  
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


facing_aabp <- facing_aabp %>% mutate(Svaba_bp = paste0(Chr_bp_front, ":", (bp_front+1)),
                                      Svaba_bp_mate = paste0(Chr, ":", (bp_back+1))) %>% 
  select(-Chr_bp_front)

# Extract contig info. from 'svaba.sv.vcf' 
contig.svaba <- sapply(str_split(svaba.sv.vcf$INFO, ";"), function(x){ x[grep("SCTG", x)] }) %>% 
  cbind(svaba.sv.vcf %>% select(`#CHROM`, POS, ALT), SCTG = .)

contig.svaba <- contig.svaba %>% 
  mutate(Svaba_bp_mate = gsub("\\[|\\]", "", contig.svaba$ALT) %>% 
           gsub('[ATCG]', "", .))
contig.svaba$`#CHROM` <- as.character(contig.svaba$`#CHROM`)

Matched_AA_Svaba <- left_join(facing_aabp,
                              contig.svaba %>% mutate(Svaba_bp = paste0(`#CHROM`, ":", POS)) %>% 
                                select(Svaba_bp, Svaba_bp_mate, SCTG),
                              by=c("Svaba_bp", "Svaba_bp_mate"))
  return(Matched_AA_Svaba)
  
}

Matched_annotated_cycle_Svaba <- matching_aa_svaba(annotated_cycle_bp)
Matched_cycle_Svaba <- matching_aa_svaba(cycle_bp)


matching_aa_svaba_unfiltered <- function(x){
  
  aa_bp <- x
  
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


facing_aabp <- facing_aabp %>% mutate(Svaba_bp = paste0(Chr_bp_front, ":", (bp_front+1)),
                                      Svaba_bp_mate = paste0(Chr, ":", (bp_back+1))) %>% 
  select(-Chr_bp_front)

# Extract contig info. from 'svaba.unfiltered.sv.vcf' 
contig.svaba <- sapply(str_split(svaba.unfiltered.sv.vcf$INFO, ";"), function(x){ x[grep("SCTG", x)] }) %>% 
  cbind(svaba.unfiltered.sv.vcf %>% select(`#CHROM`, POS, ALT), SCTG = .)

contig.svaba <- contig.svaba %>% 
  mutate(Svaba_bp_mate = gsub("\\[|\\]", "", contig.svaba$ALT) %>% 
           gsub('[ATCG]', "", .))
contig.svaba$`#CHROM` <- as.character(contig.svaba$`#CHROM`)

Matched_AA_Svaba <- left_join(facing_aabp,
                              contig.svaba %>% mutate(Svaba_bp = paste0(`#CHROM`, ":", POS)) %>% 
                                select(Svaba_bp, Svaba_bp_mate, SCTG),
                              by=c("Svaba_bp", "Svaba_bp_mate"))
  return(Matched_AA_Svaba)
  
}
Matched_annotated_cycle_Svaba_unfiltered <- matching_aa_svaba_unfiltered(annotated_cycle_bp)
Matched_cycle_Svaba_unfiltered <- matching_aa_svaba_unfiltered(cycle_bp)


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

