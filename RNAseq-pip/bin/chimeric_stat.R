#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


chimeric_path <- args[1]
htseq_path <- args[2]
sample_id <- args[3]
out_path <- args[4]

Header <- c('chr_donorA', 'brkpt_donorA', 'strand_donorA', 'chr_acceptorB', 'brkpt_acceptorB', 'strand_acceptorB', 'junction_type', 'repeat_left_lenA',
            'repeat_right_lenB', 'read_name', 'start_alnA', 'cigar_alnA', 'start_alnB', 'cigar_alnB', 'num_chim_aln', 'max_poss_aln_score', 
            'non_chim_aln_score', 'this_chim_aln_score', 'bestall_chim_aln_score', 'PEmerged_bool')
# test ####
# chimeric_path <- './adamts13_chimeric_output/HD-N85-F30.adam.chimeric.txt'
# htseq_path <- './adamts13_chimeric_output/HD-N85-F30.htseq'
# out_path <- './adamts13_chimeric_output'
# sample_id <- 'HD-N7-F12'

chimeric <- read.delim (chimeric_path, header = F)
colnames (chimeric) <- Header
chimeric <- chimeric[, c(1,2,3,4,5,6,7,10,12,14)]
chimeric$fusion_name <- paste(chimeric$chr_donorA, chimeric$brkpt_donorA, chimeric$strand_donorA, 
                              chimeric$chr_acceptorB,chimeric$brkpt_acceptorB,chimeric$strand_acceptorB, sep = '_')
chimeric_count <- data.frame(table (chimeric$fusion_name))
colnames(chimeric_count) <- c(sample_id, 'Frequency')

htseq <- read.delim(htseq_path, header = F)
count <- as.numeric(htseq[1,2])
chimeric_count$perc <- chimeric_count$Freq / count
chimeric_count <- chimeric_count [order(chimeric_count$perc, decreasing = T), ]


write.csv (chimeric_count,paste0(out_path, '/',sample_id, '.chimeric.count.csv'),row.names = F)

filter_chim <- chimeric_count [which (chimeric_count$perc > 0.001), ]
if ( dim (filter_chim) [1] > 0 ) 
  {
  x <- as.character(filter_chim[,1])
  chimeric_name <- chimeric [which(chimeric$fusion_name %in% x), 'read_name']
  write.table (chimeric_name, paste0(out_path, '/', sample_id, '.chimeric.reads.name.txt'), 
               row.names = F, col.names = F, quote = F) 
  } else {print ("Takeda pipeline: none of chimeric reads > 0.1%")}
