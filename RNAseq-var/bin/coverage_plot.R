#!/usr/bin/env Rscript

library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

con_path <- args[1]
sample_id <- args[2]
out_path <- args[3]

# test ####
# con_path <- './RNAseq_CHOHD_varscan2/HD-N85-F30/HD-N85-F30.targeted.hqaa.mpileup.consensus.vcf'
# sample_id <- 'HD-N85-F30'
# out_path <- './RNAseq_CHOHD_varscan2/HD-N85-F30'


con <- read.delim (con_path)
con_r1 <- con[, c('Position', 'Reads1', 'VarFreq')]
con_r1$Type <- 'ref'
con_r2 <-  con[, c('Position', 'Reads2', 'VarFreq')]
con_r2$Type <- 'var'
names(con_r2) <- names (con_r1)
df <- rbind (con_r1, con_r2)
p <- ggplot(data=df, aes(x=Position, y=Reads1, group=Type, color=Type)) +
  geom_line(size=1) + theme_classic()  +scale_color_manual(values = c('deepskyblue', 'red')) + 
  labs (title = sample_id) + ylab ("number of reads")

png (paste0(out_path, '/', sample_id, '.coverage_plot.png'), width = 500, height = 250)
print (p)
dev.off()