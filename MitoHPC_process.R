## For MitoHPC summary process
## Cleaned up version

# 031825
# For subset v6

## If running entire script, run this in command line first in your MitoHPC output folder
#cat mutect2.05.suspicious.tab | grep haplocheck > sus.txt

## packages
#install.packages("tidyverse")
#install.packages("data.table")

## Load
library(tidyverse)
library(data.table)

## Read in vcf
counts = read.table('out/mutect2.mutect2.05.vcf')

## Change col names for MitoHPC outputs
colnames(counts) = c('CHR','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE')

## Get info from INFO column (get AF, genotype, read depth information from vcf file)
test = unlist(strsplit(counts$INFO, split=';'))
AF = test[grepl('AF',test)]
GT = test[grepl('GT',test)]
DP = test[grepl('DP',test)]
HG = test[grepl('HG',test)]
counts$AF = as.numeric(substr(AF,4,15))
counts$Genotype = substr(GT,4,15)
counts$Read_depth = as.numeric(substr(DP,4,15))
counts$Haplogroup = NA

## Filter out het vs homo, and apply filters to het only, to not filter out homo alt variants
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]
counts = countshet

## Filter for variant level FILTER
counts = counts %>% filter(!str_detect(FILTER,"strict_strand|strand_bias|base_qual|map_qual|weak_evidence|slippage|position|Homopolymer|clustered|fragment|haplotype"))

## Filter for Read_counts < 300
counts = counts[counts$Read_depth >= 300,]

## merge
counts = rbind(counts,countshomo)
rm(countshet,countshomo,AF,DP,GT,HG,test)

## Filter for variant level INFO (remove INDEL, HP region)
counts = counts[!grepl('INDEL',counts$INFO),]
counts = counts[!grepl('Homopolymer',counts$INFO),]

## unique (the unique combination of POS_REF_ALT for each variant)
counts$unique = paste(paste(counts$POS, counts$REF, sep = '_'), counts$ALT, sep = '_')

## merge in mito annotation
ma = as.data.frame(fread('mito_genome_annotation.txt'))
countsma = merge(counts,ma,by=c('unique','POS','REF','ALT'))
counts = countsma
rm(countsma,ma)

## Sep homo and het variants
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]

## Filter out those with two ALT allele that add up to 1 (biallelic het)
countshetbi = countshet %>% group_by(SAMPLE,POS) %>% slice_max(mlc, with_ties = FALSE) %>% ungroup()
countshet = countshetbi
rm(countshetbi)

## Calculate Median, Min and Max (for het variants)
countshet = countshet %>% group_by(unique) %>%
  mutate(AoU_het_AF_Median = median(AF), AoU_het_AF_Max = max(AF), AoU_het_AF_Min = min(AF)) %>%
  as.data.frame()
countshomo$AoU_het_AF_Median = NA
countshomo$AoU_het_AF_Max = NA
countshomo$AoU_het_AF_Min = NA
counts = rbind(countshet, countshomo)
rm(countshet,countshomo)

## Non-synonymous variants
counts$mutation_nonsynonymous = NA
counts$mutation_stop = NA
countsnons1 = counts[grepl('NONSYN',counts$INFO),]
countsnons2 = counts[!grepl('NONSYN',counts$INFO),]
countsnons1$mutation_nonsynonymous = 'NONSYN'
countsnons = rbind(countsnons1,countsnons2)
counts = countsnons

## STOP
countsstop1 = counts[grepl('STOP',counts$INFO),]
countsstop2 = counts[!grepl('STOP',counts$INFO),]
countsstop1$mutation_stop = 'STOP'
countsstop = rbind(countsstop1,countsstop2)
counts = countsstop
rm(countsstop,countsnons,countsnons1,countsnons2,countsstop1,countsstop2)

## Calculate het and homo count for AoU
counts = counts %>% group_by(unique) %>%
  mutate(count_het_AoU = sum(AF != 1),
         count_homo_AoU = sum(AF == 1)) %>%
  as.data.frame()

################### Per sample count

## per sample counts
persamplecount = counts %>% group_by(SAMPLE) %>%
  summarize(count_het = sum(AF != 1),
            count_homo = sum(AF == 1),
            MSS = sum(mlc[AF != 1]),
            mMSS = sum(mlcm[AF != 1])) %>%
  as.data.frame()

## add in samples without any het/homo
everyone = as.data.frame(fread('out/cvg.tab'))
allsamplesmissed = everyone$Run[!everyone$Run %in% persamplecount$SAMPLE]
missedsampletable = data.frame(SAMPLE = allsamplesmissed, count_het = 0, count_homo = 0, MSS = 0, mMSS = 0)
persamplecount = rbind(persamplecount, missedsampletable)
rm(missedsampletable,everyone,allsamplesmissed)

## get sus samples from tab file
# Run the following command in command line first
#cat mutect2.05.suspicious.tab | grep haplocheck > sus.txt
sus = as.data.frame(fread('out/sus.txt'))
persamplecount$contamination = ifelse(persamplecount$SAMPLE %in% sus$V1,'Yes','No')
rm(sus)

## add CN
cn = as.data.frame(fread('cn_subset_v6.txt'))
cn1 = cn %>% select(Run,'mtDNA-CN')
persamplecount = merge(persamplecount,cn1,by.x='SAMPLE',by.y='Run')
colnames(persamplecount)[7] = 'CN'
rm(cn1,cn)
persamplecount$lowCN = ifelse((persamplecount$CN <= 40),'Yes','No')

## haplogroup
haplogroup = as.data.frame(fread('out/mutect2.haplogroup1.tab'))
persamplecount = merge(persamplecount,haplogroup,by.x='SAMPLE',by.y='Run')
rm(haplogroup)

## add blood vs saliva (gsutil -u $GOOGLE_PROJECT -m cp -r gs://fc-secure-0a7fb819-21e7-4e1f-b3ff-94acc1e3b4e6/srWGS_auxiliary/qc/genomic_metrics.tsv ./)
gm = as.data.frame(fread('genomic_metrics.tsv'))
gm$research_id = paste('wgs',gm$research_id,sep='_')
persamp1 = merge(persamplecount,gm,by.x='SAMPLE',by.y='research_id',all.x=TRUE)

## per complex MSS counts (het only)
# The following code is an example for getting MSS for each COMPLEX, can modify for other information interested
persamplecomplexmss = countshet %>%
  group_by(SAMPLE, COMPLEX) %>%
  summarize(MSS = sum(mito_lc_score)) %>%
  pivot_wider(names_from = COMPLEX, values_from = MSS, values_fill = 0) %>%
  as.data.frame()
# remove column NA
persamplecomplexmss = persamplecomplexmss %>% select(-'NA')


## flag file?
flag = as.data.frame(fread('flagged_samples.tsv'))
sin = persamp1[is.na(persamp1$biosample_collection_date),]$SAMPLE
sout = persamp1[!is.na(persamp1$biosample_collection_date),]$SAMPLE
flag$sample = paste('wgs',flag$s,sep='_')
table(sout %in% flag$sample, useNA = 'always')
persamp1$flagged = ifelse(persamp1$SAMPLE %in% flag$sample,'YES','NO')
flagged = persamp1 %>% filter(flagged == 'YES')
## among the 43 in subset 16k, 36 with 0 het and 7 with 1 het, no contamination really

###### Save per samp and het/homo files
countshomo = counts[counts$AF == 1,]
countshet = counts[!counts$AF == 1,]

write.table(persamp1, file='per_sample_v6_subset_022625.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(countshet, file='het_v6_subset_022525.txt', sep="\t", row.names=FALSE, quote = FALSE)
write.table(countshomo, file='homo_v6_subset_022525.txt', sep="\t", row.names=FALSE, quote = FALSE)
#write.table(persamplecomplexmss, file='complex_MSS_v6_subset_022525.txt', sep="\t", row.names=FALSE, quote = FALSE)
