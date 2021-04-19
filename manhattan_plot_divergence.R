library(data.table)
library(gridExtra)
library(tidyverse)
#library(dplyr)
###this is the output from nucmer, which I used to estimate divergence
div_mum<-fread("mummerCT2GA.filter10k.tab", header=T)
div_mum$chr_ref[div_mum$CHR1 == "chr1a" ] <- 1.1
div_mum$chr_ref[div_mum$CHR1 == "chr1b" ] <- 1.2
div_mum$chr_ref[div_mum$CHR1 == "chr1c" ] <- 1.3
div_mum$chr_ref[div_mum$CHR1 == "chr2" ] <- 2
div_mum$chr_ref[div_mum$CHR1 == "chr3" ] <- 3
div_mum$chr_ref[div_mum$CHR1 == "chr4" ] <- 4
div_mum$chr_ref[div_mum$CHR1 == "chr5" ] <- 5
div_mum$chr_ref[div_mum$CHR1 == "chr6" ] <- 6
div_mum$chr_ref[div_mum$CHR1 == "chr7" ] <- 7
div_mum$chr_ref[div_mum$CHR1 == "chr8" ] <- 8
div_mum$chr_ref[div_mum$CHR1 == "chr9" ] <- 9
div_mum$chr_ref[div_mum$CHR1 == "chr10" ] <- 10
div_mum$chr_ref[div_mum$CHR1 == "chr11" ] <- 11
div_mum$chr_ref[div_mum$CHR1 == "chr12" ] <- 12
div_mum$chr_ref[div_mum$CHR1 == "chr13" ] <- 13
div_mum$chr_ref[div_mum$CHR1 == "chr14" ] <- 14
div_mum$chr_ref[div_mum$CHR1 == "chr15" ] <- 15
div_mum$chr_ref[div_mum$CHR1 == "chr16" ] <- 16
div_mum$chr_ref[div_mum$CHR1 == "chr17" ] <- 17
div_mum$chr_ref[div_mum$CHR1 == "chr18" ] <- 18
div_mum$chr_ref[div_mum$CHR1 == "chr19" ] <- 19
div_mum$chr_ref[div_mum$CHR1 == "chr20" ] <- 20
div_mum$chr_ref[div_mum$CHR1 == "chr21" ] <- 21
div_mum$chr_ref[div_mum$CHR1 == "chr22" ] <- 22
div_mum$chr_ref[div_mum$CHR1 == "chr23" ] <- 23
div_mum$chr_ref[div_mum$CHR1 == "chr24a" ] <- 24.1
div_mum$chr_ref[div_mum$CHR1 == "chr24b" ] <- 24.2
div_mum$mb<- div_mum$S1/1000000
div_mum$DIV<-100-div_mum$IDY

### Subsample data.table by chromosome
chr8mum<-div_mum[CHR1 == 'chr8']
chr11mum<-div_mum[CHR1 == 'chr11']
chr18mum<-div_mum[CHR1 == 'chr18']
chr24mum<-div_mum[CHR1 == 'chr24a']
chr10mum<-div_mum[CHR1 == 'chr10']


### These plots are to explore the patterns of divergence across each of the 4 chromosomes with haploblocks (aka large inversions)
plot_chr8<-ggplot(chr8mum, aes(x=S1, y=DIV)) + geom_point() +ylim(c(0,13)) + ggtitle('chr8') + 
  geom_hline(yintercept=2.159173, linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(2332150,14917775), size=0.5) 
plot_chr11<-ggplot(chr11mum, aes(x=S1, y=DIV)) + geom_point() +ylim(c(0,13)) + ggtitle('chr11') + 
  geom_hline(yintercept=2.159173, linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(5918,1196215,8180880,9364918,1681812,2737423), size=0.5) 
plot_chr18<-ggplot(chr18mum, aes(x=S1, y=DIV)) + geom_point() +ylim(c(0,13)) + ggtitle('chr18') + 
  geom_hline(yintercept=2.159173, linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(2224993,9657870), size=0.5) 
plot_chr24<-ggplot(chr24mum, aes(x=S1, y=DIV)) + geom_point() +ylim(c(0,13)) + ggtitle('chr24') + 
  geom_hline(yintercept=2.159173, linetype="dashed", size=0.5) +
  geom_vline(xintercept=c(722031,10135403), size=0.5) 

### example of chromosome without big inversions
plot_chr10<-ggplot(chr10mum, aes(x=S1, y=DIV)) + geom_point() +ylim(c(0,13)) + ggtitle('chr10') + 
  geom_hline(yintercept=2.159173, linetype="dashed", size=0.5) 

grid.arrange(plot_chr8,plot_chr11,plot_chr18,plot_chr24,plot_chr10,nrow=5)

### Subsample inversion regions
chr8inv<-chr8mum[S1 >2332150 & S1 < 14917775]
weighted.mean(chr8inv$DIV,chr8inv$LEN2)

chr11inva<-chr11mum[S1 > 5918 & S1 < 1196215]
chr11invb<-chr11mum[S1 > 8180880 & S1 < 9364918]
chr11invc<-chr11mum[S1 > 1681812 & S1 < 2737423]
chr11inv<-rbind(chr11inva,chr11invb,chr11invc)
weighted.mean(chr11inv$DIV,chr11inv$LEN2)

chr18inv<-chr18mum[S1 >2224993 & S1 < 9657870]
weighted.mean(chr18inv$DIV,chr18inv$LEN2)

chr24inv<-chr24mum[S1 > 722031 & S1 < 10135403]
weighted.mean(chr24inv$DIV,chr24inv$LEN2)
### data manipulation to have all inversions together
inv<-rbind(chr8inv, chr11inv, chr18inv, chr24inv)
weighted.mean(inv$DIV,inv$LEN2)


### plot of divergence within inversions
require(plyr)
dd <- ddply(inv, .(CHR1), summarise, m=mean(DIV), wm=weighted.mean(DIV,LEN2))
chrs <- c("chr8", "chr11", "chr18", "chr24a")
inv %>% 
  mutate(CHR1 = factor(CHR1, levels = chrs)) %>% 
  ggplot(aes( CHR1,DIV)) + geom_boxplot() + geom_jitter(aes(col=CHR1), size=0.5) + geom_hline(yintercept=2.159173, linetype="dashed", size=0.5) + 
  geom_point(dd,mapping=aes(x=CHR1, y=wm, colour=CHR1), size=5) +
  theme_bw() + xlab('Inversions') + ylab("Sequence divergence (%)")

### Subsample non-inversion regions
noninv<-setdiff(div_mum, inv)
new_var <- "chr"
inv[, (new_var):=CHR1]
noninv$chr <- 'noninv'

### Combine inversion and non-inversion regions
div_mum_inv<-rbind(inv,noninv)
### set groups and colors
chrs <- c("noninv", "chr8", "chr11", "chr18", "chr24a")

dd_all <- ddply(div_mum_inv, .(chr), summarise, m=mean(DIV), wm=weighted.mean(DIV,LEN2))


### Final plot - all chromosomes
colchr<-c("#CCCCCC","#6D8470","#C96555","#AD6428","#999999")
library(ggforce)
library(tidyverse)
div_mum_inv %>% 
  mutate(chr = factor(chr, levels = chrs)) %>% 
  ggplot(aes( chr,DIV, color=chr))  + 
  scale_colour_manual(values = colchr) +
  scale_fill_manual(values = colchr) +
  geom_violin(alpha = 0) +
  geom_sina( size=0.5, alpha=0.1) +
 # geom_jitter( size=0.5, alpha=0.2) +
  geom_point(dd_all,mapping=aes(x=chr, y=wm, col=chr), size=4, shape=23, fill="white") +
  theme_classic() + xlab('Inversions') + ylab("Sequence divergence (%)")  +
  theme(legend.position = "none")

### Violin plots  comparing each chr with inversion to the rest of the genome - one plot for each chromosome
chrs_8 <- c("noninv", "chr8")
chrs_11 <- c("noninv", "chr11")
chrs_18 <- c("noninv",  "chr18")
chrs_24 <- c("noninv",  "chr24a")

dd_all_8 <- dd_all %>% 
  filter(chr %in% c('noninv', 'chr8')) 
dd_all_11 <- dd_all %>% 
  filter(chr %in% c('noninv', 'chr11')) 
dd_all_18 <- dd_all %>% 
  filter(chr %in% c('noninv', 'chr18')) 
dd_all_24 <- dd_all %>% 
  filter(chr %in% c('noninv', 'chr24a')) 

colchr_sub<-c("#CCCCCC","mediumorchid4")

div_mum_inv %>% 
  mutate(chr = factor(chr, levels = chrs_8)) %>% 
  filter(chr %in% c('noninv', 'chr8')) %>%
  ggplot(aes( chr,DIV, color=chr))  + ylim(c(0,10.5)) +
  scale_colour_manual(values = colchr_sub) +
  scale_fill_manual(values = colchr_sub) +
  geom_violin(alpha = 0) +
  geom_sina( size=0.5, alpha=0.1) +
  # geom_jitter( size=0.5, alpha=0.2) +
  geom_point(dd_all_8,mapping=aes(x=chr, y=wm, col=chr), size=4, shape=23, fill="white") +
  theme_classic()  +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks = c(0, 5, 10, 15))


div_mum_inv %>% 
  mutate(chr = factor(chr, levels = chrs_11)) %>% 
  filter(chr %in% c('noninv', 'chr11')) %>%
  ggplot(aes( chr,DIV, color=chr))  + ylim(c(0,10.5)) +
  scale_colour_manual(values = colchr_sub) +
  scale_fill_manual(values = colchr_sub) +
  geom_violin(alpha = 0) +
  geom_sina( size=0.5, alpha=0.1) +
  # geom_jitter( size=0.5, alpha=0.2) +
  geom_point(dd_all_11,mapping=aes(x=chr, y=wm, col=chr), size=4, shape=23, fill="white") +
  theme_classic()  +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())+
  scale_y_continuous(breaks = c(0, 5, 10, 15))


div_mum_inv %>% 
  mutate(chr = factor(chr, levels = chrs_18)) %>% 
  filter(chr %in% c('noninv', 'chr18')) %>%
  ggplot(aes( chr,DIV, color=chr))  +  ylim(c(0,10.5)) +
  scale_colour_manual(values = colchr_sub) +
  scale_fill_manual(values = colchr_sub) +
  geom_violin(alpha = 0) +
  geom_sina( size=0.5, alpha=0.1) +
  # geom_jitter( size=0.5, alpha=0.2) +
  geom_point(dd_all_18,mapping=aes(x=chr, y=wm, col=chr), size=4, shape=23, fill="white") +
  theme_classic()  +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank()) +
  scale_y_continuous(breaks = c(0, 5, 10, 15))


div_mum_inv %>% 
  mutate(chr = factor(chr, levels = chrs_24)) %>% 
  filter(chr %in% c('noninv', 'chr24a')) %>%
  ggplot(aes( chr,DIV, color=chr))  + ylim(c(0,10.5)) +
  scale_colour_manual(values = colchr_sub) +
  scale_fill_manual(values = colchr_sub) +
  geom_violin(alpha = 0) +
  geom_sina( size=0.5, alpha=0.1) +
  # geom_jitter( size=0.5, alpha=0.2) +
  geom_point(dd_all_24,mapping=aes(x=chr, y=wm, col=chr), size=4, shape=23, fill="white") +
  theme_classic()  +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())  +
  scale_y_continuous(breaks = c(0, 5, 10, 15))




### Manhattan plot for divergence (12 x 2.5 inch)
### All chromosomes
div_col<-rep(c("grey","mediumorchid4"), length.out = 27)
ggplot(data = div_mum,
       aes(x = S1, 
           y = DIV,
           col = as.factor(chr_ref))) +
  geom_point(size=0.2) +
  scale_color_manual(values = div_col) +
  # remove space between plot area and x axis
  scale_y_continuous(expand = c(0, 0.1)) +
  # facet by CHR
  facet_grid(cols = vars(chr_ref),
             space = "free_x",
             scales = "free_x") +
  labs(x = "Chromosome", y="Divergence (%)") + theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.1, "cm"), # adjust spacing between facets
    strip.background = element_blank(),
    strip.text.x = element_text(size=8),
    legend.position = "none",
    #axis.line.x = element_blank(),
    #axis.line.x.top = element_line(color="black")
  )

### zoom in on 4 chromosomes with haploblocks
div8<-div_mum %>%
  filter(chr_ref == 8) %>%
  ggplot(aes(mb,DIV)) + geom_point(size=0.9, col="mediumorchid4") + ylim(c(0, 11.5)) +  xlim(c(0, 17)) +
  theme_classic()  + xlab("chromosome 8 (Mb)") + ylab("Divergence (%)") +
  theme(axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  geom_vline(xintercept=c(2.332150,14.917775), size=0.5, linetype="dashed", col="gray") + geom_hline(yintercept=2.159173, col="gray", size=0.5) +
  scale_x_continuous(breaks = c(0, 5, 10, 15))

div11<-div_mum %>%
  filter(chr_ref == 11) %>%
  ggplot(aes(mb,DIV)) + geom_point(size=0.9, col="mediumorchid4") + ylim(c(0, 11.5)) +  
  theme_classic()  + xlab("chromosome 11 (Mb)") + ylab("") +
  theme(axis.title.y = element_text(size=14),
        axis.text.x = element_text(size=12, color="black"),
       axis.text.y = element_text(size=12, color="black"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))+
  geom_vline(xintercept=c(0.005918,1.196215,8.180880,9.364918,1.681812,2.737423), size=0.5, linetype="dashed", col="gray") + geom_hline(yintercept=2.159173, col="gray", size=0.5)+
  scale_x_continuous(breaks = c(0, 5, 10, 15))


div18<-div_mum %>%
  filter(chr_ref == 18) %>%
  ggplot(aes(mb,DIV)) + geom_point(size=0.9, col="mediumorchid4") + ylim(c(0, 11.5)) +  
  theme_classic()  + xlab("chromosome 18 (Mb)") + ylab("") +
  theme(axis.title.y = element_text(size=14),
   axis.text.x = element_text(size=12, color="black"),
  axis.text.y = element_text(size=12, color="black"),
  legend.text = element_text(size=12),
  legend.title = element_text(size=14))+
  geom_vline(xintercept=c(2.224993,9.657870), size=0.5, linetype="dashed", col="gray") + geom_hline(yintercept=2.159173, col="gray", size=0.5)+
  scale_x_continuous(breaks = c(0, 5, 10, 15))


div24<-div_mum %>%
  filter(chr_ref == 24.1) %>%
  ggplot(aes(mb,DIV)) + geom_point(size=0.9, col="mediumorchid4") + ylim(c(0, 11.5)) +  
  theme_classic()  + xlab("chromosome 24.1 (Mb)") + ylab("") +
  theme(axis.title.y = element_text(size=14),
   axis.text.x = element_text(size=12, color="black"),
   axis.text.y = element_text(size=12, color="black"),
  legend.text = element_text(size=12),
  legend.title = element_text(size=14))+
  geom_vline(xintercept=c(7.22031,1.0135403), size=0.5, linetype="dashed", col="gray") + geom_hline(yintercept=2.159173, col="gray", size=0.5) +
  scale_x_continuous(breaks = c(0, 5, 10, 15))

grid.arrange(div8,div11, div18,div24,nrow=1)
