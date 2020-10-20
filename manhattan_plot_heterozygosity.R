setwd("~/Dropbox/postdoc cornell/projects/manuscripts/meme_genomenote_ms")
library(data.table)
library(qqman)
var_win<-fread("shotgun_variants_windows.txt", header = FALSE)
colnames(var_win)<-c("chr","winstart", "winend","var")
var_win$chr_medaka[var_win$chr == 1 ] <- 3
var_win$chr_medaka[var_win$chr == 2 ] <- 9
var_win$chr_medaka[var_win$chr == 3 ] <- 6
var_win$chr_medaka[var_win$chr == 4 ] <- 5
var_win$chr_medaka[var_win$chr == 5 ] <- 10
var_win$chr_medaka[var_win$chr == 6 ] <- 13
var_win$chr_medaka[var_win$chr == 7 ] <- 16
var_win$chr_medaka[var_win$chr == 8 ] <- 15
var_win$chr_medaka[var_win$chr == 9 ] <- 7
var_win$chr_medaka[var_win$chr == 10 ] <- 4
var_win$chr_medaka[var_win$chr == 11 ] <- 21
var_win$chr_medaka[var_win$chr == 12 ] <- 22
var_win$chr_medaka[var_win$chr == 13 ] <- 12
var_win$chr_medaka[var_win$chr == 14 ] <- 14
var_win$chr_medaka[var_win$chr == 15 ] <- 11
var_win$chr_medaka[var_win$chr == 16 ] <- 20
var_win$chr_medaka[var_win$chr == 17 ] <- 17
var_win$chr_medaka[var_win$chr == 18 ] <- 8
var_win$chr_medaka[var_win$chr == 19 ] <- 19
var_win$chr_medaka[var_win$chr == 20 ] <- 23
var_win$chr_medaka[var_win$chr == 21 ] <- 24.1
var_win$chr_medaka[var_win$chr == 22 ] <- 1.1
var_win$chr_medaka[var_win$chr == 23 ] <- 18
var_win$chr_medaka[var_win$chr == 24 ] <- 2
var_win$chr_medaka[var_win$chr == 25 ] <- 1.2
var_win$chr_medaka[var_win$chr == 26 ] <- 1.3
var_win$chr_medaka[var_win$chr == 27 ] <- 24.2
var_win_27<-var_win[chr<28]
var_win_order<-var_win_27[order(chr_medaka, winstart)]
var_win_order$het<-var_win_order$var/50000*100

var_win_ct<-fread("10x_variants_windows.txt", header = FALSE)
colnames(var_win_ct)<-c("chr","winstart", "winend","var")
var_win_ct$chr_medaka[var_win_ct$chr == 1 ] <- 3
var_win_ct$chr_medaka[var_win_ct$chr == 2 ] <- 9
var_win_ct$chr_medaka[var_win_ct$chr == 3 ] <- 6
var_win_ct$chr_medaka[var_win_ct$chr == 4 ] <- 5
var_win_ct$chr_medaka[var_win_ct$chr == 5 ] <- 10
var_win_ct$chr_medaka[var_win_ct$chr == 6 ] <- 13
var_win_ct$chr_medaka[var_win_ct$chr == 7 ] <- 16
var_win_ct$chr_medaka[var_win_ct$chr == 8 ] <- 15
var_win_ct$chr_medaka[var_win_ct$chr == 9 ] <- 7
var_win_ct$chr_medaka[var_win_ct$chr == 10 ] <- 4
var_win_ct$chr_medaka[var_win_ct$chr == 11 ] <- 21
var_win_ct$chr_medaka[var_win_ct$chr == 12 ] <- 22
var_win_ct$chr_medaka[var_win_ct$chr == 13 ] <- 12
var_win_ct$chr_medaka[var_win_ct$chr == 14 ] <- 14
var_win_ct$chr_medaka[var_win_ct$chr == 15 ] <- 11
var_win_ct$chr_medaka[var_win_ct$chr == 16 ] <- 20
var_win_ct$chr_medaka[var_win_ct$chr == 17 ] <- 17
var_win_ct$chr_medaka[var_win_ct$chr == 18 ] <- 8
var_win_ct$chr_medaka[var_win_ct$chr == 19 ] <- 19
var_win_ct$chr_medaka[var_win_ct$chr == 20 ] <- 23
var_win_ct$chr_medaka[var_win_ct$chr == 21 ] <- 24.1
var_win_ct$chr_medaka[var_win_ct$chr == 22 ] <- 1.1
var_win_ct$chr_medaka[var_win_ct$chr == 23 ] <- 18
var_win_ct$chr_medaka[var_win_ct$chr == 24 ] <- 2
var_win_ct$chr_medaka[var_win_ct$chr == 25 ] <- 1.2
var_win_ct$chr_medaka[var_win_ct$chr == 26 ] <- 1.3
var_win_ct$chr_medaka[var_win_ct$chr == 27 ] <- 24.2
var_win_ct_27<-var_win_ct[chr<28]
var_win_ct_order<-var_win_ct_27[order(chr_medaka, winstart)]
var_win_ct_order$het<-var_win_ct_order$var/50000*100
manhattan(var_win_ct_order,chr="chr_medaka",bp="winstart",p="het",logp=FALSE,ylab="Heterozygosity (%)", ylim = c(0, 4.5), col=c("grey","lightskyblue2"), main="Connecticut", cex.axis=1.2, cex.lab=1.5, cex.main=1.5, cex=0.5)
manhattan(var_win_order,chr="chr_medaka",bp="winstart",p="het",logp=FALSE,ylab="Heterozygosity (%)", ylim = c(0, 4.5), col=c("grey","orange1"), main="Georgia", cex.axis=1.2, cex.lab=1.5, cex.main=1.5, cex=0.5)
