##############################################################################################
# Purpose: KM plots for top edge that has at least 10 patients in the two groups (hig vs. low survival)
##############################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
source("eein_cancer_util.r")
library(GenomicDataCommons)
library("survival")
library("survminer")

save_dir <- '../results/survival'
dir.create(save_dir, recursive=TRUE)

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
cpm_threshold <- 0.5
pval_thres <- 0.05
allnets <- gtools::mixedsort(list.files('../data/CRPES',full.names=TRUE))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')


for(qq in 3:length(allnets)){

    biocrpes <- readRDS(paste0('../data/',net_type[qq],"_biocrpe.Rds"))


	for(k in 1:length(cancer_type)){

		c_type <- cancer_type[k]
        out_dir1 <- paste0(save_dir,'/',c_type)
        dir.create(out_dir1)

        netn <- strsplit(basename(allnets[qq]),'[.]')[[1]][1]
        all_exr <- data.table::fread(paste0('../data/filtered_normalized_data/', c_type,'/normalized_all.txt'), header=TRUE)
        all_sr <- get_survival(c_type, 0)

        ##--- get survival results ------------
        gained_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],
            '/threshold_',cpm_threshold,'/',c_type,'_Gained_Surv_.txt'))
        lost_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
            c_type,'_Lost_Surv_.txt'))

        bioc <- igraph::as_data_frame(biocrpes[[k]])

        if(nrow(bioc) > 0){

            gc <- rep(0, length(bioc[[1]]))
            lc <- rep(0, length(bioc[[1]]))
            pval <- rep(1, length(bioc[[1]]))

            for(i in 1:length(bioc[[1]])){

                wh1 <- which(gained_surv[[1]] == bioc[[1]][i])
                wh2 <- which(gained_surv[[2]] == bioc[[2]][i])
                wha <- intersect(wh1, wh2)
                wh1 <- which(gained_surv[[2]] == bioc[[1]][i])
                wh2 <- which(gained_surv[[1]] == bioc[[2]][i])
                whb <- intersect(wh1, wh2)
                whc <- union(wha, whb)
                if(length(whc) != 0){
                    gs1 <- gained_surv[whc, ]
                    pval[i] <- gs1$pval[1]
                    gc[i] <- gs1$patient[1]
                }

                wh1 <- which(lost_surv[[1]] == bioc[[1]][i])
                wh2 <- which(lost_surv[[2]] == bioc[[2]][i])
                wha <- intersect(wh1, wh2)
                wh1 <- which(lost_surv[[2]] == bioc[[1]][i])
                wh2 <- which(lost_surv[[1]] == bioc[[2]][i])
                whb <- intersect(wh1, wh2)
                whc <- union(wha, whb)
                if(length(whc) != 0){
                    ls1 <- lost_surv[whc, ]
                    pval[i] <- ls1$pval[1]
                    lc[i] <- ls1$patient[1]
                }

            }

            biocd <- data.frame(bioc,pval,gc,lc)
            biocd <- biocd[order(biocd$pval), ]

            for(i in 1:length(biocd[[1]])){
                survival_plot_new(c_type, biocd[i,c(1,2)], cpm_threshold, all_exr, all_sr, out_dir1, paste0(net_type[qq],'_',i), pval)   
            }

        }   

	}

}


