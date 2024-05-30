##############################################################################################
# Purpose: Save excel files with different properties 
##############################################################################################

rm(list=ls())
library(data.table)
library('biomaRt')
library(Rcpp)
library(RColorBrewer)
library(ggplot2)
library(GenomicDataCommons)
source("eein_cancer_util.r")

cppFunction("List getpropa(CharacterVector ex1, CharacterVector ex2, CharacterVector exon1, CharacterVector exon2, NumericVector prop){

    int loop1 = ex1.size();
    int loop2 = exon1.size();
    NumericVector tprop;
    int demo=100000;

    for(int k=0; k<loop1; k++){

        int flag = 0;
        for(int j=0; j<loop2; j++){

            if((ex1[k] == exon1[j]) & (ex2[k] == exon2[j])){
                tprop.push_back(prop[j]);
                flag = 1;
                break;
            }

            if((ex1[k] == exon2[j]) & (ex2[k] == exon1[j])){
                tprop.push_back(prop[j]);
                flag = 1;
                break;
            }
        }

        if(flag == 0){
            tprop.push_back(demo);
        }

    }

    List L = List::create(tprop);
    return L;
  
}")


save_dir <- '../results'
dir.create(save_dir)

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
cpm_threshold <- 0.5
allnets <- gtools::mixedsort(list.files('../data/CRPES',full.names=TRUE))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')
tab_num <- c('C', 'B', 'A')

gnet <- data.table::fread(paste0('../data/PISA_survival_filt/PISA_net_final_',cpm_threshold,'.txt'), header=FALSE)
tempn <- data.table::fread('../data/final_EEINs/PISA.txt')
tempn <- tempn[,c(3,4,1,2)]
gnet <- mapProtein(gnet[[1]], gnet[[2]], tempn)

enet <- data.table::fread(paste0('../data/EPPIC_survival_filt/EPPIC_net_final_',cpm_threshold,'.txt'), header=FALSE)
tempn <- data.table::fread('../data/final_EEINs/EPPIC.txt')
tempn <- tempn[,c(3,4,1,2)]
enet <- mapProtein(enet[[1]], enet[[2]], tempn)

anet <- data.table::fread(paste0('../data/CONTACT_survival_filt/CONTACT_net_final_',cpm_threshold,'.txt'), header=FALSE)
tempn <- data.table::fread('../data/final_EEINs/CONTACT.txt')
tempn <- tempn[,c(3,4,1,2)]
anet <- mapProtein(anet[[1]], anet[[2]], tempn)

aq <- rbind(gnet, anet)
unet <- rbind(aq, enet)



##-------- for survival analysis edges -------------------------
for(qq in 1:length(allnets)){

    indir <- paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold)
    in_dir <- paste0('../data/Final_weighted_networks_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1])

    #-- save excel file ---
    wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S1_',strsplit(basename(net_type[qq]),
        '[.]')[[1]][1],'.xlsx'))
    # openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_1_',strsplit(basename(net_type[qq]),
    #     '[.]')[[1]][1],'.xlsx'), overwrite=T)

    for(k in 1:length(cancer_type)){

        c_type <- cancer_type[k]
        
        gained <- data.table::fread(paste0(indir,'/',c_type,'_Gained_Surv_.txt'))
        gained <- gained[gained$pval <= 0.05, ]
        lost <- data.table::fread(paste0(indir,'/',c_type,'_Lost_Surv_.txt'))
        lost <- lost[lost$pval <= 0.05, ]

        only_gained <- igraph::as_data_frame(igraph::difference(igraph::graph_from_data_frame(gained[,c(1,2)], directed=FALSE), igraph::graph_from_data_frame(lost[,c(1,2)], directed=FALSE)))
        only_lost <- igraph::as_data_frame(igraph::difference(igraph::graph_from_data_frame(lost[,c(1,2)], directed=FALSE), igraph::graph_from_data_frame(gained[,c(1,2)], directed=FALSE)))
        gained_lost <- igraph::as_data_frame(igraph::intersection(igraph::graph_from_data_frame(lost[,c(1,2)], directed=FALSE), igraph::graph_from_data_frame(gained[,c(1,2)], directed=FALSE)))

        og_pv <- getpropa(only_gained[[1]], only_gained[[2]], gained[[1]], gained[[2]], gained[[3]])
        only_gained$pvalue <- og_pv[[1]]
        og_pt <- getpropa(only_gained[[1]], only_gained[[2]], gained[[1]], gained[[2]], gained[[7]])
        only_gained$gained_in_patient <- og_pt[[1]]
        only_gained$lost_in_patient <- rep(0,length(only_gained[[1]]))

        og_pv <- getpropa(only_lost[[1]], only_lost[[2]], lost[[1]], lost[[2]], lost[[3]])
        only_lost$pvalue <- og_pv[[1]]
        og_pt <- getpropa(only_lost[[1]], only_lost[[2]], lost[[1]], lost[[2]], lost[[7]])
        only_lost$gained_in_patient <- rep(0,length(only_lost[[1]]))
        only_lost$lost_in_patient <- og_pt[[1]]

        og_pv <- getpropa(gained_lost[[1]], gained_lost[[2]], lost[[1]], lost[[2]], lost[[3]])
        gained_lost$pvalue <- og_pv[[1]]
        og_pt <- getpropa(gained_lost[[1]], gained_lost[[2]], gained[[1]], gained[[2]], gained[[7]])
        gained_lost$gained_in_patient <- og_pt[[1]]
        og_pt <- getpropa(gained_lost[[1]], gained_lost[[2]], lost[[1]], lost[[2]], lost[[7]])
        gained_lost$lost_in_patient <- og_pt[[1]]

        ##------ e) 
        t2 <- rbind(only_gained, only_lost)
        t3 <- rbind(t2, gained_lost)
        t5 <- t3[order(t3$pvalue), ]
        xx <- mapProtein(t5[[1]], t5[[2]], unet)
        xx$pvalue <- t5$pvalue
        xx$gained_in_patient <- t5$gained_in_patient
        xx$lost_in_patient <- t5$lost_in_patient

        # wb <- openxlsx::loadWorkbook(paste0(outdir,'/Supplementary_Table_1_',strsplit(basename(net_type[qq]),
        #     '[.]')[[1]][1],'.xlsx'))
        openxlsx::addWorksheet(wb1, sheetName = c_type)
        openxlsx::writeData(wb1, sheet = c_type, xx)
        openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S2',tab_num[qq],'.xlsx'), overwrite = T)
        
    }
}
