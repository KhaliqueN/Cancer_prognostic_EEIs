##############################################################################################
# Purpose: For a selected KEGG pathway, explore how the EEIs chnages across patients
##############################################################################################

rm(list=ls())
library(data.table)
# library(KEGGREST)
# library('biomaRt')
# library(Rcpp)
# library(RColorBrewer)
# library(ggplot2)
# library(GenomicDataCommons)
source("eein_cancer_util.r")

save_dir <- '../results_rep'
dir.create(save_dir)

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
cpm_threshold <- 0.5
pval_thres <- 0.05
allnets <- gtools::mixedsort(list.files('../data/CRPES',full.names=TRUE))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')

gnet <- data.table::fread(paste0('../data/PISA_survival_filt/PISA_net_final_',cpm_threshold,'.txt'), header=FALSE)
gnet <- mapProtein(gnet[[1]], gnet[[2]], data.table::fread('../data/final_EEINs/PISA.txt'))

enet <- data.table::fread(paste0('../data/EPPIC_survival_filt/EPPIC_net_final_',cpm_threshold,'.txt'), header=FALSE)
enet <- mapProtein(enet[[1]], enet[[2]], data.table::fread('../data/final_EEINs/EPPIC.txt'))

anet <- data.table::fread(paste0('../data/CONTACT_survival_filt/CONTACT_net_final_',cpm_threshold,'.txt'), header=FALSE)
anet <- mapProtein(anet[[1]], anet[[2]], data.table::fread('../data/final_EEINs/CONTACT.txt'))

aq <- rbind(gnet, anet)
unet <- rbind(aq, enet)


##--- download genes of "Pathways in cancer" 'hsa05200' "P53 signaling" 'hsa04115' "DNA replication" 'hsa03030' "cell cycle" 'hsa04110'
## "hedgehog signaling pathway" 'hsa04340' "VEGF signaling" 'hsa04370'
## "Sphingolipid signaling" 'hsa04071'
all_kegg <- data.table::fread('../data/all_human_pathways.txt', sep='\t')
kegg_id <- 'hsa04071'
pic <- all_kegg[all_kegg$Pathway_ID == kegg_id,]


##-------- module analysis -------------------------
for(qq in 3:length(allnets)){

    tempnet <- data.table::fread(allnets[qq], header=FALSE)
    tempnet_p <- mapProtein(tempnet[[1]], tempnet[[2]], unet)
    wh1 <- which(tempnet_p[[3]] %in% pic[[2]])
    wh2 <- which(tempnet_p[[4]] %in% pic[[2]])
    wh <- intersect(wh1, wh2)
    tempnet_f <- tempnet_p[wh,]

    ##--- For each cancer type ----
    # all_pateints <- list()
    # all_fracs <- list()

    for(kn in 1:length(cancer_type)){

        c_type <- cancer_type[kn]

        ##--- get survival results ------------
        gained_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Gained_Surv_.txt'))
        gained_surv$qvalue <- p.adjust(gained_surv$pval, 'fdr')
        gained_surv1 <- gained_surv[gained_surv$pval <= pval_thres, ]
        gs <- gained_surv1[,c(1,2,6)]
        ##-------------------------------------
        
        ##--- get survival results ------------
        lost_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Lost_Surv_.txt'))
        lost_surv1 <- lost_surv[lost_surv$pval <= pval_thres, ]
        ls <- lost_surv1[,c(1,2,6)]
        ##-------------------------------------
        gs_ls_u <- rbind(gs,ls)
        gs_ls_u <- igraph::graph_from_data_frame(gs_ls_u, directed=FALSE)
        gs_ls_u <- igraph::simplify(gs_ls_u, edge.attr.comb=list(weight="sum"))
        tempg <- igraph::as_data_frame(gs_ls_u)
        tempg$patient <- tempg[[3]]*(gained_surv1[[7]][1]/gained_surv1[[6]][1])

        ## intersection between pathway edges in the background network and the perturbed edges
        patients <- rep(0,length(tempnet_f[[1]]))
        fracs <- rep(0,length(tempnet_f[[1]]))
        counter <- 0
        for(k in 1:length(tempnet_f[[1]])){
            counter <- counter+1
            wh1 <- which(tempg[[1]] %in% tempnet_f[[1]][k])
            wh2 <- which(tempg[[2]] %in% tempnet_f[[2]][k])
            wh <- intersect(wh1, wh2)
            if(length(wh) == 0){
                wh1 <- which(tempg[[1]] %in% tempnet_f[[2]][k])
                wh2 <- which(tempg[[2]] %in% tempnet_f[[1]][k])
                wh <- intersect(wh1, wh2)
            }
            if(length(wh) == 1){
                patients[counter] <-  tempg[[4]][wh]
                fracs[counter] <-  tempg[[3]][wh]
            }else if(length(wh) == 0){
                next
            }else{
                print('Error')
                print(k)
                break
            }
        }

        tempnet_f <- cbind(tempnet_f, data.frame(ceiling(patients)))
        # tempnet_f <- cbind(tempnet_f, data.frame(fracs))


    }
    
    colnames(tempnet_f) <- c('Exon1','Exon2','Protein1','Protein2', cancer_type)
    data.table::fwrite(tempnet_f,paste0(save_dir,'/',kegg_id,'.csv'))

}
