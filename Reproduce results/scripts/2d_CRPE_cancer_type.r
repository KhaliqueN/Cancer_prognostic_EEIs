 ##############################################################################################
# Purpose: Figure for biological properties 
##############################################################################################

rm(list=ls())
library(data.table)
library('biomaRt')
library(Rcpp)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(GenomicDataCommons)
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


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


for(qq in 1:length(allnets)){

    in_dir <- paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/')
    net_file <- data.table::fread(allnets[qq], header=FALSE)

    all_pert <- list()

    for(k in 1:length(cancer_type)){
        c_type <- cancer_type[k]
        gained <- data.table::fread(paste0(in_dir,c_type,'_Gained_Surv_.txt'))
        gained_ig <- igraph::graph_from_data_frame(gained[gained$pval <= pval_thres, ][,c(1,2)], directed=FALSE)
        lost <- data.table::fread(paste0(in_dir,c_type,'_Lost_Surv_.txt'))
        lost_ig <- igraph::graph_from_data_frame(lost[lost$pval <= pval_thres, ][,c(1,2)], directed=FALSE)
        all_pert[[k]] <- igraph::union(gained_ig, lost_ig)
    }

   
    ###------------------------------------------------------------------------------

    ###--- KEGG term computations ----------------------------------------------------------------------
    net_file1 <- mapProtein(net_file[[1]], net_file[[2]], unet)
    all_kegg <- data.table::fread('../data/all_human_pathways.txt', sep='\t')
    colnames(all_kegg) <- c('ensembl_gene_id','uniprot_id','entrez_id','go_id','pathways_name')
    all_keggx <- all_kegg[all_kegg$uniprot_id %in% union(net_file1[[3]], net_file1[[4]]), ]

    attrs <- c('ensembl_gene_id', 'ensembl_exon_id')
    alld <- getBM(attributes=attrs, filters='ensembl_gene_id', values=all_keggx[[1]], mart=ensembl) 

    ## only keep the exons in my background network -------------------------------------
    alldx <- alld[alld$ensembl_exon_id %in% union(net_file[[1]], net_file[[2]]), ] 

    path_id <- c()
    path_name <- c()
    exonx <- c()
    for(k in 1:length(alldx[[1]])){
        wh <- which(all_kegg$ensembl_gene_id == alldx$ensembl_gene_id[k])
        path_name <- c(path_name, all_kegg$pathways_name[wh])
        path_id <- c(path_id, all_kegg$go_id[wh])
        exonx <- c(exonx, rep(alldx$ensembl_exon_id[k], length(wh)))
    }

    alldd <- data.frame(ensembl_exon_id=exonx, path_id=path_id, pathways_name=path_name)

    all_path <- plyr::count(alldd$path_id)
    wh_path <- unique(all_path[which(all_path$freq > 2), ]$x)
    allkegg <- alldd[alldd$path_id %in% wh_path, ] ### 285 pathways remain for NETHIGH network
    # colnames(allkegg) <- c('ensembl_gene_id','ensembl_exon_id','name_1006','go_id')
    colnames(allkegg) <- c('ensembl_exon_id','go_id','name_1006')
    kg_bgSize <- length(unique(allkegg$ensembl_exon_id))

    allKE1 <- goEnrich(allkegg, kg_bgSize, all_pert, cancer_type)
    allKE1$qval <- p.adjust(allKE1$pval, 'fdr')
    allKE <- allKE1[allKE1$qval <= pval_thres, ]

    #-- save excel file ---
    wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S3_',net_type[qq],'.xlsx'))

    for(k in 1:length(cancer_type)){

        temp <- allKE[allKE$temp_cancer == cancer_type[k], ]
        temp1a <- temp[order(temp$qval), ]

        colnames(temp1a) <- c('Cancer','KEGG ID', 'pvalue', 'KEGG pathway name', 'qvalue')
        openxlsx::addWorksheet(wb1, sheetName = cancer_type[k])
        openxlsx::writeData(wb1, sheet = cancer_type[k], temp1a)
        openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S3_',net_type[qq],'.xlsx'), overwrite = T)

    }

    
    ####----- Pairwise overlaps of enriched KEGG terms ---------------------------------------------------------------------------------
    path_bgSize <- length(unique(allkegg$go_id))

    c1 <- c()
    c2 <- c()
    value <- c()
    pvalue <- c()
    pvalueu <- c()
    loop1 <- length(cancer_type)-1
    loop2 <- length(cancer_type)

    for(k in 1:loop1){

        tempx1 <- allKE[allKE$temp_cancer == cancer_type[k], ]
        m  <- k+1
        for(j in m:loop2){

            # if(k != j){
            tempx2 <- allKE[allKE$temp_cancer == cancer_type[j], ]
            e1 <- unique(tempx1$GOterm)
            e2 <- unique(tempx2$GOterm)

            ovr <- intersect(e1, e2)
            ur <- union(e1, e2)
            pv <- phyper(length(ovr)-1, length(e1), path_bgSize-length(e1), length(e2), lower.tail=FALSE)
            pvu <- phyper(length(ovr), length(e1), path_bgSize-length(e1), length(e2), lower.tail=TRUE)

            if(length(e1) < length(e2)){
                tempo <- length(ovr)/length(e1)
            }else{
                tempo <- length(ovr)/length(e2)
            }
            value <- c(value, length(ovr))

            pvalue <- c(pvalue, pv)
            pvalueu <- c(pvalueu, pvu)
            c1 <- c(c1, paste0(cancer_type[k],' (',length(e1),')'))
            c2 <- c(c2, paste0(cancer_type[j],' (',length(e2),')'))
        # }
        }
    }

    qval <- signif(p.adjust(pvalue,'fdr'),2)
    qvalu <- signif(p.adjust(pvalueu,'fdr'),2)
    qvalx <- rep('Overlap as \nexpected by chance',length(pvalue))
    for(i in 1:length(pvalue)){
        if(qval[i] <= pval_thres){
            qvalx[i] <- 'Higher overlap than\nexpected by chance'
        }
        if(qvalu[i] <= pval_thres){
            qvalx[i] <- 'Lower overlap than\nexpected by chance'
        }
    }

    pdata <- data.frame(c1=c1, c2=c2, val=value, pval=qvalx)


    cols <- rev(c('#377eb8','#e41a1c'))

    p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
      theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
    basesize <- 8
    p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
      scale_y_discrete() +
      scale_x_discrete()+coord_fixed()+
      guides(fill=guide_legend(title="Category"))+
      geom_text(data=pdata,aes(y=c2,label=val),size=2.5)+
      theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
        axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))#+
      # guides(fill='none')
    ggsave(p,filename=paste0(save_dir,'/KEGG_overlap_',net_type[qq],'.png'),width=4.5, height=3, dpi=300)



    ###---- unique KEGG terms -------
    unique_KE <- list()
    uKEs <- c()
    for(k in 1:length(cancer_type)){
        counter <- k
        temp_uni <- unique(allKE[allKE$temp_cancer == cancer_type[k], ]$GOtermn)

        for(j in 1:length(cancer_type)){
            if(counter != j){
                temp_uni <- setdiff(temp_uni, allKE[allKE$temp_cancer == cancer_type[j], ]$GOtermn)
            }
        }
        unique_KE[[k]] <- temp_uni
        uKEs <- c(uKEs, length(temp_uni))
    }

    tempg <- data.frame(cancer=cancer_type, count=uKEs)
    p <- ggplot(tempg, aes(cancer, count)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of uniquely enriched \nKEGG pathways", limits=c(0,(max(tempg$count))+2)) +
    geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0, angle=0, size=3)+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0(save_dir,"/",net_type[qq],"_unique_KEGG.png"),width=3.5, height=3, dpi=400)


  }