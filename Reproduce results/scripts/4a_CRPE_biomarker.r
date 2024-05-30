##############################################################################################
# Purpose: numbers of BioCRPEs
##############################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
library(RColorBrewer)
source("eein_cancer_util.r")
library(GenomicDataCommons)
library("survival")
library("survminer")
library('biomaRt')

save_dir <- '../results'
dir.create(save_dir)

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
cpm_threshold <- 0.5
pval_thres <- 0.05
allnets <- gtools::mixedsort(list.files('../data/CRPES',full.names=TRUE))
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')
fig_num <- c('Supplementary_Fig_S18','Supplementary_Fig_S17','Figure_7')
tab_num <- c('C', 'B', 'A')

##--- map to survival analysis by Smith et al -----
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# attrs <- c('hgnc_symbol', 'uniprotswissprot', 'description')
survival <- data.table::fread('../data/survival_analysis_all.txt', sep='\t')
allc <- colnames(survival)[-1]
allgenes <- survival[[1]]## 20,531 genes
# alld <- getBM(attributes=attrs, filters='hgnc_symbol', values=allgenes, mart=ensembl) 
# alldd <- alld[alld$uniprotswissprot != '',]
# data.table::fwrite(alldd, '../data/gene_map.txt', row.names=FALSE, sep='\t', quote=FALSE)
alldd <- data.table::fread('../data/gene_map.txt')
alldd <- alldd[,-3]
# colnames(alldd) <- c('uniprotswissprot','hgnc_symbol')
alldd <- alldd[alldd$hgnc_symbol != '',]
survival <- survival[,-1]

save_bioCRPEs <- function(ls1, gs1, c_type){

    ls2 <- ls1[ls1$frac1 == 0, ]
    ls2 <- ls2[ls2$frac2 == 0, ]
    ls <- igraph::graph_from_data_frame(ls2[,c(1,2)], directed=FALSE)

    gs2 <- gs1[gs1$frac1 == 0, ]
    gs2 <- gs2[gs2$frac2 == 0, ]
    gs <- igraph::graph_from_data_frame(gs2[,c(1,2)], directed=FALSE)

    ps <- igraph::as_data_frame(igraph::union(ls, gs))

    ## get survival info ----
    all_exr <- data.table::fread(paste0('../data/filtered_normalized_data/', c_type,'/normalized_all.txt'), header=TRUE)
    all_sr <- get_survival(c_type, 0)

    if(nrow(ps) > 0){
        mpflag <- c()
        for(j in 1:length(ps[[1]])){
            xx <- survival_get_direction(c_type, ps[j,c(1,2)], cpm_threshold, all_exr, all_sr) 
            mpflag <- c(mpflag, xx)
        }

        pval <- c()
        for(j in 1:length(ps[[1]])){
            wh1 <- which(ls2[[1]] == ps[[1]][j])
            wh2 <- which(ls2[[2]] == ps[[2]][j])
            wha <- intersect(wh1, wh2)
            wh1 <- which(ls2[[1]] == ps[[2]][j])
            wh2 <- which(ls2[[2]] == ps[[1]][j])
            whb <- intersect(wh1, wh2)
            whx <- union(wha, whb)
            if(length(whx) == 1){
                pval <- c(pval, ls2[[3]][whx])
            }else{
                wh1 <- which(gs2[[1]] == ps[[1]][j])
                wh2 <- which(gs2[[2]] == ps[[2]][j])
                wha <- intersect(wh1, wh2)
                wh1 <- which(gs2[[1]] == ps[[2]][j])
                wh2 <- which(gs2[[2]] == ps[[1]][j])
                whb <- intersect(wh1, wh2)
                whx <- union(wha, whb)
                pval <- c(pval, gs2[[3]][whx])
            }
        }

        xps <- mapProtein(ps[[1]], ps[[2]], unet)
        xps$pval <- pval
        xps$`Survival outcome` <- mpflag
        xps <- xps[order(xps$pval),]

        ##--- add info from Smith et al. study ---
        wh <- which(allc == c_type)
        surv_pv <- as.data.frame(survival)[,wh]
        pr1 <- c()
        pr2 <- c()
        gene1 <- c()
        gene2 <- c()
        for(j in 1:length(ps[[1]])){
            p1 <- alldd[alldd$uniprotswissprot == xps$protein1[j],]
            if(nrow(p1) != 0){
                pv1 <- surv_pv[which(allgenes == p1$hgnc_symbol[1])]
                gene1 <- c(gene1, p1$hgnc_symbol[1])
                if(pv1 < -1.96){
                    pr1 <- c(pr1, 'Favorable')
                }else if(pv1 > 1.96){
                    pr1 <- c(pr1, 'Unfavorable')
                }else{
                    pr1 <- c(pr1, 'Unknown')
                }
            }else{
                pr1 <- c(pr1, 'Symbol')
                gene1 <- c(gene1, 'Symbol')
            }

            p2 <- alldd[alldd$uniprotswissprot == xps$protein2[j],]
            if(nrow(p2) != 0){
                pv2 <- surv_pv[which(allgenes == p2$hgnc_symbol[1])]
                gene2 <- c(gene2, p2$hgnc_symbol[1])
                if(pv2 < -1.96){
                    pr2 <- c(pr2, 'Favorable')
                }else if(pv2 > 1.96){
                    pr2 <- c(pr2, 'Unfavorable')
                }else{
                    pr2 <- c(pr2, 'Unknown')
                }
            }else{
                pr2 <- c(pr2, 'Symbol')
                gene2 <- c(gene2, 'Symbol')
            }
            
        }

        xps$`protein1 survival outcome` <- pr1
        xps$`protein2 survival outcome` <- pr2

        xps$`gene1` <- gene1
        xps$`gene2` <- gene2

    }else{
        xps <- ps
    }
    
    return(xps)

}


###----------------------------------------------------

qval_thres <- c(0)#, 0.01, 0.02, 0.03, 0.04, 0.05)

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


for(qq in 1:length(allnets)){

	in_dir <- paste0('../data/Final_weighted_networks_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1])
	net_file <- data.table::fread(allnets[qq], header=FALSE)
    
    pert_bio <- c()
    cancer <- c()
    q_cutoff <- c()
    fav <- c()
    unfav <- c()

    unique_bio <- list()
    counter <- 1
    wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S5',tab_num[qq],'.xlsx'))

	for(gg in 1:length(cancer_type)){
		c_type <- cancer_type[gg]
        indir <- paste0(in_dir,'/threshold_',cpm_threshold)
        gained <- data.table::fread(paste0(indir,'/',c_type,'_gained.txt'))
        lost <- data.table::fread(paste0(indir,'/',c_type,'_lost.txt'))
        ##-- a) Non-perturbed edges --------------
        nonp1 <- gained[gained$patient == 0, ]
        nonp1 <- igraph::graph_from_data_frame(nonp1[,c(1,2)], directed=FALSE)
        nonp2 <- lost[lost$patient == 0, ]
        nonp2 <- igraph::graph_from_data_frame(nonp2[,c(1,2)], directed=FALSE)
        nonp <- igraph::intersection(nonp1, nonp2)
        ##-- e) f) Patient-spcific gained and lost ----
        ps1 <- gained[gained$patient == 1, ]
        ps1 <- igraph::graph_from_data_frame(ps1[,c(1,2)], directed=FALSE)
        lost_ig <- igraph::graph_from_data_frame(lost[lost$patient > 0, ][,c(1,2)], directed=FALSE)
        ps1 <- igraph::difference(ps1, lost_ig) ## e) patient-specific gain
        ps2 <- lost[lost$patient == 1, ]
        ps2 <- igraph::graph_from_data_frame(ps2[,c(1,2)], directed=FALSE)
        gained_ig <- igraph::graph_from_data_frame(gained[gained$patient > 0, ][,c(1,2)], directed=FALSE)
        ps2 <- igraph::difference(ps2, gained_ig) ## f) patient-specific loss
        ##-- b) c) only gained/ only lost
        og <- gained[gained$patient > 1, ]
        og <- igraph::graph_from_data_frame(og[,c(1,2)], directed=FALSE)
        og <- igraph::difference(og, lost_ig) ## b) only gained --. gained in at least 2, but never lost
		ol <- lost[lost$patient > 1, ]
        ol <- igraph::graph_from_data_frame(ol[,c(1,2)], directed=FALSE)
        ol <- igraph::difference(ol, gained_ig)
        ##-- d) gained/lost
        gl <- igraph::intersection(lost_ig, gained_ig) ## d) gained/lost edges
        xx <- igraph::ecount(gl)+igraph::ecount(ol)+igraph::ecount(og)+igraph::ecount(ps1)+igraph::ecount(ps2)+igraph::ecount(nonp)

        gained_surv <- data.table::fread(paste0('../data/Final_survival_filt/',
            strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',c_type,'_Gained_Surv_.txt'))
        gained_surv1 <- gained_surv[gained_surv$pval <= pval_thres, ]
        
        lost_surv <- data.table::fread(paste0('../data/Final_survival_filt/',
            strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',c_type,'_Lost_Surv_.txt'))
        lost_surv1 <- lost_surv[lost_surv$pval <= pval_thres, ]

        # for(i in 1:length(qval_thres)){
        # if(qval_thres[i] == 0){## if the qval == 0 then save the exon pairs ---
        sbio <- save_bioCRPEs(lost_surv1, gained_surv1, cancer_type[gg])
        fav <- c(fav, length(which(sbio$`Survival outcome` == 'Favorable')))
        unfav <- c(unfav, length(which(sbio$`Survival outcome` == 'Unfavorable')))
        # }

        # pvl <- qval_thres[i]
        pvl <- qval_thres[1]
        lost_surv2 <- lost_surv1[lost_surv1$frac1 <= pvl, ]
        lost_surv2 <- lost_surv2[lost_surv2$frac2 <= pvl, ]
        ls <- igraph::graph_from_data_frame(lost_surv2[,c(1,2)], directed=FALSE)

        gained_surv2 <- gained_surv1[gained_surv1$frac1 <= pvl, ]
        gained_surv2 <- gained_surv2[gained_surv2$frac2 <= pvl, ]
        gs <- igraph::graph_from_data_frame(gained_surv2[,c(1,2)], directed=FALSE)

        ps <- igraph::union(ls, gs)

        pert_bio <- c(pert_bio, igraph::ecount(ps))
        q_cutoff <- c(q_cutoff, pvl)
        cancer <- c(cancer, c_type)

        ##---extra---
        psg <- igraph::as_data_frame(ps)
        temp_net <- mapProtein(psg[[1]], psg[[2]], unet)
        # if(i == 1){
        unique_bio[[counter]] <- ps
        counter <- counter+1
        # }
        # }
        ## save excel sheet ----
        openxlsx::addWorksheet(wb1, sheetName = cancer_type[gg])
        openxlsx::writeData(wb1, sheet = cancer_type[gg], sbio)
        openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S5',tab_num[qq],'.xlsx'), overwrite = T)
	}

    pdata <- data.frame(Cancer=cancer, cutoff=as.factor(q_cutoff), biom=pert_bio)

    ##--- number of biomarkers figure -------
    pdatax <- pdata[pdata$cutoff == 0, ]
    pdatax$Favorable <- fav
    pdatax$Unfavorable <- unfav
    pdatay <- pdatax[,-c(2,3)]
    pdataz <- reshape2::melt(pdatay)
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdataz, aes(Cancer, value, fill=variable)) + 
    geom_bar(stat="identity",position=position_stack())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of biomarker cancer-relevant \nperturbed edges (BioCRPEs)", limits=c(0,(max(pdatax$biom))+5)) +
    # geom_text(aes(label=biom), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill=guide_legend(title="Survival outcome"))
    ggsave(p,filename=paste0(save_dir,"/",fig_num[qq],"A.png"),width=5, height=3, dpi=400)





    ##---- extra analysis ----
    allg <- unique_bio[[1]]
    for(kg in 2:length(unique_bio)){
        allg <- igraph::union(allg, unique_bio[[kg]])
    }


    background_net <- length(net_file[[1]])
    perturbed_graph <- unique_bio
    ##--- figure for the overlap ---------------------------
    c1 <- c()
    c2 <- c()
    oval <- c()
    pval <- c()
    pvalu <- c()
    loop1 <- length(cancer_type)-1
    loop2 <- length(cancer_type)

    for(ov1 in 1:loop1){
        m <- ov1+1
        for(ov2 in m:loop2){
            # if(ov1 != ov2){
            e1 <- igraph::ecount(perturbed_graph[[ov1]])
            e2 <- igraph::ecount(perturbed_graph[[ov2]])
            c1 <- c(c1, paste0(cancer_type[ov1],' (', e1,')'))
            c2 <- c(c2, paste0(cancer_type[ov2],' (', e2,')'))
            inter_e <- igraph::intersection(perturbed_graph[[ov1]], perturbed_graph[[ov2]])
            oval <- c(oval, igraph::ecount(inter_e))

            hyp <- phyper(igraph::ecount(inter_e)-1,e1,background_net-e1,e2,lower.tail = FALSE)
            hyp <- signif(hyp, digits=3)
            pval <- c(pval, hyp)

            ## under representation
            hyp <- phyper(igraph::ecount(inter_e),e1,background_net-e1,e2,lower.tail = TRUE)
            hyp <- signif(hyp, digits=3)
            pvalu <- c(pvalu, hyp)
        # }
        }
    }

    qval <- p.adjust(pval, 'fdr')
    qvalu <- p.adjust(pvalu, 'fdr')
    qvalx <- rep('Expected by chance',length(pval))
    for(i in 1:length(pval)){
        if(qval[i] <= pval_thres){
            qvalx[i] <- 'Significantly \nmore overlapping'
        }
        if(qvalu[i] <= pval_thres){
            qvalx[i] <- 'Significantly \nless overlapping'
        }
    }


    pdata <- data.frame(c1=c1, c2=c2, val=oval, pval=qvalx)

    cols <- c('#377eb8','#4daf4a','#e41a1c')#rev(brewer.pal(3,"Spectral"))
    p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
      theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
    basesize <- 8
    p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
      scale_y_discrete() +
      scale_x_discrete()+coord_fixed()+
      guides(fill=guide_legend(title="Q-value"))+
      geom_text(data=pdata,aes(y=c2,label=val),size=2.5)+
      theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0.5, colour = "black"),
        axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0.5, colour = "black"))+
      guides(fill='none')
    ggsave(p,filename=paste0(save_dir,'/',fig_num[qq],'B.png'),width=2.5, height=2.5, dpi=400)

    ##-------


    #--- save the BioCRPEs -------
    saveRDS(unique_bio, file = paste0(save_dir,'/',net_type[qq],"_biocrpe.Rds"))
    # xx <- readRDS(paste0(net_type[qq],"_biocrpe.Rds"))


    ##--- category of novel biocrpes -----
    if(qq == 3){

        known <- c()
        nknown <- c()
        unknown <- c()
        knownp <- c()
        nknownp <- c()
        unknownp <- c()
        x <- c()
        for(j in 1:length(cancer_type)){
            biocrpes <- openxlsx::read.xlsx('../data/Supplementary_Table_S4_processed.xlsx',j+1)
            x <- c(x, nrow(biocrpes))
            if(nrow(biocrpes) != 0){
                wh1 <- which(biocrpes[[7]] == 'Unknown')
                wh2 <- which(biocrpes[[8]] == 'Unknown')
                wha <- length(intersect(wh1, wh2))## both unknown
                whb <- length(setdiff(union(wh1, wh2), intersect(wh1, wh2))) ## exactly one unknown
                whc <- nrow(biocrpes)-wha-whb ## both known

                known <- c(known, whc)
                nknown <- c(nknown, whb)
                unknown <- c(unknown, wha)
                knownp <- c(knownp, whc/nrow(biocrpes))
                nknownp <- c(nknownp, whb/nrow(biocrpes))
                unknownp <- c(unknownp, wha/nrow(biocrpes))
            }else{
                known <- c(known, 0)
                nknown <- c(nknown, 0)
                unknown <- c(unknown,0)
                knownp <- c(knownp, 0)
                nknownp <- c(nknownp, 0)
                unknownp <- c(unknownp, 0)
            }
        }

        ## category plot
        pData <- data.frame(A=knownp*100, B=nknownp*100, C=unknownp*100, X=x, Cancer=cancer_type)

        pData1 <- reshape2::melt(pData,id=c('Cancer','X'))
        pData1$variable <- factor(pData1$variable, levels=c("A", "B", "C"))  # setting level order
        basesize <- 8
        ppx <- ggplot(data = pData1, aes(x=Cancer, y=value, fill=variable, group=Cancer)) + geom_bar(stat="identity")+
        geom_text(aes(label=X, y=100),size=2.5, vjust=0.5, hjust=0, angle=60)+
        scale_y_continuous(limits=c(0,105), breaks = seq(0, 100, by = 10))+
        xlab("Cancer type")+ylab("% of biomarker cancer-relevant \nperturbed edges (BioCRPEs)")+
        scale_fill_manual(labels=c("A" = "Both genes", 
            "B"="Only one gene", "C"="None of the genes"), 
        values=c('#fee6ce','#fdae6b','#e6550d'))+
        theme_bw()+theme(axis.text.x = element_text(size = 1*basesize, angle = 60, vjust=1, hjust=1, colour = "black"),
            axis.text.y = element_text(size = 1*basesize, angle = 0, colour = "black"),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"))+
        guides(fill=guide_legend(title="Genes of a BioCRPE known\n to be survival biomarkers"))
        ggsave(ppx,filename=paste0(save_dir, "/Supplementary_Fig_S19.png"),width=7, height=3, dpi=400)


    }
   


}
