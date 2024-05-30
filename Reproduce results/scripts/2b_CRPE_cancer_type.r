##############################################################################################
# Purpose: Get perturbed EEIs and PPIs
## Splicing events from --> Kahles et al. Comprehensive Analysis of Alternative Splicing Across Tumors from 8,705 Patients, Cancer Cell, Volume 34, Issue 2, 2018
## Mutations from -->  Pan-cancer analysis of whole genomes (https://www.nature.com/articles/s41586-020-1969-6#MOESM3) supplementary table S1
## 
################################################################################################

rm(list=ls())
library(data.table)
library(ggplot2)
library(RColorBrewer)
source("eein_cancer_util.r")

save_dir <- '../results'
dir.create(save_dir)

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))

cpm_threshold <- 0.5
pval_thres <- 0.05
allnets <- gtools::mixedsort(list.files('../data/CRPES',full.names=TRUE))

net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')
fig_num <- c('Supplementary_Fig_S7','Supplementary_Fig_S6','Figure_3')
fig_numm <- c('Supplementary_Fig_S9','Supplementary_Fig_S8','Figure_4')
fig_nummx <- c('Supplementary_Fig_S11','Supplementary_Fig_S10','Figure_5')

gnet <- data.table::fread(paste0('../data/PISA_survival_filt/PISA_net_final_',cpm_threshold,'.txt'), header=FALSE)
gnet <- mapProtein(gnet[[1]], gnet[[2]], data.table::fread('../data/final_EEINs/PISA.txt'))

enet <- data.table::fread(paste0('../data/EPPIC_survival_filt/EPPIC_net_final_',cpm_threshold,'.txt'), header=FALSE)
enet <- mapProtein(enet[[1]], enet[[2]], data.table::fread('../data/final_EEINs/EPPIC.txt'))

anet <- data.table::fread(paste0('../data/CONTACT_survival_filt/CONTACT_net_final_',cpm_threshold,'.txt'), header=FALSE)
anet <- mapProtein(anet[[1]], anet[[2]], data.table::fread('../data/final_EEINs/CONTACT.txt'))

aq <- rbind(gnet, anet)
unet <- rbind(aq, enet)

temp_x_limit <- 5000

for(qq in 1:length(allnets)){

    if(qq==1){
        temp_bar_limit <- 1000
        temp_y_limit <- 4000
        temp_bar_limitl <- 500
    }else if(qq==2){
        temp_y_limit <- 3000
        temp_bar_limitl <- 300
        temp_bar_limit <- 1000
    }else{
        temp_bar_limitl <- 50
        temp_y_limit <- 1400
        temp_bar_limit <- 500
    }

    ppi_count <- c()
    eei_count <- c()
    ppi_max <- c()
    perturbed_graph <- list()

    xx <- c() ## number of eeis
    ff <- c() ## number of ppis
    ss <- c() ## Fraction of all PPIs
    Cancer <- c()
    eeis <- list()
    ppis <- list()
    patient_count <- as.data.frame(matrix(nrow=0,ncol=2))
    mdn <- c()
    mxp <- c()


    net_file <- igraph::graph_from_data_frame(data.table::fread(allnets[qq]), directed=FALSE)
    net_filed <- igraph::as_data_frame(net_file)
    netx <- mapProtein(net_filed[[1]], net_filed[[2]], unet)

    for(gg in 1:length(cancer_type)){

        c_type <- cancer_type[gg]
        ##--- get survival results ------------
        gained_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Gained_Surv_.txt'))
        gained_surv$qvalue <- p.adjust(gained_surv$pval, 'fdr')
        gained_surv1 <- gained_surv[gained_surv$pval <= pval_thres, ]
        gs <- gained_surv1[,c(1,2,6)]
        # mdn <- c(mdn, gained_surv1[[7]][1]/gained_surv1[[6]][1])
        # gs <- igraph::graph_from_data_frame(gained_surv1[,c(1,2,6)], directed=FALSE)
        ##-------------------------------------
        
        ##--- get survival results ------------
        lost_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Lost_Surv_.txt'))
        lost_surv1 <- lost_surv[lost_surv$pval <= pval_thres, ]
        ls <- lost_surv1[,c(1,2,6)]
        # ls <- igraph::graph_from_data_frame(lost_surv1[,c(1,2,6)], directed=FALSE)
        ##-------------------------------------
        gs_ls_u <- rbind(gs,ls)
        gs_ls_u <- igraph::graph_from_data_frame(gs_ls_u, directed=FALSE)
        gs_ls_u <- igraph::simplify(gs_ls_u, edge.attr.comb=list(weight="sum"))
        # gs_ls_u <- igraph::simplify(gs_ls_u, edge.attr.comb=list(patient="sum"))

        # igraph::union(ls, gs)
        eeis[[gg]] <- gs_ls_u
        tempg <- igraph::as_data_frame(gs_ls_u)
        mdn <- c(mdn, median(tempg[[3]]*(gained_surv1[[7]][1]/gained_surv1[[6]][1])))
        mxp <- c(mxp, max(100*tempg[[3]]))

        ## store patient counts of perturbed edges
        patient_count <- rbind(patient_count,data.frame(cancer=rep(cancer_type[gg],length(tempg[[3]])), 100*tempg[[3]]))
        # patient_count <- rbind(patient_count,data.frame(cancer=rep(cancer_type[gg],length(tempg[[3]])), tempg[[3]]))

        p1 <- c()
        p2 <- c()

        for(i in 1:length(tempg[[1]])){
            wh1 <- which(unet$exon1 == tempg[[1]][i])
            wh2 <- which(unet$exon2 == tempg[[2]][i])
            wha <- intersect(wh1, wh2)
            if(length(wha) > 0){
                tempx <- unet[wha, ]
                p1 <- c(p1, tempx$protein1[1])
                p2 <- c(p2, tempx$protein2[1])
                next
            }
            wh1 <- which(unet$exon2 == tempg[[1]][i])
            wh2 <- which(unet$exon1 == tempg[[2]][i])
            whb <- intersect(wh1, wh2)
            tempx <- unet[whb, ]
            p2 <- c(p2, tempx$protein1[1])
            p1 <- c(p1, tempx$protein2[1])
        }

        tempg$protein1 <- p1
        tempg$protein2 <- p2

        ppig <- igraph::graph_from_data_frame(tempg[,c(3,4)], directed=FALSE)
        ppi_graph <- igraph::simplify(ppig)
        ppis[[gg]] <- ppi_graph
        ppi_count <- c(ppi_count, igraph::ecount(ppi_graph))

        igraph::E(ppig)$weight <- 1
        ppig1 <- igraph::simplify(ppig, edge.attr.comb=list(weight="sum"))
        ppig2 <- igraph::as_data_frame(ppig1)
        ppig2 <- ppig2[order(-ppig2$weight), ]
        ppi_max <- c(ppi_max, ppig2$weight[1])

        eei_graph <- igraph::simplify(igraph::graph_from_data_frame(tempg[,c(1,2)], directed=FALSE))
        eei_count <- c(eei_count, igraph::ecount(eei_graph))

        ## store for plot -----
        ppig2_count <- plyr::count(ppig2$weight)
        ppig2_count$pert <- ppig2_count$freq/sum(ppig2_count$freq)
        xx <- c(xx, ppig2_count$x)
        ff <- c(ff, ppig2_count$freq)
        ss <- c(ss,ppig2_count$pert)
        Cancer <- c(Cancer, rep(cancer_type[gg], length(ppig2_count[[1]])))
        perturbed_graph[[gg]] <- gs_ls_u

    }

   
    ##----- number of cancer-relevant edges plot -----------
    pdata <- data.frame(Cancer=cancer_type, count=eei_count)
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdata, aes(Cancer, count)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,(max(pdata$count))+temp_bar_limit)) +
    geom_text(aes(label=count), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0(save_dir,"/",fig_num[qq],"A.png"),width=3.5, height=3, dpi=400)


    ### number of patients to which the perturbed edges belong -------
    pdatax <- patient_count
    pmdn <- data.frame(Cancer=cancer_type,count=signif(mdn,1),pos=mxp+8)
    colnames(pdatax) <- c('Cancer', 'count')
    p <- ggplot(pdatax, aes(Cancer, count, fill=Cancer)) + 
    geom_boxplot(alpha=0.9)+
    # geom_point(aes(fill = Cancer), shape = 21, position = position_jitterdodge()) +
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="% of patients \nwith paired samples", limits=c(0,110)) +
    geom_text(data=pmdn, aes(y=pos, x=Cancer,label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0.5, angle=0, size=3)+
    scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
        axis.text.y = element_text(size = 8, angle = 0, colour = "black"))+#,
        # panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))
    ggsave(p,filename=paste0(save_dir,"/",fig_num[qq],"B.png"),width=4.5, height=2.5, dpi=400)


    ### ---------------------- correlation with known number of cancer type-specific splicing events ------------------
    pdata <- pdata[order(pdata$Cancer), ]

    
    ###----- AS events and survival ----
    splicing <- readxl::read_excel('../data/oncosplicing.xlsx', 1)
    splicing <- as.data.frame(splicing)
    colnames(splicing) <- as.vector(splicing[1,])
    splicing <- splicing[-1,]
    splicing <- splicing[1:31,]
    splicing1 <- splicing[,c(1,5,8)]
    splicing1 <- splicing1[splicing1$`TCGA cancer` %in% cancer_type, ]
    splicing1 <- splicing1[order(splicing1$`TCGA cancer`), ]
    pdata$overall_splicing <- as.numeric(splicing1[[2]])
    pdata$splicing <- as.numeric(splicing1[[3]])
    wh <- which(!is.na(pdata$splicing))
    pdatax <- pdata[wh,]

    coral <- cor.test(x=pdatax$count, y=pdatax$splicing, method = 'spearman')
    p <- ggplot(pdatax, aes(splicing, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=12))
    basesize <- 10
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of survival associated \nalternative splicing events", limits=c(0,max(pdatax$splicing))) + 
    geom_text(aes(x=8500,y=temp_y_limit, label=paste0('Spearman correlation: ',signif(coral$estimate[[1]],3),'\np-value: ',signif(coral$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdatax$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','black','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=3))
    ggsave(p,filename=paste0(save_dir,"/",fig_numm[qq],"A.png"),width=4.5, height=2.5, dpi=400)


    ### ---------------------- correlation with known numbers of genes significantly correlated with survival ------------------
    survival <- data.table::fread('../data/survival_analysis_all.txt', sep='\t')
    survs <- survival[,-1]

    signif_can <- sapply(survs, function(x) length(which(x > 1.96 | x < -1.96))) ## 3.1 in case of the 0.001 p-value
    survival_data <- data.frame(Cancer=names(signif_can), count=signif_can)
    survival_data <- survival_data[order(survival_data$Cancer), ]
    pdata$genes <- survival_data$count


    coral1 <- cor.test(x=pdata$count, y=pdata$genes, method = 'spearman')
    # coral1 <- cor.test(x=pdata$count, y=pdata$genes, method = 'pearson')

    p <- ggplot(pdata, aes(genes, count, color=Cancer)) + 
    geom_point(size=2)+
    theme(legend.text=element_text(size=10))
    basesize <- 10
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_continuous(name="# of genes significantly \nassociated with survival", limits=c(0,max(pdata$genes))) + 
    geom_text(aes(x=5000,y=temp_y_limit, label=paste0('Spearman correlation: ',signif(coral1$estimate[[1]],3),'\np-value: ',signif(coral1$p.value,3))), size=3,show.legend = FALSE)+
    scale_y_continuous(name="# of cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,max(pdata$count))) +
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.8, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.8, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(color=guide_legend(title="Cancer type",ncol=3))
    ggsave(p,filename=paste0(save_dir,"/",fig_num[qq],"C.png"),width=4.5, height=2.5, dpi=400)

    

    # ###----------------------------------------------------------------
    # ###----- Unique CRPEs ------------------------------------------
    unique_eeis <- list()
    ueeis <- c()
    for(k in 1:length(cancer_type)){
        counter <- k
        temp_uni <- eeis[[k]]

        for(j in 1:length(cancer_type)){
            if(counter != j){
                temp_uni <- igraph::difference(temp_uni, eeis[[j]])
            }
        }
        unique_eeis[[k]] <- temp_uni
        ueeis <- c(ueeis, igraph::ecount(temp_uni))
    }

    ##----- number of cancer-relevant edges plot -----------
    pdata$ueeis <- ueeis
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdata, aes(Cancer, ueeis)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of unique cancer-relevant \nperturbed edges (CRPEs)", limits=c(0,(max(pdata$ueeis))+temp_bar_limitl)) +
    geom_text(aes(label=ueeis), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0(save_dir,"/",fig_nummx[qq],"A.png"),width=3.5, height=3, dpi=400)



    ####-------------------------------------------------------------------

    ##--- figure for the overlap of CRPEs ---------------------------
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
            c1 <- c(c1, paste0(cancer_type[ov1],' (', eei_count[ov1],')'))
            c2 <- c(c2, paste0(cancer_type[ov2],' (', eei_count[ov2],')'))
            e1 <- igraph::ecount(perturbed_graph[[ov1]])
            e2 <- igraph::ecount(perturbed_graph[[ov2]])
            inter_e <- igraph::intersection(perturbed_graph[[ov1]], perturbed_graph[[ov2]])
            if(e1 < e2){
                tempo <- igraph::ecount(inter_e)/e1
            }else{
                tempo <- igraph::ecount(inter_e)/e2
            }
            oval <- c(oval, igraph::ecount(inter_e))
            # oval <- c(oval, tempo)

            hyp <- phyper(igraph::ecount(inter_e)-1,e1,igraph::ecount(net_file)-e1,e2,lower.tail = FALSE)
            hyp <- signif(hyp, digits=3)
            pval <- c(pval, hyp)

            ## under representation
            hyp <- phyper(igraph::ecount(inter_e),e1,igraph::ecount(net_file)-e1,e2,lower.tail = TRUE)
            hyp <- signif(hyp, digits=3)
            pvalu <- c(pvalu, hyp)
        # }
        }
    }

    qval <- p.adjust(pval, 'fdr')
    qvalu <- p.adjust(pvalu, 'fdr')
    qvalx <- rep('Overlap as \nexpected by chance',length(pval))
    for(i in 1:length(pval)){
        if(qval[i] <= pval_thres){
            qvalx[i] <- 'Higher overlap than\nexpected by chance'
        }
        if(qvalu[i] <= pval_thres){
            qvalx[i] <- 'Lower overlap than\nexpected by chance'
        }
    }

    pdata <- data.frame(c1=c1, c2=c2, val=oval, pval=qvalx)


    # pdata <- data.frame(c1=c1, c2=c2, val=signif(oval*100,2), qval=signif(qval, 2))
    # pdata <- data.frame(c1=c1, c2=c2, val=oval, qval=signif(qval, 2))
    # pdata$pvalue <- cut(pdata$qval, breaks = c(0,0.01, 0.05, 1), include.lowest=TRUE)

    cols <- rev(c('#377eb8','#4daf4a','#e41a1c'))#rev(brewer.pal(3,"Spectral"))
    p <- ggplot(pdata, aes(c1, c2)) + geom_tile(aes(fill = pval),colour = "white")+
      theme(legend.text=element_text(size=8))+scale_fill_manual(values=cols,drop=FALSE)
    basesize <- 8
    p <- p + theme_grey(base_size = basesize) + labs(x = "Cancer type", y = "Cancer type") +
      scale_y_discrete() +
      scale_x_discrete()+coord_fixed()+
      guides(fill=guide_legend(title="Overlap category"))+
      geom_text(data=pdata,aes(y=c2,label=val),size=1.8)+
      theme(axis.text.x = element_text(size = basesize * 0.8,angle = 90, hjust = 0,vjust=0, colour = "black"),
        axis.text.y = element_text(size = basesize * 0.8,angle = 0, hjust = 0,vjust=0, colour = "black"))
    ggsave(p,filename=paste0(save_dir,'/',fig_nummx[qq],'B.png'),width=4.5, height=3, dpi=300)



}




