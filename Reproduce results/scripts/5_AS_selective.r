##############################################################################################
# Purpose: Figures for importance of EEI vs PPI
# Highlighting CRPEs
##############################################################################################

rm(list=ls())
library(data.table)
library('biomaRt')
library(Rcpp)
library(RColorBrewer)
library(ggplot2)
library(GenomicDataCommons)
library(dplyr)
source("eein_cancer_util.r")

cppFunction("List getprop(CharacterVector ex1, CharacterVector ex2, CharacterVector exon1, CharacterVector exon2, NumericVector prop){

    int loop1 = ex1.size();
    int loop2 = exon1.size();
    NumericVector tprop;

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
            tprop.push_back(0);
            //Rcout << k << std::endl;
            //break;
        }

    }

    List L = List::create(tprop);
    return L;
  
}")


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

alldd <- data.table::fread('../data/gene_map.txt')


save_impPPIs <- function(tempcc){

    gene1 <- c()
    gene2 <- c()
    for(j in 1:length(tempcc[[1]])){
        p1 <- alldd[alldd$uniprotswissprot == tempcc$protein1[j],]
        if(nrow(p1) != 0){
            gene1 <- c(gene1, p1$hgnc_symbol[1])
        }else{
            gene1 <- c(gene1, 'Symbol')
        }

        p2 <- alldd[alldd$uniprotswissprot == tempcc$protein2[j],]
        if(nrow(p2) != 0){
            gene2 <- c(gene2, p2$hgnc_symbol[1])
        }else{
            gene2 <- c(gene2, 'Symbol')
        }
    }

    tempcc$gene1 <- gene1
    tempcc$gene2 <- gene2
    return(tempcc)

}


for(qq in 1:length(allnets)){

    in_dir <- paste0('../data/Final_weighted_networks_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1])
    net_file <- data.table::fread(allnets[qq], header=FALSE)
    net_file <- mapProtein(net_file[[1]], net_file[[2]], unet)
    fdata <- data.frame(matrix(ncol=8, nrow=0))
    qdata <- data.frame(matrix(ncol=3, nrow=0))

    wb1 <- openxlsx::createWorkbook(paste0(save_dir,'/Supplementary_Table_S2_',net_type[qq],'.xlsx'))
    mdn <- c()

    for(k in 1:length(cancer_type)){
        c_type <- cancer_type[k]
        indir <- paste0(in_dir,'/threshold_',cpm_threshold)
        ## ----- CRPEs --------
        gained_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Gained_Surv_.txt'))
        gained_surv1 <- gained_surv[gained_surv$pval <= pval_thres, ]
        gs <- igraph::graph_from_data_frame(gained_surv1[,c(1,2)], directed=FALSE)
        lost_surv <- data.table::fread(paste0('../data/Final_survival_filt/',strsplit(basename(allnets[qq]),'[.]')[[1]][1],'/threshold_',cpm_threshold,'/',
        c_type,'_Lost_Surv_.txt'))
        lost_surv1 <- lost_surv[lost_surv$pval <= pval_thres, ]
        ls <- igraph::graph_from_data_frame(lost_surv1[,c(1,2)], directed=FALSE)
        pert <- igraph::union(gs, ls)
        allpert <- igraph::as_data_frame(pert)
        fpert <- mapProtein(allpert[[1]], allpert[[2]], unet)
        ### unique PPIs 
        tempgp <- igraph::graph_from_data_frame(fpert[,c(3,4)], directed=FALSE)
        tempgp <- igraph::simplify(tempgp)
        tempgp <- igraph::as_data_frame(tempgp)
        tempdata <- data.frame(matrix(ncol=5, nrow=0))
        counterg <- 0
        countg <- c()
        countgp <- c()
        patg <- c()
        ppic1 <- c()
        ppic2 <- c()
        tot_eeis <- c()
        pert_eeis <- c()
        patientg <- c()
        crpes <- c()
        temp_mdn <- c()
        num_patients <- gained_surv[[7]][1]/gained_surv[[6]][1]
        save_ex <- data.frame(matrix(nrow=0, ncol=5))
        ###--- both gained and lost together
        tempgl <- tempgp

        ##--- for each PPI select the EEIs 
        for(j in 1:length(tempgl[[1]])){

            wh1 <- which(net_file$protein1 == tempgl[[1]][j])
            wh2 <- which(net_file$protein2 == tempgl[[2]][j])
            wha <- intersect(wh1, wh2)
            wh1 <- which(net_file$protein1 == tempgl[[2]][j])
            wh2 <- which(net_file$protein2 == tempgl[[1]][j])
            whb <- intersect(wh1, wh2)
            wh <- union(wha, whb)
            tempxg <- net_file[wh, ] ## Among all EEIs for this PPI which of them are not perturbed 

            ## get crpe
            wh1 <- which(fpert$protein1 == tempgl[[1]][j])
            wh2 <- which(fpert$protein2 == tempgl[[2]][j])
            wha <- intersect(wh1, wh2)
            wh1 <- which(fpert$protein1 == tempgl[[2]][j])
            wh2 <- which(fpert$protein2 == tempgl[[1]][j])
            whb <- intersect(wh1, wh2)
            wh <- union(wha, whb)
            temp_crpe <- fpert[wh, ] ## Among all EEIs for this PPI which of them are not perturbed 
            tempcr_gr <- igraph::graph_from_data_frame(temp_crpe[,c(1,2)], directed=FALSE) ## make a graph of all eeis perturbed in at least one patient

            
            tempa <- getprop(tempxg[[1]], tempxg[[2]], gained_surv[[1]], gained_surv[[2]], gained_surv[[7]])
            tempb <- getprop(tempxg[[1]], tempxg[[2]], lost_surv[[1]], lost_surv[[2]], lost_surv[[7]])

            tempxg$patient <- tempa[[1]]+tempb[[1]] ## perturbed

            ##-- where are the crpes--
            crpe_pos <- c()
            for(i in 1:length(temp_crpe[[1]])){
                wh1 <- which(tempxg$exon1 == temp_crpe[[1]][i])
                wh2 <- which(tempxg$exon2 == temp_crpe[[2]][i])
                wha <- intersect(wh1, wh2)
                wh1 <- which(tempxg$exon2 == temp_crpe[[1]][i])
                wh2 <- which(tempxg$exon1 == temp_crpe[[2]][i])
                whb <- intersect(wh1, wh2)
                wh <- union(wha, whb)
                crpe_pos <- c(crpe_pos, wh)
            }
            

            # tempxg_gr <- igraph::graph_from_data_frame(tempxg[tempxg$patient != 0][,c(1,2)], directed=FALSE) ## make a graph of all eeis perturbed in at least one patient
            # diff1 <- igraph::difference(tempcr_gr, tempxg_gr)
            # diff2 <- igraph::difference(tempxg_gr,tempcr_gr)

            tempyg <- tempxg[tempxg$patient == 0, ]

            if((nrow(tempyg) != 0)){# & igraph::ecount(diff1) == 0 & igraph::ecount(diff2) == 0){ ## meaning that there are some EEIs that are not perturbed in any patient

                counterg  <- counterg+1 ## counting the number of PPIs with at least one non-perturbed EEI
                countg <- c(countg, (length(tempyg[[1]])/length(tempxg[[1]]))*100) ## % of non-perturbed
                tempzg <- tempxg[tempxg$patient != 0, ]
                countgp <- c(countgp, (length(tempzg[[1]])/length(tempxg[[1]]))*100) ## % of perturbed
                patg <- c(patg, mean(tempzg$patient))#/num_patients)*100) ## mean number of patients in which the edges are perturbed
                ppic1 <- c(ppic1, tempzg$protein1[1])
                ppic2 <- c(ppic2, tempzg$protein2[1])
                tot_eeis <- c(tot_eeis, length(tempxg[[1]]))
                pert_eeis <- c(pert_eeis, length(tempzg[[1]]))
                # crpes <- crpes+length(tempag[[1]])
                tempag <- tempxg[crpe_pos,]
                save_ex <- rbind(save_ex, save_impPPIs(tempxg))
                patientg <- c(patientg,(tempag$patient/num_patients)*100)
                temp_mdn <- c(temp_mdn, tempag$patient)
            }
            tempxg$flag <- rep('Perturbed', length(tempxg[[1]]))
            tempdata <- rbind(tempdata, tempxg)
        }

        ##-- store the fractions ---
        if(length(ppic1) != 0){

        tdatag <- data.frame(Protein1=ppic1, Protein2=ppic2, Perturbed=countgp, 
            Number_of_PPIs=rep(length(tempgp[[1]]), length(countgp)), 
            cancer=cancer_type[k], Patients=patg, Non_perturbed=countg, 
            Number_of_desired_PPIs=rep(counterg, length(countgp)), Total_EEIs=tot_eeis, Perturbed_EEIs=pert_eeis)
        fdata <- rbind(fdata, tdatag)

        tdatagg <- data.frame(Cancer=rep(cancer_type[k], length(patientg)), Patients=patientg, 
            Number_of_PPIs=rep(counterg, length(patientg)), Number_of_EEIs=rep(length(patientg), length(patientg)))
        qdata <- rbind(qdata, tdatagg)

        mdn <- c(mdn, median(temp_mdn))

        }

        ## save excel sheet ----
        openxlsx::addWorksheet(wb1, sheetName = cancer_type[k])
        openxlsx::writeData(wb1, sheet = cancer_type[k], save_ex)
        openxlsx::saveWorkbook(wb1, paste0(save_dir,'/Supplementary_Table_S2_',net_type[qq],'.xlsx'), overwrite = T)
    }


    ufl <- 400
    if(qq==2){
        ufl <- 250
    }else if(qq==3){
        ufl <- 20
    }


    ##--- number of interesting EEIs ------
    pdatax <- unique(qdata[,c(1,4)])
    # p <- ggplot(pdata, aes(Cancer, count, fill=Cancer)) + 
    p <- ggplot(pdatax, aes(Cancer, Number_of_EEIs)) + 
    geom_bar(stat="identity",position=position_dodge())+
    theme(legend.text=element_text(size=12))
    basesize <- 12
    p <- p + theme_bw(base_size = basesize * 0.8) +
    scale_x_discrete(name="Cancer type") + 
    scale_y_continuous(name="# of AS-selective cancer \nrelevant perturbed egdes (CRPEs)", limits=c(0,(max(pdatax$Number_of_EEIs))+ufl)) +
    geom_text(aes(label=Number_of_EEIs), position=position_dodge(width=0.9),hjust=0, vjust=0, angle=75, size=3)+
    # scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
    axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
    strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+
    guides(fill='none')#guide_legend(title="Cancer type",ncol=2))
    ggsave(p,filename=paste0(save_dir,"/",net_type[qq],"_intEEIs.png"),width=3.5, height=3, dpi=400)



    ##--- boxplot showing the fraction of PPIs in which at least one EEI is not perturbed 
    td <- qdata[,c(1,2)]
    mxp <-  td %>% group_by(Cancer) %>% summarise(mv = max(Patients))
    # mdn <-  td %>% group_by(Cancer) %>% summarise(mv = median(Patients))

    pmdn <- data.frame(Cancer=cancer_type,count=mdn,pos=mxp$mv+10)
    ppx <- ggplot(data = qdata, aes(y=Patients, x=Cancer, fill=Cancer, group=Cancer)) + geom_boxplot(alpha=0.9)+
    # geom_point(aes(fill = Cancer), shape = 21, position = position_jitterdodge()) +
    xlab("Cancer type")+
    ylab("% of patients \nwith paired samples")+
    scale_x_discrete()+
    scale_y_continuous(limits=c(0,100))+
    scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#ffff99','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#e31a1c','#b15928','black','#9e0142','#053061'))+
    geom_text(data=pmdn, aes(y=pos, x=Cancer,label=count), position=position_dodge(width=0.9),hjust=0.5, vjust=0.5, angle=60, size=3)+
    # geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.01, size=0.1) +
    theme_bw()+theme(axis.text.x = element_text(size = 8, angle = 60, vjust=1, hjust=1, colour = "black"),
        axis.text.y = element_text(size = 8, angle = 0, colour = "black"))+#,
        # panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill='none')#+guides(size=guide_legend(title="# of patients",ncol=1))

    # guides(color=guide_legend(title="Cancer type", ncol=2,bycol=TRUE), shape=guide_legend(title='Edge type'))
    ggsave(ppx,filename=paste0(save_dir,"/AS-selective_CRPEs_",net_type[qq],".png"),width=4.5, height=2.5, dpi=400)


}
