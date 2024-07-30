##############################################################################################--
# Purpose: Distribution of numbers of control, condition, perturbed, and survival edges with diffreent cpm thresholds
##############################################################################################--

rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(GenomicDataCommons)
source("eein_cancer_util.r")
options(warn=0)

save_dir <- '../results_rep'
dir.create(save_dir)

cancer_type <- gtools::mixedsort(c('BLCA', 'BRCA', 'KIRC', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'UCEC', 'THCA', 'COAD', 'PRAD', 'KICH', 'STAD', 'ESCA'))
cpm_threshold <- c(0.05, 0.1, 0.2, 0.5, 1, 2)
pval_thres <- 0.05

get_num_edges <- function(tempf_control, tempf_condition, exons, net_file){

	 control_graph <- list()
    condition_graph <- list()
    gained_graph <- list()
    lost_graph <- list()
    perturbed_edges <-  list()

    gain_net <- net_file
    loss_net <- net_file

    control_names <- unlist(lapply(strsplit(colnames(tempf_control), '[_]'), '[[', 1)) 
    condition_names <- unlist(lapply(strsplit(colnames(tempf_condition), '[_]'), '[[', 1)) 

    for(k in 1:length(tempf_control)){

        control_file <- tempf_control[,k]
        control_wh <- which(control_file != 0)
        control_nodes <- exons[control_wh]
        wh1 <- which(net_file$exon1 %in% control_nodes)
        wh2 <- which(net_file$exon2 %in% control_nodes)
        wh <- intersect(wh1, wh2)
        control_net <- net_file[wh, ]
        control_ig <- igraph::graph_from_data_frame(control_net, directed=FALSE)
        control_graph[[k]] <- igraph::ecount(control_ig)

        whc <- which(condition_names == control_names[k])
        condition_file <- tempf_condition[,whc]
        condition_wh <- which(condition_file != 0)
        condition_nodes <- exons[condition_wh]
        wh1 <- which(net_file$exon1 %in% condition_nodes)
        wh2 <- which(net_file$exon2 %in% condition_nodes)
        wh <- intersect(wh1, wh2)
        condition_net <- net_file[wh, ]
        condition_ig <- igraph::graph_from_data_frame(condition_net, directed=FALSE)
        condition_graph[[k]] <- igraph::ecount(condition_ig)

        lle <- igraph::difference(control_ig, condition_ig, byname=TRUE)
        lost_graph[[k]] <- lle
        if(igraph::ecount(lle) != 0){
            lled <- as.data.frame(igraph::get.edgelist(lle))
            colnames(lled) <- c('exon1', 'exon2')
            loss_net <- rbind(loss_net, lled)
        }

        gge <- igraph::difference(condition_ig, control_ig, byname=TRUE)
        gained_graph[[k]] <- gge
        if(igraph::ecount(lle) != 0){
            gged <- as.data.frame(igraph::get.edgelist(gge))
            colnames(gged) <- c('exon1', 'exon2')
            gain_net <- rbind(gain_net, gged)
        }

        perturbed_edges[[k]] <- igraph::ecount(igraph::union(lle, gge))
    }

    return(list(control_graph, condition_graph, perturbed_edges))
}

get_pert_edges <- function(tempf_control, tempf_condition, exons, net_file, allsurv){

	 control_graph <- list()
	 condition_graph <- list()
	 gained_graph <- list()
	 lost_graph <- list()
	 perturbed_graph <-  list()
	 cancer_relevant_edges <-  list()
	 cancer_relevant_graph <-  list()

	 gain_net <- net_file
	 loss_net <- net_file

	 control_names <- unlist(lapply(strsplit(colnames(tempf_control), '[_]'), '[[', 1)) 
	 condition_names <- unlist(lapply(strsplit(colnames(tempf_condition), '[_]'), '[[', 1)) 

	 for(k in 1:length(tempf_control)){

	     control_file <- tempf_control[,k]
	     control_wh <- which(control_file != 0)
	     control_nodes <- exons[control_wh]
	     wh1 <- which(net_file$exon1 %in% control_nodes)
	     wh2 <- which(net_file$exon2 %in% control_nodes)
	     wh <- intersect(wh1, wh2)
	     control_net <- net_file[wh, ]
	     control_ig <- igraph::graph_from_data_frame(control_net, directed=FALSE)
	     control_graph[[k]] <- igraph::ecount(control_ig)

	     whc <- which(condition_names == control_names[k])
	     condition_file <- tempf_condition[,whc]
	     condition_wh <- which(condition_file != 0)
	     condition_nodes <- exons[condition_wh]
	     wh1 <- which(net_file$exon1 %in% condition_nodes)
	     wh2 <- which(net_file$exon2 %in% condition_nodes)
	     wh <- intersect(wh1, wh2)
	     condition_net <- net_file[wh, ]
	     condition_ig <- igraph::graph_from_data_frame(condition_net, directed=FALSE)
	     condition_graph[[k]] <- igraph::ecount(condition_ig)

	     # if(igraph::ecount(condition_ig) == 0){
	     # 	lle <- control_ig
	     # 	lost_graph[[k]] <- lle
		 #     if(igraph::ecount(lle) != 0){
		 #         lled <- as.data.frame(igraph::get.edgelist(lle))
		 #         colnames(lled) <- c('exon1', 'exon2')
		 #         loss_net <- rbind(loss_net, lled)
		 #     }
	     # }else if(igraph::ecount(control_ig) != 0 & igraph::ecount(condition_ig) != 0){
		     lle <- igraph::difference(control_ig, condition_ig, byname=TRUE)
		     lost_graph[[k]] <- lle
		     if(igraph::ecount(lle) != 0){
		         lled <- as.data.frame(igraph::get.edgelist(lle))
		         colnames(lled) <- c('exon1', 'exon2')
		         loss_net <- rbind(loss_net, lled)
		     }
	     # }


	     gge <- igraph::difference(condition_ig, control_ig, byname=TRUE)
	     gained_graph[[k]] <- gge
	     if(igraph::ecount(lle) != 0){
	         gged <- as.data.frame(igraph::get.edgelist(gge))
	         colnames(gged) <- c('exon1', 'exon2')
	         gain_net <- rbind(gain_net, gged)
	     }

	     pert_graph <- igraph::union(lle, gge)
	     perturbed_graph[[k]] <- pert_graph
		  cancer_relevant_edges[[k]] <- igraph::ecount(igraph::intersection(pert_graph, allsurv))
		  cancer_relevant_graph[[k]] <- igraph::intersection(pert_graph, allsurv)

	 }

	 return(list(control_graph, condition_graph, cancer_relevant_edges, perturbed_graph, cancer_relevant_graph))
}


##-- read all expr data --
all_expr_data <- list()
for(j in 1:length(cancer_type)){
	all_expr_data[[j]] <- data.table::fread(paste0('../data/filtered_normalized_data/', cancer_type[j],'/normalized_all.txt'), header=TRUE)
}

##-- read all manifest data -- 
all_manifest <- list()
for(j in 1:length(cancer_type)){
	all_manifest[[j]] <- data.table::fread(paste0('../data/Manifests/',cancer_type[j],'_manifest_final.txt'))
}

##--- create networks of different confidence ----
contact_net <- data.table::fread(paste0('../data/final_EEINs/CONTACT.txt'))
contact_net <- igraph::graph_from_data_frame(contact_net[,c(1,2)], directed=FALSE)
pisa_net <- data.table::fread(paste0('../data/final_EEINs/PISA.txt'))
pisa_net <- igraph::graph_from_data_frame(pisa_net[,c(1,2)], directed=FALSE)
eppic_net <- data.table::fread('../data/final_EEINs/EPPIC.txt')
eppic_net <- igraph::graph_from_data_frame(eppic_net[,c(1,2)], directed=FALSE)
net_file1 <- igraph::as_data_frame(igraph::union(igraph::union(contact_net, pisa_net), eppic_net))

net_file2 <- igraph::union(igraph::union(igraph::intersection(pisa_net,eppic_net,keep.all.vertices=FALSE),
		igraph::intersection(contact_net,eppic_net,keep.all.vertices=FALSE)), 
		igraph::intersection(contact_net,pisa_net,keep.all.vertices=FALSE))
net_file2 <- igraph::as_data_frame(net_file2)

net_file3 <- igraph::intersection(igraph::intersection(pisa_net,eppic_net,keep.all.vertices=FALSE),
	contact_net,keep.all.vertices=FALSE) 
net_file3 <- igraph::as_data_frame(net_file3)

allnets <- list(net_file1, net_file2, net_file3)
net_type <- c('NETLOW', 'NETMEDIUM', 'NETHIGH')
fig_num <- c('Supplementary_Fig_S4','Supplementary_Fig_S3','Supplementary_Fig_S2')

for(ntype in 1:length(allnets)){

	net_file <- allnets[[ntype]]
	colnames(net_file) <- c('exon1','exon2')
	nodes <- union(net_file[[1]], net_file[[2]])

	##-- for plot -----
	control_edges <- c()
	condition_edges <- c()
	perturbed_edges <- c()
	survival_edges <- c()
	pcancer <- c()
	thresholds <- c()
	##-----------------

	for(thr in 1:length(cpm_threshold)){

		for(cn in 1:length(cancer_type)){

			all_survival <- get_survival(cancer_type[cn], 0)
			tempids <- all_manifest[[cn]]
			tempf <- as.data.frame(all_expr_data[[cn]])
			tempf <- tempf[tempf$EXON %in% nodes, ]
			sample_names_all <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 1))
			whc <- which(sample_names_all %in% union(tempids$nid,'EXON'))
			tempf <- tempf[,whc]
			exons <- tempf$EXON ## this order is constant
			tempf <- tempf[,-length(tempf)]
			##-- only keep patients with control and condition ----------------------------
			temp_brk <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 2))
			sample_names <- unlist(lapply(strsplit(colnames(tempf), '[_]'), '[[', 1))
			condition_pos_all <- which( (temp_brk == '1') | (temp_brk == '01') )
			control_pos <- which(temp_brk == '11')
			keep_samples <- sample_names[control_pos]
			condition_pos <- intersect(which(sample_names %in% keep_samples), condition_pos_all)
			tempf_control <- tempf[,control_pos]
			tempf_condition <- tempf[,condition_pos]
			tempf_control[tempf_control < cpm_threshold[thr] ] <- 0 
			tempf_condition[tempf_condition < cpm_threshold[thr] ] <- 0 

			##----for real expression--------------------------------
			out_real <- get_num_edges(tempf_control, tempf_condition, exons, net_file)
			control_edges <- c(control_edges, unlist(out_real[[1]]))
			condition_edges <- c(condition_edges, unlist(out_real[[2]]))
			perturbed_edges <- c(perturbed_edges, unlist(out_real[[3]]))
			pcancer <- c(pcancer, rep(cancer_type[cn], length(keep_samples)))
			thresholds <- c(thresholds, rep(cpm_threshold[thr], length(keep_samples)))

			##--- cancer-relevant edges ---
			pisa_g <- data.table::fread(paste0('../data/PISA_survival_filt/threshold_',cpm_threshold[thr],'/',cancer_type[cn],'_Gained_Surv_.txt'))
			pisa_g <- pisa_g[pisa_g$pval <= pval_thres, ][,c(1,2)]
			pisa_l <- data.table::fread(paste0('../data/PISA_survival_filt/threshold_',cpm_threshold[thr],'/',cancer_type[cn],'_Lost_Surv_.txt'))
			pisa_l <- pisa_l[pisa_l$pval <= pval_thres, ][,c(1,2)]
			pisa <- igraph::graph_from_data_frame(rbind(pisa_g,pisa_l), directed=FALSE)

			eppic_g <- data.table::fread(paste0('../data/EPPIC_survival_filt/threshold_',cpm_threshold[thr],'/',cancer_type[cn],'_Gained_Surv_.txt'))
			eppic_g <- eppic_g[eppic_g$pval <= pval_thres, ][,c(1,2)]
			eppic_l <- data.table::fread(paste0('../data/EPPIC_survival_filt/threshold_',cpm_threshold[thr],'/',cancer_type[cn],'_Lost_Surv_.txt'))
			eppic_l <- eppic_l[eppic_l$pval <= pval_thres, ][,c(1,2)]
			eppic <- igraph::graph_from_data_frame(rbind(eppic_g,eppic_l), directed=FALSE)

			contact_g <- data.table::fread(paste0('../data/CONTACT_survival_filt/threshold_',cpm_threshold[thr],'/',cancer_type[cn],'_Gained_Surv_.txt'))
			contact_g <- contact_g[contact_g$pval <= pval_thres, ][,c(1,2)]
			contact_l <- data.table::fread(paste0('../data/CONTACT_survival_filt/threshold_',cpm_threshold[thr],'/',cancer_type[cn],'_Lost_Surv_.txt'))
			contact_l <- contact_l[contact_l$pval <= pval_thres, ][,c(1,2)]
			contact <- igraph::graph_from_data_frame(rbind(contact_g,contact_l), directed=FALSE)

			allsurv <- igraph::union(igraph::union(pisa, eppic), contact)
			out_surv_real <- get_pert_edges(tempf_control, tempf_condition, exons, net_file, allsurv)
			survival_edges <- c(survival_edges, unlist(out_surv_real[[3]]))

		}
	}

	pdata <- data.frame(Cancer=pcancer, Control=control_edges, Condition=condition_edges, 
		Perturbed=perturbed_edges, Survival=survival_edges, Threshold=thresholds)

	pdata1 <- pdata %>%
	   group_by(Cancer, Threshold) %>% 
	   summarise_at(vars("Control", "Condition", "Perturbed", "Survival"), mean)

	colnames(pdata1) <- c('Cancer','Threshold',"Edges in control samples", "Edges in condition samples", "Perturbed edges", "Cancer-relevant perturbed edges")
	pdata2 <- reshape2::melt(as.data.frame(pdata1), id=c('Cancer','Threshold'))

	## get max ----------
	maxv <- c()
	maxvp <- c()

	for(k in 1:length(cancer_type)){
		tempq <- pdata2[pdata2$Cancer == cancer_type[k], ]
		temp3 <- tempq[tempq$variable == 'Cancer-relevant perturbed edges', ]
		maxv <- c(maxv, temp3$Threshold[which(temp3$value == max(temp3$value))[1]])
		temp3 <- tempq[tempq$variable == 'Perturbed edges', ]
		maxvp <- c(maxvp, temp3$Threshold[which(temp3$value == max(temp3$value))[1]])
	}
	thres_count_surv <- plyr::count(maxv)
	thres_count_surv$flag <- rep('Cancer-relevant perturbed edges', length(thres_count_surv[[1]]))

	thres_count_pert <- plyr::count(maxvp)
	thres_count_pert$flag <- rep('Perturbed edges', length(thres_count_pert[[1]]))

	thres_count <- rbind(thres_count_pert, thres_count_surv)

	# p <- ggplot(pdata2, aes(Threshold, value, color=variable, group=variable)) + geom_point()+geom_line()+
	# theme(legend.text=element_text(size=12))
	# basesize <- 15
	# p <- p + theme_grey(base_size = basesize * 0.8) + 
	# scale_x_continuous(name="Exon expression threshold (normalized CPM value)",breaks = seq(0, 2, by = 0.25)) + 
	# scale_y_continuous(name="# of edges",breaks = seq(0, length(net_file[[1]]), by = 2000)) +
	# scale_color_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3'))+
	# guides(fill=guide_legend(title="Threshold"))+
	# facet_wrap(~Cancer, ncol=3, scale='free')+
	# theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
	# axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
	# strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+ 
	# guides(color=guide_legend(title="Edge type"))
	# ggsave(p,filename=paste0(save_dir,'/',net_type[ntype],'_EEI_CPM_variation.png'),width=14, height=12, dpi=400)


	pdatax <- pdata[,-1]
	pdata1x <- pdatax %>%
	   group_by(Threshold) %>% 
	   summarise_at(vars("Control", "Condition", "Perturbed", "Survival"), mean)

	pdata1xsd <- pdatax %>%
	   group_by(Threshold) %>% 
	   summarise_at(vars("Control", "Condition", "Perturbed", "Survival"), sd)


	colnames(pdata1x) <- c('Threshold',"Edges in control samples", "Edges in condition samples", "Perturbed edges", "Cancer-relevant perturbed edges")
	colnames(pdata1xsd) <- c('Threshold',"Edges in control samples", "Edges in condition samples", "Perturbed edges", "Cancer-relevant perturbed edges")
	  

	pdata2x <- reshape2::melt(as.data.frame(pdata1x), id=c('Threshold'))
	pdata2xsd <- reshape2::melt(as.data.frame(pdata1xsd), id=c('Threshold'))
	pdata2x$err <- pdata2xsd$value

	p <- ggplot(pdata2x, aes(Threshold, value, color=variable, group=variable)) + geom_point()+geom_line()+
	theme(legend.text=element_text(size=12))
	basesize <- 10
	p <- p + theme_bw(base_size = basesize * 0.8) + 
	scale_x_continuous(name="Exon expression threshold \n(normalized CPM value)",breaks = seq(0, 2, by = 0.25)) + 
	scale_y_continuous(name="# of edges",breaks = seq(0, length(net_file[[1]]), by = 2000)) +
	scale_color_manual(values=c('#e41a1c','#377eb8','#4daf4a','#984ea3'))+
	geom_errorbar(aes(ymin=value-err, ymax=value+err), width=0.05 )+
	guides(fill=guide_legend(title="Threshold"))+
	theme(axis.text.x = element_text(size = basesize * 0.6, angle = 60, hjust = 0.5,vjust=0.5, colour = "black"),
	axis.text.y = element_text(size = basesize * 0.6, angle = 0, hjust = 0.5,vjust=0.5, colour = "black"), 
	strip.text = element_text(size = basesize * 0.8), axis.title=element_text(basesize * 0.8))+ 
	guides(color=guide_legend(title="Edge type"))
	ggsave(p,filename=paste0(save_dir,'/',fig_num[ntype],'.png'),width=7, height=3.5, dpi=400)


}

