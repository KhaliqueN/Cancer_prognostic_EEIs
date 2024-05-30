
get_ppis <- function(tempx, nfile){
    ##--get PPIs --
    why <- c()
    for(b in 1:length(tempx[[1]])){
        wh1 <- which(nfile[[1]] == tempx[[1]][b])
        wh2 <- which(nfile[[2]] == tempx[[2]][b])
        wha <- intersect(wh1, wh2)
        wh1 <- which(nfile[[1]] == tempx[[2]][b])
        wh2 <- which(nfile[[2]] == tempx[[1]][b])
        whb <- intersect(wh1, wh2)
        wh <- union(wha, whb)
        why <- c(why, wh[1])
    }
    tempy <- nfile[why,]
    tempy <- igraph::simplify(igraph::graph_from_data_frame(tempy[,c(5,6)], directed=FALSE), remove.multiple=TRUE, remove.loop=TRUE)
    return(tempy)
}

## CALLED IN: 
# 9a_map_EEIs.r
select_edges <- function(int_tab_filt, relevant_chains, input_dir, maps){ 
    tEXON1 <- c()
    tEXON2 <- c()
    tPROTEIN1 <- c()
    tPROTEIN2 <- c()
    tINTFAA1 <- c()
    tINTFAA2 <- c()
    tFREE_ENERGY <- c()
    tBURIED_AREA <- c()
    tHYDROGEN <- c()
    tDISULPHIDE <- c()
    tSALTBRIDGE <- c()
    tCOVALENT <- c()
    tBURIED_AREA_ABS <- c()
    tASA_AREA_ABS <- c()
    tCHAIN1 <- c()
    tCHAIN2 <- c()

    wh1 <- which(int_tab_filt$X2 %in% relevant_chains)
    wh2 <- which(int_tab_filt$X5 %in% relevant_chains)
    wh <- intersect(wh1, wh2)

    temp_int_tab <- int_tab_filt[wh,]
    int_tab_filt <- unique(temp_int_tab[,c(2,5)])
    res_files <- gtools::mixedsort(list.files(input_dir, pattern='^residue', full.names=TRUE))

    if(nrow(int_tab_filt) != 0){
        for(i in 1:length(int_tab_filt[[1]])){

            chainA <- int_tab_filt[[1]][i]
            chainB <- int_tab_filt[[2]][i]
            ochains <- c(chainA, chainB)

            ## temp int file
            temp_int_tabx <- temp_int_tab[(temp_int_tab$X2 == chainA) & (temp_int_tab$X5 == chainB), ]
            temp_int_tabx$res1 <- trimws(substr(unlist(lapply(strsplit(temp_int_tabx[[3]],'[[]'),'[[',1)),5,50))
            temp_int_tabx$res2 <- trimws(substr(unlist(lapply(strsplit(temp_int_tabx[[6]],'[[]'),'[[',1)),5,50))

            for(j in 1:length(res_files)){
                
                res_tab <- data.table::fread(res_files[j])
                res_tab <- res_tab[res_tab$Solvation_energy != 0, ] ##-- get the interface residues --> with Solvation energy != 0

                allchains <- unlist(lapply(strsplit(res_tab[[2]],'[:]'), '[[', 1))
                ichains <- c(intersect(allchains, ochains[1]), intersect(allchains, ochains[2]))
                # other_chains <- setdiff(allchains, ochains)

                # if((length(ichains) == 2) & (length(other_chains) == 0)){ ## if both the interacting chains are present
                # if(length(ichains) == 2){ ## if both the interacting chains are present
                if((length(ichains) == 2) & (ichains[1] != ichains[2])){ ## if both the interacting chains are present and the interacting chnains are different

                    wh <- which(allchains %in% ichains)
                    res_tab1 <- res_tab[wh, ]
                    allchains <- unlist(lapply(strsplit(res_tab1[[2]],'[:]'), '[[', 1))

                    ##-- Get the exons for each of the two chains
                    ##-- CHAIN 1
                    wh1 <- which(relevant_chains == ichains[1])
                    res1 <- trimws(substr(unlist(lapply(strsplit(res_tab1[which(allchains == ichains[1]), ][[2]],'[:]'),'[[',2)),4,50))

                    ##-- could be multiple exon maps --> corresponding to different proteins --> some parts of a pdb file maps to one protein while others to the other
                    ex1 <- c()
                    intf1 <- c()
                    sol_e1 <- c()
                    asa_e1 <- c()
                    bsa_e1 <- c()
                    res_e1 <- list()
                    counter <- 1
                    unipr1 <- c()
                    for(vv in 1:length(wh1)){
                        tempm <- data.table::fread(maps[wh1[vv]])
                        tp1 <- strsplit(basename(maps[wh1[vv]]),'[_]')[[1]][1]
                        whr <- which(tempm$PDBResNumAuthor %in% res1)
                        tempm1 <- tempm[whr, ]
                        exons1 <- unique(tempm1$EXON)

                        if(!identical(exons1, character(0))){
                            for(m in 1:length(exons1)){
                                ex1 <- c(ex1, exons1[m])
                                unipr1 <- c(unipr1, tp1)
                                ## interfacing residues count as per the Uniprot PDB mapping...
                                intf1 <- c(intf1, length(tempm1[tempm1$EXON == exons1[m], ][[1]]))
                                ## take the sovation energy and solvent accessible area
                                ex1_res <- tempm1[tempm1$EXON == exons1[m], ]$PDBResNumAuthor
                                res_e1[[counter]] <- ex1_res
                                res_tab2 <- res_tab1[which(res1 %in% ex1_res), ]
                                sol_e1 <- c(sol_e1, sum(res_tab2$Solvation_energy))
                                asa_e1 <- c(asa_e1, sum(res_tab2$Solvent_accessible_area))
                                bsa_e1 <- c(bsa_e1, sum(res_tab2$Buried_area))
                                counter <- counter+1
                            }
                        }
                    }
                    
                    ##-- CHAIN 2
                    wh2 <- which(relevant_chains == ichains[2])
                    res2 <- trimws(substr(unlist(lapply(strsplit(res_tab1[which(allchains == ichains[2]), ][[2]],'[:]'),'[[',2)),4,50))

                    ##-- could be multiple maps --> for different proteins
                    ex2 <- c()
                    intf2 <- c()
                    sol_e2 <- c()
                    asa_e2 <- c()
                    bsa_e2 <- c()
                    res_e2 <- list()
                    counter <- 1
                    unipr2 <- c()
                    for(vv in 1:length(wh2)){
                        tempm <- data.table::fread(maps[wh2[vv]])
                        tp2 <- strsplit(basename(maps[wh2[vv]]),'[_]')[[1]][1]
                        whr <- which(tempm$PDBResNumAuthor %in% res2)
                        tempm1 <- tempm[whr, ]
                        exons2 <- unique(tempm1$EXON)

                        if(!identical(exons2, character(0))){
                            for(m in 1:length(exons2)){
                                ex2 <- c(ex2, exons2[m])
                                unipr2 <- c(unipr2, tp2)
                                intf2 <- c(intf2, length(tempm1[tempm1$EXON == exons2[m], ][[1]]))
                                ## take the sovation energy and solvent accessible area
                                ex2_res <- tempm1[tempm1$EXON == exons2[m], ]$PDBResNumAuthor
                                res_e2[[counter]] <- ex2_res
                                res_tab2 <- res_tab1[which(res2 %in% ex2_res), ]
                                sol_e2 <- c(sol_e2, sum(res_tab2$Solvation_energy))
                                asa_e2 <- c(asa_e2, sum(res_tab2$Solvent_accessible_area))
                                bsa_e2 <- c(bsa_e2, sum(res_tab2$Buried_area))
                                counter <- counter+1
                            }
                        }
                    }
                    
                    ## save the exon-exon interactions
                    if((length(ex1)!=0) & (length(ex2)!= 0)){
                        for(mm in 1:length(ex1)){
                            for(nn in 1:length(ex2)){
                                tPROTEIN1 <- c(tPROTEIN1, unipr1[mm])
                                tPROTEIN2 <- c(tPROTEIN2, unipr2[nn])
                                tCHAIN1 <- c(tCHAIN1, ichains[1])
                                tCHAIN2 <- c(tCHAIN2, ichains[2])
                                tEXON1 <- c(tEXON1, ex1[mm])
                                tEXON2 <- c(tEXON2, ex2[nn])
                                tINTFAA1 <- c(tINTFAA1, intf1[mm])
                                tINTFAA2 <- c(tINTFAA2, intf2[nn])
                                tFREE_ENERGY <- c(tFREE_ENERGY, -1*(sol_e1[mm]+sol_e2[nn]))
                                tBURIED_AREA <- c(tBURIED_AREA, (bsa_e1[mm]+bsa_e2[nn])/(asa_e1[mm]+asa_e2[nn]))
                                tBURIED_AREA_ABS <- c(tBURIED_AREA_ABS, (bsa_e1[mm]+bsa_e2[nn]))
                                tASA_AREA_ABS <- c(tASA_AREA_ABS, (asa_e1[mm]+asa_e2[nn]))
                                wha <- which(temp_int_tabx$res1 %in% res_e1[[mm]])
                                whb <- which(temp_int_tabx$res2 %in% res_e2[[nn]])
                                wh <- intersect(wha, whb)
                                temp_int_taby <- temp_int_tabx[wh, ]
                                if(nrow(temp_int_taby) != 0){
                                    tHYDROGEN <- c(tHYDROGEN, length(which(temp_int_taby$X7 %like% 'Hydrogen')))
                                    tSALTBRIDGE <- c(tSALTBRIDGE, length(which(temp_int_taby$X7 %like% 'Salt')))
                                    tDISULPHIDE <- c(tDISULPHIDE, length(which(temp_int_taby$X7 %like% 'Disulphide')))
                                    tCOVALENT <- c(tCOVALENT, length(which(temp_int_taby$X7 %like% 'Covalent')))
                                }else{
                                    tHYDROGEN <- c(tHYDROGEN, 0)
                                    tSALTBRIDGE <- c(tSALTBRIDGE, 0)
                                    tDISULPHIDE <- c(tDISULPHIDE, 0)
                                    tCOVALENT <- c(tCOVALENT, 0)
                                }
                            }
                        }
                    }
                    
                    break ## if the match found, then no need to compare other residue files
                }
            }
        }
    }

    return(list(tPROTEIN1, tPROTEIN2, tEXON1, tEXON2, tINTFAA1, tINTFAA2, tFREE_ENERGY, tBURIED_AREA, tHYDROGEN,
        tDISULPHIDE, tSALTBRIDGE, tCOVALENT, tBURIED_AREA_ABS, tASA_AREA_ABS, tCHAIN1, tCHAIN2)) 
}

## CALLED IN: 
# 1_weighted_network_contact.r
# 2_select_edges_contact.r
# 1_weighted_network_pisa.r
# 2_select_edges_pisa.r
# 1_weighted_network_EPPIC.r
# 2_select_edges_EPPIC.r
# 2_analysis_individual_final.r
# 1_cpm_variation.r
# 2b_CRPE_cancer_type.r
get_survival <- function(c_type, survp){ 
    tempids <- data.table::fread(paste0('../data/Manifests/',c_type,'_manifest_final.txt'))
    tempids_oc <- data.table::fread(paste0('../data/Manifests/',c_type,'_manifest_final_allCancer.txt'))
    tempids_srv <- rbind(tempids, tempids_oc) ## for which patients to save the survival info

    ## get clinical data
    case_ids = cases() %>%
        filter(~ project.project_id == paste0('TCGA-',c_type)) %>%
        ids()
    clindat = gdc_clinical(case_ids)
    all_sub_cases <- tempids_srv$nid
    diag_data1 <- data.frame(clindat$diagnoses)
    diag_data2 <- data.frame(clindat$demographic)
    diag_data1$submitter_id <- unlist(lapply(strsplit(diag_data1$submitter_id, '[_]'), '[[', 1))
    diag_data2$submitter_id <- unlist(lapply(strsplit(diag_data2$submitter_id, '[_]'), '[[', 1))

    diag_data_final1 <- diag_data1[diag_data1$submitter_id %in% all_sub_cases, ]
    diag_data_final1 <- diag_data_final1[c('submitter_id','days_to_last_follow_up')]

    diag_data_final2 <- diag_data2[diag_data2$submitter_id %in% all_sub_cases, ]
    diag_data_final2 <- diag_data_final2[c('submitter_id','days_to_death','gender','vital_status')]
    diag_data_final3 <- merge(diag_data_final1, diag_data_final2, by='submitter_id')
    diag_data_final3$deceased <- diag_data_final3$vital_status == 'Dead'

    ovrs <- c()
    for(jj in 1:length(diag_data_final3[[1]])){
        if(is.na(diag_data_final3$days_to_death[jj])){
            ovrs <- c(ovrs, diag_data_final3$days_to_last_follow_up[jj])
        }else{
            ovrs <- c(ovrs, diag_data_final3$days_to_death[jj])
        }
    }
    diag_data_final3$overall_survival <- ovrs
    survival_data <- diag_data_final3
    ##--filter out patients with no overall_survival data
    survival_data <- survival_data[!is.na(survival_data$overall_survival), ]
    ##-- take all deceased patients
    dsdata <- survival_data[which(survival_data$vital_status == 'Dead'), ]
    ##-- take all non-deceased patients
    ndsdata <- survival_data[which(survival_data$vital_status == 'Alive'), ]
    ##--filter out those patients that are non-deceased but less than 5 years
    ndsdata1 <- ndsdata[ndsdata$overall_survival >= survp, ]
    ##-- all data
    surv_data <- rbind(dsdata, ndsdata1)
    return(surv_data) 
}

## CALLED IN: 
# 1_weighted_network_contact.r
# 2_select_edges_contact.r
# 1_weighted_network_pisa.r
# 2_select_edges_pisa.r
# 1_weighted_network_EPPIC.r
# 2_select_edges_EPPIC.r
# 2_analysis_individual_final.r
get_perturbed <- function(tempfd_control, tempfd_condition, exons, net_file){ 

    gain_net <- net_file
    loss_net <- net_file

    control_names <- unlist(lapply(strsplit(colnames(tempfd_control), '[_]'), '[[', 1)) 
    condition_names <- unlist(lapply(strsplit(colnames(tempfd_condition), '[_]'), '[[', 1)) 

    control_graph <- list()
    condition_graph <- list()

    for(k in 1:length(tempfd_control)){

        control_file <- tempfd_control[,k]
        control_wh <- which(control_file != 0)
        control_nodes <- exons[control_wh]
        wh1 <- which(net_file$exon1 %in% control_nodes)
        wh2 <- which(net_file$exon2 %in% control_nodes)
        wh <- intersect(wh1, wh2)
        control_net <- net_file[wh, ]
        control_ig <- igraph::graph_from_data_frame(control_net, directed=FALSE)
        control_graph[[k]] <- control_ig

        whc <- which(condition_names == control_names[k])
        condition_file <- tempfd_condition[,whc]
        condition_wh <- which(condition_file != 0)
        condition_nodes <- exons[condition_wh]
        wh1 <- which(net_file$exon1 %in% condition_nodes)
        wh2 <- which(net_file$exon2 %in% condition_nodes)
        wh <- intersect(wh1, wh2)
        condition_net <- net_file[wh, ]
        condition_ig <- igraph::graph_from_data_frame(condition_net, directed=FALSE)
        condition_graph[[k]] <- condition_ig

        lle <- igraph::difference(control_ig, condition_ig, byname=TRUE)
        if(igraph::ecount(lle) != 0){
            lled <- as.data.frame(igraph::get.edgelist(lle))
            colnames(lled) <- c('exon1', 'exon2')
            loss_net <- rbind(loss_net, lled)
        }

        gge <- igraph::difference(condition_ig, control_ig, byname=TRUE)
        if(igraph::ecount(lle) != 0){
            gged <- as.data.frame(igraph::get.edgelist(gge))
            colnames(gged) <- c('exon1', 'exon2')
            gain_net <- rbind(gain_net, gged)
        }

    }

    all_gain_net <- igraph::graph_from_data_frame(gain_net, directed=FALSE)
    igraph::E(all_gain_net)$weight <- 1
    all_gain_net1 <- igraph::simplify(all_gain_net, edge.attr.comb=list(weight="sum"))
    all_gained <- igraph::as_data_frame(all_gain_net1)
    all_gained <- all_gained[order(-all_gained$weight), ]
    all_gained$weight <- all_gained$weight-1 ##-- subtracting 1 to denote that edge was not present in any patient
    all_gained$patient <- all_gained$weight
    all_gained$weight <- all_gained$weight/length(control_names) ##-- fraction of patients in which the edge is gained/lost

    all_loss_net <- igraph::graph_from_data_frame(loss_net, directed=FALSE)
    igraph::E(all_loss_net)$weight <- 1
    all_loss_net1 <- igraph::simplify(all_loss_net, edge.attr.comb=list(weight="sum"))
    all_lost <- igraph::as_data_frame(all_loss_net1)
    all_lost <- all_lost[order(-all_lost$weight), ]
    all_lost$weight <- all_lost$weight-1 ##-- subtracting 1 to denote that edge was not present in any patient
    all_lost$patient <- all_lost$weight
    all_lost$weight <- all_lost$weight/length(control_names) ##-- fraction of patients in which the edge is gained/lost

    return(list(all_gained,all_lost))
} 

random_surv <- function(overall_survival, tempw, rnd_run){ ###---checked
    tts <- c()
    for(gg in 1:rnd_run){
        tempw$overall_survival <- sample(tempw$overall_survival) ## randomize overall survival
        sdf <- survdiff(Surv(overall_survival, deceased) ~ flag, data=tempw)
        tt <- signif(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
        tts <- c(tts, tt)
    }
    return(tts) 
}

## CALLED IN: 
# 2_select_edges_contact.r
# 2_select_edges_pisa.r
# 2_select_edges_EPPIC.r
# 2_analysis_individual_final.r
survival_qval_num <- function(gl_data, cthres, tempfd_condition, surv_data, exons, rnd_run){ 

        temp_net <- igraph::as_data_frame(gl_data)
        samplenames <- unlist(lapply(strsplit(colnames(tempfd_condition), '[_]'), '[[', 1)) 
        wh1 <- setdiff(samplenames, surv_data$submitter_id)
        wh2 <- setdiff(surv_data$submitter_id, samplenames)
        if(length(wh1) > 0){
            wh3 <- which(samplenames %in% wh1)
            samplenames <- samplenames[-wh3]
            tempfd_condition <- tempfd_condition[,-wh3]
        }
        if(length(wh2) > 0){
            wh3 <- which(surv_data$submitter_id %in% wh2)
            surv_data <- surv_data[-wh3, ]
        }

        colnames(tempfd_condition) <- samplenames
        tempx <- as.data.frame(t(tempfd_condition))
        tempx <- tempx[order(rownames(tempx)), ]
        colnames(tempx) <- exons

        temp5 <- exon2eei(temp_net, tempx, cthres)
        # vars <- colnames(temp5)
        surv_data <- surv_data[order(surv_data$submitter_id), c(1,6,7)]

        surv_data_all <- cbind(surv_data, temp5)

        sgg <- colnames(surv_data_all)
        sgg <- sgg[-c(1:3)]

        ##--- for each selected edge ----------------------------------------------------
        pvs <- c()
        s_edges <- c()

        ppvs <- list()
        counter <- 1

        if(length(sgg) != 0){
            for(i in 1:length(sgg)){
                tempe <- sgg[i]
                wh <- which(colnames(surv_data_all) == tempe)
                temps <- surv_data_all[,c(1,2,3,wh)]
                flag <- rep('',length(temps[[4]]))
                wh <- which(temps[[4]] > 0)
                flag[wh] <- 'Present'
                wh <- which(temps[[4]] <= 0)
                flag[wh] <- 'Missing'
                temps$flag <- flag
                ## check which edges have at least 10 patients in each category
                wh1 <- which(temps$flag == 'Missing')
                wh2 <- which(temps$flag == 'Present')
                s_edges <- c(s_edges, sgg[i])

                if((length(wh1) >= 10) & (length(wh2) >= 10)){
                    ##--- survival plot
                    sdf <- survdiff(Surv(overall_survival, deceased) ~ flag, data=temps)
                    tt <- signif(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
                    pvs <- c(pvs, tt)

                    if(rnd_run > 0){
                        temp_rnd <- random_surv(overall_survival, temps, rnd_run)
                        ppvs[[counter]] <- temp_rnd
                        counter <- counter+1
                    }

                }else{
                    pvs <- c(pvs, 1)
                    ppvs[[counter]] <- rep(1, rnd_run)
                    counter <- counter+1
                }
            }
        }
        if(rnd_run > 0){
            return(list(s_edges, pvs, ppvs))
        }else{
            return(list(s_edges, pvs))
        }   
}

## CALLED IN: 
# 2_select_edges_contact.r
# 2_select_edges_pisa.r
# 2_select_edges_EPPIC.r
# 2_analysis_individual_final.r
save_GL_EEI <- function(sgg,nflag,outdir3,cthres,c_type,tempw,ntf){ 
    tempnet <- igraph::as_data_frame(sgg)
    e1 <- tempnet[[1]]
    e2 <- tempnet[[2]]
    tpp <- mapProtein(e1,e2,ntf)
    tppx <- mapProtein(tpp[[1]], tpp[[2]], tempw)
    tpp$weight <- tppx$weight
    tpp$patient <- tppx$patient
    tpp <- tpp[order(-tpp$patient), ]
    otd <- paste0(outdir3,'/threshold_',cthres)
    if(!dir.exists(otd)){
        dir.create(otd)
    }
    data.table::fwrite(tpp, paste0(otd,'/',c_type,'_',nflag,'.txt'), row.names=FALSE, quote=FALSE, sep='\t') ###--- checked
}

## CALLED IN: 
# 2_select_edges_contact.r
# 2_select_edges_pisa.r
# 2_select_edges_EPPIC.r
# 2_analysis_individual_final.r
survival_qval_random <- function(gl_data, cthres, tempfd_condition, surv_data, exons, rnd_run){ 

        temp_net <- igraph::as_data_frame(gl_data)
        samplenames <- unlist(lapply(strsplit(colnames(tempfd_condition), '[_]'), '[[', 1)) 
        wh1 <- setdiff(samplenames, surv_data$submitter_id)
        wh2 <- setdiff(surv_data$submitter_id, samplenames)
        if(length(wh1) > 0){
            wh3 <- which(samplenames %in% wh1)
            samplenames <- samplenames[-wh3]
            tempfd_condition <- tempfd_condition[,-wh3]
        }
        if(length(wh2) > 0){
            wh3 <- which(surv_data$submitter_id %in% wh2)
            surv_data <- surv_data[-wh3, ]
        }

        colnames(tempfd_condition) <- samplenames
        tempx <- as.data.frame(t(tempfd_condition))
        tempx <- tempx[order(rownames(tempx)), ]
        colnames(tempx) <- exons

        temp5x <- exon2eei(temp_net, tempx, cthres)
        surv_dataxx <- surv_data

        ppvs <- list()
        for(rrn in 1:rnd_run){

            ##-- randomize the edge assignments ---
            cnames <- colnames(temp5x)
            rnames <- rownames(temp5x)
            tempfr <- as.data.table(t(temp5x))
            tempfr <- tempfr[, lapply(.SD, sample)]
            temp5 <- as.data.frame(t(tempfr))
            rownames(temp5) <- rnames
            colnames(temp5) <- cnames
            ##------------------------------------

            surv_data <- surv_dataxx[order(surv_dataxx$submitter_id), c(1,6,7)]
            surv_data_all <- cbind(surv_data, temp5)

            sgg <- colnames(surv_data_all)
            sgg <- sgg[-c(1:3)]
            ##-- for each selected edges ---
            pvs <- c()
            s_edges <- c()

            if(length(sgg) != 0){
                for(i in 1:length(sgg)){
                    tempe <- sgg[i]
                    wh <- which(colnames(surv_data_all) == tempe)
                    temps <- surv_data_all[,c(1,2,3,wh)]
                    flag <- rep('',length(temps[[4]]))
                    wh <- which(temps[[4]] > 0)
                    flag[wh] <- 'Present'
                    wh <- which(temps[[4]] <= 0)
                    flag[wh] <- 'Missing'
                    temps$flag <- flag
                    ## check which edges have at least 10 patients in each category
                    wh1 <- which(temps$flag == 'Missing')
                    wh2 <- which(temps$flag == 'Present')

                    if((length(wh1) >= 10) & (length(wh2) >= 10)){
                        ##--- survival plot -----
                        sdf <- survdiff(Surv(overall_survival, deceased) ~ flag, data=temps)
                        tt <- signif(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)
                        pvs <- c(pvs, tt)
                        s_edges <- c(s_edges, sgg[i])
                    }else{
                        pvs <- c(pvs, 1)
                        s_edges <- c(s_edges, sgg[i])
                    }
                }
            }

            ppvs[[rrn]] <- pvs
        }

        return(list(s_edges, ppvs))
}

mapProtein <- function(e1,e2,netf){ ###---checked
    pos <- c()
    for(nn in 1:length(e1)){
        wh1 <- which(netf[[1]] == e1[nn])
        wh2 <- which(netf[[2]] == e2[nn])
        wh <- intersect(wh1, wh2)
        if(length(wh) == 0){
            wh1 <- which(netf[[2]] == e1[nn])
            wh2 <- which(netf[[1]] == e2[nn])
            wh <- intersect(wh1, wh2) 
        }
        pos <- c(pos, wh[1])
    }

    tpp <- netf[pos,]
    return(tpp)
}

## CALLED IN: 
# 2_select_edges_contact.r
# 2_select_edges_pisa.r
# 2_select_edges_EPPIC.r
# 2_analysis_individual_final.r
save_selected_EEI <- function(sgg,nflag,outdir3,cthres,c_type,pvs,tempw,ntf,frac1,frac2){ 
    e1 <- unlist(lapply(strsplit(sgg,'[_]'), '[[', 1))
    e2 <- unlist(lapply(strsplit(sgg,'[_]'), '[[', 2))
    tpp <- mapProtein(e1,e2,ntf)
    tpp$pval <- pvs
    tpp$frac1 <- frac1
    tpp$frac2 <- frac2
    tppx <- mapProtein(tpp[[1]], tpp[[2]], tempw)
    tpp <- cbind(tpp, tppx[,c(3,4)])
    tpp <- tpp[order(tpp$pval), ]
    otd <- paste0(outdir3,'/threshold_',cthres)
    if(!dir.exists(otd)){
        dir.create(otd)
    }
    data.table::fwrite(tpp, paste0(otd,'/',c_type,'_',nflag,'.txt'), sep='\t', row.names=FALSE, quote=FALSE)
}




barPlotk <- function(alldata, flag, cthres, odir, cancer_type){ ##--- checked

    if(!dir.exists(odir)){dir.create(odir)}
    edge_type <- unique(alldata$pt)
    pData <- data.frame(matrix(ncol=0,nrow=0))

    for(k in 1:length(cancer_type)){
        tempd1 <- alldata[alldata$temp_cancer == cancer_type[k], ]
        for(j in 1:length(edge_type)){
            tempd3 <- tempd1[tempd1$pt == edge_type[j], ]
            if(nrow(tempd3) > 10){
                tempd3 <- tempd3[order(tempd3$qval), ][1:10,]
            }
            pData <- rbind(pData, tempd3)
        }
    }

    pData$val <- -log(pData$qval)

    outdir <- paste0(odir,'/Enrichment_KEGG','_',flag,'_',cthres)
    if(!dir.exists(outdir)){dir.create(outdir, recursive=TRUE)}
    ppData <- pData
        basesize <- 10
        p <- ggplot(ppData, aes(GOtermn, val, fill=pt)) +
        geom_bar(stat="identity", color="black", position=position_dodge()) +
          theme_grey(base_size = basesize * 0.8) + labs(x = "Cancer type", y = "# of edges") +
          scale_y_continuous()+
          scale_x_discrete()+
          guides(fill=guide_legend(title='-log(q-value)'))+
          theme(axis.text.x = element_text(size = basesize * 0.8,
                                                          angle = 90, hjust = 1,vjust=0.5, colour = "grey50"))+
          facet_wrap(~temp_cancer, scale='free', ncol=4)+
          coord_flip()+
          theme(legend.text=element_text(size=basesize * 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        ggsave(p,filename=paste0(outdir,'/KEGG_enrich.png'),width=16, height=8, dpi=400)
}

barPlot <- function(alldata, flag, cthres, odir, cancer_type){ ##--- checked 

    if(!dir.exists(odir)){dir.create(odir)}
    edge_type <- unique(alldata$pt)
    gtype <- unique(alldata$flag)
    pData <- data.frame(matrix(ncol=0,nrow=0))

    for(k in 1:length(cancer_type)){
        tempd1 <- alldata[alldata$temp_cancer == cancer_type[k], ]
        for(i in 1:length(gtype)){
            tempd2 <- tempd1[tempd1$flag == gtype[i], ]
            for(j in 1:length(edge_type)){
                tempd3 <- tempd2[tempd2$pt == edge_type[j], ]
                if(nrow(tempd3) > 20){
                    tempd3 <- tempd3[order(tempd3$qval), ][1:20,]
                }
                pData <- rbind(pData, tempd3)
            }
        }
    }

    pData$val <- -log(pData$qval)

    for(k in 1:length(gtype)){
        outdir <- paste0(odir,'/Enrichment_',gtype[k],'_',flag,'_',cthres)
        if(!dir.exists(outdir)){dir.create(outdir, recursive=TRUE)}
        ppData1 <- pData[pData$flag == gtype[k], ]
        for(j in 1:length(cancer_type)){
            ppData <- ppData1[ppData1$temp_cancer == cancer_type[j], ]
            basesize <- 10
            p <- ggplot(ppData, aes(GOtermn, val, fill=pt)) +
            geom_bar(stat="identity", color="black", position=position_dodge()) +
              theme_grey(base_size = basesize * 0.8) + labs(x = "GO term", y = "-log(q-value)") +
              scale_y_continuous()+
              scale_x_discrete()+
              guides(fill=guide_legend(title='Edge type'))+
              theme(axis.text.x = element_text(size = basesize * 0.8,
                                                              angle = 90, hjust = 1,vjust=0.5, colour = "grey50"))+
              facet_wrap(~temp_cancer, scale='free', ncol=4)+
              coord_flip()+
              theme(legend.text=element_text(size=basesize * 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
            ggsave(p,filename=paste0(outdir,'/',cancer_type[j],'.png'),width=7, height=5, dpi=400)
        
        }
    } 
}

goEnrich <- function(all_go_t, bgSize, flag, cancer_type){ ##---checked
    GOterm <- c()
    pval <- c()
    GOtermn <- c()
    temp_cancer <- c()
    for(k in 1:length(flag)){
        tempg <- igraph::as_data_frame(flag[k][[1]])
        uu <- union(tempg[[1]], tempg[[2]])
        tempd <- all_go_t[all_go_t$ensembl_exon_id %in% uu, ]
        sampleSize <- length(unique(tempd$ensembl_exon_id))
        goterms <- unique(tempd$go_id)

        for(j in 1:length(goterms)){
            setA <- unique(all_go_t[all_go_t$go_id == goterms[j], ]$ensembl_exon_id)
            setB <- unique(tempd[tempd$go_id == goterms[j], ]$ensembl_exon_id)
            ## hypergeometric test
            hyp <- phyper(length(setB)-1,length(setA),bgSize-length(setA), sampleSize,lower.tail = FALSE)
            pval <- c(pval, hyp)
            GOterm <- c(GOterm, goterms[j])
            GOtermn <- c(GOtermn, unique(tempd[tempd$go_id == goterms[j], ]$name_1006))
        }
        temp_cancer <- c(temp_cancer, rep(cancer_type[k], length(goterms)))
    }
    td <- data.frame(temp_cancer, GOterm, pval, GOtermn)
    return(td)
}


goEnrich_uniset <- function(all_go_t, bgSize, flag){ ##---checked
    GOterm <- c()
    pval <- c()
    GOtermn <- c()
    temp_cancer <- c()

    tempg <- igraph::as_data_frame(flag)
    uu <- union(tempg[[1]], tempg[[2]])
    tempd <- all_go_t[all_go_t$ensembl_exon_id %in% uu, ]
    sampleSize <- length(unique(tempd$ensembl_exon_id))
    goterms <- unique(tempd$go_id)

    for(j in 1:length(goterms)){
        setA <- unique(all_go_t[all_go_t$go_id == goterms[j], ]$ensembl_exon_id)
        setB <- unique(tempd[tempd$go_id == goterms[j], ]$ensembl_exon_id)
        ## hypergeometric test
        hyp <- phyper(length(setB)-1,length(setA),bgSize-length(setA), sampleSize,lower.tail = FALSE)
        pval <- c(pval, hyp)
        GOterm <- c(GOterm, goterms[j])
        GOtermn <- c(GOtermn, unique(tempd[tempd$go_id == goterms[j], ]$name_1006))
    }

    td <- data.frame(GOterm, pval, GOtermn)
    return(td)
}

survival_get_direction <- function(c_type, gl_data, cthres, all_exr, all_sr){ ##--- checked

        ename <- paste0(gl_data[[1]],'\n',gl_data[[2]])
        surv_data <- all_sr
        temp <- all_exr
        tnet <- gl_data
        temp1 <- temp[temp$EXON %in% union(tnet[[1]], tnet[[2]]), ]
        exons <- temp1$EXON
        temp2 <- temp1[,EXON:=NULL]
        temp2 <- as.data.frame(temp2)
        temp_brk <- unlist(lapply(strsplit(colnames(temp2), '[_]'), '[[', 2))
        sample_names <- unlist(lapply(strsplit(colnames(temp2), '[_]'), '[[', 1))
        control_pos <- which(temp_brk == '11')
        temp3 <- temp2[,-control_pos]
        samplenames <- sample_names[-control_pos]

        wh1 <- setdiff(samplenames, surv_data$submitter_id)
        wh2 <- setdiff(surv_data$submitter_id, samplenames)
        if(length(wh1) > 0){
            wh3 <- which(samplenames %in% wh1)
            samplenames <- samplenames[-wh3]
            temp3 <- temp3[,-wh3]
        }
        if(length(wh2) > 0){
            wh3 <- which(surv_data$submitter_id %in% wh2)
            surv_data <- surv_data[-wh3, ]
        }

        colnames(temp3) <- samplenames

        temp4 <- as.data.frame(t(temp3))
        temp4 <- temp4[order(rownames(temp4)), ]
        colnames(temp4) <- exons
        temp5 <- exon2eei(tnet, temp4, cthres)
        vars <- colnames(temp5)
        surv_data <- surv_data[order(surv_data$submitter_id), c(1,6,7)]
        temps <- cbind(surv_data, temp5)

        flag <- rep('',length(temps[[4]]))
        wh <- which(temps[[4]] == 1)
        flag[wh] <- 'Present'
        wh <- which(temps[[4]] == 0)
        flag[wh] <- 'Missing'
        temps$flag <- flag

        ##--- survival plot
        fit <- survfit(Surv(overall_survival, deceased) ~ flag, data=temps)
        sdf <- survdiff(Surv(overall_survival, deceased) ~ flag, data=temps)
        tt <- signif(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)

        missing <- summary(fit)$table[, "median"][1][[1]]
        present <- summary(fit)$table[, "median"][2][[1]]

        if(is.na(present)){
            present <- 100000
        }

        if(is.na(missing)){
            missing <- 100000
        }

        ## whether missing or presence of the EEI is related to lower survival
        if(missing < present){
            flag <- 'Favorable'
        }else{
            flag <- 'Unfavorable'
        }

        return(flag)
 
}


survival_plot_new <- function(c_type, gl_data, cthres, all_exr, all_sr, outdir, nflag, frac){ ##--- checked

        ename <- paste0(gl_data[[1]],'\n',gl_data[[2]])
        surv_data <- all_sr
        temp <- all_exr
        tnet <- gl_data
        temp1 <- temp[temp$EXON %in% union(tnet[[1]], tnet[[2]]), ]
        exons <- temp1$EXON
        temp2 <- temp1[,EXON:=NULL]
        temp2 <- as.data.frame(temp2)
        temp_brk <- unlist(lapply(strsplit(colnames(temp2), '[_]'), '[[', 2))
        sample_names <- unlist(lapply(strsplit(colnames(temp2), '[_]'), '[[', 1))
        control_pos <- which(temp_brk == '11')
        temp3 <- temp2[,-control_pos]
        samplenames <- sample_names[-control_pos]

        wh1 <- setdiff(samplenames, surv_data$submitter_id)
        wh2 <- setdiff(surv_data$submitter_id, samplenames)
        if(length(wh1) > 0){
            wh3 <- which(samplenames %in% wh1)
            samplenames <- samplenames[-wh3]
            temp3 <- temp3[,-wh3]
        }
        if(length(wh2) > 0){
            wh3 <- which(surv_data$submitter_id %in% wh2)
            surv_data <- surv_data[-wh3, ]
        }

        colnames(temp3) <- samplenames

        temp4 <- as.data.frame(t(temp3))
        temp4 <- temp4[order(rownames(temp4)), ]
        colnames(temp4) <- exons
        temp5 <- exon2eei(tnet, temp4, cthres)
        vars <- colnames(temp5)
        surv_data <- surv_data[order(surv_data$submitter_id), c(1,6,7)]
        temps <- cbind(surv_data, temp5)

        flag <- rep('',length(temps[[4]]))
        wh <- which(temps[[4]] == 1)
        flag[wh] <- 'Present'
        wh <- which(temps[[4]] == 0)
        flag[wh] <- 'Missing'
        temps$flag <- flag

        ##--- survival plot
        fit <- survfit(Surv(overall_survival, deceased) ~ flag, data=temps)
        sdf <- survdiff(Surv(overall_survival, deceased) ~ flag, data=temps)
        tt <- signif(1 - pchisq(sdf$chisq, length(sdf$n) - 1), 3)

        png(filename=paste0(outdir,'/',nflag,'.png'),
             width     = 3.5,
              height    = 3.5,
              units     = "in",
              res       = 600,
              pointsize = 2
              )
        p <- ggsurvplot(
               fit,                     # survfit object with calculated statistics.
               data = temps,             # data used to fit survival curves.
               risk.table = FALSE,       # show risk table.
               legend.title=ename,
               pval = FALSE,             # show p-value of log-rank test.
               xlab = "Time in days",   # customize X axis label.
               break.time.by = 500,     # break X axis in time intervals by 500.
               ggtheme = theme_light(), # customize plot and risk table with a theme.
              risk.table.y.text.col = T,# colour risk table text annotations.
              risk.table.height = 0.25, # the height of the risk table
              risk.table.y.text = FALSE,# show bars instead of names in text annotations
                                        # in legend of risk table.
              tables.theme = theme_cleantable(),
              surv.median.line = "hv",  # add the median survival pointer.
              legend.labs = c("Missing", "Present")    # change legend labels.
            )
        p$plot <- p$plot+ ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.5, vjust = 2,label = paste0("p-value: ",signif(tt,3)), size = 3)+
        ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

        print(p)
        dev.off()
}



exon2eei <- function(temp_net, tempx, cthres){ ###--- checked
    temp_mat  <- matrix(ncol=length(temp_net[[1]]), nrow=length(tempx[[1]]), 0)
    for(k in 1:length(tempx[[1]])){
        temp <- tempx[k, ]
        for(j in 1:length(temp_net[[1]])){
            t1 <- temp_net[[1]][j]
            t2 <- temp_net[[2]][j]
            expr1 <- temp[,which(colnames(temp) == t1)]
            expr2 <- temp[,which(colnames(temp) == t2)]
            if((expr1 >= cthres) & (expr2 >= cthres)){temp_mat[k,j] <- 1}
        }
    }
    temp_mat <- as.data.frame(temp_mat)
    colnames(temp_mat) <- paste0(temp_net[[1]],'_',temp_net[[2]])
    rownames(temp_mat) <- rownames(tempx)
    return(temp_mat) 
}

