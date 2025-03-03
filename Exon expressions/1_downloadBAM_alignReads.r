##############################################################################################
# Purpose: select the samples from a TCGA cancer type to process
##############################################################################################

args=commandArgs(TRUE)

library(edgeR)
library(data.table)
library(GenomicDataCommons)
library(plyr)
library(Rsubread)
library(Rsamtools) 

token <- '' ### add path to your GDC token to download .bam files
annotation_file <- '../data/Homo_sapiens.GRCh38.105_exons.gtf'
bam_files_dir <- '../data/all-bams'
survperiod <- 0 # consider all patients


TCGAtranslateID <- function(file_ids, legacy = FALSE) {
    info <- files(legacy = legacy) %>%
        filter( ~ file_id %in% file_ids) %>%
        select('cases.samples.submitter_id') %>%
        results_all()
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list <- lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file <- sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list)))
}


TCGAconsideredCancerTypes <- function(){

    ##----- Number of cases per project ------------------------------------
    res <- cases() %>% facet("project.project_id") %>% aggregations()
    res1 <- res$project.project_id
    res1$temp <- substr(res1[[2]], 1,4)
    res2 <- res1[res1$temp == 'TCGA', ]
    temp_type <- unlist(lapply(strsplit(res2$key, '[-]'),'[[',2))

    ##---- to store the cancer-types with at least 10 samples with both normal and cancer tissue
    c_type <- c()

    ##--- check for number of samples in both solid cancer and solid normal --------
    for(k in 1:length(temp_type)){

        ##----- Get manifest -----------------------------------------------
        ge_manifest = files() %>%
            filter( cases.project.project_id == paste0('TCGA-',temp_type[k])) %>% 
            filter( type == 'aligned_reads' ) %>%
            manifest()
        ge_manifest1 <- as.data.frame(ge_manifest)

        ##---- only RNA seq aligned to genome data
        temp_bam <- ge_manifest1
        wh3 <- which(temp_bam$filename %like% 'rna_seq.genomic')
        temp_bam1 <- temp_bam[wh3, ]

        #colnames(temp_bam1) <- c('file_id','filename','md5','size', 'state')
        temp_bam1 <- temp_bam1[,c(1,19,12)]
        colnames(temp_bam1) <- c('file_id','filename','md5')
        sids <- TCGAtranslateID(temp_bam1[[1]])
        temp_bam2 <- merge(temp_bam1, sids, by='file_id')

        temp_bam2$nid <- substr(temp_bam2$submitter_id, 1, 12)
        # only keep normal solid tissue (11) and primary cancer (01)
        temp_bam2$sample_id <- substr(temp_bam2$submitter_id, 14,15)
        temp_bam3 <- temp_bam2[temp_bam2$sample_id %in% c('01','11'), ]

        ## choose only those submitter ids that have at least one in each category 01 and 11
        tempx1 <- temp_bam3[temp_bam3$sample_id == '01', ]$nid
        tempx2 <- temp_bam3[temp_bam3$sample_id == '11', ]$nid
        tempx3 <- intersect(tempx1, tempx2)
        temp_bam4 <- temp_bam3[temp_bam3$nid %in% tempx3, ]

        # only keep the patients with one sample in each of normal solid tissue and primary cancer
        allnids <- unique(temp_bam4$nid)
        whc <- c()
        for(j in 1:length(allnids)){
            temp_wh <- which(temp_bam4$nid == allnids[j])
            temp <- temp_bam4[temp_wh, ]
            if(length(temp_wh) == 2){
                whc <- c(whc, temp_wh)
            }else{
                temp_wh1 <- which(temp$sample_id == '01')
                whc <- c(whc, temp_wh[temp_wh1[1]])

                temp_wh1 <- which(temp$sample_id == '11')
                whc <- c(whc, temp_wh[temp_wh1[1]])
            }
        }

        temp_final <- temp_bam4[whc, ]
        if(length(temp_final[[1]]) >= 20){
            c_type <- c(c_type, temp_type[k])
        }
    }

    return(c_type)
}


TCGAClinical <- function(ids,ct){

    case_ids = cases() %>%
        filter(~ project.project_id == paste0('TCGA-',ct)) %>%
        ids()
    clindat = gdc_clinical(case_ids)
    all_sub_cases <- ids
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
    return(diag_data_final3)
}

##-- which cancer types to concentrate on ----
c_s <-  args[1] ## 'BLCA'#
c_select <- intersect(c_s, TCGAconsideredCancerTypes())

if(length(c_select) > 0){

    out_dir <- paste0('../data/featureCounts_exons/', c_select)
    if(!dir.exists(out_dir)){
        dir.create(out_dir, recursive=TRUE)
    }

    ##----- Get manifest -----------------------------------------------
    ge_manifest = files() %>%
        filter( cases.project.project_id == paste0('TCGA-',c_select)) %>% 
        filter( type == 'aligned_reads' ) %>%
        manifest()
    ge_manifest1 <- as.data.frame(ge_manifest)

    ##---- only RNA seq aligned to genome data
    temp_bam <- ge_manifest1
    wh3 <- which(temp_bam$filename %like% 'rna_seq.genomic')
    temp_bam1 <- temp_bam[wh3, ]

    # colnames(temp_bam1) <- c('file_id','filename','md5','size', 'state')
    temp_bam1 <- temp_bam1[,c(1,19,12)]
    colnames(temp_bam1) <- c('file_id','filename','md5')
    sids <- TCGAtranslateID(temp_bam1[[1]])
    temp_bam2 <- merge(temp_bam1, sids, by='file_id')

    temp_bam2$nid <- substr(temp_bam2$submitter_id, 1, 12)
    # only keep normal solid tissue (11) and primary cancer (01)
    temp_bam2$sample_id <- substr(temp_bam2$submitter_id, 14,15)
    temp_bam3 <- temp_bam2[temp_bam2$sample_id %in% c('01','11'), ]

    ## choose only those submitter ids that have at least one in each category 01 and 11
    tempx1 <- temp_bam3[temp_bam3$sample_id == '01', ]$nid
    tempx2 <- temp_bam3[temp_bam3$sample_id == '11', ]$nid
    tempx3 <- intersect(tempx1, tempx2) ##-- patinets ids with both cancer and normal tissue

    ## choose those submitter ids that have at least one in 01 (cancer) category
    tempx1 <- temp_bam3[temp_bam3$sample_id == '01', ]$nid ##--- patient ids with cancer tissues
    tempx4 <- setdiff(tempx1, tempx3) ##--- patient ids with only cancer tissue

    temp_bam4 <- temp_bam3[temp_bam3$nid %in% tempx3, ]
    temp_bam5 <- temp_bam3[temp_bam3$nid %in% tempx4, ]

    # only keep the patients with one sample in each of normal solid tissue and primary cancer
    allnids <- unique(temp_bam4$nid)
    whc <- c()
    for(k in 1:length(allnids)){
        temp_wh <- which(temp_bam4$nid == allnids[k])
        temp <- temp_bam4[temp_wh, ]
        if(length(temp_wh) == 2){
            whc <- c(whc, temp_wh)
        }else{
            temp_wh1 <- which(temp$sample_id == '01')
            whc <- c(whc, temp_wh[temp_wh1[1]])

            temp_wh1 <- which(temp$sample_id == '11')
            whc <- c(whc, temp_wh[temp_wh1[1]])
        }
    }
    temp_final1 <- temp_bam4[whc, ]
    fwrite(temp_final1, paste0('../data/Manifests/',c_select,'_manifest_final.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    ##---- only keep the patients with one sample in primary cancer and no normal tissue sample
    allnids <- unique(temp_bam5$nid)
    whc <- c()
    for(k in 1:length(allnids)){
        temp_wh <- which(temp_bam5$nid == allnids[k])
        whc <- c(whc, temp_wh[1])
    }
    temp_final2x <- temp_bam5[whc, ]
    temp_finalx <- rbind(temp_final1, temp_final2x)

    ##----------------------------------------------------------------------------------------------------
    ##---- get clinical data
    survival_data <- TCGAClinical(temp_finalx$nid, c_select)
    ##--- filter survival data to include only those patients that have "complete" survival info
    survival_data <- survival_data[!is.na(survival_data$overall_survival), ]
    ##-- take all deceased patients
    dsdata <- survival_data[which(survival_data$vital_status == 'Dead'), ]
    ##-- take all non-deceased patients
    ndsdata <- survival_data[which(survival_data$vital_status == 'Alive'), ]
    ##--filter out those patients that are non-deceased but less than 5 years
    ndsdata1 <- ndsdata[ndsdata$overall_survival >= survperiod, ]
    ##-- all data
    surv_data <- rbind(dsdata, ndsdata1)

    ##--- filter to retain patients with only disease sample
    temp_final2 <- temp_final2x[temp_final2x$nid %in% surv_data$submitter_id, ]
    fwrite(temp_final2, paste0('../data/Manifests/',c_select,'_manifest_final_allCancer.txt'), sep='\t', row.names=FALSE, quote=FALSE)

    ##----- download bam files for some samples -------------------------
    temp_final <- rbind(temp_final1, temp_final2)

    for(j in 1:length(temp_final[[1]])){ ## for each sample, dowload it and count the features

        temp_xx <- temp_final[j, ]

        ##--- check whether normalized file is already present
        count_file <- paste0(out_dir,'/',temp_xx$nid, '_', temp_xx$sample_id,'.tsv')

        if(file.exists(count_file)){
            next
        }else{
            ##--- download the sample
            colnames(temp_xx) <- c('id','filename','md5','size','state','submitter_id','nid','sample_id')
            temp_file <- paste0('../data/',c_select,'_manifest_xx.txt')
            if(file.exists(temp_file)){
                file.remove(temp_file)
            }
            fwrite(temp_xx, temp_file, sep='\t', row.names=FALSE, quote=FALSE)
            temp_cmd <- paste0('./gdc-client download -n 5 -d ../data/all-bams -m ', temp_file,' -t ',token)
            system(temp_cmd)

            ##--- count the exon expression
            temp_name <- paste0(temp_xx$nid,'_',temp_xx$sample_id)
            bam_file <- paste0(bam_files_dir,'/',temp_xx$id,'/', temp_xx$filename)

            ##-- if the bam file was downloaded successfully ---
            ### some bam files are not able to be downloaded and are giving ""Division by Zero" error--

            if(file.exists(bam_file)){

                fl <- testPairedEndBam(bam_file)
                if(fl){flag <- 1}else{flag <- 0}
                temp <- paste0('sh ./bam2count.sh -a ',annotation_file,' -o ',out_dir,
                    ' -l ',flag,' -f ',temp_name,
                ' --removetemp --verbose -p 10 ', bam_file)
                system(temp)

                ##--- remove the bam file to save space
                temp_rm <- paste0('rm -r ',bam_files_dir,'/',temp_xx$id)
                system(temp_rm)

            }
            
        }

        cat('File', j, 'of', length(temp_final[[1]]), 'done\n')
    }
}


