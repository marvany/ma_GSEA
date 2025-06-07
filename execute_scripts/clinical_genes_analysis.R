source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")
library(ggplot2, lib.loc = "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS")
library(ggplot2)
library(scales)

setwd("/sc/arion/projects/va-biobank/marios/Clinical.Significance")



#################################################################
#################################################################
################### FUNCTIONS ###################################


# THIS IS RUN AT THE 'PREPARE VALIDATION GENESETS' PART
update.geneset.list_2 <- function(genes, original.geneset.list = NULL, genes.name = NULL){
  if(is.atomic(genes)){
    geneset.name <- genes.name
    if(is.null(original.geneset.list)) original.geneset.list <- list()
    list.of.genesets <- original.geneset.list
    list.of.genesets <- update.geneset.list(original.geneset.list = list.of.genesets, 
                                            genes.vector = genes, 
                                            geneset.name = geneset.name)
  }
  
  if(is.list(genes)){
    if(is.null(original.geneset.list)){
      list.of.genesets <- list()
      for(i in 1:length(genes)){
        geneset.name <- names(genes)[i] # each gene vector should be named!
        list.of.genesets <- update.geneset.list(original.geneset.list = list.of.genesets, 
                                                genes.vector = genes[[i]], 
                                                geneset.name = geneset.name)
        names(list.of.genesets)[i] <- geneset.name
      }
    }else{
      list.of.genesets <- original.geneset.list
      for(i in seq_along(genes)){
        geneset.name <- names(genes)[i]
        list.of.genesets <- update.geneset.list(original.geneset.list = list.of.genesets, 
                                                genes.vector = genes[[i]], 
                                                geneset.name = geneset.name)
        names(list.of.genesets)[length(list.of.genesets)] <- geneset.name
      }
    }
  }
  return(list.of.genesets)
}

###################################################
# THIS FUNCTION IS RUN AT THE END OF PRESENT SCRIPT
prepare.for.GSEA <- function(
    mydf,
    filtering.criteria = list(c("p", 0.05), c("bonferroni", 0.05), c("p", 0.01), c("fdr", 0.05), c("fdr", 0.01)),
    save.path,
    #only All Class, All subclass, Bulk, pseudobulk, all class + subclass
    only.aggregates = T,
    Focus= F){
  
  # This list will be exported as the final list that contains all the TWAS combinations
  all.TWAS.combs <- list()
  # In each iteration we want mydf will be loaded from mydf.backup so that each iteration gets the original file.
  mydf.backup = mydf
  all.TWAS.combs <- lapply(list(T,F), function(protein){
    if(protein) message("Setting argument 'protein' to TRUE") else message("Setting argument 'protein' to FALSE")
    for(i in seq_along(filtering.criteria)){
      
      mydf <- mydf.backup
      
      ######################## Filter ##########################
      
      input <- filtering.criteria[[i]]
      
      message("Running for ", input)
      
      if(protein) output <- paste(save.path, input[1], input[2], "PROTEIN.RData", sep = "_") else output <- paste(save.path, input[1], input[2], "MIXED.RData", sep = "_")
      
      mydf <- lapply(mydf, function(df){
        df <- df[df[[input[1]]] < input[[2]],]
        return(df)
      })
      
      ###################### prepare data for Sanan's project ######################
      # This is to avoid the Pseudobulk / Bulk overlap
      mydf <- lapply(mydf, function(df){
        df <- as.data.table(df)
        df[tissue == "EUR Pseudobulk Bulk", tissue := "Pseudobulk"]
      })
      
      if(!Focus) keys <- c("\\bSubclass\\b", "\\bClass\\b", "Pseudobulk", "Bulk")
      if(Focus) keys <- c("\\bSubclass\\b", "\\bClass\\b", "Pseudobulk")
      ### final list contains everything and it is readily availabe to be used as input in the GSEA Fisher pipeline.
      semifinal.list <- lapply(mydf, function(df){
        df <- as.data.table(df)
        mylist <- lapply(keys, function(name){
          patterns <- grep(name, unique(mydf[[1]]$tissue), value = T)
          genes <- sapply(patterns, function(i){
            df[tissue == i]$gene
          })
          if(length(patterns) == 1) {
            temp <- list()
            temp[[name]] <- genes
            return(temp)
          }
          return(genes)
        })
        mylist <- unlist(mylist, recursive = F)
        starting.length = length(mylist)
        mylist[["All EUR Class"]] <- unname(unlist(mylist[grep("\\bClass", names(mylist), ignore.case = T)]))
        mylist[["All EUR Subclass"]] <- unname(unlist(mylist[grep("subclass", names(mylist), ignore.case = T)]))
        mylist[["All EUR Class + Subclass"]] <- c(mylist[["All EUR Class"]], mylist[["All EUR Subclass"]])
        mylist[["sn-Pseudohomogenate"]] <- unname(unlist(mylist[grep("Pseudobulk", names(mylist), ignore.case = T)]))
        
        # this is to avoid the bulk overlap
        mylist[["All EUR Subclass"]] <- c(mylist[["All EUR Subclass"]], mylist$`EUR Class-OPCs`, mylist$`EUR Class-Astro`, mylist$`EUR Class-Endo`, mylist$`EUR Class-Oligo`)
        if(!Focus) mylist[["Bulk"]] <- unname(unlist(mylist[grep("Bulk", names(mylist), ignore.case = F)]))
        
        if(only.aggregates) mylist <- mylist[-c(1:(starting.length - 1))]
        
        return(mylist)
      })
      #names(five.tissues_coding.only_TWAS$BD_2021_Mullins_postImp)
      #View(semifinal.list)
      ######################################################################### 
      ######################### This last step is to create the final format.
      final.list <- list()
      if(protein){
        for(i in seq_along(semifinal.list)){
          #name <- names(semifinal.list)[i]
          final.list[[i]] <- unique(only.protein.coding(semifinal.list[[i]], "ENSEMBL", "ENSEMBL"))
        }
      }else{ # list will contain coding only
        for(i in seq_along(semifinal.list)){
          final.list[[i]] <- unique(semifinal.list[[i]])
        }
      }
      names(final.list) <- names(semifinal.list)
      
      
      if(protein){
        final.list <- lapply(semifinal.list, function(trait){
          sapply(trait, function(tissue){
            tissue <- unique(only.protein.coding(tissue, "ENSEMBL", "ENSEMBL"))
            tissue
          })
        })
      }else{
        final.list <- lapply(semifinal.list, function(trait){
          sapply(trait, function(tissue){
            tissue <- unique(tissue)
            tissue
          })
        })
      }
      
      
      #View(final.list)
      
      # Export
      mydata <- final.list
      #output
      save(mydata, file = output)
      #save(as.symbol(file.name), file = output)
      cat(output, "\n")
      
      # Save results in a list, on which GSEA_fisher will iterate
      name <- gsub(".RData", "", str_extract(output, "[^/]*.RData"))
      all.TWAS.combs[[name]] <- final.list
      
    }
    return(all.TWAS.combs)
  })
  # c(T,F) in lapply produce two lists, having 1 big list is what you need
  all.TWAS.combs <- unlist(all.TWAS.combs, recursive = F)
  
  # This exports and returns the list that contains all the combinations of TWAS criteria
  output = paste(save.path, "all.TWAS.combs.RData", sep = "_")
  save(all.TWAS.combs, file = output)
  cat("all.TWAS.combs can be found in this directory:", output)
  
  return(all.TWAS.combs)
}


# THIS IS RUN AT THE "MODIFICATION OF BULK TISSUE TWAS RESULTS"
read.dfs <- function(path, acr = c("gz", "csv", "tsv"), return.list = T, recursive = F){
  # extract the paths for all the items in the directory
  # these will be used for loading
  files <- list.files(path, recursive = recursive, full.names = T)
  
  # keep only the acronyms of the above paths
  # based on user's acr argument, these will be used to extract the indices
  # of the paths that will be loaded
  files.acr <- tools::file_ext(files) # the anchor $ tells R to start searching in the string backwards
  
  # if it encounters only directories it will exit.
  if(all(files.acr %in% "")) stop("No files ending in ", paste(acr), " were found in path \n", path)
  
  
  # we will load the data for each of the acr inserted
  mylist <- lapply(acr, FUN = function(thisacr){
    
    # find indices that match acr elements
    ind <- grep(thisacr, files.acr)
    
    # That is to avoid errors
    if(length(ind) == 0) return()
    
    # these.files is a list with the objects (tables for now) you want
    these.files <- pbmclapply(ind, function(i){
      fread(files[i])
    },
    mc.cores = parallel::detectCores() - 2)
    
    # the names of the variables match names of files without the acronym
    names(these.files) <- list.files(path, recursive = recursive, full.names = F)[ind]
    return(these.files)
  })
  #browser()
  if(!is.list(mylist)) return()
  names(mylist) <- acr
  #browser()
  return(mylist)
}

#################################################################
#################################################################
##################### PREPARE VALIDATION GENESETS

omim.cs <- 	readxl::read_excel(paste0(getwd(), "/Resources/clinvar.genes.omim.clinicalSynopsis.xlsx"))

neuro <- unlist(omim.cs[!is.na(unlist(omim.cs[,"neurologicCentralNervousSystem"])),1])
psych <- unlist(omim.cs[!is.na(unlist(omim.cs[,"neurologicBehavioralPsychiatricManifestations"])),1])

neuro.psych <- c(neuro, psych)


neuro.psych <- unique(only.protein.coding(neuro.psych, "Gene Symbol", "ENSEMBL"))

validation.geneset <- update.geneset.list_2(list(protein.coding.neuro.psych = neuro.psych))

######################################################
######################################################
####### MODIFICATION OF BULK TISSUE TWAS RESULTS

p.bulk <- "/sc/arion/projects/roussp01a/sanan/230924_Various_TWAS/230926_PEC_DLPFC/output/240122_V4/data.TWAS/"
bulk <- read.dfs(p.bulk, recursive = T)
bulk <- bulk[-1]
bulk <- unlist(bulk, recursive = F)

bulk <- lapply(bulk, function(temp){
  temp <- as.data.table(temp)
  temp <- temp[temp$pred_perf_r2 >= 0.01 & temp$pred_perf_pval <= 0.05 & temp$n_snps_used > 0,]
  indexMHC <- filterGenesMHC(temp$gene)
  if (length(indexMHC) > 0) {
    temp <- temp[-indexMHC,]
  }
  
  df <- temp
  df$tissue <- "Bulk"
  df$p <- df$pvalue
  df$bonferroni <- p.adjust(df$pvalue, method = "bonferroni")
  df$fdr <- p.adjust(df$pvalue, method = "fdr")
  return(df)
})

names(bulk) <- str_extract(names(bulk), ".*?(postImp)")
names(bulk) <- gsub("csv.", "", names(bulk))

## Diastayrwneis ta onomata twn TWAS pou sou edwse o Sanan me ta TWAS pou exeis.
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/local.resources/coding_only_TWAS.genes.ID_V2.RData")
bulk <- bulk[names(bulk) %in% names(coding_only_TWAS.genes.ID_V2)]
ready_for_GSEA_preparation_BULK <- bulk

#save(ready_for_GSEA_preparation_BULK, file = paste0(getwd(), "/ready_for_GSEA_preparation_BULK.RData"))
#load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/ready_for_GSEA_preparation_BULK.RData")
bulk <- ready_for_GSEA_preparation_BULK



##################################################################
##################################################################
##################### MODIFICATION OF TWAS RESULTS FOR ALL TISSUES.

# These are the raw TWAS for all tissues besides Bulk
main.dir <- "/sc/arion/projects/roussp01a/sanan/230420_TWAS_companionPipelines/combinationAnalyses/stackTWAS/output/240404_allEUR/"
raw_TWAS <- read.dfs(main.dir, acr = "tsv", recursive = T)
names(raw_TWAS$tsv) <- str_extract(names(raw_TWAS$tsv), ".*?(postImp)")

## Diastayrwneis ta onomata twn TWAS pou sou edwse o Sanan me ta TWAS pou exeis.
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/local.resources/coding_only_TWAS.genes.ID_V2.RData")

raw_TWAS <- raw_TWAS$tsv[names(raw_TWAS$tsv) %in% names(coding_only_TWAS.genes.ID_V2)]

ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK <- raw_TWAS
#save(ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK, file = paste0(getwd(), "/Resources/ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK.RData"))
#load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK.RData")

##########################################################################
##########################################################################
####################### MERGE BULK AND ALL TISSUES #######################

cat("Total rows should be", sum(nrow(bulk$ADHD_2023_Demontis_postImp), nrow(raw_TWAS$ADHD_2023_Demontis_postImp)))

coding_only_TWAS.genes.ID_V2 <- coding_only_TWAS.genes.ID_V2[order(names(coding_only_TWAS.genes.ID_V2))]
bulk <- bulk[order(names(bulk))]

ma.raw_TWAS_2 <- lapply(names(bulk), function(name){
  bind_rows(bulk[[name]], raw_TWAS[[name]])
})

names(ma.raw_TWAS_2) <- names(coding_only_TWAS.genes.ID_V2)

cat("Total rows are", nrow(ma.raw_TWAS_2$ADHD_2023_Demontis_postImp))

everything_ready_TWAS_for_GSEA_preparation <- ma.raw_TWAS_2
#save(everything_ready_TWAS_for_GSEA_preparation, file = "/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/everything_ready_TWAS_for_GSEA_preparation.RData")


##########################################################################
##################################################################
##################### MODIFICATION OF TWAS RESULTS FOR ALL TISSUES.

# These are the raw TWAS for all tissues besides Bulk
main.dir <- "/sc/arion/projects/roussp01a/sanan/230420_TWAS_companionPipelines/combinationAnalyses/stackTWAS/output/240404_allEUR/"
raw_TWAS <- read.dfs(main.dir, acr = "tsv", recursive = T)
names(raw_TWAS$tsv) <- str_extract(names(raw_TWAS$tsv), ".*?(postImp)")

## Diastayrwneis ta onomata twn TWAS pou sou edwse o Sanan me ta TWAS pou exeis.
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/local.resources/coding_only_TWAS.genes.ID_V2.RData")

raw_TWAS <- raw_TWAS$tsv[names(raw_TWAS$tsv) %in% names(coding_only_TWAS.genes.ID_V2)]

ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK <- raw_TWAS
#save(ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK, file = paste0(getwd(), "/Resources/ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK.RData"))
#load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK.RData")

##########################################################################
##########################################################################
####################### MERGE BULK AND ALL TISSUES #######################
cat("Total rows should be", sum(nrow(bulk$ADHD_2023_Demontis_postImp), nrow(raw_TWAS$ADHD_2023_Demontis_postImp)))

coding_only_TWAS.genes.ID_V2 <- coding_only_TWAS.genes.ID_V2[order(names(coding_only_TWAS.genes.ID_V2))]
bulk <- bulk[order(names(bulk))]

ma.raw_TWAS_2 <- lapply(names(bulk), function(name){
  bind_rows(bulk[[name]], raw_TWAS[[name]])
})

names(ma.raw_TWAS_2) <- names(coding_only_TWAS.genes.ID_V2)

cat("Total rows are", nrow(ma.raw_TWAS_2$ADHD_2023_Demontis_postImp))

everything_ready_TWAS_for_GSEA_preparation <- ma.raw_TWAS_2
#save(everything_ready_TWAS_for_GSEA_preparation, file = "/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/everything_ready_TWAS_for_GSEA_preparation.RData")


######################################
############# Start here #############
# This is the raw list that contains all traits/tissues.
# load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/everything_ready_TWAS_for_GSEA_preparation.RData")


all.TWAS.forGSEA.fisher <- prepare.for.GSEA(
  mydf = everything_ready_TWAS_for_GSEA_preparation,
  filtering.criteria = list(c("p", 0.05), c("bonferroni", 0.05), c("p", 0.01), c("fdr", 0.05), c("fdr", 0.01)),
  save.path = "/sc/arion/projects/va-biobank/marios/Clinical.Significance/temp/forsanan",  #/multiple.TWAS.Combinations/ma.TWAS_corrected_BULK_V1", # has the last results 
  only.aggregates = T
)

ready_for_Fisher_analysis_TWAS_list <- all.TWAS.forGSEA.fisher

#save(ready_for_Fisher_analysis_TWAS_list, file = paste0(getwd(), "/Resources/ready_for_Fisher_analysis_TWAS_list.RData"))





##########################################




