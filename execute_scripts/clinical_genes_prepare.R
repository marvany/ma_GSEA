source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")
library(ggplot2, lib.loc = "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS")
library(ggplot2)
library(scales)
library(pbmcapply)
setwd('/sc/arion/projects/va-biobank/marios/Clinical.Significance')
base.outdir = '/sc/arion/projects/va-biobank/PROJECTS/ma_gsea'

# we ran
# cp -r /sc/arion/projects/roussp01a/sanan/230420_TWAS_companionPipelines/combinationAnalyses/stackTWAS/output/250402_allEUR_revisions/ Resources/
# cp -r /sc/arion/projects/roussp01a/sanan/230420_TWAS_companionPipelines/combinationAnalyses/stackTWAS/output/250402_allEUR_revisions_CODING/ Resources/


# These are the raw TWAS for all tissues besides Bulk
main.dir <- "/sc/arion/projects/roussp01a/sanan/230420_TWAS_companionPipelines/combinationAnalyses/stackTWAS/output/240404_allEUR/"
outdir = 'classic'
bulk.filter.coding = TRUE


# revisions
main.dir <- "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/250402_allEUR_revisions/"
outdir = 'revisions'
bulk.filter.coding = FALSE


# revisions coding
main.dir <- "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/250402_allEUR_revisions_CODING/"
outdir = 'revisions_coding'
bulk.filter.coding = TRUE


# DEPENDENT VARIABLES
lists.outdir = paste0(base.outdir, "/Resources/", outdir, '/data_lists')
if(!dir.exists(lists.outdir)) dir.create(lists.outdir, recursive = TRUE)

gsea.path <- paste0(base.outdir, "/output/clinical_genes/", outdir)
if(!dir.exists(gsea.path)) dir.create(gsea.path, recursive = TRUE)

lists.outdir
gsea.path

##################################################################
##################### PREPARE VALIDATION GENESETS
if(FALSE){
  omim.cs <- 	readxl::read_excel(paste0(getwd(), "/Resources/clinvar.genes.omim.clinicalSynopsis.xlsx"))

  neuro <- unlist(omim.cs[!is.na(unlist(omim.cs[,"neurologicCentralNervousSystem"])),1])
  psych <- unlist(omim.cs[!is.na(unlist(omim.cs[,"neurologicBehavioralPsychiatricManifestations"])),1])

  neuro.psych <- c(neuro, psych)

  neuro.psych <- unique(only.protein.coding(neuro.psych, "Gene Symbol", "ENSEMBL"))


  validation.geneset <- update.geneset.list_2(list(protein.coding.neuro.psych = neuro.psych))
}

#View(validation.geneset)
######################################################
####### MODIFICATION OF BULK TISSUE TWAS RESULTS

p.bulk <- "/sc/arion/projects/roussp01a/sanan/230924_Various_TWAS/230926_PEC_DLPFC/output/240122_V4/data.TWAS/"

bulk.list <- read.dfs(p.bulk, acr = "csv", recursive = T)
bulk.list <- unlist(bulk.list, recursive = F)

bulk <- lapply(seq_along(bulk.list), function(i){
  temp <- bulk.list[[i]]
  temp <- as.data.table(temp)
  temp <- temp[temp$pred_perf_r2 >= 0.01 & temp$pred_perf_pval <= 0.05 & temp$n_snps_used > 0,]
  
  temp <- temp[pred_perf_qval <= 0.05,]
  indexMHC <- filterGenesMHC(temp$gene)
  if (length(indexMHC) > 0) {
    temp <- temp[-indexMHC,]
  }
  
  df <- temp
  df$tissue <- "Bulk"
  tryCatch({
  df$p <- df$pvalue
  df$bonferroni <- p.adjust(df$p, method = "bonferroni")
  df$fdr <- p.adjust(df$p, method = "fdr")
  },
  error = function(e){
    browser()
  })
  return(df)
})

# for one of the revisions we use an additional filter
coding.genes = geneAnnotation_v104_ensembl[Gene.Type == 'protein_coding',]$Ensembl.Gene.ID
if(bulk.filter.coding){
  bulk <- lapply(seq_along(bulk), function(i){
    
    dt <- as.data.table(bulk[[i]])
    dt = dt[gene %in% coding.genes]    
    
    dt
  })
}
names(bulk) <- names(bulk.list)
names(bulk) <- str_extract(names(bulk), ".*?(postImp)")
names(bulk) <- gsub("csv.", "", names(bulk))
str(bulk)
## Diastayrwneis ta onomata twn TWAS pou sou edwse o Sanan me ta TWAS pou exeis.
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/local.resources/coding_only_TWAS.genes.ID_V2.RData")
bulk <- bulk[names(bulk) %in% names(coding_only_TWAS.genes.ID_V2)]
names(bulk)
ready_for_GSEA_preparation_BULK <- bulk

save(ready_for_GSEA_preparation_BULK, file = paste0(lists.outdir, "/ready_for_GSEA_preparation_BULK.RData"))


##################################################################
##################################################################
##################### MODIFICATION OF TWAS RESULTS FOR ALL TISSUES.

raw_TWAS <- read.dfs(main.dir, acr = "tsv", recursive = T)
names(raw_TWAS$tsv) <- str_extract(names(raw_TWAS$tsv), ".*?(postImp)")

## Diastayrwneis ta onomata twn TWAS pou sou edwse o Sanan me ta TWAS pou exeis.
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/local.resources/coding_only_TWAS.genes.ID_V2.RData")

raw_TWAS <- raw_TWAS$tsv[names(raw_TWAS$tsv) %in% names(coding_only_TWAS.genes.ID_V2)]

temp <- lapply(seq_along(raw_TWAS), function(i){
  df = raw_TWAS[[i]]
  df$p <- 10^(-(df$neg_log10p))
  df$fdr <- df$pIB_fdr
  df$bonferroni <- df$pIB_bonferroni
  df
})
str(temp)
names(temp)
names(raw_TWAS)
names(temp) <- names(raw_TWAS)
ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK <- temp

save(ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK, file = paste0(lists.outdir, "/ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK.RData"))


##########################################################################
##########################################################################
####################### MERGE BULK AND ALL TISSUES #######################

# LOAD BULK
load(paste0(lists.outdir, "/ready_for_GSEA_preparation_BULK.RData"))
bulk <- ready_for_GSEA_preparation_BULK
str(bulk)

# LOAD ALL TISSUES
load(paste0(lists.outdir, "/ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK.RData"))
raw_TWAS <- ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK

cat("Total rows should be", sum(nrow(bulk$ADHD_2023_Demontis_postImp), nrow(raw_TWAS$ADHD_2023_Demontis_postImp)))

#coding_only_TWAS.genes.ID_V2 <- coding_only_TWAS.genes.ID_V2[order(names(coding_only_TWAS.genes.ID_V2))]

bulk <- bulk[order(names(bulk))]
raw_TWAS <- raw_TWAS[order(names(raw_TWAS))]

all(names(bulk) == names(raw_TWAS))

ma.raw_TWAS_2 <- lapply(names(bulk), function(name){
  bind_rows(bulk[[name]], raw_TWAS[[name]])
})
names(ma.raw_TWAS_2) <- names(bulk)
str(ma.raw_TWAS_2)

# Check names
all(names(ma.raw_TWAS_2) == names(raw_TWAS))
cat("Total rows are", nrow(ma.raw_TWAS_2$ADHD_2023_Demontis_postImp))

everything_ready_TWAS_for_GSEA_preparation <- ma.raw_TWAS_2
save(
  everything_ready_TWAS_for_GSEA_preparation,
  file = paste0(lists.outdir, "/everything_ready_TWAS_for_GSEA_preparation.RData")
  )
str(everything_ready_TWAS_for_GSEA_preparation)

######################################
############# Start here #############
# This is the raw list that contains all traits/tissues.
load(paste0(lists.outdir, "/everything_ready_TWAS_for_GSEA_preparation.RData"))

all.TWAS.forGSEA.fisher <- prepare.for.GSEA(
  mydf = everything_ready_TWAS_for_GSEA_preparation,
  filtering.criteria = list(c("p", 0.05), c("bonferroni", 0.05), c("p", 0.01), c("fdr", 0.05), c("fdr", 0.01)),
  save.path = paste0(lists.outdir, "TWAS_Fisher_input"),  #/multiple.TWAS.Combinations/ma.TWAS_corrected_BULK_V1", # has the last results 
  only.aggregates = T
)

ready_for_Fisher_analysis_TWAS_list <- all.TWAS.forGSEA.fisher


save(ready_for_Fisher_analysis_TWAS_list, file = paste0(lists.outdir, "/ready_for_Fisher_analysis_TWAS_list.RData"))

########################################
############# For NO_FOCUS #############

load(paste0(lists.outdir, "/ready_for_Fisher_analysis_TWAS_list.RData"))
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/Imputable.Genes.List_V4.RData")
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/Validation.Gene.List.R")


# Make sure you only have unique entries in each list # don't freak out if you don't, I believe the pipeline handles that, so just check
sapply(ready_for_Fisher_analysis_TWAS_list, function(x){
  all(sapply(x, function(y){
    sapply(y, function(z){
      length(z) == length(unique(z))
    })
  }))
})

if(!dir.exists(gsea.path)) dir.create(gsea.path, recursive = TRUE)


multi3_wrapper(
  twas.list = ready_for_Fisher_analysis_TWAS_list[3], #for NO_focus   #####all.TWAS.forGSEA.fisher #for FOCUS,
  backgrounds = list(NULL), #for NO_focus ####list(Imputable.Genes.NO.BULK, NULL) #for FOCUS,   # please do not name the elements included in this list
  gsea.path = gsea.path, # all_combinations_TWAS_GSEA_fisher_V3 uses protein specific background for Bulk and has improved aesthetics
  validation = Validation.Gene.List,
  defined.genesets =  4, # NULL for generalised run, 
  tissues_to_examine = c("/Bulk", "/sn-Pseudohomogenate", "/All EUR Class", "/All EUR Subclass", "/All EUR Class + Subclass"),
  show.title = F
)


multi3_wrapper <- function(twas.list, backgrounds, validation, gsea.path, defined.genesets = NULL, mylist = list(), tissues_to_examine, show.title){
  # Because backgrounds does NOT include named elements, we extract names computationally under the convention 
  # that these have not been named beforehand.
  # We take the string that includes all the items in the list.
  all.backgrounds.string<- gsub("list\\(|\\)", "", deparse(substitute(backgrounds)))
  
  # We create a vector that contains the different names
  all.background.names <- unlist(strsplit(all.backgrounds.string, ","))
  all.background.names <- gsub(" ", "", all.background.names)
  gsea.path.backup <- gsea.path
  for(i in seq_along(twas.list)){
    for(j in seq_along(backgrounds)){
      
      
      # same path for each iteration
      gsea.path <- gsea.path.backup
      
      # Load TWAS and background from respective lists
      twas <- twas.list[[i]]
      background <- backgrounds[[j]]
      
      # Extract twas and background names. These will be used for names of files and dirs
      # TWAS
      TWAS.name <- names(twas.list)[i]
      
      #background
      background.name <- all.background.names[j]
      # this is the single background that will be used.
      if(background.name == "NULL"){       # This relies on the fact that inside GV's function myGreatGenes will be loaded
        background = NA                    # There is room for upgrade in how you import the greatGenes, it is counterintuitive to insert NULL to get greatGenes...
        background.name = "GreatGenes"     # and why would the convention be null, instead of NA from the first place... so baaaaaaaaaaaaaaaaad.........
      }
      
      
      # Create directories
      p.all.plots <- paste0(gsea.path, "/All.Plots")
      if(!dir.exists(p.all.plots)) dir.create(p.all.plots, recursive = T)
      plot.title = paste(TWAS.name, background.name, sep = "_")
      gsea.path <- paste0(gsea.path, "/", plot.title)
      if(!dir.exists(gsea.path)) dir.create(gsea.path)
      
      
      # Prepare Validation Genesets
      # This should be done computationally using grep("protein", names(Validation.Gene.List))
      if(is.null(defined.genesets)){
        if(any(grepl("PROTEIN", TWAS.name))){
          genesets <- c(4, 5, 7)
        }else{
          genesets <- c(1, 2, 3)
        }
      }else genesets = defined.genesets
      
      
      
      # Choose between Bulk or Protein_Bulk. Either way your background genes will contain the name Bulk
      # if the twas is filtered only for protein      # if the background includes a "protein background"
      if(grepl("PROTEIN", TWAS.name, ignore.case = TRUE) && any(grepl("protein", names(background), ignore.case = TRUE))){
        
        # Remove Bulk from the list.
        background <- background[-grep("\\bBulk\\b", names(background))]
        
        # rename from Protein_Bulk to Bulk
        names(background)[grep("Protein_Bulk", names(background))] <- "Bulk"
        
        # if protein is found only in the background genes, then remove the protein backgrounds 
      }else if(any(grepl("protein", names(background), ignore.case = TRUE))){ 
        background <- background[-grep("protein", names(background), ignore.case = TRUE)]
      }
      # for debugging. first condition to skip NA. # then stop if you have 2 "bulk" entries in your background genes.
      if(length(background) != 1) if(length(grep("bulk", names(background), ignore.case = TRUE)) >= 2) stop("MA: Before you enter the fisher GSEA analysis in the 'background' argument you have >= 2 Bulk genes.")
      
      
      
      multi.tissues_multi.backgrounds_fisherGSEA(
        gsea.outdir = gsea.path,
        TWAS.genes.ID = twas, #ma.TWASp0.05PROTEIN, #five.tissues_coding.only_TWAS, #new_TWAS
        validation.geneset = validation[genesets],         #[c(1, 2, 3)],     #[c(1, 4, 5, 7)],
        background.genes =  background
      )
      
      total_sets <- length(genesets)
      
      
      only.psych.coding.only.gsea.paths <- create_tissue_paths(
        gsea.outdir = gsea.path,
        specific.tissue.dir = tissues_to_examine,# c("/Bulk", "/sn-Pseudohomogenate", "/All EUR Class", "/All EUR Subclass", "/All EUR Class + Subclass"), # c("/Bulk", "/sn-Pseudohomogenate", "/EUR Class-EN", "/EUR Class-Astro", "/EUR Class-Endo", "/EUR Class-IN", "/EUR Class-Immune", "/EUR Class-Mural", "/EUR Class-OPCs", "/EUR Class-Oligo"),   # , #c("/EUR_Bulk", "/All EUR Subclass", "/All EUR Class", "/All EUR Class&Subclass", "/230719_MegaAnalysis_EUR_superClass_bulk")
        tissues.to.trait.subdirs = "/result_textFiles_individual",     
        number_of_sets = total_sets   # remove this for christs' sake
        
        #names.of.sets = c("neuro.psych", "only.psych") 
        # they should be ordered according to their actual order of the exported tsv files.
      )
      
      #return(only.psych.coding.only.gsea.paths)
      
      require(ggplot2)
      
      
      multi.enrich.dfs <- multi_tissue_trait_OR_sig_v2(
        
        gsea.outdir = gsea.path, 
        multi.traits.multi.paths = only.psych.coding.only.gsea.paths,  # multiple.paths is a vector with all the paths for all tisues for a single trait
        
        plot.per.trait = F, # this draws the OR plot for each trait separately
        plot.all.traits = T,
        plot.height = 8,
        min.intersection  = 2,
        returnplot = T,
        saveplot = T,
        return.multi.enrich.dfs = F,
        title.of.plot = plot.title,
        all.plots.path = p.all.plots,
        show.title = show.title
      )
    }
  }
}



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
      semifinal.list <- pbmcapply::pbmclapply(mydf, function(df){
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
      }, mc.cores = parallel::detectCores() -2)

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
