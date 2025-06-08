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
tw.fdr = FALSE

# revisions
main.dir <- "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/250402_allEUR_revisions/"
outdir = 'revisions'
bulk.filter.coding = FALSE
tw.fdr = FALSE

# revisions coding
main.dir <- "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/250402_allEUR_revisions_CODING/"
outdir = 'revisions_coding'
bulk.filter.coding = TRUE
tw.fdr = FALSE

# 7 june 2025 revision
main.dir = "/sc/arion/projects/roussp01a/sanan/230420_TWAS_companionPipelines/combinationAnalyses/stackTWAS_V2/output/250602_allEUR_revisions_PLOT/"
outdir = '060725_revision'
bulk.filter.coding = FALSE
tw.fdr = TRUE

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
str(ready_for_GSEA_preparation_BULK)
save(ready_for_GSEA_preparation_BULK, file = paste0(lists.outdir, "/ready_for_GSEA_preparation_BULK.RData"))


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

if(tw.fdr){
  temp <- lapply(seq_along(temp), function(i){
    df = temp[[i]]
    df$p <- 10^(-(df$neg_log10p))
    df$fdr <- df$TW_fdr
    df$bonferroni <- df$TW_bonferroni
    df
})
}
for(i in seq_along(temp)){
  print(temp[[i]]$p[1:10])
}

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

a = c()
for(i in seq_along(bulk)){
  a = c(a, bulk[[i]]$p)
}
length(a)
length(unique(a))


# LOAD ALL TISSUES
load(paste0(lists.outdir, "/ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK.RData"))
raw_TWAS <- ready_for_GSEA_preparation_ALL_TISSUES_NO_BULK

a = c()
for(i in seq_along(raw_TWAS)){
  a = c(a, bulk[[i]]$p)
}
length(a)
length(unique(a))


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

# Necessary manual check of names
if(FALSE){
  names(everything_ready_TWAS_for_GSEA_preparation[[1]])
  df = everything_ready_TWAS_for_GSEA_preparation[[1]]
  unique(df$tissue)
}

### THIS STEP IS DATA-SPECIFIC BASED ON THE ABOVE (ALWAYS-FALSE) CONDITIONAL STATEMENT
# Mostly we check column pseudobulk, which should be

temp = lapply(everything_ready_TWAS_for_GSEA_preparation, function(x){
  x[tissue == 'EUR snBulk', tissue := "EUR Pseudobulk Bulk"]
})

names(temp) = names(everything_ready_TWAS_for_GSEA_preparation)
everything_ready_TWAS_for_GSEA_preparation = temp
str(everything_ready_TWAS_for_GSEA_preparation)

all.TWAS.forGSEA.fisher <- prepare.for.GSEA(
  mydf = everything_ready_TWAS_for_GSEA_preparation,
  filtering.criteria = list(c("p", 0.05), c("bonferroni", 0.05), c("p", 0.01), c("fdr", 0.05), c("fdr", 0.01)),
  save.path = paste0(lists.outdir, "TWAS_Fisher_input"),  #/multiple.TWAS.Combinations/ma.TWAS_corrected_BULK_V1", # has the last results 
  only.aggregates = T
)

str(all.TWAS.forGSEA.fisher)

ready_for_Fisher_analysis_TWAS_list <- all.TWAS.forGSEA.fisher

save(ready_for_Fisher_analysis_TWAS_list, file = paste0(lists.outdir, "/ready_for_Fisher_analysis_TWAS_list.RData"))

####################################
############# Analysis #############

load(paste0(lists.outdir, "/ready_for_Fisher_analysis_TWAS_list.RData"))
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/Imputable.Genes.List_V4.RData")
load("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Resources/Validation.Gene.List.R")
source('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/helper_clinical_GSEA.R')


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

