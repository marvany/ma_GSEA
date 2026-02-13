# 1. As of 261125 ; FOTIS results is the primary way to develop the input table

#############################################################################
# THIS SCRIPT IS THE FINAL PIPELINE, ITS PATH IS AT THE HOME/GLOBAL.SCRIPTS #
source("/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/060725_helper_Global.Source.R")
source("/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/060725_helper_GSEA_ORplot_Functions.R")
library(ggplot2)

meta.vs.nometa = "no_meta"
gwas.paths.list <- list(SCZ = "/sc/arion/projects/roussp01b/rachel/psychAD/work_data/disease_specific_without_meta_analysis_with_RUSH/RUSH_incl_no_meta_scz_SCZ_PRS_QDA",
                        AD = "/sc/arion/projects/roussp01b/rachel/psychAD/work_data/disease_specific_without_meta_analysis_with_RUSH/RUSH_incl_no_meta_ad_AD_PRS_QDA",
                        PD = "/sc/arion/projects/roussp01b/rachel/psychAD/work_data/disease_specific_without_meta_analysis_with_RUSH/RUSH_incl_no_meta_pd_PD_PRS_QDA",
                        BD = "/sc/arion/projects/roussp01b/rachel/psychAD/work_data/disease_specific_without_meta_analysis_with_RUSH/RUSH_incl_no_meta_bip_BIP_PRS_QDA")

# For Rachel initial
#################################################################################
#######################       NO META ANALYSIS       ############################
#######################
# FOR RACHEL FINAL
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/rachel_GSEA')
path_to_RUSH = paste0(getwd(), '/Resources/260122_genes')

meta.vs.nometa = "no_meta"
gwas.paths.list <- list(SCZ = paste0(path_to_RUSH, '/SCZ_case_SCZ_PGS'),
                        AD = paste0(path_to_RUSH, '/AD_case_AD_PGS'),
                        AD = paste0(path_to_RUSH, '/BD_case_BD_PGS'),
                        PD = paste0(path_to_RUSH, '/PD_case_PD_PGS'),
                        BD = paste0(path_to_RUSH, '/BD_case_BD_PGS'))

gene.symbol.column = "gene_id"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "p_value" #"mssm_p_value",
z.score.column = "z_score"
rachel = TRUE
p.export = '/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/rachel_GSEA/Results_260212'
# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))



########################
# FOR FOTIS
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/fotis_GSEA')
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/fotis_GSEA_v2')
path_to_genes = paste0(getwd(), '/Resources/genes/')
gwas.paths.list = list(WBS = path_to_genes)
p.export = "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/fotis_GSEA_v2/Results"

# Combine Fotis' genes (do once) - Fotis' specific
if(FALSE){
positive <- fread(list.files(path_to_genes, full.names = TRUE)[1])
negative <- fread(list.files(path_to_genes, full.names = TRUE)[2])
combined <- rbind(positive, negative)
fwrite(combined, file.path(path_to_genes, "combined_twas_genes.tsv"))
}

gene.symbol.column = "gene_name"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "pvalue" #"mssm_p_value",
z.score.column = "zscore"

# Skip first step
# you have to bring each table in a list and have it in this format
# df[, c('gene', 'feature', 'zscore', 'pvalue', 'fdr', 'bonferroni', 'gwas', 'model_ID')]
files <- list.files(path_to_genes, full.names = TRUE)

gwaslist <- lapply(files, function(f) {
  dt <- fread(f)

  # 1. Use 'gene' as symbol and match what ma.prepare.for.camera_4 does
  setnames(dt, "mol_name", "Approved.Symbol")

  # 2. Merge with annotation to get Ensembl IDs
  dt <- merge(
    dt,
    geneAnnotation_v104_ensembl,  # must exist: has Approved.Symbol + Ensembl.Gene.ID
    by = "Approved.Symbol",
    all.x = TRUE
  )

  # 3. Rename to camera-style columns
  setnames(dt, "Approved.Symbol", "gene_name")
  setnames(dt, "Ensembl.Gene.ID", "feature")

  # 4. Make 'gene' the Ensembl ID (like ma.prepare.for.camera_4)
  dt[, gene := feature]

  # 5. Add your engineered columns
  dt[, `:=`(
    pvalue     = 1e-50,
    zscore     = 10,
    fdr        = 1e-50,
    bonferroni = 1e-50,
    gwas       = "NLR",
    model_ID        = tools::file_path_sans_ext(basename(f))
  )]

  # 6. Keep only the desired columns in the desired order
  dt[, .(gene, feature, zscore, pvalue, fdr, bonferroni, gwas, model_ID)]
})
names(gwaslist) <- tools::file_path_sans_ext(basename(files))
str(gwaslist)



# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))


temp <- lapply(gwaslist, function(mydf){
  mydf$tissue <- gsub("_no_meta_analysis|_postImp","",mydf$model_ID)
  mydf$p <- mydf$pvalue
  return(mydf)
})
names(temp) <- names(gwaslist)
gwaslist <- temp

if(!dir.exists(p.export)) dir.create(p.export)

save(gwaslist, file = paste0(p.export, "/gwaslist.RData"))
# YOU CAN NOW GO STRAIGHT TO CHECKPOINT 1

########################
# FOR EIRINI
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/eirini_GSEA')
path_to_genes = paste0(getwd(), '/Resources/genes/')
dir.create(path_to_genes, recursive = T)
all_files = list.files('/sc/arion/projects/roussp01b/Eirini/Marios', full.names = T)

a = fread(all_files[1])

gwas.paths.list = list(WBS = path_to_genes)


gene.symbol.column = "gene_name"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "pvalue" #"mssm_p_value",
z.score.column = "zscore"

# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))

########################
# FOR SANAN
# THIS SCRIPT IS THE FINAL PIPELINE, ITS PATH IS AT THE HOME/GLOBAL.SCRIPTS #
source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")

library(ggplot2)

setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/sanan_mgi_GSEA')

path_to_snTIMs = paste0(getwd(), '/Resources/genes/v2')
list.files(path_to_snTIMs)
path_to_snTIMs <- "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/250402_allEUR_revisions"
list.files(path_to_snTIMs)
full.paths = list.files(path_to_snTIMs, full.names = TRUE, recursive = TRUE)
full.paths = full.paths[grep('\\.tsv', full.paths)]
gwas.paths.list = full.paths
names(gwas.paths.list) = basename(dirname(dirname(full.paths)))
#gwas.paths.list = lapply(full.paths, function(x) x)
#names(gwas.paths.list) = gsub('_.*','' , basename(full.paths))

# for ma.prepare.for.camera_v3
gene.symbol.column = "gene_name_noName"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "p" #"mssm_p_value",
z.score.column = "zscore"
feature.column = "gene"
fdr.column = "pIB_fdr"
bonferroni.column = "pIB_bonferroni"
# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))

p.export = '/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/sanan_mgi_GSEA'

#{
#  outdir = "/sc/arion/projects/roussp01b/Marios/GSEA_Rachel/Results/Fisher_V1"
#  outdir.backup <- outdir
#}

##############################################################
############## This is to normalise column names #############
source('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/060725_helper_camera.for.all.R')

# LATEST RELEASE FOR FASTER GENE ENSEMBL LOADING
gwaslist <- ma.prepare.for.camera_4(paths.list = gwas.paths.list, # expects list as input
                                    gene_symbol_col = gene.symbol.column,
                                    p.value.col = p.value.column, 
                                    z.score.col = z.score.column,
                                    return_dataframe = F)


str(gwaslist)


gwaslist <- ma.prepare.for.camera_2(paths.list = gwas.paths.list, # expects list as input
                                    gene_symbol_col = gene.symbol.column,
                                    p.value.col = p.value.column, 
                                    z.score.col = z.score.column,
                                    return_dataframe = F)
# THIS ONE WAS ESPECIALLY CREATED FOR SANAN
gwaslist <- ma.prepare.for.camera_3(paths.list = gwas.paths.list, # expects list as input
                                    gene_symbol_col = gene.symbol.column,
                                    p.value.col = p.value.column, 
                                    z.score.col = z.score.column,
                                    fdr.column = fdr.column,
                                    bonferroni.column = bonferroni.column,
                                    feature.column = feature.column,
                                    return_dataframe = F,
                                    convert.to.ensembl = F)



temp <- lapply(gwaslist, function(mydf){
  mydf$tissue <- gsub("_no_meta_analysis|_postImp","",mydf$model_ID)
  mydf$p <- mydf$pvalue
  return(mydf)
})
names(temp) <- names(gwaslist)
gwaslist <- temp


if(FALSE) rachels.gwaslist = gwaslist
if(FALSE) gwaslist = rachels.gwaslist
if(FALSE) sanans.gwaslist = gwaslist
if(FALSE) gwaslist = sanans.gwaslist

str(rachels.gwaslist[1])
str(sanans.gwaslist[1])

#gwaslist <- temp

if(!dir.exists(p.export)) dir.create(p.export)

save(gwaslist, file = paste0(p.export, "/gwaslist.RData"))



##########################################
############# CHECKPOINT 1 ###############
load(paste0(p.export, "/gwaslist.RData"))
source("/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/060725_helper_GSEA_ORplot_Functions.R")
source("/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/060725_helper_Global.Source.R")
library(ggplot2)


######################################################################
######################     ITERATION PART    #########################
# store dfs for different p-values
mylist <- lapply(1:nrow(thresh_annot), function(i){  
  # update pvalues per df
  p_val_gwaslist <- lapply(gwaslist, function(mydf){
    
    criteria_row <- as.vector(unlist(thresh_annot[i,]))
    # WHEN ALL SET UP DELETE THIS
    criteria_row <- rev(criteria_row)
    
    temp <- as.data.table(mydf[mydf[[criteria_row[1]]] < as.numeric(criteria_row[2]),])
  })
  
  # this is the list to be ran by multi3_wrapper()
  return(p_val_gwaslist)
})

names(mylist) <- sapply(1:nrow(thresh_annot), function(i) paste(rev(unlist(thresh_annot[i,])), collapse = "_"))
str(mylist)

# keep 0.05 , 0.01
mylist <- mylist[c(1,2)]

# This part checks for NA values
invisible(
  lapply(names(mylist), function(this.pval){
    thislist = mylist[[this.pval]] 
    lapply(names(thislist), function(this.name){
      thisdf = thislist[[this.name]]
      if(all(is.na(thisdf$gene))) stop('In ', this.pval, ' ', this.name,' there are NAs') else message('All good here')
    })
  })
)



filter_empty_traits <- function(mylist) {
  # Create a new list to store the filtered results
  filtered_list <- list()
  
  # Loop through each threshold type (fdr_0.05, bonferroni_0.05, etc.)
  for (threshold_name in names(mylist)) {
    threshold_list <- mylist[[threshold_name]]
    filtered_threshold <- list()
    
    # Loop through each trait (AD_2022_Bellenguez_postImp, ADHD_2023_Demontis_postImp, etc.)
    for (trait_name in names(threshold_list)) {
      trait_df <- threshold_list[[trait_name]]
      
      # Check if the data frame has features (rows)
      if (length(trait_df$feature) > 0) {
        filtered_threshold[[trait_name]] <- trait_df
      }
    }
    
    # Add the filtered threshold list to the result
    filtered_list[[threshold_name]] <- filtered_threshold
  }
  return(filtered_list)
}

mylist = filter_empty_traits(mylist)


str(mylist)
backup <- mylist

# CHECKPOINT
mylist <- backup

# Extract coding genes only from each table
temp <- pbmcapply::pbmclapply(names(mylist), function(name){
  thispval <- mylist[[name]]
  thispval <- pbmcapply::pbmclapply(thispval, function(thisdf){
    #browser()
    toreturn <- sapply(unique(thisdf$model_ID), function(tis){
      thesegenes <- thisdf[grep(tis, thisdf$model_ID), "gene"]$gene
      thesegenes <- unique(only.protein.coding(thesegenes, "ENSEMBL", "ENSEMBL")) # Keep only protein coding
    })
  }, mc.cores = parallel::detectCores() - 2)
  return(thispval)
}, mc.cores = parallel::detectCores() - 2)

names(temp) = names(mylist)
mylist = temp
####################################################
##################  Run Fisher  ####################

source("/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/060725_helper_GSEA_ORplot_Functions.R")
source("/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/060725_helper_Global.Source.R")

##############################################
################### FISHER ###################
twas.list = mylist #for NO_focus   #####all.TWAS.forGSEA.fisher #for FOCUS,
if(FALSE) rachels.twas.list = twas.list
if(FALSE) twas.list = rachels.twas.list

if(FALSE) sanan.twas.list = twas.list
if(FALSE) twas.list = sanan.twas.list


background = MultiWAS::greatGenes #for NO_focus ####list(Imputable.Genes.NO.BULK, NULL) #for FOCUS,   # please do not name the elements included in this list
#gsea.path = paste0(p.export, "/all_traits_V6") # all_combinations_TWAS_GSEA_fisher_V3 uses protein specific background for Bulk and has improved aesthetics
#gsea.path = paste0(p.export,"/fisher/MGIPhenotypes_v5_protein_coding")
gsea.path = file.path(p.export)
Validation_Genesets = mySets$standardGeneSets[c("syngoAll", "msigdbSetsPruned")] # Validation.Gene.List
#Validation_Genesets = list(MGIPhenotypes = MultiWAS::standardGeneSets$MGIPhenotypes)


defined.genesets =  NULL # NULL for generalised run, 
tissues_to_examine = NULL
show.title = T

p.all.plots <- paste0(gsea.path, "/All.Plots")
if(!dir.exists(p.all.plots)) dir.create(p.all.plots, recursive = T)

library(pbmcapply)
library(parallel)


if(FALSE){
# Top level - pval_name
lapply(
  names(twas.list), 
  function(pval_name){
    # pvalue level
    pval_level <- twas.list[[pval_name]]
    this_gsea_path <- paste0(gsea.path, "/", pval_name)
    
    # Second level - trait_name
    lapply(
      names(pval_level),
      function(trait_name){
        # trait level
        thistrait <- pval_level[[trait_name]]
        trait_gsea_path <- paste0(this_gsea_path, "/", trait_name)
        
        
        # Third level - geneset_name
        lapply(
          names(Validation_Genesets), 
          function(geneset_name) {
            Validation_Geneset <- Validation_Genesets[geneset_name]
            final_gsea_path <- paste0(trait_gsea_path, "/", geneset_name)
            # Run the analysis function
            fisherGsea_2_MA(
              testGenes = thistrait, 
              geneMetaSets = Validation_Geneset, 
              myGenes = background,
              outDir = final_gsea_path,
              tissueName = "sentinel",
              name.of.TWAS.genes = trait_name
            )
            
            return(NULL)  # Explicit return to avoid issues with parallel processing
        })
        
        return(NULL)  # Explicit return
      })
      
      return(NULL)  # Explicit return
})
} else{

### SIMPLIFIED FORMAT OF THE ABOVE : BEWARE IT ASSUMES "COLNAMES" AT SOME POINT (EXPECTING A MATRIX OR SOMETHING)
## Build a job grid: one row per (pval_name, trait_idx, trait_name, geneset_name)
job_grid_list <- lapply(names(twas.list), function(pval_name) {
  pval_level <- twas.list[[pval_name]]      # list of matrices
  
  if (!length(pval_level)) return(NULL)
  
  trait_idxs  <- seq_along(pval_level)
  
  # trait_name = column name of each matrix (negative_twas_genes / positive_twas_genes)
  trait_names <- vapply(trait_idxs, function(j) {
    cn <-names(pval_level)[j]
    if (!is.null(cn) && length(cn) >= 1 && nzchar(cn[1])) cn[1] else paste0("trait_", j)
  }, character(1))

 # For each trait, pair with all genesets
  rows <- lapply(trait_idxs, function(j) {
    data.frame(
      pval_name    = pval_name,
      trait_idx    = j,
      trait_name   = trait_names[j],
      geneset_name = names(Validation_Genesets),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
})

job_grid <- do.call(rbind, job_grid_list)
## Single-layer apply over the grid
lapply(seq_len(nrow(job_grid)), function(i) {
  pval_name    <- job_grid$pval_name[i]
  trait_idx    <- job_grid$trait_idx[i]
  trait_name   <- job_grid$trait_name[i]    # e.g. "negative_twas_genes"
  geneset_name <- job_grid$geneset_name[i]
  
  # Grab this trait's gene set (matrix -> vector)
  thesetraits <- twas.list[[pval_name]][[trait_idx]]
  if(rachel) these_model_IDs <- twas.list[[pval_name]][[trait_idx]]
  #thistrait     <- as.character(thistrait_mat[, 1])

#pbmcapply::pbmc
if(!rachel){
  lapply(names(thesetraits),
  # mc.cores = parallel::detectCores()-2,
  function(trait_name){  
    
  thistrait <- as.character(thesetraits[[trait_name]])
  # One geneset for this run
  Validation_Geneset <- Validation_Genesets[geneset_name]
  
  # Build paths: gsea.path / pval_name / trait_name / geneset_name
  this_gsea_path  <- file.path(gsea.path, pval_name)
  trait_gsea_path <- file.path(this_gsea_path, trait_name)
  final_gsea_path <- file.path(trait_gsea_path, geneset_name)


  fisherGsea_2_MA(
    testGenes          = thistrait,
    geneMetaSets       = Validation_Geneset,
    myGenes            = background,
    outDir             = final_gsea_path,
    tissueName         = trait_name,
    name.of.TWAS.genes = trait_name
  )
  })
  } else if(rachel){

    pbmcapply::pbmclapply(names(these_model_IDs),
    #mc.cores = parallel::detectCores()-2,
    function(this_model_ID){  
      
    thismodelID <- as.character(thesetraits[[this_model_ID]])
    # One geneset for this run
    Validation_Geneset <- Validation_Genesets[geneset_name]
    
    # Build paths: gsea.path / pval_name / trait_name / geneset_name
    final_gsea_path  <- file.path(gsea.path, pval_name, trait_name, this_model_ID, geneset_name)

    fisherGsea_2_MA(
      testGenes          = thismodelID,
      geneMetaSets       = Validation_Geneset,
      myGenes            = background,
      outDir             = final_gsea_path,
      tissueName         = trait_name,
      name.of.TWAS.genes = trait_name
    )
    })
  }
  invisible(NULL)
})
}





######################################################
################## PREPARE FOR PLOT ################## 

# Now we add the names manually...
syngo_name_of_pathway <- names(mySets$standardGeneSets["syngoAll"]$syngoAll$sets)
msigdb_name_of_pathway <- names(mySets$standardGeneSets$msigdbSetsPruned$sets)
mgi_phenotypes <- names(MultiWAS::standardGeneSets$MGIPhenotypes$sets)


### AS OF 112625 we skip this part.
if(FALSE){## You need to also include 0.01
mylist <- read.dfs.2(gsea.path)
str(mylist)
# this was not needed for Sanan's analysis; you only have 1 dataset
mylist <- mylist[[1]]
names(mylist)
mylist <- access_list(mylist,
                      condition = "any(grepl('individ', names(thisitem)))",
                      executable_code = "return(thisitem[2])"
                      )
str(mylist)


########
# SKIP #

mylist <- access_list(mylist,
                      condition = "any(grepl('Astro', names(thisitem)))",
                      executable_code = "if(any(grepl('ORplot', names(thisitem)))){
  thisitem <- thisitem[!grepl('ORplot', names(thisitem))]
  } else return(thisitem)"
)
names(mylist)
names(mylist[[1]][[1]])
names(mylist[[1]][[1]][[1]])
names(mylist[[1]][[1]][[1]][[1]])

mylist <- access_list(mylist,
                      condition = "any(grepl('tsv', names(thisitem)))",
                      executable_code = "
                    names(thisitem) <- gsub('syngoAll__|.genesets.tsv|msigdbSetsPruned__|MGIPhenotypes_', '', names(thisitem))
                    return(thisitem)")

names(mylist[[1]][[1]][[1]][[1]])


### NOT SURE WHERE THIS FITS
  # At this point you merge msigdb and synGO  
      mylist <- access_list(mylist,
                            condition = "any(grepl('AD', names(thisitem)))",
                            executable_code = "
                            short.list = lapply(names(thisitem), function(thisname){
                                thislist = thisitem[[thisname]]
                                thisdf = do.call('rbind', thislist)
                                thisdf$gwas = thisname
                                thisdf
                            })
                            names(short.list) = names(thisitem)
                            return(short.list)")
str(mylist)
################################################
### MAKE A FUNCTION OUT OF THIS !!!!!!!!!!!!!!!!


############# CHANGE NAME FOR 'search_for_this' VARIABLE
sanan = FALSE
if(sanan){
  # For MGIPhenotypes
  if(FALSE) mylist <- change_depth(thislist = mylist, target_depth = find_depth(mylist, 'MGIPhen') + 1, search_for_this = 'Class.Astro') # add 1 because 1 layer is lost



  if(TRUE){
        # merge MGIPhenotypes
        # At this point you merge msigdb and synGO
        mylist <- access_list(mylist,
                              condition = "any(grepl('MGIPheno', names(thisitem)))",
                              executable_code = "
                            thisitem <- do.call('rbind', thisitem)
                            return(thisitem)")
  }

mylist = mylist[[1]]

}

rachel = TRUE
if(rachel){# For syngoAll
mylist <- change_depth(thislist = mylist, target_depth = find_depth(mylist, 'syngoAll') + 1, search_for_this = 'Astro') # add 1 because 1 layer is lost

save(mylist, file = paste0(gsea.path, '/temp_mylist.RData'))

load(paste0(gsea.path, '/temp_mylist.RData'))

names(mylist)
names(mylist[[1]])
names(mylist[[1]][[1]])
names(mylist[[1]][[1]][[1]])
names(mylist[[1]][[1]][[1]][[1]])
    # At this point you merge msigdb and synGO
    mylist <- access_list(mylist,
                          condition = "any(grepl('individ', names(thisitem)))",
                          executable_code = "
                        thisitem <- do.call('rbind', thisitem)
                        return(thisitem)")

    mylist <- access_list(mylist,
                          condition = "any(grepl('syngo', names(thisitem)))",
                          executable_code = "
                        thisitem <- do.call('rbind', thisitem)
                        return(thisitem)")
    
#    mylist[[1]][[1]][['gwas']] <- NULL

names(mylist)
names(mylist[[1]])
names(mylist[[1]][[1]])
names(mylist[[1]][[1]][[1]])
names(mylist[[1]][[1]][[1]][[1]])

}

temp1 <- lapply(names(mylist), function(pval){
      
      thispval <- mylist[[pval]]
      temp2 <- lapply(names(thispval), function(trait){

          thistrait <- thispval[[trait]]
          thistrait[['gwas']] <- trait
          enrich.dfs <- thistrait
          enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = 'BH') #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
              enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = 'BH') #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
              enrich.dfs$star  <- as.character('')                                          # he creates an empty column and assigns an empty string everywhere
              enrich.dfs$star  <- ifelse(enrich.dfs$qvalue <= 0.001, '***',                 # then he accordingly annotates '***' based on the level of significance OF THE q-value
                                        ifelse(enrich.dfs$qvalue <= 0.01, '**',            # maybe change this to the FDR ?
                                                ifelse(enrich.dfs$qvalue <= 0.05, '*', '')))
              
              
              enrich.dfs <- enrich.dfs[order(enrich.dfs$pval, decreasing = F),]             # he organizes stuff based on pval order.. (maybe change this to FDR)
              '%!in%' <- function(x,y)!('%in%'(x,y))
              top_10_pathways <- unique(enrich.dfs$Reference)[1:10]
              
              enrich.dfs <- enrich.dfs[enrich.dfs$Reference %in% top_10_pathways]
              enrich.dfs$model_ID <- enrich.dfs$Set
              enrich.dfs      
      })
          
      names(temp2) <- names(thispval)
      return(temp2)
})

names(temp1) <- names(mylist)
mylist <- temp1
names(mylist)
names(mylist[[1]])
names(mylist[[1]][[1]])
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

} else { # YOU DEVELOPED THIS FOR FOTIS LAST REQUEST LET'S HOPE IT WORKS###########
    prep_gsea_for_plot <- function(mylist, top_n = 10) {
      
  # mylist: list(pvalue_0.01 = list(df1, df2, ...), pvalue_0.05 = list(...))

  out <- lapply(names(mylist), function(pval_name) {
    thispval <- mylist[[pval_name]]  # list of dfs for this pval

    inner <- lapply(names(thispval), function(trait_name) {
      df <- as.data.frame(thispval[[trait_name]])

      # If the original df is empty, just return it as-is
      if (nrow(df) == 0L) return(df)

      # 1. metadata
      df$gwas     <- trait_name       # e.g. "negative_twas_genes__syngoAll"
      df$model_ID <- df$Set
      if(rachel) df$model_ID <- trait_name # is actually the pvalue
      # 2. BH-FDR (qvalue) on raw pval
      df$qvalue <- p.adjust(df$pval, method = "BH")

      # 3. stars based on qvalue
      df$star <- ""
      df$star[df$qvalue <= 0.05]  <- "*"
      df$star[df$qvalue <= 0.01]  <- "**"
      df$star[df$qvalue <= 0.001] <- "***"

      # 4. sort by pval and keep top N unique pathways
      df <- df[order(df$pval, decreasing = FALSE), , drop = FALSE]

      uniq_refs <- unique(df$Reference)
      if (length(uniq_refs) > 0L) {
        top_refs <- uniq_refs[seq_len(min(top_n, length(uniq_refs)))]
        df <- df[df$Reference %in% top_refs, , drop = FALSE]
      }

      df
    })

    names(inner) <- names(thispval)
    inner
  })

  names(out) <- names(mylist)
  out
}



  # 1. Find all individual-level result files
indiv_files <- list.files(
  gsea.path,
  pattern = "__myTestGenes\\.genesets\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)



# 2. Build mylist: top level = pvalue threshold, second level = gwas+geneset combo
mylist <- list()

for (f in indiv_files) {
 # message(f)
  # Walk up the directory tree
  d0 <- dirname(f)    # .../geneset/result_textFiles_individual
  d1 <- dirname(d0)   # .../geneset
  d2 <- dirname(d1)   # .../negative_twas_genes or positive_twas_genes
  d3 <- dirname(d2)   # .../pvalue_0.05
  d4 <- dirname(d3)
  
  geneset_name <- basename(d1)  # syngoAll / msigdbSetsPruned
  gwas_sign    <- basename(d2)  # negative_twas_genes / positive_twas_genes
  if (rachel) model_ID <- basename(d2)
  gwas_name    <- basename(d3)  # pvalue_0.01 / pvalue_0.05
  pval_name    <- basename(d4)  # pvalue_0.01 / pvalue_0.05
  
  # Read the individual-level GSEA table
  df <- fread(f)
  df$gwas = gwas_name
  df$Set = model_ID
  # Make sure it has the expected columns (your example already does)
  # Required for your plotting function: Set, Reference, pval, odds.ratio, intersection,
  # union, FDR_AdjP, Jaccard, Bonf_AdjP, BH_AdjP, name_full, group, OR.ci95.down, OR.ci95.up

  # Construct a name for this df inside the pval layer
  df_name <- paste(gwas_name, model_ID, sep = "__")
  # e.g. "negative_twas_genes__syngoAll"

  if (is.null(mylist[[pval_name]])) {
    mylist[[pval_name]] <- list()
  }
  message(pval_name, " ", df_name)
  mylist[[pval_name]][[df_name]] <- df
}

# Optional: look at structure
str(mylist, max.level = 2)

mylist_trimmed <- prep_gsea_for_plot(mylist, top_n = 10)
mylist <- mylist_trimmed

}
#df <- mylist$pvalue_0.01_backgrounds

#########3 make similar with actual thing
#loadma("/sc/arion/projects/roussp01b/Marios/GSEA_Rachel/Results/V3_meta/0.05_pvalue/ready_to_plot_df.RData")
#for(i in names(df)){
#browser()
#  cat(i, df[[i]][1:5], "\n\n")
#}
#names(multi.enrich.dfs)

#################################################################
######################    CREATE PLOT   #########################
require(scales)
require(ggrepel)
require(ggplot2, lib.loc = "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS")
require(scales, lib.loc = "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS")
all.plot.ggpath <- paste0(p.export, "/All.Plots")
#all.plot.ggpath <- paste0(p.export, "/Results", "/all.plots_protein_only")

# when this lapply is ready you should include the heatmap thing and delete the df <- mylist.. above
# df_name: character scalar containing the line you want to append (e.g., "EUR_Astro")
no_enrichment_file <- file.path(gsea.path, "no_enrichment.tsv")
# Open list with p_value thresholds
graphs_per_pvalue <- vector("list", length(names(mylist)))
names(graphs_per_pvalue) <- names(mylist)

for (name_p in names(mylist)) {

  thislist <- mylist[[name_p]] # thislist contains all the traits-dfs

  # for conveniency
  p.value.thresh <- name_p

  # directory to store per p-value
  spec.ggpath <- paste0(all.plot.ggpath, "/", p.value.thresh)

  # define directory to store
  graph <- vector("list", length(names(thislist)))
  names(graph) <- names(thislist)
browser()
  for (df_name in names(thislist)) {

    df <- thislist[[df_name]]
    message(df_name)
    #if(df_name == 'AD__OPC_no_meta_analysis') browser()
    # title of plot
    plot_title <- paste0(name_p, "_", df_name)

    # name of file
    ggname.png <- paste0(plot_title, ".png")
    ggname.pdf <- paste0(plot_title, ".pdf")


    df$pathway <- df$Set
    df$Pvalue <- df$pval
    df$FDR <- df$FDR_AdjP
    df$neg_log10p <- -log10(df$Pvalue)
    df$model_ID <- df$name_full
    df$star <- ifelse(df$FDR < 0.05, "*", "")

    if(sum(df$intersection <2) >0) df[df$intersection<2,]$neg_log10p <- NA
    if (sum(df$intersection < 2) > (nrow(df) - 2)) {
      cat(df_name, "\n", file = no_enrichment_file, append = TRUE)
      next
    }

    colorVals <- c(min(df$neg_log10p,na.rm=T),
                   max(df$neg_log10p,na.rm=T))
    colorVals2 <- rescale(colorVals,to=c(0,1))
    colorLims <- c(0,ceiling(colorVals[2]))
    colorBreaks <- quantile(df$neg_log10p,na.rm=T,probs=seq(0,1,by=0.25))
    colorBreaks[1] <- 0

    df$model_ID <- df$Set


    # The palette with grey:
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # The palette with black:
    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


    df$Set <- gsub("_no_meta_analysis", "", df$Set)
    if(FALSE){
      df[df$Set == "Astro",]$Set <- "Astrocytes"
      df[df$Set == "EN",]$Set <- "Excitatory Neurons"
      df[df$Set == "Endo",]$Set <- "Endothelial"
      df[df$Set == "IN",]$Set <- "Inhibitory Neurons"
      df[df$Set == "Oligo",]$Set <- "Oligodendrocytes"
    }
    df$Reference <- gsub("_", " ", df$Reference)
    df$Reference <- ifelse(str_extract(df$Reference, "^.{1,2}") == "GO",
                           paste0(gsub("GO ", "", df$Reference), " (GO)"),
                           paste0(gsub("REACTOME", "", df$Reference), " (RE)")
    )

    df$Reference <- stringr::str_to_title(df$Reference)
    df$Reference <- gsub("\\(Go\\)", "(GO)", df$Reference)

    df$Reference <- gsub("Gtpase", "GTPase", df$Reference)
    df$Reference <- gsub("C Terminal", "C-Terminal", df$Reference)
    df$Reference <- gsub("Sema3a", "Sema3A", df$Reference)
    df$Reference <- gsub('Wnt5a', 'WNT5A', df$Reference)
    df$Reference <- gsub('Pak Dependent', 'PAK-Dependent', df$Reference)
    df$Reference <- gsub('Ras', 'RAS', df$Reference)

    # This is a temporary fix, this should be drawn automatically by the names of gwas.paths.list
    final_title <- gsub('_no_meta_analysis', '', plot_title)
    if(rachel) final_title <- str_extract(plot_title, '[^_]+$')
    if(final_title == 'AD') final_title <- "Alzheimer's Disease"
    if(final_title == 'BD') final_title <- "Bipolar Disorder"
    if(final_title == 'PD') final_title <- "Parkinson's Disease"
    if(final_title == 'SCZ') final_title <- "Schizophrenia"

    plot_dfs_dir <- paste0(all.plot.ggpath, "/plot_dataframes/")
    if(!dir.exists(plot_dfs_dir)) dir.create(plot_dfs_dir, recursive = T)
    fwrite(df, file = paste0(plot_dfs_dir, gsub('.pdf', '', ggname.pdf), '.csv'))

    # replace infinite values
    max.value = max(df$odds.ratio)

    message(ggname.pdf)
    message(all(is.na(df$neg_log10p)))
    if(!all(is.na(df$neg_log10p))){
      graph[[df_name]] <- ggplot(
        data = df,
        aes(
          x = Set,
          y = Reference,
          fill = neg_log10p
        )
      ) +
        geom_tile(color = "black", linewidth = 0.3) +
        scale_fill_gradientn(
          colors = c("#F5F5F5", "#F0E442", "#E69F00"),
          name   = bquote(-log[10](italic(P))),
          na.value = "gray50",
          values   = rescale(colorVals),
          breaks   = signif(as.vector(colorBreaks), digits = 2),
          guide    = "legend",
          limits   = colorLims
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid  = element_blank(),  # keep this so only tile borders act as grid
          plot.margin = unit(c(0, 0, 0, 0), units = "cm")
        ) +
        labs(title = df_name) +
        ylab("Pathway") +
        xlab("Cell-Type") +
        geom_text(aes(label = star), size = 10, fontface = "bold", color = "black") +
        coord_equal()

      if(FALSE){
        clusterOrder3_2 <- unique(df$Set)
        pathwaysPlot2 <- unique(df$Reference)

        for (a in 1:(length(clusterOrder3_2)+1)){
          graph[[df_name]] <- graph[[df_name]] + geom_segment(
            x=seq(0.5,length(clusterOrder3_2)+1,by=1)[a],
            xend=seq(0.5,length(clusterOrder3_2)+1,by=1)[a],
            y=0.5,
            yend=length(pathwaysPlot2)+0.5,
            color="black",
            linewidth=.5
          )
        }

        for (a in 1:(length(pathwaysPlot2)+1)) {
          graph[[df_name]] <- graph[[df_name]] + geom_segment(
            y=seq(0.5,length(pathwaysPlot2)+1,by=1)[a],
            yend=seq(0.5,length(pathwaysPlot2)+1,by=1)[a],
            x=0.5,
            xend=length(clusterOrder3_2)+0.5,
            color="black",
            linewidth=.5
          )
        }
      }

      if(!dir.exists(spec.ggpath)) dir.create(spec.ggpath, recursive = T)
      #ggsave(filename = ggname.png, plot = graph[[df_name]], path = spec.ggpath, width = 8, height = 8, units = "in")
      ggsave(filename = ggname.pdf, plot = graph[[df_name]], path = spec.ggpath, width = 8, height = 8, units = "in")
      ggsave(filename = ggname.png, plot = graph[[df_name]], path = spec.ggpath, width = 8, height = 8, units = "in")

      # keep identical return semantics: store the ggplot object (or NULL if skipped)
      # (this is effectively what lapply would have returned)
    } else {
      graph[[df_name]] <- NULL
    }
  }

  graphs_per_pvalue[[name_p]] <- graph
}


# STORE RESULTS IN GRID
names(graphs_per_pvalue) <- names(mylist)


require(gridExtra)
pbmcapply::pbmclapply(names(graphs_per_pvalue), function(name){
  graphs <- graphs_per_pvalue[[name]]
  
  all.traits.1.plot <- do.call('grid.arrange', graphs)
  
  
  ggsave(filename = paste0(name, ".pdf"), plot = all.traits.1.plot, path = all.plot.ggpath, width = 30, height = 30, units = "in")
  
}, mc.cores = parallel::detectCores() - 2)



