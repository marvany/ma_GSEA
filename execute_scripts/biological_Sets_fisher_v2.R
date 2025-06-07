#############################################################################
# THIS SCRIPT IS THE FINAL PIPELINE, ITS PATH IS AT THE HOME/GLOBAL.SCRIPTS #
source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")

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
path_to_RUSH = paste0(getwd(), '/Resources/genes')

meta.vs.nometa = "no_meta"
gwas.paths.list <- list(SCZ = paste0(path_to_RUSH, '/RUSH_incl_no_meta_control_SCZ_PRS_QDA'),
                        AD = paste0(path_to_RUSH, '/RUSH_incl_no_meta_control_AD_PRS_QDA'),
                        #PD = paste0(path_to_RUSH, '/RUSH_incl_no_meta_control_SCZ_PRS_QDA'),
                        BD = paste0(path_to_RUSH, '/RUSH_incl_no_meta_control_BIP_PRS_QDA'))

gene.symbol.column = "gene_id"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "p_value" #"mssm_p_value",
z.score.column = "z_score"

# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))



########################
# FOR FOTIS
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/fotis_GSEA')
path_to_genes = paste0(getwd(), '/Resources/genes/')
gwas.paths.list = list(WBS = path_to_genes)


gene.symbol.column = "gene_name"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "pvalue" #"mssm_p_value",
z.score.column = "zscore"

# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))



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
source(paste0(path.to.scripts, "/source_camera.for.all.R"))

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


ma.prepare.for.camera_3 <- function(paths.list, 
                                    gene_symbol_col, 
                                    ensembl.ID.col, 
                                    p.value.col, 
                                    z.score.col,
                                    fdr.column = NA,
                                    bonferroni.column = NA,
                                    feature.column = NA,
                                    return_dataframe = TRUE,
                                    convert.to.ensembl) {
  #browser()
  # Load required libraries if not already loaded
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }
  if (!requireNamespace("pbmcapply", quietly = TRUE)) {
    install.packages("pbmcapply")
  }
  library(data.table)
  library(pbmcapply)
  names(gwas.paths.list)
  # For each gwas path, create a dataframe for all tissues ready to be used for camera
  mylist <- pbmcapply::pbmclapply(names(paths.list), function(name) {
    
    # Get paths for this GWAS
    file_path <- paths.list[[name]]
    
    # Import files using fread instead of read.dfs
    myfile <- fread(file_path)
    
    # Change names of the dfs that correspond to each tissue
    #names(myfile) <- gsub(".csv|csv.EUR_|csv.class_|_work_files_.*", "", basename(file_path))
    ########### SPASE TO LAPPLY()/////ASDGFASDF
    # Prepare data for camera, ready_files = list of ready dfs, each df = a tissue
    dt <- as.data.table(myfile)
    

   all.dfs <- rbindlist(pbmcapply::pbmclapply(unique(dt$tissue), function(this.tissue) {

      df <- dt[tissue == this.tissue]

      # Create required columns
      df$pvalue <- df[[p.value.col]]
      df$gene_name <- df[[gene_symbol_col]]
      
      # Calculate FDR if not provided
      if (!is.na(fdr.column)) {
        df$fdr <- df[[fdr.column]] 
      } else {
        df$fdr <- p.adjust(df$pvalue, method = 'BH')
      }
      
      # Calculate Bonferroni if not provided
      if (!is.na(bonferroni.column)) {
        df$bonferroni <- df[[bonferroni.column]] 
      } else {
        df$bonferroni <- p.adjust(df$pvalue, method = 'bonferroni')
      }
      
      # Convert gene symbols to Ensembl IDs
      if(is.na(feature.column)){
          if(!is.na(gene.symbol.column)){
            df$feature <- unname(unlist(sapply(df[[gene.symbol.column]], function(x) {
            x <- toupper(x)
            ensID <- gene.symbol.to.ensembl(x, conversion = "symbol->ENSEMBL", GV.table = FALSE)
            if (ensID == 'character(0)') return(NA)
            return(ensID)
          })))
        } else stop("Both gene and feature columns are NA, exiting")
      } else{
        df$feature <- df[[feature.column]]
      }
      
      df$gene <- df$feature
      
      # Set gene column to feature
      
      df$zscore <- df[[z.score.col]]
      df$gwas <- name
      df$model_ID <- this.tissue
      
      # Select only the columns we need
      df <- df[, c('gene', 'feature', 'zscore', 'pvalue', 'fdr', 'bonferroni', 'gwas', 'model_ID')]
      
      return(df)
    }, mc.cores = parallel::detectCores() - 2))
    

    return(all.dfs)
  }, mc.cores = parallel::detectCores() -2)
  
  names(mylist) <- names(paths.list)
  
  if (return_dataframe) {
    # Using dplyr's bind_rows for consistency with original function
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      install.packages("dplyr")
    }
    library(dplyr)
    mydf <- dplyr::bind_rows(mylist)
    return(mydf)
  } else {
    return(mylist)
  }
}


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

if(!dir.exists(paste0(p.export, '/Results'))) dir.create(paste0(p.export, "/Results"))

save(gwaslist, file = paste0(p.export, "/Results/gwaslist.RData"))



##########################################
############# CHECKPOINT 1 ###############

load(paste0(p.export, "/Results/gwaslist.RData"))
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")
source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
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
#mylist <- mylist[1]

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

mylist = filter_empty_traits(mylist)

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

source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")
source(paste0(path.to.scripts, "/source_camera.for.all.R"))

##############################################
################### FISHER ###################
twas.list = mylist #for NO_focus   #####all.TWAS.forGSEA.fisher #for FOCUS,
if(FALSE) rachels.twas.list = twas.list
if(FALSE) twas.list = rachels.twas.list

if(FALSE) sanan.twas.list = twas.list
if(FALSE) twas.list = sanan.twas.list


background = MultiWAS::greatGenes #for NO_focus ####list(Imputable.Genes.NO.BULK, NULL) #for FOCUS,   # please do not name the elements included in this list
gsea.path = paste0(p.export, "/Results/all_traits_V6") # all_combinations_TWAS_GSEA_fisher_V3 uses protein specific background for Bulk and has improved aesthetics
gsea.path = paste0(p.export,"/Results/fisher/MGIPhenotypes_v5_protein_coding")
Validation_Genesets = mySets$standardGeneSets[c("syngoAll", "msigdbSetsPruned")] # Validation.Gene.List
Validation_Genesets = list(MGIPhenotypes = MultiWAS::standardGeneSets$MGIPhenotypes)


defined.genesets =  NULL # NULL for generalised run, 
tissues_to_examine = NULL
show.title = T

p.all.plots <- paste0(gsea.path, "/All.Plots")
if(!dir.exists(p.all.plots)) dir.create(p.all.plots, recursive = T)


# Top level - pval_name
pbmcapply::pbmclapply(names(twas.list), function(pval_name) {
  # pvalue level
  pval_level <- twas.list[[pval_name]]
  this_gsea_path <- paste0(gsea.path, "/", pval_name)
  
  # Second level - trait_name
  pbmcapply::pbmclapply(names(pval_level), function(trait_name) {
    # trait level
    thistrait <- pval_level[[trait_name]]
    trait_gsea_path <- paste0(this_gsea_path, "/", trait_name)
    
    # Third level - geneset_name
    pbmcapply::pbmclapply(names(Validation_Genesets), function(geneset_name) {
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
    }, mc.cores = parallel::detectCores() - 2)
    
    return(NULL)  # Explicit return
  }, mc.cores = parallel::detectCores() - 2)
  
  return(NULL)  # Explicit return
}, mc.cores = parallel::detectCores() - 2)


######################################################
################## PREPARE FOR PLOT ################## 

# Now we add the names manually...
syngo_name_of_pathway <- names(mySets$standardGeneSets["syngoAll"]$syngoAll$sets)
msigdb_name_of_pathway <- names(mySets$standardGeneSets$msigdbSetsPruned$sets)
mgi_phenotypes <- names(MultiWAS::standardGeneSets$MGIPhenotypes$sets)

## You need to also include 0.01
mylist <- read.dfs.2(gsea.path)
str(mylist)
# this was not needed for Sanan's analysis; you only have 1 dataset
mylist <- mylist[[1]]

mylist <- access_list(mylist,
                      condition = "any(grepl('individ', names(thisitem)))",
                      executable_code = "return(thisitem[2])")



########
# SKIP #

mylist <- access_list(mylist,
                      condition = "any(grepl('Astro', names(thisitem)))",
                      executable_code = "if(any(grepl('ORplot', names(thisitem)))){
  thisitem <- thisitem[!grepl('ORplot', names(thisitem))]
  } else return(thisitem)"
)


mylist <- access_list(mylist,
                      condition = "any(grepl('tsv', names(thisitem)))",
                      executable_code = "
                    names(thisitem) <- gsub('syngoAll__|.genesets.tsv|msigdbSetsPruned__|MGIPhenotypes_', '', names(thisitem))
                    return(thisitem)")


################################################
### MAKE A FUNCTION OUT OF THIS !!!!!!!!!!!!!!!!


############# CHANGE NAME FOR 'search_for_this' VARIABLE
# For syngoAll
if(FALSE) mylist <- change_depth(thislist = mylist, target_depth = find_depth(mylist, 'syngoAll') + 1, search_for_this = 'Astro') # add 1 because 1 layer is lost

# For MGIPhenotypes
if(TRUE) mylist <- change_depth(thislist = mylist, target_depth = find_depth(mylist, 'MGIPhen') + 1, search_for_this = 'Class.Astro') # add 1 because 1 layer is lost



save(mylist, file = paste0(gsea.path, '/temp_mylist.RData'))

if(TRUE){
      # merge MGIPhenotypes
      # At this point you merge msigdb and synGO
      mylist <- access_list(mylist,
                            condition = "any(grepl('MGIPheno', names(thisitem)))",
                            executable_code = "
                          thisitem <- do.call('rbind', thisitem)
                          return(thisitem)")

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
}


# choose top 10
if(FALSE){
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

}

mylist = mylist[[1]]
temp1 <- lapply(names(mylist), function(pval){
  
  thispval <- mylist[[pval]]

  temp2 <- lapply(thispval, function(enrich.dfs){
  
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
    return(enrich.dfs)
  })
  
  names(temp2) <- names(thispval)
  return(temp2)
})

names(temp1) <- names(mylist)
mylist <- temp1


df <- mylist$pvalue_0.05_backgrounds
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
all.plot.ggpath <- paste0(p.export, "/Results", "/all.plots")
all.plot.ggpath <- paste0(p.export, "/Results", "/all.plots_protein_only")

# when this lapply is ready you should include the heatmap thing and delete the df <- mylist.. above


# Open list with p_value thresholds
graphs_per_pvalue <- lapply(names(mylist), function(name_p){
  
  thislist <- mylist[[name_p]] # thislist contains all the traits-dfs
  
  # for conveniency
  p.value.thresh <- name_p
  
  # directory to store per p-value
  spec.ggpath <- paste0(all.plot.ggpath, "/", p.value.thresh)

  # define directory to store
  graph <- pbmcapply::pbmclapply(names(thislist), function(df_name){
    
    df <- thislist[[df_name]]
    
    
    
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
    
    df[df$intersection<2,]$neg_log10p <- NA
    
    
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
    
    df[df$Set == "Astro",]$Set <- "Astrocytes"
    df[df$Set == "EN",]$Set <- "Excitatory Neurons"
    df[df$Set == "Endo",]$Set <- "Endothelial"
    df[df$Set == "IN",]$Set <- "Inhibitory Neurons"
    df[df$Set == "Oligo",]$Set <- "Oligodendrocytes"
    
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
    final_title <- str_extract(plot_title, '[^_]+$')
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
    
    if(!all(is.na(df$neg_log10p))){
            graph <- ggplot(
            data = df,
            aes(
              x = Set,
              y = Reference,
              fill = neg_log10p
            )
          )+ geom_tile()+
            scale_fill_gradientn(colors=c("#F5F5F5","#F0E442","#E69F00"),
                                name=bquote(-log[10](italic(P))),na.value="gray50",values= rescale(colorVals), # rescale() rescales your data so that they linearly match a certain range of values
                                breaks=signif(as.vector(colorBreaks),digits=2),guide="legend",limits=colorLims)+  #rescale(colorVals)
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                  panel.grid = element_blank(), # Your plot was correct, but you needed to hide the panel grid!
                  plot.margin=unit(c(0,0,0,0),units="cm"))+
            labs(title = df_name)+
            ylab("Pathway")+
            xlab("Cell-Type")+
            geom_text(aes(label = star), size = 10, fontface = "bold", color = "black")+
            coord_equal()
          
          
          clusterOrder3_2 <- unique(df$Set)
          pathwaysPlot2 <- unique(df$Reference)
          
          
          for (a in 1:(length(clusterOrder3_2)+1)){
            graph <- graph + geom_segment(x=seq(0.5,length(clusterOrder3_2)+1,by=1)[a],
                                          xend=seq(0.5,length(clusterOrder3_2)+1,by=1)[a],
                                          y=0.5,
                                          yend=length(pathwaysPlot2)+0.5,
                                          color="black",
                                          linewidth=.5)
          }
          
          for (a in 1:(length(pathwaysPlot2)+1)) {
            graph <- graph + geom_segment(y=seq(0.5,length(pathwaysPlot2)+1,by=1)[a],yend=seq(0.5,length(pathwaysPlot2)+1,by=1)[a],
                                          x=0.5,xend=length(clusterOrder3_2)+0.5,color="black",linewidth=.5)
          }
          
          if(!dir.exists(spec.ggpath)) dir.create(spec.ggpath, recursive = T) 
          #ggsave(filename = ggname.png, plot = graph, path = spec.ggpath, width = 8, height = 8, units = "in")
          ggsave(filename = ggname.pdf, plot = graph, path = spec.ggpath, width = 8, height = 8, units = "in")
          ggsave(filename = ggname.png, plot = graph, path = spec.ggpath, width = 8, height = 8, units = "in")
          
          return(graph)
      }
    
  }, mc.cores = parallel::detectCores() - 2)
  names(graph) <- names(thislist)
  return(graph)
})


# STORE RESULTS IN GRID
names(graphs_per_pvalue) <- names(mylist)


require(gridExtra)
pbmcapply::pbmclapply(names(graphs_per_pvalue), function(name){
  graphs <- graphs_per_pvalue[[name]]
  
  all.traits.1.plot <- do.call('grid.arrange', graphs)
  
  
  ggsave(filename = paste0(name, ".pdf"), plot = all.traits.1.plot, path = all.plot.ggpath, width = 30, height = 30, units = "in")
  
}, mc.cores = parallel::detectCores() - 2)



