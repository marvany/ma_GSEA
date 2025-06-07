source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")
source(paste0(path.to.scripts, "/source_camera.for.all.R"))
library(data.table)

##############
# DEFINE PATHS
p.export = '/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/sanan_mgi_GSEA/Resources/genes/v2'
p.export = '/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/sanan_mgi_GSEA/Resources/genes/v3'

###############
# LOAD EXAMPLES
if(FALSE){
    p.exampleAD = '/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/rachel_GSEA/Resources/genes/RUSH_incl_no_meta_control_AD_PRS_QDA'
    p.exampleSCZ = '/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/rachel_GSEA/Resources/genes/RUSH_incl_no_meta_control_SCZ_PRS_QDA'
    p.ex.df = '/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/rachel_GSEA/Resources/genes/RUSH_incl_no_meta_control_SCZ_PRS_QDA/class_Astro_no_meta_analysis_work_files_list2.csv'

    list.files(p.exampleAD)
    list.files(p.exampleSCZ)

    ex.df = fread(p.ex.df)
    head(ex.df)
}


#############
# LOAD snTIMS
p.twas.dir = "/sc/arion/projects/roussp01a/sanan/230420_TWAS_companionPipelines/combinationAnalyses/stackTWAS/output/240406_allEUR/"
p.twas.dir = "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/250402_allEUR_revisions/"
p.all.twas = list.files(p.twas.dir, full.names = TRUE, recursive = T, pattern = '.tsv')
df = rbindlist(lapply(p.all.twas, function(p.thistwas){
    gwas = gsub('.*allEUR//|/acronyms.*', '', p.thistwas)
    thistwas = fread(p.thistwas)
    thistwas$gwas = gwas
    thistwas
}), fill = TRUE)

df = df[pIB_fdr < 0.05] # filter for FDR

setnames(df, old = c("p", "zscore", "gene_name_noName", "beta", 'tissue'), new = c("p_value", "z_score", "gene_id", "estimate", 'model_ID'))
df = df[, c("estimate", "p_value", "gene_id", "z_score", 'gwas', 'model_ID'), with = TRUE]

combo.grid = df[, .SD[1,], by = .(gwas, model_ID)][,.(gwas,model_ID)]


invisible(lapply(seq(nrow(combo.grid)), function(i){
    thisgwas = combo.grid[i,]$gwas
    thismodel = combo.grid[i,]$model_ID
    temp = df[gwas == thisgwas & model_ID == thismodel, c('estimate', 'p_value', 'gene_id', 'z_score'), with = TRUE]

    # Prepare for export
    thisdir = paste0(p.export, '/', thisgwas)
    if(!dir.exists(thisdir)) dir.create(thisdir, recursive = TRUE)
    fwrite(temp, file = paste0(thisdir, '/', gsub('-| ', '_', thismodel), '.csv'))
}))


###############################################
str(mylist)

names(gwaslist$SCZ)
nrow(gwaslist$SCZ)
nrow(mylist$pvalue_0.05$SCZ)

nrow(mylist$pvalue_0.05$SCZ)

unique(df$gwas)
colnames(gwaslist$SCZ)
temp2 = df[gwas == 'SCZ_2022_Trubetskoy_postImp']
all(temp2$p_value < 0.05)

unique(temp2$p_value)
temp2 = unlist(temp2$gene_id)
length(temp2)
length(totest)


test = rbindlist(mylist[['pvalue_0.05']]$SCZ)
totest = unlist(mylist[['pvalue_0.05']]$SCZ)
df[gwas == 'SCZ']
nrow(test)
View(mylist)



###########################################


##################
# FOR RACHEL FINAL
source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")
source(paste0(path.to.scripts, "/source_camera.for.all.R"))
library(ggplot2)

# FOR SANAN
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/sanan_mgi_GSEA')
path_to_snTIMs = paste0(getwd(), '/Resources/genes/v1')
full.paths = list.files(path_to_snTIMs, full.names = TRUE)
gwas.paths.list = lapply(full.paths, function(x) x)
names(gwas.paths.list) = gsub('_.*','' , basename(full.paths))

gene.symbol.column = "gene_id"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "p_value" #"mssm_p_value",
z.score.column = "z_score"

# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))



gwaslist <- ma.prepare.for.camera_2(paths.list = gwas.paths.list, # expects list as input
                                    gene_symbol_col = gene.symbol.column,
                                    p.value.col = p.value.column, 
                                    z.score.col = z.score.column,
                                    return_dataframe = F)





#########################
# function troubleshoot
                                    
ma.prepare.for.camera_2 <- function(paths.list, gene_symbol_col, 
                                    ensembl.ID.col, 
                                    p.value.col, 
                                    z.score.col,
                                    fdr_column = NA,
                                    bonferroni_column = NA,
                                    return_dataframe = T){
  # for each gwas path, the list has dataframe with all the tissues, which are ready to be used for camera.
  mylist <- lapply(names(paths.list), function(name){
                # take path x
                x <- paths.list[[name]]
                
                # Import files
                myfiles <- read.dfs(x)
                
                # modify list
                myfiles <- unlist(myfiles, recursive = F)
                # Change names of the dfs that correspond to each tissue
                
                names(myfiles) <-  gsub(".csv|csv.EUR_|csv.class_|_work_files_.*", "", names(myfiles))
                # prepare data from camera, ready_files = list of ready dfs, each df = a tissue
                ready_files <- pbmclapply(seq_along(myfiles), function(i){
                                df <- myfiles[[i]]
                                tissue <- names(myfiles)[i]
                                df <- as.data.frame(df)
                                df$pvalue <- df[[p.value.col]]
                                df$gene_name <- df[[gene_symbol_col]] # RUN WITH GV TABLE FOR SYMBOLS
                                if(!is.na(fdr_column)) df$fdr <- df[[fdr_column]] else df$fdr = p.adjust(df$pvalue, method = 'BH')
                                if(!is.na(bonferroni_column)) df$bonferroni <- df[[bonferroni_column]] else df$bonferroni <- p.adjust(df$pvalue, method = 'bonferroni')
                                            df$feature <- unname(unlist(sapply(df$gene_name, function(x){
                                                x <- toupper(x)
                                                ensID <- gene.symbol.to.ensembl(x, conversion = "symbol->ENSEMBL", GV.table = F)
                                                if(ensID == 'character(0)') return(NA)
                                                return(ensID)
                                })))
                                
                                # not sure if this should be gene_name
                                df$gene <- df$feature
                                df$zscore <- df[[z.score.col]]
                                df$gwas <- name
                                
                                tissue <- names(myfiles)[i]
                                df$model_ID <- tissue
                                
                                df = df[,c('gene', 'feature', 'zscore', 'pvalue', 'fdr', 'bonferroni', 'gwas', 'model_ID')]
                                
                                return(df)
                })
                
                        # combine tissue dfs
                all.dfs <- do.call("rbind", ready_files)
  })
  
  names(mylist) = names(paths.list)
  
  if(return_dataframe) {
    mydf <- do.call("bind_rows", mylist)
    return(mydf)
  }else return(mylist)
}









##################################################
# TROUBLE SHOOT

#######################
# FOR RACHEL FINAL
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/rachel_GSEA')
path_to_snTIMs = paste0(getwd(), '/Resources/genes')

meta.vs.nometa = "no_meta"
gwas.paths.list <- list(SCZ = paste0(path_to_snTIMs, '/RUSH_incl_no_meta_control_SCZ_PRS_QDA'),
                        AD = paste0(path_to_snTIMs, '/RUSH_incl_no_meta_control_AD_PRS_QDA'),
                        #PD = paste0(path_to_snTIMs, '/RUSH_incl_no_meta_control_SCZ_PRS_QDA'),
                        BD = paste0(path_to_snTIMs, '/RUSH_incl_no_meta_control_BIP_PRS_QDA'))

gene.symbol.column = "gene_id"  # MODIFY IF THE TABLE DOESN'T CONTAIN SYMBOL NAMES
p.value.column = "p_value" #"mssm_p_value",
z.score.column = "z_score"

# p_value is defined by user, the rest are default because you create them inside ma.prepare.for.camera_2
thresh_annot <- expand.grid(threshold = c("0.05", "0.01"), metric = c("pvalue", "fdr", "bonferroni"))
