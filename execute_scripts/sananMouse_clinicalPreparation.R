# IN THIS SCRIPT YOU ARE USING THE PREPARATION DONE
# in execute_scripts/clinical_genes_prepare.R

source("/hpc/users/anyfam01/Global.Scripts/Global.Source.R")
source("/sc/arion/projects/va-biobank/marios/Clinical.Significance/Scripts/GSEA_ORplot_Functions.R")
source(paste0(path.to.scripts, "/source_camera.for.all.R"))
library(ggplot2, lib.loc = "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS")
library(ggplot2)
library(scales)
library(pbmcapply)
setwd('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/sanan_mgi_GSEA')
base.outdir = '/sc/arion/projects/va-biobank/PROJECTS/ma_gsea'

# revisions
main.dir <- "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/250402_allEUR_revisions/" # contains the twas
outdir = 'revisions'
bulk.filter.coding = FALSE
background = MultiWAS::greatGenes #for NO_focus ####list(Imputable.Genes.NO.BULK, NULL) #for FOCUS,   # please do not name the elements included in this list

Validation_Genesets = mySets$standardGeneSets[c("syngoAll", "msigdbSetsPruned")] # Validation.Gene.List
Validation_Genesets = list(MGIPhenotypes = MultiWAS::standardGeneSets$MGIPhenotypes)

# DEPENDENT VARIABLES
lists.outdir = paste0(base.outdir, "/Resources/", outdir, '/data_lists')
gsea.path <- paste0(base.outdir, "/output/biological_pathways/", outdir)

if(!dir.exists(lists.outdir)) dir.create(lists.outdir, recursive = TRUE)
if(!dir.exists(gsea.path)) dir.create(gsea.path, recursive = TRUE)

# LOAD THE CREATED execute_scripts/clinical_genes_prepare.R
load(paste0(lists.outdir, "/ready_for_Fisher_analysis_TWAS_list.RData"))
twas.list = ready_for_Fisher_analysis_TWAS_list

defined.genesets =  NULL # NULL for generalised run, 
tissues_to_examine = NULL
show.title = T

p.all.plots <- paste0(gsea.path, "/All.Plots")
if(!dir.exists(p.all.plots)) dir.create(p.all.plots, recursive = T)

lapply(names(twas.list), function(pval_name){
  # pvalue level
  pval_level <- twas.list[[pval_name]]
  gsea.path <- paste0(gsea.path, "/", pval_name)
  
  lapply(names(pval_level), function(trait_name){
    # trait level
    thistrait <- pval_level[[trait_name]]
    gsea.path <- paste0(gsea.path, "/", trait_name)
    
    lapply(names(Validation_Genesets), function(geneset_name){
      Validation_Geneset <- Validation_Genesets[geneset_name]
      gsea.path <- paste0(gsea.path, "/", geneset_name)
      
      #if(is.null(defined.genesets)){
      #  if(any(grepl("PROTEIN", TWAS.name))){
      #    genesets <- c(4, 5, 7)
      #  }else{
      #    genesets <- c(1, 2, 3)
      #  }
      #}else genesets = defined.genesets
      # αυτο κανονικά είναι απο πάνω αλλα χρειαζεται debug γτ εχει το TWAS.name
      
      fisherGsea_2_MA(testGenes = thistrait, 
                      geneMetaSets = Validation_Geneset, 
                      myGenes = background,          #MultiWAS::greatGenes,
                      outDir = gsea.path,
                      tissueName = "sentinel",
                      name.of.TWAS.genes = trait_name)
      
    })
  })
})

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
if(TRUE) mylist <- change_depth(thislist = mylist, target_depth = find_depth(mylist, 'MGIPhen') + 1, search_for_this = 'Class_Astro') # add 1 because 1 layer is lost
#if(TRUE) mylist <- change_depth(thislist = mylist, target_depth = find_depth(mylist, 'MGIPhen') + 1, search_for_this = 'MGIPhenotypes__') # add 1 because 1 layer is lost

save(mylist, file = paste0(gsea.path, '/temp_mylist.RData'))
str(mylist)
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
all.plot.ggpath <- paste0(getwd(), "/Results", "/all.plots")
all.plot.ggpath <- paste0(getwd(), "/Results", "/all.plots_protein_only")

# when this lapply is ready you should include the heatmap thing and delete the df <- mylist.. above


# Open list with p_value thresholds
graphs_per_pvalue <- lapply(names(mylist), function(name_p){
  
  thislist <- mylist[[name_p]] # thislist contains all the traits-dfs
  
  # for conveniency
  p.value.thresh <- name_p
  
  # directory to store per p-value
  spec.ggpath <- paste0(all.plot.ggpath, "/", p.value.thresh)

  # define directory to store
  graph <- lapply(names(thislist), function(df_name){
    
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
            labs(title = final_title)+
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
    
  })
  names(graph) <- names(thislist)
  return(graph)
})

names(graphs_per_pvalue) <- names(mylist)


require(gridExtra)
lapply(names(graphs_per_pvalue), function(name){
  graphs <- graphs_per_pvalue[[name]]
  
  all.traits.1.plot <- do.call('grid.arrange', graphs)
  
  
  ggsave(filename = paste0(name, ".pdf"), plot = all.traits.1.plot, path = all.plot.ggpath, width = 30, height = 30, units = "in")
  
})



