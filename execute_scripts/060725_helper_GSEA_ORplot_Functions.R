#This is the list of functions used to do the GSEA using the Fisher exact test and to draw the OR plot for the results.
#It contains the following:
# 1. multi_tissues_fisherGSEA (a function that iterates across the validation genesets and each tissue)
# 2. multi.tissues_multi.backgrounds_fisherGSEA (a special form of multi_tissues_fisherGSEA function, that iterates over multiple background genesets)
# 2. fisherGSEA_2_MA (the function that performs the Fisher exact test)
# 3. create_tissue_paths (this function creates a list of the directory paths to draw the plots from. You need multiple paths that each corresponds to a specific trait across multiple tissues.)
# 4. multi_tissue_trait_OR_sig (this function creates the plot)

####################################################################################
####################################################################################
####################################################################################
####################################################################################

# 1. This function iterates over each validation geneset and over each tissue.
# This way each validation geneset is tested against the TWAS results of each tissue.
multi_tissues_fisherGSEA <- function(
    gsea.outdir,
    TWAS.genes.ID,
    validation.geneset){
  for(i in 1:length(validation.geneset)){ # we will iterate over the validation.geneset list to apply the enrichment list by list,
    for(tissue in 1:length(TWAS.genes.ID[[1]])){              #otherwise we get overlapping results.
      
      name.of.testGenes <- deparse(substitute(TWAS.genes.ID))
      
      # make a tissue specific directory path (eg. /EUR_Bulk) 
      tissue.outDir <- paste0(gsea.outdir,"/", names(TWAS.genes.ID[[1]])[tissue]) 
      
      # store the name of the tissue, this will be assigned as name_full and will be used in the OR plot
      name.of.tissue <- names(TWAS.genes.ID[[1]])[tissue]
      
      # This list will carry tissue-specific genes for each trait.
      tissue.genes <- list()
      
      for(trait in names(TWAS.genes.ID)){  
        
        # tissue.genes will be a list that for each trait (AD, BD, SCZ...), there will be a 'tissue-specific' vector of genes
        tissue.genes[[trait]] <- TWAS.genes.ID[[trait]][[tissue]]
      }
      fisherGsea_2_MA(testGenes = tissue.genes, 
                      geneMetaSets = validation.geneset[i], 
                      myGenes = MultiWAS::greatGenes,
                      outDir = tissue.outDir,
                      tissueName = name.of.tissue,
                      name.of.TWAS.genes = name.of.testGenes)
    }
  }
}


####################################################################################
####################################################################################
####################################################################################
####################################################################################


multi.tissues_multi.backgrounds_fisherGSEA <- function(
    gsea.outdir,
    TWAS.genes.ID,
    validation.geneset,
    background.genes){
  for(i in 1:length(validation.geneset)){ # we will iterate over the validation.geneset list to apply the enrichment list by list,
    for(tissue in 1:length(TWAS.genes.ID[[1]])){              #otherwise we get overlapping results.
      
      name.of.testGenes <- deparse(substitute(TWAS.genes.ID))
      # make a tissue specific directory path (eg. /EUR_Bulk) 
      tissue.outDir <- paste0(gsea.outdir,"/", names(TWAS.genes.ID[[1]])[tissue]) 

      # store the name of the tissue, this will be assigned as name_full and will be used in the OR plot
      name.of.tissue <- names(TWAS.genes.ID[[1]])[tissue]
      
      # This list will carry tissue-specific genes for each trait.
      tissue.genes <- list()
      #if(a > 1) 
      
      for(trait in names(TWAS.genes.ID)){  
        
        # tissue.genes will be a list that for each trait (AD, BD, SCZ...), there will be a 'tissue-specific' vector of genes
        tissue.genes[[trait]] <- TWAS.genes.ID[[trait]][[tissue]]
      }
      # Match the name of the tissue with the tissue specific background genes
      # background.genes is expected to be a list of vectors, where each vector has the genes that will be 
      # used as background. 
      # The names of the vectors should 100% match the names of tissues

      
      # length > 1 implies multiple backgrounds, if it is NA then tissue.spec.background gets NA and then loads myGreatGenes
      if(length(background.genes) > 1){
        tissue.spec.background <- background.genes[names(background.genes) %in% name.of.tissue]
        tissue.spec.background <- unlist(tissue.spec.background)
      }else if(is.na(background.genes)){
        tissue.spec.background <- NA
      }else tissue.spec.background <- background.genes[names(background.genes) %in% name.of.tissue]
       
      
      # This will throw an error if your background genes are incompatible with the names you have on your testgenes
      if(is.null(tissue.spec.background)) stop("tissue.spec.background is NULL, consider making sure that names of TWAS proteins and names of multiple background match")
      fisherGsea_2_MA(testGenes = tissue.genes, 
                      geneMetaSets = validation.geneset[i], 
                      myGenes = tissue.spec.background,          #MultiWAS::greatGenes,
                      outDir = tissue.outDir,
                      tissueName = name.of.tissue,
                      name.of.TWAS.genes = name.of.testGenes)
    }
  }
}

#



MultiWAS::fisherGsea_2

###################################################################
###################################################################
###################################################################
###################################################################
###################################################################

# 2. This is the main function that performs the Fisher Exact test.
fisherGsea_2_MA <- function(
    
  testGenes,       # Genes that will be tested for enrichment
  geneMetaSets,   # Enrichment set.
  myGenes, # Genes used as reference
  
  outDir,
  
  tissueName,
  name.of.TWAS.genes,
  
  provideORCI95 = T, # change this to TRUE!
  multipleTissues = T, # created by Marios
  # The below are constant for the most part.
  
  noFigs = T,
  onlySig = T,
  onlyIndividual = F, 
  addExtension = ".genesets",
  forceNoOutDirCheck = T,
  returnData = F,
  cores = 1,
  lambda = 0){
  
  list.of.packages <- c("reshape2", "fdrtool", "parallel", 
                        "vegan", "qvalue")
  suppressMessages(lapply(list.of.packages, require, character.only = TRUE))
  require("GeneOverlap")
  require(parallel)
  
  # (This gets skipped) Modify the testGenes into a list format. 
  if (!is.list(testGenes))testGenes = list(myTestGenes = testGenes)
  
  # No duplicates in testGenes
  testGenes = lapply(testGenes, unique)

  if(is.na(myGenes)[1]) myGenes = MultiWAS::greatGenes
  # This gets skipped (myGenes is defined)
  if (is.na(myGenes)[1]) {
    myGenes = read.delim(myGenesFallback, stringsAsFactors = F, 
                         col.names = c("chr", "tss", "strand", "ensembl"))
    myGenes = myGenes[!myGenes$chr %in% gseaExcludeChr, ]
    myGenes = myGenes$ensembl
  }
  
  # Filter testGenes through myGenes
  testGenes = sapply(names(testGenes), function(x) testGenes[[x]][testGenes[[x]] %in% 
                                                                    myGenes], simplify = F) # this return a list
  
  # Filter geneMetaSets through myGenes
  for (y in names(geneMetaSets)) {
    geneMetaSets[[y]]$sets = sapply(names(geneMetaSets[[y]]$sets), 
                                    function(x) {
                                      geneMetaSets[[y]]$sets[[x]][geneMetaSets[[y]]$sets[[x]] %in% 
                                                                    myGenes]
                                    }, simplify = F)
  }
  
  # Filter out "empty Sets" (doesn't do anything in present example)
  for (y in names(geneMetaSets)) {
    geneMetaSets[[y]]$sets = geneMetaSets[[y]]$sets[lapply(geneMetaSets[[y]]$sets, 
                                                           length) > 0]
  }
  
  # Analysis
  myGseas = list()
  
  
  # metaSets are the 'Sources' you used for your genes (eg. Psychiatric, Neurologic + Psychiatric) Alternatively it would be OMIM vs. ClinVar.
  for (metaSet in names(geneMetaSets)){
    print(metaSet)
    myGseas[[metaSet]] = list()
    
    # The following gets skipped
    if (!is.null(geneMetaSets[[metaSet]]$eligible_genes)) { 
      myGenes.backup <- myGenes
      myGenes <- geneMetaSets[[metaSet]]$eligible_genes
      testGenes.backup <- testGenes
      testGenes = sapply(names(testGenes), function(x) testGenes[[x]][testGenes[[x]] %in% 
                                                                        myGenes], simplify = F)
    }
    
    
    # gom object is the object in which the GeneSet Enrichment analysis is run.
    gom.obj = GeneOverlap::newGOM(testGenes, geneMetaSets[[metaSet]]$sets, 
                                  length(myGenes))
    
    # This gets skipped.
    if (provideORCI95) {
      OR.table <- expand.grid(names(geneMetaSets[[metaSet]]$sets), 
                              names(testGenes))
      names(OR.table) <- c("Reference", "Set")
      OR.table$OR.ci95.down <- numeric(nrow(OR.table))
      OR.table$OR.ci95.up <- numeric(nrow(OR.table))
      for (submetaset in names(geneMetaSets[[metaSet]]$sets)) {
        for (testgene in names(testGenes)) {
          SE.log.OR = sqrt(sum(1/as.vector(gom.obj@go.nested.list[[submetaset]][[testgene]]@cont.tbl)))
          OR.table[OR.table$Reference == submetaset & 
                     OR.table$Set == testgene, "OR.ci95.down"] <- exp(log(gom.obj@go.nested.list[[submetaset]][[testgene]]@odds.ratio) - 
                                                                        1.96 * SE.log.OR)
          OR.table[OR.table$Reference == submetaset & 
                     OR.table$Set == testgene, "OR.ci95.up"] <- exp(log(gom.obj@go.nested.list[[submetaset]][[testgene]]@odds.ratio) + 
                                                                      1.96 * SE.log.OR)
        }
      }
    }
    
    
    # Create the z object. I think z extracts information from the gom object.
    # This probably gets skipped considering that testGenes will be a list of multiple vectors.
    if (length(names(testGenes)) == 1) {
      z = cbind(names(testGenes), names(gom.obj@gsetB), 
                melt(getMatrix(gom.obj, name = "pval")), melt(getMatrix(gom.obj, 
                                                                        name = "odds.ratio"))[, 1], melt(getMatrix(gom.obj, 
                                                                                                                   name = "intersection"))[, 1], melt(getMatrix(gom.obj, 
                                                                                                                                                                name = "union"))[, 1], melt(getMatrix(gom.obj, 
                                                                                                                                                                                                      name = "Jaccard"))[, 1])
      
    }else{
      z = cbind(melt(getMatrix(gom.obj, name = "pval")), 
                melt(getMatrix(gom.obj, name = "odds.ratio"))[, 
                                                              3], melt(getMatrix(gom.obj, name = "intersection"))[, 
                                                                                                                  3], melt(getMatrix(gom.obj, name = "union"))[, 
                                                                                                                                                               3], melt(getMatrix(gom.obj, name = "Jaccard"))[, 
                                                                                                                                                                                                              3])
    }
    colnames(z) = c("Set", "Reference", "pval", "odds.ratio", 
                    "intersection", "union", "Jaccard")
    
    if (provideORCI95){
      require(dplyr)
      z <- inner_join(z, OR.table)
    }
    
    #The following add meaningful variables to the z object.
    z$Bonf_AdjP = p.adjust(z$pval, method = "bonferroni")
    z$BH_AdjP = p.adjust(z$pval, method = "BH")
    z$Z = qnorm(1 - (z$pval)/2)
    z$LogP = -log10(z$pval)
    z$FDR_AdjP = fdrtool(z$pval, statistic = "pvalue", plot = F, 
                         cutoff.method = "fndr", verbose = F)$qval
    z$FDR_AdjP = p.adjust(z$pval, method = "fdr")
 # apo edw mexri to mc lapply to fdr alla3e


    # MA: z$group is used by wrap() in OR.plot. By assigning the validation geneset, we produce forest-plot where the OR plots are per geneset,
    # and thus we can make comparisons of enrichments between the genesets.
    z$group <- metaSet
    if(multipleTissues) z$name_full <- tissueName
    
    
    z = z[order(z$pval),]
    
    # in the original format: match(z$reference, geneMetaSets[[metaSet]]$metadata$name), ]) 
    # the 'comma' led to 'uncorrect dimentions error' ; probably because geneMetaSets is a list, I erased the comma
    
    geneMetaSets[[metaSet]]$metadata = as.data.table(geneMetaSets[[metaSet]]$metadata)
    z = cbind(z, geneMetaSets[[metaSet]]$metadata[match(z$Reference, 
                                                        geneMetaSets[[metaSet]]$metadata$name)])

    z$name = NULL
    
    z$TWASgenes <- name.of.TWAS.genes
    myGseas[[metaSet]] = z[order(z$pval),]
    
    # This is skipped.
    if (!is.null(geneMetaSets[[metaSet]]$eligible_genes)){
      myGenes <- myGenes.backup
      testGenes <- testGenes.backup
    }
  } # End of the for() loop.
  
  
  
  # This is run in simple_GSEA_fisher
  if (onlySig) z <- z[z$pval <= 0.05, ]
  #View(z)
  
  # This is run in simple_GSEA_fisher
  
  if (!is.null(outDir)) {
    
    # If the directory exists stop the function.
    if (!forceNoOutDirCheck & file.exists(outDir)) stop("the provided output dir already exists. please provide a non-existant dir. Script aborted.")
    
    # If the directory doesn't exist, create a directory. This step is done.
    if (!file.exists(outDir)) dir.create(outDir, recursive = T, showWarnings = F)
    
    # This creates a 'result_textFiles' directory'.
    if (!onlyIndividual) {
      tsvDir = paste0(outDir, "/result_textFiles")
      if (!file.exists(tsvDir)) 
        dir.create(tsvDir, showWarnings = F)
    }
    
    # This creates 'result_textFiles_individual
    tsvIndiDir = paste0(outDir, "/result_textFiles_individual")
    if (!file.exists(tsvIndiDir)) dir.create(tsvIndiDir, showWarnings = F)
    
    
    niceColumns = c("Set", "Reference", "pval", "odds.ratio", "intersection", 
                    "union", "FDR_AdjP", "Jaccard", "Bonf_AdjP", "BH_AdjP", 
                    "name_full", "group")
    
    # MA: in myGseas the supposed "name_full" column, was under the name "name".
    #names(myGseas[[metaSet]])[15] <- "name_full"
    
    if (provideORCI95) niceColumns <- c(niceColumns, "OR.ci95.down", "OR.ci95.up")
    
    
    # This creates error because mtsv is not found.
    mclapply(names(myGseas), function(metaSet) {
      if (!onlyIndividual) {
        
        MultiWAS:::mtsv(myGseas[[metaSet]][, niceColumns], tsvDir, # View(myGseas[["protein.coding.neuro.psych"]][, niceColumns])
                        gz = FALSE, fileBaseName = make.names(metaSet))
      }
      lapply(unique(myGseas[[metaSet]]$Set), function(mySet) MultiWAS:::mtsv(myGseas[[metaSet]][myGseas[[metaSet]]$Set == 
                                                                                                  mySet, niceColumns], tsvIndiDir, gz = F, fileBaseName = make.names(paste0(metaSet, 
                                                                                                                                                                            "__", mySet, addExtension))))
    }, mc.cores = detectCores())
    
    
    ##############################################################################################
    # This is skipped in simple_GSEA_fisher.
    # In other words if you have reached at this point, you can straightforwardly export your data.
    if (!noFigs) {
      z = list()
      z$aggGsea = sapply(myGseas, function(x) {
        x$regDomEnrichment = x$odds.ratio
        x$peaksFullName = x$Set
        x$name = x$Reference
        x$binomTest = x$pval
        x$plotLabel = ""
        x$plotLabel[x$pval < 0.05] = "<U+00B7>"
        x$plotLabel[x$BH_AdjP < 0.05] = "#"
        x
      }, simplify = F)
      
      
      outDir0 = outDir
      # Changing outDir to the /heatmaps outDir (I am not sure I like this step)
      outDir = paste0(outDir, "/heatmaps")
      message(paste("Note: setting mc.cores to", cores))
      options(mc.cores = cores)
      
      # Create the heatmaps directory
      if (!file.exists(outDir)) 
        dir.create(outDir, showWarnings = F, recursive = T)
      
      print("About to mcmapply")
      aggGseaPlotter <- MultiWAS:::aggGseaPlotter
      mcmapply(function(i) {
        if (length(unique(z$aggGsea[[i]]$name)) < 101) {
          mySetName = z$aggGsea[[i]]$geneMetaSets[1]
          subDir = paste0(outDir, "/", mySetName)
          if (!file.exists(subDir)) 
            dir.create(subDir, showWarnings = F, recursive = T)
          aggGseaPlotter(df = z$aggGsea[[i]], annoName = i, 
                         outDir = subDir)
          aggGseaPlotter(df = z$aggGsea[[i]], annoName = i, 
                         outDir = subDir, doCluster = T, geneSets = geneMetaSets[[i]]$sets)
        }
      }, names(z$aggGsea)) # names(z$aggGsea) 
      
      #GV: heatmaps: Top something best pathways from each peakSet #NOTE: could be parallized, but take care with z$topAggGsea
      z$topAggGsea = list()
      for (i in names(z$aggGsea)) {
        w = z$aggGsea[[i]]
        w = w[order(w$binomTest), ]
        for (myTop in c(5, 10, 15, 20, 25)) {
          myName = paste0(i, "_top", myTop)
          if (length(unique(w$name)) >= myTop) {
            myPathways = unique(unlist(lapply(unique(w$peaksFullName), 
                                              function(x) {
                                                w$name_full[w$peaksFullName == x][1:myTop]
                                              })))
            z$topAggGsea[[myName]] = w[w$name_full %in% 
                                         myPathways, ]
            mySetName = z$topAggGsea[[myName]]$geneMetaSets[1]
            subDir = paste0(outDir, "/", mySetName)
            if (!file.exists(subDir)) 
              dir.create(subDir, showWarnings = F, recursive = T)
            aggGseaPlotter(df = z$topAggGsea[[myName]], 
                           annoName = myName, outDir = subDir)
            aggGseaPlotter(df = z$topAggGsea[[myName]], 
                           annoName = myName, outDir = subDir, doCluster = T, 
                           geneSets = geneMetaSets[[i]]$sets)
            if (length(unique(z$topAggGsea[[myName]]$peaksFullName)) > 
                1) 
              plotGseaBiclust(df = z$topAggGsea[[myName]], 
                              annoName = myName, outDir = subDir)
          }
        }
      }
      rm(i, myTop, w, myPathways, myName)
      outDir = outDir0
      barplotDir = paste0(outDir, "/barplots")
      lapply(names(z$aggGsea), function(metaSet) singleGseaPlotter(df = z$aggGsea[[metaSet]], 
                                                                   annoName = metaSet, barplotDir))
    } ###################################################################################
    # End of the !noFigs if statement. (this is skipped in simple_GSEA_fisher)
    # Export the data.
    if (!onlyIndividual){ 
      save(list = ls(all = T), file = paste0(outDir, "/gsea.Rdata"), 
           envir = environment())
    }
  } # end of the outDir if() statement.
  if (returnData){ 
    return(myGseas)
  }
}


##################################################################################
##################################################################################
##################################################################################
##################################################################################

# 3. create_tissue_paths creates a list for paths for the analysis for each tissue.
# maybe create an option "in multiple.tissues = FALSE", in which case, each "geneset" list (eg. neuropsych) will contain a single vector.
create_tissue_paths <- function(
    gsea.outdir, 
    specific.tissue.dir = NULL, #example c("/EUR_Bulk", "/230719_MegaAnalysis_EUR_superClass_bulk", "/All EUR Class") #if this is NULL then OR plots will be drawn for all the tissues  
    tissues.to.trait.subdirs,     # these are the directories from tissue to trait; example "/result
    number_of_sets,   # number of validation genesets that were used (eg. "biological", "clinical")
    
    names.of.sets = NA # = c("example_name_1", "example_name_2")  #this will only work if the name.of.sets are of equal number with the number of sets.
    # they should be ordered according to their actual order of the exported tsv files.
){                             
  if(is.null(specific.tissue.dir)) tissue.directories <- list.files(gsea.outdir, full.names = T) # this returns the directories for all the tissues with the full path name.
  if(!is.null(specific.tissue.dir)) tissue.directories <- paste0(gsea.outdir, specific.tissue.dir) # this returns the specific tissue directories with the full path names.
  
  path.to.traits <-paste0(tissue.directories, tissues.to.trait.subdirs)
  
  #This is for debugging / error prevention.
  debugger <- 0
  for(ind in seq_along(path.to.traits)){
    if(!dir.exists(path.to.traits[ind])){
      cat(path.to.traits[ind], "does not exist. \n\n")
      debugger <- 1
    }
  }
  
  if (debugger == 1) stop("The paths seen above are wrong. Consider checking the tissues.to.trait.subdirs argument.")
  print("All paths to traits exist.")
  
  # all.traits is a vector for all the traits tsv file names
  all.traits <- c()
  for(ind in seq_along(path.to.traits)){
    
    files <- list.files(path.to.traits[ind])
    
    # if the path.to.traits directory has files in it and all.traits has not been already assigned then assign
    # the file names of the path.to.traits directory to all.files
    if(length(files) > 0 && length(all.traits) == 0){
      all.traits <- files
      cat(path.to.traits[ind], " was used to extract the names of the final directory names were the OR plots will be stored.")
      cat("The following files were used for the names: \n\n", paste(all.traits, collapse = "\n"))
    }else if(length(files) == 0){
      cat(path.to.traits[ind], " is empty.\n\n")
      print("Pausing for a bit...")
      Sys.sleep(5)
    }  
    
    if(length(all.traits) > 0 && any(all.traits != files)) cat(path.to.traits[ind], " is different from all.traits.\n\n")
  }
  
  # the number of the tsv files per geneset used, is the traits.number.per.set
  traits.number.per.set <- length(all.traits) / number_of_sets 
  
  # we create a vector of indices which mark where is the first trait.tsv file for each geneset used
  new.set.index <- seq(from = 1, to = length(all.traits), by = traits.number.per.set)
  
  
  # this returns a list of  sequences that each corresponds to a set of traits
  traits.seq.per.set.list <- lapply(new.set.index, function(starting.number){
    ending.number <- starting.number + (traits.number.per.set - 1)
    return(starting.number:ending.number)
  })
  if(!is.na(names.of.sets[1])) if(length(names.of.sets) == length(traits.seq.per.set.list)) names(traits.seq.per.set.list) = names.of.sets
  
  # now we change the sequence numbers to the actual names of the tsv files based on the "all.traits" vector
  trait.tsv.files.per.set.list <- lapply(traits.seq.per.set.list, function(seq.number) all.traits[seq.number])
  geneset.trait.path.list <- lapply(trait.tsv.files.per.set.list, function(traits){
    lapply(traits, function(single.trait) paste0(path.to.traits, "/", single.trait))
  })
  
  # create a list that has the names of each trait for each geneset.
  names.of.traits.list <- lapply(traits.seq.per.set.list, function(geneset){
    sapply(geneset, function(seq){
      all.traits[seq]               # all.traits contains all the traits for all the genesets
    })
  })
  
  for(i in 1: length(geneset.trait.path.list)){
    names(geneset.trait.path.list[[i]]) <- names.of.traits.list[[i]]
  }
  
  return(geneset.trait.path.list)
}


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#4. This is the plotting function.
multi_tissue_trait_OR_sig <- function(
    
  gsea.outdir,                        # in the context of the MultiWAS package, this will also contain the name of the validation geneset, 
  # thus there will be no overlap of the same trait across multiple validation genesets
  multi.traits.multi.paths,            # multiple.paths is a vector with all the paths for all tisues for a single trait
  
  all.plots.path = NULL,
  show.title = T,
  
  plot.per.trait = T, # this draws the OR plot for each trait separately
  plot.all.traits = T, # this draws a graph for all the traits
  title.of.plot = NULL,
  # limit.tos.sets = NA, # limit to specific sets
  min.intersection      = 2, # do not consider results when intersection is lower
  file.pattern          = NA, # file pattern to not take all files
  # file.pattern     ="B2\\.Respiratory\\.\\.lung|B2\\.Adipose\\.\\.visceral\\.\\.STARNET|B2\\.Blood\\.\\.STARNET|B2\\.GI\\.\\.pancreas\\.\\.GTEx|B2\\.Artery\\.\\.Mammary\\.\\.STARNET",
  # !!! removes the pathway even if the intersection is less than min.intersection in one of the GWASs !!!
  returnplot            = F,
  saveplot              = T,
  facet                 = T,
  clinical              = F, # MGI, OMIM, clinical pooled / better to update these datasets before running. # MA: You will be using the OMIM, so maybe focus on the clinical
  miRNA                 = T,  # this is good
  pLI                   = T,  # this is good
  SynGO                 = T,  # This is good for brain-associated traits
  brain.location        = F,  # Could work well but only for homogenate
  cell.type             = F,  # we usually do partitioned heritability for this
  pruned                = T,  # only keep pruned when there are multiple
  top                   = 10, # top pathways which will be plotted for all GWASs if it is NA, it shows qvalue significant.
  
  custom_parent_traits  = NA, # here you can give a custom list of parent traits        #        gsea.outdir is your input !!!!!     #
  plot.width            = 8.5,                                                          ##############################################
  plot.height           = 5,
  return.multi.enrich.dfs = F
) {
  ###################################################
  # OR plot like the Figure 5d in the EpiXcan paper #
  ###################################################
  
  # 03.30.2024; synthesize a final function and create the plots, with the other validation genesets and while accounting for the protein coding only areas
  #03.30.2024; figure out what the middle section where you have commented a lot of stuff out do
  
  # These will be the arguements
  # For the manual.OR.plot.2 we will be using the draft.gsea.paths
  #multi.traits.multi.paths <- draft.gsea.paths$neuro.psych
  #all.traits.1.graph <- T
  #gsea.outdir <- paste0(getwd(), "/Results/DRAFT_GSEA")
  ###### OR plot ######
  
  # We rename to all.files for conveniency.
  all.files <- multi.traits.multi.paths
  
  # Bcause the OR plots directory is stored in the Results directory, which also contains
  # the results of the Fisher exact test. That's why we have to remove from all.files the directory that corresponds to the ORplots.
  if(length(grep("ORplots", all.files)) > 0) all.files <- all.files[-grep("ORplots", all.files)] # this needs to be lapplied
  
  ## Prepare df                                                                 
  # library(data.table) # for fread                                                                            
  all.enrich.dfs <- lapply(all.files, function(all.genesets){
    lapply(all.genesets, function(all.traits){
      all.tibles <- lapply(all.traits, FUN = fread)
      contegent.tibles <- do.call("rbind", all.tibles)                            # do.call unlist the elements of each list and passes them as arguements
      return(contegent.tibles)
    })
  })                                
  
  ############## First we define the palletes
  
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")  
  # The palette with black:                                                                               
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # First we create the directories.
  
  # This if(plot.per.trait) will create a separate directory for all traits
  if(plot.per.trait){
    for(i in 1:length(all.enrich.dfs)){
      all.files <- multi.traits.multi.paths
      #i <- 1
      # Call the trait-specific table.
      enrich.dfs <- all.enrich.dfs[[i]]
      
      # name.of.trait to create the name-specific directory
      name.of.trait <- names(all.enrich.dfs)[i]
      
      # These are the trait-specific directory paths, they will be stored with writeLines in the appropriate directory
      all.files <- all.files[i]
      
      # This creates the specific directories for each trait
      gsea.outdir.or = paste0(gsea.outdir, "/ORplots/", name.of.trait)
      if (!dir.exists(gsea.outdir.or)) dir.create(gsea.outdir.or, recursive = TRUE)
      
      #print(class(all.files))
      #print(all.files)
      #print(class(unlist(all.files)))
      # This saves the paths used to access the tissues in a txt.
      writeLines(unlist(all.files), paste0(gsea.outdir.or, "/result.files.considered.txt"))
      
      enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = "BH") #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
      
      
      
      enrich.dfs$star  <- as.character("")                                          # he creates an empty column and assigns an empty string everywhere
      enrich.dfs$star  <- ifelse(enrich.dfs$qvalue <= 0.001, "***",                 # then he accordingly annotates '***' based on the level of significance OF THE q-value
                                 ifelse(enrich.dfs$qvalue <= 0.01, "**",            # maybe change this to the FDR ?
                                        ifelse(enrich.dfs$qvalue <= 0.05, "*", "")))
      
      
      # enrich.dfs <- enrich.dfs[order(enrich.dfs$qvalue, decreasing = F),]
      enrich.dfs <- enrich.dfs[order(enrich.dfs$pval, decreasing = F),]             # he organizes stuff based on pval order.. (maybe change this to FDR)
      '%!in%' <- function(x,y)!('%in%'(x,y))                                        # this new operator checks if elements in x are NOT present in Y
      
      # 03.30.2024 figure out if these should get removed.
      # low.intersections.to.remove <- enrich.dfs[intersection < min.intersection]$name_full
      # enrich.dfs <- enrich.dfs[name_full %!in% low.intersections.to.remove]
      
      enrich.dfs <- enrich.dfs[intersection >= min.intersection] # in the lapply universe that is good
      fwrite(enrich.dfs, paste0(gsea.outdir.or, "/all.enrichments.with.min.intersection.", min.intersection,".csv")) # this is to export the table, after it has been filtered for minimum intersection.
      
      if (is.numeric(top)) { # if top then get top values based on q value
        top.genesets <- unique(enrich.dfs$name_full)[1:top]                              # the top will effectively filter out the genesets that are not included in the top 10 
        
        enrich.dfs   <- enrich.dfs[name_full %in% top.genesets]                     
      } else { # if not top then use qvalue filtering                               # this gets activated if top = NA
        top.genesets <- unique(enrich.dfs[qvalue <= 0.05, ]$name_full)              # This assigns some names to the top.genesets.
        enrich.dfs   <- enrich.dfs[name_full %in% top.genesets] }   
      
      enrich.dfs <- enrich.dfs[order(enrich.dfs$odds.ratio, decreasing = T),]
      
      # 03.30.2024 this still needs reviewing
      # the following need reviewing, this is the fwrite stuff
      #custom.sets <- unique(enrich.dfs$Set) # this list will be used for the UpSet plot                      #UpSet plots are a visualization technique used to display intersecting sets and their overlaps, 
      #if (length(custom.sets) > 5) custom.sets <- custom.sets[1:5] # up to five are plotted                    # often utilized in genomic studies to explore the relationships between different gene sets 
      #fwrite(enrich.dfs, paste0(gsea.outdir.or, "/significantly_associated_traits_asplotted.csv"))             # or pathways identified in analyses like gene set enrichment analysis (GSEA).
      #enrich.dfs$name_full <- stringr::str_wrap(enrich.dfs$name_full, width = 40) # Wrap every 40 characters   # he makes the text smaller in width to fit it in the plot
      #enrich.dfs <- as.data.frame(enrich.dfs)
      #enrich.dfs$name_full <- factor(enrich.dfs$name_full, rev(unique(enrich.dfs$name_full))) # Label wrap only works if you do it like factors / ordering as well.
      
      # he does reverse ordering
      
      
      ###### MA: 03.29.24This is supposed to change the names of the "group" column
      # the question is does this fit here?? or does it go to the first lapply statement??
      # ask Georgio what the group is supposed to represent
      
      # I am not sure under what circumstances it helps and that's why I commented it out.
      # enrich.dfs$group <-  unlist(lapply( # no need to use multicore for this        
      #  seq(nrow(enrich.dfs)),                                                  # this is a way to created indices.
      #  FUN = function(i){
      
      ## helper function to find if there is a better name                      
      #    give_new_group_name <- function(x) {                                      # this function is used to sort of "standarize" the names used in groups
      
      # This is required because id x is NA then group.conv[x,] returns a vector with NAs that has length of nrow(group.conv)
      #     ifelse (is.na(x), NA, group.conv[x,]) }                                          
      
      ## Find names                                                                      
      #    ifelse(
      #      is.na(give_new_group_name(enrich.dfs$group[i])),
      #      ifelse(is.na(enrich.dfs$group[i]), enrich.dfs$Set[i], enrich.dfs$group[i]),
      #      as.character(give_new_group_name(enrich.dfs$group[i])) ) } ) )
      #enrich.dfs <- enrich.dfs[order(enrich.dfs$Set),]
      # It is not easy to keep the order of everything
      
      
      if (nrow(enrich.dfs)>0) { # if there are not genes, there is nothing to plot           
        
        
        # library(ggplot2, lib.loc = "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"), ; library(grid)
        #View(enrich.dfs)
        #### ACTUAL PLOT ####  # the error bars (by geom_errorbars) are the CI bars.
        p.OR.CADsig <- ggplot( # reference line is a dashed line corresponds to OR = 1 (by geom_hline)                                                      
          enrich.dfs,          # text labels (by geom_text) will be added in each point for the stars
          aes(                 # there will be a legend to explain the colors and shapes, at the top right of the plot.
            x=name_full,       
            y=odds.ratio,      # the plot will be using a classic theme (theme_classic) with a base font size of 10
            label=star)) +
          geom_hline(aes(yintercept = 1), size = .25, linetype ="dashed") +
          geom_errorbar(aes(ymax = OR.ci95.up, ymin = OR.ci95.down, group =Set), size = 0.5,
                        width = 0.25, color = "gray50", position = position_dodge(width = 0.65)) +
          geom_point(aes(
            color = Set,
            group = Set, # this is required so that the labels can be aligned ;)
            shape = Set),
            size = 2.5, position = position_dodge(width = 0.65)) +
          #geom_point(aes(color = Set), size = 2.5, position = position_dodge(width = 0.65)) +
          geom_text(
            aes(group = Set), # to align the labels
            vjust=0,
            position = position_dodge(width = 0.65)) +
          scale_color_manual(
            name = "Set",
            values=rep(c(cbbPalette, cbPalette[1]), 2)[1:length(unique(enrich.dfs$Set))],
            labels = unique(enrich.dfs$Set),
            breaks = unique(enrich.dfs$Set)) +
          scale_shape_manual(
            name = "Set",
            values=c(rep(16, 9), rep(17,9))[1:length(unique(enrich.dfs$Set))],
            labels = unique(enrich.dfs$Set),
            breaks = unique(enrich.dfs$Set)) +
          theme_classic(base_size = 10) +
          xlab("") +
          ylab(expression(paste("Enrichment odds ratio for pathways"))) +
          scale_x_discrete(limits = levels(enrich.dfs$name_full)) +
          scale_y_continuous(trans='log10') +
          theme(legend.position = "top") + # c(0.77, 0.27)) +
          guides(color=guide_legend(title = paste0("TWAS"),
                                    ncol=2, byrow=T, keywidth=0.1, keyheight=0.1,
                                    default.unit="inch"),
                 shape=guide_legend(title = paste0("TWAS"),
                                    ncol=2, byrow=T, keywidth=0.1, keyheight=0.1,
                                    default.unit="inch")) +
          coord_flip()+
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        #   if (facet) p.OR.CADsig <- p.OR.CADsig + facet_grid(cols = vars(group)) # end of difficulties
        
        # You are exporting (and thus creating) only ONE plot (for a single trait)
        if (saveplot == T) {
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.pdf"), p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.png"), p.OR.CADsig, width = plot.width, height = plot.height)
        }
        # MA 03.29.2024 ; this doesn't make sense in the context of "drawing individual plots"
        # if (returnplot == T) return(p.OR.CADsig)
      } else message(paste0("For ", name.of.trait, ": there is nothing to plot"))
    }
  } #end of for loop for per.trait.plot  
  
  
  # This if(plot.all.traits) will create the all.traits.1.graph directory
  if(plot.all.traits){
    
    # all.enrich.dfs.backup <- all.enrich.dfs
    gsea.outdir.or = paste0(gsea.outdir, "/ORplots/", "all.traits.1.graph")
    if (!dir.exists(gsea.outdir.or)) dir.create(gsea.outdir.or, recursive = TRUE)
    
    # This saves the paths used to access the tissues in a txt.
    writeLines(unlist(all.files), paste0(gsea.outdir.or, "/result.files.considered.txt"))
    
    all.enrich.dfs <- lapply(all.enrich.dfs, function(all.genesets){
      
      lapply(seq_along(all.genesets), function(i){
        
        enrich.dfs <- all.genesets[[i]]
        
        #enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = "BH") #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
        #enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = "BH") #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
        
        
        enrich.dfs$qvalue <- enrich.dfs$BH_AdjP
        
        
        enrich.dfs$star  <- as.character("")                                          # he creates an empty column and assigns an empty string everywhere
        enrich.dfs$star  <- ifelse(enrich.dfs$qvalue <= 0.001, "***",                 # then he accordingly annotates '***' based on the level of significance OF THE q-value
                                   ifelse(enrich.dfs$qvalue <= 0.01, "**",            # maybe change this to the FDR ?
                                          ifelse(enrich.dfs$qvalue <= 0.05, "*", "")))
        
        
        # enrich.dfs <- enrich.dfs[order(enrich.dfs$qvalue, decreasing = F),]
        enrich.dfs <- enrich.dfs[order(enrich.dfs$pval, decreasing = F),]             # he organizes stuff based on pval order.. (maybe change this to FDR)
        '%!in%' <- function(x,y)!('%in%'(x,y))                                        # this new operator checks if elements in x are NOT present in Y
        
        # 03.30.2024 figure out if these should get removed.
        # low.intersections.to.remove <- enrich.dfs[intersection < min.intersection]$name_full
        # enrich.dfs <- enrich.dfs[name_full %!in% low.intersections.to.remove]
        
        enrich.dfs <- enrich.dfs[intersection >= min.intersection]
        
      })
    })
    
    
    
    # lapply(seq_along(multi.enrich.dfs), function(i){
    #   geneset_name <- paste0("geneset_", i) # this name is kind of sloppy, but will do for now
    #   fwrite(multi.enrich.dfs[i], paste0("/all.enrichments.with.min.intersection.", min.intersection,"_for_", geneset_name, ".csv")) # this is to export the table, after it has been filtered for minimum intersection.
    #})
    
    multi.enrich.dfs <- lapply(all.enrich.dfs, function(all.genesets){
      
      lapply(seq_along(all.genesets), function(i){
        
        enrich.dfs <- all.genesets[[i]]
        
        if (is.numeric(top)) { # if top then get top values based on q value
          top.genesets <- unique(enrich.dfs$name_full)[1:top]                              # the top will effectively filter out the genesets that are not included in the top 10 
          
          enrich.dfs   <- enrich.dfs[name_full %in% top.genesets]                     
        } else { # if not top then use qvalue filtering                               # this gets activated if top = NA
          top.genesets <- unique(enrich.dfs[qvalue <= 0.05, ]$name_full)              # This assigns some names to the top.genesets.
          enrich.dfs   <- enrich.dfs[name_full %in% top.genesets] }   
        
        enrich.dfs <- enrich.dfs[order(enrich.dfs$odds.ratio, decreasing = T),]
        
        
        # 03.30.2024 this still needs reviewing
        # the following need reviewing, this is the fwrite stuff
        #custom.sets <- unique(enrich.dfs$Set) # this list will be used for the UpSet plot                      #UpSet plots are a visualization technique used to display intersecting sets and their overlaps, 
        #if (length(custom.sets) > 5) custom.sets <- custom.sets[1:5] # up to five are plotted                    # often utilized in genomic studies to explore the relationships between different gene sets 
        #fwrite(enrich.dfs, paste0(gsea.outdir.or, "/significantly_associated_traits_asplotted.csv"))             # or pathways identified in analyses like gene set enrichment analysis (GSEA).
        #enrich.dfs$name_full <- stringr::str_wrap(enrich.dfs$name_full, width = 40) # Wrap every 40 characters   # he makes the text smaller in width to fit it in the plot
        #enrich.dfs <- as.data.frame(enrich.dfs)
        #enrich.dfs$name_full <- factor(enrich.dfs$name_full, rev(unique(enrich.dfs$name_full))) # Label wrap only works if you do it like factors / ordering as well.
        
        
        # I am not sure under what circumstances it helps and that's why I commented it out.
        # enrich.dfs$group <-  unlist(lapply( # no need to use multicore for this        
        #  seq(nrow(enrich.dfs)),                                                  # this is a way to created indices.
        #  FUN = function(i){
        
        ## helper function to find if there is a better name                      
        #    give_new_group_name <- function(x) {                                      # this function is used to sort of "standarize" the names used in groups
        
        # This is required because id x is NA then group.conv[x,] returns a vector with NAs that has length of nrow(group.conv)
        #     ifelse (is.na(x), NA, group.conv[x,]) }                                          
        
        ## Find names                                                                      
        #    ifelse(
        #      is.na(give_new_group_name(enrich.dfs$group[i])),
        #      ifelse(is.na(enrich.dfs$group[i]), enrich.dfs$Set[i], enrich.dfs$group[i]),
        #      as.character(give_new_group_name(enrich.dfs$group[i])) ) } ) )
        #enrich.dfs <- enrich.dfs[order(enrich.dfs$Set),]
        # It is not easy to keep the order of everything
      })
    })
    
    
    multi.enrich.dfs <- lapply(all.enrich.dfs, function(all.genesets){
      do.call("rbind", all.genesets)
    })
    multi.enrich.dfs <- do.call("rbind", multi.enrich.dfs)
    
    # title.of.plot <- name.of.TWAS
    if(return.multi.enrich.dfs == T) return(multi.enrich.dfs)
    
    # multi.enrich.dfs.backup <- multi.enrich.dfs
    
    if(nrow(multi.enrich.dfs) > 0){
      
      # 062624 FOR THE FINAL PLOT ENTER BREAKPOINT AT THIS POINT AND LOG THE
      # /sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/finalPlot_hotfix.R

      browser()
      ######## SPECIFIC FOR SANNAN MODIFICATIONS###############
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "All EUR Class"] <- "EUR Class-Aggregate"
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "All EUR Subclass"] <- "EUR Subclass-Aggregate"
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "Bulk"] <- "EUR Bulk" 
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "sn-Pseudohomogenate"] <- "EUR Pseudobulk Bulk" 
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "All EUR Class + Subclass"] <- "EUR Class+Subclass-Aggregate"
      multi.enrich.dfs$name_full <- factor(multi.enrich.dfs$name_full,levels=c("EUR Subclass-Aggregate","EUR Class-Aggregate","EUR Class+Subclass-Aggregate","EUR Pseudobulk Bulk", "EUR Bulk"))
      
      multi.enrich.dfs$legend <- factor(multi.enrich.dfs$name_full,levels=c("EUR Subclass-Aggregate","EUR Class-Aggregate","EUR Class+Subclass-Aggregate","EUR Pseudobulk Bulk", "EUR Bulk"))
      multi.enrich.dfs$Set <- unlist(lapply(multi.enrich.dfs$Set, function(x) strsplit(x,split="_")[[1]][1]))
      multi.enrich.dfs$Set[multi.enrich.dfs$Set == "Alcoholism"] <- "AUD"
      multi.enrich.dfs$Set <- factor(multi.enrich.dfs$Set,levels=rev(c("BD","SCZ","AD","MS","AUD","ALS","PD","Anorexia","Migraines","MDD","ADHD","Insomnia")))
      
      # this is used in scale_color_manual in values = ... This is to more easily manipulate the colors and order in your plot.
      colors_to_include <- rep(c(cbbPalette, cbPalette[1]), 2)[1:(length(unique(multi.enrich.dfs$name_full))) + 1]
      
      all.traits.p.OR.CADsig <-  ggplot( # reference line is a dashed line corresponds to OR = 1 (by geom_hline)                                                      
        
        multi.enrich.dfs,          # text labels (by geom_text) will be added in each point for the stars
        aes(                      # there will be a legend to explain the colors and shapes, at the top right of the plot.
          x=Set,       # there will be multiple panels (by facet_grid) 
          y=odds.ratio, # the plot will be using a classic theme (theme_classic) with a base font size of 10
          label=star)) +
        geom_hline(aes(yintercept = 1), size = .25, linetype ="dashed") +
        geom_errorbar(aes(ymax = OR.ci95.up, ymin = OR.ci95.down, group = name_full, color = name_full), size = 0.5,
                      width = 0.25, position = position_dodge(width = 0.65)) +
        
        geom_point(aes(
          # color = Set,
          group = name_full, 
          # shape = Set),
          size = intersection), 
          position = position_dodge(width = 0.65)) +
        scale_size_continuous(range = c(0.5,2))+
        #geom_point(aes(color = Set), size = 2.5, position = position_dodge(width = 0.65)) +
        geom_text(
          aes(group = name_full, y = OR.ci95.up), # to align the labels
          vjust = 0.7,
          hjust = -0.5,
          position = position_dodge(width = 0.65)) +
        
        # The following is optional to paint the error bars according to the pallete
        scale_color_manual( # color in the middle of the errorbar (modifies ge)
          name = "Set", # name of the legend
          values = colors_to_include, # edw praktika enwnei thn cbbPalette me to prwto element tou cbPalette kai kanei repeat 2 fores to neo palette
          labels = rev(levels(multi.enrich.dfs$name_full)),                                              # sth synexeia xrhsimopoiei tous deiktes gia na vrei ta xrwmata pou tha valei
          breaks = rev(levels(multi.enrich.dfs$name_full))) +
        scale_shape_manual( # shape in middle errorbar
          name = "Set",  # define the name in the legend
          values=c(rep(16, 9), rep(17,9))[1:length(unique(multi.enrich.dfs$Set))],  # defines the various shapes
          labels = unique(multi.enrich.dfs$Set),  # specifies the labels in the legend
          breaks = unique(multi.enrich.dfs$Set)) + # 
        theme_classic(base_size = 10) +
        xlab("") +
        ylab(expression(paste("Enrichment odds ratio for pathways"))) +
        labs(size = "Intersection"
        )+
        # scale_x_discrete(limits = levels(multi.enrich.dfs$name_full)) +
        scale_y_continuous(trans='log10') +
        theme(legend.position = "right",
              legend.text = element_text(size = 8)) + # c(0.77, 0.27)) +
        guides(color=guide_legend(title = paste0("TWAS"),
                                  ncol=1, byrow=F, keywidth=0.1, keyheight=0.1,
                                  default.unit="inch"
        ),
        shape=guide_legend(title = paste0("TWAS"),
                           ncol=1, byrow=F, keywidth=0.1, keyheight=0.1,
                           default.unit="inch")) +
        coord_flip()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      if(length(unique(multi.enrich.dfs$group)) > 1) all.traits.p.OR.CADsig <- all.traits.p.OR.CADsig + facet_wrap(~group) # this should be included if you are doing comparisons accross validation genesets
      if(show.title) all.traits.p.OR.CADsig <- all.traits.p.OR.CADsig + labs(title = title.of.plot)
      #remove bracket!!
      #   list.files(getwd())
      #list.files(paste0(getwd(), "/Scripts/GSEA_ORplot_Functions.R"))
      if (saveplot == T) {
        if(!is.null(title.of.plot)) {
          file.name <- gsub(" ", "_", title.of.plot)
          ggsave(paste0(gsea.outdir.or, "/", file.name, ".pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
          ggsave(paste0(gsea.outdir.or, "/", file.name, ".png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
          # all.plots.path <- gsub(basename(gsea.path), "", gsea.path)
          if(!is.null(all.plots.path)) {
            ggsave(paste0(all.plots.path, "/", title.of.plot, ".pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats = FALSE)
            ggsave(paste0(all.plots.path, "/", title.of.plot, ".png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)}
        }else{
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
        }
      }
      if (returnplot == T) return(all.traits.p.OR.CADsig)        #This excludes the subsequent line, fix this with a list
      
    } #else message(paste0("For ", name.of.trait, ": there is nothing to plot"))
  }
}



#################################################################
#################################################################
# 4.5 update of previous function
#################################################################
#################################################################
#################################################################
#4. This is the plotting function.
multi_tissue_trait_OR_sig_v2 <- function(
    
  gsea.outdir,                        # in the context of the MultiWAS package, this will also contain the name of the validation geneset, 
  # thus there will be no overlap of the same trait across multiple validation genesets
  multi.traits.multi.paths,            # multiple.paths is a vector with all the paths for all tisues for a single trait
  
  all.plots.path = NULL,
  show.title = T,
  
  plot.per.trait = T, # this draws the OR plot for each trait separately
  plot.all.traits = T, # this draws a graph for all the traits
  title.of.plot = NULL,
  # limit.tos.sets = NA, # limit to specific sets
  min.intersection      = 2, # do not consider results when intersection is lower
  file.pattern          = NA, # file pattern to not take all files
  # file.pattern     ="B2\\.Respiratory\\.\\.lung|B2\\.Adipose\\.\\.visceral\\.\\.STARNET|B2\\.Blood\\.\\.STARNET|B2\\.GI\\.\\.pancreas\\.\\.GTEx|B2\\.Artery\\.\\.Mammary\\.\\.STARNET",
  # !!! removes the pathway even if the intersection is less than min.intersection in one of the GWASs !!!
  returnplot            = F,
  saveplot              = T,
  facet                 = T,
  clinical              = F, # MGI, OMIM, clinical pooled / better to update these datasets before running. # MA: You will be using the OMIM, so maybe focus on the clinical
  miRNA                 = T,  # this is good
  pLI                   = T,  # this is good
  SynGO                 = T,  # This is good for brain-associated traits
  brain.location        = F,  # Could work well but only for homogenate
  cell.type             = F,  # we usually do partitioned heritability for this
  pruned                = T,  # only keep pruned when there are multiple
  top                   = 10, # top pathways which will be plotted for all GWASs if it is NA, it shows qvalue significant.
  
  custom_parent_traits  = NA, # here you can give a custom list of parent traits        #        gsea.outdir is your input !!!!!     #
  plot.width            = 8.5,                                                          ##############################################
  plot.height           = 5,
  return.multi.enrich.dfs = F
) {
  ###################################################
  # OR plot like the Figure 5d in the EpiXcan paper #
  ###################################################
  
  # 03.30.2024; synthesize a final function and create the plots, with the other validation genesets and while accounting for the protein coding only areas
  #03.30.2024; figure out what the middle section where you have commented a lot of stuff out do
  
  # These will be the arguements
  # For the manual.OR.plot.2 we will be using the draft.gsea.paths
  #multi.traits.multi.paths <- draft.gsea.paths$neuro.psych
  #all.traits.1.graph <- T
  #gsea.outdir <- paste0(getwd(), "/Results/DRAFT_GSEA")
  ###### OR plot ######
  
  # We rename to all.files for conveniency.
  all.files <- multi.traits.multi.paths
  
  # Bcause the OR plots directory is stored in the Results directory, which also contains
  # the results of the Fisher exact test. That's why we have to remove from all.files the directory that corresponds to the ORplots.
  if(length(grep("ORplots", all.files)) > 0) all.files <- all.files[-grep("ORplots", all.files)] # this needs to be lapplied
  
  ## Prepare df                                                                 
  # library(data.table) # for fread                                                                            
  all.enrich.dfs <- lapply(all.files, function(all.genesets){
    lapply(all.genesets, function(all.traits){
      all.tibles <- lapply(all.traits, FUN = fread)
      contegent.tibles <- do.call("rbind", all.tibles)                            # do.call unlist the elements of each list and passes them as arguements
      return(contegent.tibles)
    })
  })                                
  
  ############## First we define the palletes
  
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")  
  # The palette with black:                                                                               
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # First we create the directories.
  
  # This if(plot.per.trait) will create a separate directory for all traits
  if(plot.per.trait){
    for(i in 1:length(all.enrich.dfs)){
      all.files <- multi.traits.multi.paths
      #i <- 1
      # Call the trait-specific table.
      enrich.dfs <- all.enrich.dfs[[i]]
      
      # name.of.trait to create the name-specific directory
      name.of.trait <- names(all.enrich.dfs)[i]
      
      # These are the trait-specific directory paths, they will be stored with writeLines in the appropriate directory
      all.files <- all.files[i]
      
      # This creates the specific directories for each trait
      gsea.outdir.or = paste0(gsea.outdir, "/ORplots/", name.of.trait)
      if (!dir.exists(gsea.outdir.or)) dir.create(gsea.outdir.or, recursive = TRUE)
      
      #print(class(all.files))
      #print(all.files)
      #print(class(unlist(all.files)))
      # This saves the paths used to access the tissues in a txt.
      writeLines(unlist(all.files), paste0(gsea.outdir.or, "/result.files.considered.txt"))
      
      enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = "BH") #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
      
      
      
      enrich.dfs$star  <- as.character("")                                          # he creates an empty column and assigns an empty string everywhere
      enrich.dfs$star  <- ifelse(enrich.dfs$qvalue <= 0.001, "***",                 # then he accordingly annotates '***' based on the level of significance OF THE q-value
                                 ifelse(enrich.dfs$qvalue <= 0.01, "**",            # maybe change this to the FDR ?
                                        ifelse(enrich.dfs$qvalue <= 0.05, "*", "")))
      
      
      # enrich.dfs <- enrich.dfs[order(enrich.dfs$qvalue, decreasing = F),]
      enrich.dfs <- enrich.dfs[order(enrich.dfs$pval, decreasing = F),]             # he organizes stuff based on pval order.. (maybe change this to FDR)
      '%!in%' <- function(x,y)!('%in%'(x,y))                                        # this new operator checks if elements in x are NOT present in Y
      
      # 03.30.2024 figure out if these should get removed.
      # low.intersections.to.remove <- enrich.dfs[intersection < min.intersection]$name_full
      # enrich.dfs <- enrich.dfs[name_full %!in% low.intersections.to.remove]
      
      enrich.dfs <- enrich.dfs[intersection >= min.intersection] # in the lapply universe that is good
      fwrite(enrich.dfs, paste0(gsea.outdir.or, "/all.enrichments.with.min.intersection.", min.intersection,".csv")) # this is to export the table, after it has been filtered for minimum intersection.
      
      if (is.numeric(top)) { # if top then get top values based on q value
        top.genesets <- unique(enrich.dfs$name_full)[1:top]                              # the top will effectively filter out the genesets that are not included in the top 10 
        
        enrich.dfs   <- enrich.dfs[name_full %in% top.genesets]                     
      } else { # if not top then use qvalue filtering                               # this gets activated if top = NA
        top.genesets <- unique(enrich.dfs[qvalue <= 0.05, ]$name_full)              # This assigns some names to the top.genesets.
        enrich.dfs   <- enrich.dfs[name_full %in% top.genesets] }   
      
      enrich.dfs <- enrich.dfs[order(enrich.dfs$odds.ratio, decreasing = T),]
      
      # 03.30.2024 this still needs reviewing
      # the following need reviewing, this is the fwrite stuff
      #custom.sets <- unique(enrich.dfs$Set) # this list will be used for the UpSet plot                      #UpSet plots are a visualization technique used to display intersecting sets and their overlaps, 
      #if (length(custom.sets) > 5) custom.sets <- custom.sets[1:5] # up to five are plotted                    # often utilized in genomic studies to explore the relationships between different gene sets 
      #fwrite(enrich.dfs, paste0(gsea.outdir.or, "/significantly_associated_traits_asplotted.csv"))             # or pathways identified in analyses like gene set enrichment analysis (GSEA).
      #enrich.dfs$name_full <- stringr::str_wrap(enrich.dfs$name_full, width = 40) # Wrap every 40 characters   # he makes the text smaller in width to fit it in the plot
      #enrich.dfs <- as.data.frame(enrich.dfs)
      #enrich.dfs$name_full <- factor(enrich.dfs$name_full, rev(unique(enrich.dfs$name_full))) # Label wrap only works if you do it like factors / ordering as well.
      
      # he does reverse ordering
      
      
      ###### MA: 03.29.24This is supposed to change the names of the "group" column
      # the question is does this fit here?? or does it go to the first lapply statement??
      # ask Georgio what the group is supposed to represent
      
      # I am not sure under what circumstances it helps and that's why I commented it out.
      # enrich.dfs$group <-  unlist(lapply( # no need to use multicore for this        
      #  seq(nrow(enrich.dfs)),                                                  # this is a way to created indices.
      #  FUN = function(i){
      
      ## helper function to find if there is a better name                      
      #    give_new_group_name <- function(x) {                                      # this function is used to sort of "standarize" the names used in groups
      
      # This is required because id x is NA then group.conv[x,] returns a vector with NAs that has length of nrow(group.conv)
      #     ifelse (is.na(x), NA, group.conv[x,]) }                                          
      
      ## Find names                                                                      
      #    ifelse(
      #      is.na(give_new_group_name(enrich.dfs$group[i])),
      #      ifelse(is.na(enrich.dfs$group[i]), enrich.dfs$Set[i], enrich.dfs$group[i]),
      #      as.character(give_new_group_name(enrich.dfs$group[i])) ) } ) )
      #enrich.dfs <- enrich.dfs[order(enrich.dfs$Set),]
      # It is not easy to keep the order of everything
      
      
      if (nrow(enrich.dfs)>0) { # if there are not genes, there is nothing to plot           
        
        
        # library(ggplot2, lib.loc = "/sc/arion/projects/roussp01a/sanan/Rlibs/231221_R_4.2.0_MultiWAS"), ; library(grid)
        #View(enrich.dfs)
        #### ACTUAL PLOT ####  # the error bars (by geom_errorbars) are the CI bars.
        p.OR.CADsig <- ggplot( # reference line is a dashed line corresponds to OR = 1 (by geom_hline)                                                      
          enrich.dfs,          # text labels (by geom_text) will be added in each point for the stars
          aes(                 # there will be a legend to explain the colors and shapes, at the top right of the plot.
            x=name_full,       
            y=odds.ratio,      # the plot will be using a classic theme (theme_classic) with a base font size of 10
            label=star)) +
          geom_hline(aes(yintercept = 1), size = .25, linetype ="dashed") +
          geom_errorbar(aes(ymax = OR.ci95.up, ymin = OR.ci95.down, group =Set), size = 0.5,
                        width = 0.25, color = "gray50", position = position_dodge(width = 0.65)) +
          geom_point(aes(
            color = Set,
            group = Set, # this is required so that the labels can be aligned ;)
            shape = Set),
            size = 2.5, position = position_dodge(width = 0.65)) +
          #geom_point(aes(color = Set), size = 2.5, position = position_dodge(width = 0.65)) +
          geom_text(
            aes(group = Set), # to align the labels
            vjust=0,
            position = position_dodge(width = 0.65)) +
          scale_color_manual(
            name = "Set",
            values=rep(c(cbbPalette, cbPalette[1]), 2)[1:length(unique(enrich.dfs$Set))],
            labels = unique(enrich.dfs$Set),
            breaks = unique(enrich.dfs$Set)) +
          scale_shape_manual(
            name = "Set",
            values=c(rep(16, 9), rep(17,9))[1:length(unique(enrich.dfs$Set))],
            labels = unique(enrich.dfs$Set),
            breaks = unique(enrich.dfs$Set)) +
          theme_classic(base_size = 10) +
          xlab("") +
          ylab(expression(paste("Enrichment odds ratio for pathways"))) +
          scale_x_discrete(limits = levels(enrich.dfs$name_full)) +
          scale_y_continuous(trans='log10') +
          theme(legend.position = "top") + # c(0.77, 0.27)) +
          guides(color=guide_legend(title = paste0("TWAS"),
                                    ncol=2, byrow=T, keywidth=0.1, keyheight=0.1,
                                    default.unit="inch"),
                 shape=guide_legend(title = paste0("TWAS"),
                                    ncol=2, byrow=T, keywidth=0.1, keyheight=0.1,
                                    default.unit="inch")) +
          coord_flip()+
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        #   if (facet) p.OR.CADsig <- p.OR.CADsig + facet_grid(cols = vars(group)) # end of difficulties
        
        # You are exporting (and thus creating) only ONE plot (for a single trait)
        if (saveplot == T) {
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.pdf"), p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.png"), p.OR.CADsig, width = plot.width, height = plot.height)
        }
        # MA 03.29.2024 ; this doesn't make sense in the context of "drawing individual plots"
        # if (returnplot == T) return(p.OR.CADsig)
      } else message(paste0("For ", name.of.trait, ": there is nothing to plot"))
    }
  } #end of for loop for per.trait.plot  
  
  
  # This if(plot.all.traits) will create the all.traits.1.graph directory
  if(plot.all.traits){
    
    # all.enrich.dfs.backup <- all.enrich.dfs
    gsea.outdir.or = paste0(gsea.outdir, "/ORplots/", "all.traits.1.graph")
    if (!dir.exists(gsea.outdir.or)) dir.create(gsea.outdir.or, recursive = TRUE)
    
    # This saves the paths used to access the tissues in a txt.
    writeLines(unlist(all.files), paste0(gsea.outdir.or, "/result.files.considered.txt"))
    
    all.enrich.dfs <- lapply(all.enrich.dfs, function(all.genesets){
      
      lapply(seq_along(all.genesets), function(i){
        
        enrich.dfs <- all.genesets[[i]]
        
        #enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = "BH") #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
        #enrich.dfs$qvalue <- p.adjust(enrich.dfs$pval, method = "BH") #qvalue::qvalue(enrich.dfs$pval)$qvalues  # this line effectively creates a a q-value and assigns it to the corresponding column of enriched.dfs
        
        
        enrich.dfs$qvalue <- enrich.dfs$BH_AdjP
        
        
        enrich.dfs$star  <- as.character("")                                          # he creates an empty column and assigns an empty string everywhere
        enrich.dfs$star  <- ifelse(enrich.dfs$qvalue <= 0.001, "***",                 # then he accordingly annotates '***' based on the level of significance OF THE q-value
                                   ifelse(enrich.dfs$qvalue <= 0.01, "**",            # maybe change this to the FDR ?
                                          ifelse(enrich.dfs$qvalue <= 0.05, "*", "")))
        
        
        # enrich.dfs <- enrich.dfs[order(enrich.dfs$qvalue, decreasing = F),]
        enrich.dfs <- enrich.dfs[order(enrich.dfs$pval, decreasing = F),]             # he organizes stuff based on pval order.. (maybe change this to FDR)
        '%!in%' <- function(x,y)!('%in%'(x,y))                                        # this new operator checks if elements in x are NOT present in Y
        
        # 03.30.2024 figure out if these should get removed.
        # low.intersections.to.remove <- enrich.dfs[intersection < min.intersection]$name_full
        # enrich.dfs <- enrich.dfs[name_full %!in% low.intersections.to.remove]
        
        enrich.dfs <- enrich.dfs[intersection >= min.intersection]
        
      })
    })
    
    
    
    # lapply(seq_along(multi.enrich.dfs), function(i){
    #   geneset_name <- paste0("geneset_", i) # this name is kind of sloppy, but will do for now
    #   fwrite(multi.enrich.dfs[i], paste0("/all.enrichments.with.min.intersection.", min.intersection,"_for_", geneset_name, ".csv")) # this is to export the table, after it has been filtered for minimum intersection.
    #})
    
    multi.enrich.dfs <- lapply(all.enrich.dfs, function(all.genesets){
      
      lapply(seq_along(all.genesets), function(i){
        
        enrich.dfs <- all.genesets[[i]]
        
        if (is.numeric(top)) { # if top then get top values based on q value
          top.genesets <- unique(enrich.dfs$name_full)[1:top]                              # the top will effectively filter out the genesets that are not included in the top 10 
          
          enrich.dfs   <- enrich.dfs[name_full %in% top.genesets]                     
        } else { # if not top then use qvalue filtering                               # this gets activated if top = NA
          top.genesets <- unique(enrich.dfs[qvalue <= 0.05, ]$name_full)              # This assigns some names to the top.genesets.
          enrich.dfs   <- enrich.dfs[name_full %in% top.genesets] }   
        
        enrich.dfs <- enrich.dfs[order(enrich.dfs$odds.ratio, decreasing = T),]
        
        
        # 03.30.2024 this still needs reviewing
        # the following need reviewing, this is the fwrite stuff
        #custom.sets <- unique(enrich.dfs$Set) # this list will be used for the UpSet plot                      #UpSet plots are a visualization technique used to display intersecting sets and their overlaps, 
        #if (length(custom.sets) > 5) custom.sets <- custom.sets[1:5] # up to five are plotted                    # often utilized in genomic studies to explore the relationships between different gene sets 
        #fwrite(enrich.dfs, paste0(gsea.outdir.or, "/significantly_associated_traits_asplotted.csv"))             # or pathways identified in analyses like gene set enrichment analysis (GSEA).
        #enrich.dfs$name_full <- stringr::str_wrap(enrich.dfs$name_full, width = 40) # Wrap every 40 characters   # he makes the text smaller in width to fit it in the plot
        #enrich.dfs <- as.data.frame(enrich.dfs)
        #enrich.dfs$name_full <- factor(enrich.dfs$name_full, rev(unique(enrich.dfs$name_full))) # Label wrap only works if you do it like factors / ordering as well.
        
        
        # I am not sure under what circumstances it helps and that's why I commented it out.
        # enrich.dfs$group <-  unlist(lapply( # no need to use multicore for this        
        #  seq(nrow(enrich.dfs)),                                                  # this is a way to created indices.
        #  FUN = function(i){
        
        ## helper function to find if there is a better name                      
        #    give_new_group_name <- function(x) {                                      # this function is used to sort of "standarize" the names used in groups
        
        # This is required because id x is NA then group.conv[x,] returns a vector with NAs that has length of nrow(group.conv)
        #     ifelse (is.na(x), NA, group.conv[x,]) }                                          
        
        ## Find names                                                                      
        #    ifelse(
        #      is.na(give_new_group_name(enrich.dfs$group[i])),
        #      ifelse(is.na(enrich.dfs$group[i]), enrich.dfs$Set[i], enrich.dfs$group[i]),
        #      as.character(give_new_group_name(enrich.dfs$group[i])) ) } ) )
        #enrich.dfs <- enrich.dfs[order(enrich.dfs$Set),]
        # It is not easy to keep the order of everything
      })
    })
    
    
    multi.enrich.dfs <- lapply(all.enrich.dfs, function(all.genesets){
      do.call("rbind", all.genesets)
    })
    multi.enrich.dfs <- do.call("rbind", multi.enrich.dfs)
    
    # title.of.plot <- name.of.TWAS
    if(return.multi.enrich.dfs == T) return(multi.enrich.dfs)
    
    # multi.enrich.dfs.backup <- multi.enrich.dfs
    
    if(nrow(multi.enrich.dfs) > 0){
      
      # 062624 FOR THE FINAL PLOT ENTER BREAKPOINT AT THIS POINT AND LOG THE
      # /sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/execute_scripts/finalPlot_hotfix.R
      
      
      ######## SPECIFIC FOR SANNAN MODIFICATIONS###############
      multi.enrich.dfs <- multi.enrich.dfs[!grepl("All EUR Class \\+ Subclass", multi.enrich.dfs$name_full),]

      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "All EUR Class"] <- "Class-Aggregate"
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "All EUR Subclass"] <- "Subclass-Aggregate"
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "Bulk"] <- "Bulk" 
      multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "sn-Pseudohomogenate"] <- "Pseudobulk Bulk" 
      multi.enrich.dfs$name_full <- factor(multi.enrich.dfs$name_full,levels=c("Subclass-Aggregate","Class-Aggregate","Pseudobulk Bulk", "Bulk"))

      multi.enrich.dfs$legend <- factor(multi.enrich.dfs$name_full,levels=c("Subclass-Aggregate","Class-Aggregate","Pseudobulk Bulk", "Bulk"))
      multi.enrich.dfs$Set <- unlist(lapply(multi.enrich.dfs$Set, function(x) strsplit(x,split="_")[[1]][1]))
      multi.enrich.dfs$Set[multi.enrich.dfs$Set == "Alcoholism"] <- "AUD"
      multi.enrich.dfs$Set <- factor(multi.enrich.dfs$Set,levels=rev(c("BD","SCZ","AD","MS","AUD","ALS","PD","Anorexia","Migraines","MDD","ADHD","Insomnia")))


      # this is used in scale_color_manual in values = ... This is to more easily manipulate the colors and order in your plot.
      colors_to_include <- rep(c(cbbPalette, cbPalette[1]), 2)[1:(length(unique(multi.enrich.dfs$name_full))) + 1]



      #multi.enrich.dfs <- df

      all.traits.p.OR.CADsig <-  ggplot( # reference line is a dashed line corresponds to OR = 1 (by geom_hline)                                                      
        
        multi.enrich.dfs,          # text labels (by geom_text) will be added in each point for the stars
        aes(                      # there will be a legend to explain the colors and shapes, at the top right of the plot.
          x=Set,       # there will be multiple panels (by facet_grid) 
          y=odds.ratio, # the plot will be using a classic theme (theme_classic) with a base font size of 10
          label=star)) +
        geom_hline(aes(yintercept = 1), size = .25, linetype ="dashed") +
        geom_errorbar(aes(ymax = OR.ci95.up, ymin = OR.ci95.down, group = name_full, color = name_full), size = 0.5,
                      width = 0.25, position = position_dodge(width = 0.65)) +
        
        geom_point(aes(
          # color = Set,
          group = name_full, 
          # shape = Set),
          size = intersection), 
          position = position_dodge(width = 0.65)) +
        scale_size_continuous(range = c(0.5,2))+
        #geom_point(aes(color = Set), size = 2.5, position = position_dodge(width = 0.65)) +
        geom_text(
          aes(group = name_full, y = OR.ci95.up), # to align the labels
          vjust = 0.7,
          hjust = -0.5,
          position = position_dodge(width = 0.65)) +
        # The following is optional to paint the error bars according to the pallete
        scale_color_manual( # color in the middle of the errorbar (modifies ge)
          name = "Set", # name of the legend
          values = colors_to_include, # edw praktika enwnei thn cbbPalette me to prwto element tou cbPalette kai kanei repeat 2 fores to neo palette
          labels = rev(levels(multi.enrich.dfs$name_full)),                                              # sth synexeia xrhsimopoiei tous deiktes gia na vrei ta xrwmata pou tha valei
          breaks = rev(levels(multi.enrich.dfs$name_full))) +
        scale_shape_manual( # shape in middle errorbar
          name = "Set",  # define the name in the legend
          values=c(rep(16, 9), rep(17,9))[1:length(unique(multi.enrich.dfs$Set))],  # defines the various shapes
          labels = unique(multi.enrich.dfs$Set),  # specifies the labels in the legend
          breaks = unique(multi.enrich.dfs$Set)) + # 
        theme_classic(base_size = 10) +
        xlab("")+
        ylab(expression(paste("Enrichment odds ratio for pathways"))) +
        labs(size = "Intersection"
        )+
        # scale_x_discrete(limits = levels(multi.enrich.dfs$name_full)) +
        scale_y_continuous(trans='log10',
                          expand = c(0.2, 0)) +
        theme(legend.position = c(0.18, 0.82),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12)) + # c(0.77, 0.27)) +
        guides(color=guide_legend(title = paste0("TWAS"),
                                  ncol=1, byrow=F, keywidth=0.1, keyheight=0.1,
                                  default.unit="inch"
        ),
        shape=guide_legend(title = paste0("TWAS"),
                          ncol=1, byrow=F, keywidth=0.1, keyheight=0.1,
                          default.unit="inch")) +
        coord_flip()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

      if(length(unique(multi.enrich.dfs$group)) > 1) all.traits.p.OR.CADsig <- all.traits.p.OR.CADsig + facet_wrap(~group) # this should be included if you are doing comparisons accross validation genesets
      if(show.title) all.traits.p.OR.CADsig <- all.traits.p.OR.CADsig + labs(title = title.of.plot)
      #remove bracket!!
      #   list.files(getwd())
      #list.files(paste0(getwd(), "/Scripts/GSEA_ORplot_Functions.R"))
      if (saveplot == T) {
        if(!is.null(title.of.plot)) {
          file.name <- gsub(" ", "_", title.of.plot)
          ggsave(paste0(gsea.outdir.or, "/", file.name, ".pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
          ggsave(paste0(gsea.outdir.or, "/", file.name, ".png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
          # all.plots.path <- gsub(basename(gsea.path), "", gsea.path)
          if(!is.null(all.plots.path)) {
            ggsave(paste0(all.plots.path, "/", title.of.plot, ".pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats = FALSE)
            ggsave(paste0(all.plots.path, "/", title.of.plot, ".png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
            }
        }else{
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
          ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
        }
      }
      if (returnplot == T) return(all.traits.p.OR.CADsig)        #This excludes the subsequent line, fix this with a list
      
    } #else message(paste0("For ", name.of.trait, ": there is nothing to plot"))
  }
}


#################################################################
#################################################################
#5. This is the function that 



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

#################


FOCUS_transformTIMnames2_update <- function(x) {
  ref <- fread("/sc/arion/projects/roussp01a/sanan/230706_NPSAD_filePrep/psychAD_functions/ref/FOCUS_TIMtransform2.tsv")
  ref$Normal[2] <- "EUR Class-Immune"
  ref$Normal[31] <- "EUR Subclass-Adaptive"
  
  index <- grep(T,ref$FOCUS == x)
  if (length(index) > 0) {
    return(ref$Normal[index])
  } else {
    index <- grep(T,ref$Normal == x)
    if (length(index) > 0) {
      return(ref$FOCUS[index])
    } else {
      return(NA)
    }
  }
}

###############################################################################################
########################### THESE ARE THE FUNCTIONS FOR CREATING FISHER COMPATIBLE GENESET LIST



################################################################
#############            Functions               ###############


# creates validation list of validation genesets
# input may be either a vector of genes or a list of multiple vectors of genes (the list should contain multiple vectors of genes!)
# if the original.geneset.list is left
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



# This is the function that creates an appropriate validation geneset list. For further explanation, check below.





# original.geneset.list is the list that contains a collection of lists, this function adds new genesets to that list
#if left NULL a new list will be created.
# genes.vector is a vector that contains the genes we want our new validation geneset to include.
# geneset.name is the name we want to assign to the new validation geneset.

update.geneset.list <- function(original.geneset.list = NULL, genes.vector, geneset.name){
  if(is.null(original.geneset.list))   original.geneset.list <- list()
  
  # proper.geneset.format() converts the genes.vector to a validation geneset with the appropriate format for our pipeline 
  geneset.to.add <- proper.geneset.format(genes.vector, geneset.name = geneset.name)
  
  # add_validation_geneset() adds a new validation geneset to the "list of validation genesets"
  updated.geneset.list <- add_validation_geneset(validation.geneset.list = original.geneset.list, geneset.to.add = geneset.to.add)
  return(updated.geneset.list)
}





# Inputs: 1. a vector of genes, 2. the name of the gene vector, 3. the name of the geneset
# Output a standard format of a list that is recognized by our pipeline.
# The output can be readily be used as input in the "add_to_validation_geneset function".

proper.geneset.format <- function(genes, geneset.name){ #genes must be a vector # names must be specified
  
  # This is the proper format for the GSEA pipeline
  proper.list = list(
    sets = list(name.to.replace = unname(genes)),
    metadata = list(name = geneset.name)
  )
  
  names(proper.list$sets) <- geneset.name
  return(proper.list) 
}


# Input is an original list of validation genesets + a geneset which has the "proper format" (of $sets / $metadata etc)
add_validation_geneset <- function(validation.geneset.list, geneset.to.add){
  
  # this adds a new sublist to the list
  validation.geneset.list[[length(validation.geneset.list) + 1]] <- geneset.to.add
  
  # this adds the name of the new geneset to the new geneset
  names(validation.geneset.list)[length(validation.geneset.list)] <- names(geneset.to.add$sets)
  return(validation.geneset.list)
}

#



