gene    gene_name       zscore  beta    p       var_g   pred_perf_r2    pred_perf_pval  pred_perf_qval  n_snps_used     n_snps_in_cov   n_snps_in_model chr     gene_start      neg_log10p      fdr     bonferroni       
neg_log10fdr    neg_log10bonferroni     beta_dir        beta_mag        
gene_name_noName        tissue  pIB_bonferroni  pIB_fdr neg_log10pIB_bonferroni neg_log10pIB_fdr        BPcum


gene	gene_name	zscore	beta	p	var_g	pred_perf_r2	pred_perf_pval	pred_perf_qval	n_snps_used	n_snps_in_cov	n_snps_in_model	chr	gene_start	neg_log10p	
fdr	bonferroni	neg_log10fdr	neg_log10bonferroni	beta_dir	beta_mag	gene_name_noName	tissue	pIB_bonferroni	pIB_fdr	neg_log10pIB_bonferroni	neg_log10pIB_fdr	
BPcum	TW_fdr	TW_bonferroni	neg_log10TW_fdr	neg_log10TW_bonferroni


fisherGsea_2_MA
function(
    
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