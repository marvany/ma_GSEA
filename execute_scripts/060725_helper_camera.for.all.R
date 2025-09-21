#MultiWAS::camera_for_all



# PREPARE FOR CAMERA
ma.prepare.for.camera_4 <- function(paths.list, 
                                    gene_symbol_col, 
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
 #   ready_files <- lapply(seq_along(myfiles), function(i){
  
      df <- myfiles[[i]]
      tissue <- names(myfiles)[i]
      df <- as.data.frame(df)
      df$pvalue <- df[[p.value.col]]
      df$gene_name <- df[[gene_symbol_col]] # RUN WITH GV TABLE FOR SYMBOLS
      if(!is.na(fdr_column)) df$fdr <- df[[fdr_column]] else df$fdr = p.adjust(df$pvalue, method = 'BH')
      if(!is.na(bonferroni_column)) df$bonferroni <- df[[bonferroni_column]] else df$bonferroni <- p.adjust(df$pvalue, method = 'bonferroni')
      
      setnames(df, 'gene_name', 'Approved.Symbol')
      dt = merge(df, geneAnnotation_v104_ensembl, all.x = TRUE, by = 'Approved.Symbol')
      
      setnames(dt, 'Approved.Symbol', 'gene_name')
      setnames(dt, 'Ensembl.Gene.ID', 'feature')
      df = dt
      # Check if there is anything else other than ENSG in your gene column
      #df$feature <- unname(unlist(sapply(df$gene_name, function(x){
      #  x <- toupper(x)
      #  ensID <- gene.symbol.to.ensembl(x, conversion = "symbol->ENSEMBL", GV.table = F)
      #  if(ensID == 'character(0)') return(NA)
      #  return(ensID)
      #})))
      
      # not sure if this should be gene_name
      df$gene <- df$feature
      df$zscore <- df[[z.score.col]]
      df$gwas <- name
      
      tissue <- names(myfiles)[i]
      df$model_ID <- tissue
      
      df = df[,c('gene', 'feature', 'zscore', 'pvalue', 'fdr', 'bonferroni', 'gwas', 'model_ID')]
      
      return(df)
   }, mc.cores = parallel::detectCores() - 2)
   #})
    
    # combine tissue dfs
    all.dfs <- do.call("rbind", ready_files)
  })
  
  names(mylist) = names(paths.list)
  
  if(return_dataframe) {
    mydf <- do.call("bind_rows", mylist)
    return(mydf)
  }else return(mylist)
}

# PREPARE FOR CAMERA
ma.prepare.for.camera_2 <- function(paths.list, 
                                    gene_symbol_col, 
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
  browser()
    # prepare data from camera, ready_files = list of ready dfs, each df = a tissue
 #ready_files <- pbmc
 lapply(seq_along(myfiles), function(i){
 #   ready_files <- lapply(seq_along(myfiles), function(i){
      browser()
      df <- myfiles[[i]]
      tissue <- names(myfiles)[i]
      df <- as.data.frame(df)
      df$pvalue <- df[[p.value.col]]
      df$gene_name <- df[[gene_symbol_col]] # RUN WITH GV TABLE FOR SYMBOLS
      if(!is.na(fdr_column)) df$fdr <- df[[fdr_column]] else df$fdr = p.adjust(df$pvalue, method = 'BH')
      if(!is.na(bonferroni_column)) df$bonferroni <- df[[bonferroni_column]] else df$bonferroni <- p.adjust(df$pvalue, method = 'bonferroni')
      
      # Check if there is anything else other than ENSG in your gene column
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
   }) #mc.cores = parallel::detectCores() - 2)
   #})
    
    # combine tissue dfs
    all.dfs <- do.call("rbind", ready_files)
  })
  
  names(mylist) = names(paths.list)
  
  if(return_dataframe) {
    mydf <- do.call("bind_rows", mylist)
    return(mydf)
  }else return(mylist)
}

######################################################


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
    
    z = cbind(z, geneMetaSets[[metaSet]]$metadata[match(z$Reference, 
                                                        geneMetaSets[[metaSet]]$metadata$name), ])
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



######################################################
#path.to.saved.images <- "/sc/arion/projects/va-biobank/marios/GSEA_Rachel/Saved.images"
#if(!dir.exists(path.to.saved.images)) dir.create(path.to.saved.images)

camera_for_all_MA_2 <- function(
    thisdf,
    genesetlist = mySets$standardGeneSets[c("syngoAll", "msigdbSetsPruned")],
    limit.analysis = "protein coding genes",
    limit.model = NA,
    input.geneset.format = "fisher",
    outdir = NULL,
    file.name
    ){
  get_enrich_order = function(res, inter.gene.cor = 0.05) {
    
    res$qvalue = marios.qvalue(res$pvalue)$qvalue # qvalue::qvalue(res$pvalue)$qvalue
    tstat = res$zscore
    names(tstat) = rownames(res)
    tstat = tstat[!is.na(names(tstat))]
    index = limma::ids2indices(geneSetsCombined, names(tstat)) # limma::ids2indices
    limma::cameraPR(tstat, index, inter.gene.cor = inter.gene.cor) # limma::cameraPR
  }
  
  camera_wrapper <- function(df, inter.gene.cor = 0.05) {
    
    x <- as.data.frame(df)
    x <- x[order(x$gene_name, -abs(x$zscore)), ]
    x = x[!duplicated(x$gene_name), ]
    x = x[!x$gene_name == "", ]
    rownames(x) = x$gene_name
    camera_deg <- get_enrich_order(x, inter.gene.cor = inter.gene.cor)
    return(camera_deg)
  }
  if (limit.analysis == "protein coding genes") {
    #    if ("model_ID" %!in% names(thisdf)) {
    #      message("Filtering for protein coding genes works by looking into the column 'model_ID' for the keyword 'genes' (case insensitive)")
    #    }
    thisdf <- thisdf[grep("^ENSG", thisdf$feature),]
  }
  protein.coding.reliably.imputed.ensembl <- unique(only_protein_coding(thisdf$gene, 
                                                                        input.type = "ENSEMBL", output.type = "ENSEMBL", return.only = "protein_coding"))
  protein.coding.reliably.imputed.symbol <- unique(only_protein_coding(thisdf$gene, 
                                                                       input.type = "ENSEMBL", output.type = "Gene Symbol", 
                                                                       return.only = "protein_coding"))
  if ((input.geneset.format == "camera") | purrr::vec_depth(genesetlist) == 
      2) {
    message("Geneset is camera compatible")
    genesetlist <- lapply(genesetlist, FUN = function(x) {
      x <- x[x %in% protein.coding.reliably.imputed.symbol]
      x <- x[unlist(lapply(x, FUN = function(y) length(y) > 
                             1))]
      return(x)
    })
  }
  else {
     if (input.geneset.format == "fisher" & purrr::vec_depth(genesetlist) == 4) {
      message("This is a fisher geneset, will convert to camera...")
      
      # The following function is incompatible with the format of your geneset list.
      genesetlist <- Gsea_gene_set_to_camera.ma(gsea.geneset.file = genesetlist, 
                                             mygenes = protein.coding.reliably.imputed.ensembl, 
                                             return.only.type = "protein_coding")
      # View(genesetlist)
    }
    else stop("I don't know how to parse this specific geneset, look into function camera_for_all() in R/GSEA_camera_2.R to implement handling additional geneset formats")
  }
  
  assign("geneSetsCombined", genesetlist, envir = globalenv())
  comb.grid.selection = unique(thisdf[, c("gwas", "model_ID")])
  
  output <- do.call(rbind, lapply(1:nrow(comb.grid.selection), FUN = function(i) {
    
    thisdf <- as.data.table(thisdf)
    df <- thisdf[gwas == as.character(comb.grid.selection$gwas[i]) & model_ID == as.character(comb.grid.selection$model_ID[i])]
    
    if (limit.analysis == "protein coding genes") {
      df <- df[gene %in% unique(only_protein_coding(thisdf$gene, 
                                                    input.type = "ENSEMBL", output.type = "ENSEMBL", 
                                                    return.only = "protein_coding"))]
    }
    setDF(thisdf)
    setDF(df)
    annotation <- unique(df[, c("gwas", "model_ID")])
    
    results <- camera_wrapper(df)     # this is the function created at the beginning of the function
    as.data.table(cbind(annotation, data.frame(pathway = row.names(results)), results))
  })) #, mc.cores = detectCores() - 1))
  if(!is.null(outdir)) if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
  #save(output, file = paste0(outdir, "/", file.name))
  if(!is.null(outdir)) fwrite(output, paste0(outdir, "/", file.name))
  return(output)
}

#

Gsea_gene_set_to_camera.ma <- function (gsea.geneset.file = mySets$standardGeneSets[c("gtex100", 
                                                                                   "descartesDiff", "syngoAll", "hpoPruned", "jaxPhenoPruned", 
                                                                                   "msigdbSetsPruned")], mygenes = NA, return.only.type = NA) 
{
  genes <- lapply(gsea.geneset.file, FUN = function(x) {
    desc <- x$metadata
    row.names(desc) <- desc$name

    # that's how he merges passes information from one dataframe to another..
    names(x$sets) <- desc[names(x$sets), "name_full"]
    x$metadata <- NULL
    x <- x$sets
    return(x)
  })
  genes <- unlist(genes, recursive = F)
  if (!is.na(fix_length_condition(mygenes))) 
    genes <- lapply(genes, function(x) {
      x[x %in% mygenes]
    })
  genes.symbols <- only_protein_coding(genes, input.type = "ENSEMBL", 
                                       output.type = "Gene Symbol", return.only = return.only.type)
  genes.symbols <- genes.symbols[unlist(lapply(genes.symbols, 
                                               FUN = function(x) length(x) > 2))]
  return(genes.symbols)
}

#





marios.qvalue <- function (p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE,
                           pi0 = NULL, ...) 
{
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  }
  else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
                                   1)) {
    save.image(paste0(path.to.saved.images, "/fdr.level.RData"))
    stop("'fdr.level' must be in (0, 1].")
  }
  if (is.null(pi0)) {
    pi0s <- marios.pi0est(p, ...) # pi0est(p, ...)
  }
  else {
    if (pi0 > 0 && pi0 <= 1) {
      pi0s = list()
      pi0s$pi0 = pi0
    }
    else {
      stop("pi0 is not (0,1]")
    }
  }
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m/(i * (1 - 
                                                        (1 - p[o])^m))))[ro]
  }
  else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m/i))[ro]
  }
  qvals_out[rm_na] <- qvals
  if (lfdr.out) {
    lfdr <- qvalue::lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  }
  else {
    lfdr_out <- NULL
  }
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out, 
                   pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level, 
                   significant = (qvals <= fdr.level), pi0.lambda = pi0s$pi0.lambda, 
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out, 
                   pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda, 
                   lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}








marios.pi0est <- function (p, lambda = 0, pi0.method = c("smoother", "bootstrap"), smooth.df = 3, smooth.log.pi0 = FALSE, ...){ 
# function (p, lambda = seq(0.05, 0.95, 0.05), pi0.method = c("smoother", "bootstrap"), smooth.df = 3, smooth.log.pi0 = FALSE, ...) 

  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method = match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda)
  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  }
  else if (ll > 1 && ll < 4) {
    stop(sprintf(paste("ERROR:", paste("length(lambda)=", 
                                       ll, ".", sep = ""), "If length of lambda greater than 1,", 
                       "you need at least 4 values.")))
  }
  else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  if (ll == 1) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  }
  else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec = lambda))[ind])/(length(p) * 
                                                                   (1 - lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      }
      else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W/(m^2 * (1 - lambda)^2)) * (1 - W/m) + (pi0 - 
                                                         minpi0)^2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    }
    else {
      stop("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
    }
  }
  if (pi0 <= 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
  }
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda, lambda = lambda, 
              pi0.smooth = pi0Smooth))
}



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
