
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
      
      #require(ggplot2)
      
      
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
    Focus= F,
    alternative_tissues = FALSE){
  
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
      if(alternative_tissues) keys <- c("_subclass_", "_class_", "superClass_bulk", "Bulk")
      ### final list contains everything and it is readily availabe to be used as input in the GSEA Fisher pipeline.
      semifinal.list <- pbmcapply::pbmclapply(
        mydf, 
        mc.cores = parallel::detectCores() -2,
        function(df){
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
            if(!alternative_tissues){
              mylist[["All EUR Class"]] <- unname(unlist(mylist[grep("\\bClass", names(mylist), ignore.case = T)]))
              mylist[["All EUR Subclass"]] <- unname(unlist(mylist[grep("subclass", names(mylist), ignore.case = T)]))
              mylist[["All EUR Class + Subclass"]] <- c(mylist[["All EUR Class"]], mylist[["All EUR Subclass"]])
              mylist[["sn-Pseudohomogenate"]] <- unname(unlist(mylist[grep("Pseudobulk", names(mylist), ignore.case = T)]))
            } else{
              mylist[["All EUR Class"]] <- unname(unlist(mylist[grep("_class_", names(mylist), ignore.case = T)]))
              mylist[["All EUR Subclass"]] <- unname(unlist(mylist[grep("_subclass_", names(mylist), ignore.case = T)]))
              mylist[["All EUR Class + Subclass"]] <- c(mylist[["All EUR Class"]], mylist[["All EUR Subclass"]])
              mylist[["sn-Pseudohomogenate"]] <- unname(unlist(mylist[grep("superClass_bulk", names(mylist), ignore.case = T)]))
            }
            # this is to avoid the bulk overlap
            if(!alternative_tissues){
              mylist[["All EUR Subclass"]] <- c(mylist[["All EUR Subclass"]], mylist$`EUR Class-OPC`, mylist$`EUR Class-Astro`, mylist$`EUR Class-Endo`, mylist$`EUR Class-Oligo`)
            } else {
            mylist[["All EUR Subclass"]] <- c(mylist[["All EUR Subclass"]], 
                mylist[["230721_MegaAnalysis_EUR_class_OPC_prediXcan_noprior_alpha0.5_window1e6_filtered.db"]], 
                mylist[["230721_MegaAnalysis_EUR_class_Astro_prediXcan_noprior_alpha0.5_window1e6_filtered.db"]], 
                mylist[["230721_MegaAnalysis_EUR_class_Endo_prediXcan_noprior_alpha0.5_window1e6_filtered.db"]], 
                mylist[["230721_MegaAnalysis_EUR_class_Oligo_prediXcan_noprior_alpha0.5_window1e6_filtered.db"]])    
            }

            if(!Focus) mylist[["Bulk"]] <- unname(unlist(mylist[grep("Bulk", names(mylist), ignore.case = F)]))
            
            if(only.aggregates) mylist <- mylist[-c(1:(starting.length - 1))]
            
            return(mylist)
      })

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
      
      # parse pieces safely (fixes "0.05" being character)
      col <- as.character(if (is.list(input)) input[[1]] else input[1])
      thr <- as.numeric(if (is.list(input)) input[[2]] else input[2])

      # build output name with numeric threshold
      if (protein) {
        output <- paste(save.path, col, thr, "PROTEIN.RData", sep = "_")
      } else {
        output <- paste(save.path, col, thr, "MIXED.RData",   sep = "_")
      }

      # filter each df using data.table syntax
      mydf <- lapply(mydf, function(df) {
        df <- data.table::as.data.table(df)
        df[get(col) < thr]
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
          patterns <- grep(name, unique(df$tissue), value = TRUE)
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
        mylist[["All EUR Subclass"]] <- c(mylist[["All EUR Subclass"]], mylist$`EUR Class-OPC`, mylist$`EUR Class-Astro`, mylist$`EUR Class-Endo`, mylist$`EUR Class-Oligo`)
        if(!Focus) mylist[["Bulk"]] <- unname(unlist(mylist[grep("Bulk", names(mylist), ignore.case = F)]))
        
        if(only.aggregates) mylist <- mylist[-c(1:(starting.length - 1))]
        
        return(mylist)
      }, mc.cores = parallel::detectCores() -2)

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
