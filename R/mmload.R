#' mmload
#'
#' @param assembly 
#' @param coverage 
#' @param essential_genes 
#' @param kmer_pca 
#' @param kmer_BH_tSNE 
#' @param taxonomy 
#' @param additional 
#' @param verbose 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
mmload <- function(assembly,
                   coverage, 
                   essential_genes = NULL,
                   kmer_pca = FALSE,
                   kmer_BH_tSNE = FALSE,
                   taxonomy = NULL,
                   additional = NULL,
                   verbose = TRUE,
                   ...
) {
  ##### Internal functions #####
  #merge function
  mmmerge <- function(x, y, type) {
    #must be a data frame, named atomic vector, or a list of data frames and/or named vectors
    if(any(class(y) %in% c("list", "data.frame", "tbl", "tbl_df")) | is.atomic(y) | is.factor(y)) {
      #wrap non-lists in a list to work with the for-loop
      if(any(!class(y) %in% "list"))
        y <- list(y)
      #merge mm with each element in the provided list
      for(i in 1:length(y)) {
        string <- ifelse(length(y) == 1, paste0("'", type, "'"), paste0("'", type,"'", " element ", i))
        if(is.factor(y[[i]]))
          y[[i]] <- as.character(y[[i]])
        #first column must be the sequence names
        if(any(class(y[[i]]) %in% c("data.frame", "tbl", "tbl_df")) & length(y[[i]]) < 2)
          stop(paste0(string, " not accepted: Data frames must contain at least 2 columns where the first column contains the sequence names exactly matching those of the assembly."))
        
        #vectors must be named to be able to merge with mm
        if(is.atomic(y[[i]]) & is.null(names(y[[i]])))
          stop(paste0(string, " not accepted: The vector is not a named vector. The vector elements must be named by sequence names exactly matching those of the assembly."))
        
        #column names are preserved from data frames, but not from vectors. Use names of the provided list, or else create a dummy name
        if(is.atomic(y[[i]]) & !is.null(names(y[[i]]))) {
          y[[i]] <- tibble::enframe(y[[i]], name = "scaffold", value = ifelse((is.null(names(y)) | names(y)[[i]] == ""), paste0(type, i), names(y)[[i]]))
        }
        
        #merge x and y[[i]] by scaffold
        colnames(y[[i]])[1] <- "scaffold" #first columns must have same name
        y[[i]][1] <- lapply(y[[i]][1], as.character) #and must be character
        sharedScaffolds <- dplyr::intersect(x$scaffold, y[[i]][["scaffold"]]) #which scaffolds are shared between x and y[[i]]
        
        #print missing or excess scaffolds between x and y[[i]]
        if(!all(x$scaffold %in% y[[i]][["scaffold"]])) {
          missingScaffolds <- filter(x, !scaffold %in% sharedScaffolds)[[1]]
          warning(paste0("Only ", length(sharedScaffolds), " of all ", length(x$scaffold), " scaffolds in the assembly match in ", string,  ". The following ", length(missingScaffolds), " scaffolds are missing:\n\"", paste(missingScaffolds, collapse = "\", \""), "\""))
        } else if(!all(y[[i]][["scaffold"]] %in% x$scaffold)) {
          excessScaffolds <- filter(y, !scaffold %in% sharedScaffolds)[[1]]
          warning(paste0(string, " contains more scaffolds than the assembly. The following ", length(excessScaffolds), " scaffolds have not been loaded:\n\"", paste(excessScaffolds, collapse = "\", \""), "\""))
        } else if (!any(x$scaffold %in% y[[i]][["scaffold"]]))
          #no match sucks
          stop("No scaffold names match between the assembly and ", string, ". ")
        x <- dplyr::left_join(x, 
                              y[[i]], 
                              by = "scaffold")
      }
    } else
      stop("Data must be provided as a data frame, named vector, or a list of multiple data frames and/or named vectors.")
    return(x)
  }
  
  ##### Assembly #####
  #Load assembly sequences from the provided file path or object
  if(isTRUE(verbose))
    message("\nLoading assembly...")
  if(is.character(assembly)) {
    assembly <- Biostrings::readDNAStringSet(assembly, format = "fasta")
  } else if (!class(assembly) == "DNAStringSet") {
    stop("The assembly must either be an object of class \"DNAStringSet\" loaded with readDNAStringSet() from the Biostrings package, or a file path to the assembly FASTA file to be loaded.")
  }
  
  #check if a different assembly already exists in global environment
  MD5 <- digest::digest(assembly)
  if(!exists("assembly", where = .GlobalEnv, envir = .GlobalEnv)) {
    assign("assembly", assembly, envir = .GlobalEnv)
  } else if (exists("assembly", where = .GlobalEnv, envir = .GlobalEnv) & !identical(MD5, digest::digest(get("assembly", envir = .GlobalEnv)))) {
    userChoice <- readline(prompt = "An object named \"assembly\" already exists in the global environment. Do you want to overwrite it? (y/n): ")
    if(tolower(userChoice) == "n") {
      stop("Aborted by user.")
    } else if(tolower(userChoice) == "y") {
      assign("assembly", assembly, envir = .GlobalEnv)
    }
  }
  
  #attributes
  attributes <- list()
  attributes[["assemblyMD5"]] <- paste0(".mmID_", MD5)
  
  #duplicate sequence names are not allowed
  if(any(duplicated(names(assembly)))) 
    stop("The assembly contains duplicate sequence names")
  
  ##### Create base mm object #####
  if(isTRUE(verbose))
    message("Calculating GC content...")
  mm <- tibble::tibble(scaffold = names(assembly),
                       length = BiocGenerics::width(assembly),
                       gc_pct = round(as.numeric(Biostrings::letterFrequency(assembly, letters = c("CG"), as.prob=T))*100, digits = 2)
  )
  mm[1] <- lapply(mm[1], as.character)
  
  ##### Coverage #####
  if(isTRUE(verbose))
    message("Loading coverage data...")
  beforeMerge <- ncol(mm)
  mm <- mmmerge(x = mm,
                y = coverage,
                type = "coverage")
  attributes[["coverageCols"]] <- c((beforeMerge+1):ncol(mm))
  
  
  ##### Essential genes #####
  if(!is.null(essential_genes)) {
    if(is.data.frame(essential_genes) & ncol(essential_genes) > 1) {
      if(isTRUE(verbose))
        message("Loading essential genes...")
      attributes[["total_Ess.genes"]] <- nrow(essential_genes)
      if(ncol(essential_genes) > 2) {
        attributes[["unique_Ess.genes"]] <- length(unique(essential_genes$hmm.id))
      } else if(ncol(essential_genes) == 2) {
        attributes[["unique_Ess.genes"]] <- length(unique(essential_genes[,2]))
      }
      
      #replace all values in all character and factor columns to only contain alpha-numerics and dots "."
      if(any(sapply(essential_genes[,-1], class) %in% c("character", "factor"))) {
        essential_genes[,-1] <- dplyr::mutate_all(essential_genes[,-1], dplyr::funs(stringr::str_replace_all(., "[^[:alnum:].]", "")))
      }
      
      colnames(essential_genes)[1] <- "scaffold"
      essential_genes[1] <- lapply(essential_genes[1], as.character)
      essential_genes <- essential_genes %>% 
        dplyr::group_by(scaffold) %>% 
        dplyr::summarise_all(dplyr::funs(paste(., collapse = ", ")))
      
      mm <- dplyr::left_join(mm, 
                             essential_genes, 
                             by = "scaffold")
    } else if(!is.data.frame(essential_genes))
      stop("Essential genes must be provided as a data frame with at least 2 columns where the first column contains the sequence names exactly matching those of the assembly.")
  }
  
  ##### calculate tetranucleotides frequencies #####
  if(isTRUE(kmer_pca) | isTRUE(kmer_BH_tSNE)) {
    if(isTRUE(verbose))
      message("Calculating tetranucleotide frequencies in the assembly sequences...")
    kmer_fwd <- Biostrings::oligonucleotideFrequency(assembly, width = 4, as.prob = TRUE, with.labels = TRUE)
    kmer_revC <- Biostrings::oligonucleotideFrequency(Biostrings::reverseComplement(assembly), width = 4, as.prob = TRUE)
    kmer <- (kmer_fwd + kmer_revC)/2*100
  }
  
  ##### PCA of tetranucleotides #####
  if(isTRUE(kmer_pca)) {
    if(isTRUE(verbose))
      message("Calculating principal components of tetranucleotide frequencies...")
    PCA_res <- kmer %>%
      vegan::rda() %>%
      vegan::scores(choices = 1:3, display = "sites") %>%
      tibble::as.tibble()
    mm <- tibble::add_column(mm,
                             PC1 = PCA_res[[1]],
                             PC2 = PCA_res[[2]],
                             PC3 = PCA_res[[3]])
  }
  
  ##### BH tSNE of tetranucleotides #####
  if(isTRUE(kmer_BH_tSNE)){
    if(isTRUE(verbose))
      message("Calculating Barnes-Hut t-Distributed Stochastic Neighbor Embedding representations of tetranucleotide frequencies...")
    set.seed(42) # Sets seed for reproducibility
    tSNE_res <- Rtsne.multicore::Rtsne.multicore(kmer, verbose = verbose, check_duplicates = F, ...)[["Y"]] %>%
      tibble::as.tibble()
    mm <- tibble::add_column(mm,
                             tSNE1 = tSNE_res[[1]],
                             tSNE2 = tSNE_res[[2]])
  }
  
  ##### Taxonomy #####
  if(!is.null(taxonomy)) {
    if(isTRUE(verbose))
      message("Loading taxonomy...")
    if(any(class(taxonomy) %in% c("data.frame", "tbl", "tbl_df")) | is.atomic(taxonomy)) {
      its_a_vector <- ifelse(is.atomic(taxonomy), TRUE, FALSE)
      replacements <- c(" <phylum>" = "",
                        "unclassified Bacteria" = "Unclassified Bacteria",
                        "Fibrobacteres/Acidobacteria group" = "Acidobacteria",
                        "Bacteroidetes/Chlorobi group" = "Bacteroidetes",
                        "delta/epsilon subdivisions" = "Deltaproteobacteria",
                        "Chlamydiae/Verrucomicrobia group" = "Verrucomicrobia")
      for(i in 1:length(replacements)) {
        taxonomy <- lapply(taxonomy, gsub, pattern = paste(names(replacements)[i]), replacement = replacements[i])
      }
      if(!isTRUE(its_a_vector)) {
        taxonomy <- tibble::as.tibble(taxonomy)
      } else if (isTRUE(its_a_vector))
        taxonomy <- unlist(taxonomy)
    }
    mm <- mmmerge(x = mm,
                  y = taxonomy,
                  type = "taxonomy")
  }
  
  ##### Additional data #####
  if(!is.null(additional)) {
    if(isTRUE(verbose))
      message("Loading additional data...")
    mm <- mmmerge(x = mm,
                  y = additional,
                  type = "additional")
  }
  
  ##### Return #####
  if(isTRUE(verbose))
    message("Done!")
  
  #add class mm and the assembly MD5 as id
  class(mm) <- append(class(mm), c("mm", attributes[["assemblyMD5"]]))
  
  #some functions drop attributes set with attr(), so the solution is a hidden object named by an MD5 of the assembly
  assign(attributes[["assemblyMD5"]], attributes, envir = .GlobalEnv) 
  return(mm)
}
