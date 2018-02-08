#' Title
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
#' @export
#' 
#' @import tibble
#' @import digest
#' @import Biostrings
#' @import dplyr
#' @import vegan
#' @import Rtsne.multicore
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
  if(!exists("assembly", where = .GlobalEnv, envir = .GlobalEnv)) {
    assign("assembly", assembly, envir = .GlobalEnv)
  } else if (exists("assembly", where = .GlobalEnv, envir = .GlobalEnv) & !identical(digest::digest(assembly), digest::digest(get("assembly", envir = .GlobalEnv)))) {
    userChoice <- readline(prompt = "A different object named \"assembly\" already exists in the global environment. Do you want to overwrite it? (y/n or ENTER/ESC): ")
    if(any(tolower(userChoice) %in% c("y", "yes", ""))) {
      assign("assembly", assembly, envir = .GlobalEnv)
    } else
      stop("Aborted by user.")
  }
  
  #duplicate sequence names are not allowed
  if(any(duplicated(names(assembly)))) 
    stop("The assembly contains duplicate sequence names")
  
  ##### Create base mm object #####
  if(isTRUE(verbose))
    message("Calculating GC content...")
  mm <- tibble::tibble(scaffold = names(assembly),
                       length = BiocGenerics::width(assembly),
                       gc = round(as.numeric(Biostrings::letterFrequency(assembly, letters = c("CG"), as.prob=T))*100, digits = 2)
  )
  mm[1] <- lapply(mm[1], as.character)
  
  ##### Coverage #####
  if(isTRUE(verbose))
    message("Loading coverage data...")
  beforeMerge <- ncol(mm)
  mm <- mmmerge(x = mm,
                y = coverage,
                type = "coverage")
  colnames(mm)[c((beforeMerge+1):ncol(mm))] <- paste0("cov_", colnames(mm)[c((beforeMerge+1):ncol(mm))])
  
  
  ##### Essential genes #####
  if(!is.null(essential_genes)) {
    if(is.data.frame(essential_genes) & ncol(essential_genes) == 2) {
      if(isTRUE(verbose))
        message("Loading essential genes...")
      
      essential_genes[,1] <- as.character(essential_genes[,1])
      essential_genes[,2] <- as.character(essential_genes[,2])
      
      #replace all values to only contain alpha-numerics and dots "."
      essential_genes[,-1] <- dplyr::mutate_all(essential_genes[,-1, drop = FALSE], dplyr::funs(stringr::str_replace_all(., "[^[:alnum:].]", "")))
      
      colnames(essential_genes) <- c("scaffold", "geneID")
      essential_genes <- essential_genes %>% 
        dplyr::group_by(scaffold) %>% 
        dplyr::summarise_all(dplyr::funs(paste(., collapse = ",")))
      
      mm <- dplyr::left_join(mm, 
                             essential_genes, 
                             by = "scaffold")
    } else
      stop("Essential genes must be provided as a 2 column data frame, where the first column contains the sequence names exactly matching those of the assembly, and the second column the gene names/IDs.")
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
  return(mm)
}
