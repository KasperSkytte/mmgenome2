#' @title Load and combine data for use with other mmgenome2 functions
#' 
#' @description Loads, validates and combines multiple aspects of metagenome data into one dataframe for use with all mmgenome2 functions, including scaffold assembly sequences, scaffold coverage, essential genes, taxonomy, and more.
#'
#' @param assembly (\emph{required}) A character string with the path to the assembly FASTA file, or the assembly as already loaded with \code{\link{readDNAStringSet}}.
#' @param coverage (\emph{required}) A \code{vector}, \code{dataframe}, or a \code{list} hereof containing coverage of each scaffold. The prefix \code{"cov_"} will be appended to all coverage column names in the output.
#' \describe{
#'   \item{\code{vector}}{If provided as a vector, the elements of the vector must be named by the scaffold names exactly matching those of the assembly.}
#'   \item{\code{dataframe}}{If provided as a dataframe, the first column must contain the scaffold names exactly matching those of the assembly, and any additional column(s) contain coverage of each scaffold.}
#'   \item{\code{list}}{If provided as a list, it must contain any number of \code{vector} or \code{dataframe} objects as described above.}
#' }
#' @param essential_genes A 2-column dataframe with scaffold names in the first column and gene ID's in the second. Can contain duplicates. (\emph{Default: } \code{NULL}) 
#' @param taxonomy A dataframe containing taxonomy assigned to the scaffolds. The first column must contain the scaffold names. (\emph{Default: } \code{NULL})
#' @param additional A dataframe containing any additional data. The first column must contain the scaffold names. (\emph{Default: } \code{NULL})
#' @param kmer_pca (\emph{Logical}) Perform Principal Components Analysis of tetranucleotide frequencies of each scaffold and merge the scores of the 3 most significant axes. (\emph{Default: } \code{FALSE}) 
#' @param kmer_BH_tSNE (\emph{Logical}) Calculate Barnes-Hut t-Distributed Stochastic Neighbor Embedding (B-H t-SNE) representations of tetranucleotide frequencies using \code{\link[Rtsne.multicore]{Rtsne.multicore}} and merge the result. Additional arguments may be required for success, refer to the documentation of \code{\link[Rtsne.multicore]{Rtsne.multicore}}. (\emph{Default: } \code{FALSE}) 
#' @param verbose (\emph{Logical}) Whether to print status messages during the loading process. (\emph{Default: } \code{TRUE}) 
#' @param ... Additional arguments are passed on to \code{\link[Rtsne.multicore]{Rtsne.multicore}}. 
#'
#' @export
#' 
#' @return A dataframe (tibble) compatible with other mmgenome2 functions.
#' 
#' @importFrom tibble add_column as.tibble tibble
#' @importFrom digest digest
#' @importFrom Biostrings readDNAStringSet letterFrequency oligonucleotideFrequency reverseComplement
#' @importFrom BiocGenerics width
#' @importFrom dplyr mutate_all funs group_by left_join summarise_all
#' @importFrom stringr str_replace_all
#' @importFrom vegan rda scores
#' @import Rtsne.multicore
#' 
#' @examples 
#' \dontrun{
#'   library(mmgenome2)
#'   mm <- mmload(
#'     assembly = "path/to/assembly.fa",
#'     coverage = list(read.csv("path/to/coveragetable1.csv", col.names = TRUE),
#'                     read.csv("path/to/coveragetable2.csv", col.names = TRUE)),
#'     essential_genes = read.csv("path/to/ess_genes.txt", col.names = TRUE),
#'     verbose = TRUE
#'   )
#' }
#' 
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
mmload <- function(assembly,
                   coverage, 
                   essential_genes = NULL,
                   taxonomy = NULL,
                   additional = NULL,
                   kmer_pca = FALSE,
                   kmer_BH_tSNE = FALSE,
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
  mm <- tibble::tibble(scaffold = as.character(names(assembly)),
                       length = as.numeric(BiocGenerics::width(assembly)),
                       gc = round(as.numeric(Biostrings::letterFrequency(assembly, letters = c("CG"), as.prob=T))*100, digits = 2)
  )
  
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
      
      essential_genes[[1]] <- as.character(essential_genes[[1]])
      essential_genes[[2]] <- as.character(essential_genes[[2]])
      
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
