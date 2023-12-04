#' @title Load and combine data for use with other mmgenome2 functions
#'
#' @description Loads, validates and combines multiple aspects of metagenome data into one dataframe for use with all mmgenome2 functions, including scaffold assembly sequences, scaffold coverage, essential genes, taxonomy, and more.
#'
#' @param assembly (\emph{required}) A character string with the path to the assembly FASTA file, or the assembly as already loaded with \code{\link{readDNAStringSet}}.
#' @param coverage (\emph{required}) A path to a folder to scan for coverage files, or otherwise a named \code{vector}, \code{data.frame}, or a \code{list} hereof containing coverage of each scaffold. The prefix \code{"cov_"} will be appended to all coverage column names in the output so that \code{\link{mmstats}} and \code{\link{mmplot_cov_profiles}} know which columns are coverage columns.
#' \describe{
#'   \item{\code{vector}}{If provided as a vector, the elements of the vector must be named by the scaffold names exactly matching those of the assembly.}
#'   \item{\code{data.frame}}{If provided as a dataframe, the first column must contain the scaffold names exactly matching those of the assembly, and any additional column(s) contain coverage of each scaffold.}
#'   \item{\code{list}}{If provided as a list, it must contain any number of \code{vector}'s or \code{data.frame}'s as described above. If names are assigned to the objects in the list, then they will be used as column names in the output (does not apply to any dataframes that may have more than 2 columns, however).}
#'   \item{\code{path}}{If a path to a folder is provided, then all files with filenames ending with \code{"_cov"} will be loaded (by the \code{\link[data.table]{fread}} function) into a list of \code{data.frame}'s and treated as if a \code{list} of \code{data.frame}'s were provided. The filenames (stripped from extension and \code{"_cov"}) will then be used as column names in the output. \strong{Note: only the first 2 columns will be used in the loaded files!}}
#' }
#' @param essential_genes Either a path to a CSV file (comma-delimited ",") containing the essential genes, or a 2-column dataframe with scaffold names in the first column and gene ID's in the second. Can contain duplicates. (\emph{Default: } \code{NULL})
#' @param taxonomy A dataframe containing taxonomy assigned to the scaffolds. The first column must contain the scaffold names. (\emph{Default: } \code{NULL})
#' @param additional A dataframe containing any additional data. The first column must contain the scaffold names. (\emph{Default: } \code{NULL})
#' @param kmer_pca (\emph{Logical}) Perform Principal Components Analysis of kmer nucleotide frequencies (kmer size defined by \code{kmer_size}) of each scaffold and merge the scores of the 3 most significant axes. (\emph{Default: } \code{FALSE})
#' @param kmer_BH_tSNE (\emph{Logical}) Calculate Barnes-Hut t-Distributed Stochastic Neighbor Embedding (B-H t-SNE) representations of kmer nucleotide frequencies (kmer size defined by \code{kmer_size}) using \code{\link[Rtsne]{Rtsne}} and merge the result. Additional arguments may be required for success (passed on through \code{...}), refer to the documentation of \code{\link[Rtsne]{Rtsne}}. This is done in parallel, thus setting the \code{num_threads} to the number of available cores may greatly increase the calculation time of large data. (\emph{Default: } \code{FALSE})
#' @param kmer_uwot (\emph{Logical}) Calculating UWOT representations of kmer nucleotide frequencies (kmer size defined by \code{kmer_size}) using \code{\link[uwot]{uwot}} and merge the result. Additional arguments may be required for success (passed on through \code{...}), refer to the documentation of \code{\link[Rtsne]{uwot}}. This is done in parallel, thus setting the \code{n_threads} to the number of available cores may greatly increase the calculation time of large data. (\emph{Default: } \code{FALSE})

#' @param kmer_size The kmer frequency size (k) used when \code{kmer_pca = TRUE} or \code{kmer_BH_tSNE = TRUE}. The default is tetramers (\code{k = 4}). (\emph{Default: } \code{4})
#' @param verbose (\emph{Logical}) Whether to print status messages during the loading process. (\emph{Default: } \code{TRUE})
#' @param ... Additional arguments are passed on to \code{\link[Rtsne]{Rtsne}}.
#'
#' @export
#'
#' @return A dataframe (tibble) compatible with other mmgenome2 functions.
#'
#' @importFrom tibble add_column as.tibble tibble
#' @importFrom magrittr %>%
#' @importFrom Biostrings width readDNAStringSet letterFrequency oligonucleotideFrequency reverseComplement
#' @importFrom dplyr mutate_all funs group_by left_join summarise_all
#' @importFrom stringr str_replace_all str_remove
#' @importFrom data.table fread
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' \dontrun{
#' library(mmgenome2)
#' mm <- mmload(
#'   assembly = "path/to/assembly.fa",
#'   coverage = list(
#'     nameofcoverage1 = read.csv("path/to/coveragetable1.csv", col.names = TRUE),
#'     nameofcoverage2 = read.csv("path/to/coveragetable2.csv", col.names = TRUE)
#'   ),
#'   essential_genes = "path/to/ess_genes.txt",
#'   verbose = TRUE
#' )
#' mm
#' }
#'
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
mmload <- function(assembly,
                   coverage = NULL,
                   essential_genes = NULL,
                   taxonomy = NULL,
                   additional = NULL,
                   kmer_pca = FALSE,
                   kmer_BH_tSNE = FALSE,
                   kmer_uwot = FALSE, 
                   kmer_size = 4L,
                   verbose = TRUE,
                   ...) {
  ##### Assembly #####
  # Load assembly sequences from the provided file path or object
  if (isTRUE(verbose)) {
    message("\nLoading assembly...")
  }
  if (is.character(assembly)) {
    assembly <- Biostrings::readDNAStringSet(assembly, format = "fasta")
  } else if (!class(assembly) == "DNAStringSet") {
    stop("The assembly must either be an object of class \"DNAStringSet\" loaded with readDNAStringSet() from the Biostrings package, or a file path to the assembly FASTA file to be loaded.", call. = FALSE)
  }

  # check if a different assembly already exists in global environment
  if (!exists("assembly", where = .GlobalEnv, envir = .GlobalEnv)) {
    assign("assembly", assembly, envir = .GlobalEnv)
  } else if (exists("assembly", where = .GlobalEnv, envir = .GlobalEnv)) {
    checkReqPkgs("digest")
    if (!identical(digest::digest(assembly), digest::digest(get("assembly", envir = .GlobalEnv)))) {
      userChoice <- readline(prompt = "A different object named \"assembly\" already exists in the global environment. Do you want to overwrite it? (y/n or ENTER/ESC): ")
      if (any(tolower(userChoice) %in% c("y", "yes", ""))) {
        assign("assembly", assembly, envir = .GlobalEnv)
      } else {
        stop("Aborted by user.", call. = FALSE)
      }
    }
  }

  # duplicate sequence names are not allowed
  if (any(duplicated(names(assembly)))) {
    stop("The assembly contains duplicate sequence names", call. = FALSE)
  }

  ##### Create base mm object #####
  if (isTRUE(verbose)) {
    message("Calculating GC content...")
  }
  mm <- tibble::tibble(
    scaffold = as.character(names(assembly)),
    length = as.integer(Biostrings::width(assembly)),
    gc = round(as.numeric(Biostrings::letterFrequency(assembly, letters = c("CG"), as.prob = T)) * 100, digits = 2)
  )

  ##### Coverage #####
  if (!is.null(coverage)) {
    if (isTRUE(verbose)) {
      message("Loading coverage data...")
    }
    if (is.character(coverage)) {
      filepaths <- list.files(
        path = coverage,
        full.names = TRUE,
        all.files = FALSE,
        recursive = FALSE,
        ignore.case = TRUE
      )
      filepaths <- filepaths[grepl("*._cov$", tools::file_path_sans_ext(filepaths))]
      path <- coverage
      if (length(filepaths) > 0) {
        filenames <- basename(filepaths)
        coverage <- list()
        for (i in 1:length(filenames)) {
          coverage[[stringr::str_remove(tools::file_path_sans_ext(filenames)[i], "_cov$")]] <- data.table::fread(filepaths[i], data.table = FALSE)[, 1:2]
        }
        if (isTRUE(verbose)) {
          message(paste0("  Found the following ", length(filenames), " coverage files in the folder \"", path, "\":\n    ", paste0(filenames, collapse = "\n    ")))
        }
      } else {
        stop("No files with a filename ending with \"_cov\" were found in the folder \"", coverage, "\"", call. = FALSE)
      }
    }
    beforeMerge <- ncol(mm)
    mm <- mmmerge(
      x = mm,
      y = coverage,
      type = "coverage"
    )
    colnames(mm)[c((beforeMerge + 1):ncol(mm))] <- paste0("cov_", colnames(mm)[c((beforeMerge + 1):ncol(mm))])
  } else {
    warning("No coverage data loaded, this may cause trouble with other mmgenome2 functions", call. = FALSE)
  }

  ##### Essential genes #####
  if (!is.null(essential_genes)) {
    if (isTRUE(verbose)) {
      message("Loading essential genes...")
    }
    if (is.character(essential_genes)) {
      if (length(essential_genes) == 1) {
        essential_genes <- read.csv(essential_genes,
          comment.char = "#",
          header = TRUE,
          colClasses = "character"
        )
        if (any(tolower(colnames(essential_genes)) == "orf")) {
          essential_genes <- essential_genes[, -which(colnames(essential_genes) == "orf")]
        }
      }
    }
    if (is.data.frame(essential_genes) & ncol(essential_genes) == 2) {
      essential_genes[[1]] <- as.character(essential_genes[[1]])
      essential_genes[[2]] <- as.character(essential_genes[[2]])

      # replace all values to only contain alpha-numerics and dots "."
      essential_genes[, -1] <- dplyr::mutate_all(essential_genes[, -1, drop = FALSE], dplyr::funs(stringr::str_replace_all(., "[^[:alnum:].]", "")))

      colnames(essential_genes) <- c("scaffold", "geneID")
      essential_genes <- essential_genes %>%
        dplyr::group_by(scaffold) %>%
        dplyr::summarise_all(dplyr::funs(paste(., collapse = ",")))

      mm <- dplyr::left_join(mm,
        essential_genes,
        by = "scaffold"
      )
    } else {
      stop("Essential genes must be a 2 column table, where the first column contains the sequence names exactly matching those of the assembly, and the second column the gene names/IDs.", call. = FALSE)
    }
  }

  ##### calculate kmer nucleotide frequencies #####
  if (isTRUE(kmer_pca) || isTRUE(kmer_BH_tSNE) || isTRUE(kmer_uwot)) {
    if (is.numeric(kmer_size)) {
      if (isTRUE(verbose)) {
        message(paste0(
          "Calculating kmer (k=",
          as.integer(kmer_size),
          ") nucleotide frequencies in the assembly sequences..."
        ))
      }
      kmer_fwd <- Biostrings::oligonucleotideFrequency(assembly,
        width = as.integer(kmer_size),
        as.prob = TRUE,
        with.labels = TRUE
      )
      kmer_revC <- Biostrings::oligonucleotideFrequency(Biostrings::reverseComplement(assembly),
        width = as.integer(kmer_size),
        as.prob = TRUE
      )
      kmer <- (kmer_fwd + kmer_revC) / 2 * 100
      write.csv(file="kmer_frequencies.csv", kmer)
    } else {
      stop("kmer_size must be a positive integer larger than 0", call. = FALSE)
    }
    
  }

  ##### PCA of tetranucleotides #####
  if (isTRUE(kmer_pca)) {
    checkReqPkg("vegan")
    
    if (isTRUE(verbose)) {
      message(paste0(
        "Calculating principal components of kmer (k=",
        as.integer(kmer_size),
        ") nucleotide frequencies..."
      ))
    }
    PCA_res <- kmer %>%
      vegan::rda() %>%
      vegan::scores(choices = 1:3, display = "sites") %>%
      tibble::as.tibble()
    mm <- tibble::add_column(mm,
      PC1 = PCA_res[[1]],
      PC2 = PCA_res[[2]],
      PC3 = PCA_res[[3]]
    )
  }

  ##### BH tSNE of tetranucleotides #####
  if (isTRUE(kmer_BH_tSNE)) {
    checkReqPkg("Rtsne", "To install with support for multithreading run:\n  remotes::install_github('kasperskytte/Rtsne@openmp')\notherwise just install from CRAN.")
    if (isTRUE(verbose)) {
      message(paste0(
        "Calculating Barnes-Hut t-Distributed Stochastic Neighbor Embedding representations of kmer (k=",
        as.integer(kmer_size),
        ") nucleotide frequencies..."
      ))
    }
    set.seed(42) # Sets seed for reproducibility
    tSNE_res <- Rtsne::Rtsne(kmer,
      verbose = verbose,
      check_duplicates = F,
      ...
    )[["Y"]] %>%
      tibble::as.tibble()
    mm <- tibble::add_column(mm,
      tSNE1 = tSNE_res[[1]],
      tSNE2 = tSNE_res[[2]]
    )
    umap_res <- uwot::umap(kmer,
      n_neighbors = 15, 
      learning_rate = 0.5, 
      init = "random", 
      n_epochs = 20, 
      n_threads = 10,
      ...
    ) %>%
      tibble::as.tibble()
    mm <- tibble::add_column(mm,
      tSNE1 = umap_res[[1]],
      tSNE2 = umap_res[[2]]
    )
  }


  ##### UWOT/UMAP of tetranucleotides
  if (isTRUE(kmer_uwot)) {
    checkReqPkg("uwot", "To install uwot/umap run:\n install.packages('uwot')\notherwise just install from CRAN.")
    if (isTRUE(verbose)) {
      message(paste0(
        "Calculating  UMAP representations of kmer (k=",
        as.integer(kmer_size),
        ") nucleotide frequencies..."
      ))
    }
    umap_res <- uwot::umap(kmer,
      n_neighbors = 15, 
      learning_rate = 0.5, 
      init = "random", 
      n_epochs = 20, 
      n_threads = 10,
      ...
    ) %>%
      tibble::as.tibble()
    mm <- tibble::add_column(mm,
      umap1 = umap_res[[1]],
      umap2 = umap_res[[2]]
    )
  }
  ##### Taxonomy #####
  if (!is.null(taxonomy)) {
    if (isTRUE(verbose)) {
      message("Loading taxonomy...")
    }
    if (is.character(taxonomy)) {
      if (length(taxonomy) == 1) {
        taxonomy <- read.csv(taxonomy,
          comment.char = "#",
          header = TRUE,
          colClasses = "character"
        )
      }
    }
    if (any(class(taxonomy) %in% c("data.frame", "tbl", "tbl_df")) | is.atomic(taxonomy)) {
      its_a_vector <- if (is.atomic(taxonomy)) TRUE else FALSE
      replacements <- c(
        " <phylum>" = "",
        "unclassified Bacteria" = "Unclassified Bacteria",
        "Fibrobacteres/Acidobacteria group" = "Acidobacteria",
        "Bacteroidetes/Chlorobi group" = "Bacteroidetes",
        "delta/epsilon subdivisions" = "Deltaproteobacteria",
        "Chlamydiae/Verrucomicrobia group" = "Verrucomicrobia"
      )
      for (i in 1:length(replacements)) {
        taxonomy <- lapply(taxonomy, gsub, pattern = paste(names(replacements)[i]), replacement = replacements[i])
      }
      if (!isTRUE(its_a_vector)) {
        taxonomy <- tibble::as.tibble(taxonomy)
      } else if (isTRUE(its_a_vector)) {
        taxonomy <- unlist(taxonomy)
      }
    }
    mm <- mmmerge(
      x = mm,
      y = taxonomy,
      type = "taxonomy"
    )
  }

  ##### Additional data #####
  if (!is.null(additional)) {
    if (isTRUE(verbose)) {
      message("Loading additional data...")
    }
    mm <- mmmerge(
      x = mm,
      y = additional,
      type = "additional"
    )
  }

  ##### Fix colnames and return #####
  newColnames <- stringr::str_replace_all(colnames(mm), "[^[:alnum:]_.]", "") # only alpha-numerics, dots and underscores are allowed
  if (any(duplicated(newColnames))) {
    stop("Only alpha-numeric characters, dots and underscores are allowed in column names", call. = FALSE)
  }
  colnames(mm) <- newColnames
  if (isTRUE(verbose)) {
    message("Done!")
  }
  return(mm)
}
