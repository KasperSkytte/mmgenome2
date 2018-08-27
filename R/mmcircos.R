#' Title
#'
#' @param bin1 
#' @param bin2 
#'
#' @export
#' 
#' @importFrom RSQLite dbConnect SQLite dbDisconnect
#' @importFrom DECIPHER FindSynteny Seqs2DB
mmsynteny <- function(bin1, bin2) {
  # Create sequence database
  dbConn <- RSQLite::dbConnect(RSQLite::SQLite(), ":memory:")
  DECIPHER::Seqs2DB(bin1, "XStringSet", dbConn, "bin1")
  DECIPHER::Seqs2DB(bin2, "XStringSet", dbConn, "bin2")
  # Align sequences
  synteny <- DECIPHER::FindSynteny(dbConn, maskRepeats = F, storage = 0.5)
  RSQLite::dbDisconnect(dbConn)
  # Map scaffold names
  scaffold_names <- data.frame(index1 = names(bin1),
                               factor1 = unique(sort(synteny[[2]][,"index1"])),
                               index2 = names(bin2),
                               factor2 = unique(sort(synteny[[2]][,"index2"])))
  synteny <- c(synteny, list(scaffold_names))
  return(synteny)
}

#' Title
#'
#' @param synteny 
#' @param min_score 
#' @param min_length 
#'
#' @export
#' 
#' @importFrom IRanges IRanges findOverlaps reduce
#' @importFrom S4Vectors subjectHits
#' @importFrom dplyr mutate filter select arrange group_by ungroup bind_rows
#' @importFrom tidyr gather spread
mmsynteny_filter <- function(synteny, min_score, min_length = 1000) {
  # Store scaffold names
  sn <- synteny[[5]]
  # Filter data to remove self alingments and 
  synteny_filter1 <- synteny[[2]] %>%
    as.data.frame() %>%
    dplyr::mutate(cnr = 1:nrow(.)) %>%
    # Remove complete self alignments
    dplyr::filter(!(index1 == index2 & start1 == start2 & end1 == end2)) %>%
    # Remove alignments with low scores og 
    dplyr::filter(score >= min_score) %>%
    # Remove short alignments
    dplyr::filter(end1 - start1 >= min_length | end2 - start2 >= min_length)
  
  synteny_filter2 <- synteny_filter1 %>%
    # Sorts data rowwise values according to scaffold number
    dplyr::select(index1, index2, start1, start2, end1, end2, strand) %>%
    dplyr::mutate(ct = 1:nrow(.)) %>%
    tidyr::gather("type", "value", index1, index2, start1, start2, end1, end2) %>%
    dplyr::mutate(pos = gsub("[a-zA-Z]", "", type), type = gsub("[0-9]", "", type)) %>%
    tidyr::spread(type, value) %>%
    dplyr::arrange(ct, index) %>%
    dplyr::group_by(ct) %>%
    dplyr::mutate(pos = row_number()) %>%
    tidyr::gather("type", "value", index, start, end) %>%
    dplyr::mutate(type = paste(type, pos, sep = "")) %>%
    dplyr::select(-pos) %>%
    tidyr::spread(type, value) %>%
    dplyr::arrange(index1, index2) %>%
    dplyr::ungroup() %>%
    dplyr::select(index1, index2, start1, start2, end1, end2, strand)
  
  synteny_merge <- synteny_filter2 %>%
    # Merges similar or overlapping alignments
    tidyr::mutate(con_type = paste(index1, index2, sep = "_")) %>%
    split(.$con_type) %>%
    lapply(.,
           function(rgn){
             pos1 <- IRanges::IRanges(rgn$start1, rgn$end1)
             pos2 <- IRanges::IRanges(rgn$start2, rgn$end2)
             rgng <- rgn %>%
               tidyr::mutate(grp1 = S4Vectors::subjectHits(IRanges::findOverlaps(pos1, IRanges::reduce(pos1))),
                             grp2 = S4Vectors::subjectHits(IRanges::findOverlaps(pos2, IRanges::reduce(pos2))),
                             grp = paste(grp1, grp2, sep = "_"))
             rgn_merge <- rgng %>%
               split(.$grp) %>%
               lapply(., function(mergelist){
                 merged <- mergelist[1,] %>%
                   tidyr::mutate(start1 = min(mergelist$start1),
                                 start2 = min(mergelist$start2),
                                 end1 = max(mergelist$end1),
                                 end2 = max(mergelist$end2))
                 return(merged)
               }) %>% 
               dplyr::bind_rows()
             return(rgn_merge)}) %>%
    dplyr::bind_rows() %>%
    dplyr::select(scaffold1 = index1,
                  scaffold2 = index2,
                  start1,
                  start2,
                  end1,
                  end2,
                  strand) %>%
    # Map correct scaffold names
    dplyr::mutate(scaffold1 = sn$index1[match(scaffold1, sn$factor1)],
                  scaffold2 = sn$index2[match(scaffold2, sn$factor2)]) %>%
    as.data.frame()
  return(synteny_merge)
}

#' Title
#'
#' @param bin1 
#' @param bin2 
#' @param synteny_filtered 
#' @param window_size 
#' @param out_type 
#' @param min_avg_id 
#' @param min_length 
#'
#' @export
#' 
#' @importFrom dplyr row_number bind_rows group_by mutate ungroup transmute filter
#' @importFrom tidyr separate
#' @importFrom IRanges IRangesList IRanges width
#' @importFrom Biostrings extractAt reverseComplement DNAStringSet letterFrequency
#' @importFrom BiocGenerics Map
#' @importFrom DECIPHER AlignSeqs ConsensusSequence
mmsynteny_align <- function(bin1,
                            bin2,
                            synteny_filtered,
                            window_size = 500,
                            out_type = "id_windows",
                            min_avg_id = 0.7, 
                            min_length = 1000) {
  # Define Iranges for each bin
  rgn1 <- IRanges::IRangesList(apply(synteny_filtered[c("start1", "end1")],
                            1, function(x){
                              IRanges::IRanges(x[[1]], x[[2]])
                            }))
  rgn2 <- IRanges::IRangesList(apply(synteny_filtered[c("start2", "end2")],
                            1, function(x){
                              IRanges::IRanges(x[[1]], x[[2]])
                            }))
  # Extract target sequence regions from each bin
  seq1 <- bin1[synteny_filtered$scaffold1]
  seq2 <- bin2[synteny_filtered$scaffold2]
  seq1ext <- Biostrings::extractAt(seq1, rgn1) %>% unlist()
  seq2ext <- Biostrings::extractAt(seq2, rgn2) %>%
    unlist() %>%
    # Reverse complement based on synteny strand info
    BiocGenerics::Map(function(seq, strand){
      if (strand == 1) {
        seq_out <- Biostrings::reverseComplement(seq)
      } else {
        seq_out <- seq
      }
    }, ., synteny_filtered$strand) %>%
    Biostrings::DNAStringSet()
  
  # Align sequences
  alignments <- BiocGenerics::Map(function(seq_list1, seq_list2, start1, start2){
    seq_pair <- Biostrings::DNAStringSet(list(seq_list1, seq_list2))
    DECIPHER::AlignSeqs(seq_pair)
  }, seq1ext, seq2ext)
  # Name intervals for tracking
  names(alignments) <- paste(names(seq1ext), names(seq2ext), sep = "_") %>%
    data.frame(n = .) %>%
    dplyr::group_by(n) %>%
    dplyr::mutate(grp = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(grp = paste(n, grp, sep = "_")) %>%
    unlist()
  
  # Calculate identity in windows
  id_windows <- BiocGenerics::Map(function(a, start1, start2, winsize, strand){
    # Generate consensus sequence
    cs <- DECIPHER::ConsensusSequence(a)
    # Divide consensus in windows
    ws <- ifelse(length(cs[[1]]) < winsize, length(cs[[1]]), winsize)
    start <- seq(1, length(cs[[1]]), ws)
    end <- unique(c(seq(ws, length(cs[[1]]), ws), length(cs[[1]])))
    # Calculate identity in windows
    cs_wid1 <- Biostrings::extractAt(cs[[1]], IRanges::IRanges(start, end)) %>%
      {rowSums(Biostrings::letterFrequency(., c("A", "T", "C", "G"))/IRanges::width(.))}
    # Orientate id values according to strand
    cs_wid2 <- ifelse(rep(strand == 1, length(cs_wid1)), rev(cs_wid1), cs_wid1)
    # Function to map windows to sequence coordinates
    mappos_fun <- function(seq, pos){
      mapped_pos <- Biostrings::extractAt(seq, IRanges::IRanges(1, pos)) %>%
        unlist() %>%
        Biostrings::letterFrequency(., c("A", "T", "C", "G")) %>%
        rowSums()
    }
    position1 <- start1 + mappos_fun(a[1], start) + (mappos_fun(a[1], end) - mappos_fun(a[1], start))/2
    position2 <- start2 + mappos_fun(a[2], start) + (mappos_fun(a[2], end) - mappos_fun(a[2], start))/2
    # Compile output data
    wid <- data.frame(position1,
                      position2,
                      id1 = cs_wid1,
                      id2 = cs_wid2)
    return(wid)
  }, alignments,
  synteny_filtered$start1,
  synteny_filtered$start2,
  window_size,
  synteny_filtered$strand) %>%
    # Convert list to dataframe
    dplyr::bind_rows(.id = "scaffold") %>%
    tidyr::separate(scaffold, c("scaffold1", "scaffold2", "grp"), sep = "_")
  
  # Convert to circos format
  id_graph1 <- id_windows[c("scaffold1", "position1", "id1", "scaffold2", "grp")] %>%
    setNames(c("scaffold", "position", "id_value", "match", "grp"))
  id_graph2 <- id_windows[c("scaffold2", "position2", "id2", "scaffold1", "grp")] %>%
    setNames(c("scaffold", "position", "id_value", "match", "grp"))
  id_graph <- dplyr::bind_rows(id_graph1, id_graph2) %>%
    dplyr::mutate(scaffold = as.numeric(scaffold))
  
  # Identity filtering
  id_graph <- id_graph %>%
    dplyr::group_by(scaffold, match, grp) %>%
    dplyr::filter(mean(id_value) > min_avg_id) %>%
    dplyr::ungroup()
  
  # Generate output type
  if (out_type == "id_windows"){
    out <- id_graph
  } else if (out_type == "alignment") {
    out <- alignments
  } else if (out_type == "all") {
    out <- list(id_graph, alignments)
  } else {stop("Unknown output type. Options are: 'id_windows', 'alignment' or 'all'")}
  
  return(out)
}

#' Title
#'
#' @param links 
#' @param cov_sample 
#' @param out_type 
#' @param window_size 
#' @param min_link 
#' @param min_cov 
#' @param cov_fp 
#' @param end_length 
#' @param remove_regions 
#'
#' @export
#' @importFrom dplyr mutate group_by summarise ungroup left_join filter select
#' @importFrom BiocGenerics Map
mmlink_filter <- function(links,
                          cov_sample,
                          out_type,
                          window_size = 2000,
                          min_link = 3,
                          min_cov = 5,
                          cov_fp = 1/50,
                          end_length = NULL,
                          remove_regions = NULL) {
  # Calculate average coverage for window size
  c <- cov_sample %>%
    dplyr::mutate(position = floor(position/window_size) * window_size) %>%
    dplyr::group_by(scaffold, position) %>%
    dplyr::summarise(coverage = mean(coverage)) %>%
    dplyr::ungroup()
  
  # Calculate sum of links for window size
  l <- links %>%
    dplyr::mutate(position1 = floor(position1/window_size) * window_size, 
                  position2 = floor(position2/window_size) * window_size) %>%
    dplyr::group_by(scaffold1, scaffold2, position1, position2) %>%
    dplyr::summarise(connections = sum(connections)) %>%
    dplyr::ungroup()
  
  # Filter links based on coverage and expected noisy links
  lf <- l %>%
    dplyr::left_join(c, by = c("scaffold1" = "scaffold", "position1" = "position")) %>%
    dplyr::left_join(c, by = c("scaffold2" = "scaffold", "position2" = "position"),
                     suffix = c("1", "2")) %>%
    dplyr::filter(connections >= min_link) %>%
    dplyr::filter(coverage1 > min_cov & coverage2 > min_cov) %>%
    dplyr::filter(connections/coverage1 > cov_fp & connections/coverage2 > cov_fp) %>%
    dplyr::select(-coverage1, -coverage2)
  
  #Filter links based on end connections
  if (is.null(end_length)){
    lff <- lf
  } else {
    end_fun <- function(position, position_max, end_length){
      if(position_max > 2 * end_length){
        sel <- !(position %in% end_length:(position_max - end_length))
      } else {sel <- T}
      return(sel)
    }
    lff <- lf %>%
      dplyr::group_by(scaffold1) %>%
      dplyr::filter(end_fun(position1, max(position1), end_length)) %>%
      dplyr::group_by(scaffold2) %>%
      dplyr::filter(end_fun(position2, max(position2), end_length))
  }
  
  # Filter links in coupled regions
  if (is.null(remove_regions)){
    lfff <- lff
  } else {
    remove_fun <- function(scaffold1, scaffold2, position1, position2, regions){
      rf <- dplyr::filter(regions, regions$scaffold1 == scaffold1 & regions$scaffold2 == scaffold2)
      range1 <- BiocGenerics::Map(function(x1, x2){x1:x2}, rf$start1, rf$end1) %>% unlist %>% unique()
      range2 <- BiocGenerics::Map(function(x1, x2){x1:x2}, rf$start2, rf$end2) %>% unlist %>% unique()
      match <- !(position1 %in% range1 & position2 %in% range2)
      return(match)
    }
    lfff <- lff %>%
      dplyr::filter(remove_fun(scaffold1, scaffold2, position1, position2, remove_regions))
  }
  
  
  # Prepare output based on type
  if (out_type == "network_links"){
    out <- lfff %>%
      dplyr::group_by(scaffold1, scaffold2) %>%
      dplyr::summarise(connections = sum(connections)) %>%
      dplyr::ungroup()
  } else if (out_type == "circos_links"){
    out <- lfff
  } else stop("Unknown output type. Options are 'network_links' or 'circos_links'")
  return(out)
}

#' Title
#'
#' @param assembly 
#' @param coverage 
#' @param profile 
#' @param links 
#' @param assembly_col 
#' @param gc 
#' @param line_width 
#' @param window_size 
#' @param print_pdf 
#'
#' @export
#' 
#' @importFrom dplyr bind_rows mutate filter group_by transmute summarise if_else select ends_with
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Biostrings letterFrequency extractAt
#' @importFrom IRanges IRanges width
#' @importFrom data.table fread
#' @import circlize
mmplot_circos <- function(assembly = get("assembly", envir = globalenv()),
                          coverage = NULL,
                          profile = NULL,
                          links = NULL,
                          assembly_col = NULL,
                          gc = TRUE,
                          line_width = 0.5,
                          window_size = 5000,
                          print_pdf = NULL) {
  if(length(assembly) > 350L)
    stop("Too many scaffolds to plot (>350). Make a subset and try again", call. = FALSE)
  # Convert input to lists
  profile <- if(!is.null(profile)) list(profile)
  links <- if(!is.null(links)) list(links)
  
  # load coverage files
  if(!is.null(coverage)) {
    filepaths <- coverage
    coverage <- list()
    #read files and name elements in the list based on the filename
    for (i in 1:length(filepaths)) {
      coverage[[tools::file_path_sans_ext(basename(filepaths[i]))]] <- data.table::fread(filepaths[i])
    }
  }
  
  # Status messages with how many tracks will be plotted
  message("Generating circos plot...")
  if(!is.null(assembly))
    message("Assembly tracks: 1")
  else
    stop("Assembly track must be provided", call. = FALSE)
  
  if(!is.null(coverage))
    message(paste0("Coverage tracks: ", length(coverage), sep = " "))
  else
    message("Coverage tracks: 0")
  
  if(!is.null(profile))
    message(paste0("Profile tracks: ", length(profile), sep = " "))
  else
    message("Profile tracks: 0")
  
  if(!is.null(links))
    message(paste0("Links tracks: ", length(links), sep = " "))
  else
    message("Links tracks: 0")
  
  #Prepare assembly data
  asmb <- assembly
  if (is.null(assembly_col) ){
    asmb_col <- rep(RColorBrewer::brewer.pal(12, "Paired"), ceiling(length(asmb)/12))[1:length(asmb)]
  } else
    asmb_col <- assembly_col
  
  ad <- data.frame(
    scaffold = names(asmb),
    start = 0,
    end = IRanges::width(asmb),
    color = asmb_col, 
    stringsAsFactors = F)
  
  # Determine plot output
  if(!is.null(print_pdf)){
    pdf(file = print_pdf[1],
        width = ifelse(is.na(print_pdf[2]), 7, print_pdf[2]),
        height = ifelse(is.na(print_pdf[3]), 7, print_pdf[3]))
  }
  # Build plot
  circlize::circos.par(gap.after = 0.5,
                       points.overflow.warning = F,
                       start.degree = 90,
                       unit.circle.segments = 2000,
                       cell.padding = c(0.005, 0, 0.005, 0),
                       track.margin = c(0.005, 0.005))
  circlize::circos.genomicInitialize(ad, plotType = NULL)
  
  # Assembly labels
  circlize::circos.track(ylim = c(0,05),
                         panel.fun = function(x, y) {
                           circlize::circos.text(CELL_META$xcenter,
                                                 CELL_META$ycenter,
                                                 CELL_META$sector.index, 
                                                 facing = "clockwise",
                                                 niceFacing = T,
                                                 cex = 0.5)},
                         track.height = 0.05,
                         bg.border = NA)
  
  # Prepare GC data
  if (isTRUE(gc)) {
    # Loop over scaffolds in assembly
    gcd <- lapply(asmb, function(x){
      # GC function
      gc_fun <- function(x){
        letters <- Biostrings::letterFrequency(x, c("C", "G"))
        identity <- sum(letters)/length(x)
        return(identity)}
      # Calculate gc content in windows
      if (length(x) >= window_size) {w = window_size} else {w = length(x)}
      start <- seq(1, length(x), w)
      end <- unique(c(seq(w, length(x), w), length(x)))
      gc <- Biostrings::extractAt(x, IRanges::IRanges(start, end)) %>%
        {rowSums(Biostrings::letterFrequency(., c("C", "G"))/IRanges::width(.))}
      return(data.frame(start = start - 1 + (end - start)/2,
                        end = start - 1 + (end - start)/2,
                        gc))
    }) %>%
      # Convert results to table
      dplyr::bind_rows(.id = "scaffold")
    
    gc_track <- circlize::get.current.track.index()
    circlize::circos.genomicTrack(gcd, ylim = c(min(gcd$gc), max(gcd$gc)), track.index = gc_track,
                                  panel.fun = function(region, value, ...) {
                                    xlim = CELL_META$xlim
                                    circlize::circos.genomicLines(region,
                                                                  value,
                                                                  col = rgb(0.46, 0.77, 0.47),
                                                                  lwd = line_width)},
                                  track.height = 0.05,
                                  bg.border = NA)
    
  }
  
  # Assembly
  circlize::circos.track(ylim = c(0, 1),
                         cell.padding = c(0.02, 1.00, 0.005, 1.00),
                         panel.fun = function(x, y) {
                           index = CELL_META$sector.numeric.index
                           xlim = CELL_META$xlim
                           circlize::circos.rect(xlim[1],0, xlim[2], 1,
                                                 col = ad$color[index],
                                                 lwd = 0.01)
                           circlize::circos.axis(h = "top",
                                                 major.at = seq(0, xlim[2], 25000),
                                                 major.tick.length = 0.2,
                                                 minor.ticks = 0,
                                                 labels = NULL,
                                                 lwd = 0.01)
                         },
                         track.height = 0.05,
                         bg.border = NA)
  
  # Prepare coverage data
  if (!is.null(coverage) ) {
    cd <- lapply(coverage, function(x){
      as.data.frame(x) %>%
        dplyr::mutate(scaffold = as.character(scaffold)) %>%
        # Subset to assembly
        dplyr::filter(scaffold %in% names(asmb)) %>%
        # Calculate window coverage
        dplyr::group_by(scaffold) %>%
        dplyr::transmute(start = (1 + floor(position/window_size)) * 5000 - 2500,
                  start = dplyr::if_else(start > max(position), max(position), start),
                  end = start,
                  coverage) %>%
        dplyr::group_by(scaffold, start, end) %>%
        dplyr::summarise(coverage = mean(coverage)) %>%
        as.data.frame()
    })
    lapply(cd, function(x){
      circlize::circos.genomicTrack(x, ylim = c(min(x$coverage), max(x$coverage)),
                                    panel.fun = function(region, value, ...) {
                                      xlim = CELL_META$xlim
                                      circlize::circos.genomicLines(region,
                                                                    value,
                                                                    col = "lightgrey",
                                                                    lwd = 0.01,
                                                                    area = T)},
                                    track.height = 0.05,
                                    bg.border = NA)
    }) %>% invisible()
  }
  
  # Prepare profile data
  if(!is.null(profile)){
    pd <- lapply(profile, function(x){
      dplyr::mutate(x, scaffold = as.character(scaffold)) %>%
        dplyr::filter(scaffold %in% names(asmb)) %>%
        dplyr::group_by(scaffold) %>%
        dplyr::mutate(start = (1 + floor(position/window_size)) * 5000 - 2500,
                      start = dplyr::if_else(start > max(position), max(position), start),
                      end = start,
                      color = if (exists('color', where = .)) {
                        color 
                      } else if (exists('match', where = .)) {
                        ad$color[match(match, ad$scaffold)]
                      } else "black",
                      grp = if (exists('grp', where = .)) grp else 1) %>%
        dplyr::select(scaffold,
                      start,
                      end,
                      value = dplyr::ends_with("value"),
                      color,
                      grp
        ) %>%
        dplyr::group_by(scaffold, start, end, grp, color) %>%
        dplyr::summarise(value = mean(value)) %>%
        as.data.frame()
    })
    lapply(pd, function(x){
      profile_track <- circlize::get.current.track.index() + 1
      for (i in unique(x$color)) {
        pf_match <- dplyr::filter(x, color %in% i)
        for (j in unique(pf_match$grp)) {
          pf_grp <- dplyr::filter(pf_match, grp %in% j)
          circlize::circos.genomicTrack(pf_grp,
                                        ylim = c(min(x$value), max(x$value)),
                                        track.index = profile_track,
                                        panel.fun = function(region, value, ...) {
                                          xlim = CELL_META$xlim
                                          circlize::circos.genomicLines(region,
                                                                        value,
                                                                        col = pf_grp$color[1],
                                                                        lwd = line_width)
                                          circlize::circos.rect(xlim[1],0, xlim[2], 1,
                                                                lwd = 0.01)},
                                        track.height = 0.03,
                                        bg.border = NA)
        }
      }
    }
    )
  }
  
  # Prepare links data
  if ( !is.null(links) ){
    ld <- lapply(links, function(x){
      as.data.frame(x) %>%
        dplyr::mutate(scaffold1 = as.character(scaffold1),
                      scaffold2 = as.character(scaffold2)) %>%
        dplyr::filter(scaffold1 %in% names(asmb) & scaffold2 %in% names(asmb)) %>%
        dplyr::transmute(scaffold1,
                         scaffold2,
                         start1 = if (exists('position1', where = .)) position1 else start1,
                         start2 = if (exists('position2', where = .)) position2 else start2,
                         end1 = if (exists('position1', where = .)) position1 else end1,
                         end2 = if (exists('position2', where = .)) position2 else end2,
                         connections = if (exists('connections', where = .)) connections else 1,
                         color = if (exists('color', where = .)) color else rgb(0.42, 0.68, 0.84),
                         line_width = if (exists('line_width', where = .)) line_width else 1,
                         arc_height = if (exists('arc_height', where = .)) arc_height else 0.4)
    })
    lapply(ld, function(x){
      circlize::circos.genomicLink(x[,c("scaffold1", "start1", "end1")],
                                   x[,c("scaffold2", "start2", "end2")],
                                   col = x$color,
                                   lwd = x$line_width,
                                   h.ratio = x$arc_height[1], 
                                   border = NA)
    })
  }

  # Finish
  circlize::circos.clear()
  
  # Stop printing to pdf
  if(!is.null(print_pdf)){
    invisible(dev.off())
  }
  invisible()
}
#' Title
#'
#' @param data 
#'
#' @export
#' 
#' @importFrom data.table fread
mmload_links <- function(data) {
  link_files <- list.files(path = data, pattern = ".*_link.csv", full.names = T)
  names(link_files) <- basename(link_files) %>%
    gsub(".csv", "", .) %>%
    gsub("-", "_", .)
  links <- lapply(link_files, function(x){
    data.table::fread(x)
  })
  return(links)
}