#' @title Plot the most abundant bins as a heatmap plot
#'
#' @description Plots a heatmap of the coverage profiles for the bins.
#'
#' @param mm (\emph{required}) A dataframe loaded with \code{\link{mmload}}.
#' @param BIN_COL Column to group scaffolds by
#' @param CLASSIFICATION Optional: Column with gtdb classification in the format d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Tissierellales;f__Sedimentibacteraceae;g__Sedimentibacter;s__
#' @param TOPN Optional Number of "bins" to display
#' @param tax_add Optional: Taxonomic levels to add from the gtdb column c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
#'
#' @export
#' #' @import mmgenome2
#'
#' @return A ggplot object. #'
#'
#' @examples
#' library(mmgenome2)
#' data(mmgenome2)
#' mmgenome2
#' mm_heatmap(mmgenome2,
#'   BIN_COL="taxonomy",
#'   TOPN=20
#' )
#' @author Kasper Skytte Andersen \email{ksa@@bio.aau.dk}
#' @author Rasmus Kirkegaard \email{rhk@@bio.aau.dk}


mm_heatmap <-function(mm=NULL,CLASSIFICATION=NULL,BIN_COL=NULL,TOPN=20,tax_add=NULL) {
  BIN_COL<-as.name(BIN_COL)
  ### Bins: Calculate abundance for each bin
  bins_abundance <- mm %>%
    pivot_longer(cols = starts_with("cov_"), names_to = "SeqID", values_to = "Coverage") %>%
    mutate(bp_mapped = length * Coverage) %>%
    group_by_("SeqID", BIN_COL) %>%
    summarise(Mapped = sum(bp_mapped)) %>% 
    select_("SeqID", BIN_COL, "Mapped") %>%
    group_by(SeqID) %>%
    mutate(Abundance=Mapped/sum(Mapped)*100) %>%
    select_("SeqID",BIN_COL, "Abundance") %>%
    pivot_wider(names_from = SeqID, values_from = Abundance)
  bins_abundance<-as.data.frame(bins_abundance %>%
                  select(as.name(BIN_COL),starts_with("cov_")))
  if (!is.null(CLASSIFICATION)) {
  ### Bins: Get taxonomy for each bin (gtdb)
  bins_tax<-mm %>% select_(BIN_COL, CLASSIFICATION) %>%
    distinct() %>%
    separate_(col = CLASSIFICATION, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>% 
    select_(BIN_COL, "Kingdom","Phylum","Class","Order","Family","Genus","Species") %>% filter(!is.na(eval(BIN_COL)))
  # Combine TAX+abundance
  bins_abundance<- bins_abundance %>%
    left_join(bins_tax)
  }
  bins_abundance<-bins_abundance %>%
    dplyr::rename_(Display=BIN_COL)
  bins_abundance$Display[which(is.na(bins_abundance$Display))]<-"Unbinned"
  
  
  ## Make a name variable that can be used instead of tax_aggregate to display multiple levels
  suppressWarnings(
    if (!is.null(tax_add)) {
      bins_abundance <- data.frame(bins_abundance, DisplayMod = apply(bins_abundance[, c(tax_add, "Display")], 1, paste, collapse = "; "))    
      bins_abundance<- bins_abundance %>% mutate(Display=DisplayMod)
      }
  )
  top_bins<-bins_abundance %>% select(Display,starts_with("cov_")) %>% pivot_longer(cols = starts_with("cov_")) %>% group_by(Display) %>% summarise(mean_cov=mean(value)) %>% arrange(desc(mean_cov)) %>% top_n(TOPN)

  heatmap<-bins_abundance %>% filter(Display %in% top_bins$Display) %>% mutate(Display=factor(Display,levels=rev(top_bins$Display))) %>% melt %>% ggplot(aes(x = variable,y = Display))+geom_tile(aes(fill=value))+
    theme(
      axis.text.y = element_text(size = 12, color = "black", vjust = 0.4),
      axis.text.x = element_text(size = 10, color = "black", vjust = 0.5, angle = 90, hjust = 1),
      axis.title = element_blank(),
      text = element_text(size = 8, color = "black"),
      axis.line = element_blank(),
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      title = element_text(size = 8),
      panel.background = element_blank()
    )+
    scale_fill_gradientn(colours = c("white","red"), na.value = "white")+
    geom_text(aes(label=round(value,digits = 2)),size = 4, colour = "grey10", check_overlap = TRUE)
  return(heatmap)
}
