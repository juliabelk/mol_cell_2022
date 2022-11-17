root <- paste0("./fig6_outputs/")

genome_id <- "mm10"

input_files <- data.frame(
  files=c(
    "../../fragments/GSM5393356_BMDM_CTR_fragments.tsv.gz",
    "../../fragments/GSM5393357_BMDM_IFNG_fragments.tsv.gz",
    #"../../fragments/GSM5393358_BMDM_IL4_fragments.tsv.gz"#,
    "../../fragments/GSM5393359_BMDM_preIL4_fragments.tsv.gz",
    "../../fragments/GSM5393360_BMDM_preIL4_IFNG_fragments.tsv.gz"
  ),
  names=c(
    "CTR",
    "IFNG",
    "preIL4",
    "preIL4_IFNG"
  ),
  stringsAsFactors=FALSE
)

main_dr <- "All_10_LSI"
main_clusters <- "All_10_400_LSI"

sub_prefix <- "Subset"
sub_dr <- "Subset_10_LSI"
sub_clusters <- "Subset_10_500_LSI"

valid_bc <- read.delim(paste0("barcode_info/fig6_cellColData_trim.csv"),sep=",")[[1]]
sub_sel <- function(x) {
  return( rownames(x) %in% valid_bc )
}

get_constraints <- function(proj,scrna) {
  groupList <- SimpleList(
    preIL4 = SimpleList(
      ATAC = proj$cellNames[which(proj$Sample == "preIL4")],
      RNA = rownames(scrna@meta.data)[which(scrna@meta.data$orig.ident == "preIL4")]
    ),  
    preIL4_IFNG = SimpleList(
      ATAC = proj$cellNames[which(proj$Sample == "preIL4_IFNG")],
      RNA = rownames(scrna@meta.data)[which(scrna@meta.data$orig.ident == "preIL4_IFNG")]
    ),  
    IFNG = SimpleList(
      ATAC = proj$cellNames[which(proj$Sample == "IFNG")],
      RNA = rownames(scrna@meta.data)[which(scrna@meta.data$orig.ident == "IFNG")]
    ),  
    CTR = SimpleList(
      ATAC = proj$cellNames[which(proj$Sample == "CTR")],
      RNA=rownames(scrna@meta.data)[which(scrna@meta.data$orig.ident == "CTR")]
    )   
  )
  return(groupList)
}




