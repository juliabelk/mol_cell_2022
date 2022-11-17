suppressMessages({
  library(ArchR)
  library(Seurat)
  library(dplyr)
})


scrna_integrate <- function(root,scrna_pth, proj, dr, clusters,getGL) {
  if (!("predCellRNA" %in% colnames(proj@cellColData))) {
    print("performing rna integration...")
    scrna <- readRDS(scrna_pth)
    groupList <- getGL(proj,scrna)

    proj@reducedDims[["IterativeLSI"]] <- proj@reducedDims[[dr]]

    proj <- addGeneIntegrationMatrix(
      ArchRProj=proj,
      useMatrix="GeneScoreMatrix",
      matrixName="GeneIntegrationMatrix",
      reducedDims="IterativeLSI",
      seRNA=scrna,
      addToArrow=TRUE,
      groupList=groupList,
      groupRNA="Phase",
      nameCell="predCellRNA",
      nameGroup="predPhaseRNA",
      nameScore="predScoreRNA",
      force=TRUE
    )

    proj@reducedDims[["IterativeLSI"]] <- NULL

    proj_pth <- paste0(root, "r_objects/", clusters, "_proj.rds")
    saveRDS(proj, proj_pth)
  }
  return(proj)
}



