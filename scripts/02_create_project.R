suppressMessages({
  library(ArchR)
  library(dplyr)
  library(parallel)
  library(hexbin)
  #library(Signac)

  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BSgenome.Mmusculus.UCSC.mm10)

})

getGenomeInfo <- function(genome_id) {
  if (genome_id == "mm10") {
    data("geneAnnoMm10")
    data("genomeAnnoMm10")
    geneAnno <- geneAnnoMm10
    genomeAnno <- genomeAnnoMm10
  } else if (genome_id == "hg38") {
    data("geneAnnoHg38")
    data("genomeAnnoHg38")
    geneAnno <- geneAnnoHg38
    genomeAnno <- genomeAnnoHg38
  } else if (genome_id == "hg19") {
    data("geneAnnoHg19")
    data("genomeAnnoHg19")
    geneAnno <- geneAnnoHg19
    genomeAnno <- genomeAnnoHg19
  } else {
    stop("genome not found")
  }
  return(list(geneAnno=geneAnno, genomeAnno=genomeAnno))
}

subsetProj <- function(origProj, origClus, sel, newPrefix, clustering=NA) {

  proj_pth <- paste0("r_objects/", newPrefix, "_proj.rds") 
  if (file.exists(proj_pth)) {
    proj <- readRDS(proj_pth)
  } else {

    df <- data.frame(origProj@cellColData)
    cellNames <- rownames(df[which(sel(df)),])

    print(length(cellNames))

    proj <- subsetArchRProject(origProj, cells=cellNames, outputDirectory=paste0("ArchRSubset_",newPrefix),dropCells=FALSE)
    proj <- addClustering(proj, newPrefix, clustering=clustering)

    saveRDS(proj, proj_pth)
  }

  if (!is.na(clustering) & !(clustering %in% colnames(proj@cellColData))) {
    print("updating clustering")
    proj <- addClustering(proj, newPrefix, clustering=clustering)
    saveRDS(proj, proj_pth)
  } 

  return(proj)

}

createProj <- function(root, genome_id, names, clusters=NA) {


  print("create proj")
  print(getwd())

  proj_pth <- "r_objects/All_proj.rds"
  if (file.exists(proj_pth)) {
    proj <- readRDS(proj_pth)
  } else {

    g <- getGenomeInfo(genome_id)

    ArrowFiles <- paste0("Arrow/",names,".arrow")
    print(ArrowFiles)

    proj <- ArchRProject(
      ArrowFiles = ArrowFiles, geneAnnotation = g$geneAnno, 
      genomeAnnotation = g$genomeAnno, outputDirectory = "ArchRProject"
    )
    #proj <- filterDoublets(proj)

    proj <- addClustering(proj, "All",clustering=clusters)
    proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Sample")#,threads=1)
    getGroupBW(proj,groupBy="Sample",threads=1)
    saveRDS(proj, proj_pth)
  }

  return(proj)
}

addClustering <- function(proj, prefix, clustering) {

  tmp <- strsplit(clustering,"_")[[1]]
  dims <- c(as.numeric(tmp[2]))
  res_opt <- c(as.numeric(tmp[3])/1000)

  print(res_opt)
  print(dims)

  for (d in dims) {

    nm <- paste0(prefix,"_",d)

    if (!(paste0(nm,"_LSI_UMAP") %in% names(proj@embeddings))) {
      proj <- addIterativeLSI(proj,name=paste0(nm,"_LSI"),dimsToUse=1:d)
      proj <- addUMAP(proj,reducedDims=paste0(nm,"_LSI"),name=paste0(nm,"_LSI_UMAP"))
    }

    for (res in res_opt) {

      nm2 <- paste0(nm,"_",res*1000)
      nm3 <- paste0(nm2,"_LSI")

      if (!(nm3 %in% colnames(proj@cellColData))) { 
        print(paste0("clustering ",nm2))
        proj <- addClusters(proj,reducedDims=paste0(nm,"_LSI"),name=paste0(nm2,"_LSI"),resolution=res)
      }

    }
  }

  return(proj)
}

addMatrices <- function(proj, dr, clustering) {


  print("adding matrices")
  print(getwd())

  proj_pth <- paste0("r_objects/", clustering, "_proj.rds")
  if (file.exists(proj_pth)) {
    proj <- readRDS(proj_pth)
  } else {

    getGroupBW(proj,groupBy=clustering)
    proj <- addImputeWeights(proj, reducedDims=dr)
    proj <- addGroupCoverages(proj, groupBy=clustering)#,threads=1)
    proj <- addReproduciblePeakSet(proj, groupBy=clustering)
    proj <- addPeakMatrix(proj)

    proj <- addMotifAnnotations(proj,force=TRUE)
    proj <- addDeviationsMatrix(proj,force=TRUE)

    proj <- addCoAccessibility(ArchRProj=proj,reducedDims=dr,maxDist=500000)

    saveRDS(proj,proj_pth)
  }
  return(proj)
}


initProj <- function(input_files, genome_id) {

  g <- getGenomeInfo(genome_id) 

  print(input_files$files)
  print(input_files$names)

  ArrowFiles <- createArrowFiles(
    inputFiles = input_files$files,
    sampleNames = input_files$names,
    geneAnno = g$geneAnno,
    genomeAnno = g$genomeAnno
  )

  doubScores <- addDoubletScores(ArrowFiles)
}

addArchRThreads(threads = 1)

args = commandArgs(trailingOnly=TRUE)

source(paste0("scripts/",args[1],".R"))

addArchRThreads(threads = 1)

print(root)
if (!dir.exists(root)) {
  dir.create(root)
}

if (!dir.exists(paste0(root,"Arrow"))) {
  print("creating arrow files...")
  dir.create(paste0(root,"Arrow"))
  wd <- paste0(root, "Arrow")
  print(wd)
  setwd(wd)
  initProj(input_files,genome_id) 
  setwd("..")
} else {
  setwd(root)
}

if (!dir.exists("r_objects")) { dir.create("r_objects") }

main_proj <- createProj(root, genome_id, input_files$names, clusters=main_clusters)
print(colnames(main_proj@cellColData))

if (!is.na(main_clusters)) {
  main_proj <- addMatrices(main_proj, main_dr, main_clusters)
}

if (!is.na(sub_prefix)) {
  sub_proj <- subsetProj(main_proj, main_clusters, sub_sel, sub_prefix, sub_clusters)
  if (!is.na(sub_clusters)) {
    sub_proj <- addMatrices(sub_proj, sub_dr, sub_clusters)
  }
}







