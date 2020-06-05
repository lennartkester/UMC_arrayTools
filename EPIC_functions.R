
parameters <- list()
parameters$dataFolder <- "T:/pathologie/KMBP/EPIC/data/"
parameters$controlFolder <- "T:/pathologie/KMBP/EPIC/data/control/"

packages <- c("shiny","BiocManager")
bioPackages <- c("conumee","minfi","minfiData","minfiDataEPIC","IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

if (length(setdiff(packages, rownames(installed.packages()))) > 0){
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

if (length(setdiff(bioPackages, rownames(installed.packages()))) > 0){
  BiocManager::install(setdiff(bioPackages, rownames(installed.packages())))  
}

library("minfi") 
library("conumee") 
library("minfiData") 
library('minfiDataEPIC')

loadEPICFolders <- function(){
  inputChoicesEPIC <- list.dirs(parameters$dataFolder,recursive = F,full.names = F)
}

loadAndProcessSamples <- function(sampleName){
  patient<-read.metharray.exp(paste0(parameters$dataFolder,sampleName))
  gmSet = mapToGenome(preprocessRaw(patient))
  patient <- preprocessIllumina(patient)
  return(patient)
}

loadAndProcessControls <- function(control = parameters$controlFolder){
  control<-read.metharray.exp('/Users/lennartkester/Documents/EPIC/control/')
  gmSet = mapToGenome(preprocessRaw(control))
  control <- preprocessIllumina(control)
  data("exclude_regions")
  data("detail_regions")
  anno <- CNV.create_anno(array_type = "EPIC", exclude_regions = exclude_regions,detail_regions = detail_regions, chrXY = TRUE)
  refData <- list("control" = control, "anno" = anno)
  return(refData)
}

segmentData <- function(patient,refData){
  minfi.data <- CNV.load(patient) 
  minfi.controls <- CNV.load(refData$control)
  patient <- mapToGenome(patient)
  refData$anno@probes <- subsetByOverlaps(refData$anno@probes, granges(patient))
  x <- CNV.fit(minfi.data, minfi.controls, refData$anno) 
  x <- CNV.bin(x) 
  x <- CNV.detail(x) 
  x <- CNV.segment(x)
  igvData <- as.data.frame(cbind(as.character(anno@bins@seqnames),anno@bins@ranges@start,anno@bins@ranges@width+anno@bins@ranges@start,names(x@bin$ratio),x@bin$ratio))
  headerLine <- paste0("#track type=bedGraph name=",sampleNames(patient)," description=center_label visibility=full autoScale=off graphType=points viewLimits=-2:2 windowingFunction=none smoothingWindow=off")
  file <- paste0(parameters$dataFolder,sampleNames(patient),"/",sampleNames(patient),".igv")
  cat(headerLine, '\n',  file = file)
  write.table(igvData, file, append = T, col.names = F,row.names = F,quote = F,sep="\t")
  return(x)
}


