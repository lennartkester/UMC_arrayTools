
parameters <- list()
parameters$dataFolder <- "T:/pathologie/Moleculair/Diagnostiek/Methylatie array/Uitslagen Runfolder 2020/"
parameters$controlFolder <- "T:/pathologie/KMBP/EPIC/data/controlUMC/"

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

loadEPICsamples <- function(folder){
  files <- list.files(folder,recursive = F,full.names=F,pattern = ".idat")
  files <- sub("_Grn.idat","",files)
  files <- sub("_Red.idat","",files)
  return(unique(files))
}

loadAndProcessSamples <- function(dataFolder,sampleName){
  patient<-read.metharray(file.path(dataFolder,paste0(sampleName,"_Grn.idat")))
  gmSet = mapToGenome(preprocessRaw(patient))
  patient <- preprocessIllumina(patient)
  return(patient)
}

loadAndProcessControls <- function(controlFolder = parameters$controlFolder){
  control<-read.metharray.exp(controlFolder)
  gmSet = mapToGenome(preprocessRaw(control))
  control <- preprocessIllumina(control)
  data("exclude_regions")
  data("detail_regions")
  anno <- CNV.create_anno(array_type = "EPIC", exclude_regions = exclude_regions,detail_regions = detail_regions, chrXY = TRUE)
  refData <- list("control" = control, "anno" = anno)
  return(refData)
}

segmentData <- function(patient,refData,dataFolder){
  minfi.data <- CNV.load(patient) 
  minfi.controls <- CNV.load(refData$control)
  patient <- mapToGenome(patient)
  refData$anno@probes <- subsetByOverlaps(refData$anno@probes, granges(patient))
  x <- CNV.fit(minfi.data, minfi.controls, refData$anno) 
  x <- CNV.bin(x) 
  x <- CNV.detail(x) 
  x <- CNV.segment(x)
  igvData <- as.data.frame(cbind(as.character(refData$anno@bins@seqnames),refData$anno@bins@ranges@start,refData$anno@bins@ranges@width+refData$anno@bins@ranges@start,names(x@bin$ratio),x@bin$ratio))
  headerLine <- paste0("#track type=bedGraph name=",sampleNames(patient)," description=center_label visibility=full autoScale=off graphType=points viewLimits=-1:1 windowingFunction=none smoothingWindow=off")
  file <- paste0(dataFolder,"\\",sampleNames(patient),".igv")
  cat(headerLine, '\n',  file = file)
  write.table(igvData, file, append = T, col.names = F,row.names = F,quote = F,sep="\t")
  return(x)
}


