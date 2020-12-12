library("ggplot2")
library(caTools)

wus="Solyc02g083950.3.1"   # WUS
wox9="Solyc02g077390.2.1"   # WOX9
j2="Solyc12g038510.2.1" # j2

#meribroad = read.table("../peakData/ATAC-AN-meristem-broad_peaks_relative.bed")
meriPE = read.table("../peakData/ATAC-AN-meristem-PE_peaks_relative.bed")
#leafbroad = read.table("../peakData/ATAC-Leaf-broad_peaks_relative.bed")
#leafPE = read.table("../peakData/ATAC-Leaf-PE_peaks_relative.bed")

lbd29="Solyc09g066270.3.1"


#peaks = list(meribroad, meriPE, leafbroad, leafPE)


filelist = dir(pattern = "*.mat$", path = ".")
genelist = substr(filelist,1,18)
#################################
#  loadConservationMat
#
loadConservationMat <- function(locusfile, distance_correction = FALSE) {
  conse = as.matrix(read.csv(locusfile, row.names = 2, header=FALSE,as.is=TRUE))
  genomes = as.matrix(read.csv("~/workspace/solgenomes/allgenomes/genome_database.csv", header=TRUE, row.names=1, as.is=TRUE))
  
  len = nchar(conse[1,3])
  conse_mat = matrix(nrow=nrow(conse), ncol=len, data=0)
  rownames(conse_mat) = rownames(conse)
  
  for(i in rownames(conse)) {
    if(distance_correction) {
      conse_mat[i,] = as.numeric(unlist(strsplit(conse[i,3],""))[1:len])*as.numeric(genomes[i,"Distance"])
    } else {
      conse_mat[i,] = as.numeric(unlist(strsplit(conse[i,3],""))[1:len])
    }
  }
  
  mode(conse_mat)="numeric"
  conse_mat[is.na(conse_mat)]=0
  
  return(conse_mat)
}

#################################
# PlotConservation
#
# conse_mat - matrix of conservation
# peaks - list
plotConservation <- function(file, locus, peaks, distance_correction = TRUE) {
  mat <- suppressWarnings(loadConservationMat(file, distance_correction))
  
  runningaverage = runmean(colSums(mat),10)
  
  min_start = min(which(runningaverage>0))
  
  #trim 10%
#  min_start = min_start + round(length(runningaverage)*0.1)
  
#  plot(runningaverage,type="l", main = locus, xlim=c(min_start, length(runningaverage)), col="black",lwd=1, ylim=c(0,32))
  plot(smooth(runningaverage, twiceit=TRUE),type="l", main = locus, xlim=c(min_start, length(runningaverage)), col="black",lwd=1, ylim=c(0,32))
  
#  runningaverage_mean = runmean(colSums(mat),50)
#  lines(smooth(runningaverage_mean, twiceit = TRUE), col="black", lwd=1)
  
  level =max(colSums(mat)) - length(peaks)
  for(i in peaks) {
    
    locuspeaks = ncol(mat)+ i[i[,1]==locus,3:4]
    for(curpeak in 1:nrow(locuspeaks)) {
      lines(locuspeaks[curpeak,], c(level, level), lwd=5, col="red")
    }
    level = level+1
  }
}

plotConservationDecay <- function(locus,peaks) {
  mat <- loadConservationMat(locus)
  
  for(i in peaks) {
    locuspeaks = ncol(mat)+ i[i[,1]==locus,3:4]
    for(curpeak in 1:nrow(locuspeaks)) {

    }
  }

}


#############################################################
####
### gene list to peak conservation matrix

generatePeakConservationMatrix <- function(peaks, flank=1000) {
	genedirectory = "~/workspace/solgenomes/genes/"
	peakconservation=list()
	size_mid=300
	size_flank=1000
	
	for(i in 1:nrow(peaks)) {
		genename = peaks[i,1]
		genefilename = paste0(genedirectory,genename,".mat")
		if(file.exists(genefilename)) {
			geneconservemat = loadConservationMat(genefilename, TRUE)
			if(nrow(geneconservemat)>4) {
				geneconservemat = colSums(geneconservemat)
				
				geneconservemat = geneconservemat / max(geneconservemat)
			
				peak_coords = 30000+ c(peaks[i,4], peaks[i,3])
				peak_coords[peak_coords>30000]=30000
				peak_coords[peak_coords<1]=1
				
				if(abs(peak_coords[2]-peak_coords[1]) > 100) {
#					print(peaks_coords)
					peak = approx(geneconservemat[ (peak_coords[1]):(peak_coords[2]) ], n=size_mid, method="constant")$y 
					peak5 = approx(geneconservemat[ (max(1,peak_coords[1]-flank)):(peak_coords[1])  ], n=size_flank, method="constant")$y
					peak3 = approx(geneconservemat[ (peak_coords[2]):min(30000,peak_coords[2]+flank)], n=size_flank, method="constant")$y
				
					peakconservation[[peaks[i,2]]] = c(peak5,peak,peak3)
					cat(".")
				}
			}
		}
	}
	do.call(rbind,peakconservation)
}


generatePeakConservationMatrixNoCompress <- function(peaks, size=1000) {
	genedirectory = "~/workspace/solgenomes/genes/"
	peakconservation=list()
	
	for(i in 1:nrow(peaks)) {
		genename = peaks[i,1]
		genefilename = paste0(genedirectory,genename,".mat")
		if(file.exists(genefilename)) {
			geneconservemat = loadConservationMat(genefilename, TRUE)
			if(nrow(geneconservemat)>4) {
				geneconservemat = colSums(geneconservemat)
				geneconservemat = geneconservemat / max(geneconservemat)
				geneconservemat = c(rep(0, size), geneconservemat, rep(0,size))
			
				peak_coords = 30000+ c(peaks[i,4], peaks[i,3])+size
				peak_coords[peak_coords>(size+30000)]=30000+size
				peak_coords=sort(peak_coords, decreasing=FALSE)
				
				if(abs(peak_coords[2]-peak_coords[1]) > 100 && abs(peak_coords[2]-peak_coords[1])<size) {

					midpeak = mean(peak_coords)
					
					peak = geneconservemat[(midpeak-size/2):(midpeak+size/2) ]

					if(sum(peak)==0) {
						print(genefilename)
						print(peak_coords)
						print(midpeak)
						print(peaks[i,])
						print(length(peak))
					}
				
					peakconservation[[as.character(peaks[i,2])]] = peak
					cat(".")
				}
			}
		}
	}
	
	peaksmat = do.call(rbind,peakconservation)
	rownames(peaksmat) = names(peakconservation)
	peaksmat

}

generatePeakConservationForPosition <- function(peaks, start=27000, end=30000) {
	genedirectory = "~/workspace/solgenomes/genes/"
	peakconservation=list()
	
	for(i in 1:nrow(peaks)) {
		genename = peaks[i,1]
		genefilename = paste0(genedirectory,genename,".mat")
		if(file.exists(genefilename)) {
			geneconservemat = loadConservationMat(genefilename, TRUE)
			if(nrow(geneconservemat)>4) {
				geneconservemat = colSums(geneconservemat)
				geneconservemat = geneconservemat / max(geneconservemat)
		
				peak = geneconservemat[start:end]
			
				peakconservation[[peaks[i,2]]] = peak
					cat(".")
			}
		}
	}
	do.call(rbind,peakconservation)
}
