#Used librarys
library("gplots") #function rich.colors
library("preprocessCore") #function normalize.quantiles
	
###############
###Functions###				
###############

###+++###
#Function requires data frame to be normalized
#1. All 0 are set to NA, to exclude them from quatile normalization
#2. Data are quantile normalized
#3. 0 values (the one set to NA) are set back to 0
fQN <- function(x) #
{
	x[x==0] <- NA
	x_m <- as.matrix(x)
	x <- normalize.quantiles(x_m)
	x[is.na(x)] <- 0
	return(data.frame(x))
}	
###***###***###

###+++###	
#Function require a vector with expression of one gene in different tissues.
#Mean is calculated taking in account tissues with 0 expression. 2+0+4=2
fmean <- function(x)
	{
		if(!all(is.na(x)))
	 	{
	 		res <- mean(x, na.rm=TRUE)
	 	} else {
	 		res <- NA
	 	}
	 	return(res)
	}
###***###***###	

###+++###	
#Function require a vector with expression of one gene in different tissues.
#Max is calculated taking in account tissues with 0 expression. 2+0+4=2
fmax <- function(x)
	{
		if(!all(is.na(x)))
	 	{
	 		res <- max(x, na.rm=TRUE)
	 	} else {
	 		res <- NA
	 	}
	 	return(res)
	}
###***###***###	

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Minimum 2 tissues
fTau <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(max(x)!=0)
 			{
 				x <- (1-(x/max(x)))
 				res <- sum(x, na.rm=TRUE)
 				res <- res/(length(x)-1)
 			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		} 
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	} 
 	return(res)
}
###***###***###

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fGini <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x!=0))
 			{
 				res <- gini(x)*(length(x)/(length(x)-1))
  			} else {
 				res <- 0
 			}
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}	 		
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fTsi <- function(x)
{
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x!=0))
 			{
 				res <- max(x) / sum(x)
 			} else {
 				res <- 0
 			} 	
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}	
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###	

###+++###
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Function requires setting of a treshold (rpkm)	
fCounts <- function(x, rpkm)
{
	if(all(!is.na(x)))
 	{
 		res <- length(which(x > rpkm))	
 		if (res > 0)
 		{
 			res <- (1 - res/length(x))*(length(x)/(length(x)-1))  #Modification: To bring to normalized scale
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###		

###+++###
#Function require a data frame with expression data, and give back a vector with EEi values for each gene
#If expression for one tissue is not known, gene specificity for this gene is NA
fEe <- function(x)
{
	if(!all(is.na(x)))
 	{
 		x <- as.matrix(x)
 		x[x<0] <- NA
 		x <- cbind(x, r=rowSums(x, na.rm=FALSE))
		x <- rbind(x, c=colSums(x, na.rm=TRUE))	
 		x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)])
 						
 		res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max)
 		res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale
 	} else {
 		res <- NA
 		print("No data avalable.")
 	}
  	return(res)
}
###***###***###

###+++###
#Function require a data frame with expression data, and give back a vector with PEM scores
#If expression for one tissue is not known, gene specificity for this gene is NA
fPem <- function(x)
{
	if(!all(is.na(x)))
 	{
 		x <- as.matrix(x)
 		x[x<0] <- NA
 		x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
		x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
 		x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score
 		
 		x[x<1] <- 1
 		x <- log10(x)
 		 		
 		x<- abs(x)				
 		res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max) #choose only the maximal score for each gene
 		res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
 	} else {
 		res <- NA
 		print("No data avalable.")
 	}
  	return(res)
}
###***###***###

###+++###
#Hg entropy
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fHg <- function(x)
{	
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{
 			if(sum(x) !=0)
 			{
 				p <- x / sum(x)
 				res <- -sum(p*log2(p), na.rm=TRUE)
 				res <- 1 - (res/log2(length(p))) #Modification: To bring to normalized scale
 			} else {
 				res <- 0
 			} 		
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###	

###+++###
#Z-score
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fZ <- function(x)
{	
	if(all(!is.na(x)))
 	{
 	res <-  apply(scale(t(x), center=TRUE, scale=TRUE),2,max)/((length(x[1,])-1)/sqrt(length(x[1,])))
 	res[is.na(res)] <- 0
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###	

###+++###
#SPM score from TISGED
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
fSpm <- function(x)
{	
	if(all(!is.na(x)))
 	{
 		if(min(x, na.rm=TRUE) >= 0)
		{	
 			if(sum(x) !=0)
 			{
 				spm <- x^2/(x%*%x)
 				res <- max(spm) #Modification:To bring to normalized scale. Choose max
 			} else {
 				res <- 0
 			}	 		
 		} else {
 		res <- NA
 		#print("Expression values have to be positive!")
 		}
 	} else {
 		res <- NA
 		#print("No data for this gene avalable.")
 	}
 	return(res)
}
###***###***###

###+++###
#Calculate and save tissue specificity parameters
#orgExpression = data set, rpkm = cutt off, add = number of tissues, tNames = tissues to use, tNamesNew = tissues to name, RNAseq = how to normalise (log_QN, QN, log, NA)
#Only genes with Ensembl IDs are used, or for Drosophila
#Normalization is done on all tissues, not dependent which tissues are selected later
#1. Data are normalized 
#2. All expression under rpkm is set to 0
#3. Replicates mean is calculated (fReplicateMean)
#4. Genes that not expressed in any tissue are removed
#5. Tissue specificity parameters are calculated
fTS <- function(orgExpression, rpkm, add, tNames, tNamesNew, RNAseq)
{	
	orgExpression <- orgExpression[regexpr("ENS", orgExpression$Ensembl.Gene.ID)>0 | regexpr("FBgn", orgExpression$Ensembl.Gene.ID)>0 | regexpr("PPAG", orgExpression$Ensembl.Gene.ID)>0, ]
	orgExpression <- na.omit(orgExpression)
	print(summary(orgExpression))
	if(RNAseq == "log_QN"){
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 1
		x <- log2(x)
		rpkm <- log2(rpkm)
		orgExpression[,c(-1)] <- fQN(x)
	} else if (RNAseq == "QN")  {
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 0
		orgExpression[,c(-1)] <- fQN(x)
	} else if (RNAseq == "log")  {
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 1
		orgExpression[,c(-1)] <- log2(x)
		rpkm <- log2(rpkm)
	} else {
		x <- orgExpression[,c(-1)]
		x[x < rpkm] <- 0
		orgExpression[,c(-1)] <- x
	}
	
	orgExpression <- fReplicateMean(orgExpression, expDataSource, organism, paste("Averaged.RPKM.",tissuesNames, sep=""))
	orgExpression$Max <- apply(orgExpression[,c(-1)], c(1), fmax)
	orgExpression <- orgExpression[orgExpression$Max > rpkm,]
	orgExpression <- orgExpression[,c(-length(colnames(orgExpression)))]
	print(summary(orgExpression))
	fPlotExpression(orgExpression[,-1], paste("Normalized expression (cutoff", 2^rpkm, "RPKM)", sep=" "), paste("NormalizedQN_", 2^rpkm,"RPKM", sep=""), tissuesPrintNames)
	
	orgExpression <- orgExpression[,c("Ensembl.Gene.ID", paste("Averaged.RPKM.", tNames,sep="")) ]
	colnames(orgExpression) <- c("Ensembl.Gene.ID", paste("Averaged.RPKM.", tNamesNew, sep=""))
	nTissues <- length(tNamesNew)
	tissuesNames <- tNamesNew
	print(paste("Analysis done on", nTissues, "tissue:", sep=" "))
	print(tissuesNames)

	orgExpression$Tau <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fTau)
	orgExpression$Gini <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fGini)
	orgExpression$Tsi <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fTsi)
	orgExpression$Counts <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, function(x){x <- fCounts(x, rpkm)})
	orgExpression$Hg <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fHg)
	orgExpression$Zscore <- fZ(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))])
	orgExpression$Spm <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fSpm)
	orgExpression$Ee <- fEe(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))])
	orgExpression$Pem <- fPem(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))])
	
	orgExpression$Mean <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fmean)
	orgExpression$Max <- apply(orgExpression[,c(paste("Averaged.RPKM.", tissuesNames[1:nTissues], sep=""))], 1, fmax)

	print(summary(orgExpression))
	
	p <- c("Tau", "Gini", "Tsi", "Counts", "Ee", "Hg", "Zscore", "Spm", "Pem")
	x <- as.matrix(orgExpression[,p])
	xs <- cor(x, method="spearman")
	xp <- cor(x, method="pearson")
	capture.output(c("Spearman correlation"),file=paste(folder, organism, expDataSource,"CorrelationTS_", add, ".txt", sep="")) 
	capture.output(xs, append=TRUE, file=paste(folder, organism, expDataSource, "CorrelationTS_", add, ".txt", sep="")) 
	capture.output(c("Pearson correlation"), append=TRUE, file=paste(folder, organism, expDataSource,"CorrelationTS_", add, ".txt", sep="")) 
	capture.output(xp, append=TRUE, file=paste(folder, organism, expDataSource,"CorrelationTS_", add, ".txt", sep="")) 

	dev.new(height=9, width=12)
		par(cex.main=0.95, bg=my.col[1], fg=my.col[2], col.axis=my.col[2], col.lab=my.col[2], col.main=my.col[2])
		palette(rev(rich.colors(10)))
		#palette(rev(blues9))
	
		plot(density(orgExpression[,"Tau"],n=1000), main = " ", xlab="Tissue specificity",col=(1), lwd=4, lty=1
		,ylim=c(0,8), xlim=c(-0.1,1.1)
		)
		lines(density(orgExpression[,"Gini"],n = 1000), col=(2), lwd=4, lty=2)
		lines(density(orgExpression[,"Tsi"],n = 1000), col=(3), lwd=4, lty=1)
		lines(density(orgExpression[,"Counts"],n = 1000), col=(4), lwd=4, lty=2)
		lines(density(orgExpression[,"Ee"],n = 1000), col=(5), lwd=4, lty=1)
		lines(density(orgExpression[,"Hg"],n = 1000), col=(6), lwd=4, lty=2)
		lines(density(orgExpression[,"Zscore"],n = 1000), col=(7), lwd=4, lty=1)
		lines(density(orgExpression[,"Spm"],n = 1000), col=(8), lwd=4, lty=2)
		lines(density(orgExpression[,"Pem"],n = 1000), col=(9), lwd=4, lty=1)
				
		legend("topright",c("Tau", "Gini", "TSI", "Counts", "EE", "Hg", "Zscore", "SPM", "PEM"),col=(1:11), lwd=4, lty=c(1,2), bty="n", seg.len=4)
		
		dev.copy2pdf(device=quartz, file=paste(folder, organism, expDataSource, "TScomparison_9_",  add,".pdf", sep=""),onefile=TRUE)#,paper="A4r"
		#dev.off()

	write.table(orgExpression, file=paste(folder, organism, expDataSource,"TScomparisonTable_9_",  add,".txt",sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE)
		
	fScatPlot(orgExpression, "Tau", c("Gini", "Tsi", "Counts", "Hg", "Zscore","Spm", "Ee", "Pem"), add, c(0, 0.94, 4, 0.5))
	fScatPlot2(orgExpression, "Mean", c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem"), add, c(ceiling(max(orgExpression$Mean)), 0.95, 2, 0.5), ceiling(max(orgExpression$Mean)))
	fScatPlot2(orgExpression, "Max", c("Tau", "Gini", "Tsi", "Counts", "Hg", "Zscore", "Spm", "Ee", "Pem"), add, c(ceiling(max(orgExpression$Max)), 0.95, 2, 0.5), ceiling(max(orgExpression$Max)))

	return()
}
###***###***###

