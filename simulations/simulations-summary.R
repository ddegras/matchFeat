

# This code produces Table 1, Figure 1, and Figure 2 of
# "Scalable Feature Matching Across Large Data Collections" 
# by D. Degras (2021)

# The  file simulations.R should be run before simulations_summary.R
# and both should be in the same folder. 



########################
# Create summary tables
########################


## Clear workspace
rm(list=ls())

## Get names of files containing simulation results
fname <- list.files()
fname <- grep("sim_results", fname, value=TRUE)
load(fname[1])

## Summary dimensions
nsetup <- length(fname)
nval <- integer(nsetup)
method <- dimnames(result)[[1]]
nmethod <- length(method)
stat <- c("mean","median","sd","mad")
nstat <- length(stat)
error <- array(dim=c(nmethod,2,nstat,nsetup)) 
dimnames(error) <- list(method=method, 
	type=c("absolute","relative"), stat=stat)
timing <- array(dim=c(nmethod, nstat, nsetup))
dimnames(timing) <- list(method=method, stat=stat)
Rand <- array(dim=c(nmethod, nstat, nsetup))
dimnames(Rand) <- list(method=method, stat=stat)

for (i in 1:nsetup)
{
	## Load results
	load(fname[i])
	nval[i] <- n
	nrep <- dim(result)[3]
	tmp <- array(dim=c(dim(result)[1],4,nrep))
	
	## Calculate absolute and relative errors
	for (j in 1:nrep) {
		m <- min(result[,1,j], na.rm=T)
		tmp[,1,j] <- result[,1,j] - m 
		tmp[,2,j] <- tmp[,1,j] / m }
	tmp[tmp < 1e-10] <- 0
	
	## Add Rand index and timing info
	tmp[,3:4,] <- result[,2:3,]

	## Calculate summary stats
	out <- array(dim=c(dim(result)[[1]],4,nstat))
	summary.fun <- list(mean,median,sd,mad)
	for (j in 1:4)
	for (k in 1:nstat)
	out[,j,k] <- apply(tmp[,j,],1,summary.fun[[k]],na.rm=TRUE)	
	
	## Store
	method.tmp <- dimnames(result)[[1]]
	method.tmp[method.tmp == "kur.nqp"] <- "kur.iqp"
	error[method.tmp,,,i] <- out[,1:2,]
	Rand[method.tmp,,i] <- out[,3,]
	timing[method.tmp,,i] <- out[,4,]
}

dimnames(error) <- c(dimnames(error)[1:3], list(n=nval))
dimnames(Rand) <- c(dimnames(Rand)[1:2], list(n=nval))
dimnames(timing) <- c(dimnames(timing)[1:2], list(n=nval))

## Reorder by increasing n
error <- error[,,,order(nval)]
Rand <- Rand[,,order(nval)]
timing <- timing[,,order(nval)]
nval <- sort(nval)






########################################################################





#############################
# Create LaTeX table for
# relative error (mean + sd)
#############################

if (!require(xtable)) install.packages("xtable")
library(xtable)

## Create LaTeX table 
o <- order(rowMeans(error[,"relative","mean",], na.rm=TRUE))
tab <- matrix(,nmethod,2*nsetup)
tab[,seq.int(1,2* nsetup,2)] <- error[o,"relative","mean",]
tab[,seq.int(2,2* nsetup,2)] <- error[o,"relative","sd",]
rownames(tab) <- method[o]
print(xtable(tab,digits=0, display=rep("E", ncol(tab)+1)), file="rel_error.tex")
txt <- readLines("rel_error.tex")

## Detect data rows
string <- paste("|",method[-1],collapse="", sep="")
string <- paste0(method[1],string)
idx <- grep(string, txt)

## Clean up 
txt <- txt[idx]
nline <- length(idx)
out <- character(nline)
for (i in 1:nline) {
	test <- gsub(" ", "", txt[i])
	test <- unlist(strsplit(test,"&",fixed=T))
	first <- toupper(test[1])
	last <- test[length(test)]
	last <- substr(last, 1,nchar(last)-2)
	test <- paste0(test[-c(1,length(test))], c(" (",")& "), collapse="")
	out[i] <- paste0(first, " &", test, last,")\\\\", collapse="")
}
out <- gsub(".","-",out,fixed=TRUE)
out <- gsub("LABEL","LBL",out)
out <- gsub("COPKM","COP-KM",out)
# out <- gsub("KMEANS","KMM",out)
# out <- gsub("r100","R100",out)
out <- gsub("()","",out,fixed=TRUE)
out <- gsub("E+0","E+",out,fixed=TRUE)
out <- gsub("E-0","E-",out)

## Embed in LaTeX table environment
top <- c("\\begin{table}[ht]", "\\centering", 
	paste0("\\begin{tabular}{r *{",nsetup,"}{c}}"),                                                                                                                                                             
	"\\hline", paste("", paste(" &",nval, collapse=""), "\\\\"), 
	"\\hline")
bottom <- c("\\hline", "\\end{tabular}", "\\end{table}")

## Save as .tex file
writeLines(c(top,out,bottom), "rel_error.tex")




########################################################################




##################
# Plot Rand index
##################

pdf("rand_index.pdf")

## Select and sort methods to plot
o <- order(rowMeans(Rand[,"mean",], na.rm=TRUE), decreasing=TRUE)
idx <- grep("fw|label", method[o], value = TRUE, invert = TRUE)

## Main plot
rand <- Rand[idx,"mean",]
rand[is.nan(rand)] <- NA
par(mar=c(3, 4, 1, 2) + 0.1, mgp=c(3, 0.75, 0))
plot(0, 0, xlim=c(0,1050), ylim=c(0.8,1), type="n", xlab="", ylab="")
mtext(side=1,text=expression(italic(n)),line=2, cex=1.15)	
mtext(side=2,text="Rand index",line=2.5, cex=1.05)	
dims <- par("usr")
rect(dims[1],dims[3],dims[2],dims[4],col="grey92")
matlines(nval, t(rand), lwd=1.25)
grid(col="white")

## Display text labels 
idx <- match(c("r100.bca", "r100.kmm", "rec", "lcvqe", "hub", 
	"kur.ilp", "id", "r100", "gaussmix", "copkm", "mpckm", "ccls"), 
	rownames(rand))
pos <- apply(rand[idx,], 1, function(r) max(which(!is.na(r))))
xtext <- nval[pos]
ytext <- rand[cbind(idx,pos)]
# Adjust offsets as needed to avoid label overlap
xoffset <- rep(0,length(xtext))
xoffset[length(xtext)] <- -100
yoffset <- 0.002 * c(1,-1,0,0.25,.75,-1,1,-1,0,0,0,45) 
labels <- c("BCA", "KMM", "REC", "LCVQE", "HUB", 
	"KUR-ILP", "ID", "R100", "EM", "COP-KM", "MPC-KM", "CCLS")
text(x = xtext + xoffset, y = ytext + yoffset, cex = 0.56, 
	labels = labels, pos = 4, offset = 0.2)

dev.off()




########################################################################




####################
# Plot running time
####################

pdf("timing.pdf")

## Select methods to plot 
rmv <- c(grep("2x|fw", method, value = TRUE), "ilp", "kur.iqp",
	"rec", "hub", "id", "label", "r100")
idx <- setdiff(method,rmv)

## Set up timing info 
tmng <- timing[idx,"mean",]
idx <- order(rowMeans(tmng, na.rm=TRUE))
tmng <- tmng[idx,]
for (i in 1:nrow(tmng))
{
	idx <- which(tmng[i,] >= 300)
	if (length(idx) > 0) tmng[i,idx[-1]] <- NA # replace timeouts by NA
}

## Main plot 
par(mar=c(3, 4, 1, 2) + 0.1, mgp=c(3, 0.75, 0))
plot(0,0, xlim=c(0,1075), ylim=range(tmng,na.rm=T), 
	type="n", xlab="", ylab="")
mtext(side=1,text=expression(italic(n)),line=2, cex=1.15)	
mtext(side=2,text="Running time (s)",line=2.5, cex=1.05)	
dims <- par("usr")
rect(dims[1],dims[3],dims[2],dims[4],col="grey92")
matlines(nval, t(tmng), lwd=1.25)
grid(col="white")

## Display text labels
idx <- match(c("r100.bca", "r100.kmm", "lcvqe", "kur.ilp", 
	"hub.bca", "hub.kmm", "gaussmix", "copkm", "mpckm", "ccls"), 
	rownames(tmng))
pos <- apply(tmng[idx,], 1, function(r) max(which(!is.na(r))))
xtext <- nval[pos]
xoffset <- 0
ytext <- tmng[cbind(idx,pos)]
yoffset <- 3 * c(0,0,0,0,1,-1,1,1,-1,-1) 
# adjust as needed to avoid label overlap
labels <- c("R100-BCA", "R100-KMM", "LCVQE", "KUR-ILP", "BCA", 
	 "KMM", "EM", "COP-KM", "MPC-KM", "CCLS")
text(x = xtext + xoffset, y = ytext + yoffset, cex = 0.56, labels = labels, 
	pos = 4, offset = 0.2)

dev.off()



