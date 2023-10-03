if (!require(HDInterval, quietly = T)) install.packages("HDInterval")
library(HDInterval)
if (!require(graphics, quietly = T)) install.packages("graphics")
library(graphics)
if (!require(shape, quietly = T)) install.packages("shape")
library(shape)
if (!require(emdbook, quietly = T)) install.packages("emdbook")
library(emdbook)


HILL = FALSE


SF = 3

# Get the 95% credible interval of a parameter, returned as a string ready to print
get.95.str <- function(a.vector, units = "") {

	if (a.vector[1] == "-"){
		return("-")
	}

	res = hdi(a.vector, cresMass=.95)
	res = c(median(a.vector), res[[1]], res[[2]])
	res = signif(res, SF)
	if (nchar(units) > 0){
		units = paste0(" ", units)
	}
	return(paste0(res[1], units, " (", res[2], ", ", res[3], ")"))
}



BURNIN = 0.2
SHOW.INFERENCE = TRUE

HETERO.LINE.COL = "#1E88E5"
LINE.COL = "#FFC107"

DASHED.LINE.COL = "#696969"
LOGX = T

# Read in log file
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1){
	stop("Use: Rscript plotMM.R <log file> <burnin percentage>")
}


log.df = read.table(args[1], header=T, sep="\t", comment="#")
if (length(args) > 1){
	BURNIN = as.numeric(args[2]) / 100
}
log.df = log.df[ (floor(nrow(log.df) * BURNIN)):nrow(log.df), ]
outfile = paste0(args[1], ".png")


cat(paste("Saving to", outfile, "\n"))
png(outfile, width = 800, height = 800, res=130)
par(family = "Helvetica Neue Light")
par(mar = c(4,4,4,1))




# Get reactant concentrations
a = log.df[1,grep("a[.]", colnames(log.df))]



# Get reaction velocities
hasVobs = length(grep("v[.]obs[.]", colnames(log.df))) > 0
if (hasVobs){
	vobs = log.df[1,grep("v[.]obs[.]", colnames(log.df))]
}
bpred = log.df[,grep("v[.]pred[.]", colnames(log.df))]
bsample= log.df[,grep("v[.]sample[.]", colnames(log.df))]
if (hasVobs){
	b = bpred
}else{
	b = bsample
}

modelIndicator = names(sort(table(log.df[,"modelIndicator"]), decreasing=T))[1]
modelIndicator = ifelse(modelIndicator == 1, "Heterogeneous", "Homogeneous")

hillIndicator = names(sort(table(log.df[,"hillIndicator"]), decreasing=T))[1]
hillIndicator = ifelse(hillIndicator == 1, "Hill", "")

modelDesc = paste(modelIndicator, hillIndicator)


# Do not plot duplicates in b
is.duplicate = sapply(1:length(a), function(i){


	if (i == 1){
		FALSE
	}else{
		ele = as.numeric(a[i])
		#print(ele)
		if (any(a[1:(i-1)] == ele)){
			TRUE
		}else{
			FALSE
		}

	}
	

})


nobs = length(a)
nsamples = nrow(bpred)







xmax = max(a)
ymax = max(bpred)*1.2


xmax = max(a) * 1.5
xmin = 0
if (LOGX){
	xmin = min(a)/2
}


# Plot
if (LOGX){
	plot(1,1, log="x", type="n", xlim = c(xmin, xmax), ylim = c(0, ymax), xlab = "a (log space)", ylab = "r", main = "HetMM", axes=F, xaxs="i", yaxs="i")
}else{
	plot(0,0, type="n", xlim = c(xmin, xmax), ylim = c(0, ymax), xlab = "a (log space)", ylab = "r", main = "HetMM", axes=F, xaxs="i", yaxs="i")
}





legend("bottomright", c("Homogeneous", "Heterogeneous"), col=c(LINE.COL, HETERO.LINE.COL), lwd=2, bty="n")



if (SHOW.INFERENCE){

	mtext(paste("Best model:", modelDesc), cex=0.7)



	# Plot predicted
	for (i in 1:nsamples){
		m = log.df[i,"modelIndicator"]
		col = ifelse(m == 1, paste0(HETERO.LINE.COL, "55"), paste0(LINE.COL, "33"))
		lines(a[!is.duplicate], b[i,!is.duplicate], col = col)
	}

}

# Plot observed
if (hasVobs){
	points(a, vobs, cex=1.2, pch=21, col="black", bg="#dddddd")

	# R2
	lm = lm(as.numeric(vobs) ~ as.numeric(apply(b, 2, mean)))
	x = summary(lm)
	r2 = x$r.squared
	r2 = signif(r2, 3)



	#legend("bottomright", c(paste0("R2=", r2), paste0("DG=", gibbs)), bty="n")
	legend("topleft", c(paste0("R²=", r2)), bty="n")
}





axis(1)
axis(2, las=2)


dev.off()


# Print parameters
vmax = log.df[,"vmax"]
errorSD = log.df[,"errorSD"]
km = log.df[,"km"]
kmSD = log.df[,"kmSD"]
hill = log.df[,"hill"]
phetero = signif(mean(log.df[,"modelIndicator"]), SF)
phill = signif(mean(log.df[,"hillIndicator"]), SF)

if (phetero > 0.95 | phetero < 0.05){
	phetero = paste0("\\textbf{", phetero, "}")
}
if (phill > 0.95 | phill < 0.05){
	phill = paste0("\\textbf{", phill, "}")
}


if (modelIndicator == "Homogeneous"){
	kmSD = "-"
}else{
	kmSD = kmSD[log.df$modelIndicator == 1]
}


if (hillIndicator != "Hill"){
	hill = "-"
}else{
	hill = hill[log.df$hillIndicator == 1]
}


cat(paste("Summarising posterior distribution:\n"))
cat(paste("median (2.5 percentile, 97.5 percentile)\n"))
cat(paste("\t      Vmax =", get.95.str(vmax), "\n"))
cat(paste("\t        Km =",   get.95.str(km), "\n"))
cat(paste("\t         ε =",get.95.str(errorSD), "\n"))
if (kmSD != "-"){
	cat(paste("\t      kmSD =", get.95.str(kmSD), "\n"))
}
if (hill != "-"){
	cat(paste("\t      Hill =", get.95.str(hill), "\n"))
}
cat(paste("\t p(hetero) =", phetero, "\n"))
cat(paste("\t   p(Hill) =", phill, "\n"))
cat(paste("\tBest model =", modelDesc, "\n"))




