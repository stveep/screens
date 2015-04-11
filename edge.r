require("edgeR")

edgeets <- function(x,dmsocols,rexp="L.+_(.+)_.+_.+_.+") {
	if (names(x)[1] == "barcode") {
		x <- data.frame(x[,-1],row.names=x[,1])
	}
	if(missing(dmsocols)) {
		dmsocols = grep("DMSO",names(x))
	}
#	group = seq(from=2,to=(ncol(x) + 1))
#	group[dmsocols] <- 1
	group = gsub(rexp,"\\1",names(x))
	y <- DGEList(counts=x,group=group)
	y <- calcNormFactors(y)
	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	y
	lst <- list()
	for(i in group) {
		pr <- c("DMSO",as.character(i))
		et <- exactTest(y,pair=pr)
		lst[[i]] <- et
	}
	lst
}

dgeobjs <- function(x,dmsocols,rexp="L.+_(.+)_.+_.+_.+") {
        if (names(x)[1] == "barcode") {
                x <- data.frame(x[,-1],row.names=x[,1])
        }
        if(missing(dmsocols)) {
                dmsocols = grep("DMSO",names(x))
        }
#       group = seq(from=2,to=(ncol(x) + 1))
#       group[dmsocols] <- 1
        group = gsub(rexp,"\\1",names(x))
        y <- DGEList(counts=x,group=group)
        y <- calcNormFactors(y)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        y
}

edgean <- function(x,dmsocols,rexp="L.+_(.+)_.+_.+_.+") {
	if (names(x)[1] == "barcode") {
		x <- data.frame(x[,-1],row.names=x[,1])
	}
	if(missing(dmsocols)) {
		dmsocols = grep("DMSO",names(x))
	}
#	group = seq(from=2,to=(ncol(x) + 1))
#	group[dmsocols] <- 1
	group = gsub(rexp,"\\1",names(x))
	y <- DGEList(counts=x,group=group)
	y <- calcNormFactors(y)
	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	y
	df <- y$pseudo.counts
	for(i in group) {
		pr <- c("DMSO",as.character(i))
		et <- exactTest(y,pair=pr)
		et <- et$table[,c("logFC","PValue")]
		names(et) <- sapply(names(et),function(a,j) paste(j,a,sep="_"),j=i)
		df <- cbind(df,et)
	}
	df
}

smearSens <- function(lib,drug,ntag) {
plotSmear(edgelist[[lib]][[drug]],de.tags=row.names(topTags(sens[[lib]][[drug]],ntag)),cex=0.5,main=paste(lib,drug,sep=" "))
}

dhits <- function(etlst,drug,pt=0.001,fct=-1,sens=TRUE,vthr=1000) {
	# Get library names that have results for the drug
	sublist <- drugInLib(etlst,drug)
	# Discount those with high variance
	sublist <- sublist[sapply(sublist, function(x,etlst) {tbl <- etlst[[x]][[drug]]$table; var(tbl[tbl$logCPM > 4,"logFC"]) < vthr},etlst=etlst)]	
	sublist <- etlst[sublist]
	lapply(sublist,function(x,d) hits(x[[d]],pthr=pt,fcthr=fct,sens=sens)$table,d=drug)
}

dzhits <- function(etlst,drug,zt=-2,cpmt=4,sens=TRUE,vthr=1000) {
	# Get library names that have results for the drug
	sublist <- drugInLib(etlst,drug)
	# Discount those with high variance
	sublist <- sublist[sapply(sublist, function(x,etlst) {tbl <- etlst[[x]][[drug]]$table; var(tbl[tbl$logCPM > 4,"logFC"]) < vthr},etlst=etlst)]	
	sublist <- etlst[sublist]
	lapply(sublist,function(x,d) zhits(x[[d]],zt,cpmt,sens=sens)$table,d=drug)
}

dhitgenes <- function(etlst,maplst,drug,fn=dhits,mincount=1,sens=TRUE,vthr=1000,...) {
	unmappedhits <- fn(etlst,drug,sens=sens,vthr=vthr,...)
	table(unlist(lapply(names(unmappedhits),function(x) as.character(maplst[[x]][maplst[[x]]$barcode %in% row.names(unmappedhits[[x]]),"gene"]))))[
	table(unlist(lapply(names(unmappedhits),function(x) as.character(maplst[[x]][maplst[[x]]$barcode %in% row.names(unmappedhits[[x]]),"gene"])))) > mincount]
}

hits <- function(et,pthr=0.001,fcthr=-1,sens=TRUE) {
	if (sens==TRUE) {
	et[et$table$logFC < fcthr & et$table$PValue < pthr,]
	} else {
	et[et$table$logFC > fcthr & et$table$PValue < pthr,]
	}
}

zhits <- function(et,zt=-2,cpmt=4,sens=TRUE) {
	if (sens==TRUE) {
	et[et$table$z < zt & et$table$logCPM > cpmt,]
	} else {
	et[et$table$z > zt & et$table$logCPM > cpmt,]
	}
}
drugInLib <- function(lst,drug) {
names(lst[unlist(lapply(lst,function(x) drug %in% names(x)))])
}

# Function for volcano plots

volcs <- function(df,drug,end="") {
	pcol <- paste(drug,"_PValue",sep="")
	fccol <- paste(drug,"_logFC",sep="")
	plot(df[,fccol],log10(df[,pcol]),pch=19,main=paste(drug,end),xlab="Fold enrichment/depletion",ylab="P-value (log10)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

}

vplots <- function(dfs,dfw,xs,xw,filen="volcano.pdf",rexp="L.+_(.+)_.+_.+_.+") {
	groups = unique(gsub(rexp,"\\1",names(xs[,-1])))
#	groups = gsub("^(\\d)","X\\1",groups)
	groupw = unique(gsub(rexp,"\\1",names(xw[,-1])))
#	groupw = gsub("^(\\d)","X\\1",groupw)
	# Second sub was needed because an X is added to col. names of drugs that start with a number (e.g. 17AAG), apparently now no longer nec.
	pdf(filen)
	for(i in groups) {
		volcs(dfs,i,"sims")
	}
	for(i in groupw) {
		volcs(dfw,i,"west")
	}
	dev.off()
}
