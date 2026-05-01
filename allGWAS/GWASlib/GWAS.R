args<-commandArgs(trailingOnly=TRUE);
library(foreach)
library(doParallel)


getcolname <- function(o, x, y, z) {
	outtmp <- o
	for (i in 1:length(x)) {
		if (z == "CDS") {
			iname <- c( paste(x[i],"Coefficient",sep="_"), paste(x[i],"Zvalue",sep="_"), paste(x[i],"Pvalue",sep="_") )
		} else if (z == "SNP") {
			iname <- c( paste(x[i],"Coefficient",sep="_"), paste(x[i],"Pvalue",sep="_"), paste(x[i],"Zvalue",sep="_") )
		} else if (z == "nonCDS") {
			iname <- c( paste(x[i],"Pvalue",sep="_"), paste(x[i],"Coefficient",sep="_"), paste(x[i],"Zvalue",sep="_") )
		}
		outtmp <- c(outtmp, iname)
	}
	anovatmp <- paste(y,rep("anova_Pvalue",length(y)),sep="_")
	outtmp <- c(outtmp, anovatmp)	
	outtmp <- c(outtmp, c("R2_McFadden","R2_Cox.Snell","R2_Nagelkerke","R2_McKelvey.Zavoina","R2_Effron","R2_Count"))
	return(outtmp)
}

ptype <- args[4];
dat <- read.table(args[1],as.is=TRUE,header=FALSE,sep="\t");
if (ptype == "CDS") {
	mat <- dat[,2:ncol(dat)]
} else if (ptype == "SNP") {
	mat <- dat[,3:ncol(dat)]
} else if (ptype == "nonCDS") {
	mat <- dat[,4:ncol(dat)]
}

rowNum <- nrow(mat);
cols <- unlist(strsplit(args[2], split=','));
ncols <- length(cols); 
types <- unlist(strsplit(args[3], split=','));

mattmp <- data.frame(matrix(as.vector(unlist(mat[1, ])), nrow=ncol(mat)/ncols, ncol=ncols, byrow=T), stringsAsFactors=FALSE);
ColNum <- 2;
keyname <- c("Allele1","Allele2");
keyname2 <- cols[1:(ncols-1)];
if (types[1] != "N") {
	for (j in 2:(ncols-1)) {
		if (types[j-1] == "f") {
			mattmp[, j] <- as.factor(as.character(mattmp[, j]))
			ColNum <- ColNum + length(levels(mattmp[, j]))
			keyname = c(keyname, paste(rep(cols[j],length(levels(mattmp[, j]))),levels(mattmp[, j]),sep=""))
		} else if (types[j-1] == "n") {
			ColNum <- ColNum + 1
			keyname = c(keyname, cols[j])
		} else {
			stop("Please provide right data types for the covariates!")
		}
	}
}
rm(mattmp);

if (ptype == "CDS") {
	ColNum <- 1 + 3 + ColNum * 3 + ncols - 1 + 6
	tmp <- c("Site","Intercept_Pvalue","Null_deviance","Residual_deviance")
	Colnames <- getcolname(tmp, keyname, keyname2, args[4])
} else if (ptype == "SNP") {
	ColNum <- 2 + 3 + ColNum * 3 + ncols - 1 + 6
	tmp <- c("RefName","Site","Intercept_Pvalue","Null_deviance","Residual_deviance")
	Colnames <- getcolname(tmp, keyname, keyname2, args[4])
} else if (ptype == "nonCDS") {
	ColNum <- 3 + 3 + ColNum * 3 + ncols - 1 + 6
	tmp <- c("RefName","Site","Sub","Intercept_Pvalue","Null_deviance","Residual_deviance")
	Colnames <- getcolname(tmp, keyname, keyname2, args[4])
}

res_all <- data.frame(matrix(0, rowNum, ColNum), stringsAsFactors=FALSE);

LogisRre <- function(i) {
	mat_ii <- data.frame(matrix(as.vector(unlist(mat[i,])), nrow=ncol(mat)/ncols, ncol=ncols, byrow=T), stringsAsFactors=FALSE);
	check <- FALSE;
	mat_ii[,ncols] <- as.factor(as.character(mat_ii[,ncols])); 
	phenof <- factor(mat_ii[,ncols]);
	mat_ii[,1] <- as.factor(as.character(mat_ii[,1]));
	freq <- table(phenof,mat_ii[,1]);
	leveltmp <- length(levels(mat_ii[,1]));
	if ( (ncol(freq)==1) | ( (ncol(as.matrix(freq[,which(freq[1,]==0)])) >= (leveltmp-1)) | (ncol(as.matrix(freq[,which(freq[2,]==0)])) >= (leveltmp-1)) ) ) {check = TRUE};
	if (types[1] != "N") {
		for (j in 2:(ncol(mat_ii)-1)) {
			if (types[j-1] == "f") {
				mat_ii[,j] <- as.factor(as.character(mat_ii[,j]))
				freq <- table(phenof,mat_ii[,j])
				leveltmp <- length(levels(mat_ii[,j]))
				if ( (ncol(freq)==1) | ( (ncol(as.matrix(freq[,which(freq[1,]==0)])) >= (leveltmp-1)) | (ncol(as.matrix(freq[,which(freq[2,]==0)])) >= (leveltmp-1)) ) ) {check = TRUE} 
			} else if (types[j-1] == "n") {
				mat_ii[,j] <- as.numeric(mat_ii[,j])
			} else {
				stop("Please provide right data types for the covariates!")
			}
		}
}

	colnames(mat_ii) <- cols;
	if ( check == TRUE ) {
		model_ii <- glm(formula = Phenotype ~ .,family=binomial(link='logit'),data=mat_ii,control=glm.control(maxit=1))
	} else {
		model_ii <- glm(formula = Phenotype ~ .,family=binomial(link='logit'),data=mat_ii)
	}

	if (ptype == "CDS") {
		out_i <- c(dat[i,1])
	} else if (ptype == "SNP") {
		out_i <- c(dat[i,1], dat[i,2])
	} else if (ptype == "nonCDS") {
		out_i <- c(dat[i,1], dat[i,2], as.numeric(as.character(dat[i,3])))
	}

	summary_ii <- summary(model_ii);
	intercept_pVal <- summary_ii$coefficients[1,4];
	intercept_pBH=p.adjust(intercept_pVal, method = "BH")
	Null_deviance <- summary_ii$null.deviance;
	Residual_deviance <- summary_ii$deviance;
	out_i <- c(out_i, as.numeric(as.character(intercept_pVal)),as.numeric(as.character(intercept_pBH)), as.numeric(as.character(Null_deviance)), as.numeric(as.character(Residual_deviance)));
	coeffi <- as.data.frame(summary_ii$coefficients);
	coeffi$rownames <- rownames(coeffi);
	for (j in 1:length(keyname)) {
		key <- keyname[j]
		dattmp = coeffi[which(coeffi$rownames == key),]
		if (nrow(dattmp) == 0) {
			out_i <- c(out_i, NA, NA, NA)
		} else if (nrow(dattmp) == 1) {
			esiti_tmp <- as.numeric(as.character(dattmp[1,1]))
			zVal_tmp <- as.numeric(as.character(dattmp[1,3]))
			pVal_tmp <- as.numeric(as.character(dattmp[1,4]))
			if (ptype == "CDS") {
				out_i <- c(out_i, esiti_tmp, zVal_tmp, pVal_tmp)
			} else if (ptype == "SNP") {
				out_i <- c(out_i, esiti_tmp, pVal_tmp, zVal_tmp)
			} else if (ptype == "nonCDS") {
				out_i <- c(out_i, pVal_tmp, esiti_tmp, zVal_tmp)
			}
		} else {
			cat("There is wrong with:\t",key,"\n")
			stop()
		}
	}

	anova_ii <- anova(model_ii,test="LRT");
	anova_ii <- as.data.frame(anova_ii);
	anova_ii$rownames = rownames(anova_ii);
	for (j in 1:length(keyname2)) {
		key2 <- keyname2[j]
		dattmp2 = anova_ii[which(anova_ii$rownames == key2),]
		if (nrow(dattmp2) == 0) {
			out_i <- c(out_i, NA)
		} else if (nrow(dattmp) == 1) {
			anova_p <- as.numeric(as.character(dattmp2[1,5]))
			out_i <- c(out_i, anova_p)
		} else {
			cat("There is wrong with:\t",key,"\n")
			stop()
		}
	}
	
	R2 <- PseudoR2(model_ii);
	out_i <- c(out_i, as.numeric(as.character(R2[1])), as.numeric(as.character(R2[3])), as.numeric(as.character(R2[4])), as.numeric(as.character(R2[5])), as.numeric(as.character(R2[6])), as.numeric(as.character(R2[7])));
	return(out_i)
}


var <- ls(all.names=TRUE);
cl <- makeCluster(as.numeric(args[5]),type="FORK");
registerDoParallel(cl);
res <- foreach(i=1:rowNum,
              .combine=rbind,
              .export = var,
              .packages = c("BaylorEdPsych")) %dopar% {
       LogisRre(i)
}
stopImplicitCluster()


res_all <- data.frame(res,stringsAsFactors=FALSE);
colnames(res_all) <- Colnames;

if (ptype == "CDS") {
	for (i in 2:ncol(res_all)) {res_all[,i] <- as.numeric(res_all[,i])}
} else if (ptype == "SNP") {
	for (i in 3:ncol(res_all)) {res_all[,i] <- as.numeric(res_all[,i])}
} else {
	for (i in 4:ncol(res_all)) {res_all[,i] <- as.numeric(res_all[,i])}
}

res_all <- cbind( res_all[,1:(ncol(res_all)-6)][,colSums(is.na(res_all[,1:(ncol(res_all)-6)]))<length(res_all)] , res_all[,(ncol(res_all)-5):ncol(res_all)] );

if (ptype == "CDS") {
	for (i in 2:ncol(res_all)) {res_all[,i] <- signif(res_all[,i],3)}
} else if (ptype == "SNP") {
	for (i in 3:ncol(res_all)) {res_all[,i] <- signif(res_all[,i],3)}
} else {
	for (i in 4:ncol(res_all)) {res_all[,i] <- signif(res_all[,i],3)}
}

options(digits = 2,max.print=.Machine$integer.max,width=9999);
res_all <- res_all[order(res_all[,7]),];
res_filter <- res_all[-which(res_all[,7] > as.numeric(args[6])),];

write.table(res_all,args[7],sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE);
write.table(res_filter,args[8],sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE);

