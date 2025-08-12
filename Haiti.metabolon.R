
library(ggplot2)
library(ggrepel)
library(ape)
library(plyr)
library(reshape2)
library(cluster)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gplots)
library(vegan)
library(irr)
library(useful)
library(lubridate)
library(pscl)
library(parallel)
library(doParallel)
library(igraph)
library(randomForest)
library(ranger)
library(ROCR)
library(stringi)
library(mixOmics)
library(ggfortify)
library(ggsci)
library(ggforce)
library(ggbeeswarm)
library(ggvenn)
library(lme4)
library(lmerTest)
library(emmeans)
library(tableone)
library(abind)
library(limma)
library(psych)
library(glmnet)
library(caret)
library(gpboost)
library(SHAPforxgboost)
library(MetaboAnalystR)
library(umap)
library(WGCNA)


source("utils.R")
source("mcc.R")

cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")


get_color_list <- function(mvar) {
	missing_colors <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Dark2"), brewer.pal(8, "Pastel1"), brewer.pal(8, "Pastel2"), brewer.pal(8, "Accent"), brewer.pal(12, "Set3"))

	cols.cd <- c("#E71111", "#808080"); names(cols.cd) <- c("1", "0")
	cols.hivstatus <- c("#808080", "#E71111"); names(cols.hivstatus) <- c("Negative", "Positive")

	color_list <- list(cd=cols.cd, HIVStatus=cols.hivstatus)

	retval <- {}
	if (mvar %in% names(color_list)) {
		retval <- color_list[[mvar]]
	} else {
		retval <- missing_colors
	}
	return(retval)
}

siglevel <- 0.1
ntop <- 20
dircolors <- c("blue", "red", "grey", "black"); names(dircolors) <- c("down", "up", "NS", "manual")
shapes.sig <- c(19, 1); names(shapes.sig) <- c("sig", "NS")

sig_lookup <- function(p) {
	retval <- "ns"
	if (is.na(p)) {
		return("ns")
	}
	else if (p < 0.001) {
		retval <- "***"
	} else if (p < 0.01) {
		retval <- "**"
	} else if (p < 0.1) {
		retval <- "*"
	}
	return(retval)
}

# connect mapping file and Metabolon metadata file
metadata_fn <- "data/Haiti/SampleMetaData_UNMERGE.txt"
metadata <- read.table(metadata_fn, header=T, as.is=T, sep="\t")
metadata$SAMPLEID <- ifelse(metadata$SAMPLE_TYPE=="CMTRX", metadata$PARENT_SAMPLE_NAME, sprintf("%s.BMK", metadata$CLIENT_SAMPLE_ID))
rownames(metadata) <- metadata$SAMPLEID

mapping_fn <- "data/Haiti/Haiti_Mapping.metabolon.082724.txt"
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t", colClasses="character", comment.char="", quote="")
colnames(mapping)[1] <- "SampleID"
rownames(mapping) <- mapping$SampleID
sel <- intersect(rownames(mapping), metadata$SAMPLEID)
mapping <- mapping[sel,]
mapping$PARENT_SAMPLE_NAME <- metadata[rownames(mapping), "PARENT_SAMPLE_NAME"]


metadata_variables <- read.table("/Lab_Share/Haiti/metadata_variables.112823.txt", header=T, as.is=T, sep="\t", row.names=1)
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping[,mvar] <- factor(mapping[,mvar])
		if (metadata_variables[mvar, "baseline"] != "") {
			mapping[,mvar] <- relevel(mapping[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "ordered") {
		if (metadata_variables[mvar, "baseline"] == "") {
			lvls <- unique(mapping[,mvar])
		} else {
			lvls <- unlist(strsplit(metadata_variables[mvar, "baseline"], ","))
		}
		mapping[,mvar] <- ordered(mapping[,mvar], levels=lvls)
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping[,mvar] <- as.numeric(as.character(mapping[,mvar]))
	}
}

# set demographic variables for Table 1
demo_vars <- c("CD4Count", "CD4CountBinned", "ViralLoad", "ViralLoadBinned", "MomBMI", "MomAntibiotics", "Parity", "Delivery", "MomAge", "BreastProblems", "BBAge", "BabyGender", "ExclusiveBF", "ZWEI", "ZLEN")

# set some metabolites to always label
always_label <- c("kynurenine", "tryptophan", "9,10-DiHOME", "12,13-DiHOME", "X-12127", "X-12100", "dimethylarginine (SDMA + ADMA)", "cytosine")
always_label.extended <- c("kynurenine", "tryptophan", "serotonin", "kynurenate", "9,10-DiHOME", "12,13-DiHOME", "X-12127", "X-12100", "dimethylarginine (SDMA + ADMA)", "cytosine", "quinolinate", "indolelactate", "indolepropionate")

# enable weights for grouped analyses
weights.list <- {}
enable_weights <- FALSE

output_dir <- "output"
out_pdf <- sprintf("%s/metabolon_analysis.Haiti_%s.pdf", output_dir, format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)


#########################################################################################################
### read in Metabolon data
strip_quotes <- function(cvec) {
	return(gsub("\"", "", cvec))
}

metabolite_levels <- c("BIOCHEMICAL")
## ChemicalAnnotation.txt
chem_annot <- read.table("data/ZEBS/ChemicalAnnotation_UNMERGE.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
rownames(chem_annot) <- chem_annot$CHEM_ID
chem_annot$CHEMICAL_NAME <- strip_quotes(chem_annot$CHEMICAL_NAME)
chem_annot$SUB_PATHWAY <- strip_quotes(chem_annot$SUB_PATHWAY)
metabolon_map <- chem_annot[, c("CHEM_ID", "CHEMICAL_NAME", "SUB_PATHWAY", "SUPER_PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")]
colnames(metabolon_map) <- c("CHEM.ID", "CHEMICAL.NAME", "SUB.PATHWAY", "SUPER.PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")
rownames(metabolon_map) <- metabolon_map$CHEMICAL.NAME
## LogTransformedData.txt
df.metabolon <- list()
metabolon <- read.table(sprintf("data/ZEBS/LogTransformedData_UNMERGE.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping)[match(rownames(metabolon), mapping$PARENT_SAMPLE_NAME)]
# remove Drug - Antibiotics and Drug - Antiviral
to_remove <- rownames(metabolon_map)[which(metabolon_map[colnames(metabolon), "SUB.PATHWAY"] %in% c("Drug - Antibiotic", "Drug - Antiviral"))]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon[["BIOCHEMICAL"]] <- metabolon
# BatchNormalizedData.txt
metabolon <- read.table(sprintf("data/ZEBS/BatchNormalizedData_UNMERGE.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping)[match(rownames(metabolon), mapping$PARENT_SAMPLE_NAME)]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon[["BIOCHEMICAL_raw"]] <- metabolon
# PeakAreaData.txt
metabolon <- read.table(sprintf("data/ZEBS/PeakAreaData_UNMERGE.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping)[match(rownames(metabolon), mapping$PARENT_SAMPLE_NAME)]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon[["BIOCHEMICAL_peak"]] <- metabolon
cols.superpathway <- c(brewer.pal(9, "Set1"), "#bbbbbb", "#111111")[1:(length(unique(metabolon_map$SUPER.PATHWAY))+1)]; names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")


#########################################################################################################
## HUU vs HEU (by HIVStatus)

mvar <- "HIVStatus"
mapping.sel <- subset(mapping, HIVStatus %in% c("Negative", "Positive"))
mapping.sel[, mvar] <- droplevels(mapping.sel[,mvar])


## PCA
data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
pca <- prcomp(data.sel, center=F, scale=F)
eigs <- pca$sdev^2
pvar <- 100*(eigs / sum(eigs))
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID=rownames(pca$x))
for (mv in c("SampleID", "HIVStatus", demo_vars)) {
	df[, mv] <- mapping.sel[rownames(df), mv]
	if (metadata_variables[mv, "type"] == "factor") {
		p <- ggplot(df, aes_string(x="PC1", y="PC2", group=mv, color=mv)) + geom_point() + theme_classic() + ggtitle(sprintf("PCA (%s)", mv)) + scale_color_manual(values=get_color_list(mv), drop=T)
		print(p)
	} else {
		p <- ggplot(df, aes_string(x="PC1", y="PC2", group=mv, color=mv)) + geom_point() + theme_classic() + ggtitle(sprintf("PCA (%s)", mv))
		print(p)
	}
}


## boxplots of each metabolite over time
agg.melt <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt) <- c("SampleID", "metabolite", "value")
data.raw <- df.metabolon[["BIOCHEMICAL_raw"]][rownames(mapping.sel),]
agg.raw <- melt(as.matrix(data.raw), as.is=T); colnames(agg.raw) <- c("SampleID", "metabolite", "value")
agg.raw$Detected <- factor(!is.na(agg.raw$value))
data.peak <- df.metabolon[["BIOCHEMICAL_peak"]][rownames(mapping.sel),]
agg.peak <- melt(as.matrix(data.peak), as.is=T); colnames(agg.peak) <- c("SampleID", "metabolite", "value")
agg.peak$Detected <- factor(!is.na(agg.peak$value))
for (m in c("SampleID", "HIVStatus", "BBAge", "ViralLoadBinned")) {
	agg.melt[,m] <- mapping.sel[agg.melt$SampleID, m]
	agg.raw[,m] <- mapping.sel[agg.raw$SampleID, m]
	agg.peak[,m] <- mapping.sel[agg.peak$SampleID, m]
}

pl <- list()
for (m in sort(unique(agg.melt$metabolite))) {
	m.original <- m
	df <- subset(agg.melt, metabolite==m)
	df.raw <- subset(agg.raw, metabolite==m.original); rownames(df.raw) <- df.raw$SampleID
	df$Detected <- df.raw[df$SampleID, "Detected"]
	# compute detection statistics
	tab <- table(df[,c("HIVStatus", "Detected")])
	df.detection <- data.frame(x=rownames(tab), str=sprintf("%d/%d (%.1f%%)", tab[, "TRUE"], rowSums(tab), 100*(tab[, "TRUE"]/rowSums(tab)))); colnames(df.detection) <- c(mvar, "str")
	ypos <- min(df$value)*1.03
	p <- ggplot(df, aes_string(x=mvar, y="value", group=mvar, color=mvar)) + geom_boxplot(outlier.shape=NA) + geom_beeswarm(cex=2) + geom_text(data=df.detection, aes_string(x=mvar, label="str"), y=-Inf, size=2, inherit.aes=F) + theme_classic() + ggtitle(sprintf("%s (%d/%d - %.1f%% detected)", m.original, sum(tab[,"TRUE"]), sum(tab), 100*sum(tab[,"TRUE"])/sum(tab))) + scale_color_manual(values=get_color_list(mvar)) + theme(legend.position="none", plot.title=element_text(size=6), axis.text=element_text(size=6), axis.title=element_text(size=8))
	pl[[length(pl)+1]] <- p
	p <- ggplot(df, aes_string(x=mvar, y="value", group=mvar, color="ViralLoadBinned")) + geom_boxplot(outlier.shape=NA) + geom_beeswarm(cex=2) + geom_text(data=df.detection, aes_string(x=mvar, label="str"), y=-Inf, size=2, inherit.aes=F) + theme_classic() + ggtitle(sprintf("%s (%d/%d - %.1f%% detected)", m.original, sum(tab[,"TRUE"]), sum(tab), 100*sum(tab[,"TRUE"])/sum(tab))) + scale_color_manual(values=get_color_list("ViralLoadBinned")) + theme(plot.title=element_text(size=6), axis.text=element_text(size=6), axis.title=element_text(size=8))
	pl[[length(pl)+1]] <- p
	p <- ggplot(df, aes_string(x="BBAge", y="value", group="SampleID")) + stat_smooth(data=df, aes_string(x="BBAge", y="value"), inherit.aes=F, method="loess", se=T) + theme_classic() + ggtitle(sprintf("%s by %s", m.original, "BBAge")) + theme(plot.title=element_text(size=8), axis.text=element_text(size=6), axis.title=element_text(size=6), legend.position="none") + scale_color_manual(values=get_color_list(mvar)) + scale_fill_manual(values=get_color_list(mvar))
	pl[[length(pl)+1]] <- p
	p <- ggplot(df, aes_string(x="BBAge", y="value", group="SampleID", color=mvar)) + stat_smooth(data=df, aes_string(x="BBAge", y="value", group=mvar, color=mvar, fill=mvar), inherit.aes=F, method="loess", se=T) + theme_classic() + ggtitle(sprintf("%s by %s", m.original, mvar)) + theme(plot.title=element_text(size=8), axis.text=element_text(size=6), axis.title=element_text(size=6), legend.position="none") + scale_color_manual(values=get_color_list(mvar)) + scale_fill_manual(values=get_color_list(mvar))
	pl[[length(pl)+1]] <- p
}
multiplot(plotlist=pl, cols=2, rows=2)


## linear regression
mlevel <- "BIOCHEMICAL"; data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds==0)); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites with zero variance
df <- merge(data.sel, mapping.sel, by="row.names"); 
res <- {}
for (mvar in c("HIVStatus")) {
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
}
res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
res$padj <- p.adjust(res$p.value, method="fdr")
res <- res[order(res$p.value, decreasing=F),]
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
res$exp_estimate <- exp(res$estimate)
for (addvar in c("COMP_ID", "HMDB", "KEGG", "PUBCHEM", "PLATFORM", "SUB.PATHWAY", "SUPER.PATHWAY")) {
	res[, addvar] <- metabolon_map[as.character(res$metabolite), addvar]
}
write.table(res, file=sprintf("%s/emmeans.%s.%s.txt", output_dir, "Haiti", mvar), quote=F, sep="\t", row.names=F, col.names=T)
# volcano plots and circular dendrogram (Circos-style) plots
res <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, "Haiti", mvar), header=T, as.is=T, sep="\t", quote="", comment.char="")
res$contrast <- factor(res$contrast)
for (mv in unique(res$metadata_variable)) {
	df <- subset(res, metadata_variable==mv); df$contrast <- droplevels(df$contrast)
	for (ct in levels(df$contrast)) {
		df2 <- subset(df, contrast==ct & !is.na(estimate))
		lims <- max(abs(df2$estimate), na.rm=T)
		pl <- list()
		
		df2$siglabel <- ifelse(df2$dir=="NS", ifelse(df2$metabolite %in% always_label, df2$metabolite, NA), df2$metabolite)
		df2$dir <- ifelse(df2$dir == "NS", ifelse(df2$metabolite %in% always_label, "manual", "NS"), df2$dir)
		df2$padj <- pmax(df2$padj, 0.001) # censor at padj < 0.001 to improve plotting
		p <- ggplot(df2, aes(x=estimate, y=-log10(padj), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=3, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s %s", mv, ct)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
		p <- ggplot(df2, aes(x=estimate, y=-log10(p.value), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=3, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s %s (unadjusted p)", mv, ct)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
	}
}


## linear regression (with BBAge as a covariate)
mlevel <- "BIOCHEMICAL"; data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds==0)); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites with zero variance
df <- merge(data.sel, mapping.sel, by="row.names"); 
res <- {}
for (mvar in c("HIVStatus")) {
	enable_weights <- mvar %in% weights.list
	for (metabolite in colnames(data.sel)) {
		mod <- lm(as.formula(sprintf("%s ~ %s + BBAge", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr
		res <- rbind(res, tmp)
	}
}
res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
res$padj <- p.adjust(res$p.value, method="fdr")
res <- res[order(res$p.value, decreasing=F),]
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
res$exp_estimate <- exp(res$estimate)
for (addvar in c("COMP_ID", "HMDB", "KEGG", "PUBCHEM", "PLATFORM", "SUB.PATHWAY", "SUPER.PATHWAY")) {
	res[, addvar] <- metabolon_map[as.character(res$metabolite), addvar]
}
write.table(res, file=sprintf("%s/emmeans.%s.%s.with_BBAge.txt", output_dir, "Haiti", mvar), quote=F, sep="\t", row.names=F, col.names=T)
# volcano plots and circular dendrogram (Circos-style) plots
res <- read.table(sprintf("%s/emmeans.%s.%s.with_BBAge.txt", output_dir, "Haiti", mvar), header=T, as.is=T, sep="\t", quote="", comment.char="")
res$contrast <- factor(res$contrast)
for (mv in unique(res$metadata_variable)) {
	df <- subset(res, metadata_variable==mv); df$contrast <- droplevels(df$contrast)
	for (ct in levels(df$contrast)) {
		df2 <- subset(df, contrast==ct & !is.na(estimate))
		lims <- max(abs(df2$estimate), na.rm=T)
		pl <- list()
		
		df2$siglabel <- ifelse(df2$dir=="NS", ifelse(df2$metabolite %in% always_label, df2$metabolite, NA), df2$metabolite)
		df2$dir <- ifelse(df2$dir == "NS", ifelse(df2$metabolite %in% always_label, "manual", "NS"), df2$dir)
		df2$padj <- pmax(df2$padj, 0.001) # censor at padj < 0.001 to improve plotting
		p <- ggplot(df2, aes(x=estimate, y=-log10(padj), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=3, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s %s (+BBAge)", mv, ct)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
		p <- ggplot(df2, aes(x=estimate, y=-log10(p.value), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=3, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s %s (+BBAge, unadjusted p)", mv, ct)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
	}
}


## comparison of linear regression coefficients from Haiti versus ZEBS
visits.selected <- c("1 Wk", "1 Mo", "4 Mo", "4.5 Mo", "6 Mo", "9 Mo", "12 Mo", "15 Mo", "18 Mo")
res.zebs <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, "All_vs_HUU", "all_visits"), header=T, as.is=T, sep="\t", quote="", comment.char="")
res.zebs$contrast <- factor(res.zebs$contrast)
res <- read.table(sprintf("%s/emmeans.%s.%s.with_BBAge.txt", output_dir, "Haiti", mvar), header=T, as.is=T, sep="\t", quote="", comment.char="")
res$contrast <- factor(res$contrast)

sel.metabolites <- intersect(res.zebs$metabolite, res$metabolite)
for (visit in visits.selected) {
	df <- subset(res.zebs, Visit==visit & metabolite %in% sel.metabolites)
	rownames(df) <- df$metabolite
	df <- df[sel.metabolites,]
	df$estimate.haiti <- res[match(sel.metabolites, res$metabolite), "estimate"]
	df$estimate.haiti <- pmin(df$estimate.haiti, 1.0) # censor Haiti estimates at 1.0 to improve plotting
	df$siglabel <- ifelse(df$dir=="NS", ifelse(df$metabolite %in% always_label, df$metabolite, NA), df$metabolite)
	lims <- max(abs(c(df$estimate, df$estimate.haiti)), na.rm=T)
	
	test <- cor.test(~ estimate.haiti + estimate, data=df, method="pearson")
	test2 <- cor.test(~ estimate.haiti + estimate, data=df, method="spearman")
	p <- ggplot(df, aes(x=estimate.haiti, y=estimate)) + geom_point(aes(color=dir)) + geom_text_repel(aes(color=dir, label=siglabel), size=3, max.overlaps=15) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + stat_smooth(method="lm", se=F) + theme_classic() + ggtitle(sprintf("Correlation in estimates (Haiti vs ZEBS, %s, r=%.4g, p=%.4g, rho=%.4g, p=%.4g)", visit, test$estimate, test$p.value, test2$estimate, test2$p.value)) + xlim(c(-lims, lims)) + ylim(c(-lims, lims)) + scale_color_manual(values=dircolors)
	print(p)
	
	df2 <- subset(df, dir!="NS")
	test <- cor.test(~ estimate.haiti + estimate, data=df2, method="pearson")
	test2 <- cor.test(~ estimate.haiti + estimate, data=df2, method="spearman")
	p <- ggplot(df2, aes(x=estimate.haiti, y=estimate)) + geom_point(aes(color=dir)) + geom_text_repel(aes(color=dir, label=siglabel), size=3, max.overlaps=15) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + stat_smooth(method="lm", se=F) + theme_classic() + ggtitle(sprintf("Correlation in estimates sighits only (Haiti vs ZEBS, %s, r=%.4g, p=%.4g, rho=%.4g, p=%.4g)", visit, test$estimate, test$p.value, test2$estimate, test2$p.value)) + xlim(c(-lims, lims)) + ylim(c(-lims, lims)) + scale_color_manual(values=dircolors)
	print(p)
}


## kynurenine to tryptophan ratio
data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
df <- merge(data.sel, mapping.sel, by="row.names"); 
df$KTratio <- exp(df[, "kynurenine"]) / exp(df[, "tryptophan"])
df$KTlogratio <- log10(df$KTratio)
# by emmeans
res <- {}
for (mvar in c("HIVStatus")) {
	mod <- lm(as.formula(sprintf("%s ~ %s", "KTlogratio", mvar)), data=df); modelstr <- "LM"
	emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s", mvar)), adjust="none")
	tmp <- as.data.frame(emm$contrasts)
	tmp$metabolite <- "KTlogratio"; tmp$metadata_variable <- mvar; tmp$model <- modelstr
	tmp <- subset(tmp, !is.na(estimate))
	tmp <- tmp[,c("metabolite", setdiff(colnames(tmp), "metabolite"))]
	tmp$padj <- p.adjust(tmp$p.value, method="fdr")
	tmp <- tmp[order(tmp$p.value, decreasing=F),]
	tmp$dir <- ifelse(tmp$padj < siglevel, ifelse(sign(tmp$estimate)==1, "up", "down"), "NS")
	res <- rbind(res, tmp)
	# boxplot
	p <- ggplot(df, aes_string(x=mvar, y="KTlogratio", color=mvar)) + geom_boxplot(outlier.shape=NA) + geom_beeswarm(size=3, cex=1.5) + theme_classic() + ggtitle(sprintf("%s by %s", "KT log-ratio", mvar)) + scale_color_manual(values=get_color_list(mvar))
	print(p)
}
write.table(res, file=sprintf("%s/KTRatio.%s.%s.txt", output_dir, "Haiti", mvar), quote=F, sep="\t", row.names=F, col.names=T)



dev.off()


