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
library(stringr)
library(mixOmics)
library(ggfortify)
library(ggsci)
library(ggforce)
library(ggbeeswarm)
library(ggvenn)
library(ggraph)
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
	cols.currentbf <- pal_startrek("uniform")(7)[1:2]; names(cols.currentbf) <- c("No", "Yes")
	cols.hivexposure <- c("#808080", "#E71111"); names(cols.hivexposure) <- c("Unexposed", "Exposed")
	
	color_list <- list(cd=cols.cd, CurrentBF=cols.currentbf, HIVExposure=cols.hivexposure)
	
	retval <- {}
	if (mvar %in% names(color_list)) {
		retval <- color_list[[mvar]]
	} else {
		retval <- missing_colors
	}
	return(retval)
}

siglevel <- 0.05
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
	} else if (p < 0.05) {
		retval <- "*"
	}
	return(retval)
}

# connect mapping file and Metabolon metadata file
metadata_fn <- "data/ZEBS/SampleMetaData.txt"
metadata <- read.table(metadata_fn, header=T, as.is=T, sep="\t")
metadata$datestr <- unlist(lapply(metadata$CLIENT_IDENTIFIER, function(x) unlist(strsplit(x, " "))[2]))
metadata$Collection.Date <- as.Date(parse_date_time(metadata$datestr, orders="dmy"))
metadata$SAMPLEID <- sprintf("%s.%s", metadata$SUBJECT_OR_ANIMAL_ID, format(metadata$Collection.Date, "%m%d%Y"))
rownames(metadata) <- metadata$SAMPLEID

mapping_fn <- "data/ZEBS/ZEBS_BMK_1410_Mapping.082624.txt"
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t", colClasses="character")
rownames(mapping) <- mapping$SampleID
sel <- intersect(rownames(mapping), metadata$SAMPLEID)
mapping <- mapping[sel,] # filter to samples with valid Metabolon data (n=1597)
mapping$PARENT_SAMPLE_NAME <- metadata[rownames(mapping), "PARENT_SAMPLE_NAME"]

mapping.full <- mapping

# filter to selected Visits
visits.selected <- c("1 Wk", "1 Mo", "4 Mo", "4.5 Mo", "6 Mo", "9 Mo", "12 Mo", "15 Mo", "18 Mo")
visit_to_age_in_months <- c(0.25, 1, 4, 4.5, 6, 9, 12, 15, 18); names(visit_to_age_in_months) <- visits.selected
mapping <- subset(mapping, Visit %in% visits.selected)

metadata_variables <- read.table("data/ZEBS/metadata_variables.082624.txt", header=T, as.is=T, sep="\t", row.names=1)
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
# set specific levels for Visit
visits <- c("0 Wk", "1 Wk", "1 Mo", "2 Mo", "3 Mo", "4 Mo", "4.5 Mo", "5 Mo", "6 Mo", "9 Mo", "12 Mo", "15 Mo", "18 Mo", "21 Mo", "24 Mo")
visits.filt <- c("1 Wk", "1 Mo", "4 Mo", "4.5 Mo", "6 Mo", "9 Mo", "12 Mo", "15 Mo", "18 Mo")
mapping$Visit <- factor(as.character(mapping$Visit), levels=visits.filt)
mapping$lnvis <- factor(as.character(mapping$lnvis), levels=visits)
mapping$fpvis <- factor(as.character(mapping$fpvis), levels=visits)

# NA out log10vlrna for HUU
mapping[which(mapping$Group3=="HUU"), "log10vlrna"] <- NA

# remove all CurrentBF=No samples
mapping <- subset(mapping, CurrentBF=="Yes")

# set demographic variables for Table 1
demo_vars <- c("HIVExposure", "MOMAGE", "log10vlrna", "delgage", "MODE", "CD4C", "CD8C", "CD3C", "BWT", "SEX", "CD4_12M", "CD8_12M", "CD3_12M", "cd", "md", "timec", "timem")

# set some metabolites to always label
always_label <- c("kynurenine", "tryptophan", "9,10-DiHOME", "12,13-DiHOME", "X-12127", "X-12100", "dimethylarginine (SDMA + ADMA)", "cytosine")
always_label.extended <- c("kynurenine", "tryptophan", "serotonin", "kynurenate", "9,10-DiHOME", "12,13-DiHOME", "X-12127", "X-12100", "dimethylarginine (SDMA + ADMA)", "cytosine", "quinolinate", "indolelactate", "indolepropionate")
always_label.1mo <- c("kynurenine", "tryptophan", "9,10-DiHOME", "12,13-DiHOME", "X-12127", "X-12100", "dimethylarginine (SDMA + ADMA)", "cytosine", "N-acetylalanine", "N-acetylserine", "X-07765")

weights.list <- {}
enable_weights <- FALSE

# list of emmeans comparisons for boxplot labeling
emmeans_files <- list("emmeans.All_vs_HUU.all_visits.txt"=c("Exposed - Unexposed"), "emmeans.HEU_vs_HUU.all_visits.txt"=c("HEU - HUU"))

# list of correlation variables
correlation_vars <- c("CD4C", "log10vlrna", "bm_vl_log10")

output_dir <- "output/"
out_pdf <- sprintf("%s/metabolon_analysis.%s.pdf", output_dir, format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)


#########################################################################################################
### read in Metabolon data
metabolite_levels <- c("BIOCHEMICAL")
## ChemicalAnnotation.txt
chem_annot <- read.table("data/ZEBS/ChemicalAnnotation.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
rownames(chem_annot) <- chem_annot$CHEM_ID
metabolon_map <- chem_annot[, c("CHEM_ID", "CHEMICAL_NAME", "SUB_PATHWAY", "SUPER_PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")]
colnames(metabolon_map) <- c("CHEM.ID", "CHEMICAL.NAME", "SUB.PATHWAY", "SUPER.PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")
rownames(metabolon_map) <- metabolon_map$CHEMICAL.NAME
## LogTransformedData.txt
df.metabolon <- list()
metabolon <- read.table(sprintf("data/ZEBS/LogTransformedDataWithUnnamed.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping)[match(rownames(metabolon), mapping$PARENT_SAMPLE_NAME)]
# remove Drug - Antibiotics and Drug - Antiviral
to_remove <- rownames(metabolon_map)[which(metabolon_map[colnames(metabolon), "SUB.PATHWAY"] %in% c("Drug - Antibiotic", "Drug - Antiviral"))]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon[["BIOCHEMICAL"]] <- metabolon
# BatchNormalizedData.txt
metabolon <- read.table(sprintf("data/ZEBS/BatchNormalizedData.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping)[match(rownames(metabolon), mapping$PARENT_SAMPLE_NAME)]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon[["BIOCHEMICAL_raw"]] <- metabolon
# PeakAreaData.txt
metabolon <- read.table(sprintf("data/ZEBS/PeakAreaData.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping)[match(rownames(metabolon), mapping$PARENT_SAMPLE_NAME)]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon[["BIOCHEMICAL_peak"]] <- metabolon
cols.superpathway <- c(brewer.pal(9, "Set1"), "#bbbbbb", "#111111")[1:(length(unique(metabolon_map$SUPER.PATHWAY))+1)]; names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")


#########################################################################################################
### read in Metabolon data (n=1597 full mapping set)
metabolite_levels <- c("BIOCHEMICAL")
## ChemicalAnnotation.txt
chem_annot <- read.table("data/ZEBS/ChemicalAnnotation.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
rownames(chem_annot) <- chem_annot$CHEM_ID
metabolon_map <- chem_annot[, c("CHEM_ID", "CHEMICAL_NAME", "SUB_PATHWAY", "SUPER_PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")]
colnames(metabolon_map) <- c("CHEM.ID", "CHEMICAL.NAME", "SUB.PATHWAY", "SUPER.PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")
rownames(metabolon_map) <- metabolon_map$CHEMICAL.NAME
## LogTransformedData.txt
df.metabolon.full <- list()
metabolon <- read.table(sprintf("data/ZEBS/LogTransformedDataWithUnnamed.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping.full$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping.full)[match(rownames(metabolon), mapping.full$PARENT_SAMPLE_NAME)]
# remove Drug - Antibiotics and Drug - Antiviral
to_remove <- rownames(metabolon_map)[which(metabolon_map[colnames(metabolon), "SUB.PATHWAY"] %in% c("Drug - Antibiotic", "Drug - Antiviral"))]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon.full[["BIOCHEMICAL"]] <- metabolon
# BatchNormalizedData.txt
metabolon <- read.table(sprintf("data/ZEBS/BatchNormalizedData.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping.full$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping.full)[match(rownames(metabolon), mapping.full$PARENT_SAMPLE_NAME)]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon.full[["BIOCHEMICAL_raw"]] <- metabolon
# PeakAreaData.txt
metabolon <- read.table(sprintf("data/ZEBS/PeakAreaData.txt"), header=T, as.is=T, sep="\t", quote="", comment.char="", row.names=1)
colnames(metabolon) <- chem_annot[gsub("^X", "", colnames(metabolon)), "CHEMICAL_NAME"]
sel <- which(rownames(metabolon) %in% mapping.full$PARENT_SAMPLE_NAME); metabolon <- metabolon[sel,]
rownames(metabolon) <- rownames(mapping.full)[match(rownames(metabolon), mapping.full$PARENT_SAMPLE_NAME)]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]
df.metabolon.full[["BIOCHEMICAL_peak"]] <- metabolon



##########################################################################################################
### analysis by HIVExposure

mapping.sel <- mapping
mvar <- "HIVExposure"

dm <- list()
distance_metrics <- c("euclidean")

for (distance_metric in distance_metrics) {
	dm[[length(dm)+1]] <- as.matrix(vegdist(df.metabolon[["BIOCHEMICAL"]], method=distance_metric))
}
names(dm) <- distance_metrics

## sample diagram
for (mvar in c("HIVExposure")) {
	agg <- aggregate(Visit ~ ID1, mapping.sel, function(x) max(as.numeric(x))); agg <- agg[order(agg$Visit),]
	mapping.sel$ID1 <- factor(mapping.sel$ID1, levels=agg$ID1)
	tmp <- table(unique(mapping.sel[,c("ID1", mvar)])[mvar])
	lut <- sprintf("%s (n=%d)", names(tmp),tmp); names(lut) <- names(tmp)
	df <- mapping.sel
	df$mvarstr <- lut[df[,mvar]]
	p <- ggplot(df, aes_string(x="Visit", y="ID1", group="ID1", color=mvar)) + geom_point() + geom_line() + facet_wrap(as.formula(sprintf("~%s", mvar)), ncol=1, scales="free_y") + theme_classic() + ggtitle(sprintf("Overall sample diagram (n=%d samples from %d PIDs)", nrow(df), length(unique(df$ID1)))) + scale_color_manual(values=get_color_list(mvar), drop=T)
	print(p)
	p <- ggplot(df, aes_string(x="Visit", y="ID1", group="ID1", color=mvar)) + geom_point() + geom_line() + facet_wrap(as.formula(sprintf("~%s", "mvarstr")), ncol=1, scales="free_y") + theme_classic() + ggtitle(sprintf("Overall sample diagram (n=%d samples from %d PIDs)", nrow(df), length(unique(df$ID1)))) + scale_color_manual(values=get_color_list(mvar), drop=T)
	print(p)
	tab <- table(mapping.sel[,c(mvar, "Visit")])
	df <- melt(tab)
	p <- ggplot(df, aes_string(x="Visit", y=mvar, size="value", color=mvar)) + geom_point(aes(size=2*value), shape=1) + geom_text(aes(label=value), hjust=0.5, size=2) + scale_size_area() + theme_classic() + ggtitle(sprintf("Overall sample diagram (n=%d samples from %d PIDs)", nrow(mapping.sel), length(unique(mapping.sel$ID1)))) + scale_color_manual(values=get_color_list(mvar), drop=T) + scale_size_continuous(range = c(1,16))
	print(p)
	tab <- table(mapping.sel[,mvar])
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Sample counts by %s", mvar)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(melt(tab)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
	print(p)
}


## tableOne
demo <- unique(mapping.sel[, c("ID1", demo_vars)])
demo$CD4_gt_350 <- ifelse(demo$CD4C > 350, ">350", "<=350")
demo$CD4_gt_300 <- ifelse(demo$CD4C > 300, ">300", "<=300")
demo$CD4_gt_200 <- ifelse(demo$CD4C > 200, ">200", "<=200")
for (strata_var in c("HIVExposure")) {
	tab1 <- CreateTableOne(vars=c(demo_vars, "CD4_gt_350", "CD4_gt_300", "CD4_gt_200"), strata=c(strata_var), data=demo, smd=T)
	write.table(print(tab1, noSpaces=T), file=sprintf("%s/Table_1.%s.txt", output_dir, strata_var), quote=F, sep="\t", row.names=T, col.names=T)
	write.table(print(tab1, noSpaces=T, nonnormal=demo_vars), file=sprintf("%s/Table_1_nonnormal.%s.txt", output_dir, strata_var), quote=F, sep="\t", row.names=T, col.names=T)
}


## PCA
pca <- prcomp(df.metabolon[["BIOCHEMICAL"]], center=F, scale=F)
eigs <- pca$sdev^2
pvar <- 100*(eigs / sum(eigs))
df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], SampleID=rownames(pca$x))
mvar <- "HIVExposure"
p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar, label="SampleID")) + geom_point() + geom_text_repel(size=4) + theme_classic() + ggtitle(sprintf("%s PCA (Euclidean distance)", mvar)) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + stat_ellipse(mapping=aes_string(x="PC1", y="PC2", colour=mvar), inherit.aes=F, type="t") + scale_color_manual(values=get_color_list(mvar))
print(p)
# plot loadings
for (pc in 1:3) {
	d <- pca$rotation[,pc]
	d <- data.frame(biomarker=names(d), loading=d)
	d <- d[order(abs(d$loading), decreasing=T),]
	d <- d[1:20,]
	d$biomarker <- factor(as.character(d$biomarker), levels=as.character(d$biomarker))
	p <- ggplot(d, aes(x=biomarker, y=loading)) + geom_bar(stat="identity") + scale_fill_manual(values="blue") + coord_flip() + theme_classic() + ggtitle(sprintf("PC %d", pc))
	print(p)
}


## PERMANOVA on metabolomics
permanova_vars <- c("Visit", "HIVExposure", "CD4C", "SEX")
sink(sprintf("%s/PERMANOVA.%s.txt", output_dir, "metabolites"), append=F)
for (distance_metric in distance_metrics) {
	print(distance_metric)

	inds_to_remove <- which(is.na(mapping[,permanova_vars]), arr.ind=T)[, "row"]
	mapping.permanova <- mapping[setdiff(1:nrow(mapping), inds_to_remove),]
	dm.permanova <- dm[[distance_metric]][rownames(mapping.permanova), rownames(mapping.permanova)]
	
	form <- as.formula(sprintf("as.dist(dm.permanova) ~ %s", paste(permanova_vars, collapse="+")))
	res <- adonis2(form , data=mapping.permanova, permutations=999, by="margin")
	res$R2 <- res$SumOfSqs / sum(res$SumOfSqs)
	print(res)
}
sink()


## boxplots of each metabolite over time
mvars.plot <- c("HIVExposure")
mlevel <- "BIOCHEMICAL"; data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
agg.melt <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt) <- c("SampleID", "metabolite", "value")
data.raw <- df.metabolon[["BIOCHEMICAL_raw"]][rownames(mapping.sel),]
agg.raw <- melt(as.matrix(data.raw), as.is=T); colnames(agg.raw) <- c("SampleID", "metabolite", "value")
agg.raw$Detected <- factor(!is.na(agg.raw$value))
data.peak <- df.metabolon[["BIOCHEMICAL_peak"]][rownames(mapping.sel),]
agg.peak <- melt(as.matrix(data.peak), as.is=T); colnames(agg.peak) <- c("SampleID", "metabolite", "value")
agg.peak$Detected <- factor(!is.na(agg.peak$value))
for (m in c("ID1", "Group2", "Group3", "HIVExposure", "Visit", "log10vlrna", "CD4C", "bm_vl_log10")) {
	agg.melt[,m] <- mapping.sel[agg.melt$SampleID, m]
	agg.raw[,m] <- mapping.sel[agg.raw$SampleID, m]
	agg.peak[,m] <- mapping.sel[agg.peak$SampleID, m]
}

tmp <- read.table(sprintf("%s/%s", output_dir, "emmeans.All_vs_HUU.all_visits.txt"), header=T, as.is=T, sep="\t", comment.char="", quote="")
rfres <- read.table(sprintf("%s/ranger_METABOLITE.%s.all_visits.features_detailed.txt", output_dir, "Exposed_vs_Unexposed"), header=T, as.is=T, sep="\t", comment.char="", quote="")
res.corr <- {}
for (m in sort(unique(agg.melt$metabolite))) {
	pl <- list()
	lmres <- {}
	for (emmeans_file in names(emmeans_files)) {
		selected_contrasts <- emmeans_files[[emmeans_file]]
		tmp <- subset(read.table(sprintf("%s/%s", output_dir, emmeans_file), header=T, as.is=T, sep="\t", comment.char="", quote=""), metabolite==m)
		tmp <- subset(tmp, contrast %in% selected_contrasts)
		tmp <- tmp[, setdiff(colnames(tmp), "weighted")]
		lmres <- rbind(lmres, tmp)
	}
	mvar <- "HIVExposure"
	pl <- list()
	for (mvar in mvars.plot) {
		df <- subset(agg.melt, metabolite==m)
		df.siglab <- subset(lmres, metabolite==m & metadata_variable==mvar)
		df.siglab2 <- aggregate(dir ~ Visit, df.siglab, function(x) ifelse(any(x!="NS"), "*", ""))
		df.siglab2$Visit <- factor(df.siglab2$Visit, levels=levels(df$Visit))
		df.rf <- subset(rfres, metabolite==m & metadata_variable==mvar)
		if (nrow(df.rf) > 0) {
			df.rf$dir <- "*" # only features selected by RF are in this data frame to begin with
			df.rf$Visit <- factor(df.rf$Visit, levels=levels(df$Visit)) 
		}
		p <- ggplot(df, aes_string(x="Visit", y="value", group="ID1", color=mvar)) + stat_smooth(data=df, aes_string(x="Visit", y="value", group=mvar, color=mvar), inherit.aes=F, method="loess", se=T, n=nlevels(df$Visit)) + theme_classic() + ggtitle(sprintf("%s by %s", m, mvar)) + theme(plot.title=element_text(size=8), axis.text=element_text(size=6), axis.title=element_text(size=6), legend.position="none") + scale_color_manual(values=get_color_list(mvar))
		pl[[length(pl)+1]] <- p
		p <- ggplot(df, aes_string(x="Visit", y="value", group="ID1", color=mvar)) + stat_smooth(data=df, aes_string(x="Visit", y="value", group=mvar, color=mvar), inherit.aes=F, method="loess", se=T, n=nlevels(df$Visit)) + geom_text(data=df.siglab2, aes_string(x="Visit", y=Inf, label="dir"), size=8, color="red", vjust="inward", inherit.aes=F) + theme_classic() + ggtitle(sprintf("%s by %s", m, mvar)) + theme(plot.title=element_text(size=8), axis.text=element_text(size=6), axis.title=element_text(size=6), legend.position="none") + scale_color_manual(values=get_color_list(mvar))
		pl[[length(pl)+1]] <- p
		p <- ggplot(df, aes_string(x="Visit", y="value", group="ID1", color=mvar)) + stat_smooth(data=df, aes_string(x="Visit", y="value", group=mvar, color=mvar), inherit.aes=F, method="loess", se=T, n=nlevels(df$Visit)) + theme_classic() + ggtitle(sprintf("%s by %s", m, mvar)) + theme(plot.title=element_text(size=8), axis.text=element_text(size=6), axis.title=element_text(size=6), legend.position="none") + scale_color_manual(values=get_color_list(mvar))
		df.siglab3 <- subset(df.siglab2, dir=="*")
		if (nrow(df.siglab3) > 0) {
			p <- p + geom_rect(data=df.siglab3, aes_string(xmin="as.numeric(Visit)-0.25", xmax="as.numeric(Visit)+0.2", ymin=-Inf, ymax=Inf), size=2, alpha=0.25, fill="red", color=NA, inherit.aes=F)
		}
		pl[[length(pl)+1]] <- p
		p <- ggplot(df, aes_string(x="Visit", y="value", group="ID1", color=mvar)) + stat_smooth(data=df, aes_string(x="Visit", y="value", group=mvar, color=mvar), inherit.aes=F, method="loess", se=T, n=nlevels(df$Visit)) + geom_text(data=df.siglab2, aes_string(x="Visit", y=Inf, label="dir"), size=8, color="red", vjust="inward", inherit.aes=F) + theme_classic() + ggtitle(sprintf("%s by %s", m, mvar)) + theme(plot.title=element_text(size=8), axis.text=element_text(size=6), axis.title=element_text(size=6), legend.position="none") + scale_color_manual(values=get_color_list(mvar))
		if (nrow(df.rf) > 0) {
			p <- p + geom_text(data=df.rf, aes_string(x="Visit", y=-Inf, label="dir"), size=8, color="black", vjust="inward", inherit.aes=F)
		}
		pl[[length(pl)+1]] <- p
		df.rf2 <- subset(df.rf, dir=="*")
		p <- ggplot(df, aes_string(x="Visit", y="value", group="ID1", color=mvar)) + stat_smooth(data=df, aes_string(x="Visit", y="value", group=mvar, color=mvar), inherit.aes=F, method="loess", se=T, n=nlevels(df$Visit)) + geom_text(data=df.siglab2, aes_string(x="Visit", y=Inf, label="dir"), size=8, color="red", vjust="inward", inherit.aes=F) + theme_classic() + ggtitle(sprintf("%s by %s", m, mvar)) + theme(plot.title=element_text(size=8), axis.text=element_text(size=6), axis.title=element_text(size=6), legend.position="none") + scale_color_manual(values=get_color_list(mvar))
		if (nrow(df.rf2) > 0) {
			p <- p + geom_rect(data=df.rf2, aes_string(xmin="as.numeric(Visit)-0.25", xmax="as.numeric(Visit)+0.2", ymin=-Inf, ymax=Inf), size=2, alpha=0.25, fill="red", color=NA, inherit.aes=F)
		}
		pl[[length(pl)+1]] <- p
		
		for (cvar in correlation_vars) {
			test <- cor.test(as.formula(sprintf("~ value + %s", cvar)), data=df, method="spearman")
			res.corr <- rbind(res.corr, c(m, cvar, test$statistic, test$p.value, test$estimate, test$method))
			p <- ggplot(df, aes_string(x=cvar, y="value")) + geom_point(size=0.8) + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("%s vs %s (rho=%.4g, p=%.4g)", cvar, m, test$estimate, test$p.value)) + theme(plot.title=element_text(size=6), text=element_text(size=6))
			pl[[length(pl)+1]] <- p			
			p <- ggplot(df, aes_string(x=cvar, y="value", group=mvar, color=mvar)) + geom_point(size=0.8) + stat_smooth(method="lm") + scale_color_manual(values=get_color_list(mvar)) + theme_classic() + ggtitle(sprintf("%s vs %s by %s", cvar, m, mvar)) + theme(plot.title=element_text(size=6), text=element_text(size=6))
			pl[[length(pl)+1]] <- p
		}
	}
	multiplot(plotlist=pl, cols=3, rows=3)
}
res.corr <- as.data.frame(res.corr)
colnames(res.corr) <- c("metabolite", "cvar", "statistic", "p.value", "estimate", "method")
res.corr$p.value <- as.numeric(res.corr$p.value)
res.corr$padj <- p.adjust(res.corr$p.value, method="fdr")
res.corr$estimate <- as.numeric(res.corr$estimate)
res.corr <- res.corr[order(res.corr$p.value, decreasing=F),]
res.corr$dir <- ifelse(res.corr$padj < siglevel, ifelse(sign(res.corr$estimate)==1, "up", "down"), "NS")
res.corr$exp_estimate <- exp(res.corr$estimate)
for (addvar in c("COMP_ID", "HMDB", "KEGG", "PUBCHEM", "PLATFORM", "SUB.PATHWAY", "SUPER.PATHWAY")) {
	res.corr[, addvar] <- metabolon_map[as.character(res.corr$metabolite), addvar]
}
write.table(res.corr, file=sprintf("%s/metabolite_correlation_CD4_log10vlrna.%s.txt", output_dir, "all_samples"), quote=F, sep="\t", row.names=F, col.names=T)


## linear regression (stratified by Visit)
mlevel <- "BIOCHEMICAL"; data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds==0)); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites with zero variance
df <- merge(data.sel, mapping.sel, by="row.names"); 
res <- {}
for (mvar in c("HIVExposure")) {
	for (metabolite in colnames(data.sel)) {
		mod <- lmer(as.formula(sprintf("%s ~ %s*Visit + (1 | ID1)", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | Visit", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr; tmp$weighted <- enable_weights
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
write.table(res, file=sprintf("%s/emmeans.%s.%s.txt", output_dir, "All_vs_HUU", "all_visits"), quote=F, sep="\t", row.names=F, col.names=T)
# volcano plots and rank-order plots
res <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, "All_vs_HUU", "all_visits"), header=T, as.is=T, sep="\t", quote="", comment.char="")
res$contrast <- factor(res$contrast); res$Visit <- droplevels(factor(res$Visit, levels=visits))
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS") # updated dir based on padj=0.05 cutoff
for (mv in unique(res$metadata_variable)) {
	df <- subset(res, metadata_variable==mv); df$contrast <- droplevels(df$contrast)
	for (ct in levels(df$contrast)) {
		df2 <- subset(df, contrast==ct & !is.na(estimate))
		df2$Visit <- droplevels(df2$Visit)
		lims <- max(abs(df2$estimate), na.rm=T)
		pl <- list()
		for (visit in levels(df2$Visit)) {
			df3 <- subset(df2, Visit==visit)
			df3$siglabel <- ifelse(df3$dir=="NS", ifelse(df3$metabolite %in% always_label, df3$metabolite, NA), df3$metabolite)
			df3$dir <- ifelse(df3$dir == "NS", ifelse(df3$metabolite %in% always_label, "manual", "NS"), df3$dir)
			p <- ggplot(df3, aes(x=estimate, y=-log10(padj), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=3, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s, %s", ct, visit)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
			pl[[length(pl)+1]] <- p
			df3 <- df3[order(df3$estimate),]
			df3$rank <- 1:nrow(df3)
			p <- ggplot(df3, aes(x=rank, y=estimate, color=dir)) + geom_point(size=2) + geom_text_repel(aes(label=siglabel), size=1.5, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s, %s", ct, visit)) + scale_color_manual(values=dircolors) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
			pl[[length(pl)+1]] <- p
		}
		multiplot(plotlist=pl, cols=2, rows=2)
	}
}


## linear regression (averaged over Visit)
mlevel <- "BIOCHEMICAL"; data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds==0)); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites with zero variance
df <- merge(data.sel, mapping.sel, by="row.names"); 
res <- {}
for (mvar in c("HIVExposure")) {
	for (metabolite in colnames(data.sel)) {
		mod <- lmer(as.formula(sprintf("%s ~ %s*Visit + (1 | ID1)", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts) 
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr; tmp$weighted <- enable_weights
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
write.table(res, file=sprintf("%s/emmeans.%s.%s.txt", output_dir, "All_vs_HUU", "averaged_over_Visit"), quote=F, sep="\t", row.names=F, col.names=T)
# volcano plots and rank-order plots
res <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, "All_vs_HUU", "averaged_over_Visit"), header=T, as.is=T, sep="\t", quote="", comment.char="")
res$contrast <- factor(res$contrast)
for (mv in unique(res$metadata_variable)) {
	df <- subset(res, metadata_variable==mv); df$contrast <- droplevels(df$contrast)
	for (ct in levels(df$contrast)) {
		df2 <- subset(df, contrast==ct & !is.na(estimate))
		lims <- max(abs(df2$estimate), na.rm=T)
		df3 <- df2
		df3$siglabel <- ifelse(df3$dir=="NS", ifelse(df3$metabolite %in% always_label, df3$metabolite, NA), df3$metabolite)
		df3$dir <- ifelse(df3$dir == "NS", ifelse(df3$metabolite %in% always_label, "manual", "NS"), df3$dir)
		df3$padj <- pmax(df3$padj, 0.001) # censor at padj < 0.001 to improve plotting
		p <- ggplot(df3, aes(x=estimate, y=-log10(padj), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=1.5, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s, averaged over Visit", ct)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
		df3 <- df3[order(df3$estimate),]
		df3$rank <- 1:nrow(df3)
		p <- ggplot(df3, aes(x=rank, y=estimate, color=dir)) + geom_point(size=2) + geom_text_repel(aes(label=siglabel), size=1.0, max.overlaps=35) + theme_classic() + ggtitle(sprintf("%s, averaged over Visit", ct)) + scale_color_manual(values=dircolors) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
	}
}


## kynurenine to tryptophan ratio
data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
df <- merge(data.sel, mapping.sel, by="row.names"); 
df$KTratio <- exp(df[, "kynurenine"]) / exp(df[, "tryptophan"])
df$KTlogratio <- log10(df$KTratio)
# by emmeans
res <- {}
for (mvar in c("HIVExposure")) {
	mod <- lmer(as.formula(sprintf("%s ~ %s*Visit + (1 | ID1)", "KTlogratio", mvar)), data=df); modelstr <- "LM"
	emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | Visit", mvar)), adjust="none")
	tmp <- as.data.frame(emm$contrasts)
	tmp$metabolite <- "KTlogratio"; tmp$metadata_variable <- mvar; tmp$model <- modelstr
	tmp <- subset(tmp, !is.na(estimate))
	tmp <- tmp[,c("metabolite", setdiff(colnames(tmp), "metabolite"))]
	tmp$padj <- p.adjust(tmp$p.value, method="fdr")
	tmp <- tmp[order(tmp$p.value, decreasing=F),]
	tmp$dir <- ifelse(tmp$padj < siglevel, ifelse(sign(tmp$estimate)==1, "up", "down"), "NS")
	res <- rbind(res, tmp)
	# boxplot
	pd <- position_dodge(width=0.8)
	p <- ggplot(df, aes_string(x="Visit", y="KTlogratio", fill=mvar)) + geom_boxplot(position=pd) + theme_classic() + ggtitle(sprintf("%s by %s", "KT log-ratio", mvar)) + scale_fill_manual(values=get_color_list(mvar))
	print(p)
	# forest plot
	tmp <- tmp[order(tmp$estimate, decreasing=T),]
	lims <- max(abs(tmp$estimate) + abs(tmp$SE))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(tmp, aes(x=Visit, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=Visit, ymin=estimate-SE, max=estimate+SE, group=metabolite), width=0.2, position=pd) + geom_hline(yintercept=0) + facet_wrap(~contrast) + theme_classic() + ggtitle(sprintf("%s by %s*Visit", "KT log-ratio", mvar)) + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
	print(p)
	p <- ggplot(tmp, aes(x=Visit, y=estimate, group=contrast, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=Visit, ymin=estimate-SE, max=estimate+SE, group=contrast), width=0.2, position=pd) + geom_text(aes(x=Visit, y=-lims, group=contrast, label=contrast), color="black", hjust=0, size=1, position=pd) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("%s by %s*Visit", "KT log-ratio", mvar)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
	print(p)
}
write.table(res, file=sprintf("%s/KTRatio.%s.%s.txt", output_dir, "emmeans", "by_Visit"), quote=F, sep="\t", row.names=F, col.names=T)


## random forest classifying by HIVExposure, separately by Visit
# set weights for ranger's case.weights argument
tab <- table(mapping.sel$HIVExposure)
wt <- 1-tab/sum(tab)
mapping.sel$ipw.HIVExposure <- wt[as.character(mapping.sel$HIVExposure)]
# run ranger
mlevel <- "BIOCHEMICAL"
mvar <- "HIVExposure"
for (visit in visits.selected) {
	# subset data
	mapping.visit <- subset(mapping.sel, Visit==visit)
	data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.visit),]
	samples.visit <- rownames(mapping.visit)
		
	# RF
	response <- mapping.visit[, mvar]; names(response) <- rownames(mapping.visit)
	agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	num_iter <- 1000 # 1000 for full run, 20 for testing
	ncores <- 20
	ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
	out <- mclapply(1:num_iter, function (dummy) {
			importance(ranger(x=data.sel, y=response, num.trees=10000, importance="permutation", seed=ranger.seeds[dummy], num.threads=1))
	}, mc.cores=ncores )	
	collated.importance <- do.call(cbind, out)
	out <- mclapply(1:num_iter, function (dummy) {
			rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
		}, mc.cores=ncores )
	collated.cv <- do.call(cbind, out)

	write.table(collated.importance, file=sprintf("%s/ranger_METABOLITE.%s.%s.importance.txt", output_dir, "Exposed_vs_Unexposed", visit), quote=F, sep="\t", row.names=T, col.names=F)
	write.table(collated.cv, file=sprintf("%s/ranger_METABOLITE.%s.%s.cv.txt", output_dir, "Exposed_vs_Unexposed", visit), quote=F, sep="\t", row.names=T, col.names=F)
	## END BLOCK TO COMMENT ##

	collated.importance <- read.table(sprintf("%s/ranger_METABOLITE.%s.%s.importance.txt", output_dir, "Exposed_vs_Unexposed", visit), header=F, as.is=T, sep="\t", row.names=1, quote="")
	collated.cv <- read.table(sprintf("%s/ranger_METABOLITE.%s.%s.cv.txt", output_dir, "Exposed_vs_Unexposed", visit), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	inds <- order(importance.mean, decreasing=T)
	inds <- inds[1:min(20, as.numeric(names(cv.mean)[which.min(cv.mean)]))] # minimum CVE, capped at 20 features
	write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger_METABOLITE.%s.%s.features.txt", output_dir, "Exposed_vs_Unexposed", visit), quote=F, sep="\t", row.names=T, col.names=F)

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
	sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, importance="permutation", seed=sample(1:num_iter,1), probability=F)
	save(sparseRanger, file=sprintf("%s/ranger_METABOLITE.%s.%s.model", output_dir, "Exposed_vs_Unexposed", visit))
	load(sprintf("%s/ranger_METABOLITE.%s.%s.model", output_dir, "Exposed_vs_Unexposed", visit))
	# accuracy of final sparseRF model
	pred <- predictions(sparseRanger)
	pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
	pred_df_out <- merge(pred_df, data.sel, by="row.names")
	write.table(pred_df_out, file=sprintf("%s/ranger_METABOLITE.%s.%s.predictions.txt", output_dir, "Exposed_vs_Unexposed", visit), quote=F, sep="\t", row.names=F, col.names=T)
	confusion_matrix <- table(pred_df[, c("true", "predicted")])
	class_errors <- unlist(lapply(levels(mapping.visit[, mvar]), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.visit[,mvar])
	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
	mccvalue <- mcc(vec.pred, vec.true)
	cf <- confusionMatrix(confusion_matrix, positive=sprintf("Exposed"))
	cflabel <- sprintf("Positive: %s Sens: %.4g  Spec: %.4g\n PPV: %.4g  NPV: %.4g", cf$positive, cf$byClass[["Sensitivity"]], cf$byClass[["Specificity"]], cf$byClass[["Pos Pred Value"]], cf$byClass[["Neg Pred Value"]])
	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s) (accuracy = %.2f%%, MCC = %.4f)", "Exposed_vs_Unexposed", visit, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + annotate(geom="text", x=2.5, y=10, label=cflabel)
	print(p)

	write.table(confusion_matrix, file=sprintf("%s/ranger_METABOLITE.%s.%s.confusion_matrix.txt", output_dir, "Exposed_vs_Unexposed", visit), quote=F, sep="\t", row.names=T, col.names=T)
	## END BLOCK TO COMMENT ##

	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", "Exposed_vs_Unexposed", visit)))
	# plotting - per-group variables
	df <- data.frame(OTU=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	df$metabolite_name <- as.character(df$OTU)
	if (mlevel == "BIOCHEMICAL") {
		df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
		df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
		df$OTU_string <- sprintf("%s (%s; %s)", df$OTU, df$subpathway, df$superpathway)
	} else if (mlevel == "SUB.PATHWAY") {
		df$subpathway <- df$metabolite_name
		df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
		df$OTU_string <- sprintf("%s", df$OTU)
	} else {
		df$superpathway <- df$metabolite_name
		df$OTU_string <- sprintf("%s", df$OTU)
	}
	df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
	# load effect sizes from linear regression
	lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, "All_vs_HUU", "all_visits"), header=T, as.is=T, sep="\t", quote=""); lmres <- subset(lmres, contrast=="Exposed - Unexposed" & Visit==visit)
	rownames(lmres) <- lmres$metabolite
	for (lmvar in c("estimate", "SE", "padj", "dir")) {
		df[,lmvar] <- lmres[df$metabolite_name, lmvar]
	}
	p <- ggplot(df, aes(x=OTU, y=importance, label=OTU, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=OTU, y=0, label=OTU_string), size=3, hjust=0) + ggtitle(sprintf("%s features (%s)", "Exposed_vs_Unexposed", visit)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
	print(p)
	lims <- max(abs(df$estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=OTU, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=OTU, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_tile(aes(x=OTU, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s", "Exposed_vs_Unexposed", visit)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
	print(p)
	# store detailed RF feature data
	feature.detailed <- data.frame(metabolite=factor(names(importance.mean)[inds], levels=names(importance.mean)[inds]), importance=importance.mean[inds], sd=importance.sd[inds])
	feature.detailed <- merge(feature.detailed, lmres, by="metabolite", sort=F)
	write.table(feature.detailed, file=sprintf("%s/ranger_METABOLITE.%s.%s.features_detailed.txt", output_dir, "Exposed_vs_Unexposed", visit), quote=F, sep="\t", row.names=F, col.names=T)
	# shading rectangles of importance values
	df.rect <- df
	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	p <- ggplot(df.rect, aes(x=x, y=OTU, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("RF importance (%s, %s)", "Exposed_vs_Unexposed", visit)) + scale_fill_gradient(low="white", high="black")
	print(p)
	# violin plots of metabolite values
	agg.melt <- agg.melt.stored
	agg.melt[,mvar] <- mapping.visit[agg.melt$SampleID, mvar]
	agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt[,mvar]))
	agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
	agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "Exposed_vs_Unexposed", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "Exposed_vs_Unexposed", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	npages <- n_pages(p)
	for (ip in 1:npages) {
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "Exposed_vs_Unexposed", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
		print(p)
	}
	agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
	p <- ggplot(agg.melt2, aes_string(x="metabolite", y="value", color=mvar)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "Exposed_vs_Unexposed", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)
	agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
	agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
	p <- ggplot(agg.melt2, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "Exposed_vs_Unexposed", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)
	p <- ggplot(agg.melt2, aes_string(x="metabolite", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "Exposed_vs_Unexposed", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)
}


#########################################################################################################
### compare HEU to HUU
mapping.sel <- subset(mapping, Group3 %in% c("HEU", "HUU"))
mapping.sel$Group3 <- droplevels(mapping.sel$Group3)

## sample diagram
mvar <- "Group3"
agg <- aggregate(Visit ~ ID1, mapping.sel, function(x) max(as.numeric(x))); agg <- agg[order(agg$Visit),]
mapping.sel$ID1 <- factor(mapping.sel$ID1, levels=agg$ID1)
tab <- table(mapping.sel$Visit)
p <- ggplot(mapping.sel, aes_string(x="Visit", y="ID1", group="ID1", color=mvar)) + geom_point() + geom_line() + facet_wrap(as.formula(sprintf("~%s", mvar)), ncol=1, scales="free_y") + theme_classic() + ggtitle(sprintf("HEU vs HUU (n=%d samples from %d PIDs)", nrow(mapping.sel), length(unique(mapping.sel$ID1)))) + scale_color_manual(values=get_color_list(mvar))
print(p)
tab <- table(mapping.sel[,c(mvar, "Visit")])
df <- melt(tab)
p <- ggplot(df, aes_string(x="Visit", y=mvar, size="value", color=mvar)) + geom_point(aes(size=2*value), shape=1) + geom_text(aes(label=value), hjust=0.5, size=2) + scale_size_area() + theme_classic() + ggtitle(sprintf("HEU vs HUU (n=%d samples from %d PIDs)", nrow(mapping.sel), length(unique(mapping.sel$ID1)))) + scale_color_manual(values=get_color_list(mvar)) + scale_size_continuous(range = c(1,16))
print(p)


## tableOne
demo <- unique(mapping.sel[, c("ID1", demo_vars)])
tab1 <- CreateTableOne(vars=demo_vars, strata=mvar, data=demo, smd=T)
write.table(print(tab1, noSpaces=T), file=sprintf("%s/Table_1.%s.txt", output_dir, "HEU_vs_HUU"), quote=F, sep="\t", row.names=T, col.names=T)


## linear regression (stratified by Visit)
mlevel <- "BIOCHEMICAL"; mvar <- "Group3"; data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
tmp <- apply(data.sel, 2, function(x) length(which(x > min(x, na.rm=T)))); to_remove <- names(which(tmp <= 1)); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites that are detected in < 2 samples
df <- merge(data.sel, mapping.sel, by="row.names"); 
res <- {}
for (metabolite in colnames(data.sel)) {
	mod <- lmer(as.formula(sprintf("%s ~ %s*Visit + (1 | ID1)", metabolite, mvar)), data=df); modelstr <- "LM"
	emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | Visit", mvar)), adjust="none")
	tmp <- as.data.frame(emm$contrasts) 
	tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr; tmp$weighted <- enable_weights
	res <- rbind(res, tmp)
}
res <- res[,c("metabolite", setdiff(colnames(res), "metabolite"))]
res$padj <- p.adjust(res$p.value, method="fdr")
res <- res[order(res$p.value, decreasing=F),]
res$dir <- ifelse(res$padj < siglevel, ifelse(sign(res$estimate)==1, "up", "down"), "NS")
res$exp_estimate <- exp(res$estimate)
for (addvar in c("COMP_ID", "HMDB", "KEGG", "PUBCHEM", "PLATFORM", "SUB.PATHWAY", "SUPER.PATHWAY")) {
	res[, addvar] <- metabolon_map[as.character(res$metabolite), addvar]
}
write.table(res, file=sprintf("%s/emmeans.%s.%s.txt", output_dir, "HEU_vs_HUU", "all_visits"), quote=F, sep="\t", row.names=F, col.names=T)
# volcano plots
res <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, "HEU_vs_HUU", "all_visits"), header=T, as.is=T, sep="\t", quote="", comment.char="")
res$contrast <- factor(res$contrast); res$Visit <- droplevels(factor(res$Visit, levels=visits))
for (mv in unique(res$metadata_variable)) {
	df <- subset(res, metadata_variable==mv)
	for (ct in levels(df$contrast)) {
		df2 <- subset(df, contrast==ct & !is.na(estimate))
		lims <- max(abs(df2$estimate), na.rm=T)
		pl <- list()
		for (visit in levels(df2$Visit)) {
			df3 <- subset(df2, Visit==visit)
			df3$siglabel <- ifelse(df3$dir=="NS", ifelse(df3$metabolite %in% always_label, df3$metabolite, NA), df3$metabolite)
			df3$dir <- ifelse(df3$dir == "NS", ifelse(df3$metabolite %in% always_label, "manual", "NS"), df3$dir)
			p <- ggplot(df3, aes(x=estimate, y=-log10(padj), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=2, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s, %s", ct, visit)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
			pl[[length(pl)+1]] <- p
		}
		multiplot(plotlist=pl, cols=2, rows=2)
	}
}


## linear regression (averaged over Visit)
mlevel <- "BIOCHEMICAL"; data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
name_map <- data.frame(original=colnames(data.sel), valid=make.names(colnames(data.sel))); rownames(name_map) <- name_map$valid; colnames(data.sel) <- make.names(colnames(data.sel))
sds <- apply(data.sel, 2, sd); to_remove <- names(which(sds==0)); data.sel <- data.sel[, setdiff(colnames(data.sel), to_remove)] # remove metabolites with zero variance
df <- merge(data.sel, mapping.sel, by="row.names"); 
res <- {}
for (mvar in c("Group3")) {
	for (metabolite in colnames(data.sel)) {
		mod <- lmer(as.formula(sprintf("%s ~ %s*Visit + (1 | ID1)", metabolite, mvar)), data=df); modelstr <- "LM"
		emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s", mvar)), adjust="none")
		tmp <- as.data.frame(emm$contrasts)
		tmp$metabolite <- name_map[metabolite, "original"]; tmp$metadata_variable <- mvar; tmp$model <- modelstr; tmp$weighted <- enable_weights
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
write.table(res, file=sprintf("%s/emmeans.%s.%s.txt", output_dir, "HEU_vs_HUU", "averaged_over_Visit"), quote=F, sep="\t", row.names=F, col.names=T)
# volcano plots and rank-order plots
res <- read.table(sprintf("%s/emmeans.%s.%s.txt", output_dir, "HEU_vs_HUU", "averaged_over_Visit"), header=T, as.is=T, sep="\t", quote="", comment.char="")
res$contrast <- factor(res$contrast)
for (mv in unique(res$metadata_variable)) {
	df <- subset(res, metadata_variable==mv); df$contrast <- droplevels(df$contrast)
	for (ct in levels(df$contrast)) {
		df2 <- subset(df, contrast==ct & !is.na(estimate))
		lims <- max(abs(df2$estimate), na.rm=T)
		df3 <- df2
		df3$siglabel <- ifelse(df3$dir=="NS", ifelse(df3$metabolite %in% always_label, df3$metabolite, NA), df3$metabolite)
		df3$dir <- ifelse(df3$dir == "NS", ifelse(df3$metabolite %in% always_label, "manual", "NS"), df3$dir)
		df3$padj <- pmax(df3$padj, 0.001) # censor at padj < 0.001 to improve plotting
		p <- ggplot(df3, aes(x=estimate, y=-log10(padj), color=dir)) + geom_point() + geom_text_repel(aes(label=siglabel), size=1.5, max.overlaps=15) + theme_classic() + ggtitle(sprintf("%s, averaged over Visit", ct)) + geom_hline(yintercept=-log10(siglevel)) + scale_color_manual(values=dircolors) + xlim(c(-lims, lims)) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
		df3 <- df3[order(df3$estimate),]
		df3$rank <- 1:nrow(df3)
		p <- ggplot(df3, aes(x=rank, y=estimate, color=dir)) + geom_point(size=2) + geom_text_repel(aes(label=siglabel), size=1.0, max.overlaps=35) + theme_classic() + ggtitle(sprintf("%s, averaged over Visit", ct)) + scale_color_manual(values=dircolors) + theme(title=element_text(size=10), axis.text=element_text(size=8, color="black"), axis.title=element_text(size=8))
		print(p)
	}
}


## kynurenine to tryptophan ratio
data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.sel),]
df <- merge(data.sel, mapping.sel, by="row.names"); 
df$KTratio <- exp(df[, "kynurenine"]) / exp(df[, "tryptophan"])
df$KTlogratio <- log10(df$KTratio)
# by emmeans
res <- {}
for (mvar in c("HIVExposure")) {
	mod <- lmer(as.formula(sprintf("%s ~ %s*Visit + (1 | ID1)", "KTlogratio", mvar)), data=df); modelstr <- "LM"
	emm <- emmeans(mod, as.formula(sprintf("trt.vs.ctrl ~ %s | Visit", mvar)), adjust="none")
	tmp <- as.data.frame(emm$contrasts)
	tmp$metabolite <- "KTlogratio"; tmp$metadata_variable <- mvar; tmp$model <- modelstr
	tmp <- subset(tmp, !is.na(estimate))
	tmp <- tmp[,c("metabolite", setdiff(colnames(tmp), "metabolite"))]
	tmp$padj <- p.adjust(tmp$p.value, method="fdr")
	tmp <- tmp[order(tmp$p.value, decreasing=F),]
	tmp$dir <- ifelse(tmp$padj < siglevel, ifelse(sign(tmp$estimate)==1, "up", "down"), "NS")
	res <- rbind(res, tmp)
	# boxplot
	pd <- position_dodge(width=0.8)
	p <- ggplot(df, aes_string(x="Visit", y="KTlogratio", fill=mvar)) + geom_boxplot(position=pd) + theme_classic() + ggtitle(sprintf("%s by %s %s", "KT log-ratio", mvar, "HEU_vs_HUU")) + scale_fill_manual(values=get_color_list(mvar))
	print(p)
	# forest plot
	tmp <- tmp[order(tmp$estimate, decreasing=T),]
	lims <- max(abs(tmp$estimate) + abs(tmp$SE))*1.0
	pd <- position_dodge(0.8)
	p <- ggplot(tmp, aes(x=Visit, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=Visit, ymin=estimate-SE, max=estimate+SE, group=metabolite), width=0.2, position=pd) + geom_hline(yintercept=0) + facet_wrap(~contrast) + theme_classic() + ggtitle(sprintf("%s by %s*Visit %s", "KT log-ratio", mvar, "HEU_vs_HUU")) + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
	print(p)
	p <- ggplot(tmp, aes(x=Visit, y=estimate, group=contrast, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=Visit, ymin=estimate-SE, max=estimate+SE, group=contrast), width=0.2, position=pd) + geom_text(aes(x=Visit, y=-lims, group=contrast, label=contrast), color="black", hjust=0, size=1, position=pd) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("%s by %s*Visit %s", "KT log-ratio", mvar, "HEU_vs_HUU")) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
	print(p)
}
write.table(res, file=sprintf("%s/KTRatio.%s.%s.txt", output_dir, "emmeans", "HEU_vs_HUU"), quote=F, sep="\t", row.names=F, col.names=T)


## line plots of each metabolite over time
mvars.plot <- c("HIVExposure")
agg.melt <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt) <- c("SampleID", "metabolite", "value")
data.raw <- df.metabolon[["BIOCHEMICAL_raw"]][rownames(mapping.sel),]
agg.raw <- melt(as.matrix(data.raw), as.is=T); colnames(agg.raw) <- c("SampleID", "metabolite", "value")
agg.raw$Detected <- factor(!is.na(agg.raw$value))
data.peak <- df.metabolon[["BIOCHEMICAL_peak"]][rownames(mapping.sel),]
agg.peak <- melt(as.matrix(data.peak), as.is=T); colnames(agg.peak) <- c("SampleID", "metabolite", "value")
agg.peak$Detected <- factor(!is.na(agg.peak$value))
for (m in c("ID1", "Group2", "Group3", "HIVExposure", "Visit", "log10vlrna", "CD4C", "bm_vl_log10")) {
	agg.melt[,m] <- mapping.sel[agg.melt$SampleID, m]
	agg.raw[,m] <- mapping.sel[agg.raw$SampleID, m]
	agg.peak[,m] <- mapping.sel[agg.peak$SampleID, m]
}


## for each Visit:
## - RF
mlevel <- "BIOCHEMICAL"
set.seed(nrow(mapping.sel))
for (visit in visits.selected) {
	# subset data
	mapping.visit <- subset(mapping.sel, Visit==visit)
	data.sel <- df.metabolon[["BIOCHEMICAL"]][rownames(mapping.visit),]
	samples.visit <- rownames(mapping.visit)
	
	# RF
	response <- mapping.visit[, mvar]; names(response) <- rownames(mapping.visit)
	agg.melt.stored <- melt(as.matrix(data.sel), as.is=T); colnames(agg.melt.stored) <- c("SampleID", "metabolite", "value")
	wt <- 1 - (table(response) / length(response)) # weights for ranger's sample.fraction argument

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	num_iter <- 1000 # 1000 for full run
	ncores <- 20
	ranger.seeds <- sample(1:num_iter, num_iter, replace=T) # set up a vector of seeds for the ranger mclapply
	out <- mclapply(1:num_iter, function (dummy) {
			importance(ranger(x=data.sel, y=response, num.trees=10000, sample.fraction=wt, importance="permutation", seed=ranger.seeds[dummy], num.threads=1))
	}, mc.cores=ncores )	
	collated.importance <- do.call(cbind, out)
	out <- mclapply(1:num_iter, function (dummy) {
			rgcv2(trainx=data.sel, trainy=response, cv.fold=10, step=0.9, num.threads=1)$error.cv
		}, mc.cores=ncores )
	collated.cv <- do.call(cbind, out)

	write.table(collated.importance, file=sprintf("%s/ranger_METABOLITE.%s.%s.importance.txt", output_dir, "HEU_vs_HUU", visit), quote=F, sep="\t", row.names=T, col.names=F)
	write.table(collated.cv, file=sprintf("%s/ranger_METABOLITE.%s.%s.cv.txt", "HEU_vs_HUU", output_dir, visit), quote=F, sep="\t", row.names=T, col.names=F)
	## END BLOCK TO COMMENT ##

	collated.importance <- read.table(sprintf("%s/ranger_METABOLITE.%s.%s.importance.txt", output_dir, "HEU_vs_HUU", visit), header=F, as.is=T, sep="\t", row.names=1, quote="")
	collated.cv <- read.table(sprintf("%s/ranger_METABOLITE.%s.%s.cv.txt", output_dir, "HEU_vs_HUU", visit), header=F, as.is=T, sep="\t", row.names=1)
	importance.mean <- rowMeans(collated.importance)
	importance.sd <- unlist(apply(collated.importance, 1, sd))
	cv.mean <- rowMeans(collated.cv)
	cv.sd <- unlist(apply(collated.cv, 1, sd))
	inds <- order(importance.mean, decreasing=T)
	inds <- inds[1:min(20, as.numeric(names(cv.mean)[which.min(cv.mean)]))] # minimum CVE, capped at 20 features
	write.table(melt(importance.mean[inds]), file=sprintf("%s/ranger_METABOLITE.%s.%s.features.txt", output_dir, "HEU_vs_HUU", visit), quote=F, sep="\t", row.names=T, col.names=F)

	## after running for the first time, COMMENT OUT THIS BLOCK ##
	# using a sparse model with N predictors
	sparseRanger <- ranger(x=data.sel[, names(importance.mean[inds]), drop=F], y=response, num.trees=10000, sample.fraction=wt, importance="permutation", seed=sample(1:num_iter,1), probability=F)
	save(sparseRanger, file=sprintf("%s/ranger_METABOLITE.%s.%s.model", output_dir, "HEU_vs_HUU", visit))
	load(sprintf("%s/ranger_METABOLITE.%s.%s.model", output_dir, "HEU_vs_HUU", visit))
	# accuracy of final sparseRF model
	pred <- predictions(sparseRanger)
	pred_df <- data.frame(SampleID=names(response), predicted=pred, true=response, stringsAsFactors=F); pred_df$predicted <- factor(pred_df$predicted, levels=levels(pred_df$true))
	pred_df_out <- merge(pred_df, data.sel, by="row.names")
	write.table(pred_df_out, file=sprintf("%s/ranger_METABOLITE.%s.%s.predictions.txt", output_dir, "HEU_vs_HUU", visit), quote=F, sep="\t", row.names=F, col.names=T)
	confusion_matrix <- table(pred_df[, c("true", "predicted")])
	class_errors <- unlist(lapply(levels(mapping.visit[, mvar]), function(x) 1-(confusion_matrix[x,x] / sum(confusion_matrix[x,])) )); names(class_errors) <- levels(mapping.visit[,mvar])
	accuracy <- 100*(sum(diag(confusion_matrix)) / sum(confusion_matrix))
	vec.pred <- as.numeric(pred_df$predicted)-1; vec.true <- as.numeric(pred_df$true)-1
	mccvalue <- mcc(vec.pred, vec.true)
	cf <- confusionMatrix(confusion_matrix, positive=sprintf("HUU"))
	cflabel <- sprintf("Positive: %s Sens: %.4g  Spec: %.4g\n PPV: %.4g  NPV: %.4g", cf$positive, cf$byClass[["Sensitivity"]], cf$byClass[["Specificity"]], cf$byClass[["Pos Pred Value"]], cf$byClass[["Neg Pred Value"]])
	df <- cbind(confusion_matrix, class_errors[rownames(confusion_matrix)])
	p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("confusion matrix (%s, %s) (accuracy = %.2f%%, MCC = %.4f)", "HEU_vs_HUU", visit, accuracy, mccvalue)) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + annotate(geom="text", x=2.5, y=10, label=cflabel)
	print(p)

	write.table(confusion_matrix, file=sprintf("%s/ranger_METABOLITE.%s.%s.confusion_matrix.txt", output_dir, "HEU_vs_HUU", visit), quote=F, sep="\t", row.names=T, col.names=T)
	## END BLOCK TO COMMENT ##

	# plotting - per-group sparse model
	df <- data.frame(m=cv.mean, sd=cv.sd, numvar=as.numeric(names(cv.mean)))
	colnames(df) <- c("CV_error", "CV_stddev", "num_variables")
	print(ggplot(df, aes(x=num_variables, y=CV_error)) + geom_errorbar(aes(ymin=CV_error-CV_stddev, ymax=CV_error+CV_stddev), width=.1) + geom_line() + geom_point() + ggtitle(sprintf("Model selection - %s %s", "HEU_vs_HUU", visit)))
	# plotting - per-group variables
	df <- data.frame(feature=factor(names(importance.mean)[inds], levels=rev(names(importance.mean)[inds])), importance=importance.mean[inds], sd=importance.sd[inds])
	df$metabolite_name <- as.character(df$feature)
	if (mlevel == "BIOCHEMICAL") {
		df$subpathway <- metabolon_map[df$metabolite_name, "SUB.PATHWAY"]
		df$superpathway <- metabolon_map[df$metabolite_name, "SUPER.PATHWAY"]
		df$feature_string <- sprintf("%s (%s; %s)", df$feature, df$subpathway, df$superpathway)
	} else if (mlevel == "SUB.PATHWAY") {
		df$subpathway <- df$metabolite_name
		df$superpathway <- metabolon_map[match(df$metabolite_name, metabolon_map$SUB.PATHWAY), "SUPER.PATHWAY"]
		df$feature_string <- sprintf("%s", df$feature)
	} else {
		df$superpathway <- df$metabolite_name
		df$feature_string <- sprintf("%s", df$feature)
	}
	df$superpathway[which(is.na(df$superpathway))] <- "NOT_METABOLITE"
	# load effect sizes from linear regression
	lmres <- read.table(sprintf("%s/emmeans.%s.%s.txt", "HEU_vs_HUU", "all_visits"), header=T, as.is=T, sep="\t", quote="")
	lmres <- subset(lmres, Visit==visit)
	rownames(lmres) <- lmres$metabolite
	for (lmvar in c("estimate", "SE", "padj", "dir")) {
		df[,lmvar] <- lmres[df$metabolite_name, lmvar]
	}
	p <- ggplot(df, aes(x=feature, y=importance, label=feature, fill=superpathway)) + geom_bar(position=position_dodge(), stat="identity", color=NA) + geom_errorbar(aes(ymin=importance-sd, ymax=importance+sd), width=.2, position=position_dodge(.9)) + coord_flip() + geom_text(aes(x=feature, y=0, label=feature_string), size=3, hjust=0) + ggtitle(sprintf("%s features (%s)", "HEU_vs_HUU", visit)) + scale_fill_manual(values=cols.superpathway) + theme(axis.text.y=element_blank())
	print(p)
	lims <- max(abs(df$estimate) + abs(df$SE), na.rm=T)*1.0; pd <- position_dodge(0.8)
	p <- ggplot(df, aes(x=feature, y=estimate, color=dir)) + geom_point(position=pd) + geom_errorbar(aes(x=feature, ymin=estimate-SE, max=estimate+SE), width=0.2, position=pd) + geom_tile(aes(x=feature, y=-lims*0.95, fill=importance), height=0.1, inherit.aes=F) + geom_hline(yintercept=0) + theme_classic() + ggtitle(sprintf("LM estimates %s %s", "HEU_vs_HUU", visit)) + coord_flip() + scale_color_manual(values=dircolors) + ylim(c(-lims, lims)) + scale_fill_gradient(low="white", high="black")
	print(p)
	# store detailed RF feature data
	feature.detailed <- data.frame(metabolite=factor(names(importance.mean)[inds], levels=names(importance.mean)[inds]), importance=importance.mean[inds], sd=importance.sd[inds])
	feature.detailed <- merge(feature.detailed, lmres, by="metabolite", sort=F)
	write.table(feature.detailed, file=sprintf("%s/ranger_METABOLITE.%s.%s.features_detailed.txt", output_dir, "HEU_vs_HUU", visit), quote=F, sep="\t", row.names=F, col.names=T)
	# shading rectangles of importance values
	df.rect <- df
	df.rect$x <- 1; df.rect$y <- 1:nrow(df.rect)
	p <- ggplot(df.rect, aes(x=x, y=feature, fill=importance)) + geom_tile() + theme_classic() + ggtitle(sprintf("RF importance (%s, %s)", "HEU_vs_HUU", visit)) + scale_fill_gradient(low="white", high="black")
	print(p)
	# violin plots of metabolite values
	agg.melt <- agg.melt.stored
	agg.melt[,mvar] <- mapping.visit[agg.melt$SampleID, mvar]
	agg.melt$Prediction <- ordered(as.character(pred_df[agg.melt$SampleID, "predicted"]), levels=levels(agg.melt[,mvar]))
	agg.melt <- subset(agg.melt, metabolite %in% rownames(df))
	agg.melt$metabolite <- factor(agg.melt$metabolite, levels=rownames(df))
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="free", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "HEU_vs_HUU", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)
	p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "HEU_vs_HUU", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	npages <- n_pages(p)
	for (ip in 1:npages) {
		p <- ggplot(agg.melt, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap_paginate(~metabolite, scales="free", ncol=3, nrow=4, page=ip) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "HEU_vs_HUU", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
		print(p)
	}
	agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway=="NOT_METABOLITE")))
	p <- ggplot(agg.melt2, aes_string(x="metabolite", y="value", color=mvar)) + geom_violin() + geom_point(position=position_dodge(width=0.9)) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "HEU_vs_HUU", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)
	agg.melt2 <- subset(agg.melt, metabolite %in% rownames(subset(df, superpathway!="NOT_METABOLITE")))
	agg.melt2$metabolite <- factor(as.character(agg.melt2$metabolite), levels=rev(levels(agg.melt2$metabolite)))
	p <- ggplot(agg.melt2, aes_string(x=mvar, y="value", color=mvar)) + geom_violin() + geom_point() + facet_wrap(~metabolite, scales="fixed", ncol=3) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "HEU_vs_HUU", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)
	p <- ggplot(agg.melt2, aes_string(x="metabolite", y="value", color=mvar)) + geom_violin(position=position_dodge(width=0.7)) + geom_point(position=position_dodge(width=0.7), size=1) + theme_classic() + ggtitle(sprintf("Rel. abund. of RF metabolites (%s, %s)", "HEU_vs_HUU", visit)) + coord_flip() + scale_color_manual(values=get_color_list(mvar))
	print(p)

}



#########################################################################################################
### Combined RF model plots for various comparisons

## Exposed_vs_Unexposed (by HIVExposure)
res <- {}
for (visit in visits.selected) {
	tmp <- read.table(sprintf("%s/ranger_METABOLITE.%s.%s.features_detailed.txt", output_dir, "Exposed_vs_Unexposed", visit), header=T, as.is=T, sep="\t", quote="", comment.char="")
	res <- rbind(res, tmp)
}
write.table(res, file=sprintf("%s/ranger_METABOLITE.%s.all_visits.features_detailed.txt", output_dir, "Exposed_vs_Unexposed"), quote=F, col.names=T, row.names=F, sep="\t")

cex.main <- par("cex.main")
par(cex.main=0.7)
df <- dcast(metabolite ~ Visit, value.var="importance", data=res)
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("RF importance values (%s, %s)", "Exposed_vs_Unexposed", "all visits"))

df <- dcast(metabolite ~ Visit, value.var="estimate", data=res)
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("blue", "white", "red"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("LM estimates from RF features (%s, %s)", "Exposed_vs_Unexposed", "all visits"))

# only show features detected in at least 2 visits
selected_features <- names(which(table(res$metabolite)>1))
df <- dcast(metabolite ~ Visit, value.var="importance", data=subset(res, metabolite %in% selected_features))
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("RF importance values (%s, %s)", "Exposed_vs_Unexposed", "all visits, no singleton features"))

df <- dcast(metabolite ~ Visit, value.var="estimate", data=subset(res, metabolite %in% selected_features))
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
hm <- heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("blue", "white", "red"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("LM estimates from RF features (%s, %s)", "Exposed_vs_Unexposed", "all visits, no singleton features"))
par(cex.main=cex.main)

# as a circle plot
df2 <- melt(df); colnames(df2) <- c("feature", "Visit", "estimate")
df2$feature <- factor(df2$feature, levels=labels(hm$rowDendrogram))
df2$dir <- ifelse(df2$estimate==0, "NS", ifelse(sign(df2$estimate)==1, "up", "down"))
p <- ggplot(df2, aes(x=Visit, y=feature, fill=estimate, color=dir, group=feature)) + geom_point(shape=21, size=4) + theme_classic() + ggtitle(sprintf("LM estimates from RF features (%s, %s)", "Exposed_vs_Unexposed", "all visits, no singleton features")) + scale_fill_gradient2(low="blue", mid="white", high="red") + scale_color_manual(values=dircolors)
print(p)


## HEU vs HUU
## combined plot of all RF models (across Visit)
res <- {}
for (visit in visits.selected) {
	tmp <- read.table(sprintf("%s/ranger_METABOLITE.%s.%s.features_detailed.txt", output_dir, "HEU_vs_HUU", visit), header=T, as.is=T, sep="\t", quote="", comment.char="")
	res <- rbind(res, tmp)
}

write.table(res, file=sprintf("%s/ranger_METABOLITE.%s.all_visits.features_detailed.txt", output_dir, "HEU_vs_HUU"), quote=F, col.names=T, row.names=F, sep="\t")

cex.main <- par("cex.main")
par(cex.main=0.7)

df <- dcast(metabolite ~ Visit, value.var="importance", data=res)
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("RF importance values (%s, %s)", "HEU_vs_HUU", "all visits"))

df <- dcast(metabolite ~ Visit, value.var="estimate", data=res)
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("blue", "white", "red"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("LM estimates from RF features (%s, %s)", "HEU_vs_HUU", "all visits"))

# only show features detected in at least 2 visits
selected_features <- names(which(table(res$metabolite)>1))
df <- dcast(metabolite ~ Visit, value.var="importance", data=subset(res, metabolite %in% selected_features))
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("RF importance values (%s, %s)", "HEU_vs_HUU", "all visits, no singleton features"))

df <- dcast(metabolite ~ Visit, value.var="estimate", data=subset(res, metabolite %in% selected_features))
rownames(df) <- df$metabolite; df <- df[,-1]
df <- df[, visits.selected]
df <- as.matrix(df)
df[which(is.na(df), arr.ind=T)] <- 0
hm <- heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("blue", "white", "red"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.4, main=sprintf("LM estimates from RF features (%s, %s)", "HEU_vs_HUU", "all visits, no singleton features"))

par(cex.main=cex.main)

# as a circle plot
df2 <- melt(df); colnames(df2) <- c("feature", "Visit", "estimate")
df2$feature <- factor(df2$feature, levels=labels(hm$rowDendrogram))
df2$dir <- ifelse(df2$estimate==0, "NS", ifelse(sign(df2$estimate)==1, "up", "down"))
p <- ggplot(df2, aes(x=Visit, y=feature, fill=estimate, color=dir, group=feature)) + geom_point(shape=21, size=4) + theme_classic() + ggtitle(sprintf("LM estimates from RF features (%s, %s)", "HEU_vs_HUU", "all visits, no singleton features")) + scale_fill_gradient2(low="blue", mid="white", high="red") + scale_color_manual(values=dircolors)
print(p)



##########################################################################################################
### analysis of quantitative KT data (n=101 BMK, n=121 plasma)

matrices <- c("BMK", "Plasma")
assay_types <- c("Global", "Quantitative")

## read in quantitative Plasma and BMK data from UCLA-07-23PHTASA
phtasa <- read.table("data/ZEBS/UCLA-07-23PHTASA.txt", header=T, as.is=T, sep="\t")
phtasa2 <- read.table("data/ZEBS/UCLA-02-24PHTASA.txt", header=T, as.is=T, sep="\t", colClasses="character")

## build a combined table of all values
combined <- list()
# quantitative - BMK
tmp <- subset(phtasa, Client.Matrix=="Whole Breast Milk")
sel <- intersect(tmp$SampleID, rownames(mapping))
tmp2 <- dcast(SampleID ~ Analyte.Name, data=tmp, value.var="Calculated.Concentration..ug.mL.")
colnames(tmp2) <- c("SampleID", "kynurenine", "tryptophan")
tmp2$kynurenine <- as.numeric(tmp2$kynurenine); tmp2$tryptophan <- as.numeric(tmp2$tryptophan)
tmp2$KTratio <- tmp2$kynurenine / tmp2$tryptophan
tmp2$KTlogratio <- log10(tmp2$KTratio)
tmp2$AssayType <- "Quantitative"
# global - BMK
mlevel <- "BIOCHEMICAL"
data.sel <- df.metabolon[["BIOCHEMICAL"]][sel,]
df <- merge(data.sel, mapping, by="row.names"); 
df$KTratio <- exp(df[, "kynurenine"]) / exp(df[, "tryptophan"])
df$KTlogratio <- log10(df$KTratio)
tmp3 <- df[,c("Row.names", "kynurenine", "tryptophan", "KTratio", "KTlogratio")]; colnames(tmp3) <- c("SampleID", "kynurenine", "tryptophan", "KTratio", "KTlogratio")
tmp3$AssayType <- "Global"
tostore <- rbind(tmp2, tmp3)
tostore$Matrix <- "BMK"
combined[["BMK"]] <- tostore
# quantitative - Plasma
tmp <- subset(phtasa, Client.Matrix=="Plasma")
tmp$SampleID <- gsub(".Plasma", "", tmp$SampleID)
sel <- intersect(tmp$SampleID, rownames(mapping)); sel.plasma.quant <- sel
tmp2 <- dcast(SampleID ~ Analyte.Name, data=tmp, value.var="Calculated.Concentration..ug.mL.")
colnames(tmp2) <- c("SampleID", "kynurenine", "tryptophan")
tmp2$kynurenine <- as.numeric(tmp2$kynurenine); tmp2$tryptophan <- as.numeric(tmp2$tryptophan)
tmp2$KTratio <- tmp2$kynurenine / tmp2$tryptophan
tmp2$KTlogratio <- log10(tmp2$KTratio)
tmp2$AssayType <- "Quantitative"
tmp <- subset(phtasa2, Client.Matrix=="Plasma")
tmp$SampleID <- tmp$Client.Sample.ID
sel <- intersect(tmp$SampleID, rownames(mapping))
tmp <- subset(tmp, SampleID %in% sel)
tmp2b <- dcast(SampleID ~ Analyte, data=tmp, value.var="Result")
colnames(tmp2b) <- c("SampleID", "kynurenine", "tryptophan")
tmp2b$kynurenine <- as.numeric(tmp2b$kynurenine); tmp2b$tryptophan <- as.numeric(tmp2b$tryptophan)
tmp2b$KTratio <- tmp2b$kynurenine / tmp2b$tryptophan
tmp2b$KTlogratio <- log10(tmp2b$KTratio)
tmp2b$AssayType <- "Quantitative"
tmp2.combined <- rbind(tmp2, tmp2b)

out <- do.call(rbind, combined)
out$HIVExposure <- mapping[as.character(out$SampleID), "HIVExposure"]

## compare Plasma vs BMK
p <- ggplot(mapping) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("Plasma vs BMK"), size=16)
print(p)
res <- {}
res_to_store <- {}
bmkplasmaratio_to_store <- {}
for (assay_type in "Quantitative") {
	# combine plasma and BMK data
	sel <- intersect(combined[["BMK"]]$SampleID, combined[["Plasma"]]$SampleID)
	data.bmk <- subset(combined[["BMK"]], AssayType==assay_type); data.bmk$Matrix <- "BMK"
	data.plasma <- subset(combined[["Plasma"]], AssayType==assay_type); data.plasma$Matrix <- "Plasma"
	data.sel <- rbind(data.bmk, data.plasma)
	df <- merge(data.sel, mapping)
	df$Matrix <- factor(df$Matrix, levels=c("Plasma", "BMK"))
	df$CD4C_binned <- ifelse(df$CD4C > 350, ">350", "<=350")
	df$CD4C_binned2 <- ifelse(df$CD4C > 300, ">300", "<=300")
	df$CD4C_binned3 <- ifelse(df$CD4C > 200, ">200", "<=200")
	
	# identify and remove outliers
	outliers_to_remove <- {}
	for (qvar in c("kynurenine", "tryptophan", "KTratio", "KTlogratio")) {
		df2 <- dcast(SampleID ~ Matrix, data=df, value.var=qvar)
		df2 <- merge(df2, mapping, by="SampleID")
		if (assay_type=="Quantitative") {
			df2$BMKPlasmaRatio <- df2$BMK / df2$Plasma
		} else {
			df2$BMKPlasmaRatio <- exp(df2$BMK) / exp(df2$Plasma)
		}
		# mark outliers as +/- 4 SDs from mean
		for (mtx in c("Plasma", "BMK", "BMKPlasmaRatio")) {
			df3 <- df2 %>% arrange(desc( .data[[mtx]] )) %>% group_by(HIVExposure) %>% mutate(m = mean(.data[[mtx]], na.rm=T), sd=sd(.data[[mtx]], na.rm=T), outlier=(.data[[mtx]] > (m+4*sd) | (.data[[mtx]] < (m-4*sd))), outlier_type=mtx, qvar=qvar) %>% filter(outlier) %>% select(SampleID, HIVExposure, outlier_type, qvar)
			outliers_to_remove <- rbind(outliers_to_remove, df3)
		}
	}
	write.table(outliers_to_remove, file=sprintf("%s/Quantitative_KT.outliers_removed.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)
	# remove outliers
	df <- subset(df, !(SampleID %in% outliers_to_remove$SampleID))

	for (qvar in c("kynurenine", "tryptophan", "KTratio", "KTlogratio")) {	
		# scatterplot of values
		df2 <- dcast(SampleID ~ Matrix, data=df, value.var=qvar)
		df2 <- merge(df2, mapping, by="SampleID")
		df2$CD4C_binned <- ifelse(df2$CD4C > 350, ">350", "<=350")
		df2$CD4C_binned2 <- ifelse(df2$CD4C > 300, ">300", "<=300")
		df2$CD4C_binned3 <- ifelse(df2$CD4C > 200, ">200", "<=200")

		test.spearman <- cor.test(as.formula(sprintf("~%s+%s", "BMK", "Plasma")), data=df2, method="spearman")
		test.pearson <- cor.test(as.formula(sprintf("~%s+%s", "BMK", "Plasma")), data=df2, method="pearson")
		p <- ggplot(df2, aes_string(x="BMK", y="Plasma", color=mvar)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("%s (BMK vs Plasma, %s; rho=%.4g, p=%.4g, r=%.4g, p=%.4g)", qvar, assay_type, test.spearman$estimate, test.spearman$p.value, test.pearson$estimate, test.pearson$p.value)) + scale_color_manual(values=get_color_list(mvar))
		print(p)
		# boxplot of values
		pd <- position_dodge(width=0.8)
		p <- ggplot(df, aes_string(x=mvar, y=qvar, color=mvar)) + geom_boxplot(outlier.shape=NA, position=pd, width=0.6) + geom_beeswarm(size=3, cex=2, dodge.width=0.8) + theme_classic() + facet_wrap(~Matrix) + ggtitle(sprintf("%s by %s+%s (%s)", qvar, mvar, "Matrix", assay_type)) + scale_color_manual(values=get_color_list(mvar))
		print(p)
		p <- ggplot(df, aes_string(x=mvar, y=qvar, color=mvar)) + geom_boxplot(outlier.shape=NA, position=pd, width=0.6) + geom_beeswarm(size=3, cex=2, dodge.width=0.8) + theme_classic() + facet_wrap(~Matrix, scales="free_y") + ggtitle(sprintf("%s by %s+%s (%s)", qvar, mvar, "Matrix", assay_type)) + scale_color_manual(values=get_color_list(mvar))
		print(p)
		# print summary statistics to file
		summ <- df %>% group_by(HIVExposure, Matrix) %>% summarise_at(qvar, list(min=min, median=median, mean=mean, max=max, sd=sd, var=var), na.rm=T)
		summ$qvar <- qvar; summ$AssayType <- assay_type
		res_to_store <- rbind(res_to_store, summ)
		# ratio (BMK/Plasma)
		if (assay_type=="Quantitative") {
			df2$BMKPlasmaRatio <- df2$BMK / df2$Plasma
		} else {
			df2$BMKPlasmaRatio <- exp(df2$BMK) / exp(df2$Plasma)
		}
		df2$qvar <- qvar; df2$AssayType <- assay_type
		bmkplasmaratio_to_store <- rbind(bmkplasmaratio_to_store, df2)
		test <- wilcox.test(as.formula(sprintf("%s ~ %s", "BMKPlasmaRatio", mvar)), data=df2)
		p <- ggplot(df2, aes_string(x=mvar, y="BMKPlasmaRatio", color=mvar)) + geom_boxplot(outlier.shape=NA, position=pd, width=0.6) + geom_beeswarm(size=3, cex=1.5, dodge.width=0.8) + theme_classic() + ggtitle(sprintf("%s BMK/Plasma ratio by %s+%s (%s, wilcoxon p=%.4g)", qvar, mvar, "Matrix", assay_type, test$p.value)) + scale_color_manual(values=get_color_list(mvar))
		print(p)
		df3 <- df2 %>% arrange(desc(BMKPlasmaRatio)) %>% group_by(HIVExposure) %>% mutate(m = mean(BMKPlasmaRatio, na.rm=T), sd=sd(BMKPlasmaRatio, na.rm=T), outlier=(BMKPlasmaRatio > (m+4*sd) | (BMKPlasmaRatio < (m-4*sd)))) # color outliers as +/- 4 SDs from mean
		df$outlier <- df$SampleID %in% subset(df3, outlier)$SampleID
		p <- ggplot(df3, aes_string(x=mvar, y="BMKPlasmaRatio", group=mvar, color="outlier")) + geom_boxplot(outlier.shape=NA, position=pd, width=0.6) + geom_beeswarm(size=3, cex=1.5, dodge.width=0.8) + theme_classic() + ggtitle(sprintf("%s BMK/Plasma ratio by %s+%s (%s, wilcoxon p=%.4g)", qvar, mvar, "Matrix", assay_type, test$p.value)) + scale_color_manual(values=c("black", "red"))
		print(p)
		p <- ggplot(df, aes_string(x="Matrix", y=qvar, color=mvar)) + geom_point() + geom_text_repel(aes(label=SampleID), size=2) + geom_line(aes(group=SampleID), size=0.1, linetype="dotted") + theme_classic() + ggtitle(sprintf("%s by %s+%s (%s)", qvar, mvar, "Matrix", assay_type)) + scale_color_manual(values=get_color_list(mvar))
		print(p)
		p <- ggplot(df, aes_string(x="Matrix", y=qvar, color="outlier")) + geom_point() + geom_text_repel(aes(label=SampleID), size=2) + geom_line(aes(group=SampleID), size=0.1, linetype="dotted") + theme_classic() + ggtitle(sprintf("%s by %s+%s (%s)", qvar, mvar, "Matrix", assay_type)) + scale_color_manual(values=get_color_list(mvar)) + scale_color_manual(values=c("black", "red"))
		print(p)
		p <- ggplot(df, aes_string(x="Matrix", y=qvar)) + geom_point() + geom_text_repel(aes(label=SampleID), size=2) + geom_segment(x=1, xend=2, aes_string(y="Plasma", yend="BMK", group="SampleID", color="BMKPlasmaRatio"), inherit.aes=F, data=df2, size=0.2) + theme_classic() + ggtitle(sprintf("%s by %s+%s (%s)", qvar, mvar, "Matrix", assay_type)) + scale_color_gradient(low="blue", high="red")
		print(p)

		summ <- df2 %>% group_by(HIVExposure) %>% summarise_at("BMKPlasmaRatio", list(min=min, median=median, mean=mean, max=max, sd=sd, var=var), na.rm=T)
		summ$Matrix <- "BMKPlasmaRatio"; summ$qvar <- qvar; summ$AssayType <- assay_type
		res_to_store <- rbind(res_to_store, summ)
	}
}
write.table(res_to_store, file=sprintf("%s/Quantitative_KT.summary_statistics.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)
write.table(bmkplasmaratio_to_store, file=sprintf("%s/Quantitative_KT.BMKPlasmaRatio.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)


## compare global vs quantitative
p <- ggplot(mapping) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("Global vs Quantitative"), size=16)
print(p)
res <- {}
for (mvar in c("HIVExposure")) {
	for (matrix in matrices) {
		data.sel <- combined[[matrix]]
		df <- merge(data.sel, mapping)
		outliers_to_remove <- read.table(sprintf("/Lab_Share/ZEBS/metabolon/Quantitative_KT.outliers_removed.txt"), header=T, colClasses="character", as.is=T, sep="\t")
	# remove outliers
		df <- subset(df, !(SampleID %in% outliers_to_remove$SampleID))
		
		for (qvar in c("kynurenine", "tryptophan", "KTratio", "KTlogratio")) {	
			# scatterplot of values
			df2 <- dcast(SampleID ~ AssayType, data=df, value.var=qvar)
			df2 <- merge(df2, mapping, by="SampleID")
			test.spearman <- cor.test(as.formula(sprintf("~%s+%s", "Global", "Quantitative")), data=df2, method="spearman")
			test.pearson <- cor.test(as.formula(sprintf("~%s+%s", "Global", "Quantitative")), data=df2, method="pearson")
			p <- ggplot(df2, aes_string(x="Global", y="Quantitative", color=mvar)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("%s (global vs quantitative, %s; rho=%.4g, p=%.4g, r=%.4g, p=%.4g)", qvar, matrix, test.spearman$estimate, test.spearman$p.value, test.pearson$estimate, test.pearson$p.value)) + scale_color_manual(values=get_color_list(mvar))
			print(p)
			p <- ggplot(df2, aes_string(x="Global", y="Quantitative")) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("%s (global vs quantitative, %s; rho=%.4g, p=%.4g, r=%.4g, p=%.4g)", qvar, matrix, test.spearman$estimate, test.spearman$p.value, test.pearson$estimate, test.pearson$p.value)) + scale_color_manual(values=get_color_list(mvar))
			print(p)
			# boxplot of values
			pd <- position_dodge(width=0.8)
			p <- ggplot(df, aes_string(x=mvar, y=qvar, color=mvar)) + geom_boxplot(outlier.shape=NA, position=pd, width=0.6) + geom_beeswarm(dodge.width=0.8, size=3, cex=2) + theme_classic() + facet_wrap(~AssayType) + ggtitle(sprintf("%s by %s+%s (%s)", qvar, mvar, "AssayType", matrix)) + scale_color_manual(values=get_color_list(mvar))
			print(p)
			p <- ggplot(df, aes_string(x=mvar, y=qvar, color=mvar)) + geom_boxplot(outlier.shape=NA, position=pd, width=0.6) + geom_beeswarm(dodge.width=0.8, size=3, cex=2) + theme_classic() + facet_wrap(~AssayType, scales="free_y") + ggtitle(sprintf("%s by %s+%s (%s)", qvar, mvar, "AssayType", matrix)) + scale_color_manual(values=get_color_list(mvar))
			print(p)
			# t-test on quantitative data
			if (nlevels(df2[,mvar])==2) {
				test <- t.test(as.formula(sprintf("%s ~ %s", "Quantitative", mvar)), data=df2)
			} else {
				test <- kruskal.test(as.formula(sprintf("%s ~ %s", "Quantitative", mvar)), data=df2)
			}
			res <- rbind(res, c(qvar, mvar, matrix, test$statistic, test$p.value, test$method))
		}
	}
}
res <- as.data.frame(res)
colnames(res) <- c("metabolite", "metadata_variable", "matrix", "statistic", "p.value", "method")
res$p.value <- as.numeric(as.character(res$p.value))
res$padj <- p.adjust(res$p.value, method="fdr")
write.table(res, file=sprintf("%s/Quantitative_KT.all_stats.txt", output_dir), quote=F, sep="\t", row.names=F, col.names=T)


## compare vs CD4, log10vlrna, bm_vl_log10
p <- ggplot(mapping) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("CD4 and log10vlrna"), size=16)
print(p)
group_var <- "HIVExposure"
for (qvar in c("kynurenine", "tryptophan", "KTratio", "KTlogratio")) {
	for (matrix in matrices) {
		data.sel <- combined[[matrix]]
		df <- merge(data.sel, mapping)
		df[which(df$Group2=="HUU"), "log10vlrna"] <- NA # NA out HUU viral loads
		for (assay_type in assay_types) {
			df2 <- subset(df, AssayType==assay_type)
			pl <- list()
			# scatterplot of values
			for (mvar in c("CD4C", "log10vlrna", "bm_vl_log10")) {
				test.spearman <- cor.test(as.formula(sprintf("~%s+%s", qvar, mvar)), data=df2, method="spearman")
				p <- ggplot(df2, aes_string(x=qvar, y=mvar)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("%s vs. %s, %s %s (rho=%.4g, p=%.4g)", qvar, mvar, matrix, assay_type, test.spearman$estimate, test.spearman$p.value))
				pl[[length(pl)+1]] <- p
				p <- ggplot(df2, aes_string(x=qvar, y=mvar, color=group_var)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("%s vs. %s, %s %s (rho=%.4g, p=%.4g)", qvar, mvar, matrix, assay_type, test.spearman$estimate, test.spearman$p.value)) + scale_color_manual(values=get_color_list(group_var))
				pl[[length(pl)+1]] <- p
			}
			print(multiplot(plotlist=pl, cols=2, rows=2))
		}
	}
}



	
dev.off()
	
	
