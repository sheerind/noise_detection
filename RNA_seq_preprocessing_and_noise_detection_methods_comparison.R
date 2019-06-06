library(limma)
library(edgeR)
library(biomaRt)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rowr)
library(scales)
library(UpSetR)
library(plyr)
library(cowplot)
library(DESeq2)
library(geneplotter)
library(Hmisc)
library(gplots)
library(grid)
library(gridExtra)

#read in sample info and noise thresholds
sample_info <- read.csv("/home/dylan/Documents/pipeline_comparison/EUCLIDS/sample_info.csv", header = T, row.names = 1)
min_noise <- 2^(sample_info$min_noise)

#load stringtie data
hisat_string_raw <- read.csv("/home/dylan/Documents/pipeline_comparison/EUCLIDS/Hisat2/stringtie/gene_count_matrix.csv", header=T, row.names=1)
colnames(hisat_string_raw) <- gsub("X", "", colnames(hisat_string_raw))
colnames(hisat_string_raw) <- gsub(colnames(hisat_string_raw), pattern = "\\.", replacement = "-")
hisat_string_raw <- hisat_string_raw[-c(17,18)] #sex discrepancy
dim(hisat_string_raw)
#58395 60
hi_str_exp <- hisat_string_raw[rowSums(hisat_string_raw) != 0,]
dim(hi_str_exp)
#54482 60
hi_str_exp <- rbind(hi_str_exp,min_noise)

star_string_raw <- read.csv("/home/dylan/Documents/pipeline_comparison/EUCLIDS/STAR/stringtie/gene_count_matrix.csv", header=T, row.names=1)
colnames(star_string_raw) <- gsub("X", "", colnames(star_string_raw))
colnames(star_string_raw) <- gsub(colnames(star_string_raw), pattern = "\\.", replacement = "-")
star_string_raw <- star_string_raw[-c(17,18)] #sex discrepancy
dim(star_string_raw)
#58395 60
star_str_exp <- star_string_raw[rowSums(star_string_raw) != 0,]
dim(star_str_exp)
#54074 60
star_str_exp <- rbind(star_str_exp,min_noise)

#load htseq data
hisat_htseq <- list.files(path="/home/dylan/Documents/pipeline_comparison/EUCLIDS/Hisat2/htseq", pattern="*-1.txt", full.names=T)
hisat_htseq_files <- lapply(hisat_htseq, read.table)
hisat_htseq_raw <- as.data.frame(sapply(hisat_htseq_files, function(x) x[,2]))
hisat_htseq <- gsub("-1.txt", "", hisat_htseq)
hisat_htseq <- gsub("/home/dylan/Documents/pipeline_comparison/EUCLIDS/Hisat2/htseq/result_", "", hisat_htseq)
colnames(hisat_htseq_raw) <- hisat_htseq
row.names(hisat_htseq_raw) <- hisat_htseq_files[[1]]$V1
hisat_htseq_raw <- hisat_htseq_raw[,-c(17,18)] #sex discrepancy
hisat_htseq_raw <- hisat_htseq_raw[-nrow(hisat_htseq_raw),]
hisat_htseq_raw <- hisat_htseq_raw[-nrow(hisat_htseq_raw),]
hisat_htseq_raw <- hisat_htseq_raw[-nrow(hisat_htseq_raw),]
hisat_htseq_raw <- hisat_htseq_raw[-nrow(hisat_htseq_raw),]
hisat_htseq_raw <- hisat_htseq_raw[-nrow(hisat_htseq_raw),]
dim(hisat_htseq_raw)
#58395 60
hi_ht_exp <- hisat_htseq_raw[rowSums(hisat_htseq_raw) != 0,]
dim(hi_ht_exp)
#53862 60
hi_ht_exp <- rbind(hi_ht_exp,min_noise)

star_htseq <- list.files(path="/home/dylan/Documents/pipeline_comparison/EUCLIDS/STAR/htseq", pattern="*.txt", full.names=T)
star_htseq_files <- lapply(star_htseq, read.table)
star_htseq_raw <- as.data.frame(sapply(star_htseq_files, function(x) x[,2]))
star_htseq <- gsub(".txt", "", star_htseq)
star_htseq <- gsub("/home/dylan/Documents/pipeline_comparison/EUCLIDS/STAR/htseq/result_", "", star_htseq)
colnames(star_htseq_raw) <- star_htseq
row.names(star_htseq_raw) <- star_htseq_files[[1]]$V1
star_htseq_raw <- star_htseq_raw[,-c(17,18)] #sex discrepancy
star_htseq_raw <- star_htseq_raw[-nrow(star_htseq_raw),]
star_htseq_raw <- star_htseq_raw[-nrow(star_htseq_raw),]
star_htseq_raw <- star_htseq_raw[-nrow(star_htseq_raw),]
star_htseq_raw <- star_htseq_raw[-nrow(star_htseq_raw),]
star_htseq_raw <- star_htseq_raw[-nrow(star_htseq_raw),]
dim(star_htseq_raw)
#58395 60
star_ht_exp <- star_htseq_raw[rowSums(star_htseq_raw) != 0,]
dim(star_ht_exp)
#53855 60
star_ht_exp <- rbind(star_ht_exp,min_noise)

#load featurecounts data
hisat_feature_raw <- read.table("/home/dylan/Documents/pipeline_comparison/EUCLIDS/Hisat2/featurecounts/featurecounts_hisat2.txt", header=T, row.names=1)
hisat_feature_raw <- hisat_feature_raw[,-c(1:5)]
colnames(hisat_feature_raw) <- gsub("X", "", colnames(hisat_feature_raw))
colnames(hisat_feature_raw) <- gsub(".1.bam", "", colnames(hisat_feature_raw))
colnames(hisat_feature_raw) <- gsub(colnames(hisat_feature_raw), pattern = "\\.", replacement = "-")
hisat_feature_raw <- hisat_feature_raw[c(9,14,19,24,29,34,39,44,49,54,59,1,5,10,15,20,25,30,35,40,45,50,55,60,2,6,11,16,21,26,31,36,41,46,51,56,61,3,7,12,17,22,32,37,42,47,52,57,62,4,8,13,18,23,28,33,38,43,48,53,58,63)]
hisat_feature_raw <- hisat_feature_raw[,-c(17,18)]#sex discrepancy
dim(hisat_feature_raw)
#58395 60
hi_fc_exp <- hisat_feature_raw[rowSums(hisat_feature_raw) != 0,]
dim(hi_fc_exp)
#54281 60
hi_fc_exp <- rbind(hi_fc_exp,min_noise)

star_feature_raw <- read.table("/home/dylan/Documents/pipeline_comparison/EUCLIDS/STAR/featurecounts/featurecounts_STAR.txt", header=T, row.names=1)
star_feature_raw <- star_feature_raw[,-c(1:5)]
colnames(star_feature_raw) <- gsub("X", "", colnames(star_feature_raw))
colnames(star_feature_raw) <- gsub(".bam", "", colnames(star_feature_raw))
colnames(star_feature_raw) <- gsub(colnames(star_feature_raw), pattern = "\\.", replacement = "-")
star_feature_raw <- star_feature_raw[c(9,14,19,24,28,33,38,43,48,53,58,1,5,10,15,20,25,29,34,39,44,49,54,59,2,6,11,16,21,26,30,35,40,45,50,55,60,3,7,12,17,22,31,36,41,46,51,56,61,4,8,13,18,23,27,32,37,42,47,52,57,62)]
star_feature_raw <- star_feature_raw[,-c(17,18)]#sex discrepancy
dim(star_feature_raw)
#58395 60
star_fc_exp <- star_feature_raw[rowSums(star_feature_raw) != 0,]
dim(star_fc_exp)
#54086 60
star_fc_exp <- rbind(star_fc_exp,min_noise)

#load STAR quantMode data
reads_per_gene <- list.files(path="/home/dylan/Documents/pipeline_comparison/EUCLIDS/STAR/reads_per_gene", pattern="*ReadsPerGene.out.tab$", full.names=T)
reads_per_gene_files <- lapply(reads_per_gene, read.table, skip = 4 )
reads_per_gene_raw <- as.data.frame(sapply(reads_per_gene_files, function(x) x[,4]))
reads_per_gene <- gsub("ReadsPerGene[.]out[.]tab", "", reads_per_gene)
reads_per_gene <- gsub("/home/dylan/Documents/pipeline_comparison/EUCLIDS/STAR/reads_per_gene/", "", reads_per_gene)
colnames(reads_per_gene_raw) <- reads_per_gene
row.names(reads_per_gene_raw) <- reads_per_gene_files[[1]]$V1
reads_per_gene_raw <- reads_per_gene_raw[,-c(17,18,43,44)]
dim(reads_per_gene_raw)
#58299 60
star_star_exp <- reads_per_gene_raw[rowSums(reads_per_gene_raw) != 0,]
dim(star_star_exp)
#51586 60
star_star_exp <- rbind(star_star_exp,min_noise)

#load kallisto data
ensdb <- EnsDb.Hsapiens.v86
txk <- keys(ensdb, keytype = "TXNAME")
tx2gene <- select(ensdb, txk, "GENEID", "TXNAME")
tx2gene <- tx2gene[,-3]
kallisto <- list.files(path="/home/dylan/Documents/pipeline_comparison/EUCLIDS/kallisto/counts", pattern="*.tsv", full.names=T)
txi <- tximport(kallisto, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
kallisto_counts <- txi$counts
kallisto <- gsub(".tsv", "", kallisto)
kallisto <- gsub("/home/dylan/Documents/pipeline_comparison/EUCLIDS/kallisto/counts/", "", kallisto)
colnames(kallisto_counts) <- kallisto
kallisto_counts <- kallisto_counts[,-c(17,18)] #sex discrepancy
dim(kallisto_counts)
#39429 60
kall_exp <- kallisto_counts[rowSums(kallisto_counts) != 0,]
kall_exp <- as.data.frame(kall_exp)
kall_exp <- round(kall_exp, digits = 0)
dim(kall_exp)
#37907 60
kall_exp <- rbind(kall_exp,min_noise)

#find intersection of expressed genes
#int_all_nz <- Reduce(intersect,list(rownames(x),rownames(x2),rownames(x3)...)
common_x <- x[int_all_nz,]

#determine number of unique genes
unique_x <- x[-which(rownames(x) %in% rownames(x)),]
dim(unique_x)
#nz_unique <- rowr::cbind.fill(unique_x[,1],unique_x2[,1],unique_x3[,1]..., fill = NA)
l_unique <- log2(nz_unique+1)
gg_unique_dens <- reshape2::melt(as.matrix(l_unique))
gg_uni_dens <- ggplot(gg_unique_dens, aes(x = Expression, colour = Sample)) +
  geom_density(alpha = 0.2, size = 1.5) +
  scale_x_continuous(limits = c(-5,15)) +
  scale_y_continuous(limits = c(0,0.6)) +
  labs(title = "Unique genes from each method combination", x = (expression(log[2](count+1))), y = "Density") +
  theme(plot.title = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  geom_vline(xintercept = median(gg_unique_dens$Expression,na.rm = T), colour = "blue") +
  theme_gray()
dev.off()

#MA plots of genes in the intersection
#nz_intersect <- cbind(common_x[,1],common_x2[,1],common_x3[,1]...)
for (i in 1:ncol(nz_intersect)){
  for (j in 1:ncol(nz_intersect)){
    m <- log2((nz_intersect[,i]+1)/(nz_intersect[,j]+1))
    a <- log2(((nz_intersect[,i]+1) + (nz_intersect[,j]+1))/2)
    a.bins <- as.factor(floor(a*10))
    boxplot(m~a.bins, main = paste(colnames(nz_intersect)[i], " vs. ", colnames(nz_intersect)[j]),  xlab = "(A) mean average", ylab = "(M) log ratio", xlim = c(0,160), ylim = c(-20,20), xaxt='n')
    abline(h=0,col="red");abline(h=-0.5,col="blue");abline(h=0.5,col="blue")
    axis(side = 1)
  }
}
dev.off()

#total genes
#genes <- rowr::cbind.fill(rownames(x),rownames(x2),rownames(x3)...fill = NA)
#select genes in the intersection for two samples
int_genes_s1 <- rowr::cbind.fill(common_x$sample1,common_x2$sample1,common_x3$sample1, fill = NA)
int_genes_s2 <- rowr::cbind.fill(common_x$sample2,common_x2$sample2,common_x3$sample2, fill = NA)
#select 100 most variable genes in the intersection for two samples
log_int_genes_s1 <- cpm(int_genes_s1, log = T)
log_int_genes_s2 <- cpm(int_genes_s2, log = T)
log_int_genes_s1 <- as.matrix(log_int_genes_s1)
log_int_genes_s2 <- as.matrix(log_int_genes_s2)
var_int_genes_s1 <- apply(log_int_genes_s1, 1, var)
var_int_genes_s2 <- apply(log_int_genes_s2, 1, var)
select_var_int_genes_s1 <- names(sort(var_int_genes_s1, decreasing=TRUE))[1:100]
select_var_int_genes_s2 <- names(sort(var_int_genes_s2, decreasing=TRUE))[1:100]
high_var_int_genes_s1 <- log_int_genes_s1[select_var_int_genes_s1,]
high_var_int_genes_s2 <- log_int_genes_s2[select_var_int_genes_s2,]
#plot heatmap
par(cex.main=1)
heatmap.2(high_var_int_genes_s1, main = "100 most variable genes common to all method combinations (sample 1)", srtCol = 45, cexCol = 1, margins = c(8, 4), labRow = F, dendrogram = "column", keysize=1, key.par = list(cex=0.5))
dev.off()
par(cex.main=1)
heatmap.2(high_var_int_genes_s2, main = "100 most variable genes common all method combinations (sample 2)", srtCol = 45, cexCol = 1, margins = c(8, 4), labRow = F, dendrogram = "column", keysize=1, key.par = list(cex=0.5))
dev.off()

#create DGEList object
group <- factor(paste(sample_info$group,sample_info$time_point,sep="_"))
dge_x <- DGEList(x, group = group)

#pre-normalisation density plots
cpm_dens_plot <- function(x){
  cpm_x <- edgeR::cpm(x, log = T)
  x_dens <- reshape2::melt(as.matrix(cpm_x))
  colnames(x_dens) <- c("Gene","Sample","Expression")
  ggplot(x_dens, aes(x = Expression, colour = Sample)) +
    geom_density(alpha = 0.2, size = 1.5) +
    scale_x_continuous(limits = c(-10,20)) +
    scale_y_continuous(limits = c(0,0.6)) +
    geom_vline(xintercept = median(x_dens$Expression,na.rm = T), colour = "blue") +
    theme_gray()
}

#crude filter for low abundance genes
keep <- function(x){
  keep <- rowSums(edgeR::cpm(x)>1)>=15
  x[keep,, keep.lib.sizes=F]
}
dge_x_filt <- keep(dge_x)

#calculate normalisation factors for limma/voom
dge_x <- calcNormFactors(x, method = "TMM")
dge_x_filt <- calcNormFactors(hi_str_filt, method = "TMM")

#create design matrix
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))
contr.matrix <- makeContrasts(controls = Control_D1 - Control_D0, tests = Test_D1 - Test_D0, test_v_control = Test_D1 - Control_D1, levels=colnames(design))

#voom transform
v_x <- voom(dge_x, design, plot = T)
v_x_filt <- voom(dge_x, design, plot = T)

#post-voom density plots
post_dens_plot <- function(x){
  x_dens <- reshape2::melt(as.matrix(x$E))
  colnames(x_dens) <- c("Gene","Sample","Expression")
  ggplot(x_dens, aes(x = Expression, colour = Sample)) +
    geom_density(alpha = 0.2, size = 1.5) +
    scale_x_continuous(limits = c(-10,20)) +
    scale_y_continuous(limits = c(0,0.6)) +
    geom_vline(xintercept = median(x_dens$Expression,na.rm = T), colour = "blue") +
    theme_gray()
}

#apply noise thresholds to voom transformed data
v_nc_x <- v_x
v_nc_x_filt <- v_x_filt
for(i in 1:ncol(v_nc_x$E)){
  v_nc_x$E[,i] <- ifelse(v_nc_x$E[,i] < v_nc_x$E[nrow(v_nc_x),i], v_nc_x$E[nrow(v_nc_x),i], v_nc_x$E[,i])
}
for(i in 1:ncol(v_nc_x_filt$E)){
  v_nc_x_filt$E[,i] <- ifelse(v_nc_x_filt$E[,i] < v_nc_x_filt$E[nrow(v_nc_x_filt),i], v_nc_x_filt$E[nrow(v_nc_x_filt),i], v_nc_x_filt$E[,i])
}
  
#remove rows corresponding to sum of noise threshold
rownames(v_nc_x$weights) <- rownames(v_nc_x$E)
v_nc_x$weights <- v_nc_x$weights[rowSums(v_nc_x$E) != sum(v_nc_x$E[nrow(v_nc_x),]),]
v_nc_x$E <- v_nc_x$E[rowSums(v_nc_x$E) != sum(v_nc_x$E[nrow(v_nc_x),]),]
rownames(v_nc_x_filt$weights) <- rownames(v_nc_x_filt$E)
v_nc_x_filt$weights <- v_nc_x_filt$weights[rowSums(v_nc_x_filt$E) != sum(v_nc_x_filt$E[nrow(v_nc_x_filt),]),]
v_nc_x_filt$E <- v_nc_x_filt$E[rowSums(v_nc_x_filt$E) != sum(v_nc_x_filt$E[nrow(v_nc_x_filt),]),]

#fit linear models and contrasts matrix
fit_x <- lmFit(v_x, design)
fit_x <- contrasts.fit(fit_x, contrasts = contr.matrix)
fit_nc_x <- lmFit(v_nc_x, design)
fit_nc_x <- contrasts.fit(fit_nc_x, contrasts = contr.matrix)
fit_x_filt <- lmFit(v_x_filt, design)
fit_x_filt <- contrasts.fit(fit_x_filt, contrasts = contr.matrix)
fit_nc_x_filt <- lmFit(v_nc_x_filt, design)
fit_nc_x_filt <- contrasts.fit(fit_nc_x_filt, contrasts = contr.matrix)

#empirical Bayes
efit_x <- eBayes(fit_x)
efit_nc_x <- eBayes(fit_nc_x)
efit_x_filt <- eBayes(fit_x_filt)
efit_nc_x_filt <- eBayes(fit_nc_x_filt)

#create results tables
tt_contr <- function(x){
  tt <- topTable(x, 
                 coef="controls",
                 number=Inf, adjust.method="BH",
                 sort.by="logFC")
  tt[-1,]
}
tt_test <- function(x){
  tt <- topTable(x, 
                 coef="tests",
                 number=Inf, adjust.method="BH",
                 sort.by="logFC")
  tt[-1,]
}
padj.cutoff <- 0.01
lfc.cutoff <- 1
res_x_contr <- tt_contr(efit_x)
res_x_contr$threshold <- res_x_contr$adj.P.Val < padj.cutoff & abs(res_x_contr$logFC) > lfc.cutoff
res_x_test <- tt_test(efit_x)
res_x_test$threshold <- res_x_test$adj.P.Val < padj.cutoff & abs(res_x_test$logFC) > lfc.cutoff
#etc...

#MA plots
ma_plot <- function(x){
  ggplot(data=x, aes(x=AveExpr, y=logFC, colour=threshold)) + 
    geom_point(alpha=0.4, size=0.5) + 
    geom_hline(aes(yintercept = 0), colour = "blue", size = 0.5) +
    geom_hline(aes(yintercept = 1), colour = "red", size = 0.5) +
    geom_hline(aes(yintercept = -1), colour = "red", size = 0.5) +
    ylim(c(-6,6)) + 
    xlim(c(-5,15)) +
    theme(axis.title.x = element_text(face = "bold", size = 15),
          axis.text.x = element_text(face = "bold", size = 12)) +
    theme(axis.title.y = element_text(face = "bold", size = 15),
          axis.text.y = element_text(face = "bold", size = 12)) +
    scale_colour_manual(labels = c("Not significant", "Significant"), values = c("lightblue2","red")) +
    theme(legend.text = element_text(size = 12)) +
    theme_gray()
}

