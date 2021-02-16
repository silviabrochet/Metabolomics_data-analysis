# Technical reproducibility ------
# The intensity values corresponding to the two injections for each sample are plotted against each other. 

full_data <- read.table("techrep_sheet1.txt",h=T)
full_data_log <- log(full_data)
sample_names <- colnames(full_data)
pdf("techrep_sheet1.pdf")
par(mfrow=c(2,2))
sample_index <- seq(1,by=2,len=122) #122 sample pairs in sheet 1
for (i in sample_index) {
  name <- substring(sample_names[i],2)
  plot(full_data_log[,i],full_data_log[,(i+1)], xlab=("rep1:log intensity"),ylab=("rep2:log intensity"), main=name)
}   
dev.off()

full_data2 <- read.table("techrep_sheet2.txt",h=T)
temp <- cbind(full_data[,245],full_data2)
new_name <- "X245_F5_261_PE_T16_4_F5_261_PE"
temp_names <- colnames(temp)
colnames(temp)[colnames(temp)==temp_names[1]] <- new_name
full_data2 <- temp
full_data2_log <- log(full_data2)
sample_names <- colnames(full_data2)

pdf("techrep_sheet2.pdf")
par(mfrow=c(2,2))
sample_index <- seq(1,by=2,len=87) #87 sample pairs in sheet 2
for (i in sample_index) {
  name <- substring(sample_names[i],2)
  plot(full_data2_log[,i],full_data2_log[,(i+1)], xlab=("rep1:log intensity"),ylab=("rep2:log intensity"), main=name)
}   
dev.off()


# Determine which ions are pollen-derived ------

# Create dataframe with pollen-dilution data = pollen_dilution_data.txt

data_pollen <- read.table("pollen_dilution_data.txt", h=T, check.names = FALSE)

data_pollen$Sample <- factor(data_pollen$Sample, levels = c("1X", "2X", "4X", "8X", "16X", "32X", "64X", "128X", "256X", "512X", "1024X", "2048X"))
data_pollen_levels <- as.numeric(data_pollen$Sample)

# Plot for each ion its intensity values along the dilution serie

pdf("pollen_dilution_serie.pdf")
par(mfrow=c(2,2))
for (i in 2:657) {
  plot(data_pollen_levels, data_pollen[,i], pch=16, main=i-1, xaxt='n',xlab="dilution_serie", ylab="ion_int")
  axis(side=1, cex.axis=0.7, at=1:12, labels=c("1X", "2X", "4X", "8X", "16X", "32X", "64X", "128X", "256X", "512X", "1024X", "2048X"))
  abline(mod <- lm(data_pollen[,i]~data_pollen_levels, data=data_pollen))
  legend("topright", bty="n", legend=paste("R2 is", format(summary(mod)$adj.r.squared, digits=4)))
}
dev.off()


# Retrieve R2 values from each plot

r_sq <- function(dataset) {
  dataset_levels1 <- as.numeric(dataset$Sample)
  r.sq.mod <- vector(mode="numeric", length=0)
  for (i in c(2:length(dataset[1,]))){
    mod <- lm(dataset[,i]~dataset_levels1, data=dataset)
    sum.mod <- summary(mod)
    r.sq.mod[i-1] <- sum.mod$adj.r.squared
  }
  return(r.sq.mod)
}

r_sq <- data.frame(r_sq(data_pollen))

# Calculate log2FC between averages of non-diluted pollen samples (PE1X, 3 replicates) and H2O (average of 6 injections giving 3 replicates) 
# This was done in Excel and then a dataframe coupling the log2FC between PE1X and H2O and the R2 calculated for the dilution curve was created
# log2FC_rsquared.txt

log2FC_rsq <- read.table("log2FC_rsquared.txt", h=T, check.names = FALSE)

# These data were used to produce a Volcano plot and set the threshold for the pollen-derived ions

pdf("volcano_log2FC_rsquared.pdf")
with(log2FC_rsq, plot(log2FC_PE1X_H20, r_sq.data_pollen., pch=20, cex=0.5, main="VP:pollen_vs_H2O", xlim=c(-2,9)))
abline(v=0)
with(subset(log2FC_rsq, log2FC_PE1X_H20<0), points(log2FC_PE1X_H20, r_sq.data_pollen., pch=20, cex=0.5, col="orange"))

abline(h=0.75, col="red")
abline(v=2, col="red")

dev.off()

# Cutoffs to define a ion as pollen-derived = log2FC>=2 and rsq>0.75 
# We have 406 pollen-derived ions out of 657 ions!

log2FC_rsq_subset = subset(log2FC_rsq, log2FC_rsq$log2FC>=2 & log2FC_rsq$r_sq>0.75)
write.table(log2FC_rsq_subset$ion, "pollen_derived_ions.txt")

# Calculate the log2FC (and significance) between samples per each time-point -----
# We the mean_intensities datafile
# T-tests, fold-change and Volcano plots for individual strains
# Will create one file per strain = F5_183_FC, F5_184_FC, F5_185_FC, F5_186_FC
# Combined in all_log2FC.xlsx 

# Read in mean intensities 

data <- read.table("mean_intensities.txt",h=T)

# Lists of samples corresponding to the two time-points measured (not controls)
T0 <-c("F5_183_PE_T0","F5_185_PE_T0","F5_184_PE_T0","F5_186_PE_T0")

T16 <-c("F5_183_PE_T16","F5_185_PE_T16","F5_184_PE_T16","F5_186_PE_T16")

t_test <- function(data, group1, group2) {
  group1_data <- subset(data, data$Group==group1)
  group2_data <- subset(data, data$Group==group2)
  sample_data <- rbind(group1_data,group2_data)
  sample_data <- droplevels(sample_data)
  nb_ions <- length(colnames(sample_data))
  t.test.list <- lapply(sample_data[,c(3:nb_ions)],function(x) t.test(x~sample_data$Group))
  group1_log2 <- log2(group1_data[,c(3:nb_ions)])
  group2_log2 <- log2(group2_data[,c(3:nb_ions)])
  group1_sd <- apply(group1_log2,2,function(x){sd(x)})
  group2_sd <-  apply(group2_log2,2,function(x){sd(x)})
  final_se <- sqrt((group1_sd^2) + (group2_sd^2))
  t.test.list <- list(t.test.list,final_se)
  return(t.test.list)
}

# Loop for generating volcano plots, and t-test data-tables

nb_samples <- length(T0)
index <- c(1:nb_samples)
pdf("volcano_plots.pdf")
par(mfrow=c(2,2))
for (i in index) {
  T0_sample <- T0[i]
  T16_sample <- T16[i]
  strain <- substr(T0_sample,1,6)
  table_name <- paste0(strain,"_FC.txt")
  t.test <- t_test(data, T16_sample, T0_sample)
  t.test.list <- t.test[[1]]
  SE <- t.test[[2]]
  mean.values.list <- lapply(t.test.list, function(x) x$estimate)
  mean.values.T0 <- unlist(lapply(mean.values.list, function(x) x[1]))
  mean.values.T16 <- unlist(lapply(mean.values.list, function(x) x[2]))
  p.values.list <- lapply(t.test.list,function(x) x$p.value)
  p.values <- unlist(p.values.list)
  p.values.adjust <- p.adjust(p.values, method="BH")
  p.values.log <- -log10(p.values.adjust)
  fold_change <- mean.values.T16/mean.values.T0
  FC_log2 <- log2(fold_change)
  ions <- c(1:length(p.values))
  out.table <- data.frame(ions, mean.values.T0, mean.values.T16, fold_change, FC_log2, SE, p.values.adjust)
  names(out.table) <- c("ions",T0_sample,T16_sample,"FC","log2_FC","SE","p.value.adjusted")
  write.table(out.table,file=table_name,row.names=FALSE,quote=FALSE,sep="\t")
  small_index <- which(FC_log2 > -1 & FC_log2 < 1)
  large_index <- which(FC_log2 < -1 | FC_log2 > 1)
  small_x <- FC_log2[small_index]
  small_y <- p.values.log[small_index]
  large_x <- FC_log2[large_index]
  large_y <- p.values.log[large_index]
  plot(large_x,large_y, pch=19, col="blue", main=strain, ylab="-log10 p-value", xlab="log2(Fold-change)", ylim=c(0,10), xlim=c(-6,6), cex=0.7)
  points(small_x, small_y, pch=19, col="grey", cex=0.7)
  abline(h=2,col="red",lty=4)
  grid()
} 
dev.off()


# Identify shared/strainspecific ions (a posteriori) ----
# Start from a table (sharing.txt) where we gathered the ions that change significantly for each strain

sharing <- read.table("sharing.txt", h=T, check.names = FALSE)
vector <- unique(c(sharing$ESL183, sharing$ESL184, sharing$ESL185, sharing$ESL186))
share_res <- c()
i <- 1
for (metabolite in vector) {
  share_res[i] <- sum(c(metabolite %in% sharing$ESL183, metabolite %in% sharing$ESL184, metabolite %in% sharing$ESL185, metabolite %in% sharing$ESL186)) 
  share_res[i]
  i=i+1
}

shared_res <- data.frame("vector"=vector, "share_res"=share_res)

write.table(shared_res, "shared_res.txt")

shared4_subset <- subset(shared_res, shared_res$share_res=="4")
shared4_transpose <- t(shared4_subset$vector)
shared3_subset <- subset(shared_res, shared_res$share_res=="3")
shared3_transpose <- t(shared3_subset$vector)
shared2_subset <- subset(shared_res, shared_res$share_res=="2")
shared2_transpose <- t(shared2_subset$vector)
shared1_subset <- subset(shared_res, shared_res$share_res=="1")
shared1_transpose <- t(shared1_subset$vector)

# Loop to generate volcanos with sharing info

nb_samples <- length(T0)
index <- c(1:nb_samples)
pdf("volcano_plots_shared.pdf")
par(mfrow=c(2,2))
for (i in index) {
  T0_sample <- T0[i]
  T16_sample <- T16[i]
  strain <- substr(T0_sample,1,6)
  table_name <- paste0(strain,"_FC.txt")
  t.test <- t_test(data, T16_sample, T0_sample)
  t.test.list <- t.test[[1]]
  SE <- t.test[[2]]
  mean.values.list <- lapply(t.test.list, function(x) x$estimate)
  mean.values.T0 <- unlist(lapply(mean.values.list, function(x) x[1]))
  mean.values.T16 <- unlist(lapply(mean.values.list, function(x) x[2]))
  p.values.list <- lapply(t.test.list,function(x) x$p.value)
  p.values <- unlist(p.values.list)
  p.values.adjust <- p.adjust(p.values, method="BH")
  p.values.log <- -log10(p.values.adjust)
  fold_change <- mean.values.T16/mean.values.T0
  FC_log2 <- log2(fold_change)
  ions <- c(1:length(p.values))
  out.table <- data.frame(ions, mean.values.T0, mean.values.T16, fold_change, FC_log2, SE, p.values.adjust)
  names(out.table) <- c("ions",T0_sample,T16_sample,"FC","log2_FC","SE","p.value.adjusted")
  write.table(out.table,file=table_name,row.names=FALSE,quote=FALSE,sep="\t")
  small_index <- which(FC_log2 > -1 & FC_log2 < 1)
  large_index <- which(FC_log2 < -1 | FC_log2 > 1)
  shared4 <- out.table[out.table$ions %in% shared4_transpose, ]
  shared3 <- out.table[out.table$ions %in% shared3_transpose, ]
  shared2 <- out.table[out.table$ions %in% shared2_transpose, ]
  shared1 <- out.table[out.table$ions %in% shared1_transpose, ]
  small_x <- FC_log2[small_index]
  small_y <- p.values.log[small_index]
  large_x <- FC_log2[large_index]
  large_y <- p.values.log[large_index]
  plot(large_x,large_y, pch=19, main=strain, ylab="-log10 p-value", xlab="log2(Fold-change)", ylim=c(0,10), xlim=c(-6,6), cex=0.7)
  points(shared4$log2_FC, (-log10(shared4$p.value.adjusted)), pch=19, col="darkblue", cex=0.8)
  points(shared3$log2_FC, (-log10(shared3$p.value.adjusted)), pch=19, col="deepskyblue3", cex=0.8)
  points(shared2$log2_FC, (-log10(shared2$p.value.adjusted)), pch=19, col="deepskyblue1", cex=0.8)
  points(shared1$log2_FC, (-log10(shared1$p.value.adjusted)), pch=19, col="darkslategray1", cex=0.8)
  points(small_x, small_y, pch=19, col=adjustcolor("grey", alpha.f = 0.5), cex=0.7)
  abline(v=c(-1,1),col="red",lty=4)
  abline(h=2,col="red",lty=4)
  grid()
} 

dev.off()

# Produce dotplots for each ion -----
# Create dataset merging files obtained before - selecting only the ions significantly
# changing in at least one strain

dotplot_data <- read.table("dotplot.txt", check.names = FALSE, h=T)

library(gridExtra)
library(ggplot2)

attach(dotplot_data)
dotplot_data$Ions <- as.factor(dotplot_data$Ions)
plot = list()
n=1
a=1
pdf("dotplots.pdf", width=10, height=10)
for (ii in levels(dotplot_data$Ions)) {
  print(ii)
  p=ggplot(dotplot_data[dotplot_data$Ions==ii,], aes(Strain, Values)) + 
    coord_flip() + 
    geom_point(aes(colour = Strain, size = 10, shape = Strain)) + theme_bw() + 
    labs(title = dotplot_data$Top_annotation_name[a]) + 
    geom_hline(yintercept = 0) + theme(legend.position = "none") + ylab(label="log2FC") +xlab(label="") +
    geom_errorbar(aes(ymin = Values - SE, ymax = Values + SE), width=0.3) +
    scale_y_continuous(breaks=c(-7,-6,-5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5,6,7), limits = c(-7,7)) +
    scale_x_discrete(limits = rev(levels(Strain))) +
    scale_color_manual(values = c("Lapi" = "#f2962c", "Lhel" = "#2027e8", "Lmel" = "#fa52de", "Lkul" = "#12e636")) +
    scale_shape_manual(values = c("Lapi" = 15, "Lhel" = 17, "Lmel" = 18, "Lkul" = 16)) +
    theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"), plot.title = element_text(size=20))
  plot[[n]]=p
  n=n+1
  a=a+4
  print(p)
}
dev.off()


# PCA ------

PCA2 <- read.table("PCA_R.txt", h=T, check.names = FALSE)
attach(PCA2)
library(vegan)
library(ggfortify)
library(ggplot2)

rnames <- PCA2[,1]
PCA2$Strains <- NULL
rownames(PCA2) = make.names(rnames, unique=TRUE)
PCA_data1 <- rda(PCA2[,2:657])

pdf("PCA.pdf")

par(mfrow=c(1,1)) 
shrink <- FALSE 
cols <- rep(c("orange", "green", "blue", "magenta"), each = 5) 
spp.sc <- scores(PCA_data1, display = "species", shrink = shrink) 
site.sc <- scores(PCA_data1, display = "sites") 
ylim <- range(spp.sc[,2], site.sc[,2]) 
xlim <- range(spp.sc[,2], site.sc[,2]) 
# Obtain axis labels by looking at PC1 and PC2 values in PCA_data1 
ordiplot(PCA_data1,type = "n", xlab="PC1 = 34.31%", ylab="PC2 = 16.40")
points(site.sc, col = cols, labels = rownames(PCA2), cex = 2, pch = c(15, 20, 17, 18)[as.numeric(PCA2$Group)]) 
legend("bottomright", legend = c("Lb_185","Lb_183","Lb_184", "Lb_186"), pch = c(15, 17, 18, 20), 
       col = c("orange", "blue", "magenta", "green"), cex=0.8)
dev.off()

PCA_data2 <- prcomp(PCA2[,2:657])
cols <- rep(c("orange", "green", "blue", "magenta"), each = 5) 
theme_set(theme_bw())
autoplot(PCA_data2, colour =cols, loadings = TRUE, 
         loadings.colour = 'white',
         loadings.label = TRUE, loadings.label.size = 3) 

