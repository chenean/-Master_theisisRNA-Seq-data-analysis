#limma-voom differential expression analysis on normalized continuous reads.

setwd("~/Documents/uni/MISB2_semester4/master thesis intership/computional part/RNA_seq training/GSE6024_dataset_normalizedread")

##find confounding factors: e.g. sex, age, etenic,..
library("GEOquery") #version 3.8
gsm<-getGEO(filename = "GSE60424_series_matrix.txt.gz")
table<-gsm@phenoData@data # extract the data.frame


library(ggplot2)
library(ggfortify)
##select whole blood samples
#consider only data set with at least 4 CTRL and 4 Disease for further analysis
index_blood<-which(as.character(table$`celltype:ch1`)=="Whole Blood")
table_whole_blood <- table[index_blood,]

index1 <- which(table_whole_blood$`diseasestatus:ch1`=="Healthy Control")
index2 <- which(table_whole_blood$`diseasestatus:ch1`=="Type 1 Diabetes")
index3 <- c(index1, index2)
table_wb_diabetes <- table_whole_blood[index3,]
#the other diseases do not have the same dim() vs Ctrl
# t.test

#age
index_HC<-which(table_wb_diabetes$`diseasestatus:ch1`=="Healthy Control")
index_Diab<-which(table_wb_diabetes$`diseasestatus:ch1`=="Type 1 Diabetes")
age_HC<-as.numeric(table_wb_diabetes$`age:ch1`[index_HC])
age_Diab<-as.numeric(table_wb_diabetes$`age:ch1`[index_Diab])
t.test(age_HC,age_Diab)
#p-value = 0.9635

#gender
#For the categorical test, only 2 groups have at lest 4 samples.
g_HC <- as.character(table_wb_diabetes$`gender:ch1`["index_HC"])
g_Diab <- as.character(table_wb_diabetes$`gender:ch1`["index_Diab"])
g_matrix <- t(rbind(g_HC, g_Diab))
HC_Male <- length(subset(g_matrix[,1] ,g_HC == "Male"))
HC_Female <- length(subset(g_matrix[,1] ,g_HC == "Female"))
Diab_Male <- length(subset(g_matrix[,1] ,g_Diab == "Male"))
Diab_Female <- length(subset(g_matrix[,1] ,g_Diab == "Female"))
gender_row1 <- c(HC_Male, HC_Female)
gender_row2 <- c(Diab_Male,Diab_Female)
gender_matrix <- rbind(gender_row1, gender_row2)
colnames(gender_matrix ) <-c("Male", "Female") 
rownames(gender_matrix) <- c("HC", "Diab")
#If you just check the gender matrix, you see that there are no males in the HC.
#In this case no need for fisher test!
fisher.test(gender_matrix)
#p-value = 1

#ethenic
r_HC <- as.character(table_wb_diabetes$`race:ch1`["index_HC"])
r_Diab <- as.character(table_wb_diabetes$`race:ch1`["Tindex_Diab"])
r_matrix <- t(rbind(r_HC, r_Diab))
HC_Hisp <- length(subset(r_matrix[,1] ,r_HC == "Hispanic"))
HC_White <- length(subset(r_matrix[,1] ,r_HC == "White"))
Diab_Hisp <- length(subset(r_matrix[,1] ,r_Diab == "Hispanic"))
Diab_White <- length(subset(r_matrix[,1] ,r_Diab == "White"))
race_row1 <- c(HC_Hisp, HC_White)
race_row2 <- c(Diab_Hisp,Diab_White)
race_matrix <- rbind(race_row1, race_row2)
colnames(race_matrix ) <-c("Hispanic", "White") 
rownames(race_matrix) <- c("HC", "Diab")
#If you check the race matrix, you see that there are no Hispanic with diabetes.
fisher.test(race_matrix)
#p-value = 1


##select neutrophils
index_neutrophils<-which(as.character(table$`celltype:ch1`)=="Neutrophils")
table_neutrophils <- table[index_neutrophils,]
index1n <- which(table_neutrophils$`diseasestatus:ch1`=="Healthy Control")
index2n <- which(table_neutrophils$`diseasestatus:ch1`=="Type 1 Diabetes")
index3n <- c(index1n, index2n)
table_neutro_diabetes <- table_neutrophils[index3n,]
# t.test (the age is the same cause cell types from the same patients CTRL vs DISEASE)
#index_HC_N<-which(table_neutro_diabetes$`diseasestatus:ch1`=="Healthy Control")
#index_Diab_N<-which(table_neutro_diabetes$`diseasestatus:ch1`=="Type 1 Diabetes")
#age_HC_N<-as.numeric(table_neutro_diabetes$`age:ch1`[index_HC_N])
#age_Diab_N<-as.numeric(table_neutro_diabetes$`age:ch1`[index_Diab_N])
#t.test(age_HC_N,age_Diab_N)
 

##select monocytes
index_mono<-which(as.character(table$`celltype:ch1`)=="Monocytes")
table_mono <- table[index_mono,]
index1m <- which(table_mono$`diseasestatus:ch1`=="Healthy Control")
index2m <- which(table_mono$`diseasestatus:ch1`=="Type 1 Diabetes")
index3m <- c(index1m, index2m)
table_mono_diabetes <- table_mono[index3m,]

##select B-cells
index_B<-which(as.character(table$`celltype:ch1`)=="B-cells")
table_B <- table[index_B,]
index1B <- which(table_B$`diseasestatus:ch1`=="Healthy Control")
index2B <- which(table_B$`diseasestatus:ch1`=="Type 1 Diabetes")
index3B <- c(index1B, index2B)
table_B_diabetes <- table_B[index3B,]

##select CD4
index_CD4<-which(as.character(table$`celltype:ch1`)=="CD4")
table_CD4 <- table[index_CD4,]
index1CD4 <- which(table_CD4$`diseasestatus:ch1`=="Healthy Control")
index2CD4 <- which(table_CD4$`diseasestatus:ch1`=="Type 1 Diabetes")
index3CD4 <- c(index1CD4, index2CD4)
table_CD4_diabetes <- table_CD4[index3CD4,]


##select CD8
index_CD8<-which(as.character(table$`celltype:ch1`)=="CD8")
table_CD8 <- table[index_CD8,]
index1CD8 <- which(table_CD8$`diseasestatus:ch1`=="Healthy Control")
index2CD8 <- which(table_CD8$`diseasestatus:ch1`=="Type 1 Diabetes")
index3CD8 <- c(index1CD8, index2CD8)
table_CD8_diabetes <- table_CD8[index3CD8,]


#select NK (less than 4 samples)
#index_NK<-which(as.character(table$`celltype:ch1`)=="NK")
#table_NK <- table[index_NK,]
#index1NK <- which(table_whole_blood$`diseasestatus:ch1`=="Healthy Control")
#index2NK <- which(table_whole_blood$`diseasestatus:ch1`=="Type 1 Diabetes")
#index3NK <- c(index1NK, index2NK)
#table_NK_diabetes <- table_neutrophils[index3NK,]


library(dplyr)


#voom differential analysis
#limma approach is based on "dge" which is based on raw_counts
Norma_countData <- read.table("GSE60424_GEOSubmit_FC1to11_normalized_counts.txt", header = TRUE, sep="\t")
Norma_data <- data.matrix(Norma_countData[,c(2:ncol(Norma_countData))])
rownames(Norma_data) <- Norma_countData$genenames
colnames(Norma_data) <- table$title


#Select the samples from whole blood
wb_samples_diabetes<-table_wb_diabetes$title
colindex<-vector(mode="integer",length=length(wb_samples_diabetes))
for (i in 1:length(wb_samples_diabetes)){
  colindex[i]<-which(colnames(Norma_data)==wb_samples_diabetes[i])
}
Norma_data_wb_diabetes<-Norma_data[,colindex]


#remove rows where colums=0 (for all rows at least one is bigger than 0 -> filter out all rows with all 0)
final_Norma_data <- Norma_data_wb_diabetes[rowSums(Norma_data_wb_diabetes > 0) >= 1,]


design<-model.matrix(~0+table_wb_diabetes$`diseasestatus:ch1`)
colnames(design) <- c("Ctrl", "Diabetes")
#design the model to be fitted. before using voom because voom uses variances of the model residuals (observed-fitted)
#age = as.factor(table$`age:ch1`) as numeric


library(limma)
# command : BiocManager::install("limma", version = "3.8")
v <- voom(final_Norma_data, design, plot=TRUE )
table.v <- v$E #the log values from here for heatmap (without row Z score, st) and PCA (with row z score: Zscore gives the correlation matrix)
colnames(table.v) <- c("Ctrl","Ctrl","Ctrl","Ctrl","Diabetes","Diabetes","Diabetes","Diabetes")

table.v_id<- cbind(rownames(table.v), table.v)

#####################
# Top 10 gene ontology terms between 5 & 500 according to q-values
#GO1 ribosome biogenesis
index <- integer(270)
for (i in 1:270)
  index[i] <- which(table.v_id[,1] == gene_ribosome[i])

ribosome_table<- table.v[index,]
#GO1 average
GO1 <- t(colMeans(ribosome_table))
rownames(GO1) <- c("GO1")

#GO2 ribonucleoprotein complex 
index <- integer(425)
for (i in 1:425)
  index[i] <- which(table.v_id[,1] == gene_ribonucleoprotein[i])

ribonucleoprotein_table<- table.v[index,]
#GO2 average
GO2 <- t(colMeans(ribonucleoprotein_table))
rownames(GO2) <- c("GO2")

#GO3 mitochondrial gene expression 
index <- integer(154)
for (i in 1:154)
  index[i] <- which(table.v_id[,1] == gene_mitochondrial[i])

mitochondrial_table<- table.v[index,]
#GO3 average
GO3 <- t(colMeans(mitochondrial_table))
rownames(GO3) <- c("GO3")


#GO4 translational initiation
index <- integer(184)
for (i in 1:184)
  index[i] <- which(table.v_id[,1] == gene_translational[i])

translational_table<- table.v[index,]
#GO4 average
GO4 <- t(colMeans(translational_table))
rownames(GO4) <- c("GO4")

#GO5 protein localization to membran
index <- integer(301)
for (i in 1:301)
  index[i] <- which(table.v_id[,1] == gene_proteinlocalization[i])

proteinlocalization_table<- table.v[index,]
#GO5 average
GO5 <- t(colMeans(proteinlocalization_table))
rownames(GO5) <- c("GO5")

#GO6 oxidatice phosphorylation
index <- integer(119)
for (i in 1:119)
  index[i] <- which(table.v_id[,1] == gene_oxidative[i])

oxidative_table<- table.v[index,]
#GO6 average
GO6 <- t(colMeans(oxidative_table))
rownames(GO6) <- c("GO6")


#GO7 nucleoside triphosphate
index <- integer(276)
for (i in 1:276)
  index[i] <- which(table.v_id[,1] == gene_nucleosidetri[i])

nucleosidetri_table<- table.v[index,]
#GO7 average
GO7 <- t(colMeans(nucleosidetri_table))
rownames(GO7) <- c("GO7")


#GO8 nucleoside monophosphate 
index <- integer(288)
for (i in 1:288)
  index[i] <- which(table.v_id[,1] == gene_nucleosidemono[i])

nucleosidemono_table<- table.v[index,]
#GO8 average
GO8 <- t(colMeans(nucleosidemono_table))
rownames(GO8) <- c("GO8")


#GO9 oxidatice phosphorylation
index <- integer(69)
for (i in 1:69)
  index[i] <- which(table.v_id[,1] == gene_ribosomallarge[i])

ribosomallarge_table<- table.v[index,]
#GO9 average
GO9 <- t(colMeans(ribosomallarge_table))
rownames(GO9) <- c("GO9")


#GO10 ATP metabolic process 
index <- integer(225)
for (i in 1:225)
  index[i] <- which(table.v_id[,1] == gene_ATP[i])

ATP_table<- table.v[index,]
#GO10 average
GO10 <- t(colMeans(ATP_table))
rownames(GO10) <- c("GO10")

GO_10<- rbind(GO1, GO2, GO3, GO4,GO5, GO6, GO7, GO8,GO9, GO10)

GO_10_mean <- apply(GO_10, 1, mean)
GO_10_sd <- apply(GO_10, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
GO_10_zscore = (GO_10 - GO_10_mean)/GO_10_sd
Heatmap(GO_10_zscore, name = "Row Z-score", column_title = "Whole blood", cluster_columns = FALSE)

#######################################
# Top 10 enrichment pathway between 5 & 500 according to q-values
#GO1 Translation
index <- integer(291)
for (i in 1:291)
  index[i] <- which(table.v_id[,1] == gene_Translation[i])

Translation_table<- table.v[index,]
#GO1 average
P1 <- t(colMeans(Translation_table))
rownames(P1) <- c("P1")

#GO2 Ribosome-Homo sapiens (human)
index <- integer(127)
for (i in 1:127)
  index[i] <- which(table.v_id[,1] == gene_Ribosomehuman[i])

Ribosomehuman_table<- table.v[index,]
#GO2 average
P2 <- t(colMeans(Ribosomehuman_table))
rownames(P2) <- c("P2")

#GO3 Cap-dependent Translation initiation 
index <- integer(120)
for (i in 1:120)
  index[i] <- which(table.v_id[,1] == gene_Cap[i])

Cap_table<- table.v[index,]
#GO3 average
P3 <- t(colMeans(Cap_table))
rownames(P3) <- c("P3")


#GO4 Eukaryotic Translation Initiation
index <- integer(120)
for (i in 1:120)
  index[i] <- which(table.v_id[,1] == gene_Initiation[i])

Initiation_table<- table.v[index,]
#GO4 average
P4 <- t(colMeans(Initiation_table))
rownames(P4) <- c("P4")

#GO5 GTP hydrolysis and joining of the 60S ribosomal subunit
index <- integer(113)
for (i in 1:113)
  index[i] <- which(table.v_id[,1] == gene_GTP[i])

GTP_table<- table.v[index,]
#GO5 average
P5 <- t(colMeans(GTP_table))
rownames(P5) <- c("P5")

#GO6 L13a-mediated translational silencing of Ceruloplasmin expression
index <- integer(112)
for (i in 1:112)
  index[i] <- which(table.v_id[,1] == gene_L13a[i])

L13a_table<- table.v[index,]
#GO6 average
P6 <- t(colMeans(L13a_table))
rownames(P6) <- c("P6")


#GO7 SRP-dependent cotranslational protein targeting to membrane
index <- integer(113)
for (i in 1:113)
  index[i] <- which(table.v_id[,1] == gene_SRP[i])

SRP_table<- table.v[index,]
#GO7 average
P7 <- t(colMeans(SRP_table))
rownames(P7) <- c("P7")


#GO8 Formation of a pool of free 40S subunits
index <- integer(102)
for (i in 1:102)
  index[i] <- which(table.v_id[,1] == gene_Formation[i])

Formation_table<- table.v[index,]
#GO8 average
P8 <- t(colMeans(Formation_table))
rownames(P8) <- c("P8")


#GO9 Eukaryotic Translation Elongation
index <- integer(96)
for (i in 1:96)
  index[i] <- which(table.v_id[,1] == gene_Eukaryotic[i])

Eukaryotic_table<- table.v[index,]
#GO9 average
P9 <- t(colMeans(Eukaryotic_table))
rownames(P9) <- c("P9")


#GO10 Oxidative phosphorylation- Homo sapiens(human)
index <- integer(119)
for (i in 1:119)
  index[i] <- which(table.v_id[,1] == gene_Oxidative[i])

Oxidative_table<- table.v[index,]
#GO10 average
P10 <- t(colMeans(Oxidative_table))
rownames(P10) <- c("P10")

GO_10<- rbind(P1, P2, P3, P4,P5, P6, P7, P8,P9, P10)

GO_10_mean <- apply(GO_10, 1, mean)
GO_10_sd <- apply(GO_10, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
GO_10_zscore = (GO_10 - GO_10_mean)/GO_10_sd
Heatmap(GO_10_zscore, name = "Row Z-score", column_title = "Whole blood",cluster_columns = FALSE)


# #######################################
#calculate mean of each row
mean <- apply(table.v, 1, mean)
sd <- apply(table.v, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
zscore = (table.v - mean)/sd
colnames(zscore) <- c("Ctrl","Ctrl","Ctrl","Ctrl","Diabetes","Diabetes","Diabetes","Diabetes")


#ComplexHeatmap on zscore
###########################
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")

library(ComplexHeatmap)
#remove all rownames
rownames(zscore) <- c()
zscorefinal <- zscore
Heatmap(zscore, name = "Row Z-score", column_title = "Whole blood", cluster_columns = FALSE)
# 
# 
# #PCA 
# autoplot(prcomp(table.v))
# 
# 
# ######################################
# #gene list from Gene ontology based on Row Z-score 
# zscorefinale <- cbind (table.v[, 0], zscore)
# 
# zscorefinale_id<- cbind(rownames(zscorefinale), zscorefinale)
# rownames(zscorefinale_id) <- c()
# 
# #select the row Zscore for the corresponding genes names from zscorefinale matrix
# 
# #Ribosome biogenesis
# index <- integer(270)
# for (i in 1:270)
#   index[i] <- which(zscorefinale_id[,1] == gene_ribosome[i])
# 
# ribosome_table<- zscorefinale[index,]
# rownames(ribosome_table) <- c()
# Heatmap(ribosome_table, name = "Row Z-score", column_title = "Ribosome biogenesis")
# 
# autoplot(prcomp(ribosome_table))
# 
# 
# 
# #Ribonucleoprotein complex biogenesis
# index <- integer(425)
# for (i in 1:425)
#   index[i] <- which(zscorefinale_id[,1] == gene_ribonucleoprotein[i])
# 
# ribonucleoprotein_table<- zscorefinale[index,]
# rownames(ribonucleoprotein_table) <- c()
# Heatmap(ribonucleoprotein_table, name = "Row Z-score", column_title = "Ribonucleoprotein complex biogenesis")
# 
# 
# #mitochondrial gene expression
# index <- integer(154)
# for (i in 1:154)
#   index[i] <- which(zscorefinale_id[,1] == gene_mitochondrial[i])
# 
# mitochondiral_table<- zscorefinale[index,]
# rownames(mitochondiral_table) <- c()
# Heatmap(mitochondiral_table, name = "Row Z-score", column_title = "Mitochondrial gene expression")
# 
# 
# #translational initiation
# index <- integer(184)
# for (i in 1:184)
#   index[i] <- which(zscorefinale_id[,1] == gene_translational[i])
# 
# translational_table<- zscorefinale[index,]
# rownames(translational_table) <- c()
# Heatmap(translational_table, name = "Row Z-score", column_title = "Translational initiation")
# 
# 
# #establishment of protein localization to membrane
# index <- integer(301)
# for (i in 1:301)
#   index[i] <- which(zscorefinale_id[,1] == gene_proteinlocalization[i])
# 
# proteinlocalization_table<- zscorefinale[index,]
# rownames(proteinlocalization_table) <- c()
# Heatmap(proteinlocalization_table, name = "Row Z-score", column_title = "Establishment of protein localization to membrane")








#PCA
PCA<- prcomp(table.v)
biplot(PCA, scale = 0)
################################################################
#enrichment analysis: rowMean table as input in ConsensusPathDB
Ctrl_mean <- rowMeans(table.v[, 1:4])
Diabetes_mean <- rowMeans(table.v[, 5:8])
Enrichment_table <- cbind(Ctrl_mean, Diabetes_mean)
write.table(Enrichment_table, file="GSE60424_Enrichmentblood.txt", sep="\t", col.names = FALSE, quote=FALSE)
###############################################################



#differential expression (example:compare normal against obese)
contrast.matrix <- makeContrasts(Diff = Ctrl-Diabetes, levels = design) #take the difference between 2 groups
#ImFit fits a linear model using weighted least squares for each gene
fit <- lmFit(v, design)
#Estimate contrast for each gene
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)

#What genes are most differentially expressed?
Result_Table <- topTable(fit3, sort.by = "p", n= Inf)
head(Result_Table, 10)
#or Result_Table <- topTable(fit.3, coef = 'Diff') # this gives the top genes which give the difference
View(Result_Table)

Top_10 <- Result_Table[1:10, ]
write.table(Top_10, file="GSE60424_top10nonsignificant.csv", sep="\t", col.names = TRUE, quote=FALSE)



length(which(Result_Table$adj.P.Val< 0.05))
#[1] 0

sign_table <- Result_Table[which(Result_Table$adj.P.Val< 0.05),]


write.csv(sign_table, file="GSE60424_sign.txt", sep="\t", row.names = TRUE, quote=FALSE)
#write output into txt file
write.table(Result_Table, file="GSE60424_resultblood.txt", sep="\t", col.names = TRUE, quote=FALSE)




######################################################################################
#Select the samples from neutrophils
neutro_samples_diabetes<-table_neutro_diabetes$title
colindexneutro<-vector(mode="integer",length=length(neutro_samples_diabetes))
for (i in 1:length(neutro_samples_diabetes)){
  colindexneutro[i]<-which(colnames(Norma_data)==neutro_samples_diabetes[i])
}
Norma_data_neutro_diabetes<-Norma_data[,colindexneutro]


#remove rows where colums=0 (for all rows at least one is bigger than 0 -> filter out all rows with all 0)
finalneutro_Norma_data <- Norma_data_neutro_diabetes[rowSums(Norma_data_neutro_diabetes > 0) >= 1,]



#design the model to be fitted. before using voom because voom uses variances of the model residuals (observed-fitted)
#age = as.factor(table$`age:ch1`) as numeric
designneutro<-model.matrix(~0+table_neutro_diabetes$`diseasestatus:ch1`)
colnames(designneutro) <- c("Ctrl", "Diabetes")
#voom transformation into LinearModel
library(limma)
# command : BiocManager::install("limma", version = "3.8")
vneutro <- voom(finalneutro_Norma_data, designneutro, plot=TRUE )
table.vneutro <- vneutro$E #the log values from here for heatmap (without row Z score, st) and PCA (without row z score)

#differential expression (example:compare normal against obese)
contrast.matrixneutro <- makeContrasts(Diff = Ctrl-Diabetes, levels = designneutro) #take the difference between 2 groups
#ImFit fits a linear model using weighted least squares for each gene
fitneutro <- lmFit(vneutro, designneutro)
#Estimate contrast for each gene
fit2neutro <- contrasts.fit(fitneutro, contrast.matrixneutro)
fit3neutro <- eBayes(fit2neutro)

#What genes are most differentially expressed?
Result_Tableneutro <- topTable(fit3neutro, sort.by = "p", n= Inf)
#head(Result_Table, 10)
#or Result_Table <- topTable(fit.3, coef = 'Diff') # this gives the top genes which give the difference
#View(Result_Table)

length(which(Result_Tableneutro$adj.P.Val< 0.05))
#[1] 0

#sign_table <- Result_Table[which(Result_Table$adj.P.Val< 0.05),]

#write.csv(sign_table, file="GSE60424_sign.txt", sep="\t", row.names = TRUE, quote=FALSE)
#write output into txt file
write.table(Result_Tableneutro, file="GSE60424_resultneutro.txt", sep="\t", col.names = TRUE, quote=FALSE)






#Select the samples from monocytes
mono_samples_diabetes<-table_mono_diabetes$title
colindexmono<-vector(mode="integer",length=length(mono_samples_diabetes))
for (i in 1:length(mono_samples_diabetes)){
  colindexmono[i]<-which(colnames(Norma_data)==mono_samples_diabetes[i])
}
Norma_data_mono_diabetes<-Norma_data[,colindexmono]


#remove rows where colums=0 (for all rows at least one is bigger than 0 -> filter out all rows with all 0)
finalmono_Norma_data <- Norma_data_mono_diabetes[rowSums(Norma_data_mono_diabetes > 0) >= 1,]



#design the model to be fitted. before using voom because voom uses variances of the model residuals (observed-fitted)
#age = as.factor(table$`age:ch1`) as numeric
designmono<-model.matrix(~0+table_mono_diabetes$`diseasestatus:ch1`)
colnames(designmono) <- c("Ctrl", "Diabetes")
#voom transformation into LinearModel
library(limma)
# command : BiocManager::install("limma", version = "3.8")
vmono <- voom(finalmono_Norma_data, designmono, plot=TRUE )
table.vmono <- vmono$E #the log values from here for heatmap (without row Z score, st) and PCA (without row z score)

#differential expression (example:compare normal against obese)
contrast.matrixmono <- makeContrasts(Diff = Ctrl-Diabetes, levels = designmono) #take the difference between 2 groups
#ImFit fits a linear model using weighted least squares for each gene
fitmono <- lmFit(vmono, designmono)
#Estimate contrast for each gene
fit2mono <- contrasts.fit(fitmono, contrast.matrixmono)
fit3mono <- eBayes(fit2mono)

#What genes are most differentially expressed?
Result_Tablemono <- topTable(fit3mono, sort.by = "p", n= Inf)
#head(Result_Table, 10)
#or Result_Table <- topTable(fit.3, coef = 'Diff') # this gives the top genes which give the difference
#View(Result_Table)

length(which(Result_Tablemono$adj.P.Val< 0.05))
#[1] 0

#sign_table <- Result_Table[which(Result_Table$adj.P.Val< 0.05),]

#write.csv(sign_table, file="GSE60424_sign.txt", sep="\t", row.names = TRUE, quote=FALSE)
#write output into txt file
write.table(Result_Tablemono, file="GSE60424_resultmono.txt", sep="\t", col.names = TRUE, quote=FALSE)





#############################################################################################
#Select the samples from B-cells
B_samples_diabetes<-table_B_diabetes$title
colindexB<-vector(mode="integer",length=length(B_samples_diabetes))
for (i in 1:length(B_samples_diabetes)){
  colindexB[i]<-which(colnames(Norma_data)==B_samples_diabetes[i])
}
Norma_data_B_diabetes<-Norma_data[,colindexB]


#remove rows where colums=0 (for all rows at least one is bigger than 0 -> filter out all rows with all 0)
finalB_Norma_data <- Norma_data_B_diabetes[rowSums(Norma_data_B_diabetes > 0) >= 1,]



#design the model to be fitted. before using voom because voom uses variances of the model residuals (observed-fitted)
#age = as.factor(table$`age:ch1`) as numeric
designB<-model.matrix(~0+table_B_diabetes$`diseasestatus:ch1`)
colnames(designB) <- c("Ctrl", "Diabetes")
#voom transformation into LinearModel
library(limma)
# command : BiocManager::install("limma", version = "3.8")
vB <- voom(finalB_Norma_data, designB, plot=TRUE )
table.vB <- vB$E #the log values from here for heatmap (without row Z score, st) and PCA (without row z score)

#differential expression (example:compare normal against obese)
contrast.matrixB <- makeContrasts(Diff = Ctrl-Diabetes, levels = designB) #take the difference between 2 groups
#ImFit fits a linear model using weighted least squares for each gene
fitB <- lmFit(vB, designB)
#Estimate contrast for each gene
fit2B <- contrasts.fit(fitB, contrast.matrixB)
fit3B <- eBayes(fit2B)

#What genes are most differentially expressed?
Result_TableB <- topTable(fit3B, sort.by = "p", n= Inf)
#head(Result_Table, 10)
#or Result_Table <- topTable(fit.3, coef = 'Diff') # this gives the top genes which give the difference
#View(Result_Table)

length(which(Result_TableB$adj.P.Val< 0.05))
#[1] 0

#sign_table <- Result_Table[which(Result_Table$adj.P.Val< 0.05),]

#write.csv(sign_table, file="GSE60424_sign.txt", sep="\t", row.names = TRUE, quote=FALSE)
#write output into txt file
write.table(Result_TableB, file="GSE60424_resultB.txt", sep="\t", col.names = TRUE, quote=FALSE)





###########################################################################################
#Select the samples from CD4 cells
CD4_samples_diabetes<-table_CD4_diabetes$title
colindexCD4<-vector(mode="integer",length=length(CD4_samples_diabetes))
for (i in 1:length(CD4_samples_diabetes)){
  colindexCD4[i]<-which(colnames(Norma_data)==CD4_samples_diabetes[i])
}
Norma_data_CD4_diabetes<-Norma_data[,colindexCD4]


#remove rows where colums=0 (for all rows at least one is bigger than 0 -> filter out all rows with all 0)
finalCD4_Norma_data <- Norma_data_CD4_diabetes[rowSums(Norma_data_CD4_diabetes > 0) >= 1,]



#design the model to be fitted. before using voom because voom uses variances of the model residuals (observed-fitted)
#age = as.factor(table$`age:ch1`) as numeric
designCD4<-model.matrix(~0+table_CD4_diabetes$`diseasestatus:ch1`)
colnames(designCD4) <- c("Ctrl", "Diabetes")
#voom transformation into LinearModel
library(limma)
# command : BiocManager::install("limma", version = "3.8")
vCD4 <- voom(finalCD4_Norma_data, designCD4, plot=TRUE )
table.vCD4 <- vCD4$E #the log values from here for heatmap (without row Z score, st) and PCA (without row z score)
colnames(table.vCD4) <- c("Ctrl","Ctrl","Ctrl","Ctrl","Diabetes","Diabetes","Diabetes","Diabetes")
table.v_idCD4<- cbind(rownames(table.vCD4), table.vCD4)
##########################################################################
# Top 10 gene ontology terms between 5 & 500 according to q-values
#GO1 ribosome biogenesis
index <- integer(173)
for (i in 1:173)
  index[i] <- which(table.v_idCD4[,1] == gene_chromosome[i])

chromosome_table<- table.vCD4[index,]
#GO1 average
GO1 <- t(colMeans(chromosome_table))
rownames(GO1) <- c("GO1")

#GO2 nuclear division
index <- integer(351)
for (i in 1:351)
  index[i] <- which(table.v_idCD4[,1] == gene_nucleardivision[i])

nucleardivision_table<- table.vCD4[index,]
#GO2 average
GO2 <- t(colMeans(nucleardivision_table))
rownames(GO2) <- c("GO2")

#GO3 organelle fission 
index <- integer(386)
for (i in 1:386)
  index[i] <- which(table.v_idCD4[,1] == gene_organllefission[i])

organellefission_table<- table.vCD4[index,]
#GO3 average
GO3 <- t(colMeans(organellefission_table))
rownames(GO3) <- c("GO3")


#GO4 translational initiation
index <- integer(227)
for (i in 1:227)
  index[i] <- which(table.v_idCD4[,1] == gene_nuclearchromosome[i])

nuclearchromosome_table<- table.vCD4[index,]
#GO4 average
GO4 <- t(colMeans(nuclearchromosome_table))
rownames(GO4) <- c("GO4")

#GO5 DNA replication
index <- integer(251)
for (i in 1:251)
  index[i] <- which(table.v_idCD4[,1] == gene_DNAreplication[i])

DNAreplication_table<- table.vCD4[index,]
#GO5 average
GO5 <- t(colMeans(DNAreplication_table))
rownames(GO5) <- c("GO5")

#GO6 small molecule
index <- integer(341)
for (i in 1:341)
  index[i] <- which(table.v_idCD4[,1] == gene_smallmolecule[i])

smallmolecule_table<- table.vCD4[index,]
#GO6 average
GO6 <- t(colMeans(smallmolecule_table))
rownames(GO6) <- c("GO6")


#GO7 mitotic
index <- integer(259)
for (i in 1:259)
  index[i] <- which(table.v_idCD4[,1] == gene_mitotic[i])

mitotic_table<- table.vCD4[index,]
#GO7 average
GO7 <- t(colMeans(mitotic_table))
rownames(GO7) <- c("GO7")


#GO8 cell cycle 
index <- integer(206)
for (i in 1:206)
  index[i] <- which(table.v_idCD4[,1] == gene_cellcycle[i])

cellcycle_table<- table.vCD4[index,]
#GO8 average
GO8 <- t(colMeans(cellcycle_table))
rownames(GO8) <- c("GO8")


#GO9 nucleic acid
index <- integer(256)
for (i in 1:256)
  index[i] <- which(table.v_idCD4[,1] == gene_nucleicacid[i])

nucleicacid_table<- table.vCD4[index,]
#GO9 average
GO9 <- t(colMeans(nucleicacid_table))
rownames(GO9) <- c("GO9")


#GO10 sister chromatid
index <- integer(173)
for (i in 1:173)
  index[i] <- which(table.v_idCD4[,1] == gene_sister[i])

sister_table<- table.vCD4[index,]
#GO10 average
GO10 <- t(colMeans(sister_table))
rownames(GO10) <- c("GO10")

GO_10<- rbind(GO1, GO2, GO3, GO4,GO5, GO6, GO7, GO8, GO9, GO10)

GO_10_mean <- apply(GO_10, 1, mean)
GO_10_sd <- apply(GO_10, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
GO_10_zscore = (GO_10 - GO_10_mean)/GO_10_sd
Heatmap(GO_10_zscore, name = "Row Z-score", column_title = "CD4 T cells", cluster_columns = FALSE)

##############################################################
# Top 10 enrichment pathway between 5 & 500 according to q-values
#GO1 Cell cycle mitotic
index <- integer(419)
for (i in 1:419)
  index[i] <- which(table.v_idCD4[,1] == gene_Cell[i])

Cell_table<- table.vCD4[index,]
#GO1 average
P1 <- t(colMeans(Cell_table))
rownames(P1) <- c("P1")

#GO2 DNA repair
index <- integer(280)
for (i in 1:280)
  index[i] <- which(table.v_idCD4[,1] == gene_DNArepair[i])

DNArepair_table<- table.vCD4[index,]
#GO2 average
P2 <- t(colMeans(DNArepair_table))
rownames(P2) <- c("P2")

#GO3 Retinoblastoma Gene in Cancer 
index <- integer(86)
for (i in 1:86)
  index[i] <- which(table.v_idCD4[,1] == gene_Retinoblastoma[i])

Retinoblastoma_table<- table.vCD4[index,]
#GO3 average
P3 <- t(colMeans(Retinoblastoma_table))
rownames(P3) <- c("P3")


#GO4 Cell cycle checkpoints
index <- integer(215)
for (i in 1:215)
  index[i] <- which(table.v_idCD4[,1] == gene_Checkpoints[i])

Checkpoints_table<- table.vCD4[index,]
#GO4 average
P4 <- t(colMeans(Checkpoints_table))
rownames(P4) <- c("P4")

#GO5 Neddylation
index <- integer(212)
for (i in 1:212)
  index[i] <- which(table.v_idCD4[,1] == gene_Neddylation[i])

Neddylation_table<- table.vCD4[index,]
#GO5 average
P5 <- t(colMeans(Neddylation_table))
rownames(P5) <- c("P5")

#GO6 DNA Double strand break repair
index <- integer(135)
for (i in 1:135)
  index[i] <- which(table.v_idCD4[,1] == gene_DNAdouble[i])

DNAdouble_table<- table.vCD4[index,]
#GO6 average
P6 <- t(colMeans(DNAdouble_table))
rownames(P6) <- c("P6")


#GO7 Mphase
index <- integer(282)
for (i in 1:282)
  index[i] <- which(table.v_idCD4[,1] == gene_Mphase[i])

Mphase_table<- table.vCD4[index,]
#GO7 average
P7 <- t(colMeans(Mphase_table))
rownames(P7) <- c("P7")


#GO8 G2/M checkpoints
index <- integer(93)
for (i in 1:93)
  index[i] <- which(table.v_idCD4[,1] == gene_G2[i])

G2_table<- table.vCD4[index,]
#GO8 average
P8 <- t(colMeans(G2_table))
rownames(P8) <- c("P8")


#GO9 DNA replication
index <- integer(77)
for (i in 1:77)
  index[i] <- which(table.v_idCD4[,1] == gene_DNAreplication2[i])

DNAreplication2_table<- table.vCD4[index,]
#GO9 average
P9 <- t(colMeans(DNAreplication2_table))
rownames(P9) <- c("P9")


#GO10 Homology directed repair
index <- integer(110)
for (i in 1:110)
  index[i] <- which(table.v_idCD4[,1] == gene_Homology[i])

Homology_table<- table.vCD4[index,]
#GO10 average
P10 <- t(colMeans(Homology_table))
rownames(P10) <- c("P10")

GO_10<- rbind(P1, P2, P3, P4,P5, P6, P7, P8,P9, P10)

GO_10_mean <- apply(GO_10, 1, mean)
GO_10_sd <- apply(GO_10, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
GO_10_zscore = (GO_10 - GO_10_mean)/GO_10_sd
Heatmap(GO_10_zscore, name = "Row Z-score", column_title = "CD4 T cells", cluster_columns = FALSE)

#######################################


#calculate mean of each row 
meanCD4 <- apply(table.vCD4, 1, mean)
sdCD4 <- apply(table.vCD4, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
zscoreCD4 = (table.vCD4 - meanCD4)/sdCD4
colnames(zscoreCD4) <- c("Ctrl","Ctrl","Ctrl","Ctrl","Diabetes","Diabetes","Diabetes","Diabetes")
#remove all rownames
rownames(zscoreCD4) <- c()

##############################################
#ComplexHeatmap on table.v
Heatmap(zscoreCD4, name = "Row Z-score", column_title = "CD4 T cells", cluster_columns = FALSE)
PCA(zscoreCD4)

#select the row Zscore for the corresponding genes names from zscorefinale matrix

#Chromosome segregation
index <- integer(279)
for (i in 1:279)
  index[i] <- which(zscorefinale_id[,1] == gene_chromosome[i])

chromosome_table<- zscorefinale[index,]
rownames(chromosome_table) <- c()
Heatmap(chromosome_table, name = "Row Z-score", column_title = "Chromosome segregation")



#PCA
PCA<- prcomp(table.vCD4)
biplot(PCA, scale = 0)

#####################################################################
#enrichment analysis: rowMean table as input in ConsensusPathDB
Ctrl_mean_CD4 <- rowMeans(table.vCD4[, 1:4])
Diabetes_mean_CD4 <- rowMeans(table.vCD4[, 5:8])
Enrichment_table_CD4 <- cbind(Ctrl_mean_CD4, Diabetes_mean_CD4)
write.table(Enrichment_table_CD4, file="GSE60424_EnrichmentCD4.txt", sep="\t", col.names = FALSE, quote=FALSE)
#######################################################################

#differential expression (example:compare normal against obese)
contrast.matrixCD4 <- makeContrasts(Diff = Ctrl-Diabetes, levels = designCD4) #take the difference between 2 groups
#ImFit fits a linear model using weighted least squares for each gene
fitCD4 <- lmFit(vCD4, designCD4)
#Estimate contrast for each gene
fit2CD4 <- contrasts.fit(fitCD4, contrast.matrixCD4)
fit3CD4 <- eBayes(fit2CD4)

#What genes are most differentially expressed?
Result_TableCD4 <- topTable(fit3CD4, sort.by = "p", n= Inf)
#head(Result_Table, 10)
#or Result_Table <- topTable(fit.3, coef = 'Diff') # this gives the top genes which give the difference
#View(Result_Table)

length(which(Result_TableCD4$adj.P.Val< 0.05))
#[1] 0

#sign_table <- Result_Table[which(Result_Table$adj.P.Val< 0.05),]

#write.csv(sign_table, file="GSE60424_sign.txt", sep="\t", row.names = TRUE, quote=FALSE)
#write output into txt file
write.table(Result_TableCD4, file="GSE60424_resultCD4.txt", sep="\t", col.names = TRUE, quote=FALSE)



######################################################################################
#Select the samples from CD8 cells
CD8_samples_diabetes<-table_CD8_diabetes$title
colindexCD8<-vector(mode="integer",length=length(CD8_samples_diabetes))
for (i in 1:length(CD8_samples_diabetes)){
  colindexCD8[i]<-which(colnames(Norma_data)==CD8_samples_diabetes[i])
}
Norma_data_CD8_diabetes<-Norma_data[,colindexCD8]


#remove rows where colums=0 (for all rows at least one is bigger than 0 -> filter out all rows with all 0)
finalCD8_Norma_data <- Norma_data_CD8_diabetes[rowSums(Norma_data_CD8_diabetes > 0) >= 1,]



#design the model to be fitted. before using voom because voom uses variances of the model residuals (observed-fitted)
#age = as.factor(table$`age:ch1`) as numeric
designCD8<-model.matrix(~0+table_CD8_diabetes$`diseasestatus:ch1`)
colnames(designCD8) <- c("Ctrl", "Diabetes")
#voom transformation into LinearModel
library(limma)
# command : BiocManager::install("limma", version = "3.8")
vCD8 <- voom(finalCD8_Norma_data, designCD8, plot=TRUE )
table.vCD8 <- vCD8$E #the log values from here for heatmap (without row Z score, st) and PCA (without row z score)
colnames(table.vCD8) <- c("Ctrl","Ctrl","Ctrl","Ctrl","Diabetes","Diabetes","Diabetes","Diabetes")
table.v_idCD8<- cbind(rownames(table.vCD8), table.vCD8)
###############################################################################################
# Top 18 gene ontology terms between 5 & 500 according to q-values
#GO1 chromosome segregation
index <- integer(277)
for (i in 1:277)
  index[i] <- which(table.v_idCD8[,1] == chromosome_segregation[i])

chromosome_table<- table.vCD8[index,]
#GO1 average
GO1 <- t(colMeans(chromosome_table))
rownames(GO1) <- c("GO1")

#GO2 organelle fission
index <- integer(383)
for (i in 1:383)
  index[i] <- which(table.v_idCD8[,1] == organelle_fission[i])

organelle_fission_table<- table.vCD8[index,]
#GO2 average
GO2 <- t(colMeans(organelle_fission_table))
rownames(GO2) <- c("GO2")

#GO3 mitochondrial gene expression 
index <- integer(154)
for (i in 1:154)
  index[i] <- which(table.v_idCD8[,1] == mitochondrial_gene_expression[i])

mitochondrial_gene_expression_table<- table.vCD8[index,]
#GO3 average
GO3 <- t(colMeans(mitochondrial_gene_expression_table))
rownames(GO3) <- c("GO3")


#GO4 ribonucleoprotein complex biogenesis
index <- integer(425)
for (i in 1:425)
  index[i] <- which(table.v_idCD8[,1] == ribonucleo_protein[i])

ribonucleo_protein_table<- table.vCD8[index,]
#GO4 average
GO4 <- t(colMeans(ribonucleo_protein_table))
rownames(GO4) <- c("GO4")

#GO5 nuclear division
index <- integer(348)
for (i in 1:348)
  index[i] <- which(table.v_idCD8[,1] == nuclear_division[i])

nuclear_division_table<- table.vCD8[index,]
#GO5 average
GO5 <- t(colMeans(nuclear_division_table))
rownames(GO5) <- c("GO5")

#GO6 nuclear chromosome segregation
index <- integer(225)
for (i in 1:225)
  index[i] <- which(table.v_idCD8[,1] == nuclear_chromosome_segregation[i])

nuclear_chromosome_segregation_table<- table.vCD8[index,]
#GO6 average
GO6 <- t(colMeans(nuclear_chromosome_segregation_table))
rownames(GO6) <- c("GO6")


#GO7 mitotic nuclear division
index <- integer(257)
for (i in 1:257)
  index[i] <- which(table.v_idCD8[,1] == mitotic_nuclear_division[i])

mitotic_nuclear_division_table<- table.vCD8[index,]
#GO7 average
GO7 <- t(colMeans(mitotic_nuclear_division_table))
rownames(GO7) <- c("GO7")


#GO8 translational elongation
index <- integer(128)
for (i in 1:128)
  index[i] <- which(table.v_idCD8[,1] == translational_elongation[i])

translational_elongation_table<- table.vCD8[index,]
#GO8 average
GO8 <- t(colMeans(translational_elongation_table))
rownames(GO8) <- c("GO8")


#GO9 mitotic cell cycle phase transition
index <- integer(455)
for (i in 1:455)
  index[i] <- which(table.v_idCD8[,1] == mitotic_cell_cycle[i])

mitotic_cell_cycle_table<- table.vCD8[index,]
#GO9 average
GO9 <- t(colMeans(mitotic_cell_cycle_table))
rownames(GO9) <- c("GO9")


#GO10 sister chromatid segregation
index <- integer(172)
for (i in 1:172)
  index[i] <- which(table.v_idCD8[,1] == sister_chromatid [i])

sister_chromatid_table<- table.vCD8[index,]
#GO10 average
GO10 <- t(colMeans(sister_chromatid_table))
rownames(GO10) <- c("GO10")

#GO11 ribosome biogenesis
index <- integer(269)
for (i in 1:269)
  index[i] <- which(table.v_idCD8[,1] == ribosome_biogenesis[i])

ribosome_biogenesis_table<- table.vCD8[index,]
#GO10 average
GO11 <- t(colMeans(ribosome_biogenesis_table))
rownames(GO11) <- c("GO11")


#GO12 nucleoside triphosphate metabolic process
index <- integer(270)
for (i in 1:270)
  index[i] <- which(table.v_idCD8[,1] == nucleoside_triphosphate[i])

nucleoside_triphosphate_table<- table.vCD8[index,]
#GO10 average
GO12 <- t(colMeans(nucleoside_triphosphate_table))
rownames(GO12) <- c("GO12")



#GO13 protein-DNA complex subunit organization
index <- integer(207)
for (i in 1:207)
  index[i] <- which(table.v_idCD8[,1] == protein_DNA_complex[i])

protein_DNA_complex_table<- table.vCD8[index,]
#GO10 average
GO13 <- t(colMeans(protein_DNA_complex_table))
rownames(GO13) <- c("GO13")



#GO14 protein-DNA complex assembly
index <- integer(178)
for (i in 1:178)
  index[i] <- which(table.v_idCD8[,1] == protein_DNA_assembly[i])

protein_DNA_assembly_table<- table.vCD8[index,]
#GO10 average
GO14 <- t(colMeans(protein_DNA_assembly_table))
rownames(GO14) <- c("GO14")



#GO15 DNA replication
index <- integer(253)
for (i in 1:253)
  index[i] <- which(table.v_idCD8[,1] == DNA_replication[i])

DNA_replication_table<- table.vCD8[index,]
#GO10 average
GO15 <- t(colMeans(DNA_replication_table))
rownames(GO15) <- c("GO15")



#GO16 mitotic sister chromatid segregation
index <- integer(147)
for (i in 1:147)
  index[i] <- which(table.v_idCD8[,1] == mitotic_sister_chromatid[i])

mitotic_sister_chromatid_table<- table.vCD8[index,]
#GO10 average
GO16 <- t(colMeans(mitotic_sister_chromatid_table))
rownames(GO16) <- c("GO16")


#GO17 cell cycle checkpoint 
index <- integer(205)
for (i in 1:205)
  index[i] <- which(table.v_idCD8[,1] == cell_cycle_checkpoint[i])

cell_cycle_checkpoint_table<- table.vCD8[index,]
#GO10 average
GO17 <- t(colMeans(cell_cycle_checkpoint_table))
rownames(GO17) <- c("GO17")


#GO18 oxidative phosphorylation
index <- integer(118)
for (i in 1:118)
  index[i] <- which(table.v_idCD8[,1] == oxidative_phosphorylation[i])

oxidative_phosphorylation_table<- table.vCD8[index,]
#GO10 average
GO18 <- t(colMeans(oxidative_phosphorylation_table))
rownames(GO18) <- c("GO18")




GO_sum<- rbind(GO1, GO2, GO3, GO4,GO5, GO6, GO7, GO8, GO9, GO10, GO11, GO12, GO13, GO14, GO15, GO16, GO17, GO18)

GO_sum_mean <- apply(GO_sum, 1, mean)
GO_sum_sd <- apply(GO_sum, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
GO_sum_zscore = (GO_sum - GO_sum_mean)/GO_sum_sd
Heatmap(GO_sum_zscore, name = "Row Z-score", column_title = "CD8 T cells", cluster_columns = FALSE)

##############################################################
# Top 10 enrichment pathway between 5 & 500 according to q-values
#GO1 Translation
index <- integer(291)
for (i in 1:291)
  index[i] <- which(table.v_idCD8[,1] == Translation1[i])

Cell_table<- table.vCD8[index,]
#GO1 average
P1 <- t(colMeans(Cell_table))
rownames(P1) <- c("P1")

#GO2 Cell Cycle, Mitotic
index <- integer(418)
for (i in 1:418)
  index[i] <- which(table.v_idCD8[,1] == Cell_Cycle1[i])

Cell_Cycle1_table<- table.vCD8[index,]
#GO2 average
P2 <- t(colMeans(Cell_Cycle1_table))
rownames(P2) <- c("P2")

  #GO3 Cell Cycle Checkpoints
index <- integer(215)
for (i in 1:215)
  index[i] <- which(table.v_idCD8[,1] == CellCheckpoints1[i])

CellCheckpoints1_table<- table.vCD8[index,]
#GO3 average
P3 <- t(colMeans(CellCheckpoints1_table))
rownames(P3) <- c("P3")


#GO4 Risobome
index <- integer(126)
for (i in 1:126)
  index[i] <- which(table.v_idCD8[,1] == Ribosome1[i])

Ribosome1_table<- table.vCD8[index,]
#GO4 average
P4 <- t(colMeans(Ribosome1_table))
rownames(P4) <- c("P4")

#GO5 Mitochondrial translation
index <- integer(90)
for (i in 1:90)
  index[i] <- which(table.v_idCD8[,1] == Mitochondrial1[i])

Mitochondrial1_table<- table.vCD8[index,]
#GO5 average
P5 <- t(colMeans(Mitochondrial1_table))
rownames(P5) <- c("P5")

#GO6 Neddylation
index <- integer(211)
for (i in 1:211)
  index[i] <- which(table.v_idCD8[,1] == Neddylation1[i])

Neddylation1_table<- table.vCD8[index,]
#GO6 average
P6 <- t(colMeans(Neddylation1_table))
rownames(P6) <- c("P6")


#GO7 Respiratory electron transport
index <- integer(121)
for (i in 1:121)
  index[i] <- which(table.v_idCD8[,1] == Respiratory_electron_transport1[i])

Respiratory_electron_transport1_table<- table.vCD8[index,]
#GO7 average
P7 <- t(colMeans(Respiratory_electron_transport1_table))
rownames(P7) <- c("P7")


#GO8 Metabolism of amino acids and derivatives
index <- integer(265)
for (i in 1:265)
  index[i] <- which(table.v_idCD8[,1] == Metabolism_of_amino_acids1[i])

Metabolism_of_amino_acids1_table<- table.vCD8[index,]
#GO8 average
P8 <- t(colMeans(Metabolism_of_amino_acids1_table))
rownames(P8) <- c("P8")


#GO9 Processing if CappedIntron-Containing Pre-mRNA
index <- integer(231)
for (i in 1:231)
  index[i] <- which(table.v_idCD8[,1] == Processing_of_Capped1[i])

Processing_of_Capped1_table<- table.vCD8[index,]
#GO9 average
P9 <- t(colMeans(Processing_of_Capped1_table))
rownames(P9) <- c("P9")


#GO10 Electron Transport Chain
index <- integer(99)
for (i in 1:99)
  index[i] <- which(table.v_idCD8[,1] == Electron_Transport_Chai1[i])

Electron_Transport_Chai1_table<- table.vCD8[index,]
#GO10 average
P10 <- t(colMeans(Electron_Transport_Chai1_table))
rownames(P10) <- c("P10")


#GO11 Oxidative phosphorylation
index <- integer(117)
for (i in 1:117)
  index[i] <- which(table.v_idCD8[,1] == Oxidative_phosphorylation1[i])

Oxidative_phosphorylation1_table<- table.vCD8[index,]
#GO10 average
P11 <- t(colMeans(Oxidative_phosphorylation1_table))
rownames(P11) <- c("P11")



GO_average<- rbind(P1, P2, P3, P4,P5, P6, P7, P8, P9, P10, P11)

GO_average_mean <- apply(GO_average, 1, mean)
GO_average_sd <- apply(GO_average, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
GO_average_zscore = (GO_average - GO_average_mean)/GO_average_sd
Heatmap(GO_average_zscore, name = "Row Z-score", column_title = "CD8 T cells", cluster_columns = FALSE)

#######################################


#calculate mean of each row 
meanCD8 <- apply(table.vCD8, 1, mean)
sdCD8 <- apply(table.vCD8, 1, sd)
#calculate row-Zscores for heatmap presentation from mean and sd
zscoreCD8 = (table.vCD8 - meanCD8)/sdCD8
colnames(zscoreCD8) <- c("Ctrl","Ctrl","Ctrl","Ctrl","Diabetes","Diabetes","Diabetes","Diabetes")
#remove all rownames
rownames(zscoreCD8) <- c()

##############################################
#ComplexHeatmap on table.v
Heatmap(zscoreCD8, name = "Row Z-score", column_title = "CD8 T cells", cluster_columns = FALSE)
PCA(zscoreCD4)

#select the row Zscore for the corresponding genes names from zscorefinale matrix

#Chromosome segregation
index <- integer(279)
for (i in 1:279)
  index[i] <- which(zscorefinale_id[,1] == gene_chromosome[i])

chromosome_table<- zscorefinale[index,]
rownames(chromosome_table) <- c()
Heatmap(chromosome_table, name = "Row Z-score", column_title = "Chromosome segregation")

##########################################################################################
#enrichment analysis: rowMean table as input in ConsensusPathDB
Ctrl_mean_CD8 <- rowMeans(table.vCD8[, 1:4])
Diabetes_mean_CD8 <- rowMeans(table.vCD8[, 5:8])
Enrichment_table_CD8 <- cbind(Ctrl_mean_CD8, Diabetes_mean_CD8)
write.table(Enrichment_table_CD8, file="GSE60424_EnrichmentCD8.txt", sep="\t", col.names = FALSE, quote=FALSE)
#######################################################################




#differential expression (example:compare normal against obese)
contrast.matrixCD8 <- makeContrasts(Diff = Ctrl-Diabetes, levels = designCD8) #take the difference between 2 groups
#ImFit fits a linear model using weighted least squares for each gene
fitCD8 <- lmFit(vCD8, designCD8)
#Estimate contrast for each gene
fit2CD8 <- contrasts.fit(fitCD8, contrast.matrixCD8)
fit3CD8 <- eBayes(fit2CD8)

#What genes are most differentially expressed?
Result_TableCD8 <- topTable(fit3CD8, sort.by = "p", n= Inf)
#head(Result_Table, 10)
#or Result_Table <- topTable(fit.3, coef = 'Diff') # this gives the top genes which give the difference
#View(Result_Table)

length(which(Result_TableCD8$adj.P.Val< 0.05))
#[1] 0

#sign_table <- Result_Table[which(Result_Table$adj.P.Val< 0.05),]

#write.csv(sign_table, file="GSE60424_sign.txt", sep="\t", row.names = TRUE, quote=FALSE)
#write output into txt file
write.table(Result_TableCD8, file="GSE60424_resultCD8.txt", sep="\t", col.names = TRUE, quote=FALSE)




