######################################################################
##Prepare data set
######################################################################

#For Microarray data:
##For Affymetrix, download raw (tar) data files and perform RMA-normalization
##For Agilent, Lowess normalization and log2 transformation
##If combining data from multiple datasets/experiments, consider median- or quantile-normalization

#For RNAseq count data:
##Consider removing genes with maximum 4 or fewer counts.
##Consider RNAseq-specific normalization, such as DESeq2 or EdgeR rather than count data. 

#Note: rownames are Entrez IDs, colnames are samples

#Example
library(DESeq2)
MyCountData <- read.table("countTable.txt", row.names="SYMBOL", header=TRUE)
CountData=as.matrix(MyCountData) #Raw count data with integer values and rownames are gene symbols
condition=factor(colnames(CountData))
condition_df=data.frame(condition)

#Get DESeq normalized data
dds=DESeqDataSetFromMatrix(CountData, condition_df, ~ condition)
dds <- estimateSizeFactors(dds)
DESeq2_rlog_norm = assay(rlog( dds )) #The internal rlog normalization - preferred

#OR you can do log2 transformation with a pseudocount. Also acceptable
DESeq2Norm_data <- counts(dds, normalized=TRUE)
DESeq2Norm_data_log2 <- log2(norm.counts + 1)

#Change gene symbols to more stable Entrez gene ID's
library(AnnotationDbi)
library(org.Hs.eg.db)
EntrezSymbol=toTable(org.Hs.egSYMBOL)
DESeq2_rlog_Entrez=merge(EntrezSymbol, DESeq2_rlog_norm, by.x='symbol', by.y=0, all=F)
row.names(DESeq2_rlog_Entrez)= DESeq2_rlog_Entrez$gene_id
DESeq2_rlog_Entrez_final=DESeq2_rlog_Entrez[,-c(1:2)]

data=DESeq2_rlog_Entrez_final


###################################################### 
###Import Signatures
######################################################

SignaturesUp=as.data.frame(read.table('SignatureList_Up_EntrezID.txt', sep='\t', header=TRUE))
SigsList_Up=lapply(SignaturesUp, function (x) x[!is.na(x)])

SignaturesDown=as.data.frame(read.table('SignatureList_Down_EntrezID.txt', sep='\t', header=TRUE))
SigsList_Down=lapply(SignaturesDown, function (x) x[!is.na(x)])

#Signatures with 'down' component
DownComponent=c('ACIDOSIS.21672245', 'BCAT.ONCOGENE.16273092', 'BRCA1.11823860', 'BRCA1.UP.16843262', 'CCND1.12914697', 'CCND1K112E.12914697', 'CD44POS.CD24NEG.17229949', 'E2F3.ONCOGENE.16273092', 'EIF4E.17638893', 'FIBRO.SERUM1.14737219', 'FIBRO.SERUM2.14737219', 'GENE70.11823860', 'GGI.16478745', 'GLUCOSEDEPLETION.21672245', 'HSP90.20920318', 'HYPOXIA.21672245', 'LACTICACIDOSIS.21672245', 'MAMLIN.LumProg.19648928', 'MAMLIN.MaSC.19648928', 'MAMLIN.MatLum.19648928', 'MAMLIN.STROMAL.19648928', 'MAPK.16585219', 'MYC.ONCOGENE.16273092', 'MYCSIG.19690609', 'PDGF.23583284', 'PIK3CA.MUT.20479250', 'PTENDEL.17452630', 'RAS.ACT.20591134', 'RAS.ONCOGENE.16273092', 'RBBP8.UP.16843262', 'SRC.ONCOGENE.16273092', 'TGFBETA.18394990', 'YAP.UP.18413746')


###########################################
###Calculate Expression Signatures Score###
###########################################


data=DESeq2_rlog_Entrez_final #rownames are Entrez IDs, colnames are samples
ExpressionScore=as.data.frame(colnames(data)) #Create data frame to accept scores
colnames(ExpressionScore)=c('fileName')

for(i in (1:length(SigsList_Up)))
{
  upList=ls(SigsList_Up[i])
  downList=ls(SigsList_Down[i])
  print(upList)
  print(downList)
  UpListDF=as.data.frame(SigsList_Up[[i]])
  colnames(UpListDF)=upList
  DownListDF=as.data.frame(SigsList_Down[[i]])
  colnames(DownListDF)=downList
  sigDataUp=merge(UpListDF,data,by.x=upList,by.y=0, all=FALSE)
  sigDataDown=merge(DownListDF,data,by.x=downList,by.y=0, all=FALSE)
  UpScore=as.data.frame(apply(sigDataUp[,2:length(colnames(sigDataUp))],2,mean, na.rm=T))
  DownScore=as.data.frame(apply(sigDataDown[,2:length(colnames(sigDataDown))],2,mean, na.rm=T))
  DownScore2=ifelse(is.na(DownScore[,1]),0,DownScore[,1])
  FinalScore=UpScore-DownScore2
  colnames(FinalScore)= upList
  ExpressionScore=merge(ExpressionScore,FinalScore,by.x='fileName',by.y=0, all=FALSE) #OR AgilentExpressionScore
}

write.table(ExpressionScore, file='ExpressionScore.txt', sep='\t', row.names=F)
