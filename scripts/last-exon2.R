library(ggplot2)
library(reshape2)
library(dplyr)
library(biomaRt)
library(changepoint)
library(argparse)
library(seqinr)
library(EDASeq)


###  Arguments Handler  ###
parser <- ArgumentParser(description='Build table of samples along with binary indicators of usage.')
parser$add_argument("--vectorPB", type="character", nargs="+",
                    default="depth/vector_PB.txt",
                    help="Path to Pbodies depth of coverage vectors.")
parser$add_argument("--vectorCYTO", type="character", nargs="+", default="depth/vector_CYTO.txt",
                    help="Path to Cytoplasmic depth of coverage vectors.")
parser$add_argument("--bed", type="character", nargs="+", default="data/last_exon_final.bed",
                    help="Complete Last Exon BED file.")
parser$add_argument("--fasta", type="character", nargs="+", default="data/signif_last_exon.fa",
                    help="Filtered Last Exon Fasta file.")
parser$add_argument("--output", type="character", default="results/",
                    help="Paths to output folder.")
parser$add_argument('--log', type="character", help='Path to log file.')
args <- parser$parse_args()

###  Utils  ###

logging <- function(str) {
  sink(file=paste(args$log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}

### Import and Format Datas ###
logging("Importing input datas")
#Import raw files
listup<-scan(args$vectorPB,sep = "\n",what="")
listdown<-scan(args$vectorCYTO,sep = "\n",what="")
last_exon_final <- read.delim(args$bed, header=FALSE)
enstup<-c()
list2up<-list()

for (i in  1:length(listup)){
  
  enstup<-c(enstup,strsplit(listup[i]," ")[[1]][1])
  concatup<-strsplit(listup[i]," ")[1]
  list2up[[i]]<-as.numeric(concatup[[1]][-1])
  
}
names(list2up)<-enstup

enstdown<-c()
list2down<-list()

for (i in  1:length(listdown)){
  
  enstdown<-c(enstdown,strsplit(listdown[i]," ")[[1]][1])
  concatdown<-strsplit(listdown[i]," ")[1]
  list2down[[i]]<-as.numeric(concatdown[[1]][-1])
  
}
names(list2down)<-enstdown

logging("Formating datas")
# Select Sens and Antisens genes
PosList<-last_exon_final%>%
  dplyr::filter(V6=="+")%>%
  dplyr::select(V4)
NegList<-last_exon_final%>%
  dplyr::filter(V6=="-")%>%
  dplyr::select(V4)

# Reverse Antisens genes
ListSens<-list2up[names(list2up)%in%PosList$V4]
ListAntiSens<-list2up[names(list2up)%in%NegList$V4]
for (i in 1:length(ListAntiSens)){
  ListAntiSens[[i]]<-rev(ListAntiSens[[i]])
}
list2up<-c(ListSens,ListAntiSens)

ListSens<-list2down[names(list2down)%in%PosList$V4]
ListAntiSens<-list2down[names(list2down)%in%NegList$V4]
for (i in 1:length(ListAntiSens)){
  ListAntiSens[[i]]<-rev(ListAntiSens[[i]])
}
list2down<-c(ListSens,ListAntiSens)

#Erasing non-covered 3' part of the exon
for (i in names(list2up)) {
  for (pos in length(list2up[[i]]):1){
    if (list2up[[i]][pos]+list2down[[i]][pos]==0){
      list2up[[i]]<-list2up[[i]][-pos]
      list2down[[i]]<-list2down[[i]][-pos]
    }
    else {
      break()
    }
  }
}

logging("Calculating ratio")
###  Ratio Pbodies/Presort  ###
list2final<-list()
for (x in names(list2up)) {
  list2final[[x]]<-log(list2up[[x]]/list2down[[x]])
  list2final[[x]][is.nan(list2final[[x]])]<-0
}

NormalizedList<-list()
for (j in names(list2final)){
  vec<-list2final[[j]]
  breaks<-round(seq(1,length(vec),length.out=100))
  WholeMeanExtract<-c()
  for(x in 1:(length(breaks)-1)){
    extract<-list2final[[j]][c(breaks[x]:breaks[x+1])]
    meanExtract<-mean(as.numeric(extract))
    WholeMeanExtract<-c(WholeMeanExtract,meanExtract)
  }
  NormalizedList[[j]]<-WholeMeanExtract
}
showOnPlot<-round(seq(1,99,length.out=25))
NormalizedListSOP<-list()
for (k in names(NormalizedList)){
  NormalizedListSOP[[k]]<-NormalizedList[[k]][showOnPlot]
}

MatrixNormalized<-matrix(unlist(NormalizedListSOP), ncol = 25, byrow = TRUE)
rownames(MatrixNormalized)<-names(NormalizedListSOP)
colnames(MatrixNormalized)<-showOnPlot
Melted<-reshape2::melt(as.data.frame(MatrixNormalized))

pdf(paste0(args$output,"/ratioAllGenes.pdf",collapse = ""))
ratioplot<-ggplot(Melted,aes(variable,value))+geom_boxplot()+ggtitle("Last Exon IN/OUT ratio of depth of coverage (4174 genes)")+labs(y="log(DepthIN+1/DepthOUT+1)",x="Position on last exon (lenghts normalized to 100")
dev.off()


###  Find most significant events  ###
logging("Starting Chi2 Method")

# Filter event with low 3' coverage
goin_up<-c()
filteredlist<-list()
for (j in names(NormalizedList)){
  if (length(which(is.na(list2final[[j]])))==0){
    filteredlist[[j]]<-list2final[[j]]
    if ((tail(list2down[[j]],n=1)+tail(list2up[[j]],n=1))>5){
        goin_up<-c(goin_up,j)
    }
  }
}

# Statistical test (ChiSquare)
chisqlist<-list()
oddsratio<-list()
see<-c()
for (k in goin_up) {
    len_exon<-length(list2up[[k]])
    one_fifth<-round(len_exon/100*20)
    two_fifth<-round(len_exon/100*40)
    three_fifth<-round(len_exon/100*60)
    four_fifth<-round(len_exon/100*80)
    limits<-c(two_fifth,three_fifth,four_fifth,len_exon)
    chisqlist[[k]]<-1
    for (l in 1:(length(limits)-1)){
      contingency<-as.table(rbind(c(mean(list2up[[k]][1:one_fifth]), mean(list2down[[k]][1:one_fifth])), c(mean(list2up[[k]][limits[l]:limits[l+1]]),mean(list2down[[k]][limits[l]:limits[l+1]]))))
      dimnames(contingency) <- list(position=c("Start","End"), condition=c("Pbodies","Presort"))
    test<-chisq.test(contingency)
    if (!is.na(test$p.value) && test$p.value<chisqlist[[k]]){
      chisqlist[[k]]<-test$p.value
      oddsratio[[k]]<-(contingency[2,1]/contingency[1,1])-(contingency[2,2]/contingency[1,2])
    }
  }
}
chisqfinal<-as.data.frame(t(as.data.frame(chisqlist)))
oddsfinal<-t(as.data.frame(oddsratio))
chisqadjust<-as.data.frame(p.adjust(chisqfinal[,1]))
rownames(chisqadjust)<-rownames(chisqfinal)
significant1<-rownames(chisqadjust)[which(chisqadjust[,1]<=0.05)] 

#List of last exons with no ovelaping gene / List generated by an "utils" script
significant<- c("ENST00000196371","ENST00000209884","ENST00000216297","ENST00000216832","ENST00000219097","ENST00000219905","ENST00000221130","ENST00000221200","ENST00000229854","ENST00000240316","ENST00000257787","ENST00000258341","ENST00000260641","ENST00000261858","ENST00000261875","ENST00000263239","ENST00000264220","ENST00000265333","ENST00000271850","ENST00000274242","ENST00000285379","ENST00000296411","ENST00000296503","ENST00000296701","ENST00000300870","ENST00000303319","ENST00000303436","ENST00000312377","ENST00000313899","ENST00000313949","ENST00000314191","ENST00000315717","ENST00000320848","ENST00000321442","ENST00000331738","ENST00000334701","ENST00000335007","ENST00000336032","ENST00000336053","ENST00000336314","ENST00000336812","ENST00000337288","ENST00000337620","ENST00000337859","ENST00000338179","ENST00000339399","ENST00000339697","ENST00000340513","ENST00000344034","ENST00000344756","ENST00000347630","ENST00000348943","ENST00000351217","ENST00000354185","ENST00000356053","ENST00000356221","ENST00000356840","ENST00000358794","ENST00000360541","ENST00000360663","ENST00000361078","ENST00000361436","ENST00000366922","ENST00000367047","ENST00000367075","ENST00000367142","ENST00000367241","ENST00000367874","ENST00000368321","ENST00000368690","ENST00000369577","ENST00000370339","ENST00000370355","ENST00000371370","ENST00000371429","ENST00000373232","ENST00000373548","ENST00000373647","ENST00000373719","ENST00000373734","ENST00000374171","ENST00000375370","ENST00000376499","ENST00000378455","ENST00000379160","ENST00000379535","ENST00000388711","ENST00000395913","ENST00000401878","ENST00000421503","ENST00000422000","ENST00000423485","ENST00000426480","ENST00000428425","ENST00000433176","ENST00000448750","ENST00000455997","ENST00000457657","ENST00000469257","ENST00000476371","ENST00000497289","ENST00000502553","ENST00000504154","ENST00000506538","ENST00000517671","ENST00000535419","ENST00000535601","ENST00000536438","ENST00000549190","ENST00000549855","ENST00000556029","ENST00000576409","ENST00000578681","ENST00000584379","ENST00000590758","ENST00000610854","ENST00000612630","ENST00000617365","ENST00000620073","ENST00000622999","ENST00000625464","ENST00000640418","ENST00000642050","ENST00000646102","ENST00000646647","ENST00000648638","ENST00000650070","ENST00000662716","ENST00000672724","ENST00000672937","ENST00000674920","ENST00000677213","ENST00000677446")

###  Detection of inflexion point  ###
logging("Starting Changepoint Method")
changepoints<-list()
for (x in significant){
  if(any(is.infinite(list2final[[x]]))){
  }
  else if (any(is.na(list2final[[x]]))){
  }
  else {
    lo<-loess(list2final[[x]]~seq(1:length(list2final[[x]])))
    xl <- seq(min(seq(1:length(list2final[[x]]))),max(seq(1:length(list2final[[x]]))))
    out = predict(lo,xl)
    changepoints[[x]]<-cpt.mean(out)@cpts[1]
  }
}
changepointsDT<-as.data.frame(changepoints)
changepointsDT[2,]<-colnames(changepointsDT)
finalTable<-as.data.frame(t(changepointsDT))
colnames(finalTable)<-c("changepoint","ENST")


pdf(paste0(args$output,"/inflection_est.pdf",collapse = ""))
for (i in colnames(changepointsDT)){
  plot(list2final[[i]],main=i)
  abline(v=changepointsDT[1,i])
}
dev.off()


###  Fortmat datas for output  ###
logging("Generating final output")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_transcript_id", attributes= c("ensembl_transcript_id","hgnc_symbol"),
                  values = finalTable$ENST, mart= mart)

last_exon_final$V7<-apply( last_exon_final[ , c(2:3) ] , 1 , paste , collapse = "-" )
last_exon_final$V8<-apply( last_exon_final[ , c(1,7) ] , 1 , paste , collapse = ":" )
coord_table<-last_exon_final[,c(4,8)]
colnames(coord_table)<-c("ENST","coord")

chisqfinal$ENST<-rownames(chisqfinal)
colnames(chisqfinal)<-c("chisq.padj","ENST")

#oddsfinal$ENST<-rownames(oddsfinal)
#colnames(oddsfinal)<-c("ratio","ENST")
  
exportTable1<-merge(x = finalTable,y=gene_IDs,by.x = "ENST",by.y = "ensembl_transcript_id")
exportTable2<-merge(x = exportTable1,y=coord_table,by= "ENST",all=FALSE)
exportTable3<-merge(x = exportTable2,y=chisqfinal,by= "ENST",all=FALSE)
#exportTable3bis<-merge(x = exportTable3,y=oddsfinal,by= "ENST",all=FALSE)

###  Parse fasta sequence to create fasta primers  ###

Handler<-read.fasta(file=args$fasta, seqtype = "DNA",as.string = TRUE,forceDNAtolower = FALSE)
preprimer<-list()
postprimer<-list()
for (gene in names(Handler)) {
  id<-substr(gene,1,15)
  cut<-as.numeric(exportTable3[exportTable3$ENST==id,"changepoint"])
  preprimer[[id]]<-substr(Handler[[gene]][1],1,cut)
  postprimer[[id]]<-substr(Handler[[gene]][1],cut+1,nchar(Handler[[gene]][1]))
}
preprimerDF<-as.data.frame(t(as.data.frame(preprimer)))
postprimerDF<-as.data.frame(t(as.data.frame(postprimer)))
preprimerDF$V2<-rownames(preprimerDF)
postprimerDF$V2<-rownames(postprimerDF)
colnames(preprimerDF)<-c("preprimer","ENST")
colnames(postprimerDF)<-c("postprimer","ENST")


#Add primers to the output
primerTable<-merge(x = preprimerDF,y=postprimerDF,by= "ENST")
exportTable4<-merge(x = exportTable3, y = primerTable, by = "ENST")


#Add DEG informations
DEG.PresortVPBodies <- read.delim("~/Documents/PBodies/DEG-PresortVPBodies.tsv")
ENSG_2 <- sub('\\.[0-9]*$', '', DEG.PresortVPBodies$ENSG)
DEG.PresortVPBodies$ENSG <- ENSG_2
exportTable5<-merge(x = exportTable4,y=DEG.PresortVPBodies,by.x = "hgnc_symbol",by.y = "GeneName",all=FALSE)

exportTable5$gene_len<-getGeneLengthAndGCContent(exportTable5$ENSG, "hsa")[,"length"]
gene_len<-exportTable5$gene_len
#exportTable5$gene_len<-gene_len
RPKM<-exportTable5[,c(11:16)]/exportTable5$gene_len
exportTable5[,c(11:16)]<-RPKM*1000

colnames(exportTable5)[c(11:16)]<-c("PB1","PB2","PB3","Cyto1","Cyto2","Cyto3")


###  Generate output  ###
write.table(exportTable5,paste0(args$output, "/Final_Last_Exon.tsv",collapse = ""),sep="\t",row.names = F,quote=F)
logging(paste0("Final tsv generated : ", args$output, "/Final_Last_Exon.tsv",collapse = ""))