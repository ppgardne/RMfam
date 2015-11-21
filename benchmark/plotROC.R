#R CMD BATCH --no-save plotROC.R

######################################################################
#benchmark1: RMfam

rocFS      <-read.table("ROC_data_fracSeqs.dat",          header = T, sep = "\t")
rocSB      <-read.table("ROC_data_sumBits.dat",           header = T, sep = "\t")
rocWB      <-read.table("ROC_data_weightedSumBits.dat",   header = T, sep = "\t")
rocSSSB    <-read.table("ROC_data_single_seq_sumBits.dat",header = T, sep = "\t")
rocSSSBM   <-read.table("ROC_data_single_seq_sumBits-motifPerMotif.dat",header = T, sep = "\t")
rocB2ALL   <-read.table("ROC_benchmark2-all.dat",           header = T, sep = "\t")

#METRIC	TP	TN	FP	FN	MCC	SENS	SPEC	PPV	NPV	ACC	FDR
#fracSeqs
#sumBits
#weightedSumBits

########

pdf(file="figures/ROC_accuracy.pdf", width=20, height=20)
op<-par(mfrow=c(2,2),cex=3.0)
plot(rocFS$SPEC, rocFS$SENS, col="red", lwd=8, type="l",ylim=c(0,1),xlim=c(0,1),xlab="Specificity",ylab="Sensitivity",main="ROC plot")      
lines(rocSB$SPEC, rocSB$SENS, col="blue",lwd=8)
lines(rocWB$SPEC, rocWB$SENS, col="purple",lwd=8)
lines(rocSSSB$SPEC, rocSSSB$SENS, col="cyan",lwd=8)
points(rocFS$SPEC[rocFS$MCC==max(rocFS$MCC)], rocFS$SENS[rocFS$MCC==max(rocFS$MCC)], col="red",pch="x",lwd=12)
points(rocSB$SPEC[rocSB$MCC==max(rocSB$MCC)], rocSB$SENS[rocSB$MCC==max(rocSB$MCC)], col="blue",pch="x",lwd=12)
points(rocWB$SPEC[rocWB$MCC==max(rocWB$MCC)], rocWB$SENS[rocWB$MCC==max(rocWB$MCC)], col="purple",pch="x",lwd=12)
points(rocSSSB$SPEC[rocSSSB$MCC==max(rocSSSB$MCC)], rocSSSB$SENS[rocSSSB$MCC==max(rocSSSB$MCC)], col="cyan",pch="x",lwd=12)
legend("bottomleft",c("fraction seqs.","sum bits", "weighted sum bits", "single-sequence sum bits"),col=c("red","blue","purple", "cyan"),lwd=8,cex=0.7,ncol=1)

plot(rocFS$ACC, rocFS$FDR, col="red", lwd=8, type="l",ylim=c(0.2,1),xlim=c(0.2,1),xlab="ACC",ylab="FDR",main="ROC-like plot")      
lines(rocSB$ACC, rocSB$FDR, col="blue",lwd=8)
lines(rocWB$ACC, rocWB$FDR, col="purple",lwd=8)
lines(rocSSSB$ACC, rocSSSB$FDR, col="cyan",lwd=8)
points(rocFS$ACC[rocFS$MCC==max(rocFS$MCC)], rocFS$FDR[rocFS$MCC==max(rocFS$MCC)], col="red",pch="x",lwd=12)
points(rocSB$ACC[rocSB$MCC==max(rocSB$MCC)], rocSB$FDR[rocSB$MCC==max(rocSB$MCC)], col="blue",pch="x",lwd=12)
points(rocWB$ACC[rocWB$MCC==max(rocWB$MCC)], rocWB$FDR[rocWB$MCC==max(rocWB$MCC)], col="purple",pch="x",lwd=12)
points(rocSSSB$ACC[rocSSSB$MCC==max(rocSSSB$MCC)], rocSSSB$FDR[rocSSSB$MCC==max(rocSSSB$MCC)], col="cyan",pch="x",lwd=12)

plot(rocFS$fracSeqs/max(rocFS$fracSeqs), rocFS$MCC, col="red", lwd=8, type="l",ylim=c(0,max(rocSB$MCC)),xlim=c(0,1),xlab="threshold",ylab="MCC",main="MCC vs threshold")      
lines(rocSB$sumBits/max(rocSB$sumBits), rocSB$MCC, col="blue",lwd=8)
lines(rocWB$weightedSumBits/max(rocWB$weightedSumBits), rocWB$MCC, col="purple",lwd=8)
lines(rocSSSB$sngl_seq_sumBits/max(rocSSSB$sngl_seq_sumBits), rocSSSB$MCC, col="cyan",lwd=8)
#legend("topright",c("fraction seqs.","sum bits", "weighted sum bits", "single-sequence sum bits"),col=c("red","blue","purple","cyan"),lwd=8,cex=0.7,ncol=1)

plot(rocFS$fracSeqs/max(rocFS$fracSeqs), rocFS$FDR, col="red", lwd=8, type="l",ylim=c(0,1),xlim=c(0,1),xlab="threshold",ylab="FDR",main="FDR plot")      
lines(rocSB$sumBits/max(rocSB$sumBits), rocSB$FDR, col="blue",lwd=8)
lines(rocWB$weightedSumBits/max(rocWB$weightedSumBits), rocWB$FDR, col="purple",lwd=8)
lines(rocSSSB$sngl_seq_sumBits/max(rocSSSB$sngl_seq_sumBits), rocSSSB$FDR, col="cyan",lwd=8)
#legend("topright",c("fraction seqs.","sum bits", "weighted sum bits", "single-sequence sum bits"),col=c("red","blue","purple","cyan"),lwd=8,cex=0.7,ncol=1)
dev.off()

########

pdf(file="figures/ROC.pdf", width=20, height=20)
op<-par(mfrow=c(1,1),cex=5.0)
plot(rocSB$SPEC, rocSB$SENS, col="blue", lwd=12, type="l",ylim=c(0,1),xlim=c(0,1),xlab="Specificity",ylab="Sensitivity",main="ROC plot")      
#lines(rocSB$SPEC, rocSB$SENS, col="blue",lwd=12)
#lines(rocWB$SPEC, rocWB$SENS, col="purple",lwd=12)
lines(rocSSSB$SPEC, rocSSSB$SENS, col="red",lwd=8)
lines(rocB2ALL$SPEC, rocB2ALL$SENS, col="rosybrown",lwd=8)
#points(rocFS$SPEC[rocFS$MCC==max(rocFS$MCC)], rocFS$SENS[rocFS$MCC==max(rocFS$MCC)], col="red",   pch="x",lwd=16)
points(rocSB$SPEC[rocSB$MCC==max(rocSB$MCC)], rocSB$SENS[rocSB$MCC==max(rocSB$MCC)], col="blue",  pch="x",lwd=16)
#points(rocWB$SPEC[rocWB$MCC==max(rocWB$MCC)], rocWB$SENS[rocWB$MCC==max(rocWB$MCC)], col="purple",pch="x",lwd=16)
points(rocSSSB$SPEC[rocSSSB$MCC==max(rocSSSB$MCC)], rocSSSB$SENS[rocSSSB$MCC==max(rocSSSB$MCC)], col="red",pch="x",lwd=12)
points(rocB2ALL$SPEC[rocB2ALL$MCC==max(rocB2ALL$MCC)], rocB2ALL$SENS[rocB2ALL$MCC==max(rocB2ALL$MCC)], col="rosybrown",pch="x",lwd=12)
lines(c(0.93,0.6),c(0.76,0.6), col="blue",lwd=4)
text(0.6,0.6,'sum-bits=145\nmaxMCC=0.68',pos=2,col="blue",cex=0.5)
lines(c(1,0.8),c(0.48,0.4), col="rosybrown",lwd=4)
text(0.8,0.4,'bits=18.5\nmaxMCC=0.66',pos=2,col="rosybrown",cex=0.5)
text(0.99,0.22,'bits=19\nmaxMCC=0.38',pos=2,col="red",cex=0.5)
legend("bottomleft",c("Rfam alignments (sum-bits)", "Rfam sequences (bits)", "RMfam sequences (bits)"),col=c("blue","red","rosybrown"),lwd=12,cex=0.7,ncol=1)
dev.off()

########

pdf(file="figures/ROC_thresh.pdf", width=20, height=20)
op<-par(mfrow=c(1,1),cex=5.0)
#plot(rocSB$sumBits/max(rocSB$sumBits), rocSB$MCC,                 col="blue",lwd=12, type="l",ylim=c(0,max(rocSB$MCC)),xlim=c(0,1),xlab="threshold",ylab="MCC",main="MCC vs threshold")      
plot(rocSB$sumBits, rocSB$MCC,                 col="blue",lwd=12, type="l",ylim=c(0,max(rocSB$MCC)),xlim=c(0,1000),xlab="threshold (bits)",ylab="MCC",main="MCC vs threshold")      
#lines(rocSB$sumBits/max(rocSB$sumBits), rocSB$MCC,                 col="blue",lwd=12)
#lines(rocWB$weightedSumBits/max(rocWB$weightedSumBits), rocWB$MCC, col="purple",lwd=12)
#lines(rocSSSB$sngl_seq_sumBits/max(rocSSSB$sngl_seq_sumBits), rocSSSB$MCC, col="red",lwd=8)
lines(rocSSSB$sngl_seq_sumBits, rocSSSB$MCC, col="red",lwd=12)
legend("bottomright",c("Rfam alignment (sum bits)", "Rfam sequences (bits)"),col=c("blue","red"),lwd=12,cex=0.7,ncol=1)
dev.off()

######################################################################
#benchmark2: individual motifs
#cat ROC_benchmark2.dat | cut -f 13 | sort -d | uniq | awk '{print "\42"$1"\42,"}'
motif<-c(
"ANYA",
"AUF1_binding",
"C-loop",
"CRC_binding",
"CsrA_binding",
"CUYG",
"docking_elbow",
"Domain-V",
"GNRA",
"HuR_binding",
"k-turn-1",
"k-turn-2",
"pK-turn",
"RBS_B_subtilis",
"RBS_E_coli",
"RBS_H_pylori",
"right_angle-2",
"right_angle-3",
"Roquin_binding",
"sarcin-ricin-1",
"sarcin-ricin-2",
"SRP_S_domain",
"tandem-GA",
"Terminator1",
"Terminator2",
"T-loop",
"TRIT",
"twist_up",
"UAA_GAN",
"UMAC",
"UNCG",
"U-turn",
"vapC_target",
"VTS1_binding"
);
rocB2      <-read.table("ROC_benchmark2.dat",              header = T, sep = "\t")
rocB2A     <-read.table("ROC_data_sumBits_benchmark2.dat", header = T, sep = "\t")
thresholds <-read.table("thresholds.dat",                  header = T, sep = "\t")

pdf(file="figures/ROC_benchmark2.pdf", width=20, height=20)
op<-par(mfrow=c(1,1),cex=5.0)
plot(c(1,0), c(0,1), col="black",lwd=12, type="l",ylim=c(0,1),xlim=c(0,1),xlab="Specificity",ylab="Sensitivity",main="ROC curve: RMfam CMs to seqs")      
for(i in 1:length(motif)){
      lines(rocB2$SPEC[rocB2$motif==motif[i]], rocB2$SENS[rocB2$motif==motif[i]], lty=i, col=i, lwd=12)
}
legend("bottomleft", motif, col=1:length(motif), lty=1:length(motif), cex=0.3, ncol=3)
dev.off()



emptyColnames<-motif
maxMCC<-mat.or.vec(3,length(motif))
sumMaxMCC<-mat.or.vec(1,length(motif))
for(i in 1:length(motif)){

      maxMCC[1,i]<-max(rocB2$MCC[rocB2$motif==motif[i]])

      if(max(rocSSSBM$MCC[rocSSSBM$motif==motif[i]])>0){
	maxMCC[2,i]<-max(rocSSSBM$MCC[rocSSSBM$motif==motif[i]])
      }
      else{
	maxMCC[2,i]<-0
      }

      if(max(rocB2A$MCC[rocB2A$motif==motif[i]])>0){
	maxMCC[3,i]<-max(rocB2A$MCC[rocB2A$motif==motif[i]])
      }
      else{
	maxMCC[3,i]<-0
      }
      
      sumMaxMCC[i]<-maxMCC[1,i]+maxMCC[2,i]+maxMCC[3,i]
      emptyColnames[i]<-"";
}
maxMCCs<-sort(sumMaxMCC,decreasing=T,index.return=T)

maxMCC<-array(data = maxMCC[,maxMCCs$ix], dim = c(3,length(motif)))
colnames(maxMCC)<-motif[maxMCCs$ix]	  
motifR<-motif[maxMCCs$ix]	  

pdf(file="figures/maxMCC_benchmark2.pdf", width=37, height=20)
op<-par(mfrow=c(1,2),cex=3.0,las=1)#bottom, left, top, right
plot(rocSB$SPEC, rocSB$SENS, col="blue", lwd=12, type="l",ylim=c(0,1),xlim=c(0,1),xlab="Specificity",ylab="Sensitivity",main="ROC plots for 3 RMfam benchmarks")      
lines(rocSSSB$SPEC, rocSSSB$SENS, col="red",lwd=8)
lines(rocB2ALL$SPEC, rocB2ALL$SENS, col="rosybrown",lwd=8)
points(rocSB$SPEC[rocSB$MCC==max(rocSB$MCC)], rocSB$SENS[rocSB$MCC==max(rocSB$MCC)], col="blue",  pch="x",lwd=16)
points(rocSSSB$SPEC[rocSSSB$MCC==max(rocSSSB$MCC)], rocSSSB$SENS[rocSSSB$MCC==max(rocSSSB$MCC)], col="red",pch="x",lwd=12)
points(rocB2ALL$SPEC[rocB2ALL$MCC==max(rocB2ALL$MCC)], rocB2ALL$SENS[rocB2ALL$MCC==max(rocB2ALL$MCC)], col="rosybrown",pch="x",lwd=12)
lines(c(0.93,0.6),c(0.76,0.6), col="blue",lwd=4)
text(0.6,0.6,'sum-bits=145\nmaxMCC=0.68',pos=2,col="blue",cex=1.0)
lines(c(1,0.8),c(0.48,0.4), col="rosybrown",lwd=4)
text(0.8,0.4,'bits=18.5\nmaxMCC=0.66',pos=2,col="rosybrown",cex=1.0)
lines(c(0.9,0.99),c(0.22,0.22), col="red",lwd=4)
text(0.9,0.22,'bits=19\nmaxMCC=0.38',pos=2,col="red",cex=1.0)
legend(0.025,0.2,c("Rfam alignments (sum-bits)", "Rfam sequences (bits)", "RMfam sequences (bits)"),fill=c("blue","red","rosybrown"),cex=1.0,ncol=1)
op<-par(mar=c(5, 4+2, 4, 2+1)+0.1)
barplot(maxMCC, ylab="", xlim=c(0,1), xlab="max(MCC)", main="max(MCC) values for RMfam motifs",col=c("rosybrown","red","blue"),horiz = TRUE, cex.names=0.85, beside=TRUE )
dev.off()


pdf(file="figures/maxMCC_benchmark2-2.pdf", width=40, height=20)
op<-par(mfrow=c(1,3),cex=3.0,las=1,mar=c(5, 4+4, 4, 2+1)+0.1)#bottom, left, top, right
barplot(maxMCC[1,], ylab="", xlim=c(0,1), xlab="max(MCC)", main="RMfam sequences (bits)",col=c("rosybrown"),horiz = TRUE, axisnames = FALSE)
axis(2, 1.2*(1:length(motifR))-0.5, motifR,cex=0.75,tick=FALSE)
op<-par(mar=c(5, 0+0, 4, 2+1)+0.1)#bottom, left, top, right
barplot(maxMCC[2,], ylab="", xlim=c(0,1), xlab="max(MCC)", main="Rfam sequences (bits)",col=c("red"), horiz = TRUE, axisnames = FALSE )
legend("topright",c("Rfam alignments (sum-bits)", "Rfam sequences (bits)", "RMfam sequences (bits)"),fill=c("blue","red","rosybrown"),ncol=1) #,cex=1.2
barplot(maxMCC[3,], ylab="", xlim=c(0,1), xlab="max(MCC)", main="Rfam alignments (sum-bits)",col=c("blue"),horiz = TRUE, axisnames = FALSE )
dev.off()

######################################################################
######################################################################
######################################################################
######################################################################

for(i in 1:length(motif)){

pdfName<-sprintf("figures/%s_benchmark.pdf",   motif[i])
pdf(file=pdfName, width=20, height=20)
op<-par(mfrow=c(2,2),cex=3.0)

main<-sprintf("ROC curve: %s",   motif[i])
plot(c(1,0), c(0,1), col="black",lwd=6, type="l",ylim=c(0,1),xlim=c(0,1),xlab="Specificity",ylab="Sensitivity",main=main)
lines(rocB2$SPEC[rocB2$motif==motif[i]], rocB2$SENS[rocB2$motif==motif[i]], lty=1, col="rosybrown", lwd=6)
lines(rocB2A$SPEC[rocB2A$motif==motif[i]], rocB2A$SENS[rocB2A$motif==motif[i]], lty=1, col="blue", lwd=6)
lines(rocSSSBM$SPEC[rocSSSBM$motif==motif[i]], rocSSSBM$SENS[rocSSSBM$motif==motif[i]], lty=1, col="red", lwd=6)
legend("bottomleft",c("Rfam alignments (sum-bits)", "Rfam sequences (bits)", "RMfam sequences (bits)"),col=c("blue","red","rosybrown"),lwd=12,cex=0.7,ncol=1)

main<-sprintf("ROC-like curve: %s",   motif[i])
plot(c(1,0), c(0,1), col="black",lwd=6, type="l",ylim=c(0,1),xlim=c(0,1),xlab="PPV",ylab="Sensitivity",main=main)
lines(rocB2$PPV[rocB2$motif==motif[i]], rocB2$SENS[rocB2$motif==motif[i]], lty=1, col="rosybrown", lwd=6)
lines(rocB2A$PPV[rocB2A$motif==motif[i]], rocB2A$SENS[rocB2A$motif==motif[i]], lty=1, col="blue", lwd=6)
lines(rocSSSBM$PPV[rocSSSBM$motif==motif[i]], rocSSSBM$SENS[rocSSSBM$motif==motif[i]], lty=1, col="red", lwd=6)
legend("bottomleft",c("Rfam alignments (sum-bits)", "Rfam sequences (bits)", "RMfam sequences (bits)"),col=c("blue","red","rosybrown"),lwd=12,cex=0.7,ncol=1)


###
thresh<-thresholds$GA[thresholds$NAME==motif[i]]
main<-sprintf("Score vs MCC: %s",   motif[i])
op<-par(las=2)
plot(rocB2$bits[rocB2$motif==motif[i]], rocB2$MCC[rocB2$motif==motif[i]], lty=1, col="rosybrown", lwd=6, type="l",ylim=c(-1,1),xlim=c(0,50),xlab="CM score (bits)",ylab="MCC",main=main)
lines(rocB2A$bits[rocB2A$motif==motif[i]], rocB2A$MCC[rocB2A$motif==motif[i]], lty=1, col="blue", lwd=6, type="l")
lines(rocSSSBM$sngl_seq_sumBits[rocSSSBM$motif==motif[i]], rocSSSBM$MCC[rocSSSBM$motif==motif[i]], lty=1, col="red", lwd=6, type="l")
lines(c(0,1000), c(0,0), lty=1, col="black",)
lines(c(thresh,thresh), c(-1,1), lty=1, col="black")

###
maxMCCme<-rocB2[rocB2$MCC==maxMCC[1,colnames(maxMCC)==motif[i]] & rocB2$motif==motif[i],]
maxMCCme<-maxMCCme[maxMCCme$TN==max(maxMCCme$TN),]
maxMCCme<-maxMCCme[maxMCCme$FN==min(maxMCCme$FN),6:12]

maxMCCSme<-rocSSSBM[rocSSSBM$MCC==maxMCC[2,colnames(maxMCC)==motif[i]] & rocSSSBM$motif==motif[i], ]
if(length(row.names(maxMCCSme))==0){
	maxMCCSme<-as.data.frame(array(data = 0*1:13, dim = c(1,13)))
	colnames(maxMCCSme)<-colnames(rocSSSBM)
}
maxMCCSme<-maxMCCSme[maxMCCSme$TN==max(maxMCCSme$TN),]
maxMCCSme<-maxMCCSme[maxMCCSme$FN==min(maxMCCSme$FN),6:12]

maxMCCAme<-rocB2A[rocB2A$MCC==maxMCC[3,colnames(maxMCC)==motif[i]] & rocB2A$motif==motif[i], ]
if(length(row.names(maxMCCAme))==0){
	maxMCCAme<-as.data.frame(array(data = 0*1:13, dim = c(1,13)))
	colnames(maxMCCAme)<-colnames(rocB2A)
}
maxMCCAme<-maxMCCAme[maxMCCAme$TN==max(maxMCCAme$TN),]
maxMCCAme<-maxMCCAme[maxMCCAme$FN==min(maxMCCAme$FN),6:12]

stats<-array(data = as.matrix(c(maxMCCme, maxMCCSme, maxMCCAme)), dim = c(7,3)); 

rownames(stats)<-colnames(rocB2[1,6:12])

main<-sprintf("Performance: %s",   motif[i])
op<-par(las=2)
barplot(t(stats), xlab="", ylab="Stats", ylim=c(0,3), main=main,col=c("rosybrown","red","blue"))

dev.off()

}
