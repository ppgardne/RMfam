#R CMD BATCH --no-save plotROC.R

rocFS      <-read.table("ROC_data_fracSeqs.dat",          header = T, sep = "\t")
rocSB      <-read.table("ROC_data_sumBits.dat",           header = T, sep = "\t")
rocWB      <-read.table("ROC_data_weightedSumBits.dat",   header = T, sep = "\t")
rocSSSB    <-read.table("ROC_data_single_seq_sumBits.dat",header = T, sep = "\t")

#METRIC	TP	TN	FP	FN	MCC	SENS	SPEC	PPV	NPV	ACC	FDR
#fracSeqs
#sumBits
#weightedSumBits

########

pdf(file="ROC_accuracy.pdf", width=20, height=20)
op<-par(mfrow=c(2,2),cex=3.0)
plot(rocFS$SPEC, rocFS$SENS, col="red", lwd=8, type="l",ylim=c(0,1),xlim=c(0,1),xlab="Specificity",ylab="Sensitivity",main="ROC plot")      
lines(rocSB$SPEC, rocSB$SENS, col="blue",lwd=8)
lines(rocWB$SPEC, rocWB$SENS, col="purple",lwd=8)
lines(rocSSSB$SPEC, rocSSSB$SENS, col="cyan",lwd=8)
points(rocFS$SPEC[rocFS$MCC==max(rocFS$MCC)], rocFS$SENS[rocFS$MCC==max(rocFS$MCC)], col="red",pch="x",lwd=12)
points(rocSB$SPEC[rocSB$MCC==max(rocSB$MCC)], rocSB$SENS[rocSB$MCC==max(rocSB$MCC)], col="blue",pch="x",lwd=12)
points(rocWB$SPEC[rocWB$MCC==max(rocWB$MCC)], rocWB$SENS[rocWB$MCC==max(rocWB$MCC)], col="purple",pch="x",lwd=12)
points(rocSSSB$SPEC[rocSSSB$MCC==max(rocSSSB$MCC)], rocSSSB$SENS[rocSSSB$MCC==max(rocSSSB$MCC)], col="cyan",pch="x",lwd=12)
#legend("bottomleft",c("fraction seqs.","sum bits", "weighted sum bits", "single-sequence sum bits"),col=c("red","blue","purple", "cyan"),lwd=8,cex=0.7,ncol=1)

plot(rocFS$ACC, rocFS$FDR, col="red", lwd=8, type="l",ylim=c(0.2,1),xlim=c(0.2,1),xlab="ACC",ylab="FDR",main="ROC-like plot")      
lines(rocSB$ACC, rocSB$FDR, col="blue",lwd=8)
lines(rocWB$ACC, rocWB$FDR, col="purple",lwd=8)
lines(rocSSSB$ACC, rocSSSB$FDR, col="cyan",lwd=8)
points(rocFS$ACC[rocFS$MCC==max(rocFS$MCC)], rocFS$FDR[rocFS$MCC==max(rocFS$MCC)], col="red",pch="x",lwd=12)
points(rocSB$ACC[rocSB$MCC==max(rocSB$MCC)], rocSB$FDR[rocSB$MCC==max(rocSB$MCC)], col="blue",pch="x",lwd=12)
points(rocWB$ACC[rocWB$MCC==max(rocWB$MCC)], rocWB$FDR[rocWB$MCC==max(rocWB$MCC)], col="purple",pch="x",lwd=12)
points(rocSSSB$ACC[rocSSSB$MCC==max(rocSSSB$MCC)], rocSSSB$FDR[rocSSSB$MCC==max(rocSSSB$MCC)], col="cyan",pch="x",lwd=12)
legend("bottomleft",c("fraction seqs.","sum bits", "weighted sum bits", "single-sequence sum bits"),col=c("red","blue","purple", "cyan"),lwd=8,cex=0.7,ncol=1)

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

pdf(file="ROC.pdf", width=20, height=20)
op<-par(mfrow=c(1,1),cex=5.0)
plot(rocSB$SPEC, rocSB$SENS, col="blue", lwd=12, type="l",ylim=c(0,1),xlim=c(0,1),xlab="Specificity",ylab="Sensitivity",main="ROC plot")      
#lines(rocSB$SPEC, rocSB$SENS, col="blue",lwd=12)
#lines(rocWB$SPEC, rocWB$SENS, col="purple",lwd=12)
lines(rocSSSB$SPEC, rocSSSB$SENS, col="red",lwd=8)
#points(rocFS$SPEC[rocFS$MCC==max(rocFS$MCC)], rocFS$SENS[rocFS$MCC==max(rocFS$MCC)], col="red",   pch="x",lwd=16)
points(rocSB$SPEC[rocSB$MCC==max(rocSB$MCC)], rocSB$SENS[rocSB$MCC==max(rocSB$MCC)], col="blue",  pch="x",lwd=16)
#points(rocWB$SPEC[rocWB$MCC==max(rocWB$MCC)], rocWB$SENS[rocWB$MCC==max(rocWB$MCC)], col="purple",pch="x",lwd=16)
points(rocSSSB$SPEC[rocSSSB$MCC==max(rocSSSB$MCC)], rocSSSB$SENS[rocSSSB$MCC==max(rocSSSB$MCC)], col="red",pch="x",lwd=12)
legend("bottomleft",c("alignment sum-bits", "sequence sum bits"),col=c("blue","red"),lwd=12,cex=0.7,ncol=1)
dev.off()

########

pdf(file="ROC_thresh.pdf", width=20, height=20)
op<-par(mfrow=c(1,1),cex=5.0)
#plot(rocSB$sumBits/max(rocSB$sumBits), rocSB$MCC,                 col="blue",lwd=12, type="l",ylim=c(0,max(rocSB$MCC)),xlim=c(0,1),xlab="threshold",ylab="MCC",main="MCC vs threshold")      
plot(rocSB$sumBits, rocSB$MCC,                 col="blue",lwd=12, type="l",ylim=c(0,max(rocSB$MCC)),xlim=c(0,1000),xlab="threshold (bits)",ylab="MCC",main="MCC vs threshold")      
#lines(rocSB$sumBits/max(rocSB$sumBits), rocSB$MCC,                 col="blue",lwd=12)
#lines(rocWB$weightedSumBits/max(rocWB$weightedSumBits), rocWB$MCC, col="purple",lwd=12)
#lines(rocSSSB$sngl_seq_sumBits/max(rocSSSB$sngl_seq_sumBits), rocSSSB$MCC, col="red",lwd=8)
lines(rocSSSB$sngl_seq_sumBits, rocSSSB$MCC, col="red",lwd=8)
legend("topright",c("alignment sum bits", "sequence sum bits"),col=c("blue","red"),lwd=12,cex=0.7,ncol=1)
dev.off()

