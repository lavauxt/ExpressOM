#' Fork to SPIA Plot function
#' 
#' This function generates a two-way evidence plot for SPIA results, showing the relationship between the p-values for NDE (Number of Differentially Expressed genes) and PERT (Perturbation). It also highlights significant pathways based on the specified threshold and indicates the combined p-value method used (Fisher or Normal Inverse). The plot includes lines representing the significance thresholds and labels for significant pathways.
#' @export
#' @param x A SPIA results data.frame
#' @param threshold Threshold for significance
#' @return NULL (generates a plot)
plotP_fork<-function(x,threshold=0.01){

if(!inherits(x, "data.frame") | dim(x)[1]<1 | !all(c("ID","pNDE","pPERT","pG","pGFdr","pGFWER")%in%names(x)))
{
 stop("plotP can be applied only to a dataframe produced by spia function!!!") 
}

if(threshold<x[1,"pGFdr"]){
print(paste("The threshold value was corrected to be equal to ",x[1,"pGFdr"]))
threshold <- x[1,"pGFdr"]
}

pb<-x[,"pPERT"]
ph<-x[,"pNDE"]

# FIX: Added SPIA::: to access internal 'combfunc'
combinemethod=ifelse(sum(SPIA:::combfunc(pb,ph,"fisher")==x$pG)>sum(SPIA:::combfunc(pb,ph,"norminv")==x$pG),"fisher","norminv")

okx<-(ph<1e-6)
oky<-(pb<1e-6)
ph[ph<1e-6]<-1e-6
pb[pb<1e-6]<-1e-6

plot(-log(ph),-log(pb),xlim=c(0,max(c(-log(ph),-log(pb))+1,na.rm=TRUE)),
 ylim=c(0,max(c(-log(ph),-log(pb)+1),na.rm=TRUE)),pch=19,main="SPIA two-way evidence plot",cex=1.5,
 xlab="-log(P NDE)",ylab="-log(P PERT)")
tr<-threshold/dim(na.omit(x))[1]

if(combinemethod=="fisher"){
# FIX: Added SPIA::: to access internal 'getP2'
points(c(0,-log(SPIA:::getP2(tr,"fisher")^2)),c(-log(SPIA:::getP2(tr,"fisher")^2),0),col="red",lwd=2,cex=0.7,type="l")
}else{
somep1=exp(seq(from=min(log(ph)),to=max(log(ph)),length=200))
somep2=pnorm(qnorm(tr)*sqrt(2)-qnorm(somep1))
points(-log(somep1),-log(somep2),col="red",lwd=2,cex=0.7,type="l") 
}

trold=tr
tr<-max(x[,"pG"][x[,"pGFdr"]<=threshold])
if(tr<=trold){tr=trold*1.03}

if(combinemethod=="fisher"){
# FIX: Added SPIA::: to access internal 'getP2'
points(c(0,-log(SPIA:::getP2(tr,"fisher")^2)),c(-log(SPIA:::getP2(tr,"fisher")^2),0),col="blue",lwd=2,cex=0.7,type="l")
 }else{
somep1=exp(seq(from=min(log(ph)),to=max(log(ph)),length=200))
somep2=pnorm(qnorm(tr)*sqrt(2)-qnorm(somep1))
points(-log(somep1),-log(somep2),col="blue",lwd=2,cex=0.7,type="l") 
}

oks<-x[,"pGFWER"]<=threshold
oks2<-x[,"pGFdr"]<=threshold

points(-log(ph)[oks2],-log(pb)[oks2],pch=19,col="blue",cex=1.5)
points(-log(ph)[oks],-log(pb)[oks],pch=19,col="red",cex=1.5)

if(sum(oks)>0){
 text(-log(ph)[oks]+0.70,-log(pb)[oks],labels=as.vector(x$ID)[oks],cex=0.65)
}
if(sum(oks2)>0){
 text(-log(ph)[oks2]+0.70,-log(pb)[oks2],labels=as.vector(x$ID)[oks2],cex=0.65)
}

testx <- sum(okx)
if (!(is.na(testx)) && (sum(okx)>0)) {
    points(-log(ph)[okx]-0.12,-log(pb)[okx],pch="|",col="black",cex=1.5)
}

testy <- sum(oky)
if (!(is.na(testy)) && (sum(oky)>0)) {
    points(-log(ph)[oky],-log(pb)[oky]-0.12,pch="_",col="black",cex=1.5)
}

}