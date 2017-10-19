##############Altemose et al. eLife 2017 R plotting code
###############################################################################
#executing these commands in R will reproduce the plots in Figures 1,2,3, and 5
#code was tested in R 3.3.1 [Mac GUI 1.68 Mavericks build (7238)]
#ensure the following libraries are installed: gtools, seqLogo, boot
#download and decompress all files in the zipped directories 
#    "Altemose2017_plotting_files_1of2.zip" and 
#    "Altemose2017_plotting_files_2of2.zip"
#    and move them to the active R directory
###############################################################################


###############################################
### Figure 1a and Figure 1-Figure supplement 3
###############################################

#see FindMotifs_Human.R for code used to generate motif PWMs

rm(list=ls())
#Get PWMs from all motifs, stored in intermediate file "Human_Motif_Results_Final_iter240.r"
load("Human_Motif_Results_Final_iter240.r")
library(gtools)
library(seqLogo)
mySeqLogo = seqLogo::seqLogo 
bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") | sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
body(mySeqLogo)[bad] = NULL 
motnum=length(znew$scorematdim)
plotting=TRUE
dimvec = znew$scorematdim
scorematset = znew$scoremat
alpha= znew$alpha
iter=240
compvec=c(0,1,0,0,0,0,1,1,0,0,1,1,0,0,0,1,1) #specify whether to take reverse complement of each motif for final output (0=no, 1=yes)
keep=c(4,6,8,9,10,15,17,5,11,13,7,16,12,3,2,14,1) #specify the final ordering of motif names

#plot motifs and output individual PWM text files for each motif
system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
props = rep(0,motnum)
for(i in 1:motnum){
	newstarts=c(1,cumsum(dimvec)+1)
	newends=cumsum(dimvec)
	comp=compvec[i]

	testmat=scorematset[newstarts[i]:newends[i],]
	testmat=exp(testmat)
	compmat=testmat[,c(4:1)]
	testmat=testmat/rowSums(testmat)
	compmat=compmat/rowSums(compmat)
	pwm = t(testmat)
	compmat=compmat[nrow(compmat):1,]
	grid.newpage()
	name=which(keep==i)
	if(comp==1){
		mySeqLogo(t(compmat))
		#write.table(compmat,file=paste("plots/PWM.Iter",iter,".Motif",name,".tbl",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
	}else{
		mySeqLogo(pwm)
		#write.table(t(pwm),file=paste("plots/PWM.Iter",iter,".Motif",name,".tbl",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
	}
	dev.copy2pdf(file=paste("plots/LOGO_Human_Motif",name,".pdf",sep=""),width=36,height=24)
	
	hist(znew$whichpos[znew$whichmot==i],breaks=seq(0,300,25),col='blue',xlab="",ylab="",main="")
	dev.copy2pdf(file=paste("plots/Histogram_Human_Motif",name,".pdf",sep=""),width=36,height=24)
	props[name]=mean(znew$whichpos[znew$whichmot==i]>100 & znew$whichpos[znew$whichmot==i]<200)
}
print(props)
# [1] 0.8961039 0.9119718 0.8520710 0.8634686 0.8613139 0.8543860 0.8623188
# [8] 0.7878788 0.7822581 0.7577320 0.7631579 0.7486486 0.7391304 0.6950820
#[15] 0.6071429 0.4716981 0.2398374


#############
### Figure 1b
#############

rm(list=ls())
dataH0=read.table("ForceCall.Pratto_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataH0)
dataH=subset(dataH0,covg>5 & covg<quantile(dataH0$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500))
dim(dataH)

dataD=read.table("ForceCall.Pratto_DMC1_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataD)
dataD=subset(dataD,covg>5 & covg<quantile(dataD$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500))
dim(dataD)

dataHH=read.table("ForceCall.YFP_Human_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataHH)
dataHH=subset(dataHH,covg>5 & covg<quantile(dataHH$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500))
dim(dataHH)

input1=subset(dataH,dataH$AnyH3K4me3Overlap==0)
input2=subset(dataHH,dataH$AnyH3K4me3Overlap==0)
input3=subset(dataD,dataH$AnyH3K4me3Overlap==0)

#123662 peaks informative
#78333 have covg in range
#37118 do not overlap PRDM9-independent H3K4me3

stderrprop <- function(x,n) sqrt(x*(1-x)/n)

bins1=quantile(input1$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1

results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub1 = subset(input1, enrich >=bins1[i] & enrich <bins1[i+1])
	results1[i] = mean(sub1$H.pval<=pthresh,na.rm=TRUE)
	seprop1 = stderrprop(results1[i],dim(sub1)[1])
	lower1[i] = results1[i] - 2*seprop1
	upper1[i] = results1[i] + 2*seprop1
	midpoints1[i] = median(sub1$enrich,na.rm=TRUE)

	sub2 = subset(input2, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.pval<=pthresh,na.rm=TRUE)
	seprop2 = stderrprop(results2[i],dim(sub2)[1])
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	
	sub3 = subset(input3, enrich >=bins1[i] & enrich <bins1[i+1])
	results3[i] = mean(sub3$H.pval<=pthresh,na.rm=TRUE)
	seprop3 = stderrprop(results3[i],dim(sub3)[1])
	lower3[i] = results3[i] - 2*seprop3
	upper3[i] = results3[i] + 2*seprop3
}


plot(midpoints1,results1,,col='red4',ylim=c(0.1,1),type="l",lwd=2,xlab="",ylab="")
points(midpoints1,results1,col='red4',pch=16)
points(midpoints1,results2,col='red',pch=16)
points(midpoints1,results3,col='orange',pch=16)
lines(midpoints1,results2,col='red',lwd=2)
lines(midpoints1,results3,col='orange',lwd=2)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1b.pdf",sep=""),width=5,height=5)


#############
### Figure 1c
#############

rm(list=ls())
###code used to generate intermediate files
#dataH=read.table("ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt",header=T)
#dim(dataH)
##170198 26
#dataH=subset(dataH,cov_input<quantile(dataH$cov_input,0.999) & anyENCODE_H3K4me3_overlap==0 & cov_input>5 & is.na(which_motif)==F & is.na(FIMO_score)==F & promoter_overlap==0)
#dim(dataH)
##21963 26
#bins1=quantile(dataH$enrichment,seq(0,1,0.25),na.rm=TRUE)
#bins1[length(bins1)]=bins1[length(bins1)]+1
#for(i in 1:(length(bins1)-1)){
#	sub2 = subset(dataH, enrichment >=bins1[i] & enrichment <bins1[i+1])
#	write.table(cbind(sub2[,1:3],rep("+",dim(sub2)[1])),file=paste("ChIPseq_Peak_Positions.YFP_HumanPRDM9.PRDM9quartile",i,".bed",sep=""),row.names=F,col.names=F,sep="\t",quote=F)	
#}
#system("perl MakeProfilePlot.pl ChIPseq_Peak_Positions.YFP_HumanPRDM9.PRDM9quartile1.bed HumanRecMap.HapMapCEU.CMperMb.bed ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile1.bed 10000 1 8 1")
#system("perl MakeProfilePlot.pl ChIPseq_Peak_Positions.YFP_HumanPRDM9.PRDM9quartile2.bed HumanRecMap.HapMapCEU.CMperMb.bed ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile2.bed 10000 1 8 1")
#system("perl MakeProfilePlot.pl ChIPseq_Peak_Positions.YFP_HumanPRDM9.PRDM9quartile3.bed HumanRecMap.HapMapCEU.CMperMb.bed ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile3.bed 10000 1 8 1")
#system("perl MakeProfilePlot.pl ChIPseq_Peak_Positions.YFP_HumanPRDM9.PRDM9quartile4.bed HumanRecMap.HapMapCEU.CMperMb.bed ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile4.bed 10000 1 8 1")

###plotting code given intermediate files
data1=read.table("ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile1.bed",header=T)
data2=read.table("ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile2.bed",header=T)
data3=read.table("ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile3.bed",header=T)
data4=read.table("ProfilePlot.RecRate.YFP_HumanPRDM9PeakCentres.PRDM9quartile4.bed",header=T)
plot(ksmooth(data4$Position,data4$Mean, "normal", bandwidth=25),type='l',col=colors()[497],ylim=c(1.8,11),xlab="",ylab="",main="",lwd=2)
lines(ksmooth(data3$Position,data3$Mean, "normal", bandwidth=25),type='l',col=colors()[496],lwd=2)
lines(ksmooth(data2$Position,data2$Mean, "normal", bandwidth=25),type='l',col=colors()[495],lwd=2)
lines(ksmooth(data1$Position,data1$Mean, "normal", bandwidth=25),type='l',col=colors()[494],lwd=2)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1c.pdf",sep=""),width=5,height=5)


#############
### Figure 1d
#############

#left plot
rm(list=ls())
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
motlens = c(32,31,32,50,24,23,31)

data=read.table("ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt",header=T)
dataS=subset(data, cov_input>5 & cov_input<quantile(cov_input,0.999) & is.na(HapMap_recombination_rate)==F)
dataSprom=subset(dataS, promoter_overlap==1 & DNaseHS_overlap==1 & HEK293T_H3K4me3_overlap==1 )
dataSnoprom = subset(dataS, promoter_overlap==0  & anyENCODE_H3K4me3_overlap==0 ) 

mot1noprom = subset(dataSnoprom,which_motif==1 )
mot2noprom = subset(dataSnoprom,which_motif==2 | which_motif==3 | which_motif==5)
mot3noprom = subset(dataSnoprom,which_motif==3 )
mot4noprom = subset(dataSnoprom,which_motif==4 )
mot5noprom = subset(dataSnoprom,which_motif==5 )
mot6noprom = subset(dataSnoprom,which_motif==6 )
mot7noprom = subset(dataSnoprom,which_motif==7 )

input1=mot1noprom
input2=mot2noprom
input3=mot3noprom
input4=mot4noprom
input5=mot5noprom
input6=mot6noprom
input7=mot7noprom

bins1=quantile(input1$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins2=quantile(input2$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins3=quantile(input3$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins4=quantile(input4$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins5=quantile(input5$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins6=quantile(input6$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins7=quantile(input7$enrichment,seq(0,1,0.25),na.rm=TRUE)

bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
bins3[length(bins3)]=bins3[length(bins3)]+1
bins4[length(bins4)]=bins4[length(bins4)]+1
bins5[length(bins5)]=bins5[length(bins5)]+1
bins6[length(bins6)]=bins6[length(bins6)]+1
bins7[length(bins7)]=bins7[length(bins7)]+1

results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins2)-1)
results3 = rep(0,length(bins3)-1)
results4 = rep(0,length(bins4)-1)
results5 = rep(0,length(bins5)-1)
results6 = rep(0,length(bins6)-1)
results7 = rep(0,length(bins7)-1)

midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins2)-1)
midpoints3 = rep(0,length(bins3)-1)
midpoints4 = rep(0,length(bins4)-1)
midpoints5 = rep(0,length(bins5)-1)
midpoints6 = rep(0,length(bins6)-1)
midpoints7 = rep(0,length(bins7)-1)

lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins2)-1)
upper2 = rep(NA,length(bins2)-1)
lower3 = rep(NA,length(bins3)-1)
upper3 = rep(NA,length(bins3)-1)
lower4 = rep(NA,length(bins4)-1)
upper4 = rep(NA,length(bins4)-1)
lower5 = rep(NA,length(bins5)-1)
upper5 = rep(NA,length(bins5)-1)
lower6 = rep(NA,length(bins6)-1)
upper6 = rep(NA,length(bins6)-1)
lower7 = rep(NA,length(bins7)-1)
upper7 = rep(NA,length(bins7)-1)

for(i in 1:(length(bins1)-1)){
	sub1 = subset(input1, enrichment >=bins1[i] & enrichment <bins1[i+1])
	results1[i] = mean(sub1$HapMap_recombination_rate,na.rm=TRUE)
	seprop1 = stderrmean(sub1$HapMap_recombination_rate)
	lower1[i] = results1[i] - 2*seprop1
	upper1[i] = results1[i] + 2*seprop1
	midpoints1[i] = median(sub1$enrichment,na.rm=TRUE)
	
	sub2 = subset(input2, enrichment >=bins2[i] & enrichment <bins2[i+1])
	results2[i] = mean(sub2$HapMap_recombination_rate,na.rm=TRUE)
	seprop2 = stderrmean(sub2$HapMap_recombination_rate)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints2[i] = median(sub2$enrichment,na.rm=TRUE)
	
	sub3 = subset(input3, enrichment >=bins3[i] & enrichment <bins3[i+1])
	results3[i] = mean(sub3$HapMap_recombination_rate,na.rm=TRUE)
	seprop3 = stderrmean(sub3$HapMap_recombination_rate)
	lower3[i] = results3[i] - 2*seprop3
	upper3[i] = results3[i] + 2*seprop3
	midpoints3[i] = median(sub3$enrichment,na.rm=TRUE)
	
	sub4 = subset(input4, enrichment >=bins4[i] & enrichment <bins4[i+1])
	results4[i] = mean(sub4$HapMap_recombination_rate,na.rm=TRUE)
	seprop4 = stderrmean(sub4$HapMap_recombination_rate)
	lower4[i] = results4[i] - 2*seprop4
	upper4[i] = results4[i] + 2*seprop4
	midpoints4[i] = median(sub4$enrichment,na.rm=TRUE)
	
	sub5 = subset(input5, enrichment >=bins5[i] & enrichment <bins5[i+1])
	results5[i] = mean(sub5$HapMap_recombination_rate,na.rm=TRUE)
	seprop5 = stderrmean(sub5$HapMap_recombination_rate)
	lower5[i] = results5[i] - 2*seprop5
	upper5[i] = results5[i] + 2*seprop5
	midpoints5[i] = median(sub5$enrichment,na.rm=TRUE)
	
	sub6 = subset(input6, enrichment >=bins6[i] & enrichment <bins6[i+1])
	results6[i] = mean(sub6$HapMap_recombination_rate,na.rm=TRUE)
	seprop6 = stderrmean(sub6$HapMap_recombination_rate)
	lower6[i] = results6[i] - 2*seprop6
	upper6[i] = results6[i] + 2*seprop6
	midpoints6[i] = median(sub6$enrichment,na.rm=TRUE)
	
	sub7 = subset(input7, enrichment >=bins7[i] & enrichment <bins7[i+1])
	results7[i] = mean(sub7$HapMap_recombination_rate,na.rm=TRUE)
	seprop7 = stderrmean(sub7$HapMap_recombination_rate)
	lower7[i] = results7[i] - 2*seprop7
	upper7[i] = results7[i] + 2*seprop7
	midpoints7[i] = median(sub7$enrichment,na.rm=TRUE)
}

#only plot 1,(2,3,5),4,6,7
colorprof = c(133,54,148,51,565,548,640)
plot(midpoints1,results1,xlab="HEK293T PRDM9 enrichment",ylab="Rec. Rate (cM/Mb)",main="",col=colors()[colorprof[1]],pch=16,xlim=c(min(c(midpoints1,midpoints2,midpoints4,midpoints6,midpoints7)),max(c(midpoints1,midpoints2,midpoints4,midpoints6,midpoints7))),ylim=c(min(c(lower1,lower2,lower4,lower6,lower7)),max(c(upper1,upper2,upper4,upper6,upper7))))
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[colorprof[1]])
points(midpoints2,results2,col=colors()[colorprof[2]],pch=16)
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[colorprof[2]])
#points(midpoints3,results3,col=colors()[colorprof[3]],pch=16)
#segments(midpoints3,lower3,midpoints3,upper3,col=colors()[colorprof[3]])
points(midpoints4,results4,col=colors()[colorprof[4]],pch=16)
segments(midpoints4,lower4,midpoints4,upper4,col=colors()[colorprof[4]])
#points(midpoints5,results5,col=colors()[colorprof[5]],pch=16)
#segments(midpoints5,lower5,midpoints5,upper5,col=colors()[colorprof[5]])
points(midpoints6,results6,col=colors()[colorprof[6]],pch=16)
segments(midpoints6,lower6,midpoints6,upper6,col=colors()[colorprof[6]])
points(midpoints7,results7,col=colors()[colorprof[7]],pch=16)
segments(midpoints7,lower7,midpoints7,upper7,col=colors()[colorprof[7]])
#legend("topleft",legend=c("M 1","M 2/3/5","M 4","M 6","M 7"),col=colors()[colorprof[c(1,2,4,6,7)]],pch=c(16),bty='n')

lmNot7 = lm(V2~V1,as.data.frame(cbind(c(midpoints1,midpoints2,midpoints4,midpoints6),c(results1,results2,results4,results6))))
abline(reg=lmNot7)
lm7 = lm(results7~midpoints7,as.data.frame(cbind(midpoints7,results7)))
abline(reg=lm7,col=colors()[colorprof[7]])

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1d_left.pdf",sep=""),width=3.5,height=5)


#right plot
library(boot)
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
motlens = c(32,31,32,50,24,23,31)
colorprof = c(133,54,148,51,565,548,640)

finaldata=read.table(file="PrattoPeakInfo.WithOverlappingPRDM9PeakInfo.txt",header=TRUE)
data = subset(finaldata,(AA_hotspots==1 | AB_hotspots==1) & is.na(WhichMotif)==F) #| AB_hotspots==1

ABonly = subset(finaldata,(AA_hotspots==0 & AB_hotspots==1) & is.na(WhichMotif)==F)
allA = subset(finaldata,(AA_hotspots==1) & is.na(WhichMotif)==F)
ABonly[ABonly$WhichMotif==3 | ABonly$WhichMotif==5,"WhichMotif"]=2
allA[allA$WhichMotif==3 | allA$WhichMotif==5,"WhichMotif"]=2

nummots=5
motlist=c(1,2,4,6,7)
midpoints = rep(0,nummots)
results = rep(0,nummots)
lower = rep(0,nummots)
upper = rep(0,nummots)
for(m in motlist){
	motifprops = function(x,i,j=m) mean(ABonly$WhichMotif[i]==j,na.rm=T) / mean(allA$WhichMotif[i]==j,na.rm=T)
	b1=boot(ABonly,motifprops,1000,sim="ordinary",stype="i")
	c1=boot.ci(b1,conf=0.95,type="basic")
	midpoints[m] = (mean(ABonly$WhichMotif==m)+mean(allA$WhichMotif==m))/2
	results[m] = b1$t0
	lower[m] = c1$basic[4]
	upper[m] = c1$basic[5]
}
pm = c(1,2,4,6,7)
results = results[pm]
lower = lower[pm]
upper = upper[pm]

#results
#0.7014925 0.8589744 1.0470588 0.6447368 2.0839161

plot(results,-1:-5,col=colors()[colorprof[pm]],pch=18,cex=1.5,ylab="",xlab="p(AB)/p(AA)",xlim=c(min(lower),max(upper)),yaxt='n',xaxt='n')
abline(v=1,col='gray',lty=3)
points(results,-1:-5,col=colors()[colorprof[pm]],pch=18,cex=1.5)
segments(lower,-1:-5,upper,-1:-5,col=colors()[colorprof[pm]],lwd=2)
axis(1,at=c(0,0.5,1.0,1.5,2.0))

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1d_right.pdf",sep=""),width=2,height=4)


##################################
### Figure 1- Figure supplement 1a
##################################
rm(list=ls())

##code used to compute overlap fractions
##list bed files with positions of DSB hotspots (from Pratto 2014), with various filters applied (telomeric, non-telomeric, all)
#filesa=c("PrattoDSB.Aintersect.lt3kb.autosomal.AB1heat.NonTel.bed","PrattoDSB.Aintersect.lt3kb.autosomal.AB1heat.Tel.bed","PrattoDSB.Aintersect.lt3kb.autosomal.AB1heat.bed")
##list bed files with positions of PRDM9 peak centres, applying various p-value thresholds between 1E-8 and 1E-3
#filesb=c("../1_CallPeaks/SingleBasePeaks.YFP_HumanPRDM9.p0.00000001.sep250.autosomal.posonly.bed","../1_CallPeaks/SingleBasePeaks.YFP_HumanPRDM9.p0.0000001.sep250.autosomal.posonly.bed","../1_CallPeaks/SingleBasePeaks.YFP_HumanPRDM9.p0.000001.sep250.autosomal.posonly.bed","../1_CallPeaks/SingleBasePeaks.YFP_HumanPRDM9.p0.00001.sep250.autosomal.posonly.bed","../1_CallPeaks/SingleBasePeaks.YFP_HumanPRDM9.p0.0001.sep250.autosomal.posonly.bed","../1_CallPeaks/SingleBasePeaks.YFP_HumanPRDM9.p0.001.sep250.autosomal.posonly.bed")
#
#lena=length(filesa)
#lenb=length(filesb)
#results = matrix(0,lenb,lena)
##declare a function to determine overlapping fractions between bedfiles and perform chance overlap correction
#getcorrected=function(filea,fileb){
#	pratto=read.table(filea,header=F)
#	widths = pratto[,3]-pratto[,2]
#	gsize = 2634892424
#	tot=length(widths)
#
#	peaknum=as.integer((strsplit(system(paste("wc -l",fileb),intern=TRUE)," ")[[1]])[1])
#	found = as.integer(system(paste("bedtools intersect -u -a",filea,"-b",fileb,"| wc -l"),intern=TRUE))
#	chance=sum(1-exp(peaknum*(log(gsize-widths)-log(gsize))))
#
#	corrected = (1-((1-(found/tot))/(1-(chance/tot))))
#
#	return(corrected)
#}
#
#for(i in 1:lena){
#	for(j in 1:lenb){
#		results[j,i]=getcorrected(filesa[i],filesb[j])
#	}
#}

#these are the results
threshes = c(-8,-7,-6,-5,-4,-3)
props = c(0.492,0.528,0.569,0.617,0.675,0.743)
propsN = c(0.580,0.617,0.654,0.699,0.755,0.816)
propsT= c(0.389,0.422,0.468,0.520,0.579,0.656)

plot(threshes,props,type='l',ylim=c(0,1),xlab="",ylab="")
lines(threshes,propsT,col='darkturquoise')
lines(threshes,propsN,col='orange3')
points(threshes,props,pch=16)
points(threshes,propsT,col='darkturquoise',pch=16)
points(threshes,propsN,col='orange3',pch=16)
#legend("bottomright",legend=c("Non-Telomeric","ALL","Telomeric"),col=c('orange3','black','darkturquoise'),pch=16)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S1a.pdf",sep=""),width=3.5,height=4)


##################################
### Figure 1- Figure supplement 1b
##################################
rm(list=ls())
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
stderrprop <- function(x,n) sqrt(x*(1-x)/n)

data=read.table("PrattoDSB.Aintersect.lt3kb.autosomal.AB1heat.NonTel.PRDM9e6flagged.bed",header=F)

bins1=quantile(data[,5],seq(0,1,0.1),na.rm=TRUE)
props=rep(0,10)
mids=rep(0,10)
upper=rep(0,10)
lower=rep(0,10)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(data, data[,5] >=bins1[i] & data[,5] <bins1[i+1])
	widths = sub2[,3]-sub2[,2]
	gsize = 2634892424
	tot=length(widths)

	peaknum=167066
	found = sum(sub2[,6])
	chance=sum(1-exp(peaknum*(log(gsize-widths)-log(gsize))))

	corrected = (1-((1-(found/tot))/(1-(chance/tot))))

	props[i]=corrected
	mids[i]=median(sub2[,5])
	upper[i]=props[i]+2*stderrprop(corrected,dim(sub2)[1])
	lower[i]=props[i]-2*stderrprop(corrected,dim(sub2)[1])
}

dataT=read.table("PrattoDSB.Aintersect.lt3kb.autosomal.AB1heat.Tel.PRDM9e6flagged.bed",header=F)
bins1T=quantile(dataT[,5],seq(0,1,0.1),na.rm=TRUE)
propsT=rep(0,10)
midsT=rep(0,10)
upperT=rep(0,10)
lowerT=rep(0,10)
for(i in 1:(length(bins1T)-1)){
	sub2 = subset(dataT, dataT[,5] >=bins1T[i] & dataT[,5] <bins1T[i+1])
	widths = sub2[,3]-sub2[,2]
	gsize = 2634892424
	tot=length(widths)

	peaknum=167066
	found = sum(sub2[,6])
	chance=sum(1-exp(peaknum*(log(gsize-widths)-log(gsize))))

	corrected = (1-((1-(found/tot))/(1-(chance/tot))))

	propsT[i]=corrected
	midsT[i]=median(sub2[,5])
	upperT[i]=propsT[i]+2*stderrprop(corrected,dim(sub2)[1])
	lowerT[i]=propsT[i]-2*stderrprop(corrected,dim(sub2)[1])
}

dataA=read.table("PrattoDSB.Aintersect.lt3kb.autosomal.AB1heat.PRDM9e6flagged.bed",header=F)
bins1A=quantile(dataA[,5],seq(0,1,0.1),na.rm=TRUE)
propsA=rep(0,10)
midsA=rep(0,10)
upperA=rep(0,10)
lowerA=rep(0,10)
for(i in 1:(length(bins1A)-1)){
	sub2 = subset(dataA, dataA[,5] >=bins1A[i] & dataA[,5] <bins1A[i+1])
	widths = sub2[,3]-sub2[,2]
	gsize = 2634892424
	tot=length(widths)

	peaknum=167066
	found = sum(sub2[,6])
	chance=sum(1-exp(peaknum*(log(gsize-widths)-log(gsize))))

	corrected = (1-((1-(found/tot))/(1-(chance/tot))))

	propsA[i]=corrected
	midsA[i]=median(sub2[,5])
	upperA[i]=propsA[i]+2*stderrprop(corrected,dim(sub2)[1])
	lowerA[i]=propsA[i]-2*stderrprop(corrected,dim(sub2)[1])
}


plot(mids,props,main="",xlim=c(0,max(midsT)),ylim=c(0,1),col='orange3',pch=16,xlab="",ylab="")
segments(mids,lower,mids,upper,col='orange3')
points(midsT,propsT,col='darkturquoise',pch=16)
segments(midsT,lowerT,midsT,upperT,col='darkturquoise')
points(midsA,propsA,col='black',pch=16)
segments(midsA,lowerA,midsA,upperA,col='black')
#legend("bottomright",legend=c("Non-Telomeric","ALL","Telomeric"),col=c('orange3','black','darkturquoise'),pch=16)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S1b.pdf",sep=""),width=3.5,height=4)


####################################
### Figure 1- Figure supplement 1c/d
####################################
rm(list=ls())
load("finalMeansProfile.HumanK4enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q1.all.r")
finalMeansH=finalMeans
load("finalMeansProfile.UTK4enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q1.all.r")
finalMeansU=finalMeans

colorprof=c(553,554,555,556)
plot(ksmooth(seq(-5000,5000,1),finalMeansH, "normal", bandwidth=200),type='l',col=colors()[colorprof[1]],xlab="",ylab="",main="",lwd=3,ylim=c(0.4,2.1))

for(i in 2:4){
	load(paste0("finalMeansProfile.HumanK4enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q",i,".all.r"))
	finalMeansH=finalMeans
	load(paste0("finalMeansProfile.UTK4enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q",i,".all.r"))
	finalMeansU=finalMeans

	lines(ksmooth(seq(-5000,5000,1),finalMeansH, "normal", bandwidth=200),type='l',col=colors()[colorprof[i]],xlab="",ylab="",main="",lwd=3)
}

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S1c.pdf",sep=""),width=4,height=5)

rm(list=ls())
load("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q1.all.r")
finalMeansH=finalMeans
load("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q1.all.r")
finalMeansU=finalMeans

colorprof=c(468,469,470,471)
plot(ksmooth(seq(-5000,5000,1),finalMeansH, "normal", bandwidth=200),type='l',col=colors()[colorprof[1]],xlab="",ylab="",main="",lwd=3,ylim=c(0.17,0.34))

for(i in 2:4){
	load(paste0("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q",i,".all.r"))
	finalMeansH=finalMeans
	load(paste0("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q",i,".all.r"))
	finalMeansU=finalMeans

	lines(ksmooth(seq(-5000,5000,1),finalMeansH, "normal", bandwidth=200),type='l',col=colors()[colorprof[i]],xlab="",ylab="",main="",lwd=3)
}

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S1d.pdf",sep=""),width=4,height=5)


#################################
### Figure 1- Figure supplement 2
#################################

rm(list=ls())
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

dataH0=read.table("ForceCall.Pratto_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataH0)
dataH=subset(dataH0,covg>5 & covg<quantile(dataH0$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500))
dim(dataH)

dataD=read.table("ForceCall.Pratto_DMC1_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataD)
dataD=subset(dataD,covg>5 & covg<quantile(dataD$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500))
dim(dataD)

dataHH=read.table("ForceCall.YFP_Human_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataHH)
dataHH=subset(dataHH,covg>5 & covg<quantile(dataHH$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500))
dim(dataHH)


#1_S2a compare PRDM9 and HEK H3K4me3 enrichment 

input2=subset(dataHH,dataH$H3K4me3Overlap==0)
input2nt=subset(dataHH,dataH$H3K4me3Overlap==0 & dataH$TelFlag==0)
input2t=subset(dataHH,dataH$H3K4me3Overlap==0 & dataH$TelFlag==1)
  
cor(input2$enrich,input2$H.enrichment,use="complete.obs")
#0.4777871
cor(input2nt$enrich,input2nt$H.enrichment,use="complete.obs")
#0.5119167
cor(input2t$enrich,input2t$H.enrichment,use="complete.obs")
#0.489168

bins1=quantile(input2$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}

plot(midpoints1,results2,xlab="HEK293T PRDM9",ylab="HEK293T H3K4me3",col='black',pch=16,main="",xlim=c(0,8),ylim=c(0,1.5))
segments(midpoints1,lower2,midpoints1,upper2,col='black')

bins1=quantile(input2nt$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2nt, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}
points(midpoints1,results2,xlab="HEK293T PRDM9",ylab="HEK293T H3K4me3",col='orange3',pch=16)
segments(midpoints1,lower2,midpoints1,upper2,col='orange3')


bins1=quantile(input2t$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2t, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}
points(midpoints1,results2,xlab="HEK293T PRDM9",ylab="HEK293T H3K4me3",col='darkturquoise',pch=16)
segments(midpoints1,lower2,midpoints1,upper2,col='darkturquoise')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S2a.pdf",sep=""),width=3.5,height=4)


#1_S2b compare PRDM9 and Testis H3K4me3 enrichment 

input2=subset(dataH,dataH$H3K4me3Overlap==0)
input2nt=subset(dataH,dataH$H3K4me3Overlap==0 & dataH$TelFlag==0)
input2t=subset(dataH,dataH$H3K4me3Overlap==0 & dataH$TelFlag==1)
  
cor(input2$enrich,input2$H.enrichment,use="complete.obs")
#0.4777871
cor(input2nt$enrich,input2nt$H.enrichment,use="complete.obs")
#0.5119167
cor(input2t$enrich,input2t$H.enrichment,use="complete.obs")
#0.489168

bins1=quantile(input2$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}

plot(midpoints1,results2,xlab="HEK293T PRDM9",ylab="Testis H3K4me3",col='black',pch=16,main="",xlim=c(0,8),ylim=c(0,0.8))
segments(midpoints1,lower2,midpoints1,upper2,col='black')

bins1=quantile(input2nt$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2nt, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}
points(midpoints1,results2,xlab="HEK293T PRDM9",ylab="HEK293T H3K4me3",col='orange3',pch=16)
segments(midpoints1,lower2,midpoints1,upper2,col='orange3')


bins1=quantile(input2t$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2t, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}
points(midpoints1,results2,xlab="HEK293T PRDM9",ylab="HEK293T H3K4me3",col='darkturquoise',pch=16)
segments(midpoints1,lower2,midpoints1,upper2,col='darkturquoise')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S2b.pdf",sep=""),width=3.5,height=4)


#1_S2c compare PRDM9 and Testis DMC1 enrichment 

input2=subset(dataD,dataH$H3K4me3Overlap==0)
input2nt=subset(dataD,dataH$H3K4me3Overlap==0 & dataH$TelFlag==0)
input2t=subset(dataD,dataH$H3K4me3Overlap==0 & dataH$TelFlag==1)
  
cor(input2$enrich,input2$H.enrichment,use="complete.obs")
#0.2098
cor(input2nt$enrich,input2nt$H.enrichment,use="complete.obs")
#0.3627
cor(input2t$enrich,input2t$H.enrichment,use="complete.obs")
#0.2252

bins1=quantile(input2$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}

plot(midpoints1,results2,xlab="HEK293T PRDM9",ylab="Testis DMC1",col='black',pch=16,main="",xlim=c(0,8),ylim=c(0,3.1))
segments(midpoints1,lower2,midpoints1,upper2,col='black')

bins1=quantile(input2nt$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2nt, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}
points(midpoints1,results2,xlab="HEK293T PRDM9",ylab="HEK293T H3K4me3",col='orange3',pch=16)
segments(midpoints1,lower2,midpoints1,upper2,col='orange3')


bins1=quantile(input2t$enrich,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2t, enrich >=bins1[i] & enrich <bins1[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints1[i] = median(sub2$enrich,na.rm=TRUE)
}
points(midpoints1,results2,xlab="HEK293T PRDM9",ylab="HEK293T H3K4me3",col='darkturquoise',pch=16)
segments(midpoints1,lower2,midpoints1,upper2,col='darkturquoise')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S2c.pdf",sep=""),width=3.5,height=4)


#1_S2d compare HEK293 H3K4me3 and Testis H3K4me3 enrichment 

input3=subset(dataH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataHH$H.enrichment) & H.enrichment>0)
input2=subset(dataHH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataH$H.enrichment) & dataH$H.enrichment>0)
input3t=subset(dataH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataHH$H.enrichment) & H.enrichment>0 & dataH$TelFlag==1)
input2t=subset(dataHH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataH$H.enrichment) & dataH$H.enrichment>0 & dataH$TelFlag==1)
input3nt=subset(dataH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataHH$H.enrichment) & H.enrichment>0 & dataH$TelFlag==0)
input2nt=subset(dataHH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataH$H.enrichment) & dataH$H.enrichment>0 & dataH$TelFlag==0)

cor(input2$H.enrichment,input3$H.enrichment,use="complete.obs")
#0.2281716
cor(input2t$H.enrichment,input3t$H.enrichment,use="complete.obs")
#0.1740453
cor(input2nt$H.enrichment,input3nt$H.enrichment,use="complete.obs")
#0.3008702


stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
bins1=quantile(input2$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
plot(midpoints1,results2,ylab="Testis H3K4me3",xlab="HEK293T H3K4me3",col='black',pch=16,main="",xlim=c(0,1.8),ylim=c(0,0.6))
segments(midpoints1,lower2,midpoints1,upper2,col='black')

input2=input2t
input3=input3t
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
bins1=quantile(input2$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
points(midpoints1,results2,ylab="Testis H3K4me3",xlab="HEK293T H3K4me3",col='darkturquoise',pch=16,main="",xlim=c(0,max(midpoints1)),ylim=c(0,max(upper2)))
segments(midpoints1,lower2,midpoints1,upper2,col='darkturquoise')

input2=input2nt
input3=input3nt
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
bins1=quantile(input2$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
points(midpoints1,results2,ylab="Testis H3K4me3",xlab="HEK293T H3K4me3",col='orange3',pch=16,main="",xlim=c(0,max(midpoints1)),ylim=c(0,max(upper2)))
segments(midpoints1,lower2,midpoints1,upper2,col='orange3')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S2d.pdf",sep=""),width=3.5,height=4)



#1_S2e compare HEK293 H3K4me3 and Testis DMC1 enrichment 

input2=subset(dataHH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataD$H.enrichment) & H.enrichment>0)
input3=subset(dataD,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataHH$H.enrichment) & dataHH$H.enrichment>0)
input2t=subset(dataHH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataD$H.enrichment) & H.enrichment>0 & dataH$TelFlag==1)
input2nt=subset(dataHH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataD$H.enrichment) & H.enrichment>0 & dataH$TelFlag==0)
input3t=subset(dataD,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataHH$H.enrichment) & dataHH$H.enrichment>0 & dataH$TelFlag==1)
input3nt=subset(dataD,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataHH$H.enrichment) & dataHH$H.enrichment>0 & dataH$TelFlag==0)


cor(input2$H.enrichment,input3$H.enrichment,use="complete.obs")
#0.2042315
cor(input2t$H.enrichment,input3t$H.enrichment,use="complete.obs")
#0.1776932
cor(input2nt$H.enrichment,input3nt$H.enrichment,use="complete.obs")
#0.2855162

stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
bins1=quantile(input2$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
plot(midpoints1,results2,xlab="HEK293T H3K4me3",ylab="Testis DMC1",col='black',pch=16,main="",xlim=c(0,1.8),ylim=c(0,3))
segments(midpoints1,lower2,midpoints1,upper2,col='black')

input2 =input2t
input3 =input3t
bins1=quantile(input2t$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
points(midpoints1,results2,xlab="HEK293T H3K4me3",ylab="Testis DMC1",col='darkturquoise',pch=16,main="",xlim=c(0,max(midpoints1)),ylim=c(0,max(upper2)))
segments(midpoints1,lower2,midpoints1,upper2,col='darkturquoise')

input2 =input2nt
input3 =input3nt
bins1=quantile(input2t$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
points(midpoints1,results2,xlab="Testis H3K4me3",ylab="Testis DMC1",col='orange3',pch=16,main="",xlim=c(0,max(midpoints1)),ylim=c(0,max(upper2)))
segments(midpoints1,lower2,midpoints1,upper2,col='orange3')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S2e.pdf",sep=""),width=3.5,height=4)


#1_S2f compare Testis H3K4me3 and Testis DMC1 enrichment 

input2=subset(dataH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataD$H.enrichment) & H.enrichment>0)
input3=subset(dataD,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataH$H.enrichment) & dataH$H.enrichment>0)
input2t=subset(dataH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataD$H.enrichment) & H.enrichment>0 & dataH$TelFlag==1)
input2nt=subset(dataH,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataD$H.enrichment) & H.enrichment>0 & dataH$TelFlag==0)
input3t=subset(dataD,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataH$H.enrichment) & dataH$H.enrichment>0 & dataH$TelFlag==1)
input3nt=subset(dataD,dataH$H3K4me3Overlap==0 & !is.na(H.enrichment) & !is.na(dataH$H.enrichment) & dataH$H.enrichment>0 & dataH$TelFlag==0)


cor(input2$H.enrichment,input3$H.enrichment,use="complete.obs")
#0.3079681
cor(input2t$H.enrichment,input3t$H.enrichment,use="complete.obs")
#0.3418105
cor(input2nt$H.enrichment,input3nt$H.enrichment,use="complete.obs")
#0.5508352

stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
bins1=quantile(input2$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
plot(midpoints1,results2,xlab="Testis H3K4me3",ylab="Testis DMC1",col='black',pch=16,main="",xlim=c(0,max(midpoints1)),ylim=c(0,4.5))
segments(midpoints1,lower2,midpoints1,upper2,col='black')

input2 =input2t
input3 =input3t
bins1=quantile(input2t$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
points(midpoints1,results2,xlab="Testis H3K4me3",ylab="Testis DMC1",col='darkturquoise',pch=16,main="",xlim=c(0,max(midpoints1)),ylim=c(0,max(upper2)))
segments(midpoints1,lower2,midpoints1,upper2,col='darkturquoise')

input2 =input2nt
input3 =input3nt
bins1=quantile(input2t$H.enrichment,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(input2, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	midpoints1[i] = median(sub2$H.enrichment,na.rm=TRUE)
	sub3 = subset(input3, input2$H.enrichment >=bins1[i] & input2$H.enrichment <bins1[i+1])
	results2[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub3$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
}
points(midpoints1,results2,xlab="Testis H3K4me3",ylab="Testis DMC1",col='orange3',pch=16,main="",xlim=c(0,max(midpoints1)),ylim=c(0,max(upper2)))
segments(midpoints1,lower2,midpoints1,upper2,col='orange3')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S2f.pdf",sep=""),width=3.5,height=4)



##################################
### Figure 1- Figure supplement 4a
##################################
rm(list=ls())
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
motlens = c(32,31,32,50,24,23,31)

dataH0=read.table("ForceCall.Pratto_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataH0)
dataH=subset(dataH0,covg>5 & covg<quantile(dataH0$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500) & dataH0$AnyH3K4me3Overlap==0)
dim(dataH)

dataHH=read.table("ForceCall.YFP_Human_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataHH)
dataHH=subset(dataHH,covg>5 & covg<quantile(dataHH$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500) & dataH0$AnyH3K4me3Overlap==0)
dim(dataHH)
dataHH$WhichMotif=dataH$WhichMotif

dataSnoprom = dataHH 

mot1noprom = subset(dataSnoprom,WhichMotif==1 )
mot2noprom = subset(dataSnoprom,WhichMotif==2 | WhichMotif==3 | WhichMotif==5)
mot3noprom = subset(dataSnoprom,WhichMotif==3 )
mot4noprom = subset(dataSnoprom,WhichMotif==4 )
mot5noprom = subset(dataSnoprom,WhichMotif==5 )
mot6noprom = subset(dataSnoprom,WhichMotif==6 )
mot7noprom = subset(dataSnoprom,WhichMotif==7 )

input1=mot1noprom
input2=mot2noprom
input3=mot3noprom
input4=mot4noprom
input5=mot5noprom
input6=mot6noprom
input7=mot7noprom

bins1=quantile(input1$enrich,seq(0,1,0.25),na.rm=TRUE)
bins2=quantile(input2$enrich,seq(0,1,0.25),na.rm=TRUE)
bins3=quantile(input3$enrich,seq(0,1,0.25),na.rm=TRUE)
bins4=quantile(input4$enrich,seq(0,1,0.25),na.rm=TRUE)
bins5=quantile(input5$enrich,seq(0,1,0.25),na.rm=TRUE)
bins6=quantile(input6$enrich,seq(0,1,0.25),na.rm=TRUE)
bins7=quantile(input7$enrich,seq(0,1,0.25),na.rm=TRUE)

bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
bins3[length(bins3)]=bins3[length(bins3)]+1
bins4[length(bins4)]=bins4[length(bins4)]+1
bins5[length(bins5)]=bins5[length(bins5)]+1
bins6[length(bins6)]=bins6[length(bins6)]+1
bins7[length(bins7)]=bins7[length(bins7)]+1

results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins2)-1)
results3 = rep(0,length(bins3)-1)
results4 = rep(0,length(bins4)-1)
results5 = rep(0,length(bins5)-1)
results6 = rep(0,length(bins6)-1)
results7 = rep(0,length(bins7)-1)

midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins2)-1)
midpoints3 = rep(0,length(bins3)-1)
midpoints4 = rep(0,length(bins4)-1)
midpoints5 = rep(0,length(bins5)-1)
midpoints6 = rep(0,length(bins6)-1)
midpoints7 = rep(0,length(bins7)-1)

lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins2)-1)
upper2 = rep(NA,length(bins2)-1)
lower3 = rep(NA,length(bins3)-1)
upper3 = rep(NA,length(bins3)-1)
lower4 = rep(NA,length(bins4)-1)
upper4 = rep(NA,length(bins4)-1)
lower5 = rep(NA,length(bins5)-1)
upper5 = rep(NA,length(bins5)-1)
lower6 = rep(NA,length(bins6)-1)
upper6 = rep(NA,length(bins6)-1)
lower7 = rep(NA,length(bins7)-1)
upper7 = rep(NA,length(bins7)-1)

for(i in 1:(length(bins1)-1)){
	sub1 = subset(input1, enrich >=bins1[i] & enrich <bins1[i+1])
	results1[i] = mean(sub1$H.enrichment,na.rm=TRUE)
	seprop1 = stderrmean(sub1$H.enrichment)
	lower1[i] = results1[i] - 2*seprop1
	upper1[i] = results1[i] + 2*seprop1
	midpoints1[i] = median(sub1$enrich,na.rm=TRUE)
	
	sub2 = subset(input2, enrich >=bins2[i] & enrich <bins2[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints2[i] = median(sub2$enrich,na.rm=TRUE)
	
	sub3 = subset(input3, enrich >=bins3[i] & enrich <bins3[i+1])
	results3[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop3 = stderrmean(sub3$H.enrichment)
	lower3[i] = results3[i] - 2*seprop3
	upper3[i] = results3[i] + 2*seprop3
	midpoints3[i] = median(sub3$enrich,na.rm=TRUE)
	
	sub4 = subset(input4, enrich >=bins4[i] & enrich <bins4[i+1])
	results4[i] = mean(sub4$H.enrichment,na.rm=TRUE)
	seprop4 = stderrmean(sub4$H.enrichment)
	lower4[i] = results4[i] - 2*seprop4
	upper4[i] = results4[i] + 2*seprop4
	midpoints4[i] = median(sub4$enrich,na.rm=TRUE)
	
	sub5 = subset(input5, enrich >=bins5[i] & enrich <bins5[i+1])
	results5[i] = mean(sub5$H.enrichment,na.rm=TRUE)
	seprop5 = stderrmean(sub5$H.enrichment)
	lower5[i] = results5[i] - 2*seprop5
	upper5[i] = results5[i] + 2*seprop5
	midpoints5[i] = median(sub5$enrich,na.rm=TRUE)
	
	sub6 = subset(input6, enrich >=bins6[i] & enrich <bins6[i+1])
	results6[i] = mean(sub6$H.enrichment,na.rm=TRUE)
	seprop6 = stderrmean(sub6$H.enrichment)
	lower6[i] = results6[i] - 2*seprop6
	upper6[i] = results6[i] + 2*seprop6
	midpoints6[i] = median(sub6$enrich,na.rm=TRUE)
	
	sub7 = subset(input7, enrich >=bins7[i] & enrich <bins7[i+1])
	results7[i] = mean(sub7$H.enrichment,na.rm=TRUE)
	seprop7 = stderrmean(sub7$H.enrichment)
	lower7[i] = results7[i] - 2*seprop7
	upper7[i] = results7[i] + 2*seprop7
	midpoints7[i] = median(sub7$enrich,na.rm=TRUE)
}

#only plot 1,(2,3,5),4,6,7
colorprof = c(133,54,148,51,565,548,640)
plot(midpoints1,results1,xlab="HEK293T PRDM9 enrichment",ylab="HEK293T H3K4me3 enrich.",main="",col=colors()[colorprof[1]],pch=16,xlim=c(min(c(midpoints1,midpoints2,midpoints4,midpoints6,midpoints7)),max(c(midpoints1,midpoints2,midpoints4,midpoints6,midpoints7))),ylim=c(min(c(lower1,lower2,lower4,lower6,lower7)),max(c(upper1,upper2,upper4,upper6,upper7))))
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[colorprof[1]])
points(midpoints2,results2,col=colors()[colorprof[2]],pch=16)
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[colorprof[2]])
#points(midpoints3,results3,col=colors()[colorprof[3]],pch=16)
#segments(midpoints3,lower3,midpoints3,upper3,col=colors()[colorprof[3]])
points(midpoints4,results4,col=colors()[colorprof[4]],pch=16)
segments(midpoints4,lower4,midpoints4,upper4,col=colors()[colorprof[4]])
#points(midpoints5,results5,col=colors()[colorprof[5]],pch=16)
#segments(midpoints5,lower5,midpoints5,upper5,col=colors()[colorprof[5]])
points(midpoints6,results6,col=colors()[colorprof[6]],pch=16)
segments(midpoints6,lower6,midpoints6,upper6,col=colors()[colorprof[6]])
points(midpoints7,results7,col=colors()[colorprof[7]],pch=16)
segments(midpoints7,lower7,midpoints7,upper7,col=colors()[colorprof[7]])
#legend("topleft",legend=c("M 1","M 2/3/5","M 4","M 6","M 7"),col=colors()[colorprof[c(1,2,4,6,7)]],pch=c(16),bty='n')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S4a.pdf",sep=""),width=3.5,height=4)


##################################
### Figure 1- Figure supplement 4b
##################################
rm(list=ls())
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
motlens = c(32,31,32,50,24,23,31)

dataH0=read.table("ForceCall.Pratto_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataH0)
dataH=subset(dataH0,covg>5 & covg<quantile(dataH0$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500) & dataH0$AnyH3K4me3Overlap==0)
dim(dataH)

dataHH=read.table("ForceCall.YFP_Human_H3K4me3_on_YFP_HumanPRDM9.autosomal.bed",header=TRUE)
dim(dataHH)
dataHH=subset(dataHH,covg>5 & covg<quantile(dataHH$covg,0.999) & dataH0$H.cov.g >5 & dataH0$H.cov.g <150 & (dataH0$H.cov.r1 + dataH0$H.cov.r2 <500) & dataH0$AnyH3K4me3Overlap==0)
dim(dataHH)
dataHH$WhichMotif=dataH$WhichMotif

dataSnoprom = dataH 

mot1noprom = subset(dataSnoprom,WhichMotif==1 )
mot2noprom = subset(dataSnoprom,WhichMotif==2 | WhichMotif==3 | WhichMotif==5)
mot3noprom = subset(dataSnoprom,WhichMotif==3 )
mot4noprom = subset(dataSnoprom,WhichMotif==4 )
mot5noprom = subset(dataSnoprom,WhichMotif==5 )
mot6noprom = subset(dataSnoprom,WhichMotif==6 )
mot7noprom = subset(dataSnoprom,WhichMotif==7 )

input1=mot1noprom
input2=mot2noprom
input3=mot3noprom
input4=mot4noprom
input5=mot5noprom
input6=mot6noprom
input7=mot7noprom

bins1=quantile(input1$enrich,seq(0,1,0.25),na.rm=TRUE)
bins2=quantile(input2$enrich,seq(0,1,0.25),na.rm=TRUE)
bins3=quantile(input3$enrich,seq(0,1,0.25),na.rm=TRUE)
bins4=quantile(input4$enrich,seq(0,1,0.25),na.rm=TRUE)
bins5=quantile(input5$enrich,seq(0,1,0.25),na.rm=TRUE)
bins6=quantile(input6$enrich,seq(0,1,0.25),na.rm=TRUE)
bins7=quantile(input7$enrich,seq(0,1,0.25),na.rm=TRUE)

bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
bins3[length(bins3)]=bins3[length(bins3)]+1
bins4[length(bins4)]=bins4[length(bins4)]+1
bins5[length(bins5)]=bins5[length(bins5)]+1
bins6[length(bins6)]=bins6[length(bins6)]+1
bins7[length(bins7)]=bins7[length(bins7)]+1

results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins2)-1)
results3 = rep(0,length(bins3)-1)
results4 = rep(0,length(bins4)-1)
results5 = rep(0,length(bins5)-1)
results6 = rep(0,length(bins6)-1)
results7 = rep(0,length(bins7)-1)

midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins2)-1)
midpoints3 = rep(0,length(bins3)-1)
midpoints4 = rep(0,length(bins4)-1)
midpoints5 = rep(0,length(bins5)-1)
midpoints6 = rep(0,length(bins6)-1)
midpoints7 = rep(0,length(bins7)-1)

lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins2)-1)
upper2 = rep(NA,length(bins2)-1)
lower3 = rep(NA,length(bins3)-1)
upper3 = rep(NA,length(bins3)-1)
lower4 = rep(NA,length(bins4)-1)
upper4 = rep(NA,length(bins4)-1)
lower5 = rep(NA,length(bins5)-1)
upper5 = rep(NA,length(bins5)-1)
lower6 = rep(NA,length(bins6)-1)
upper6 = rep(NA,length(bins6)-1)
lower7 = rep(NA,length(bins7)-1)
upper7 = rep(NA,length(bins7)-1)

for(i in 1:(length(bins1)-1)){
	sub1 = subset(input1, enrich >=bins1[i] & enrich <bins1[i+1])
	results1[i] = mean(sub1$H.enrichment,na.rm=TRUE)
	seprop1 = stderrmean(sub1$H.enrichment)
	lower1[i] = results1[i] - 2*seprop1
	upper1[i] = results1[i] + 2*seprop1
	midpoints1[i] = median(sub1$enrich,na.rm=TRUE)
	
	sub2 = subset(input2, enrich >=bins2[i] & enrich <bins2[i+1])
	results2[i] = mean(sub2$H.enrichment,na.rm=TRUE)
	seprop2 = stderrmean(sub2$H.enrichment)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints2[i] = median(sub2$enrich,na.rm=TRUE)
	
	sub3 = subset(input3, enrich >=bins3[i] & enrich <bins3[i+1])
	results3[i] = mean(sub3$H.enrichment,na.rm=TRUE)
	seprop3 = stderrmean(sub3$H.enrichment)
	lower3[i] = results3[i] - 2*seprop3
	upper3[i] = results3[i] + 2*seprop3
	midpoints3[i] = median(sub3$enrich,na.rm=TRUE)
	
	sub4 = subset(input4, enrich >=bins4[i] & enrich <bins4[i+1])
	results4[i] = mean(sub4$H.enrichment,na.rm=TRUE)
	seprop4 = stderrmean(sub4$H.enrichment)
	lower4[i] = results4[i] - 2*seprop4
	upper4[i] = results4[i] + 2*seprop4
	midpoints4[i] = median(sub4$enrich,na.rm=TRUE)
	
	sub5 = subset(input5, enrich >=bins5[i] & enrich <bins5[i+1])
	results5[i] = mean(sub5$H.enrichment,na.rm=TRUE)
	seprop5 = stderrmean(sub5$H.enrichment)
	lower5[i] = results5[i] - 2*seprop5
	upper5[i] = results5[i] + 2*seprop5
	midpoints5[i] = median(sub5$enrich,na.rm=TRUE)
	
	sub6 = subset(input6, enrich >=bins6[i] & enrich <bins6[i+1])
	results6[i] = mean(sub6$H.enrichment,na.rm=TRUE)
	seprop6 = stderrmean(sub6$H.enrichment)
	lower6[i] = results6[i] - 2*seprop6
	upper6[i] = results6[i] + 2*seprop6
	midpoints6[i] = median(sub6$enrich,na.rm=TRUE)
	
	sub7 = subset(input7, enrich >=bins7[i] & enrich <bins7[i+1])
	results7[i] = mean(sub7$H.enrichment,na.rm=TRUE)
	seprop7 = stderrmean(sub7$H.enrichment)
	lower7[i] = results7[i] - 2*seprop7
	upper7[i] = results7[i] + 2*seprop7
	midpoints7[i] = median(sub7$enrich,na.rm=TRUE)
}

#only plot 1,(2,3,5),4,6,7
colorprof = c(133,54,148,51,565,548,640)
plot(midpoints1,results1,xlab="HEK293T PRDM9 enrichment",ylab="Testis H3K4me3 enrich.",main="",col=colors()[colorprof[1]],pch=16,xlim=c(min(c(midpoints1,midpoints2,midpoints4,midpoints6,midpoints7)),max(c(midpoints1,midpoints2,midpoints4,midpoints6,midpoints7))),ylim=c(min(c(lower1,lower2,lower4,lower6,lower7)),max(c(upper1,upper2,upper4,upper6,upper7))))
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[colorprof[1]])
points(midpoints2,results2,col=colors()[colorprof[2]],pch=16)
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[colorprof[2]])
#points(midpoints3,results3,col=colors()[colorprof[3]],pch=16)
#segments(midpoints3,lower3,midpoints3,upper3,col=colors()[colorprof[3]])
points(midpoints4,results4,col=colors()[colorprof[4]],pch=16)
segments(midpoints4,lower4,midpoints4,upper4,col=colors()[colorprof[4]])
#points(midpoints5,results5,col=colors()[colorprof[5]],pch=16)
#segments(midpoints5,lower5,midpoints5,upper5,col=colors()[colorprof[5]])
points(midpoints6,results6,col=colors()[colorprof[6]],pch=16)
segments(midpoints6,lower6,midpoints6,upper6,col=colors()[colorprof[6]])
points(midpoints7,results7,col=colors()[colorprof[7]],pch=16)
segments(midpoints7,lower7,midpoints7,upper7,col=colors()[colorprof[7]])
#legend("topleft",legend=c("M 1","M 2/3/5","M 4","M 6","M 7"),col=colors()[colorprof[c(1,2,4,6,7)]],pch=c(16),bty='n')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S4b.pdf",sep=""),width=3.5,height=4)


##################################
### Figure 1- Figure supplement 4c
##################################
#make q-q plot of DMC1 data with error bars
rm(list=ls())
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
motlens = c(32,31,32,50,24,23,31)
library(boot)

finaldata = read.table("PrattoPeakInfo.WithOverlappingPRDM9PeakInfo.txt",header=T)
data = subset(finaldata,(AA_hotspots==1 | AB_hotspots==1) & is.na(WhichMotif)==F)
#[1] 14220    56
table(data$WhichMotif)
#   1    2    3    4    5    6    7 
#3413  663  684 2676 1070 3308 2406


data$ABenrich = (data$AB1_strength + data$AB2_strength)
data$AAenrich = (data$AA1_strength + data$AA2_strength)


input1=subset(data,WhichMotif==1 )
input2=subset(data,WhichMotif==2 | WhichMotif==3 | WhichMotif==5)
input4=subset(data,WhichMotif==4 )
input6=subset(data,WhichMotif==6 )
input7=subset(data,WhichMotif==7 )

inputlist = list(input1,input2,input4,input6,input7)

aaPoints = matrix(0,nrow=4,ncol=5)
abPoints = matrix(0,nrow=4,ncol=5)

aaLower = matrix(0,nrow=4,ncol=5)
aaUpper = matrix(0,nrow=4,ncol=5)

abLower = matrix(0,nrow=4,ncol=5)
abUpper = matrix(0,nrow=4,ncol=5)


quantilefunc = function(x,i,q) quantile(x[i],q)
quantiles = seq(0.125,0.875,0.25)
for(ip in 1:5){
	input = inputlist[[ip]]
	for(q1 in 1:4){
		b1=boot(input$AAenrich,quantilefunc,q=quantiles[q1],1000,sim="ordinary",stype="i")
		c1=boot.ci(b1,conf=0.95,type="basic")
		aaPoints[q1,ip]=b1$t0
		aaLower[q1,ip]=c1$basic[4]
		aaUpper[q1,ip]=c1$basic[5]
		b1=boot(input$ABenrich,quantilefunc,q=quantiles[q1],1000,sim="ordinary",stype="i")
		c1=boot.ci(b1,conf=0.95,type="basic")
		abPoints[q1,ip]=b1$t0
		abLower[q1,ip]=c1$basic[4]
		abUpper[q1,ip]=c1$basic[5]
	}
}

colorprof = c(133,54,51,548,640)
plot(abPoints[,1],aaPoints[,1],xlab="AB DMC1 heat",ylab="AA DMC1 heat",main="",col=colors()[colorprof[1]],pch=16,xlim=c(0,225),ylim=c(0,350))
segments(abLower[,1],aaPoints[,1],abUpper[,1],aaPoints[,1],col=colors()[colorprof[1]])
segments(abPoints[,1],aaLower[,1],abPoints[,1],aaUpper[,1],col=colors()[colorprof[1]])

for(i in 2:5){
	points(abPoints[,i],aaPoints[,i],xlab="AB DMC1 heat",ylab="AA DMC1 heat",main="",col=colors()[colorprof[i]],pch=16,xlim=c(0,250),ylim=c(0,350))
	segments(abLower[,i],aaPoints[,i],abUpper[,i],aaPoints[,i],col=colors()[colorprof[i]])
	segments(abPoints[,i],aaLower[,i],abPoints[,i],aaUpper[,i],col=colors()[colorprof[i]])
}
#legend("topleft",legend=c("M 1","M 2/3/5","M 4","M 6","M 7"),col=colors()[colorprof],pch=c(16),bty='n')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure1_S4c.pdf",sep=""),width=3.5,height=4)

#############
### Figure 2a
#############
rm(list=ls())
data=read.table("ForceCall.HumanK4onTSSpositions.autosomal.plusFPKM.plusPRDM9enrich.bed",header=T)
data=data[order(data$gene),]
data2=read.table("ForceCall.UTK4onTSSpositions.autosomal.plusFPKMv2.bed",header=T)
data2=data2[order(data2$gene),]
data2$PRDM9enrich = data$PRDM9enrich

sub = subset(data2,is.na(H.pval)==F & H.cov.g>5 & H.pval<=1)
dim(sub)
#[1] 18625    28
sub1 = subset(sub,H.pval<0.00001)
#12915
sub2 = subset(sub,H.pval>=0.5)
#4044

input1 = sub 
input0 = subset(sub,H.enrichment==0)
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
colorprof=rev(c(81,257,657,495))
colorprof=c(494,495,496,497)

bins1=quantile(input1$H.enrichment,seq(0,1,0.25),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
bins2=quantile(input1$PRDM9enrich,seq(0,1,0.25),na.rm=TRUE)
bins2[length(bins2)]=bins2[length(bins2)]+1

for(p in 1:(length(bins2)-1)){

results1 = c(0,length(bins1)-1)
midpoints1 = c(0,length(bins1)-1)
lower1 = c(NA,length(bins1)-1)
upper1 = c(NA,length(bins1)-1)

for(i in 1:(length(bins1)-1)){
	sub1 = subset(input1, H.enrichment >=bins1[i] & H.enrichment <bins1[i+1])
	found1 = sum(sub1$PRDM9enrich>=bins2[p],na.rm=TRUE) 
	tot1 = dim(sub1)[1]
	corrected = (1-((1-(found1/tot1))/(1-(0.06255194))))
	results1[i] = corrected
	seprop1 = stderrprop(results1[i],dim(sub1)[1])
	lower1[i] = results1[i] - 2*seprop1
	upper1[i] = results1[i] + 2*seprop1
	midpoints1[i] = median(sub1$H.enrichment,na.rm=TRUE)
}


if(p==1){
	plot(midpoints1,results1,col=colors()[colorprof[p]],ylim=c(0,1),type="l",lwd=2,xlab="",ylab="")
	points(midpoints1,results1,col=colors()[colorprof[p]],pch=16)
	#abline(h=0.6473557,col='gray',lty=2)
	segments(midpoints1,lower1,midpoints1,upper1,col=colors()[colorprof[p]])
}else{
	lines(midpoints1,results1,col=colors()[colorprof[p]],ylim=c(0,1),type="l",lwd=2,xlab="",ylab="")
	points(midpoints1,results1,col=colors()[colorprof[p]],pch=16)
	#abline(h=0.6473557,col='gray',lty=2)
	segments(midpoints1,lower1,midpoints1,upper1,col=colors()[colorprof[p]])
}

}

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2a.pdf",sep=""),width=3.5,height=4)


#############
### Figure 2b
#############
rm(list=ls())
merge=read.table("ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt",header=T)
data=subset(merge, cov_input>5 & cov_input<quantile(merge$cov_input,0.999))
proms=subset(data,promoter_overlap==1 & DNaseHS_overlap==1 & HEK293T_H3K4me3_overlap==1 )
nonproms=subset(data,promoter_overlap==0 & DNaseHS_overlap==0 & anyENCODE_H3K4me3_overlap==0 )
proms[is.na(proms$which_motif)==T,"which_motif"]=0
nonproms[is.na(nonproms$which_motif)==T,"which_motif"]=0

TP=table(proms$which_motif)/sum(table(proms$which_motif))
#         0          1          2          3          4          5          6          7 
#0.40162540 0.12990613 0.02167202 0.01556102 0.04088704 0.03225603 0.08681409 0.27127827 

TNP= table(nonproms$which_motif)/sum(table(nonproms$which_motif))
#         0          1          2          3          4          5          6          7 
#0.34879171 0.16003021 0.02152272 0.02783880 0.11389537 0.04184402 0.14588768 0.14018948

colorprof = c(133,54,148,51,565,548,640)
barplot(matrix(c(rev(TNP),rev(TP)),nrow=8,ncol=2),col=rev(c('gray',colors()[colorprof])))

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2b.pdf",sep=""),width=2,height=4)


###############
### Figure 2c/d
###############

rm(list=ls())
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

data=read.table("ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt",header=T)
dataS=subset(data, cov_input>5 & cov_input<quantile(cov_input,0.999)  & is.na(HapMap_recombination_rate)==F & repeat_overlap==0 & telomere_overlap==0)
dataSprom=subset(dataS, promoter_overlap==1 & DNaseHS_overlap==1 & HEK293T_H3K4me3_overlap==1 )
dataSnoprom = subset(dataS, promoter_overlap==0 & DNaseHS_overlap==0 & anyENCODE_H3K4me3_overlap==0 )
mot1noprom = dataSnoprom
mot1prom = dataSprom

input1=mot1prom
input2=mot1noprom
cor(input2$enrichment,input2$HapMap_recombination_rate,use="complete.obs")
#0.2320102

bins1=quantile(input1$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins2=quantile(input2$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins2)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins2)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins2)-1)
upper2 = rep(NA,length(bins2)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05

for(i in 1:(length(bins1)-1)){
	sub1 = subset(input1, enrichment >=bins1[i] & enrichment <bins1[i+1])
	results1[i] = mean(sub1$HapMap_recombination_rate,na.rm=TRUE)
	seprop1 = stderrmean(sub1$HapMap_recombination_rate)
	lower1[i] = results1[i] - 2*seprop1
	upper1[i] = results1[i] + 2*seprop1
	midpoints1[i] = median(sub1$enrichment,na.rm=TRUE)
	
	sub2 = subset(input2, enrichment >=bins2[i] & enrichment <bins2[i+1])
	results2[i] = mean(sub2$HapMap_recombination_rate,na.rm=TRUE)
	seprop2 = stderrmean(sub2$HapMap_recombination_rate)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints2[i] = median(sub2$enrichment,na.rm=TRUE)
}

plot(midpoints2,results2,xlab="PRDM9 ChIP-seq enrich.",ylab="Mean Rec. rate (cM/Mb)",col=colors()[188],pch=16,xlim=c(0,max(midpoints2)),ylim=c(0,max(upper2)))
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[188])
points(midpoints1,results1,col=colors()[525],pch=15)
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[525])
legend("topleft",legend=c("Non-Promoter","Promoter"),col=c(colors()[188],colors()[525]),pch=c(16,15))

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2c.pdf",sep=""),width=3.5,height=4)

#mot1promsub=subset(mot1prom,enrich>1 & enrich<2)
#mot1nopromsub=subset(mot1noprom,enrich>1 & enrich<2)
#write.table(cbind(paste(mot1promsub[,1],sep=""),mot1promsub[,2]-150+16+mot1promsub[,25]-1,mot1promsub[,2]-150+16+mot1promsub[,25],"+"),file="MotifALL.Promoters.Matched.bed",quote=F,sep="\t",row.names=F,col.names=F)
#write.table(cbind(paste(mot1nopromsub[,1],sep=""),mot1nopromsub[,2]-150+16+mot1nopromsub[,25]-1,mot1nopromsub[,2]-150+16+mot1nopromsub[,25],"+"),file="MotifALL.NonPromoters.Matched.bed",quote=F,sep="\t",row.names=F,col.names=F)
#perl MakeProfilePlot.pl MotifALL.Promoters.Matched.bed HumanRecMap.HapMapCEU.CMperMb.bed ProfilePlot.Recrate.YFP_HumanPRDM9Peaks.MotifALL.Promoter.Matched.bed 10000 1 8 1
#perl MakeProfilePlot.pl MotifALL.NonPromoters.Matched.bed HumanRecMap.HapMapCEU.CMperMb.bed ProfilePlot.Recrate.YFP_HumanPRDM9Peaks.MotifALL.NonPromoter.Matched.bed 10000 1 8 1

prof=read.table(paste("ProfilePlot.Recrate.YFP_HumanPRDM9Peaks.MotifALL.Promoter.Matched.bed",sep=""),header=T)
prof2=read.table(paste("ProfilePlot.Recrate.YFP_HumanPRDM9Peaks.MotifALL.NonPromoter.Matched.bed",sep=""),header=T)

smoothed1 = ksmooth(prof$Position,prof$Mean, "normal", bandwidth=10)
smoothed1x=smoothed1$x[seq(1,20001,500)]
smoothed1y=smoothed1$y[seq(1,20001,500)]
smoothed2 = ksmooth(prof2$Position,prof2$Mean, "normal", bandwidth=10)
smoothed2x=smoothed2$x[seq(1,20001,500)]
smoothed2y=smoothed2$y[seq(1,20001,500)]

plot(smoothed1x,smoothed1y,type='l',lty=3,col=colors()[525],ylim=c(0.7,4.9),main="",lwd=3,xlab="distance from motif center (bp)",ylab="mean rec. rate (cM/Mb)")
lines(smoothed2,col=colors()[188],ylim=c(0,5.5),main="",lwd=3,xlab="distance from motif center (bp)",ylab="mean rec. rate (cM/Mb)")
#legend("topleft",legend=c("Non-Promoter (enrich. 1-2)","Promoter (enrich. 1-2)"),col=c(colors()[188],colors()[525]),lty=c(1,3),lwd=3,cex=1,bty='n')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2d.pdf",sep=""),width=3.5,height=4)



#################
### Figure 2e/f/g
#################
rm(list=ls())
load("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q1.all.r")
finalMeansH=finalMeans
load("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q1.all.r")
finalMeansU=finalMeans

colorprof=c(468,469,470,471)
plot(ksmooth(seq(-5000,5000,1),finalMeansH/finalMeansU, "normal", bandwidth=200),type='l',col=colors()[colorprof[1]],xlab="",ylab="",main="",lwd=3,ylim=c(0.75,1.9))

for(i in 2:4){
	load(paste0("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q",i,".all.r"))
	finalMeansH=finalMeans
	load(paste0("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q",i,".all.r"))
	finalMeansU=finalMeans

	lines(ksmooth(seq(-5000,5000,1),finalMeansH/finalMeansU, "normal", bandwidth=200),type='l',col=colors()[colorprof[i]],xlab="",ylab="",main="",lwd=3)
}

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2e.pdf",sep=""),width=3.5,height=4)


load("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.Promoters.q1.all.r")
finalMeansH=finalMeans
load("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.Promoters.q1.all.r")
finalMeansU=finalMeans

colorprof=c(468,469,470,471)
plot(ksmooth(seq(-5000,5000,1),finalMeansH/finalMeansU, "normal", bandwidth=200),type='l',col=colors()[colorprof[1]],xlab="",ylab="",main="",lwd=3,ylim=c(0.5,3.1))

for(i in 2:4){
	load(paste0("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.Promoters.q",i,".all.r"))
	finalMeansH=finalMeans
	load(paste0("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.Promoters.q",i,".all.r"))
	finalMeansU=finalMeans

	lines(ksmooth(seq(-5000,5000,1),finalMeansH/finalMeansU, "normal", bandwidth=200),type='l',col=colors()[colorprof[i]],xlab="",ylab="",main="",lwd=3)
}

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2f.pdf",sep=""),width=3.5,height=4)


load("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q4.all.r")
finalMeansH=finalMeans
load("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.NonPromoters.q4.all.r")
finalMeansU=finalMeans

smoothedH = ksmooth(seq(-5000,5000,1),finalMeansH, "normal", bandwidth=200)
smoothedHx=smoothedH$x[seq(1,10001,100)]
smoothedHy=smoothedH$y[seq(1,10001,100)]
smoothedU = ksmooth(seq(-5000,5000,1),finalMeansU, "normal", bandwidth=200)
smoothedUx=smoothedU$x[seq(1,10001,100)]
smoothedUy=smoothedU$y[seq(1,10001,100)]

colorprof=c(468,469,470,471)
plot(smoothedHx,smoothedHy,type='l',col=colors()[188],xlab="",ylab="",main="",lwd=3,ylim=c(0,0.5))
lines(smoothedUx,smoothedUy,type='l',col=colors()[188],xlab="",ylab="",main="",lwd=3,lty=3)


load("finalMeansProfile.HumanK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.Promoters.q4.all.r")
finalMeansH=finalMeans
load("finalMeansProfile.UTK36enrichHumanPeakCentres.MotifCentred.MotifsOnly.Promoters.q4.all.r")
finalMeansU=finalMeans

smoothedH = ksmooth(seq(-5000,5000,1),finalMeansH, "normal", bandwidth=200)
smoothedHx=smoothedH$x[seq(1,10001,100)]
smoothedHy=smoothedH$y[seq(1,10001,100)]
smoothedU = ksmooth(seq(-5000,5000,1),finalMeansU, "normal", bandwidth=200)
smoothedUx=smoothedU$x[seq(1,10001,100)]
smoothedUy=smoothedU$y[seq(1,10001,100)]

colorprof=c(468,469,470,471)
lines(smoothedHx,smoothedHy,type='l',col=colors()[525],xlab="",ylab="",main="",lwd=3,ylim=c(0,0.35))
lines(smoothedUx,smoothedUy,type='l',col=colors()[525],xlab="",ylab="",main="",lwd=3,lty=3)

#legend("topright",legend=c("NonProm., Transf.","NonProm., Untransf.","Prom., Transf.","Prom., Untransf."),col=c(colors()[188],colors()[188],colors()[525],colors()[525]),lty=c(1,3,1,3),lwd=3,cex=0.75,bty='n')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2g.pdf",sep=""),width=3.5,height=4)


#############
### Figure 2h
#############

rm(list=ls())
data=read.table("ForceCall.HumanK4onTSSpositions.autosomal.plusFPKM.plusPRDM9enrich.bed",header=T)
data=data[order(data$gene),]
data2=read.table("ForceCall.UTK4onTSSpositions.autosomal.plusFPKMv2.bed",header=T)
data2=data2[order(data2$gene),]
data2$PRDM9enrich = data$PRDM9enrich

sub = subset(data2,is.na(H.pval)==F & H.cov.g>5 & H.pval<=1 & PRDM9enrich>2 & PRDM9enrich<10)
rownames(sub)=paste(sub[,1],sub[,2],sub[,4],sep=":")
sub[is.na(sub$PRDM9enrich)==T,"PRDM9enrich"]=0

data3 = read.table("ForceCall.UTK36onTSSpositions.autosomal.bed",header=T)
rownames(data3)=paste(data3[,1],data3[,2],data3[,4],sep=":")
data3sub = data3[rownames(sub),]

dim(data3sub)
#[1] 4641   28
dim(sub)
#[1] 4641   28

data4 = read.table("ForceCall.HumanK36onTSSpositions.autosomal.bed",header=T)
rownames(data4)=paste(data4[,1],data4[,2],data4[,4],sep=":")
data4sub = data4[rownames(sub),]

data3sub$PRDM9enrich1 = sub$PRDM9enrich
data4sub$PRDM9enrich1 = sub$PRDM9enrich
data3sub$PRDM9enrich = sub$H.enrichment
data4sub$PRDM9enrich = sub$H.enrichment

inputA = subset(sub,is.na(data3sub$H.enrichment)==F & is.na(data4sub$H.enrichment)==F)
inputB = subset(data3sub,is.na(data3sub$H.enrichment)==F & is.na(data4sub$H.enrichment)==F)
inputC = subset(data4sub,is.na(data3sub$H.enrichment)==F & is.na(data4sub$H.enrichment)==F)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
colorprof=c(468,469,470,471)

for(p in 1:2){
	input1 = inputB
	if(p==2){
		input1=inputC
	}

	bins1=quantile(input1$PRDM9enrich,seq(0,1,0.2),na.rm=TRUE)
	bins1[length(bins1)]=bins1[length(bins1)]+1
	results1 = rep(0,length(bins1)-1)
	midpoints1 = rep(0,length(bins1)-1)
	lower1 = rep(NA,length(bins1)-1)
	upper1 = rep(NA,length(bins1)-1)

	for(i in 1:(length(bins1)-1)){
		sub1 = subset(input1, PRDM9enrich >=bins1[i] & PRDM9enrich <bins1[i+1])
		print(i)
		print(mean(sub1$PRDM9enrich1))
		results1[i] = mean(sub1$H.enrichment)
		seprop1 = stderr(sub1$H.enrichment)
		lower1[i] = results1[i] - 2*seprop1
		upper1[i] = results1[i] + 2*seprop1
		midpoints1[i] = median(sub1$PRDM9enrich,na.rm=TRUE)
	}

	if(p==1){
		plot(midpoints1,results1,col=colors()[525],type="l",lwd=3,xlab="",ylab="",ylim=c(0,0.17),lty=3)
		points(midpoints1,results1,col=colors()[525],pch=16)
		segments(midpoints1,lower1,midpoints1,upper1,col=colors()[525],lwd=2)
	}else{
		lines(midpoints1,results1,col=colors()[525],type="l",lwd=3,xlab="",ylab="")
		points(midpoints1,results1,col=colors()[525],pch=16)
		segments(midpoints1,lower1,midpoints1,upper1,col=colors()[525],lwd=2)
	}
}
#legend("topright",legend=c("Transfected","Untransfected"),col=c(colors()[525]),lty=c(1,3),lwd=3,cex=1,bty='n')

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2h.pdf",sep=""),width=3.5,height=4)

##################################
### Figure 2- Figure supplement 1a
##################################
rm(list=ls())
human=read.table("ChIPseq_Peaks.HumanPRDM9_HAV5.antiHAV5.protocolC.p10e-5.sep250.txt",header=T)
chimp=read.table("ChIPseq_Peaks.ChimpPRDM9_HAV5.antiHAV5.protocolC.p10e-5.sep250.Annotated.txt",header=T)
human1= subset(human,chr=="chr1")
chimp1=subset(chimp,chr=="chr1")

hist(human1[,2],50,col=rgb(0,0,1,0.5),freq=F,xlab="position on chr1",main="Distribution of Peaks on Chr1",border=rgb(0,0,0,0))
hist(chimp1[,2],50,col=rgb(0,1,0,0.5),freq=F,add=T,border=rgb(0,0,0,0))
#legend("topright",legend=c("Human Peaks","Chimp Peaks"),col=c('blue','green'),pch=15,bty=F)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S1a.pdf",sep=""),width=6,height=4)


##################################
### Figure 2- Figure supplement 1b
##################################
rm(list=ls())
#system("awk -v OFS="\t" 'FNR>1 {print $1,$2,$3}' ChIPseq_Peaks.HumanPRDM9_HAV5.antiHAV5.protocolC.p10e-5.sep250.txt > HumanPeakCentres.bed")
#system("awk -v OFS="\t" 'FNR>1 {print $1,$2,$3}' ChIPseq_Peaks.ChimpPRDM9_HAV5.antiHAV5.protocolC.p10e-5.sep250.Annotated.txt > ChimpPeakCentres.bed")
#system("bedtools slop -b 500 -i HumanPeakCentres.bed -g hg19.chromsizes.tbl > HumanPeakCentres.slop500.bed")
#system("bedtools slop -b 500 -i ChimpPeakCentres.bed -g hg19.chromsizes.tbl > ChimpPeakCentres.slop500.bed")
#
#filea = "HumanPeakCentres.slop500.bed"
#fileb = "ChimpPeakCentres.bed"
#pratto=read.table(filea,header=F)
#widths = pratto[,3]-pratto[,2]
#gsize = 2634892424
#tot=length(widths)
#peaknum=as.integer((strsplit(system(paste("wc -l",fileb),intern=TRUE)," ")[[1]])[1])
#found = as.integer(system(paste("bedtools intersect -u -a",filea,"-b",fileb,"| wc -l"),intern=TRUE))
#chance=sum(1-exp(peaknum*(log(gsize-widths)-log(gsize))))
#tot
##213885
#peaknum
##247717
#found 
##13423
#chance
##19210.22
#found/peaknum
##0.05418683
#chance/peaknum
##0.07754906
#
#system("bedtools intersect -u -a ChimpPeakCentres.bed -b HumanPeakCentres.slop500.bed >ChimpPeakCentres.HumanOverlap.bed")

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
stderrprop <- function(x,n) sqrt(x*(1-x)/n)

data=read.table("ChIPseq_Peaks.ChimpPRDM9_HAV5.antiHAV5.protocolC.p10e-5.sep250.Annotated.txt",header=T)
rownames(data)=paste0(data[,1],":",data[,2])
data2 = read.table("ChimpPeakCentres.HumanOverlap.bed",header=T)
rownames(data2)=paste0(data2[,1],":",data2[,2])
data$Overlap=0
data[rownames(data2),"Overlap"]=1

bins1=quantile(data[,9],seq(0,1,0.1),na.rm=TRUE)
bins1[11]=bins1[11]+1
props=rep(0,10)
mids=rep(0,10)
upper=rep(0,10)
lower=rep(0,10)
for(i in 1:(length(bins1)-1)){
	sub2 = subset(data, data[,9] >=bins1[i] & data[,9] <bins1[i+1])
	widths = sub2[,3]-sub2[,2]
	gsize = 2634892424
	tot=length(widths)
	peaknum=247717
	found = sum(sub2$Overlap)
	corrected = found/tot
	props[i]=corrected
	mids[i]=median(sub2[,9])
	upper[i]=props[i]+2*stderrprop(corrected,dim(sub2)[1])
	lower[i]=props[i]-2*stderrprop(corrected,dim(sub2)[1])
}

plot(mids,props,main="",xlim=c(0,max(mids)),ylim=c(0,max(upper)),col='darkgreen',pch=16,xlab="",ylab="")
segments(mids,lower,mids,upper,col='darkgreen')
abline(h=0.07754906,col='gray',lty=2)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S1b.pdf",sep=""),width=3.5,height=4)


##################################
### Figure 2- Figure supplement 1c
##################################

#see Rscript FindMotifs_Chimp.R for code used to generate chimp motif PWM

#Get PWMs from all motifs, stored in intermediate file "Chimp_Motif_Results_Final_iter65.r"
load("Chimp_Motif_Results_Final_iter65.r")
library(gtools)
library(seqLogo)
mySeqLogo = seqLogo::seqLogo 
bad = (sapply( body(mySeqLogo), "==", "grid.newpage()") | sapply( body(mySeqLogo), "==", "par(ask = FALSE)"))
body(mySeqLogo)[bad] = NULL 
motnum=length(znew$scorematdim)
plotting=TRUE
dimvec = znew$scorematdim
scorematset = znew$scoremat
alpha= znew$alpha
iter=65
compvec=c(1,1) #specify whether to take reverse complement of each motif for final output (0=no, 1=yes)
keep=c(2,1) #specify the final ordering of motif names

#plot motifs and output individual PWM text files for each motif
system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
props = rep(0,motnum)
for(i in 1:motnum){
	newstarts=c(1,cumsum(dimvec)+1)
	newends=cumsum(dimvec)
	comp=compvec[i]

	testmat=scorematset[newstarts[i]:newends[i],]
	testmat=exp(testmat)
	compmat=testmat[,c(4:1)]
	testmat=testmat/rowSums(testmat)
	compmat=compmat/rowSums(compmat)
	pwm = t(testmat)
	compmat=compmat[nrow(compmat):1,]
	grid.newpage()
	name=which(keep==i)
	if(comp==1){
		mySeqLogo(t(compmat))
		#write.table(compmat,file=paste("plots/chimpPWM.Iter",iter,".Motif",name,".tbl",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
	}else{
		mySeqLogo(pwm)
		#write.table(t(pwm),file=paste("plots/chimpPWM.Iter",iter,".Motif",name,".tbl",sep=""),row.names=F,col.names=F,quote=F,sep="\t")
	}
	dev.copy2pdf(file=paste("plots/LOGO_Chimp_Motif",name,".pdf",sep=""),width=36,height=24)
	
	if(name==1){
		dev.copy2pdf(file="plots/Figure2_S1c.pdf",width=36,height=24)
	}
	
	hist(znew$whichpos[znew$whichmot==i],breaks=seq(0,300,25),col='blue',xlab="",ylab="",main="")
	dev.copy2pdf(file=paste("plots/Histogram_Chimp_Motif",name,".pdf",sep=""),width=36,height=24)
	props[name]=mean(znew$whichpos[znew$whichmot==i]>100 & znew$whichpos[znew$whichmot==i]<200)
}
print(props)

#other panels of Figure2-S1c generated using online ZF binding motif prediction tool
#at http://compbio.cs.princeton.edu/zf/form.html


#################################
### Figure 2-Figure supplement 2a
#################################
rm(list=ls())
data=read.table("fimo.Human1.sorted.peakcentreoverlaps.300bp.bed")
bins1=seq(12,36,1)
bins1=quantile(data[,5],seq(0,1,0.001),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
results1 = c(0,length(bins1)-1)
midpoints1 = c(0,length(bins1)-1)
for(i in 1:(length(bins1)-1)){
	sub = subset(data, V5 >=bins1[i] & V5 <bins1[i+1])
	prop = mean(sub$V7,na.rm=T)
	corrected = prop #(1-((1-(prop))/(1-(0.04373913))))
	results1[i] = corrected
	midpoints1[i] = median(sub$V5,na.rm=TRUE)	
}
plot(midpoints1,results1,xlab="",ylab="",col=rgb(0,0,1,0.5),pch=16,main="",ylim=c(0,0.5))

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S2a.pdf",sep=""),width=3.5,height=4)


#####################################
### Figure 2-Figure supplement 2b/c/d
#####################################
rm(list=ls())
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
stderrprop <- function(x,n) sqrt(x*(1-x)/n)

merge=read.table("ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt",header=T)
data=subset(merge, cov_input>5 & cov_input<quantile(cov_input,0.999) & which_motif==1 & is.na(FIMO_score)==F)
proms=subset(data, promoter_overlap==1 & DNaseHS_overlap==1 & HEK293T_H3K4me3_overlap==1 )
nonproms = subset(data, promoter_overlap==0 & DNaseHS_overlap==0 & anyENCODE_H3K4me3_overlap==0 )


bins1=quantile(proms$FIMO_score,seq(0,1,0.1),na.rm=TRUE)
bins2=quantile(nonproms$FIMO_score,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)

for(i in 1:(length(bins1)-1)){
	sub1 = subset(proms, FIMO_score >=bins1[i] & FIMO_score <bins1[i+1])
	sub2 = subset(nonproms, FIMO_score >=bins2[i] & FIMO_score <bins2[i+1])
	
	results1[i] = mean(sub1$cov_input,na.rm=TRUE)
	results2[i] = mean(sub2$cov_input,na.rm=TRUE)
	se1 = stderr(sub1$cov_input)
	se2 = stderr(sub2$cov_input)
	lower1[i] = results1[i] - 2*se1
	upper1[i] = results1[i] + 2*se1
	lower2[i] = results2[i] - 2*se2
	upper2[i] = results2[i] + 2*se2
	midpoints1[i] = median(sub1$FIMO_score,na.rm=TRUE)
	midpoints2[i] = median(sub2$FIMO_score,na.rm=TRUE)
}
plot(midpoints2,results2,xlab="",ylab="",col=colors()[188],pch=16,main="",xlim=c(12,max(midpoints2)),ylim=c(0,50))
points(midpoints1,results1,col=colors()[525],pch=16)
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[188])
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[525])
#legend("topleft",legend=c("Non-promoter","Promoter"),col=c(colors()[188],colors()[525]),pch=16)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S2b.pdf",sep=""),width=3.5,height=4)



bins1=quantile(proms$FIMO_score,seq(0,1,0.1),na.rm=TRUE)
bins2=quantile(nonproms$FIMO_score,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)

for(i in 1:(length(bins1)-1)){
	sub1 = subset(proms, FIMO_score >=bins1[i] & FIMO_score <bins1[i+1])
	sub2 = subset(nonproms, FIMO_score >=bins2[i] & FIMO_score <bins2[i+1])
	
	results1[i] = mean(sub1$cov_r1+sub1$cov_r2,na.rm=TRUE)
	results2[i] = mean(sub2$cov_r1+sub2$cov_r2,na.rm=TRUE)
	se1 = stderr(sub1$cov_r1+sub1$cov_r2)
	se2 = stderr(sub2$cov_r1+sub2$cov_r2)
	lower1[i] = results1[i] - 2*se1
	upper1[i] = results1[i] + 2*se1
	lower2[i] = results2[i] - 2*se2
	upper2[i] = results2[i] + 2*se2
	midpoints1[i] = median(sub1$FIMO_score,na.rm=TRUE)
	midpoints2[i] = median(sub2$FIMO_score,na.rm=TRUE)
}
plot(midpoints2,results2,xlab="",ylab="",col=colors()[188],pch=16,main="",xlim=c(12,max(midpoints2)),ylim=c(60,200))
points(midpoints1,results1,col=colors()[525],pch=16)
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[188])
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[525])
#legend("topleft",legend=c("Non-promoter","Promoter"),col=c(colors()[188],colors()[525]),pch=16)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S2c.pdf",sep=""),width=3.5,height=4)


bins1=quantile(proms$FIMO_score,seq(0,1,0.1),na.rm=TRUE)
bins2=quantile(nonproms$FIMO_score,seq(0,1,0.1),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins1)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins1)-1)
upper2 = rep(NA,length(bins1)-1)

for(i in 1:(length(bins1)-1)){
	sub1 = subset(proms, FIMO_score >=bins1[i] & FIMO_score <bins1[i+1])
	sub2 = subset(nonproms, FIMO_score >=bins2[i] & FIMO_score <bins2[i+1])
	
	results1[i] = mean(sub1$enrichment,na.rm=TRUE)
	results2[i] = mean(sub2$enrichment,na.rm=TRUE)
	se1 = stderr(sub1$enrichment)
	se2 = stderr(sub2$enrichment)
	lower1[i] = results1[i] - 2*se1
	upper1[i] = results1[i] + 2*se1
	lower2[i] = results2[i] - 2*se2
	upper2[i] = results2[i] + 2*se2
	midpoints1[i] = median(sub1$FIMO_score,na.rm=TRUE)
	midpoints2[i] = median(sub2$FIMO_score,na.rm=TRUE)
}

plot(midpoints2,results2,xlab="",ylab="",col=colors()[188],pch=16,main="",xlim=c(12,max(midpoints2)),ylim=c(0,max(upper2)))
points(midpoints1,results1,col=colors()[525],pch=16)
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[188])
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[525])
#legend("topleft",legend=c("Non-promoter","Promoter"),col=c(colors()[188],colors()[525]),pch=16)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S2d.pdf",sep=""),width=3.5,height=4)



###################################
### Figure 2-Figure supplement 2e/f
###################################

rm(list=ls())
stderrprop <- function(x,n) sqrt(x*(1-x)/n)
stderrmean <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

data=read.table("ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt",header=T)
dataS=subset(data, cov_input>5 & cov_input<quantile(cov_input,0.999)  & is.na(HapMap_recombination_rate)==F & repeat_overlap==0 & telomere_overlap==0)
rownames(dataS)=paste0(dataS[,1],":",dataS[,2])
dataSprom=subset(dataS, promoter_overlap==1 & DNaseHS_overlap==1 & HEK293T_H3K4me3_overlap==1 )
dataSnoprom = subset(dataS, promoter_overlap==0 & DNaseHS_overlap==0 & anyENCODE_H3K4me3_overlap==0 )
mot1noprom = dataSnoprom
mot1prom = dataSprom

dmc1data = read.table("ForceCall.Pratto_DMC1_on_YFP_HumanPRDM9.autosomal.bed",header=T)
rownames(dmc1data)=paste0(dmc1data[,1],":",dmc1data[,2])
mot1noprom$dmc1 = NA
mot1prom$dmc1 = NA

dmc1datanoprom = dmc1data[rownames(mot1noprom),]
dmc1dataprom = dmc1data[rownames(mot1prom),]

mot1noprom[rownames(dmc1datanoprom),"dmc1"]= dmc1datanoprom$H.enrichment
mot1prom[rownames(dmc1dataprom),"dmc1"]= dmc1dataprom$H.enrichment
sum(is.na(mot1prom$dmc1) == FALSE)
#6367
sum(is.na(mot1noprom$dmc1) == FALSE)
#6615

mot1noprom = subset(mot1noprom, is.na(dmc1) == FALSE)
mot1prom = subset(mot1prom, is.na(dmc1) == FALSE)

input1=mot1prom
input2=mot1noprom
cor(input2$enrich,input2$dmc1,use="complete.obs")
#0.373

bins1=quantile(input1$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins2=quantile(input2$enrichment,seq(0,1,0.25),na.rm=TRUE)
bins1[length(bins1)]=bins1[length(bins1)]+1
bins2[length(bins2)]=bins2[length(bins2)]+1
results1 = rep(0,length(bins1)-1)
results2 = rep(0,length(bins2)-1)
results3 = rep(0,length(bins1)-1)
midpoints1 = rep(0,length(bins1)-1)
midpoints2 = rep(0,length(bins2)-1)
lower1 = rep(NA,length(bins1)-1)
upper1 = rep(NA,length(bins1)-1)
lower2 = rep(NA,length(bins2)-1)
upper2 = rep(NA,length(bins2)-1)
lower3 = rep(NA,length(bins1)-1)
upper3 = rep(NA,length(bins1)-1)
pthresh=0.05
for(i in 1:(length(bins1)-1)){
	sub1 = subset(input1, enrichment >=bins1[i] & enrichment <bins1[i+1])
	results1[i] = mean(sub1$dmc1,na.rm=TRUE)
	seprop1 = stderrmean(sub1$dmc1)
	lower1[i] = results1[i] - 2*seprop1
	upper1[i] = results1[i] + 2*seprop1
	midpoints1[i] = median(sub1$enrichment,na.rm=TRUE)
	
	sub2 = subset(input2, enrichment >=bins2[i] & enrichment <bins2[i+1])
	results2[i] = mean(sub2$dmc1,na.rm=TRUE)
	seprop2 = stderrmean(sub2$dmc1)
	lower2[i] = results2[i] - 2*seprop2
	upper2[i] = results2[i] + 2*seprop2
	midpoints2[i] = median(sub2$enrichment,na.rm=TRUE)
}

plot(midpoints2,results2,xlab="PRDM9 ChIP-seq enrich.",ylab="DMC1 enrich.",col=colors()[188],pch=16,xlim=c(0,max(midpoints2)),ylim=c(0,max(upper2)))
segments(midpoints2,lower2,midpoints2,upper2,col=colors()[188])
points(midpoints1,results1,col=colors()[525],pch=15)
segments(midpoints1,lower1,midpoints1,upper1,col=colors()[525])
legend("topleft",legend=c("Non-Promoter","Promoter"),col=c(colors()[188],colors()[525]),pch=c(16,15))

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S2e.pdf",sep=""),width=3.5,height=4)

#mot1promsub=subset(mot1prom,enrich>1 & enrich<2)
#mot1nopromsub=subset(mot1noprom,enrich>1 & enrich<2)
#write.table(cbind(paste(mot1promsub[,1],sep=""),mot1promsub[,2]-150+16+mot1promsub[,25]-1,mot1promsub[,2]-150+16+mot1promsub[,25],"+"),file="MotifALL.Promoters.Matched.bed",quote=F,sep="\t",row.names=F,col.names=F)
#write.table(cbind(paste(mot1nopromsub[,1],sep=""),mot1nopromsub[,2]-150+16+mot1nopromsub[,25]-1,mot1nopromsub[,2]-150+16+mot1nopromsub[,25],"+"),file="MotifALL.NonPromoters.Matched.bed",quote=F,sep="\t",row.names=F,col.names=F)
#system("bedtools genomecov -bg -i FragPos.AB1.bed -g hg19.chromsizes.tbl > FragDepth.AB1.bedgraph")
#system("perl MakeProfilePlot.pl MotifALL.Promoters.Matched.bed FragDepth.AB1.bedgraph ProfilePlot.DMC1_AB1.OldHumanPRDM9Peaks.MotifALL.Promoter.Matched.bed 10000 1 8 1")
#system("perl MakeProfilePlot.pl MotifALL.NonPromoters.Matched.bed FragDepth.AB1.bedgraph ProfilePlot.DMC1_AB1.OldHumanPRDM9Peaks.MotifALL.NonPromoter.Matched.bed 10000 1 8 1")

prof=read.table(paste("ProfilePlot.DMC1_AB1.OldHumanPRDM9Peaks.MotifALL.Promoter.Matched.bed",sep=""),header=T)
prof2=read.table(paste("ProfilePlot.DMC1_AB1.OldHumanPRDM9Peaks.MotifALL.NonPromoter.Matched.bed",sep=""),header=T)

smoothed1 = ksmooth(prof$Position,prof$Mean, "normal", bandwidth=500)
smoothed1x=smoothed1$x[seq(1,20001,250)]
smoothed1y=smoothed1$y[seq(1,20001,250)]
smoothed2 = ksmooth(prof2$Position,prof2$Mean, "normal", bandwidth=500)
smoothed2x=smoothed2$x[seq(1,20001,250)]
smoothed2y=smoothed2$y[seq(1,20001,250)]

plot(smoothed1x,smoothed1y,type='l',lty=3,col=colors()[525],ylim=c(1.8,3.4),main="",lwd=3,xlab="distance from motif center (bp)",ylab="mean DMC1 coverage (reads/bp)")
lines(smoothed2,type='l',col=colors()[188],main="",lwd=3,xlab="distance from Motif Center (bp)")
legend("topleft",legend=c("Non-Promoter (enrich. 1-2)","Promoter (enrich. 1-2)"),col=c(colors()[188],colors()[525]),lty=c(1,3),lwd=2, cex=0.75)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S2f.pdf",sep=""),width=3.5,height=4)


################################
### Figure 2-Figure supplement 3
################################

data0=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.Humanatac0nuccentre.bed",header=TRUE)
data1=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.Humanatac1nuccentre.bed",header=TRUE)
plot(ksmooth(data0$Position,data0$Mean/mean(data0$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='gray',xlim=c(-1000,1000),main="",ylab="",lwd=3,xlab="",ylim=c(0.9,2.4))
points(ksmooth(data1$Position,data1$Mean/mean(data1$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='blue',lwd=3)
#legend("topright",legend=c("Nuc.-free","MonoNuc."),col=c('gray','blue'),lty=1,lwd=3,bty="n")
system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S3a.pdf",sep=""),width=3.5,height=4)

data0=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.ZFonlyatac0nuccentre.bed",header=TRUE)
data1=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.ZFonlyatac1nuccentre.bed",header=TRUE)
plot(ksmooth(data0$Position,data0$Mean/mean(data0$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='gray',xlim=c(-1000,1000),main="",xlab="",ylab="",lwd=3,ylim=c(0.9,2.4))
points(ksmooth(data1$Position,data1$Mean/mean(data1$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='blue',lwd=3)
#legend("topright",legend=c("Nuc.-free","MonoNuc."),col=c('gray','blue'),lty=1,lwd=3,bty="n")
system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S3b.pdf",sep=""),width=3.5,height=4)

data0=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.UTatac0nuccentre.bed",header=TRUE)
data1=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.UTatac1nuccentre.bed",header=TRUE)
plot(ksmooth(data0$Position,data0$Mean/mean(data0$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='gray',xlim=c(-1000,1000),main="",lwd=3,xlab="",ylab="",ylim=c(0.9,2.4))
points(ksmooth(data1$Position,data1$Mean/mean(data1$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='blue',lwd=3)
#legend("topright",legend=c("Nuc.-free","MonoNuc."),col=c('gray','blue'),lty=1,lwd=3,bty="n")
system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S3c.pdf",sep=""),width=3.5,height=4)

data0=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.noZFatac0nuccentre.bed",header=TRUE)
data1=read.table("ProfilePlotInfo.BoundHumanMotifs.NonDNAse.strand.noZFatac1nuccentre.bed",header=TRUE)
plot(ksmooth(data0$Position,data0$Mean/mean(data0$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='gray',xlim=c(-1000,1000),main="",xlab="",ylab="",lwd=3,ylim=c(0.9,2.4))
points(ksmooth(data1$Position,data1$Mean/mean(data1$Mean[1:1500],na.rm=T), "normal", bandwidth=50),type='l',col='blue',lwd=3)
#legend("topright",legend=c("Nuc.-free","MonoNuc."),col=c('gray','blue'),lty=1,lwd=3,bty="n")

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure2_S3d.pdf",sep=""),width=3.5,height=4)



#############
### Figure 3a
#############
rm(list=ls())
means=read.table("FPKM.means.txt",header=T,row.names=1)
lo=read.table("FPKM.lo.txt",header=T,row.names=1)
hi=read.table("FPKM.hi.txt",header=T,row.names=1)

relmeans=means
relmeans[2,]=means[2,]/means[1,]
relmeans[3,]=means[3,]/means[1,]
relmeans[4,]=means[4,]/means[1,]

rello=lo
rello[2,]=lo[2,]/hi[1,]
rello[3,]=lo[3,]/hi[1,]
rello[4,]=lo[4,]/hi[1,]

relhi=hi
relhi[2,]=hi[2,]/lo[1,]
relhi[3,]=hi[3,]/lo[1,]
relhi[4,]=hi[4,]/lo[1,]

colorprof5=c('indianred4','darkslateblue','lightsteelblue')

#bp=barplot(log(t(as.matrix(relmeans[2:4,])),base=2),density=c(15,40,-1),beside=T,col=colorprof5,legend.text=c("VCX","CTCFL","CTCF"),ylim=c(-3,8.5),ylab="log2 fold RPKM over untransfected",xlab="Transfection with",main="RNA-seq results",names=c("Human","ZFonly","Chimp"))
bp=barplot(log(t(as.matrix(relmeans[2:4,])),base=2),beside=T,col=colorprof5,ylim=c(-3,8.5),names=c("Human","ZFonly","Chimp"))
segments(as.vector(bp),log(as.vector(t(rello[2:4,])),base=2),as.vector(bp),log(as.vector(t(relhi[2:4,])),base=2))

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure3a_left.pdf",sep=""),width=3.5,height=4)


means=read.table("means.txt",header=T,row.names=1)
sdevs=read.table("sdevs.txt",header=T,row.names=1)

colorprof5=c('indianred4','darkslateblue','lightsteelblue')

bp=barplot(t(as.matrix(means[,1:3])),beside=T,col=colorprof5,ylim=c(-3,8.5),ylab="Fold expression over untransfected",xlab="Transfection with")
segments(as.vector(bp),as.vector(t(means[,1:3]))-as.vector(t(sdevs[,1:3])),as.vector(bp),as.vector(t(means[,1:3]))+as.vector(t(sdevs[,1:3])))

data = read.table("FinalValues.txt",header=T)
points(seq(1.25,1.75,0.25),data[2,7:9],col='darkgrey',cex=0.5)
points(seq(2.25,2.75,0.25),data[2,4:6],col='darkgrey',cex=0.5)
points(seq(3.25,3.75,0.25),data[2,1:3],col='darkgrey',cex=0.5)
points(seq(5.25,5.75,0.25),data[3,7:9],col='darkgrey',cex=0.5)
points(seq(6.25,6.75,0.25),data[3,4:6],col='darkgrey',cex=0.5)
points(seq(7.25,7.75,0.25),data[3,1:3],col='darkgrey',cex=0.5)
points(seq(9.25,9.75,0.25),data[1,7:9],col='darkgrey',cex=0.5)
points(seq(10.25,10.75,0.25),data[1,4:6],col='darkgrey',cex=0.5)
points(seq(11.25,11.75,0.25),data[1,1:3],col='darkgrey',cex=0.5)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure3a_right.pdf",sep=""),width=3.5,height=4)



#################################################
### Figure 5b/c and Figure 5-Figure supplement 3d
#################################################
rm(list=ls())
#coIP/input band intensities for each sample (computed using ImageLab software)
frac2 = c(0.51,0.45,0.04,0.08,0.01,0.26,0.12)
barplot(frac2,space=c(0,0,0,0,0,0,0),col=c('orange3','orange3','olivedrab','olivedrab','olivedrab','lightsteelblue3','lightsteelblue3'),horiz=F)
barplot(frac2[1:5],space=c(0.1,0.1,0.1,0.1,0.1),col=c('olivedrab','olivedrab','orange3','orange3','orange3','lightsteelblue3','lightsteelblue3'),horiz=F)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure5b.pdf",sep=""),width=4,height=3)

#band intensities of upper (chimp) and lower (human) bands in each lane (computed using ImageLab software)
upper = c(3.24,1.89,8.85,1.28)
lower = c(1,5.38,1.62,8.36)
frac = upper/lower
frac
#[1] 3.2400000 0.3513011 5.4629630 0.1531100
log(frac,2)
#[1]  1.695994 -1.509220  2.449684 -2.707359

logfrac = c(-1.695994,-2.449684,1.509220,2.707359) #swap 2 and 3
barplot(logfrac,space=c(0,0,0.5,0),col=c('deeppink4','darkslateblue','deeppink4','darkslateblue'),horiz=F)

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure5c.pdf",sep=""),width=4,height=3.5)

#band intensity values for chimp and human input/IP bands with and without benzonase treatment (computed using ImageLab software)
c.i.nb=1394354 #chimp, input, no benzonase
c.i.b=430800 #chimp, input, benzonase
c.IP.nb=2971000 #chimp, IP, no benzonase
c.IP.b=2021775 #chimp, IP, benzonase
h.i.nb=3322774 #human, input, no benzonase
h.i.b=1175544 #human, input, benzonase
h.IP.nb=12004550 #human, IP, no benzonase
h.IP.b=11351500 #human, IP, benzonase

frac = c((12004550/3322774)/(2971000/1394354),(11351500/1175544)/(2021775/430800))
logfrac = log(frac,2)
barplot(logfrac,space=c(0.5),col=c('deeppink4','deeppink4'),horiz=F,density=c(250,10))

system("mkdir plots",ignore.stdout = T, ignore.stderr = T)
dev.copy2pdf(file=paste("plots/Figure5_S3d.pdf",sep=""),width=3,height=3)

