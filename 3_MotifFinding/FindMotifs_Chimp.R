#FindMotifs_Chimp.R
#R code used to find motifs in chimp PRDM9 ChIP-seq peaks given peak file containing peak sequences (peak centre +- 150 bp)
#Nicolas Altemose
#2015
####If you use this program, please cite Altemose et al. eLife 2017
####This is free software shared in the hope it may be of use; no warranty is given or implied

library(gtools)
library(seqLogo)
source("denovoMotifFinder_functions.R")
fullpeakfile="ChIPseq_Peaks.ChimpPRDM9_HAV5.antiHAV5.protocolC.p10e-5.sep250.Annotated.txt"

#read in a list of peaks with 300bp sequences (here, in column 17, with repeats masked as lower case)
fullpeakinfo = read.table(fullpeakfile,header=T,comment.char="?")
fullpeakinfo[,17]=as.character(fullpeakinfo[,17]) #ensure column with sequence is stored as character class

#filter peaks to a stringent set with genomic input coverage between 1 and 10, no repeatmasked bases, pvalue <1e-10, enrichment >2, CI width <=50bp, and sum ChIP coverage >=30
ours=which(fullpeakinfo[,8]>=1 & fullpeakinfo[,8]<=10 & toupper(fullpeakinfo[,17])==fullpeakinfo[,17] & as.double(fullpeakinfo[,10])<=1e-10 & as.double(fullpeakinfo[,9])>2 & (fullpeakinfo[,5]-fullpeakinfo[,4])<=50 & as.double(fullpeakinfo[,6])+as.double(fullpeakinfo[,7])>=30)
ourseqs=fullpeakinfo[ours,17]

#sort peaks from largest to smallest enrichment values
ourscores=(as.double(fullpeakinfo[,9]))[ours]
ourseqs=ourseqs[order(-ourscores)]
ourscores=ourscores[order(-ourscores)]

#seed the algorithm by applying the "findamotif" function to the top 2000 peak seuences, with parameters below
#finds a sequence of length 10 enriched for exact matches among the peak sequences
#then refines this motif across all sequences
z=findamotif(ourseqs,ourscores,len=10,nits=100,ntries=1,n_for_refine=2000,outfile="motiftry1.pdf",prior=NULL)
save(z,ourscores,file="mot1try.out")

#re-run the algorithm on sequences with a posterior match probability <0.75 to any previous motif to find a new motif, up to 14 iterations or until no new matches are found
for(i in 2:20){
	if(length(attributes(z))){
		newseqs=ourseqs[z$regprobs<0.75]
		newscores=ourscores[z$regprobs<0.75]
		z=findamotif(newseqs,newscores,len=10,nits=100,ntries=1,n_for_refine=min(2000,length(newseqs)),outfile=paste("motiftry",i,".pdf",sep=""),prior=z$prior,updateprior=0)
		save(z,newscores,file=paste("mot",i,"try.out",sep=""))
	}
}

#read in final motif seeds
mot=matrix(nrow=0,ncol=4)
sizes=vector(length=0)
prior=c()
for(i in 1:20){
	outfile=paste("mot",i,"try.out",sep="")
	load(outfile)
	if(i==1) prior=z$prior
	mot=rbind(mot,z$scoremat)
	sizes=c(sizes,nrow(z$scoremat))
}

#jointly call and refine all final motifs at top 5000 peaks (most computationally intensive step)
znew=getmotifs(mot,sizes,ourseqs[1:length(ourseqs)],maxwidth=max(nchar(ourseqs)),alpha=rep(1/3,3),incprob=0.99999,maxits=2,plen=0.9,updatemot=1,updatealpha=1,updateprior=0,ourprior=prior,bg=-1,plotting=FALSE)
for(iter in 1:65){
	print(paste("Iteration Number",iter))
	mot=znew$scoremat
	sizes=znew$scorematdim
	znew=getmotifs(mot,sizes,ourseqs[1:5000],maxwidth=max(nchar(ourseqs)),alpha=znew$alpha,incprob=0.99999,maxits=2,plen=0.9,updatemot=1,updatealpha=1,updateprior=0,ourprior=znew$prior,bg=-1,plotting=FALSE)
	save(iter,znew,file="Chimp_Motif_Results.r")
}

#final results are stored in "Chimp_Motif_Results_Final_iter65.r"


###now select non-degenerate motifs to call at ALL peaks (not just top 5000)
rm(list=ls())
load("Chimp_Motif_Results_Final_iter65.r")
compvec=c(1,1) #specify whether to take reverse complement of each motif for final output (0=no, 1=yes)
keep=c(2) #specify the indices of motifs to keep (non-degenerate motifs)

prior=znew$prior
starts=c(1,cumsum(znew$scorematdim)+1)
starts=starts[1:length(znew$scorematdim)]
ends=cumsum(znew$scorematdim)
mot=matrix(nrow=0,ncol=4)
sizes=vector(length=0)
for(j in keep){
 	scoremat = znew$scoremat[starts[j]:ends[j],]
 	if(compvec[j]==1){
 		compmat=scoremat[,c(4:1)]
		compmat=compmat[nrow(compmat):1,]
		scoremat=compmat
	}
	mot=rbind(mot,scoremat)
	sizes=c(sizes,nrow(scoremat))
}
#save the final non-degenerate motif info
save(mot,sizes,prior,file="FinalMotif.Chimp.r")

rm(list=ls())
library(gtools)
library(seqLogo)
source("denovoMotifFinder_functions.R")

#load final non-degenerate motifs, use saved prior
load("FinalMotif.Chimp.r")
fullpeakfile="ChIPseq_Peaks.ChimpPRDM9_HAV5.antiHAV5.protocolC.p10e-5.sep250.Annotated.txt"

fullpeakinfo = read.table(fullpeakfile,header=T)
fullpeakinfo[,17]=as.character(fullpeakinfo[,17]) #ensure sequence column is stored as character class
ourseqs=toupper(fullpeakinfo[,17]) #revert all sequences to uppercase

#now search for non-degenerate motifs at ALL peak sequences using the same algorithm, but do not update motif PWMs
print(date())
znewALL=getmotifs(mot,sizes,ourseqs[1:length(ourseqs)],maxwidth=301,alpha=1,incprob=0.99999,maxits=2,plen=0.9,updatemot=0,updatealpha=1,updateprior=0,ourprior=prior,bg=-1,plotting=FALSE)
for(iter in 1:100){
	print(paste("Iteration Number",iter))
	znewALL=getmotifs(mot,sizes,ourseqs[1:length(ourseqs)],maxwidth=301,alpha=znewALL$alpha,incprob=0.99999,maxits=2,plen=0.9,updatemot=0,updatealpha=1,updateprior=0,ourprior=prior,bg=-1,plotting=FALSE)
	save(znewALL,file="Chimp_Motif_Results_Final_ALLpeaks.r")
	print(date())
}

#add in motif info and print out final results
fullpeakinfo = cbind(fullpeakinfo,znewALL$whichmot, znewALL$whichpos, znewALL$whichstrand)
write.table(fullpeakinfo,"ChimpPeaks.withMotifs.tbl",sep="\t",quote=F,row.names=F,col.names=F)







