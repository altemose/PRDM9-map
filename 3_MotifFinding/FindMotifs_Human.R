#FindMotifs_Human.R
#R code used to find motifs in human PRDM9 ChIP-seq peaks given peak file containing peak sequences (peak centre +- 150 bp)
#Nicolas Altemose
#2015
####If you use this program, please cite Altemose et al. eLife 2017
####This is free software shared in the hope it may be of use; no warranty is given or implied


library(gtools)
library(seqLogo)
source("denovoMotifFinder_functions.R")
fullpeakfile="ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt"

#read in a list of peaks with 300bp sequences (here, in column 26, with repeats masked as lower case)
fullpeakinfo = read.table(fullpeakfile,header=T,comment.char="?")
fullpeakinfo[,26]=as.character(fullpeakinfo[,26])

#filter peaks to a stringent set with genomic input coverage between 2 and 26 (10th and 90th %iles), no repeatmasked bases, pvalue <1e-10, enrichment >2, CI width <=50bp, and sum ChIP coverage >=50
ours=which(as.double(fullpeakinfo[,9])>2 & fullpeakinfo[,8]>2 & fullpeakinfo[,8]<=26 & as.double(fullpeakinfo[,10])<=1e-10 & toupper(fullpeakinfo[,26])==fullpeakinfo[,26] & (fullpeakinfo[,5]-fullpeakinfo[,4])<=50 & as.double(fullpeakinfo[,6])+as.double(fullpeakinfo[,7])>=50)
ourseqs=fullpeakinfo[ours,26]

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
znew=getmotifs(mot,sizes,ourseqs[1:length(ourseqs)],maxwidth=max(nchar(ourseqs)),alpha=rep(1/length(sizes),length(sizes)),incprob=0.99999,maxits=2,plen=0.9,updatemot=1,updatealpha=1,updateprior=0,ourprior=prior,bg=-1,plotting=FALSE)
for(iter in 1:240){
	print(paste("Iteration Number",iter))
	mot=znew$scoremat
	sizes=znew$scorematdim
	znew=getmotifs(mot,sizes,ourseqs[1:5000],maxwidth=max(nchar(ourseqs)),alpha=znew$alpha,incprob=0.99999,maxits=2,plen=0.9,updatemot=1,updatealpha=1,updateprior=0,ourprior=znew$prior,bg=-1,plotting=FALSE)
	save(iter,znew,file="Human_Motif_Results.r")
}

#final results from published run are stored in "Human_Motif_Results_Final_iter240.r"


###now select non-degenerate motifs to call at ALL peaks (not just top 5000)
rm(list=ls())
load("Human_Motif_Results_Final_iter240.r")
compvec=c(0,1,0,0,0,0,1,1,0,0,1,1,0,0,0,1,1) #specify whether to take reverse complement of each motif for final output (0=no, 1=yes)
keep=c(4,6,8,9,10,15,17,5,11,13,7,16,2,3,12,14,1) #specify the indices of motifs to keep (non-degenerate motifs)
keep=c(4,6,8,9,10,15,17)


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
save(mot,sizes,prior,file="FinalMotifs.Human.r")


rm(list=ls())
library(gtools)
library(seqLogo)
source("denovoMotifFinder_functions.R")

#load final non-degenerate motifs, use saved prior
load("FinalMotifs.Human.r")
fullpeakfile="ChIPseq_Peaks.YFP_HumanPRDM9.antiGFP.protocolN.p10e-5.sep250.Annotated.txt"

fullpeakinfo = read.table(fullpeakfile,header=T)
fullpeakinfo[,26]=as.character(fullpeakinfo[,26]) #ensure sequence column is stored as character class
ourseqs=toupper(fullpeakinfo[,26]) #revert all sequences to uppercase

#now search for non-degenerate motifs at ALL peak sequences using the same algorithm, but do not update motif PWMs
print(date())
znewALL=getmotifs(mot,sizes,ourseqs[1:length(ourseqs)],maxwidth=301,alpha=1,incprob=0.99999,maxits=2,plen=0.9,updatemot=0,updatealpha=1,updateprior=0,ourprior=prior,bg=-1,plotting=FALSE)
for(iter in 1:100){
	print(paste("Iteration Number",iter))
	znewALL=getmotifs(mot,sizes,ourseqs[1:length(ourseqs)],maxwidth=301,alpha=znewALL$alpha,incprob=0.99999,maxits=2,plen=0.9,updatemot=0,updatealpha=1,updateprior=0,ourprior=prior,bg=-1,plotting=FALSE)
	save(znewALL,file="Human_Motif_Results_Final_ALLpeaks.r")
	print(date())
}

#add in motif info and print out final results
fullpeakinfo = cbind(fullpeakinfo,znewALL$whichmot, znewALL$whichpos, znewALL$whichstrand)
write.table(fullpeakinfo,"HumanPeaks.withMotifs.tbl",sep="\t",quote=F,row.names=F,col.names=F)

