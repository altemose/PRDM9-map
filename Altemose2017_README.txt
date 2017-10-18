##############README for Altemose et al. eLife 2017 source code
#The following bash commands demonstrate our peak-calling and motif-finding pipelines

####################
#####1_Preprocessing
####################
##steps used to process ChIP-seq BAM files into fragment position bed files used
##for peak calling. example: YFP_Human peaks

##input files (downloadable from GEO accession GSE99407)
#	ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.hg19.bam
#	ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.hg19.bam
#	ChIPseq_Reads.YFP_HumanPRDM9.Input.ProtocolN.hg19.bam
#	ChIPseq_Reads.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.hg19.bam

##use samtools to remove likely PCR duplicate fragments
samtools rmdup ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.hg19.bam ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.hg19.rmdup.bam
samtools rmdup ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.hg19.bam ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.hg19.rmdup.bam
samtools rmdup ChIPseq_Reads.YFP_HumanPRDM9.Input.ProtocolN.hg19.bam ChIPseq_Reads.YFP_HumanPRDM9.Input.ProtocolN.hg19.rmdup.bam
samtools rmdup ChIPseq_Reads.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.hg19.bam ChIPseq_Reads.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.hg19.rmdup.bam
samtools index ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.hg19.rmdup.bam
samtools index ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.hg19.rmdup.bam
samtools index ChIPseq_Reads.YFP_HumanPRDM9.Input.ProtocolN.hg19.rmdup.bam
samtools index ChIPseq_Reads.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.hg19.rmdup.bam

##use GetFragmentDepth.pl to produce fragment position files and fragment depth files from BAM files
#USAGE: samtools view -F12 -q1 <SortedBAMfile.bam> | perl GetFragmentDepth.pl <Fragment Depth Output File> <Fragment Position Output File> <read length [default 51]> <max frag length [default 10000]> <path to bedtools executable [default ~/bin]> <chr size file [default hg19.chromsizes.tbl]>
samtools view -F12 -q1 ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.hg19.rmdup.bam | perl GetFragmentDepth.pl ChIPseq_FragmentDepth.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.q1.rmdup.hg19.bedgraph FragPos.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.bed
samtools view -F12 -q1 ChIPseq_Reads.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.hg19.rmdup.bam | perl GetFragmentDepth.pl ChIPseq_FragmentDepth.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.q1.rmdup.hg19.bedgraph FragPos.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.bed
samtools view -F12 -q1 ChIPseq_Reads.YFP_HumanPRDM9.Input.ProtocolN.hg19.rmdup.bam | perl GetFragmentDepth.pl ChIPseq_FragmentDepth.YFP_HumanPRDM9.Input.ProtocolN.q1.rmdup.hg19.bedgraph FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed
samtools view -F12 -q1 ChIPseq_Reads.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.hg19.rmdup.bam | perl GetFragmentDepth.pl ChIPseq_FragmentDepth.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.q1.rmdup.hg19.bedgraph FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed

##optional: use MakePseudoreplicates.pl to create pseudoreplicates when needed (e.g. when only 1 replicate is available, to make it compatible with the peak caller's expected input of 2 replicates)
#USAGE: perl MakePseudoreplicates.pl <input file> <Proportion in Replicate 1 [default 0.5]> <chr size file [default hg19.chromsizes.tbl]> <bedtools path [default ~/bin]>
#sample usage:
perl MakePseudoreplicates.pl FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed.sorted.bed


##################
#####2_PeakCalling
##################
##steps used to call peaks given fragment position bed files from two ChIP replicates and one
##genomic input chromatin control

##first, estimate constants alpha1, alpha2, and beta using EstimateConstants.R
#USAGE: Rscript EstimateConstants.R <path to directory containing input files> <sample name, used to generate output files> <filename of bedfile containing ChIP replicate 1 fragment positions> <filename of bedfile containing ChIP replicate 2 fragment positions> <filename of bedfile containing genomic input control fragment positions>
Rscript EstimateConstants.R ../1_Preprocessing YFP_HumanPRDM9.antiGFP FragPos.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.bed.sorted.bed FragPos.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.bed.sorted.bed FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed.sorted.bed
	#generates a file called "Constants.YFP_HumanPRDM9.antiGFP.txt"
Rscript EstimateConstants.R ../1_Preprocessing YFP_HumanPRDM9.antiH3K4me3 FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed.sorted.bed.PR1.sorted.bed FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed.sorted.bed.PR2.sorted.bed FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed.sorted.bed
	#generates a file called "Constants.YFP_HumanPRDM9.antiH3K4me3.txt"

##now, call peaks using these constant estimates
#USAGE: Rscript DeNovoPeakCalling-SingleBase.R <path to directory containing input files> <sample name, used to generate output files> <filename of bedfile containing ChIP replicate 1 fragment positions> <filename of bedfile containing ChIP replicate 2 fragment positions> <filename of bedfile containing genomic input control fragment positions> <filename of Constants file output by EstimateConstants.R> <p-value threshold> <minimum separation between peak centres> <binary flags indicating whether intermediate coverage files should be computed for each sample, default "1,1,1"> <binary flag indicating whether coverage files should be converted to binary format, default "1"> <binary flag indicating whether likelihood values should be computed, default "1"> 
Rscript DeNovoPeakCalling-SingleBase.R ../1_Preprocessing YFP_HumanPRDM9.antiGFP FragPos.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.bed.sorted.bed FragPos.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.bed.sorted.bed FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed.sorted.bed Constants.YFP_HumanPRDM9.antiGFP.txt 0.000001 250 1,1,1 1 1
	#generates a file called "SingleBasePeaks.YFP_HumanPRDM9.antiGFP.p0.000001.sep250.ALL.bed"
Rscript DeNovoPeakCalling-SingleBase.R ../1_Preprocessing YFP_HumanPRDM9.antiGFP FragPos.YFP_HumanPRDM9.antiGFP.Rep1.ProtocolN.bed.sorted.bed FragPos.YFP_HumanPRDM9.antiGFP.Rep2.ProtocolN.bed.sorted.bed FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed.sorted.bed Constants.YFP_HumanPRDM9.antiGFP.txt 0.000001 1000 0,0,0 0 0
	#generates a file called "SingleBasePeaks.YFP_HumanPRDM9.antiGFP.p0.000001.sep1000.ALL.bed"

##optional:use ForceCallPeaks.R to determine if there is significant enrichment around a fixed set of genomic positions
#USAGE: Rscript ForceCallPeaks.R <path to directory containing input files> <sample name, used to generate output files> <filename of bedfile containing ChIP replicate 1 fragment positions> <filename of bedfile containing ChIP replicate 2 fragment positions> <filename of bedfile containing genomic input control fragment positions> <filename of Constants file output by EstimateConstants.R> <bedfile containing positions of windows at which to force-call enrichment> <output filename> 
#sample usage (first generate a 3-column bed file with peak center positions only and no header):
awk -v OFS="\t" 'NR>1 {print $1,$2,$3}' SingleBasePeaks.YFP_HumanPRDM9.antiGFP.p0.000001.sep1000.ALL.bed >SingleBasePeaks.YFP_HumanPRDM9.antiGFP.p0.000001.sep1000.ALL.PeakCenters.bed
Rscript ForceCallPeaks.R ../1_Preprocessing FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed.sorted.bed.PR1.sorted.bed FragPos.YFP_HumanPRDM9.antiH3K4me3.ProtocolN.bed.sorted.bed.PR2.sorted.bed FragPos.YFP_HumanPRDM9.Input.ProtocolN.bed.sorted.bed Constants.YFP_HumanPRDM9.antiH3K4me3.txt SingleBasePeaks.YFP_HumanPRDM9.antiGFP.p0.000001.sep1000.ALL.PeakCenters.bed ChIPseq_ForceCalling.HumanH3K4me3_on_YFP_HumanPRDM9_peaks.protocolN.p10e-5.sep1000.txt
	#generates a file called "ChIPseq_ForceCalling.HumanH3K4me3_on_YFP_HumanPRDM9_peaks.protocolN.p10e-5.sep1000.txt"


###################
#####3_MotifFinding
###################

#see FindMotifs_Human.R for code used to call motifs
#with all required libraries installed, code can be executed simply as
Rscript FindMotifs_Human.R
Rscript FindMotifs_Chimp.R
#note: algorithm is stochastic, so exact results will vary run to run, though similar motifs will appear stably across runs


###############
#####4_Plotting
###############

#see Altemose2017_PlotCode.R for detailed plotting code
#data for profile plots (as in Figure 1c) were generated using MakeProfilePlot.pl
#USAGE: perl MakeProfilePlot.pl <Positions.bed> <Features.bedgraph> <Outfile.txt> <range [10000 default]> <value multiplier [default 1]> <decimal places to round to for printing means [default 4]> <motif length [default 1]> <path to bedtools [default ~/bin]>


##########
#####Other
##########

#Hapmap CEU recombination map data were downloaded from ftp://ftp.ncbi.nlm.nih.gov/hapmap/
#and converted to cM/Mb in the file HumanRecMap.HapMapCEU.CMperMb.bed
#using the perl script ParseHapMap.pl



