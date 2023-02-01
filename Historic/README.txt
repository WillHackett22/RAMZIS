Hello!
RAMZIS (Ranking Assessment of m/z Identifications by Similarity) is a tool to help researchers evaluate their data and perform better analysis.

The code has recently (6-10-19) undergone a major revision and a lot of its code became defunct leading me to purge it.

It does not currently function as a traditional R Package. Instead once you clone the repository the way to use it is as follows.

In R:
filelocation<-Location Of The Cloned Repository
source(paste0(filelocation,'RAMZISDraft.R'))

filename1<-First File of GP abundances to look at
filename2<-Second File of GP abundances to look at
plottitle<-What do you want the plot to be called

SimObject<-SimPlotFromFile(plottile,filename1,filename2)

From there the SimObject has a bunch of information in it.

Things that I'm still revamping the code for on my branch:
Functions:
Overall Data Quality Assessment
Individual Identification Data Quality Assessment
Ranking Assessment
	-Plotting
	-Output
	-Thresholding

Normalization Vector Management