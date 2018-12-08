library("seqinr")
#set path
setwd('path/to/fasta')
#set fasta file with SEQ
filename='fasta.fa'
species <- read.fasta(file = filename)
speciesseq <- species[[1]]

#FUNCTION
slidingwindowplot <- function(windowsize, inputseq, speciesname)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  
  n <- length(starts)
  chunkGCs <- numeric(n) #empty vecor of length equal to the intervals
  
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    chunkGCs[i] <- chunkGC
  }
  return(chunkGCs)
}

#Set the number of bases you want each time to count the GC content
t=slidingwindowplot(100, speciesseq, filename)
write.table(t,paste(getwd(), '/', 'fasta_bpwindow.txt', sep=''), sep="\n", row.names=FALSE, col.names=filename)

