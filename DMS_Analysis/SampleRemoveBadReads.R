  
#set folder
folders <- "../1"
  mat.sum <- as.data.frame(matrix(0,length(folders),3)) #creating summary matrix, length of samples we have, then 3 columns
  colnames(mat.sum) <- c('sample','#_merged_sequences','#_good_sequences') #will name columns

  c <- 0 #variable, initialize to 0
  for(i in folders) #i is another variable
  {
    raw.seqs <- read.table(paste(i, "/combined.fastq", sep=""), comment.char="", stringsAsFactors=F, header=F) #reads in combined.fastq
    add.cols <- as.data.frame(matrix(0,length(raw.seqs[,1]),2)) #adds two more columns to keep track of expected error and seq length
    colnames(add.cols) <- c('total_expected_error','sequence_length')
    c <- c + 1 #every time it changes a new folder it will increment c
    mat.sum[c,1] <- i #in pos matrix set
    mat.sum[c,2] <- length(raw.seqs[,1]) #for column 2 put raw seq

    for(j in 1:length(raw.seqs[,1])) #j is variable, go through fastq until end
    {
      acs <- strtoi(charToRaw(as.character(raw.seqs[j,2])),16L) #stringtointeger, convert to integer
      q <- acs - 33 #minus 33 for phred score, now its a q score
      # ss.q <- subset(q, q>30) #unsure what this is, #'ed out
      q.adj <- (-1 * q)/10 #convert q score to probability, 10 log10
      p <- 10^q.adj 
      add.cols[j,1] <- sum(p) #adding probability, total expected error, tracking this for each seq. 
      add.cols[j,2] <- nchar(as.character(raw.seqs[j,2])) #determine sequence length here
    }
    merge  <- cbind(raw.seqs,add.cols) #binding fastq with total expected error and seq length
    ss.err <- subset(merge, merge[,3] < 1) #this is the filter, if error is less than 1 keep seq.
    mat.sum[c,3] <- length(ss.err[,1]) #count seq that pass error filter
    write.csv(ss.err[,1], file=paste(i, "/good-reads.csv", sep=""), row.names=F, quote=F) #write it
    pdf(file=paste(i, "/seq_lengths.pdf", sep="")) #done be a pdf histogram that is invisible, cannot make png on cluster, only pdf
    par(mar=c(5,5,4,2))
    hist(ss.err[,4], xlab="Sequence Length (bp)", cex.axis=1.5, 
         ylab="Frequency", main="", cex.lab=1.5)
    dev.off()
    pdf(file=paste(i, "/expected_errors.pdf", sep=""))
    par(mar=c(5,5,4,2))
    hist(merge[,3], xlab="Total Expected Error", cex.axis=1.5,
         ylab="Frequency", main="", cex.lab=1.5)
    dev.off()
    pdf(file=paste(i, "/filtered_expected_errors.pdf", sep=""))
    par(mar=c(5,5,4,2))
    hist(ss.err[,3], xlab="Total Expected Error", cex.axis=1.5,
         ylab="Frequency", main="", cex.lab=1.5)
    dev.off()
  }
  
  write.csv(mat.sum, file="summary.csv", row.names=F, quote=F)

