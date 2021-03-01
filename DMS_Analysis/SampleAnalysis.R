#This script counts the number of single mutants

#set sample)
for(sample in c('1'))
{
  setwd(paste0("../", sample))
  ref <- 'QVWSGSAGGGVSVTVSQDLRFRNIWIKCANNSWNFFRTGPDGIYFIASDGGWLRFQIHSNGLGFKNIADSRSVPNAIMVENEZ' #set wt
  ref_length <- nchar(ref)
  aa_csvs <- list.files(pattern='aa_')
  
  for(tsl_file in aa_csvs) #loop across all translation files (one for each barcode)
  {
    seqs <- read.csv(paste0(tsl_file), comment.char="", stringsAsFactors=F, header=T) #load csv of sequences
    for(seq_num in 1:length(seqs[,1])) #loop across all sequences in each file
    {
      c_muts <- 0 #count of number of mutations initialized to 0
      seqs[seq_num,2] <- nchar(seqs[seq_num,1]) #total sequence length
      for(aa in 1:seqs[seq_num,2]) #loop across all positions in each sequence and compare to WT
      {
        if(aa > ref_length) #if for some reason the translated seq length exceeds the reference, stop comparing aa's
        {
          break
        }
        if(substr(ref,aa,aa)!=substr(seqs[seq_num,1],aa,aa)) #if translated does not match WT, increment mutation count
        {
          c_muts <- c_muts + 1
        }
      }
      seqs[seq_num,3] <- c_muts #record the number of mutations
    }
    colnames(seqs) <- c('aa_sequence','length','#_muts')
    ss.5muts.or.less <- subset(seqs, seqs[,3]<=5) 

    pdf(file=paste0("analysis_pdfs/",sub('\\.csv$','',tsl_file),"_number_muts_histogram.pdf")) 
    par(mar=c(5,5,4,2))
    hist(seqs[,3],main="",cex.axis=1.5,xlab="# Mutations",
         cex.lab=1.5, cex.main=1.5)
    dev.off()
    
    pdf(file=paste0("analysis_pdfs/",sub('\\.csv$','',tsl_file),"_number_muts_histogram_5_or_less.pdf")) 
    par(mar=c(5,5,4,2))
    hist(ss.5muts.or.less[,3],main="",cex.axis=1.5,xlab="# Mutations",breaks=seq(-0.5,5.5,by=1),
         cex.lab=1.5, cex.main=1.5)
    dev.off()
    
    pdf(file=paste0("analysis_pdfs/",sub('\\.csv$','',tsl_file),"_sequence_lengths.pdf")) 
    par(mar=c(5,5,4,2))
    hist(seqs[,2],main="",cex.axis=1.5,xlab="Length (AA)",
         cex.lab=1.5, cex.main=1.5)
    dev.off()
    
    write.csv(seqs, file=paste0("analysis_csvs/",sub('\\.csv$','',tsl_file),"_tsl_summary.csv"), row.names=F)
    
    ss.1mut <- subset(seqs,seqs[,3]==1)
    ss.2mut <- subset(seqs,seqs[,3]==2) 
    
    dist.1mut <- matrix(0,nchar(ref),21)
    colnames(dist.1mut) <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','Z')
    dist.2mut <- dist.1mut

    if(length(ss.1mut[,1]) != 0)
    {
    for(i in 1:length(ss.1mut[,1]))
     {
      for(j in 1:nchar(ref))
      {
        aa.temp <- substr(ss.1mut[i,1], j, j)
        for(k in 1:length(dist.1mut[1,]))
        {
          if(aa.temp == colnames(dist.1mut)[k])
          {
            dist.1mut[j,k] <- dist.1mut[j,k] + 1
            break
          }
        }
      }
     }
    }
    if(length(ss.2mut[,1]) != 0)
    {
    for(i in 1:length(ss.2mut[,1]))
     {
      for(j in 1:nchar(ref))
      {
        aa.temp <- substr(ss.2mut[i,1], j, j)
        for(k in 1:length(dist.2mut[1,]))
        {
          if(aa.temp == colnames(dist.2mut)[k])
          {
            dist.2mut[j,k] <- dist.2mut[j,k] + 1
            break
          }
        }
      }
     }
    }


    write.csv(dist.1mut, file=paste0("analysis_csvs/",sub('\\.csv$','',tsl_file),"_aa_distribution_1mut.csv"), row.names=F)
    write.csv(dist.2mut, file=paste0("analysis_csvs/",sub('\\.csv$','',tsl_file),"_aa_distribution_2mut.csv"), row.names=F)
  }
}

  
    
  

    
