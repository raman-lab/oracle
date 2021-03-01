#ATTENTION - this script will not work if your sequencing has N's
#To remove sequences with N, use sed -i '/pattern_to_match/d' ./infile_name

#Defined function dna2aa - unsuprisingly, it converts dna sequences (3bp) into an amino acid
#--------------------------------------------------------------------------------------------------------
dna2aa <- function(dna) 
{
  if(dna == "AAA")
  {
    return("K")
  }
  if(dna == "AAC")
  {
    return("N")
  }
  if(dna == "AAG")
  {
    return("K")
  }
  if(dna == "AAT")
  {
    return("N")
  }
  if(dna == "ACA")
  {
    return("T")
  }
  if(dna == "ACC")
  {
    return("T")
  }
  if(dna == "ACG")
  {
    return("T")
  }
  if(dna == "ACT")
  {
    return("T")
  }
  if(dna == "AGA")
  {
    return("R")
  }
  if(dna == "AGC")
  {
    return("S")
  }
  if(dna == "AGG")
  {
    return("R")
  }
  if(dna == "AGT")
  {
    return("S")
  }
  if(dna == "ATA")
  {
    return("I")
  }
  if(dna == "ATC")
  {
    return("I")
  }
  if(dna == "ATG")
  {
    return("M")
  }
  if(dna == "ATT")
  {
    return("I")
  }
  if(dna == "CAA")
  {
    return("Q")
  }
  if(dna == "CAC")
  {
    return("H")
  }
  if(dna == "CAG")
  {
    return("Q")
  }
  if(dna == "CAT")
  {
    return("H")
  }
  if(dna == "CCA")
  {
    return("P")
  }
  if(dna == "CCC")
  {
    return("P")
  }
  if(dna == "CCG")
  {
    return("P")
  }
  if(dna == "CCT")
  {
    return("P")
  }
  if(dna == "CGA")
  {
    return("R")
  }
  if(dna == "CGC")
  {
    return("R")
  }
  if(dna == "CGG")
  {
    return("R")
  }
  if(dna == "CGT")
  {
    return("R")
  }
  if(dna == "CTA")
  {
    return("L")
  }
  if(dna == "CTC")
  {
    return("L")
  }
  if(dna == "CTG")
  {
    return("L")
  }
  if(dna == "CTT")
  {
    return("L")
  }
  if(dna == "GAA")
  {
    return("E")
  }
  if(dna == "GAC")
  {
    return("D")
  }
  if(dna == "GAG")
  {
    return("E")
  }
  if(dna == "GAT")
  {
    return("D")
  }
  if(dna == "GCA")
  {
    return("A")
  }
  if(dna == "GCC")
  {
    return("A")
  }
  if(dna == "GCG")
  {
    return("A")
  }
  if(dna == "GCT")
  {
    return("A")
  }
  if(dna == "GGA")
  {
    return("G")
  }
  if(dna == "GGC")
  {
    return("G")
  }
  if(dna == "GGG")
  {
    return("G")
  }
  if(dna == "GGT")
  {
    return("G")
  }
  if(dna == "GTA")
  {
    return("V")
  }
  if(dna == "GTC")
  {
    return("V")
  }
  if(dna == "GTG")
  {
    return("V")
  }
  if(dna == "GTT")
  {
    return("V")
  }
  if(dna == "TAA")
  {
    return("Z")
  }
  if(dna == "TAC")
  {
    return("Y")
  }
  if(dna == "TAG")
  {
    return("Z")
  }
  if(dna == "TAT")
  {
    return("Y")
  }
  if(dna == "TCA")
  {
    return("S")
  }
  if(dna == "TCC")
  {
    return("S")
  }
  if(dna == "TCG")
  {
    return("S")
  }
  if(dna == "TCT")
  {
    return("S")
  }
  if(dna == "TGA")
  {
    return("Z")
  }
  if(dna == "TGC")
  {
    return("C")
  }
  if(dna == "TGG")
  {
    return("W")
  }
  if(dna == "TGT")
  {
    return("C")
  }
  if(dna == "TTA")
  {
    return("L")
  }
  if(dna == "TTC")
  {
    return("F")
  }
  if(dna == "TTG")
  {
    return("L")
  }
  if(dna == "TTT")
  {
    return("F")
  }  
}
#--------------------------------------------------------------------------------------------------------

#Defined function getBarcode - it scans the 5' end of a sequence and searches for a barcode
#--------------------------------------------------------------------------------------------------------
getBarcode <- function(seq) 
{
  bc <- 'none' #assumes no barcode, if it finds one it updates it
  count <- 0 #tracking how many times it finds a barcode
  
  for(nt in 1:20) #scan window from 1 to 21, takes 8 bp
  {
    bc_octet <- substr(seq, nt, (nt + 7)) #checks for this every nucleotide until pos21
    if(bc_octet == 'TAGCCAGC') 
    {
      count <- count + 1
      bc <- 'TAGCCAGC'
    } 
    if(bc_octet == 'CGGACAGC')
    {
      count <- count + 1
      bc <- 'CGGACAGC'
    }
    if(bc_octet == 'AGTACAGC')
    {
      bc <- 'AGTACAGC'
      count <- count + 1
    }
  }
  
  if(count <= 1) #if count of barcode is 1, returns barcode
  {
    return(bc)
  } else
  {
    return('unsure') #returns unsure if it has multiple, then pray
  }
}
#--------------------------------------------------------------------------------------------------------

#Defined function getReadingFrame - it scans the 5' end and searches for a user-defined start sequence
#--------------------------------------------------------------------------------------------------------
getReadingFrame <- function(seq) 
{
  start_position <- 0 #assumes cant find one at first, if it finds one below, updates updates to location
  count <- 0
  
  for(nt in 1:50)
  {
    read_octet <- substr(seq, nt, (nt + 7))
    if(read_octet == 'CGTGGACT') #checks for constant region, which is first coding region
    {
      count <- count + 1
      start_position <- nt + 8
    } 
  }
  
  if(count <= 1) #if it finds start octet, return pos of start translation, otherwise returns 0, or it found too many
  {
    return(start_position)
  } else
  {
    return(0)
  }
}
#--------------------------------------------------------------------------------------------------------

#Defined function getTranslation - it converts the dna sequence to aa sequence
#--------------------------------------------------------------------------------------------------------
getTranslation <- function(seq,start_pos) #take in seq and start pos
{
  aa.len <- floor((nchar(seq)-start_pos+1)/3) #determine AA it can translate
  aa.seq <- as.data.frame(matrix(0,1,aa.len)) # creates sequence matching max length of read
  for(aa.pos in 1:aa.len) # for every pos determine codon to translate
  {
    codon <- substr(seq,(start_pos+3*(aa.pos-1)),(start_pos+(3*aa.pos)-1))
    aa.seq[1,aa.pos] <- dna2aa(codon)
    if(aa.seq[1,aa.pos]=='Z')
    {
      break #stop translating if stop codon
    }
  }
  aa.seq[aa.seq==0] <- ""
  aa.seq <- apply(aa.seq, 1, paste, collapse = "")
  return(aa.seq) #finally returns AA seq
}
#--------------------------------------------------------------------------------------------------------

setwd("../1") #for phages, alter directory then remove _rev in file location


  all_seqs <- as.data.frame(read.csv(paste0("good-reads.csv"), header=T, stringsAsFactors=F))
  aa       <- as.data.frame(matrix(0,length(all_seqs[,1]),3)) #first column is barcode, then reading frame, then translation seq
  colnames(aa)[1:3] <- c('barcode', 'reading_frame_start','translation')

  for(seq_num in 1:length(all_seqs[,1]))
  {
    aa[seq_num,1] <- getBarcode(all_seqs[seq_num,1])
    aa[seq_num,2] <- getReadingFrame(all_seqs[seq_num,1])
    if(aa[seq_num,2] > 0)
    {
      aa[seq_num,3] <- getTranslation(all_seqs[seq_num,1],aa[seq_num,2])
    }
  }
  
  tsl.sum <- as.data.frame(matrix(0,4,3)) #this is translation summary matrix, how many times we see each barcode, % seq kept
  colnames(tsl.sum) <- c('barcode', '#_seqs', '%_total')
  c <- 0
  for(barcode in unique(aa[,1])) #for each barcode write out seq sequences in separate file
  {
    ss <- subset(aa[,3], aa[,1]==barcode & aa[,3]!=0)
    write.csv(ss, file=paste0("aa_", barcode, ".csv"), quote=F, row.names=F)
    
    c <- c + 1
    tsl.sum[c,1] <- barcode
    tsl.sum[c,2] <- length(ss)
    tsl.sum[c,3] <- 100 * length(ss) / length(all_seqs[,1])
  }
  write.csv(tsl.sum, file=paste0("tsl_summary.csv"), quote=F, row.names=F)

    
