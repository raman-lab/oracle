aa_list <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
ddg <- as.data.frame(matrix(0,1,22))
c <- 0
for(pos in 463:553)
{
  for(aa in aa_list)
  {
    c <- c + 1
    score <- read.table(paste0("../muts/pos_",pos,"/",aa,"/ddg_predictions.out"), header=T, stringsAsFactors=F)
    score[1,2] <- paste0(substr(score[1,2],1,1),pos,aa)
    ddg[c,] <- score[1,]
  }
}
colnames(ddg) <- colnames(score)
write.csv(ddg, "../analysis/ddG_breakdown.csv", row.names=F, quote=F)

