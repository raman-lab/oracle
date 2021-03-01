  for j in 1 2 3 4 5 6 7 8 9 10 11 12
  {
    cp ${j}/out.extendedFrags.fastq  ${j}/combined.fastq
    sed -i '/@M07074/d' ${j}/combined.fastq
    sed -i '/@M04241/d' ${j}/combined.fastq
    sed -i '/@M02713/d' ${j}/combined.fastq
    sed -i '2~3d' ${j}/combined.fastq
    sed -i 'N;s/\n/ /' ${j}/combined.fastq 
  } 
