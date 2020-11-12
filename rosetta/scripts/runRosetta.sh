#!/bin/bash

#Setup folders for each ddG calculation: variable positions are 462-553
#---------------------------------------------------------------------------------------------------------------------
mkdir ../muts/
for pos in {462..553}
{
  mkdir ../muts/pos_${pos}
  for aa in A C D E F G H I K L M N P Q R S T V W Y
  {
    mkdir ../muts/pos_${pos}/${aa}
  }
}
#---------------------------------------------------------------------------------------------------------------------

#Generate mutfiles that will inform Rosetta which mutations to make for each ddG calculations
#Variables pos1, pos2 and pos3 refer to individual monomers of the trimeric structure (each is mutated simultaneously)
#---------------------------------------------------------------------------------------------------------------------
fasta='FVAKSKAWTQVWSGSAGGGVSVTVSQDLRFRNIWIKCANNSWNFFRTGPDGIYFIASDGGWLRFQIHSNGLGFKNIADSRSVPNAIMVENEFVAKSKAWTQVWSGSAGGGVSVTVSQDLRFRNIWIKCANNSWNFFRTGPDGIYFIASDGGWLRFQIHSNGLGFKNIADSRSVPNAIMVENEFVAKSKAWTQVWSGSAGGGVSVTVSQDLRFRNIWIKCANNSWNFFRTGPDGIYFIASDGGWLRFQIHSNGLGFKNIADSRSVPNAIMVENE'
for pos in {1..91}
{
  pos1=$((${pos}+462))
  wt=$(echo ${fasta} | cut -c${pos}-${pos})
  for aa in A C D E F G H I K L M N P Q R S T V W Y
  {
    cp templates/temp.mutfile ../muts/pos_${pos1}/${aa}/pos_${pos1}_${aa}.mutfile
    sed -i "s|pos1|${pos}|g" ../muts/pos_${pos1}/${aa}/pos_${pos1}_${aa}.mutfile
    sed -i "s|mut1|${aa}|g" ../muts/pos_${pos1}/${aa}/pos_${pos1}_${aa}.mutfile
    sed -i "s|wt|${wt}|g" ../muts/pos_${pos1}/${aa}/pos_${pos1}_${aa}.mutfile
  }
}
for pos in {92..182}
{
  pos2=$((${pos}+371))
  wt=$(echo ${fasta} | cut -c${pos}-${pos})
  for aa in A C D E F G H I K L M N P Q R S T V W Y
  {
    sed -i "s|pos2|${pos}|g" ../muts/pos_${pos2}/${aa}/pos_${pos2}_${aa}.mutfile
    sed -i "s|mut2|${aa}|g" ../muts/pos_${pos2}/${aa}/pos_${pos2}_${aa}.mutfile
  }
}
for pos in {183..273}
{
  pos3=$((${pos}+280))
  wt=$(echo ${fasta} | cut -c${pos}-${pos})
  for j in A C D E F G H I K L M N P Q R S T V W Y
  {
    sed -i "s|pos3|${pos}|g" ../muts/pos_${pos3}/${aa}/pos_${pos3}_${aa}.mutfile
    sed -i "s|mut3|${aa}|g" ../muts/pos_${pos3}/${aa}/pos_${pos3}_${aa}.mutfile
  }
}
#---------------------------------------------------------------------------------------------------------------------

#Make and execute Condor submission files for each ddG calculation
#---------------------------------------------------------------------------------------------------------------------
mkdir ../subs/ ../shs/
for pos in {463..553}
{
  for aa in A C D E F G H I K L M N P Q R S T V W Y
  {
    cp templates/temp.sub ../subs/pos_${pos}_${aa}.sub
    sed -i "s|number|${pos}|g" ../subs/pos_${pos}_${aa}.sub
    sed -i "s|identity|${aa}|g" ../subs/pos_${pos}_${aa}.sub
    cp templates/temp.sh ../shs/pos_${pos}_${aa}.sh
    sed -i "s|number|${pos}|g" ../shs/pos_${pos}_${aa}.sh
    sed -i "s|identity|${aa}|g" ../shs/pos_${pos}_${aa}.sh
    condor_submit ../subs/pos_${pos}_${aa}.sub
  }
}
#---------------------------------------------------------------------------------------------------------------------


