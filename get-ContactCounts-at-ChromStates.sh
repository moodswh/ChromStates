#!/bin/bash

#method='CaptureC'
#days=(SC D3 D6)

method='HiC'
days=(RE100 RE40 RE20 RE10 RE5)

for (( i=0; i<${#days[@]}; i++ )); do
        day=${days[$i]}

ChromatinStatesFile=/srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/ChromatinStates/ChromStatesData/primarykeratinocyte_SC_14_states_dense.bed

#HiC:
resolution=$day
SignificantContactsFile=/srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/hic_flexible_repl_combined_SAVED/data/fit-hi-c-results/$resolution/afterICE/SC/SC.spline_pass2.pvals.txt.gz
RESitesFile=/srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/hic_flexible_repl_combined_SAVED/data/RESites/REFrags.$resolution

#CaptureC:
#SignificantContactsFile=/srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/promoters_flexible_repl_combined_OLD/data/fit-hi-c-results_50bins/RE1/afterICE/$day/$day.spline_pass2.pvals.txt.gz
#ChromatinStatesFile=/srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/ChromatinStates/ChromStatesData/primarykeratinocyte_${day}_14_states_dense.bed
#RESitesFile=/srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/promoters_flexible_repl_combined/data/RESites/REFrags.RE1

PROJDIR=/srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/ChromatinStates
python=/srv/gs1/software/python/python-2.7/bin/python
bedtools=/srv/gs1/software/bedtools/2.21.0/bin/bedtools

############################
## Get contacts with Q<=0.01 and Combine SignificantContactsFile and RESitesFile to get corresponding stat and end points of RE-sites
###########################

# write three files:
#file1: chr1 mid1 id 
#file2: chr2 mid2 id
#file3: id, NumOfContacts
zcat $SignificantContactsFile | awk '{if ($7<=0.01) print $1, $2, NR}' | sort -k1,1 -k2n,2 > $PROJDIR/file1
zcat $SignificantContactsFile | awk '{if ($7<=0.01) print $3, $4, NR}' | sort -k1,1 -k2n,2 > $PROJDIR/file2
zcat $SignificantContactsFile | awk '{if ($7<=0.01) print NR, $5}' | sort -k1n,1 > $PROJDIR/file3

cat $RESitesFile | sort -k1,1 > $PROJDIR/RESites

# write site1: chr1 start1 end1 id mid1 
#join based on two fields; id=bin index
join $PROJDIR/file1 $PROJDIR/RESites | awk '{ if ($2==$6) print $1, $4, $5, $3, $2}' | /usr/bin/perl -p -i -e 's/ /\t/g' > $PROJDIR/site1

# write site2: chr2 start2 end2 id mid2 
#join based on two fields; id=bin index
join $PROJDIR/file2 $PROJDIR/RESites | awk '{ if ($2==$6) print $1, $4, $5, $3, $2}' | /usr/bin/perl -p -i -e 's/ /\t/g' > $PROJDIR/site2

########################################
#Calculate NumOfContacts in each State-State bin
#######################################

# find intersection of RE-sites and chromatin states
# test1: id state1 LengthOfOverlap1
$bedtools intersect -wao -a $PROJDIR/site1 -b $ChromatinStatesFile >  $PROJDIR/overlap1
cat $PROJDIR/overlap1 | awk '{print $4, $9, $15}' | sort -k1n,1 -k2n,2 > $PROJDIR/test1
# test2: id state2 LengthOfOverlap2
$bedtools intersect -wao -a $PROJDIR/site2 -b $ChromatinStatesFile >  $PROJDIR/overlap2 
cat $PROJDIR/overlap2 | awk '{print $4, $9, $15}' | sort -k1n,1 -k2n,2 > $PROJDIR/test2

#Compute background probabilities
cat $PROJDIR/overlap1 > $PROJDIR/background_pairs
cat $PROJDIR/overlap2 >> $PROJDIR/background_pairs

cat $PROJDIR/background_pairs | awk '{print $1, $5, $9, $15}'| sort -k1,1 -k2n,2 -k3n,3 | awk '{print $4, $1, $2, $3}' | uniq -f 1 | awk '{print $4, $1}'| sort -k1n,1 > $PROJDIR/background_states_sorted
echo "" >> $PROJDIR/background_states_sorted
cat $PROJDIR/background_states_sorted | gawk '($1==key1) { sum+=$2} ($1 != key1) {if (NR>1){print key1, sum} key1=$1; sum=$2}'| awk '{a[NR]=$0;x+=(b[NR]=$2)}END{while(++i<=NR)print i"\t"b[i]/x}' > $PROJDIR/background_final
join -12 -22 -o 1.1,2.1 -t '
' $PROJDIR/background_final $PROJDIR/background_final | sed 'N; s/\n/ /' | awk '{print $1, $3, $2*$4}'| sort -k1n,1 -k2n,2 > $PROJDIR/background_prob

#need to sum raws
echo "" >> $PROJDIR/test1
echo "" >> $PROJDIR/test2

#write test11: id StateOvelap1(combined for a given state)
cat $PROJDIR/test1 | gawk '($1==key1 && $2==key2) { sum+=$3} ($1 != key1 || $2 != key2) {if (NR>1){print key1, key2, sum} key1=$1; key2=$2; sum=$3}' > $PROJDIR/test11
#write sum1: id TotalOvelap1(all states combined for a given id)
cat $PROJDIR/test11 | gawk '($1==key1) { sum+=$3} ($1 != key1) {if (NR>1){print key1, sum} key1=$1; sum=$3}'> $PROJDIR/sum1

#write test22: id StateOvelap2(combined for a given state)
cat $PROJDIR/test2 | gawk '($1==key1 && $2==key2) { sum+=$3} ($1 != key1 || $2 != key2) {if (NR>1){print key1, key2, sum} key1=$1; key2=$2; sum=$3}' > $PROJDIR/test22
#write sum2: id TotalOvelap2(all states combined for a given id)
cat $PROJDIR/test22 | gawk '($1==key1) { sum+=$3} ($1 != key1) {if (NR>1){print key1, sum} key1=$1; sum=$3}'> $PROJDIR/sum2

#write: id state OverlapFraction(for a given state in a bin)
join $PROJDIR/test11 $PROJDIR/sum1 | awk '{print $1, $2, $3/$4}'| sort -k1n,1 -k2n,2 > $PROJDIR/test11_sum1
join $PROJDIR/test22 $PROJDIR/sum2 | awk '{print $1, $2, $3/$4}'| sort -k1n,1 -k2n,2 > $PROJDIR/test22_sum2

#join two files by id: id state1 OverlapFraction1  state2 OverlapFraction2
join $PROJDIR/test11_sum1 $PROJDIR/test22_sum2 > $PROJDIR/pairs

#join pairs with ContactCounts by id and get ContactCounts for state pairs
join $PROJDIR/pairs $PROJDIR/file3| awk '{print $2, $4, $3*$5*$6}' | sort -k1n,1 -k2n,2 > $PROJDIR/pairs2
echo "" >> $PROJDIR/pairs2
#some ContactCounts for all uniq state pairs: state1 state2 ContactCounts
cat $PROJDIR/pairs2| gawk '($1==key1 && $2==key2) { sum+=$3} ($1 != key1 || $2 != key2) {if (NR>1){print key1, key2, sum} key1=$1; key2=$2; sum=$3}'  > $PROJDIR/${method}_${day}_StatePairs

#Normalize StatePairs
#join based on two fields
join $PROJDIR/${method}_${day}_StatePairs $PROJDIR/background_prob | awk '{ if ($2==$5) print $1, $2, $4, $5}' | /usr/bin/perl -p -i -e 's/ /\t/g' > $PROJDIR/${method}_${day}_StatePairs_norm
join -j 2 $PROJDIR/${method}_${day}_StatePairs $PROJDIR/background_prob |awk '{print $1, $2, $3/$5}' | sort -k1n,1 -k2n,2> $PROJDIR/${method}_${day}_StatePairs_norm

$python $PROJDIR/plot14by14.py $PROJDIR/${method}_${day}_StatePairs_norm 14 $method $day

rm $PROJDIR/file* $PROJDIR/test* $PROJDIR/pairs* $PROJDIR/site* $PROJDIR/RESites $PROJDIR/sum* $PROJDIR/overlap* $PROJDIR/background* $PROJDIR/${method}_${day}_StatePairs
done

#to run
#    PROJDIR=`pwd`
#    qsub -l h_vmem=20G -l h_rt=20:00:00 -m ea  -M taranova.maryna@gmail.com -o $PROJDIR/o.out -e $PROJDIR/e.error $PROJDIR/get-ContactCounts-at-ChromStates.sh 
