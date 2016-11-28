#!/usr/bin/env bash

####################################
#
# Written by Ilker Tunc
#
###################################


#REFERENCE=/home/tunci/crispr_ZouJ/Reference/mmul/Macaca_mulatta.MMUL_1.softmasked.genome
#bedtools=/v/apps/bedtools2/bedtools2.25/bin/bedtools 


reference=$1
inlist=$2


#while IFS=$'\t' read -r c1 c2 c3 c4 c5 c6 ; do
#    echo $c1,$c2,$c3
#done < $1


while read line; do
    ID=$(echo $line | awk '{print $1}')
    st=$(echo $line | awk '{print $2}')
    chrom=$(echo $line | awk '{print $3}')
    pos=$(echo $line | awk '{print $4}')
    seq0=$(echo $line | awk '{print $5}')
    L=${#seq0}

    if [[ $st == "+" ]];then
	#echo "plus"
	start=$pos 
	end=$(( $pos + $L )) 
    fi

    if [[ $st == "-" ]];then
	#echo "negative"
	start=$pos
	end=$(( $pos + $L)) 
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" $chrom $start $end "name" 0 $st > temp.bed	
    bedtools getfasta -s -fi $reference.fa  -bed temp.bed  -fo temp.out -tab
    seq=$(awk '{print $2}' temp.out) 


    if [[ $st == "+" ]];then
	#echo "plus"
	start=$(( $pos + $L )) 
	end=$(( $start + 3 )) 
    fi

    if [[ $st == "-" ]];then
	#echo "negative"
	start=$(( $pos - 3 ))
	end=$(( $pos )) 
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" $chrom $start $end "name" 0 $st > temp.bed	
    bedtools getfasta -s -fi $reference.fa  -bed temp.bed  -fo temp.out -tab
    pam=$(awk '{print $2}' temp.out) 

    
    #printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $ID $seq $pam $chrom $st $pos $seq0
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $ID $seq $pam $chrom $st $pos
    


done <$inlist


rm temp.bed temp.out
