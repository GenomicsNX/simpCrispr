#!/usr/bin/env bash

####################################
#
# Written by Ilker Tunc
#
###################################



while [ "$1" != "" ]; do
  case $1 in
    -infile) shift; infile=$1;;    
    -guide) shift; guide=$1;;
    -top ) shift; top=$1;;
    -dn_primer) shift; dn_primer=$1;;
    -up_primer) shift; up_primer=$1;;
    -out ) shift; out=$1;;
	-ref ) shift; reference=$1;;
	-n ) shift; n=$1;;
  esac
shift
done



cc=0
while read line; 
do
	guide=$(echo $line | awk  '{ print $1}' )
	awk -v id=$guide ' $1==id ' $infile  | head -n $n > temp.offtarget 

	while read line2;
	do
		offtargetId=$(echo $line2 | awk  '{ print $3}' )
		offtarget=$(echo $line2 | awk  '{ print $7}' )
		offtargetpam=$(echo $line2 | awk  '{ print $8}' )
		chrom=$(echo $line2 | awk  '{ print $4}' )
		strand=$(echo $line2 | awk  '{ print $5}' )
		pos=$(echo $line2 | awk  '{ print $6}' )
 
		start=$(($pos-dn_primer))
		end=$(($pos+up_primer))

		#echo $chrom,$pos,$strand,$start,$end
		name=$guide,$offtargetId,$offtarget,${offtargetpam},chr${chrom},${start},${end},${strand},
		printf "%s\t%s\t%s\t%s\t%s\t%s\n" $chrom $start $end $name 0 $strand > bedf
		#printf "%s\t%s\t%s\t%s\t%s\t%s\n" $chrom $start $end $name 0 $strand 

		bedtools getfasta -s -fi $reference.fa  -bed bedf  -name -fo temp.primers -tab
		cat  temp.primers >> $out
		#((cc++))
		#if [[ $cc -ge $n ]]; then
			#break
		#fi
	done< temp.offtarget
done < $guide 


rm temp.primers temp.offtarget
rm bedf
