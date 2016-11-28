#!/bin/bash


usage(){
    echo "Usage: $0 "

    echo "-fasta      input fasta (list of gRNAs 20bp) or (region/s to search for gRNAS minimum 23bp)"
    echo "-pam        NGG, NGA or both , default both  "
    echo "-chr        chromosome name to find gRNAs [ignored if fasta file given]"
    echo "-start      start position to find gRNAs [ignored if fasta file given]"
    echo "-end        end position to find gRNAs [ignored if fasta file given]"
    echo "-strand     strand  to find gRNAs [ignored if fasta file given]"
    echo "-dn_primer  downstream bp length for primer region , default 300"
    echo "-up_primer  upstream bp length for primer region , default 300"
    echo "-top_primer number of top offtargets for primer region , default 10"
    echo "-annotate   annotate offtarget sites , default 1"
    echo "-mask       remove lower case letter bases , default 0"
    echo "-unmask     convert lower case to upper , default 0"
    echo "-excel      save to excel file, default 1"
    echo "-primers    find primer regions, default 1"
    echo "-offtargets find offtarget sites, otherwise just lists the gRNAs defult 1"
    echo "-out        outname prefix, default out"
    echo "-out_dir    outname directory, default current directory"
    echo "-h          help"
    exit 1
}


HOMEPATH="path where main simCrispr script and bin directort located"
BIN=$HOMEPATH/bin

REFERENCE="path to target genome fasta and bowtie indices"
GTF="path for target genome gtf file for annotation"

# bedtools, bowtie and homer annotatePeaks.pl should be installed and paths should exported 


while [ "$1" != "" ]; do
  case $1 in
    -fasta) shift; FASTAFILE=$1;;    
    -pam) shift; PAM=$1;;
    -chr ) shift; CHROM=$1;;
    -start ) shift; START=$1;;
    -end ) shift; END=$1;;
    -strand ) shift ; STRAND=$1;;
    -dn_primer) shift; DN_PRIMER=$1;;
    -up_primer) shift; UP_PRIMER=$1;;
    -top_primer) shift; N_TOP_PRIMER=$1;;
	-offtargets) shift;offtargets=$1;;
    -annotate) shift; annotate=$1;;
    -mask) shift; mask=$1;;
    -unmask) shift; unmask=$1;;
    -excel) shift; excel=$1;;
    -primers) shift; primers=$1;;
    -out ) shift; OUTNAME=$1;;
    -outdir ) shift ; SAVEDIR=$1;;
    -h ) usage;exit;;
    * ) usage;exit ;;
  esac
shift
done


# default parameters
CHROM=${CHROM:-1}
START=${START:-0}
END=${END:-23}
DN_PRIMER=${DN_PRIMER:-300}
UP_PRIMER=${UP_PRIMER:-300}
N_TOP_PRIMER=${TOP_PRIMER:-10}
STRAND=${STRAND:-+}

annotate=${annotate:-1}
mask=${mask:-0}
unmask=${unmask:-0}
excel=${excel:-1}
primers=${primers:-1}
offtargets=${offtargets:-1}

PAM=${PAM:-both}
OUTNAME=${OUTNAME:-out}
SAVEDIR=${SAVEDIR:-$(pwd)}


# temporary files generated
TEMPDIR=$(mktemp -d -p .)
TEMP=$TEMPDIR/tmp
OUT=$TEMPDIR/out


date


# if n fasta file is given then search the region given by chr,start,end
if [[ -z $FASTAFILE ]];then
    echo "no fasta file given,extracting sequence from specified position..."
    printf "chr%s\t%s\t%s\t%s\t%s\t%s\n" $CHROM $START $END  $OUT 0 > $TEMP.bed	
    bedtools getfasta -s -name  -fi $REFERENCE.fa  -bed $TEMP.bed -fo $TEMP.fa
    python $BIN/find_gRNAs.py -f $TEMP.fa -o $OUT.gRNA.fa -p $PAM
    rm $TEMP.*
else
    echo "Given fasta file, finding gRNAs..."
    python $BIN/find_gRNAs.py -f $FASTAFILE -o $OUT.gRNA.fa -p $PAM
fi


# if option -offtargtes 0, exit
if [[ $offtargets -eq 0 ]]; then 
	mv $OUT.gRNA.fa  $SAVEDIR/$OUTNAME.gRNA.fa
	rm -rf $TEMPDIR
	exit 1 
fi

echo "running bowtie..."
bowtie  -a -v 0 -p 1 $REFERENCE -f $OUT.gRNA.fa  $TEMP.bwt
bash $BIN/get_pam_list.sh $REFERENCE $TEMP.bwt > $TEMP.all.guides
python $BIN/check_pam_seq.py  $TEMP.all.guides $OUT 	
python $BIN/insert_mutations.py -i $OUT.guides  -o $OUT.mutations 

echo "running bowtie against $REFERENCE"
# report alignments with at most 3 MMs -v 3, no quality or seed length applies
bowtie  -a -v 3 -p 1 $REFERENCE -f $OUT.mutations $TEMP.bwt
python $BIN/bowtie2bed.py -i $TEMP.bwt -o $TEMP.bed
bedtools getfasta -s -fi  $REFERENCE.fa  -bed $TEMP.bed  -fo $TEMP.flanking -name -tab
# fix one missing column on bowtie format
awk ' {if (NF==7) print $0,"on-target";else print $0}' $TEMP.bwt > $TEMP.bwt.fixed
paste $TEMP.bwt.fixed  $TEMP.flanking > $OUT.merge


echo "filter out hits with no PAM matching..."
python $BIN/filter_hits.py -i $OUT.merge -g $OUT.guides -o $TEMP.filtered 
awk -v FS="\t" -v OFS="\t" '{print $2,$3,$8,$9,$10,$4,$5,$6,$7}' $TEMP.filtered > $TEMP.filtered.new
awk -v FS="\t" ' $9<5 {print $0}' $TEMP.filtered.new  | sort | uniq > $TEMP.nodups.filter
awk -v OFS="\t" ' {if ($9==0) label="on-target";else label="offtarget_"NR;print $1,$2,label,$3,$4,$5,$6,$7,$8,$9}' $TEMP.nodups.filter > $OUT.filter


# if option -mask exists
if [[ $mask -eq 1 ]]; then
	echo "removing lowercase bases..."	
	python $BIN/mask_OT_sites.py $OUT.filter $TEMP.masked.filter
	mv $TEMP.masked.filter $OUT.filter
fi


# if option -unmask exists
if [[ $unmask -eq 1 ]];then
	echo "converting lowercase bases to uppercase..."	
	python $BIN/unmask_OT_sites.py $OUT.filter $TEMP.unmasked.filter
	mv $TEMP.unmasked.filter $OUT.filter
fi


echo "calculating scores..."
python $BIN/score_offtargets.py -i $OUT.filter -o $TEMP.offtargets 
python $BIN/append_score.py -i $TEMP.offtargets -g $OUT.guides -o $TEMP.scored.guides
sort -k1,1 -k11,11nr $TEMP.offtargets > $OUT.sorted.offtargets
sort $TEMP.scored.guides > $OUT.sorted.guides


# if option -annotate exists
if [[ $annotate -eq 1 ]];then
	echo "annotating off-target sites for exons..."
	awk -v OFS="\t" ' {print $4,$6,$6+length($2),$1"_"$3,0,$5}'  $OUT.sorted.offtargets > $TEMP.bed
	annotatePeaks.pl $TEMP.bed $REFERENCE.fa -gtf $GTF  > $TEMP.annotate
	awk -v FS="\t" '{ print $1} ' $TEMP.annotate > $TEMP.col1
	awk -v FS="\t" '{ print $8} ' $TEMP.annotate | awk ' {print $1,$2} ' > $TEMP.col2
	paste $TEMP.col1 $TEMP.col2 | awk ' $2 ~/exon/ || $2 ~/intron/  { print $1,$2,$3} ' | sed 's/[(,]//g' > $TEMP.exon
	awk ' FNR==NR {aa[$1]=$2"-"$3;next} {if (aa[$1"_"$3]>0) print $0,aa[$1"_"$3];else print $0,"--"}' $TEMP.exon $OUT.sorted.offtargets > $OUT.offtargets.exons
fi


# if option -primers exists
if [[ $primers -eq 1 ]];then
	echo "finding region for primers..."
	bash  $BIN/write_primers.sh  -infile $OUT.sorted.offtargets -ref $REFERENCE -guide $OUT.sorted.guides  -dn_primer ${DN_PRIMER} -up_primer ${UP_PRIMER} -n ${N_TOP_PRIMER} -out $OUT.top${N_TOP_PRIMER}.primers
fi


# if option -excel exists
if [[ $excel -eq 1 ]];then
	echo "saving results to excel file..."
	python $BIN/write_to_excel.py $OUT $OUT ${N_TOP_PRIMER} ${DN_PRIMER} ${UP_PRIMER}
fi

   
echo "removing temporary files..."

mkdir -p $SAVEDIR
mv $OUT.sorted.guides  $SAVEDIR/$OUTNAME.guides
mv $OUT.sorted.offtargets  $SAVEDIR/$OUTNAME.offtargets
mv $OUT.xlsx  $SAVEDIR/$OUTNAME.xlsx
mv $OUT.top${N_TOP_PRIMER}.primers  $SAVEDIR/$OUTNAME.top${N_TOP_PRIMER}.primers

date

rm -rf $TEMPDIR

