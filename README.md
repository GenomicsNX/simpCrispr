simpCrispr
=====================
**simpCrispr** is a module to analyze CRISPR offtarget sites.

It can find gRNAs for a given fasta file or given chr,start,end and for each gRNA it will find all potential offtarget sites and score these with the same formula used in (http://crispr.mit.edu/about), the score matrix can be modified if needed. Then it will rank gRNAs based on the total score of corresponding offtarget sites, will write primer region for a given number of top rated offtargets for all gRNAs, and finally it will write all the results to an excel file

It should be run on a unix machine and it mainly has bash and python scripts and it employs NGS tools, bowtie for alignment, bedtools and homer tools.
For python , pandas and numpy packages should be installed. 

##How it works:

These applications should be installed:

* It runs bowtie for alignment 
* It runs bedtools 
* It runs homer tools, annotatePeaks.pl if annotation is needed


```

example1 - input fasta given -

simpCrispr.sh -fasta examples/dbt.fa -outdir examples -out dbt-out
```
 
  
```

example2 - chr,start,end given -

bash simpCrispr.sh -chr X -start 101370000 -end 101370050 -out btk

```

 
