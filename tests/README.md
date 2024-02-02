Summary statistics were downloaded from the following links:

```bash
# Neuroticism
wget https://ctg.cncr.nl/documents/p1651/sumstats_neuroticism_ctg_format.txt.gz

# Height
wget https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
```

Following Speed's instructions:

```bash
# Neuroticism
gunzip -c sumstats_neuroticism_ctg_format.txt.gz | 
    awk '(NR>1){snp=$3":"$4;a1=$5;a2=$6;z=$9;n=$11}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T") && $12>0.95){print snp, a1, a2, z, n}' - |
    awk '!seen[$1]++' > neur.txt.old

# Height
gunzip -c GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz | awk '(NR>1){snp=$1;a1=$2;a2=$3;Z=($5/$6);n=$8}(NR==1){print "Predictor A1 A2 Z n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, Z, n}' - > height.raw

wget https://www.dropbox.com/s/xabjdu6squ6u56r/hapmap3.snps

awk '(NR==FNR){arr[$1]=$2;ars[$1]=$3$4;next}(FNR==1){print $0}($1 in arr && ($2$3==ars[$1]||$3$2==ars[$1])){$1=arr[$1];print $0}' hapmap3.snps height.raw > height.txt.old
```

Finally, a simple filter to minimize the number of lines in these files:

```bash
# Neuroticism
xsv join -d ' ' Predictor ldak.thin.genotyped.gbr.tagging Predictor neur.txt.old |
    xsv select 'Predictor,A1,A2,Z,n' |
    xsv fmt -t ' ' -o neur.txt

# Height
xsv join -d ' ' Predictor ldak.thin.genotyped.gbr.tagging Predictor height.txt.old |
    xsv select 'Predictor,A1,A2,Direction,Stat,n' |
    xsv fmt -t ' ' -o height.txt
```

The actual outputs of LDAK on these data are in `data/neur.hers` and `data/height.hers`.
