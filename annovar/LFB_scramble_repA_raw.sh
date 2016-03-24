#!/bin/bash
set -e
seqtools_annovar=/home/brb/annovar
export bdge_bcftools_PATH=/opt/SeqTools/bin/bcftools-1.3
export PATH=$bdge_bcftools_PATH:$PATH

export bdge_htslib_PATH=/opt/SeqTools/bin/samtools-1.3/htslib-1.3
export PATH=$bdge_htslib_PATH:$PATH

inputVCF="/home/brb/SeqTestdata/RNASeqFibroblast/outputhg38/LFB_scramble_repA_raw.vcf"
cosmicVCF="/home/brb/SeqTestdata/usefulvcf/hg38/CosmicCodingMuts.vcf.gz"
dbSNPVCF="/home/brb/SeqTestdata/usefulvcf/hg38/common_all_20151104.vcf.gz"
genomeRef="/home/brb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
minQual=1
minReadDepth=1
minMapQual=1
genomeVer=hg38
RcodeDir="/opt/SeqTools/bin/SeqTools/code"
outputDir="/home/brb/SeqTestdata/RNASeqFibroblast/outputhg38_cli_annovar"
tmpname="$(basename $inputVCF)"
outputFileNameShare="${tmpname::-4}"
tmpfd="tmp_$outputFileNameShare"
if [ ! -d "$outputDir/tmp/$tmpfd" ]; then mkdir -p "$outputDir/tmp/$tmpfd";  fi

echo "running sample: $inputVCF"
bcftools filter -i"QUAL >= $minQual && DP >= $minReadDepth && MQ >= $minMapQual" "$inputVCF" > "$outputDir/tmp/$tmpfd/filtered.vcf"
(bcftools norm -m-both -o "$outputDir/tmp/$tmpfd/splitted.vcf" "$outputDir/tmp/$tmpfd/filtered.vcf") 2>&1
(bcftools norm -f "$genomeRef" -o "$outputDir/tmp/$tmpfd/leftnormalized.vcf" "$outputDir/tmp/$tmpfd/splitted.vcf") 2>&1
count=$(grep -v '#' "$inputVCF" | cut -f1 | grep 'chr' | wc -l)
if (( $count > 0 ))
then
	perl -pe 's/^chr//' "$outputDir/tmp/$tmpfd/leftnormalized.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.nochr.vcf"
	bgzip -c "$outputDir/tmp/$tmpfd/leftnormalized.nochr.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.nochr.vcf.gz"
	tabix -p vcf "$outputDir/tmp/$tmpfd/leftnormalized.nochr.vcf.gz"
else
	bgzip -c "$outputDir/tmp/$tmpfd/leftnormalized.vcf" > "$outputDir/tmp/$tmpfd/leftnormalized.vcf.gz"
	tabix -p vcf "$outputDir/tmp/$tmpfd/leftnormalized.vcf.gz"
fi
if [[ $cosmicVCF == *.vcf ]]; then bgzip -c "$cosmicVCF" > "$cosmicVCF.gz"; tabix -p vcf "$cosmicVCF.gz"; cosmicVCF="$cosmicVCF.gz"; fi
if [[ $dbSNPVCF == *.vcf ]]; then bgzip -c "$dbSNPVCF" > "$dbSNPVCF.gz"; tabix -p vcf "$dbSNPVCF.gz"; dbSNPVCF="$dbSNPVCF.gz"; fi
## NOTICE: .vcf.gz files will be unzipped, bgzipped and then tabix-ed since it might not be tabix-ed correctly.
if [[ $dbSNPVCF == *.vcf.gz ]]; then tabix -f -p vcf "$dbSNPVCF"; fi
if [[ $cosmicVCF == *.vcf.gz ]]
then
  if tabix -f -p vcf "$cosmicVCF"; then
    echo "The cosmic VCF file can be tabixed."
  else
    gunzip -f -d "$cosmicVCF" > "${cosmicVCF::-3}"
    bgzip -c "${cosmicVCF::-3}" > "$cosmicVCF";
    tabix -f -p vcf "$cosmicVCF"
  fi
fi
### dbsnp annotation via bcftools
if (( $count > 0 ))
then
	bcftools annotate -c ID -a "$dbSNPVCF" "$outputDir/tmp/$tmpfd/leftnormalized.nochr.vcf.gz" > "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf"
else
	bcftools annotate -c ID -a "$dbSNPVCF" "$outputDir/tmp/$tmpfd/leftnormalized.vcf.gz" > "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf"
fi
bgzip -c "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf" > "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf.gz"
tabix -p vcf "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf.gz"
### cosmic annotation via bcftools
bcftools annotate -c ID,+GENE -a "$cosmicVCF" "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf.gz" > "$outputDir/tmp/$tmpfd/cosmic_dbsnp.vcf"
sed '/\trs[0-9]\+\t/d' "$outputDir/tmp/$tmpfd/cosmic_dbsnp.vcf" > "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf"
#bgzip -c $outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf > $outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf.gz
#tabix -p vcf $outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf.gz
##download database files required for running ANNOVAR annotation
file1=$genomeVer
file1+="_refGene.txt"
file2=$genomeVer
file2+="_refGeneMrna.fa"
file3=$genomeVer
file3+="_refGeneVersion.txt"
if [ ! -f "$seqtools_annovar/humandb/$file1" ] || [ ! -f "$seqtools_annovar/humandb/$file2" ] || [ ! -f "$seqtools_annovar/humandb/$file3" ] ; then
    (perl $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar refGene "$seqtools_annovar/humandb/") 2>&1;
fi;
file1=$genomeVer
file1+="_knownGene.txt"
file2=$genomeVer
file2+="_knownGeneMrna.fa"
if [ ! -f "$seqtools_annovar/humandb/$file1" ] || [ ! -f "$seqtools_annovar/humandb/$file2" ] ; then
    (perl $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar knownGene "$seqtools_annovar/humandb/") 2>&1;
fi
if [ "$genomeVer" != "hg38" ] ; then
  file1=$genomeVer
  file1+="_ensGene.txt"
  file2=$genomeVer
  file2+="_ensGeneMrna.fa"
  if [ ! -f "$seqtools_annovar/humandb/$file1" ] || [ ! -f "$seqtools_annovar/humandb/$file2" ] ; then
    (pel $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar ensGene "$seqtools_annovar/humandb/") 2>&1 | tee -a  "$outputDir/tmp/log.txt";
  fi
fi
file1=$genomeVer
file1+="_dbnsfp30a.txt"
file2=$genomeVer
file2+="_dbnsfp30a.txt.idx"
if [ ! -f "$seqtools_annovar/humandb/$file1" ] || [ ! -f "$seqtools_annovar/humandb/$file2" ] ; then
    (perl $seqtools_annovar/annotate_variation.pl -buildver $genomeVer -downdb -webfrom annovar dbnsfp30a "$seqtools_annovar/humandb/") 2>&1;
fi
## convert vcf file to avinput, a format used by annovar
(perl $seqtools_annovar/convert2annovar.pl -format vcf4 "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" > "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.avinput" -includeinfo -comment) 2>&1
## identify splicing variants and nonsynonymous SNPs
(perl $seqtools_annovar/variants_reduction.pl "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.avinput" "$seqtools_annovar/humandb/" -protocol nonsyn_splicing,dominant -operation g,m -out "$outputDir/tmp/$tmpfd/reduced" -buildver $genomeVer) 2>&1
## generate a vcf file from the reduced results in .avinput
cut -f8- "$outputDir/tmp/$tmpfd/reduced.step2.varlist" > "$outputDir/tmp/$tmpfd/tmpfile"
grep '#' "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" > "$outputDir/tmp/$tmpfd/metainfo"
cat "$outputDir/tmp/$tmpfd/metainfo" "$outputDir/tmp/$tmpfd/tmpfile" > "$outputDir/tmp/$tmpfd/reduced.vcf"
rm "/$outputDir/tmp/$tmpfd/tmpfile"
rm "/$outputDir/tmp/$tmpfd/metainfo"
## annotate the reduced vcf through annovar
perl "$seqtools_annovar/convert2annovar.pl" -includeinfo -allsample -withfreq -format vcf4 "$outputDir/tmp/$tmpfd/reduced.vcf" > "$outputDir/tmp/$tmpfd/test.avinput"
(perl "$seqtools_annovar/table_annovar.pl" "$outputDir/tmp/$tmpfd/reduced.vcf" "$seqtools_annovar/humandb/" -buildver $genomeVer -out "$outputDir/tmp/$tmpfd/$outputFileNameShare" -remove -protocol refGene,dbnsfp30a -operation g,f -nastring . -vcfinput) 2>&1
annoTable="$outputFileNameShare"
annoTable+='_annoTable.txt'
annoVCF="$outputFileNameShare"
annoVCF+='_annotated.vcf'
file="$outputFileNameShare.$genomeVer"
file+='_multianno'
mv "$outputDir/tmp/$tmpfd/$file.txt" "$outputDir/tmp/$tmpfd/$annoTable"
mv "$outputDir/tmp/$tmpfd/$file.vcf" "$outputDir/tmp/$tmpfd/$annoVCF"
## re-arrange the annotation table
cut -f -7,9-44,50 "$outputDir/tmp/$tmpfd/$annoTable" > "$outputDir/tmp/$tmpfd/tmp0.txt"
grep -v 'Chr' "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp.txt"
grep 'Chr' "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmph.txt"
awk '{ printf("%s\tCOSMIC ID\n", $0); }' "$outputDir/tmp/$tmpfd/tmph.txt" > "$outputDir/tmp/$tmpfd/tmph1.txt"
cat "$outputDir/tmp/$tmpfd/tmph1.txt" "$outputDir/tmp/$tmpfd/tmp.txt" > "$outputDir/tmp/$tmpfd/tmp0.txt"
cut -f 1-5 "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp1.txt"
cut -f 6-43 "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp2.txt"
cut -f 44 "$outputDir/tmp/$tmpfd/tmp0.txt" > "$outputDir/tmp/$tmpfd/tmp3.txt"
paste -d'\t' "$outputDir/tmp/$tmpfd/tmp1.txt" <(cat "$outputDir/tmp/$tmpfd/tmp3.txt") > "$outputDir/tmp/$tmpfd/tmp4.txt"
paste -d'\t' "$outputDir/tmp/$tmpfd/tmp4.txt" <(cat "$outputDir/tmp/$tmpfd/tmp2.txt") >  "$outputDir/tmp/$tmpfd/$annoTable"
## reorder the genelist file and save it as genelist.txt
grep 'Number_of_deleterious_alleles' "$outputDir/tmp/$tmpfd/reduced.step2.genelist" > "$outputDir/tmp/$tmpfd/tmpa"
grep -v 'Number_of_deleterious_alleles' "$outputDir/tmp/$tmpfd/reduced.step2.genelist"  | sort > "$outputDir/tmp/$tmpfd/tmpb"
genelistFile="$outputFileNameShare"
genelistFile+='_genelist.txt'
cat "$outputDir/tmp/$tmpfd/tmpa" "$outputDir/tmp/$tmpfd/tmpb" > "$outputDir/tmp/$tmpfd/$genelistFile"
cut -f1,3 "$outputDir/tmp/$tmpfd/$genelistFile" > "$outputDir/tmp/$tmpfd/tmp_gl.txt"

rm "/$outputDir/tmp/$tmpfd/tmpa"
rm "/$outputDir/tmp/$tmpfd/tmpb"
## gene annotation through refseq, ensembl and ucsc, for the purpose of drawing figures in R
if [ "$genomeVer" == "hg38" ]
  then
	perl $seqtools_annovar/table_annovar.pl "$outputDir/tmp/$tmpfd/leftnormalized.vcf" "$seqtools_annovar/humandb/" -buildver $genomeVer -out "$outputDir/tmp/$tmpfd/anno_all" -remove -protocol refGene,knownGene -operation g,g -nastring . -vcfinput >> "$outputDir/tmp/log.txt"
	file=anno_all.$genomeVer
  else
	perl $seqtools_annovar/table_annovar.pl "$outputDir/tmp/$tmpfd/leftnormalized.vcf" "$seqtools_annovar/humandb/" -buildver $genomeVer -out "$outputDir/tmp/$tmpfd/anno_all" -remove -protocol refGene,knownGene,ensGene -operation g,g,g -nastring . -vcfinput >> "$outputDir/tmp/log.txt"
	file=anno_all.$genomeVer
fi	
file+='_multianno'
if [ "$genomeVer" == "hg38" ]
  then
	cut -d$'\t' -f1-15 "$outputDir/tmp/$tmpfd/$file.txt" > "$outputDir/tmp/$tmpfd/anno_all.txt"
  else
	cut -d$'\t' -f1-20 "$outputDir/tmp/$tmpfd/$file.txt" > "$outputDir/tmp/$tmpfd/anno_all.txt"
fi
## run R code to generate figures
($RcodeDir/./statsVariantGeneAnno_Auto.sh "$outputDir/tmp/$tmpfd/anno_all.txt" "$outputDir/tmp/$tmpfd") 2>&1
### Generate a log file
read num0 <<< $(grep -v '#' $inputVCF | wc -l)
echo -e "$inputVCF contains $num0 lines where each line corresponds to one or multiple variants." 
read num0a <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/filtered.vcf" | wc -l)
echo -e "There are $num0a variants after pre-filtering." 
read num <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/leftnormalized.vcf" | wc -l)
echo -e "There are $num variants after decomposing and normalizing variants."
read num2 <<< $(cut -f 3 "$outputDir/tmp/$tmpfd/dbsnp_anno.vcf" | grep -v '#' | grep 'rs' | wc -l)
echo -e "There are $num2 variants (out of $num variants) that are reported by dbSNP."
read num3 <<< $(cut -f 3 "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" | grep -v '#' | grep 'COSM' | wc -l)
echo -e "There are $num3 variants (out of $num variants) that are reported by COSMIC. "
read num3a <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" | wc -l)
echo -e "There are $((num3a-num+num2)) variants (out of $num variants) that are reported by both dbSNP build 144 and COSMIC v74."
echo -e "After removing these variants reported in dbSNP database and keeping these variants reported in COSMIC database, there are $num3a variants left for analysis." 
read num4 <<< $(wc -l "$outputDir/tmp/$tmpfd/reduced.step2.varlist" | awk '{print $1}')
read num5 <<< $(grep -v 'Number_of_deleterious_alleles' "$outputDir/tmp/$tmpfd/$genelistFile" | wc -l)
echo -e "There are $num4 variants (out of $num3a variants) that are nonsynonymous or splicing ones, which are kept for further analysis."
read num6 <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/reduced.vcf" | grep "COSM" | wc -l)
echo -e "There are $num6 variants (out of $num4 variants) that are reported by COSMIC v74."
read num_spg <<< $(grep -v 'exonic;splicing' "$outputDir/tmp/$tmpfd/$annoTable" | grep 'splicing' | wc -l)
echo -e "There are $num_spg variants (out of $num4 variants) that are splicing."
read num_del <<< $(grep 'frameshift deletion' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
echo -e "There are $num_del variants (out of $num4 variants) that are frameshift deletion."
read num_ins <<< $(grep 'frameshift insertion' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
echo -e "There are $num_ins variants (out of $num4 variants) that are frameshift insersion."
read num_sl <<< $(grep 'stoploss' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
echo -e "There are $num_sl variants (out of $num4 variants) that are stoploss."
read num_sg <<< $(grep 'stopgain' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
echo -e "There are $num_sg variants (out of $num4 variants) that are stopgain."
read num_non <<< $(grep 'nonsynonymous' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
echo -e "There are $num_non variants (out of $num4 variants) that are nonsynonymous."
echo -e "There are $num5 genes that are associated with $num4 variants, and the gene list is saved in "$outputDir/tmp/$tmpfd/$genelistFile"."
### Genenrate a statistics summary table read by python to create a table in a pdf file
if [ -f "$outputDir/tmp/$tmpfd/statistics.txt" ]; then rm "$outputDir/tmp/$tmpfd/statistics.txt"; fi
echo -e "Total number of variants in the raw VCF file $num0" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants left after the filter QUAL >= $minQual, DP >=$minReadDepth, MQ >= $minMapQual $num0a" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants after decomposing and left normalization $num" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants reported in dbSNP database $num2" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants reported in COSMIC database $num3" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants reported in both dbSNP and COSMIC database $((num3a-num+num2))" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants remaining after removing variants reported in dbSNP while keeping variants in COSMIC $num3a"  >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants (out of $num3a variants) that are nonsynonymous or splicing ones $num4" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants (out of $num4 variants) that are reported in COSMIC $num6" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of genes associated with $num4 variants $num5" >> "$outputDir/tmp/$tmpfd/statistics.txt"
## Genenrate an effect type summary table read by python to create a table in a pdf file
if [ -f "$outputDir/tmp/$tmpfd/effectType.txt" ]; then rm "$outputDir/tmp/$tmpfd/effectType.txt"; fi
echo -e "Number of variants (out of $num4 variants) that are frameshift deletion $num_del" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are frameshift insertion $num_ins" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stoploss $num_sl" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stopgain $num_sg" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are mis-sense $num_non" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are splicing $num_spg" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Total number of variants that are nonsynonymous or splicing $num4" >> "$outputDir/tmp/$tmpfd/effectType.txt"
## A file saving the directories of outputted files and the raw VCF file
if [ -f "$outputDir/tmp/$tmpfd/fileDetails.txt" ]; then rm "$outputDir/tmp/$tmpfd/fileDetails.txt"; fi
echo -e "$inputVCF" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
echo -e "$outputDir/$genelistFile" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
echo -e "$outputDir/$annoTable" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
echo -e "$outputDir/$annoVCF" >> "$outputDir/tmp/$tmpfd/fileDetails.txt"
## A file saving the filtering parameters
if [ -f "$outputDir/tmp/$tmpfd/filterParam.txt" ]; then rm "$outputDir/tmp/$tmpfd/filterParam.txt"; fi
echo -e "minimum quality $minQual" >> "$outputDir/tmp/$tmpfd/filterParam.txt"
echo -e "minimum read depth $minReadDepth" >> "$outputDir/tmp/$tmpfd/filterParam.txt"
echo -e "minimum mapping quality $minMapQual" >> "$outputDir/tmp/$tmpfd/filterParam.txt"
## run python code to generate summary report
outReportName="$outputFileNameShare"
outReportName+='_summaryReport'
cd "$outputDir/tmp/$tmpfd"
R -e "rmarkdown::render('$RcodeDir/summaryreportANNOVARhtml.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
R -e "rmarkdown::render('$RcodeDir/summaryreportANNOVARpdf.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
cd -
mv "$outputDir/tmp/$tmpfd/summaryreportANNOVARpdf.pdf" "$outputDir/$outReportName.pdf"
mv "$outputDir/tmp/$tmpfd/summaryreportANNOVARhtml.html" "$outputDir/$outReportName.html"
if (( $count > 0 ))
then
	perl -pe 's/^([^#])/chr\1/' "$outputDir/tmp/$tmpfd/$annoVCF" > "$outputDir/tmp/$tmpfd/tmp.vcf"
	mv "$outputDir/tmp/$tmpfd/tmp.vcf" "$outputDir/$annoVCF"
else
	mv "$outputDir/tmp/$tmpfd/$annoVCF" "$outputDir/$annoVCF"
fi
mv "$outputDir/tmp/$tmpfd/$annoTable" "$outputDir/$annoTable"
mv "$outputDir/tmp/$tmpfd/tmp_gl.txt" "$outputDir/$genelistFile"
echo "An annotation table associated with $num4 variants is saved in $outputDir/$annoTable, and a vcf file with annotation information is saved in $outputDir/$annoVCF."
rm -r -f "$outputDir/tmp/$tmpfd"
echo completed `date +'%Y-%m-%d %T'`
