#!/bin/bash
set -e
seqtools_snpeff=/opt/SeqTools/bin/snpEff
export bdge_bcftools_PATH=/opt/SeqTools/bin/bcftools-1.3
export PATH=$bdge_bcftools_PATH:$PATH

export bdge_htslib_PATH=/opt/SeqTools/bin/samtools-1.3/htslib-1.3
export PATH=$bdge_htslib_PATH:$PATH

inputVCF="/home/brb/SeqTestdata/RNASeqFibroblast/outputhg38/LFB_scramble_repA_raw.vcf"
cosmicVCF="/home/brb/SeqTestdata/usefulvcf/hg38/CosmicCodingMuts.vcf.gz"
dbSNPVCF="/home/brb/SeqTestdata/usefulvcf/hg38/common_all_20151104.vcf.gz"
genomeRef="/home/brb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
dbnsfpFile="/home/brb/SeqTestdata/usefulvcf/hg38/dbNSFP.txt.gz"
minQual=1
minReadDepth=1
minMapQual=1
genomeVer=hg38
RcodeDir="/opt/SeqTools/bin/SeqTools/code"
outputDir="/home/brb/SeqTestdata/RNASeqFibroblast/outputhg38_cli_snpeff"
tmpname="$(basename $inputVCF)"
outputFileNameShare="${tmpname::-4}"
tmpfd="tmp_$outputFileNameShare"
if [ ! -d "$outputDir/tmp/$tmpfd" ]; then mkdir -p "$outputDir/tmp/$tmpfd";  fi

echo "running sample: $inputVCF" >> "$outputDir/tmp/log.txt"
bcftools filter -i"QUAL >= $minQual && DP >= $minReadDepth && MQ >= $minMapQual" "$inputVCF" > "$outputDir/tmp/$tmpfd/filtered.vcf"
(bcftools norm -m-both -o "$outputDir/tmp/$tmpfd/splitted.vcf" "$outputDir/tmp/$tmpfd/filtered.vcf") 2>&1
(bcftools norm -f $genomeRef -o "$outputDir/tmp/$tmpfd/leftnormalized.vcf" "$outputDir/tmp/$tmpfd/splitted.vcf") 2>&1
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
(java -Xmx4G -jar "$seqtools_snpeff/snpEff.jar" -canon -no-downstream -no-upstream -no-intergenic -no-intron -no-utr -noNextProt -noMotif $genomeVer -s "$outputDir/tmp/$tmpfd/annodbsnpRemove.html" "$outputDir/tmp/$tmpfd/cosmic_dbsnp_rem.vcf" > "$outputDir/tmp/$tmpfd/snpeff_anno.vcf") 2>&1 
(cat "$outputDir/tmp/$tmpfd/snpeff_anno.vcf" | java -jar "$seqtools_snpeff/SnpSift.jar" filter "(ANN[*].BIOTYPE = 'protein_coding') | (ANN[*].EFFECT has 'splice')"  > "$outputDir/tmp/$tmpfd/snpeff_proteincoding.vcf") 2>&1
# Remove synonymous variants
grep -Ev 'synonymous_variant|start_retained|stop_retained_variant' "$outputDir/tmp/$tmpfd/snpeff_proteincoding.vcf" > "$outputDir/tmp/$tmpfd/nonsyn_splicing.vcf"
grep -wv 'LOW' "$outputDir/tmp/$tmpfd/nonsyn_splicing.vcf" > "$outputDir/tmp/$tmpfd/nonsyn_splicing2.vcf"
if [[ $genomeVer == hg38 ]] || [[ $genomeVer == GRCh38* ]]
then
dbnsfpField=SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,CADD_raw,CADD_phred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,GERP++_NR,GERP++_RS,phyloP7way_vertebrate,phastCons7way_vertebrate,SiPhy_29way_logOdds 
else
dbnsfpField=SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,CADD_raw,CADD_phred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,GERP++_NR,GERP++_RS,phyloP100way_vertebrate,phastCons100way_vertebrate,SiPhy_29way_logOdds 
fi
java -jar "$seqtools_snpeff/SnpSift.jar" dbNSFP -f $dbnsfpField -v -db "$dbnsfpFile" "$outputDir/tmp/$tmpfd/nonsyn_splicing2.vcf" > "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf"
if [[ $genomeVer == hg38 ]] || [[ $genomeVer == GRCh38* ]]
then
java -jar "$seqtools_snpeff/SnpSift.jar" extractFields -s "," -e "." "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf" CHROM POS ID REF ALT "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].GENEID" "ANN[0].FEATURE" "ANN[0].FEATUREID" "ANN[0].BIOTYPE" "ANN[0].HGVS_C" "ANN[0].HGVS_P" dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_score dbNSFP_LRT_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_FATHMM_score dbNSFP_FATHMM_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_VEST3_score dbNSFP_CADD_raw dbNSFP_CADD_phred dbNSFP_MetaSVM_score dbNSFP_MetaSVM_pred dbNSFP_MetaLR_score dbNSFP_MetaLR_pred dbNSFP_GERP___NR dbNSFP_GERP___RS dbNSFP_phyloP7way_vertebrate dbNSFP_phastCons7way_vertebrate dbNSFP_SiPhy_29way_logOdds > "$outputDir/tmp/$tmpfd/annoTable.txt"
else
java -jar "$seqtools_snpeff/SnpSift.jar" extractFields -s "," -e "." "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf" CHROM POS ID REF ALT "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].GENEID" "ANN[0].FEATURE" "ANN[0].FEATUREID" "ANN[0].BIOTYPE" "ANN[0].HGVS_C" "ANN[0].HGVS_P" dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_score dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_score dbNSFP_LRT_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_FATHMM_score dbNSFP_FATHMM_pred dbNSFP_PROVEAN_score dbNSFP_PROVEAN_pred dbNSFP_VEST3_score dbNSFP_CADD_raw dbNSFP_CADD_phred dbNSFP_MetaSVM_score dbNSFP_MetaSVM_pred dbNSFP_MetaLR_score dbNSFP_MetaLR_pred dbNSFP_GERP___NR dbNSFP_GERP___RS dbNSFP_phyloP100way_vertebrate dbNSFP_phastCons100way_vertebrate dbNSFP_SiPhy_29way_logOdds > "$outputDir/tmp/$tmpfd/annoTable.txt"
fi
annoTable="$outputFileNameShare"
annoTable+='_annoTable.txt'
sed -r 's/(\[|\])//g' <"$outputDir/tmp/$tmpfd/annoTable.txt" > "$outputDir/tmp/$tmpfd/tmp.txt"
sed -r 's/\#//g' <"$outputDir/tmp/$tmpfd/tmp.txt" > "$outputDir/tmp/$tmpfd/tmp2.txt"
sed -r 's/dbNSFP\_//g' <"$outputDir/tmp/$tmpfd/tmp2.txt" > "$outputDir/tmp/$tmpfd/tmp3.txt"
sed -r 's/ANN0\.//g' <"$outputDir/tmp/$tmpfd/tmp3.txt" > "$outputDir/tmp/$tmpfd/$annoTable"
##run annotation for variants after pre-filtering
(java -Xmx4G -jar "$seqtools_snpeff/snpEff.jar" -canon -noNextProt -noMotif $genomeVer -s "$outputDir/tmp/$tmpfd/annoPrefilter.html" "$outputDir/tmp/$tmpfd/leftnormalized.vcf" > "$outputDir/tmp/$tmpfd/snpEffAnnoAll.vcf") 2>&1 | tee -a  $outputDir/tmp/$tmpfd/log.txt
java -jar "$seqtools_snpeff/SnpSift.jar" extractFields -s "," -e "." "$outputDir/tmp/$tmpfd/snpEffAnnoAll.vcf" "ANN[0].EFFECT" > "$outputDir/tmp/$tmpfd/annoTableAll.txt"
($RcodeDir/./createGeneListSNPEFF.sh "$outputDir/tmp/$tmpfd/$annoTable" "$outputDir/tmp/$tmpfd/annoTableAll.txt" "$outputDir/tmp/$tmpfd") 2>&1 | tee -a "$outputDir/tmp/$tmpfd/log.txt"
rm "$outputDir/tmp/$tmpfd/$annoTable"
grep -v "protein_protein_contact" "$outputDir/tmp/$tmpfd/annoTableAfterClean.txt" > "$outputDir/tmp/$tmpfd/$annoTable"
annoVCF="$outputFileNameShare"
annoVCF+='_annotated.vcf'
grep -v "protein_protein_contact" "$outputDir/tmp/$tmpfd/nonsyn_splicing_dbnsfp.vcf" > "$outputDir/tmp/$tmpfd/$annoVCF"
genelistFile="$outputFileNameShare"
genelistFile+='_genelist.txt'
mv "$outputDir/tmp/$tmpfd/genelist.txt" "$outputDir/tmp/$tmpfd/$genelistFile"
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
read num4 <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/$annoVCF" | wc -l)
read num5 <<< $(grep -v 'Mutation' "$outputDir/tmp/$tmpfd/$genelistFile" | wc -l)
echo -e "There are $num4 variants (out of $num3a variants) that are nonsynonymous or splicing ones, which are kept for further analysis."
read num6 <<< $(grep -v '#' "$outputDir/tmp/$tmpfd/$annoVCF" | grep "COSM" | wc -l)
echo -e "There are $num6 variants (out of $num4 variants) that are reported by COSMIC v74."
read num_spg <<< $(grep 'splice_region_variant' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
echo -e "There are $num_spg variants (out of $num4 variants) that are splicing."
read num_del <<< $(grep 'frameshift_variant' "$outputDir/tmp/$tmpfd/$annoTable" | wc -l)
echo -e "There are $num_del variants (out of $num4 variants) that are frameshift."
read num_sl <<< $(grep 'stop_lost' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) 
echo -e "There are $num_sl variants (out of $num4 variants) that are stoploss."
read num_sg <<< $(grep 'stop_gained' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l)
echo -e "There are $num_sg variants (out of $num4 variants) that are stopgain." 
read num_stl <<< $(grep 'start_lost' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l) >> "$outputDir/tmp/log.txt"
echo -e "There are $num_stl variants (out of $num4 variants) that are startloss." >> "$outputDir/tmp/log.txt"
read num_non <<< $(grep 'missense_variant' "$outputDir/tmp/$tmpfd/$annoTable"  | wc -l)
echo -e "There are $num_non variants (out of $num4 variants) that are mis-sense."
echo -e "There are $num5 genes that are associated with $num4 variants, and the gene list is saved in $outputDir/tmp/$tmpfd/$genelistFile."
### Genenrate a statistics summary table read by python to create a table in a pdf file
if [ -f "$outputDir/tmp/$tmpfd/statistics.txt" ]; then rm "$outputDir/tmp/$tmpfd/statistics.txt"; fi
echo -e "Total number of variants in the raw VCF file $num0" >> "$outputDir/tmp/$tmpfd/statistics.txt"
echo -e "Number of variants left after the filter by QUAL >= $minQual, DP >=$minReadDepth, MQ >= $minMapQual $num0a" >> "$outputDir/tmp/$tmpfd/statistics.txt"
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
echo -e "Number of variants (out of $num4 variants) that are frameshift $num_del" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stoploss $num_sl" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are stopgain $num_sg" >> "$outputDir/tmp/$tmpfd/effectType.txt"
echo -e "Number of variants (out of $num4 variants) that are startloss $num_stl" >> "$outputDir/tmp/$tmpfd/effectType.txt"
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
outReportName=$outputFileNameShare
outReportName+='_summaryReport'
cd "$outputDir/tmp/$tmpfd"
R -e "rmarkdown::render('$RcodeDir/summaryreportSnpEffhtml.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
R -e "rmarkdown::render('$RcodeDir/summaryreportSnpEffpdf.Rmd',intermediates_dir='$outputDir/tmp/$tmpfd/',output_dir='$outputDir/tmp/$tmpfd/')" --args "$outputDir/tmp/$tmpfd/"
cd -
mv "$outputDir/tmp/$tmpfd/summaryreportSnpEffpdf.pdf" "$outputDir/$outReportName.pdf"
mv "$outputDir/tmp/$tmpfd/summaryreportSnpEffhtml.html" "$outputDir/$outReportName.html"
if (( $count > 0 ))
then
	perl -pe 's/^([^#])/chr\1/' "$outputDir/tmp/$tmpfd/$annoVCF" > "$outputDir/tmp/$tmpfd/tmp.vcf"
	mv "$outputDir/tmp/$tmpfd/tmp.vcf" "$outputDir/$annoVCF"
else
	mv "$outputDir/tmp/$tmpfd/$annoVCF" "$outputDir/$annoVCF"
fi
mv "$outputDir/tmp/$tmpfd/$annoTable" "$outputDir/$annoTable"
mv "$outputDir/tmp/$tmpfd/$genelistFile" "$outputDir/$genelistFile"
echo "An annotation table associated with $num4 variants is saved in $outputDir/$annoTable, and a vcf file with annotation information is saved in $outputDir/$annoVCF." >> $outputDir/log.txt
rm -r -f "$outputDir/tmp/$tmpfd"
echo completed `date +'%Y-%m-%d %T'`
