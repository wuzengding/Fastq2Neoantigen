##########################################################################
# File Name:        MC38.analysis.sh
# Author:           wuzengding
# mail:             wuzengding@163.com
# Created Time:     Mon 19 Sep 2022 01:32:15 PM CST
##########################################################################

tumor_fq1=$1
tumor_fq2=$2 # if single sequence,fastq2 must with input of "None"
normal_fq1=$3
normal_fq2=$4 # if single sequence,fastq2 must with input of "None"
outpath=$5
tumorname=$6
normalname=$7
xomic=$8

if [ -z "${xomic}" ]; then
    echo "paragram exited \
          you should give ensential parameters!"
    exit 1
else
    echo tumor_fq1:"$tumor_fq1"
    echo tumor_fq2="$tumor_fq2"
    echo normal_fq1="$normal_fq1"
    echo normal_fq2="$normal_fq2"
    echo outpath="$outpath"
    echo tumorname="$tumorname"
    echo normalname="$normalname"
    echo xomic="$xomic"
    
fi

###################################################################
#########         config                              #############
###################################################################
fasp=/mnt/data2/wuzengding/03.biotools/software/fastp
STAR=/mnt/data2/wuzengding/03.biotools/miniconda3/bin/STAR
genomedir=/mnt/data2/wuzengding/00.database/15.index_mouse_genome/index_GRCm39_ensemble107.star
bwa=/mnt/data2/wuzengding/03.biotools/software/bwa/bwa
index_bwa=/mnt/data2/wuzengding/00.database/15.index_mouse_genome/index_GRCm39_ensemble107.bwa/index_Mus_musculus.GRCm39.dna.primary_assembly_ensemble107
sambamba=/mnt/data2/wuzengding/03.biotools/miniconda3/bin/sambamba
VarDict=/mnt/data2/wuzengding/03.biotools/software/VarDictJava/build/install/VarDict/lib/VarDict-1.8.3.jar
geneome=/mnt/data2/wuzengding/00.database/15.index_mouse_genome/Mus_musculus.GRCm39.dna.primary_assembly.fa
bedfile=/mnt/data2/wuzengding/00.database/16.bedfile_mouse/mm39.primary_chr.bed
testsomaticR=/mnt/data2/wuzengding/03.biotools/software/VarDictJava/build/install/VarDict/bin/testsomatic.R
var2vcf_pairedpl=/mnt/data2/wuzengding/03.biotools/software/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl
SnpSift=/mnt/data2/wuzengding/03.biotools/software/snpEff/SnpSift.jar
vepCacheDir=/mnt/data2/wuzengding/03.biotools/software/ensembl-vep/vep_data


if false;then
######### step1  fasp trim and qc for raw fastq reads #############
###################################################################
mkdir -p ${outpath}/${tumorname}/01.QC
mkdir -p ${outpath}/${normalname}/01.QC
if [ $tumor_fq2 != "None" ];then
    echo "Your inputs of tumor sample is PE sequence, start fastp for qc"
    ${fasp} --thread 10 -i ${tumor_fq1} -I ${tumor_fq2} -o  ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz -O ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz --unpaired1 ${outpath}/${tumorname}/01.QC/${tumorname}_R1.unpaired.fastq.gz --unpaired2 ${outpath}/${tumorname}/01.QC/${tumorname}_R2.unpaired.fastq.gz --html ${outpath}/${tumorname}/01.QC/${tumorname}.html --json ${outpath}/${tumorname}/01.QC/${tumorname}.json 
elif [ $tumor_fq2 = "None" ];then
    echo "Your inputs of tumor sample is SE sequence,"
    ${fasp} --thread 10 -i ${tumor_fq1} -o  ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz  --unpaired1 ${outpath}/${tumorname}/01.QC/${tumorname}_R1.unpaired.fastq.gz  --html ${outpath}/${tumorname}/01.QC/${tumorname}.html --json ${outpath}/${tumorname}/01.QC/${tumorname}.json 
fi

if [ $normal_fq2 != "None" ];then
    echo "Your inputs of normal sample is PE sequence,start fastp for qc"
    ${fasp} --thread 10 -i ${normal_fq1} -I ${normal_fq2} -o ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz -O ${outpath}/${normalname}/01.QC/${normalname}_R2.fastq.gz --unpaired1 ${outpath}/${normalname}/01.QC/${normalname}_R1.unpaired.fastq.gz --unpaired2 ${outpath}/${normalname}/01.QC/${normalname}_R2.unpaired.fastq.gz --html ${outpath}/${normalname}/01.QC/${normalname}.html --json ${outpath}/${normalname}/01.QC/${normalname}.json
elif [ $normal_fq2 = "None" ];then
    echo "Your inputs of normal sample is SE sequence, start fastp for qc"
    ${fasp} --thread 10 -i /mnt/data2/cuifenglei/data/normal_mouse_tissue/${normalname}_pass_1.fastq.gz  -o ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz  --unpaired1 ${outpath}/${normalname}/01.QC/${normalname}_R1.unpaired.fastq.gz  --html ${outpath}/${normalname}/01.QC/${normalname}.html --json ${outpath}/${normalname}/01.QC/${normalname}.json
fi
    
fi
###################################################################
###########   step2  STAR align for RN                #############
###################################################################
mkdir -p ${outpath}/${tumorname}/02.Aln
mkdir -p ${outpath}/${normalname}/02.Aln

if [ $xomic = "RNA" ];then
    if [ $tumor_fq2 != "None" ];then 
        echo "Your inputs of tumor sample is PE sequence,start alignning"
        ${STAR}  --runThreadN 18  --genomeDir ${genomedir} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${tumorname}/02.Aln/${tumorname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
    else
        echo "Your inputs of tumor sample is SE sequence,start alignning"
        ${STAR}  --runThreadN 18  --genomeDir ${genomedir} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz  --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${tumorname}/02.Aln/${tumorname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
    fi
    
    if [ $normal_fq2 != "None" ];then 
        echo "Your inputs of normal sample is PE sequence,start alignning"
        ${STAR}  --runThreadN 18  --genomeDir ${genomedir} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz ${outpath}/${normalname}/01.QC/${normalname}_R2.fastq.gz --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${normalname}/02.Aln/${normalname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
    else
        echo "Your inputs of normal sample is PE sequence,start alignning"
        ${STAR}  --runThreadN 18  --genomeDir ${genomedir} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz  --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${normalname}/02.Aln/${normalname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
    fi
    ${sambamba} index -t 18 ${outpath}/${tumorname}/02.Aln/${tumorname}Aligned.sortedByCoord.out.bam  > ${outpath}/${tumorname}/02.Aln/${tumorname}Aligned.sortedByCoord.out.bam.bai
    ${sambamba} index -t 18 ${outpath}/${normalname}/02.Aln/${normalname}Aligned.sortedByCoord.out.bam > ${outpath}/${normalname}/02.Aln/${normalname}Aligned.sortedByCoord.out.bam.bai
###################################################################
###########   step2.2  BWA align for DNA              #############
###################################################################
elif [ $xomic = "DNA" ];then
    echo "your inputs are DNA"
    if [ $tumor_fq2 != "None" ];then 
        echo ",start alignning"
        ${bwa} mem -t 20 ${index_bwa} ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sam
    else
        echo ",start alignning"
        ${bwa} mem -t 20 ${index_bwa} ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sam
    fi
    if [ $normal_fq2 != "None" ];then
        echo "start alignning"
        ${bwa} mem -t 20 ${index_bwa} ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz ${outpath}/${normalname}/01.QC/${normalname}_R2.fastq.gz > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sam
    else
        echo "start alignning"
        ${bwa} mem -t 20 ${index_bwa} ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sam
    fi
    ${sambamba} view-t 20 --format=bam ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sam > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.bam
    ${sambamba} sort -t 20 ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.bam > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sored.bam
    ${sambamba} index -t 20 ${outpath}/${tumorname}/02.Aln/${tumorname}Aligned.sortedByCoord.out.bam  > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sored.bam.bai
    
    ${sambamba} view-t 20 --format=bam ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sam > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.bam
    ${sambamba} sort -t 20 ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.bam > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sored.bam
    ${sambamba} index -t 20 ${outpath}/${normalname}/02.Aln/${normalname}Aligned.sortedByCoord.out.bam  > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sored.bam.bai
fi
###################################################################
###########   step3  index for BAM                    #############
###################################################################



###################################################################
###########   step4  SNV calling                      #############
###################################################################
mkdir -p ${outpath}/${tumorname}/03.SNV
/usr/bin/java -jar -Xmx20g ${VarDict} -G ${geneome} -f 0.002 -N ${tumorname} -b "${outpath}/${tumorname}/02.Aln/${tumorname}Aligned.sortedByCoord.out.bam|${outpath}/${normalname}/02.Aln/${normalname}Aligned.sortedByCoord.out.bam"  -c 1 -S 2 -E 3 -g 4 ${bedfile} -th 40 | ${testsomaticR} | ${var2vcf_pairedpl} -N "${tumorname}|${normalname}"  -f 0.002 >${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.vcf


###################################################################
###########   step5  VEP annotation for vcf           #############
###################################################################

docker run --rm -v ${vepCacheDir}:/opt/vep/cache -v ${outpath}/${tumorname}/03.SNV:${temp}/03.SNV  ensemblorg/ensembl-vep  vep  --fork 20 --input_file  ${temp}/03.SNV/${tumorname}.vardict.vcf --output_file  ${temp}/03.SNV/${tumorname}.vardict.vep.vcf --format vcf --vcf --symbol --species mus_musculus  --terms SO --tsl --offline --cache  --dir_cache /opt/vep/cache  --plugin Frameshift --plugin Wildtype


###################################################################
#########     step6   somatic and AF filter           #############
###################################################################
/usr/bin/java -Xmx10g -jar ${SnpSift} filter "( na FILTER ) | (FILTER = 'PASS') && (QUAL >=30) && (AF >= 0.1) && ( DP >= 5 )" --file ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.vep.vcf |grep -v Germline > ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.vep.somatic.vcf 

###################################################################
###########   step7  neoantigen prediction            #############
###################################################################

docker run --rm -v ${outpath}/${tumorname}:/tmp/${tumorname}  griffithlab/pvactools pvacseq run  --n-threads 25 --iedb-install-directory /opt/iedb --allele-specific-binding-thresholds --keep-tmp-files  --minimum-fold-change 1 /tmp/${tumorname}/03.SNV/${tumorname}.vardict.vep.somatic.vcf   ${tumorname} H-2-Kb,H-2-Db,H2-IAb  NetMHCIIpan NetMHCpan /tmp/${tumorname}/05.AntigenPrediction  -e1 8,9,10,11 -e2 15 


docker run --rm -v ${outpath}/${tumorname}/05.AntigenPrediction/combined:/tmp/combined  griffithlab/pvactools pvacseq generate_aggregated_report /tmp/combined/${tumorname}.filtered.tsv /tmp/combined/${tumorname}.filtered.temp.tsv
head -1 ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.filter.temp.tsv > ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.filtered.aggregated.tsv
tail -n+2  ${outpath}/${tumorname}/05.AntigenPrediction/netMHC/MHC_Class_I/${tumorname}.filtered.temp.tsv|cut -f3,4|while read gene_variant;
	do
		gene_variant=$(echo ${gene_variant} | tr ' ' '\t')
		echo ${gene_variant}
		grep "${gene_variant}" ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.all_epitopes.aggregated.tsv >> ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.filtered.aggregated.tsv
	done

## Fusion calling
#export PATH=/mnt/data2/wuzengding/03.biotools/software/STAR-2.7.2d/bin/Linux_x86_64:$PATH
#/mnt/data2/wuzengding/03.biotools/software/STAR-Fusion.v1.9.0/STAR-Fusion  --genome_lib_dir /mnt/data2/wuzengding/03.biotools/software/STAR-Fusion.v1.9.0/Mouse_gencode_M24_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir --left_fq ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz --right_fq ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz --output_dir ${outpath}/${tumorname}/04.Fusion


## generate protein.fasta
#/usr/bin/java -Xmx10g -jar /mnt/data2/wuzengding/03.biotools/software/snpEff/snpEff.jar -c /mnt/data2/wuzengding/03.biotools/software/snpEff/snpEff.config -noLog -v -noInteraction -nextProt mm39  -fastaProt ${outpath}/${tumorname}/${tumorname}.somatic.rot.fa ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.somatic.filter.vcf > ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.somatic.annot.vcf

## epitope pridiction
#/mnt/data1/feicaiyi/netmhc/netpan41/netMHCpan -a H-2-Kb,H-2-Db -s -f  ${outpath}/${tumorname}/03.SNV/${tumorname}.somatic.prot.fa -xls -xlsfile ${outpath}/${tumorname}/05.AntigenPrediction/${tumorname}.NetMHCpan_out.xls -inptype 0 -BA


