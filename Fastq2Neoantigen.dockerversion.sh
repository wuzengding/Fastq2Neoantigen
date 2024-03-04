##########################################################################
# File Name:        Fastq2Neoantigen.sh
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
xomic=$8  # "RNA" or "DNA"
species=$9 # "mus_musculus" or "homo_sapiens"
runmodel=${10} # "all" or model name
mhctype=${11} # example "H-2-Kb,H-2-Db,H2-IAb" 或者"HLA-A11:01,HLA-A26:01,HLA-B39:01,HLA-B57:01"等方式

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
    echo species=$species
    echo runmodel=$runmodel
    echo mhctype=${mhctype}
fi

###################################################################
#########         config                              #############
###################################################################
##这个配置文件在初次配置时，或者环境变换时，需要重新配置
fasp=/usr/local/bin/fastp
STAR=/usr/local/bin/STAR
sambamba=/usr/local/bin/sambamba
bwa=/usr/local/bin/bwa
VarDict=/usr/local/bin/VarDict.jar
testsomaticR=/usr/local/bin/testsomatic.R
var2vcf_pairedpl=/usr/local/bin/var2vcf_paired.pl
SnpSift=/usr/local/bin/SnpSift.jar
vepCacheDir=/mnt/user/wzd/00.database/21.ensemble_vep
vepPlugins=/mnt/user/wzd/00.database/21.ensemble_vep/Plugins
realpath=/mnt/user/wzd
realoutpath=$(echo $outpath|sed "s:/opt:${realpath}:g")
userid=2001
groupid=1000

if [ "${species}" == "mus_musculus" ];then
  annotdir=/mnt/user/wzd/00.database/15.index_mouse_genome
  fasta=Mus_musculus.GRCm39.dna.toplevel.fa
  index_star=/opt/00.database/15.index_mouse_genome/index_GRCm39_ensemble110.star
  index_bwa=/opt/00.database/15.index_mouse_genome/index_GRCm39_ensemble110.bwa/index_Mus_musculus.GRCm39.dna.toplevel
  genome=/opt/00.database/15.index_mouse_genome/Mus_musculus.GRCm39.dna.toplevel.fa
  bedfile=/opt/00.database/16.annot_mouse/Mus_musculus.GRCm39.107.bed
 
elif [ "${species}" == "homo_sapiens" ];then
  annotdir=/mnt/user/wzd/00.database/17.annot_human
  fasta=hg38_SIRV_ERCC_longSIRV.fa
  index_star=/opt/00.database/17.annot_human/index_hg38_SIRV_ERCC_STAR
  index_bwa=/opt/00.database/17.annot_human/index_hg38_SIRV_ERCC_BWA/hg38_SIRV_ERCC_longSIRV
  genome=/opt/00.database/17.annot_human/hg38_SIRV_ERCC_longSIRV.fa
  bedfile=/opt/00.database/bed_hg38_hg19/hg38_Gencode_V28.bed
fi




###################################################################
######### step1  fasp trim and qc for raw fastq reads #############
###################################################################
if [[ "${runmodel}" == *"fastp"* ]];then
    echo "分析包含模块fastp"
    setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "mkdir -p ${outpath}/${tumorname}/01.QC"
    setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "mkdir -p ${outpath}/${normalname}/01.QC"
    if [ $tumor_fq2 != "None" ];then
        echo "Your inputs of tumor sample is PE sequence, start fastp for qc"
      setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${fasp} --thread 10 -i ${tumor_fq1} -I ${tumor_fq2} -o  ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz -O ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz --unpaired1 ${outpath}/${tumorname}/01.QC/${tumorname}_R1.unpaired.fastq.gz --unpaired2 ${outpath}/${tumorname}/01.QC/${tumorname}_R2.unpaired.fastq.gz --html ${outpath}/${tumorname}/01.QC/${tumorname}.html --json ${outpath}/${tumorname}/01.QC/${tumorname}.json"
    elif [ $tumor_fq2 = "None" ];then
        echo "Your inputs of tumor sample is SE sequence,"
      setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c  "${fasp} --thread 10 -i ${tumor_fq1} -o  ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz  --unpaired1 ${outpath}/${tumorname}/01.QC/${tumorname}_R1.unpaired.fastq.gz  --html ${outpath}/${tumorname}/01.QC/${tumorname}.html --json ${outpath}/${tumorname}/01.QC/${tumorname}.json"
    fi
    
    if [ $normal_fq2 != "None" ];then
        echo "Your inputs of normal sample is PE sequence,start fastp for qc"
       setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${fasp} --thread 10 -i ${normal_fq1} -I ${normal_fq2} -o ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz -O ${outpath}/${normalname}/01.QC/${normalname}_R2.fastq.gz --unpaired1 ${outpath}/${normalname}/01.QC/${normalname}_R1.unpaired.fastq.gz --unpaired2 ${outpath}/${normalname}/01.QC/${normalname}_R2.unpaired.fastq.gz --html ${outpath}/${normalname}/01.QC/${normalname}.html --json ${outpath}/${normalname}/01.QC/${normalname}.json"
    elif [ $normal_fq2 = "None" ];then
        echo "Your inputs of normal sample is SE sequence, start fastp for qc"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c  "${fasp} --thread 10  -i ${normal_fq1}  -o ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz  --unpaired1 ${outpath}/${normalname}/01.QC/${normalname}_R1.unpaired.fastq.gz  --html ${outpath}/${normalname}/01.QC/${normalname}.html --json ${outpath}/${normalname}/01.QC/${normalname}.json"
    fi
fi 

###################################################################
###########   step2  align for fq              #############
###################################################################

if [[ "${runmodel}" == *"align"* ]];then
    echo "分析包含模块align"
    setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c  "mkdir -p ${outpath}/${tumorname}/02.Aln"
    setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c  "mkdir -p ${outpath}/${normalname}/02.Aln"

    ###########   step2  STAR align for RNA  ######################
    if [ $xomic = "RNA" ];then
        if [ $tumor_fq2 != "None" ];then 
            echo "Your inputs of tumor sample is PE sequence,start alignning"
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups   ${STAR}  --runThreadN 18  --genomeDir ${index_star} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${tumorname}/02.Aln/${tumorname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
        else
            echo "Your inputs of tumor sample is SE sequence,start alignning"
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups  ${STAR}  --runThreadN 18  --genomeDir ${index_star} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz  --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${tumorname}/02.Aln/${tumorname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
        fi
        
        if [ $normal_fq2 != "None" ];then 
            echo "Your inputs of normal sample is PE sequence,start alignning"
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups  ${STAR}  --runThreadN 18  --genomeDir ${index_star} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz ${outpath}/${normalname}/01.QC/${normalname}_R2.fastq.gz --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${normalname}/02.Aln/${normalname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
        else
            echo "Your inputs of normal sample is PE sequence,start alignning"
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups  ${STAR}  --runThreadN 18  --genomeDir ${index_star} --outReadsUnmapped Fastx --readFilesIn ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz  --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${outpath}/${normalname}/02.Aln/${normalname} --twopassMode Basic --chimSegmentMin 12 --outSAMstrandField intronMotif --chimJunctionOverhangMin 8 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --outSAMunmapped Within --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 --genomeLoad NoSharedMemory --alignInsertionFlush Right --alignSplicedMateMapLmin 30 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --chimScoreJunctionNonGTAG -4 --alignSplicedMateMapLminOverLmate 0 --alignSJstitchMismatchNmax 5 -1 5 5 --readFilesCommand 'gunzip -c' --chimOutType Junctions WithinBAM SoftClip
        fi

    ################   step2.2  BWA align for DNA    ###################

    elif [ $xomic = "DNA" ];then
        echo "your inputs are DNA"
        
        if [ $tumor_fq2 != "None" ];then 
            echo "start alignning"
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${bwa} mem -t 20 ${index_bwa} ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sam"
        else
            echo "start alignning"
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${bwa} mem -t 20 ${index_bwa} ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sam"
        fi
        
        if [ $normal_fq2 != "None" ];then
            echo "start alignning"
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c  "${bwa} mem -t 20 ${index_bwa} ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz ${outpath}/${normalname}/01.QC/${normalname}_R2.fastq.gz > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sam"
        else
            echo "start alignning"         
            setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${bwa} mem -t 20 ${index_bwa} ${outpath}/${normalname}/01.QC/${normalname}_R1.fastq.gz > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sam"    
        fi
    fi
fi

###################################################################
###########   step3  index for BAM                    #############
###################################################################
if [[ "${runmodel}" == *"index"* ]];then
    echo "分析包含模块index"
    if [ $xomic = "RNA" ];then
        echo "setpriv --reuid=${userid} --regid=${groupid} --clear-groups  bash -c ${sambamba} index -t 18 ${outpath}/${tumorname}/02.Aln/${tumorname}Aligned.sortedByCoord.out.bam"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${sambamba} index -t 18 ${outpath}/${tumorname}/02.Aln/${tumorname}Aligned.sortedByCoord.out.bam"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${sambamba} index -t 18 ${outpath}/${normalname}/02.Aln/${normalname}Aligned.sortedByCoord.out.bam"
    
    elif [ $xomic = "DNA" ];then
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${sambamba} view -t 20 --sam-input --format=bam ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sam >${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.bam"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${sambamba} sort -t 20 ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.bam > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sorted.bam"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${sambamba} index -t 20 ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sorted.bam  > ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sorted.bam.bai"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "rm -f ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.bam"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "rm -f ${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sam"
          
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c  "${sambamba} view -t 20 --sam-input --format=bam ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sam > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.bam"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c  "${sambamba} sort -t 20 ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.bam > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sorted.bam"
       
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "${sambamba} index -t 20 ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sorted.bam  > ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sorted.bam.bai"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "rm -f ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.bam"
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "rm -f ${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sam"
    fi
fi

###################################################################
###########   step4  SNV calling                      #############
###################################################################
if [[ "${runmodel}" == *"snvcall"* ]];then
    echo "分析包含模块SNV calling"

    setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "mkdir -p ${outpath}/${tumorname}/03.SNV"
    if [ $xomic = "RNA" ];then
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "/usr/bin/java -jar -Xmx20g ${VarDict} -G ${genome} -f 0.002 -N ${tumorname} -b \"${outpath}/${tumorname}/02.Aln/${tumorname}Aligned.sortedByCoord.out.bam|${outpath}/${normalname}/02.Aln/${normalname}Aligned.sortedByCoord.out.bam\"  -c 1 -S 2 -E 3 -g 4 ${bedfile} -th 40 | ${testsomaticR} | ${var2vcf_pairedpl} -N \"${tumorname}|${normalname}\"  -f 0.002 >${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.vcf"
    elif [ $xomic = "DNA" ];then
        setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "/usr/bin/java -jar -Xmx20g ${VarDict} -G ${genome} -f 0.002 -N ${tumorname} -b \"${outpath}/${tumorname}/02.Aln/${tumorname}.BWA.aligned.sorted.bam|${outpath}/${normalname}/02.Aln/${normalname}.BWA.aligned.sorted.bam\"  -c 1 -S 2 -E 3 -g 4 ${bedfile} -th 40 | ${testsomaticR} | ${var2vcf_pairedpl} -N \"${tumorname}|${normalname}\"  -f 0.002 >${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.vcf"
    fi
fi

###################################################################
###########   step5  VEP annotation for vcf           #############
###################################################################
if [[ "${runmodel}" == *"vepannot"* ]];then
    echo "分析包含模块SNV VEP annotation"
    
    setpriv --reuid=${userid} --regid=${groupid} --clear-groups bash -c "mkdir -p ${outpath}/${tumorname}/04.FUSION"

    echo " docker run --user ${userid}:${groupid} --rm -v ${annotdir}:/annotdir -v ${vepCacheDir}:/opt/vep/cache  -v ${vepPlugins}:/opt/vepPlugins -v ${realoutpath}/${tumorname}/03.SNV:/temp/03.SNV  ensemblorg/ensembl-vep  vep  --input_file  /temp/03.SNV/${tumorname}.vardict.vcf  --output_file  /temp/03.SNV/${tumorname}.vardict.vep.vcf  --dir_plugins /opt/vepPlugins  --cache --dir_cache /opt/vep/cache      --fork 20 --format vcf --vcf  --symbol --species ${species}  --terms SO --tsl  --biotype --assembly --fasta /annotdir/${fasta}  --offline   --plugin Frameshift --plugin Wildtype  --pick --force_overwrite"

    docker run --user ${userid}:${groupid} --rm -v ${annotdir}:/annotdir -v ${vepCacheDir}:/opt/vep/cache  -v ${vepPlugins}:/opt/vepPlugins -v ${realoutpath}/${tumorname}/03.SNV:/temp/03.SNV  ensemblorg/ensembl-vep  vep  --input_file  /temp/03.SNV/${tumorname}.vardict.vcf  --output_file  /temp/03.SNV/${tumorname}.vardict.vep.vcf  --dir_plugins /opt/vepPlugins  --cache  --dir_cache /opt/vep/cache --fork 20 --format vcf --vcf  --symbol --species ${species}  --terms SO --tsl  --biotype  --hgvs  --fasta /annotdir/${fasta} --offline   --plugin Frameshift --plugin Wildtype  --pick --force_overwrite
fi

###################################################################
#########     step6   somatic and AF filter           #############
###################################################################
if [[ "${runmodel}" == *"snvfilter"* ]];then
    echo "分析包含模块somatic and AF filter"
    setpriv --reuid=2001 --regid=1000   --clear-groups bash -c  "/usr/bin/java -Xmx10g -jar ${SnpSift} filter \"( na FILTER ) | (FILTER = 'PASS') && (QUAL >=30) && (AF >= 0.1) && ( DP >= 5 )\" --file ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.vep.vcf |grep -v Germline > ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.vep.somatic.vcf"
 fi
 
###################################################################
###########   step7  neoantigen prediction            #############
###################################################################
if [[ "${runmodel}" == *"neopred"* ]];then
    echo "分析包含模块 neoantigen prediction"

    echo "docker run  --rm -v ${realoutpath}/${tumorname}:/tmp/${tumorname} --user ${userid}:${groupid}  griffithlab/pvactools pvacseq run  --n-threads 25 \ --iedb-install-directory /opt/iedb --allele-specific-binding-thresholds --keep-tmp-files  --minimum-fold-change 1 \ /tmp/${tumorname}/03.SNV/${tumorname}.vardict.vep.somatic.vcf   ${tumorname} ${mhctype}  NetMHCIIpan NetMHCpan /tmp/${tumorname}/05.AntigenPrediction  -e1 8,9,10,11 -e2 15"

    docker run  --rm -v ${realoutpath}/${tumorname}:/tmp/${tumorname} --user ${userid}:${groupid}  griffithlab/pvactools pvacseq run  --n-threads 25 --iedb-install-directory /opt/iedb --allele-specific-binding-thresholds --keep-tmp-files  --minimum-fold-change 1 /tmp/${tumorname}/03.SNV/${tumorname}.vardict.vep.somatic.vcf   ${tumorname} ${mhctype} NetMHCpan NetMHCIIpan /tmp/${tumorname}/05.AntigenPrediction  -e1 8,9,10,11 -e2 15

    docker run --rm -v ${realoutpath}/${tumorname}/05.AntigenPrediction/combined:/tmp/combined  --user ${userid}:${groupid} griffithlab/pvactools pvacseq generate_aggregated_report /tmp/combined/${tumorname}.filtered.tsv /tmp/combined/${tumorname}.filtered.temp.tsv

    setpriv --reuid=2001 --regid=1000   --clear-groups bash -c "head -1 ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.filtered.temp.tsv > ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.filtered.aggregated.tsv"
    tail -n+2  ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.filtered.temp.tsv|cut -f4,5|while read line;
        do
            gene_variant=$(echo ${line} | tr ' ' '\t')
            echo ${gene_variant}
            setpriv --reuid=2001 --regid=1000   --clear-groups bash -c "grep \"${gene_variant}\" ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.all_epitopes.aggregated.tsv >> ${outpath}/${tumorname}/05.AntigenPrediction/combined/${tumorname}.filtered.aggregated.tsv"
        done
fi

###################################################################
###########   step7  Fusion calling                   #############
###################################################################
if [[ "${runmodel}" == *"fusioncall"* ]];then
    echo "分析包含模块 Fusion calling "
    
    #export PATH=/mnt/data2/wuzengding/03.biotools/software/STAR-2.7.2d/bin/Linux_x86_64:$PATH
    #/mnt/data2/wuzengding/03.biotools/software/STAR-Fusion.v1.9.0/STAR-Fusion  --genome_lib_dir /mnt/data2/wuzengding/03.biotools/software/STAR-Fusion.v1.9.0/Mouse_gencode_M24_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir --left_fq ${outpath}/${tumorname}/01.QC/${tumorname}_R1.fastq.gz --right_fq ${outpath}/${tumorname}/01.QC/${tumorname}_R2.fastq.gz --output_dir ${outpath}/${tumorname}/04.Fusion
    
    
    ## generate protein.fasta
    #/usr/bin/java -Xmx10g -jar /mnt/data2/wuzengding/03.biotools/software/snpEff/snpEff.jar -c /mnt/data2/wuzengding/03.biotools/software/snpEff/snpEff.config -noLog -v -noInteraction -nextProt mm39  -fastaProt ${outpath}/${tumorname}/${tumorname}.somatic.rot.fa ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.somatic.filter.vcf > ${outpath}/${tumorname}/03.SNV/${tumorname}.vardict.somatic.annot.vcf
    
    ## epitope pridiction
    #/mnt/data1/feicaiyi/netmhc/netpan41/netMHCpan -a H-2-Kb,H-2-Db -s -f  ${outpath}/${tumorname}/03.SNV/${tumorname}.somatic.prot.fa -xls -xlsfile ${outpath}/${tumorname}/05.AntigenPrediction/${tumorname}.NetMHCpan_out.xls -inptype 0 -BA
fi