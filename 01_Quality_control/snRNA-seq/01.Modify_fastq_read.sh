cat ./00.raw_fastq_path_B1_L2.txt | while read line
do
	library=`basename ${line}`
	sample=`basename ${library} | cut -d "." -f 1`
	if [[ "$library" =~ "_1.fq.gz" ]]; then
    # for read 1
    echo "#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=24:mem=20G
#PBS -l walltime=24:00:00
#PBS -N ${sample}_modify_read1
#PBS -o ${sample}_modify_read1.o
#PBS -e ${sample}_modify_read1.e
zcat ${line} | awk '/^@/ && (NR % 4 == 1) {sub(/\/1/, \" 1\", \$0)} 1' | gzip > ./$sample.modified.fastq.gz">./${sample}.modify_read.sh
elif [[ "$library" =~ "_2.fq.gz" ]]; then
    # for read 2
        echo "#!/bin/bash
#PBS -q normal
#PBS -l select=1:ncpus=24:mem=20G
#PBS -l walltime=24:00:00
#PBS -N ${sample}_modify_read2
#PBS -o ${sample}_modify_read2.o
#PBS -e ${sample}_modify_read2.e
zcat ${line} | awk '/^@/ && (NR % 4 == 1) {sub(/\/2/, \" 2\", \$0)} 1' | gzip > ./$sample.modified.fastq.gz">./${sample}.modify_read.sh
else
    echo "File does not match expected pattern (_1.fq.gz or _2.fq.gz)"
fi
done
