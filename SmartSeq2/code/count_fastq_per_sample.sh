root_dir="/gfs/work/kralbrecht/18_matthias/scrnaseq-full/data.dir"

outfile=count_fastq_per_sample.txt
rm ${outfile}

for folder in $(ls ${root_dir})
do
    # echo $folder
    n_fastq=$(ls ${root_dir}/${folder} | wc -l)
    echo -e "${folder}\t${n_fastq}" >> count_fastq_per_sample.txt
done
