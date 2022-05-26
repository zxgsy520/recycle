# recycle
Obtain the mitochondrial or chlorophyll sequence from the reads of genome sequencing, and perform assembly.

## Requirements
* [Python](https://www.python.org/)
* [Unicycler](https://github.com/zxgsy520/Unicycler) Modified part of the code on the basis of the original version to make it more suitable for the assembly of bacteriophages, mitochondria and chloroplasts
* [gfa1](https://github.com/lh3/gfa1) Processing and transformation of gfa and fastg files
* [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) Software for assembling mitochondria and chloroplasts (some species have better results)
* [Bandage](https://github.com/rrwick/Bandage) Some samples need to be manually solved with Bandage
## Installation
```
git clone https://github.com/zxgsy520/recycle.git
cd  biodb/bin    
chmod 755 *    #Linux system can be used directly
cd ../scripts
chmod 755 *   #Supports all systems but requires python and related packages
```
or
```
wget -c https://github.com/zxgsy520/recycle/archive/refs/heads/main.zip
unzip main.zip

```
## Options and usage
Extract mitochondrial chloroplast sequence
```
minimap2 -t 10 -x sr mitochondrion.genomic.fasta clean.r1.fq.gz clean.r2.fq.gz >ngs.paf
choose_ml_map --input ngs.paf -p mgi --gc 20 --base all --score 0 \
  --read1 clean.r1.fq.gz --read2 clean.r2.fq.gz --prefix mitochondrion_read    #如果是三代数据，那不需要提供read2，-p 参数设置为illumina，paf文件输入三代比对后的结果即可。
```
Assemble mitochondria from metagenomes
```
megahit -1 clean.r1.fq.gz -2 clean.r2.fq.gz  --num-cpu-threads 20 --out-dir megahit
megahit_core contig2fastg 119 megahit/intermediate_contigs/k119.contigs.fa > megahit.k119.fastg
fastg2gfa megahit.k119.fastg > megahit.k119.gfa
get_organelle_from_assembly.py -F animal_mt -g megahit.k119.gfa -t 10 -o getorganelle
```

![mit_filter_graph](https://github.com/zxgsy520/recycle/blob/main/docs/mit_filter_graph.png)
