# recycle
Obtain the mitochondrial or chlorophyll sequence from the reads of genome sequencing, and perform assembly.

## Requirements
* [Python](https://www.python.org/)
* [Unicycler](https://github.com/zxgsy520/Unicycler) Modified part of the code on the basis of the original version to make it more suitable for the assembly of bacteriophages, mitochondria and chloroplasts
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
  --read1 clean.r1.fq.gz --read2 clean.r2.fq.gz --prefix mitochondrion_read
```


![graph](https://user-images.githubusercontent.com/36355222/147458005-47ff1a65-9cca-42cf-b321-6a270b06ac95.png)
