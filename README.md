# Bioinformatics
The application of tools of computation and analysis to the capture and interpretation of biological data
1) To know what tools to use and how to use the tools
2) How to interpret data and test our hypotheses: statistic tests, null model, and visualization

## Tools
* [anaconda or miniconda](https://docs.conda.io/en/latest/miniconda.html): get python 3.9 version\
Many tools + their required environments can be easily installed by conda!\
Can you find how to install blast by conda? Hint: search `conda install blast` on Google.\
Homework 1: install conda and [install blast by conda](https://anaconda.org/bioconda/blast).
* [pip](https://anaconda.org/anaconda/pip)\
Many python packages can be easily installed pip!\
Homework 2: install pip and use pip to install pandas (`pip install pandas`)
* sequence similarity search: of course, to find sequences in your data that are similar to what have been annotated by other people\
    * Blast+: benchmark, slow but accurate
        * `makeblastdb -in your_database -dbtype nucl`: make db for DNA
        * `makeblastdb -in your_database -dbtype prot`: make db for protein
    * [Diamond](https://anaconda.org/bioconda/diamond): fast, protein database only
        * `diamond makedb --in your_database -d your_database.dmnd`
        * `diamond blastp --query your_sample.faa --db your_database.dmnd --out your_result.txt --id 50 --query-cover 80 --outfmt 6 --max-target-seqs 2 --evalue 1e-1 --threads 10`
    * [Usearch](https://www.drive5.com/usearch/download.html): fast, nucleotide database only?
        * `usearch -makeudb_usearch your_database -output your_database.udb`
    * [hmm](https://anaconda.org/bioconda/hmmer): domain based
        * `hmmsearch --tblout your_result.txt --cpu 10 -E 0.01 your_database.hmm your_sample.faa`
Homework 3: find HGT regions among 2 genomes by similarity >= 100% over 1kbp genomic region
* sequence mapping/alignment: bioinformaticians just enjoy making up words\
To align reads, usually metagenomes, to a known reference genome
    * abundance of a known species in a sample
    * SNPs, indels of a known species in a sample -> infer strains, genes under positive selection\
**Tools**
    * [Mapper](https://github.com/mathjeff/Mapper): fast, accurate, short read only, easy to use.
    * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml): slow, short read only, designed for human genomes, requires other tools such as samtools and bcftools for downstream analysis.
    * [minimap2](https://github.com/lh3/minimap2): slightly faster than bowtie2, long read too, designed for human genomes, requires other tools such as samtools and bcftools for downstream analysis.\
**Homework 4**: use mapper to compute the abundance of E.coli in a metagenome
* assembly: assemble short reads into long reads 
    * for WGS (whole genome sequencing): SPAdes
    * for metagenomes: metaSPAdes, IDBA-UD, MEGAHIT
    * for long reads: metaFlye, Raven, Canu (challenges: high indel errors)\
**binning: long reads -> draft genome**
    * MetaBAT2, CONCOCT, MaxBin, DasTOOL\
**dereplicate genomes**:
    * dRep (all genomes from all environments at 97-99% identity)
**annotate taxonomy**
    * gtdb classify_wf: also for WGS\
**check contamination**
    * checkM
* Metagenome analysis
    * taxonomy annotation
        * kraken
        * strainphylan
        * metaphylan
    * gene annotation
        * [Eggnog (COG)](http://eggnog5.embl.de/): hmm, free
        * KEGG: not free [BlastKOALA](https://www.kegg.jp/blastkoala/) blast, online only; [KofamKOALA](https://www.genome.jp/tools/kofamkoala/): hmm, online and local
        * [CAZY](http://www.cazy.org/): Carbohydrate-Active enZYmes Database, online and local
        * ARG: [arg_ranker](https://github.com/caozhichongchong/arg_ranker/tree/v2.0/arg_ranker): risk analysis, local; [CARD](https://card.mcmaster.ca/): online and local; [args_oap](https://smile.hku.hk/SARGs): online only;
        * Database of your interest
    * strain identification
        * InStrain, strainGE, strainfinder
    * MGE identification
        * [SRID](https://github.com/XiaofangJ/SRID)
        * [MGEfinder](https://github.com/bhattlab/MGEfinder)
    * homework 5: annotate ARGs in a metagenome sample
* genome analysis
    * taxonomy annotation
        * 16S: blast to 16S database, such as [greengenes](https://greengenes.secondgenome.com/),
        * gtdb classify_wf
    * gene annotation
        * prokka
        * same databases for metagenome analysis, different cutoff, why?
    * MGE identification
        * [Identify integrative and conjugative elements on genomes](https://db-mml.sjtu.edu.cn/ICEfinder/instruction.html)
        * [Differentiate chromosomes from plasmids](https://cge.cbs.dtu.dk/services/PlasmidFinder/)
        * [Identify integrons on genomes](https://github.com/caozhichongchong/I-VIP)
        * [MGEfinder](https://github.com/bhattlab/MGEfinder): requires fastq
    * SNPs and indels identification 
        * mapper
        * DeepVariant: bam files (output of bowtie2, minimap2, ...) as input
    * genome - genome comparison
        * fastANI: pairwise comparison
        * pangenome analysis: core, flexible, unique genes
            * [roary](https://sanger-pathogens.github.io/Roary/)
* 16S analysis
    * Qiime2 (install + execute in qiime2.sh.txt)

## Data interpretation
* 

## Contact
anniz44@mit.edu or caozhichongchong@gmail.com
