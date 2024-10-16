# A script to fill buffer characters denoted by "N" in a DNA sequence template while avoiding restriction enzymes and maintaining GC content

### Installation
To setup environment on Elzar (CSHL high performance computer):
```
module load EBModules BLAST+
conda activate biopython
```

To setup environment on local device, need to first install [blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), then setup conda environment:
```
conda create -n biopython
conda activate biopython
conda install -c conda-forge biopython
```

### Obtaining necessary genome reference files to avoid the Human and Mouse genomes
To get the Human genome reference:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

To get the [mouse genome reference](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/), download and find the fasta file.

The reference sequeces will need to be assembled in `fill_buffer_seq.py` for the first time by uncommenting the line `# subprocess.run(db_cmd, shell=True)`. This line can then be commented out for future runs once the ref db is assembled, which will save significant runtime.

### Running the main script
To run main script, add a single command line argument with a sequence containing any distribution of N characters to be filled:
```
python fill_buffer_seq.py GCAATGNNNNNNNNNNNNNNCTATACGTTGNNNNNNNNNNNNNCGTAGC
```


### (Optional) Running the main script from .dna file type
To get the input sequence from .dna file type, convert it to .fa file in SnapGene, then use `sequence=$(cat filename.fa | grep -v "^>" | tr -d '\n')`, and then can run the main script with `python fill_buffer_seq.py $sequence`.
