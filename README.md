To setup environment on Elzar:
```
module load EBModules BLAST+
conda activate biopython
```

To get the Human genome reference:
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

To get the [mouse genome reference](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/).

The reference sequeces will need to be assembled in `fill_buffer_seq.py` for the first time by uncommenting the line `# subprocess.run(db_cmd, shell=True)`. This line can then be commented out for future runs once the ref db is assembled, which will save significant runtime.

To run main script, add command line argument with a sequence with any distribution of N characters to be filled:
```
python fill_buffer_seq.py GCAATGNNNNNNNNNNNNNNCTATACGTTGNNNNNNNNNNNNNCGTAGC
```