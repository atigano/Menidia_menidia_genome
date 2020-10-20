# Menidia_menidia_genome
Scripts accompanying manuscript on sequence and structural variation in the Atlantic silverside (Menidia menidia)

## Genome assembly with 10X data

```
supernova run --id supernova_run_output_270mreads --fastqs /workdir/chromium/H7NF5BCX2_10388115_supernova_mkfastq_output_24Jan2018/outs/fastq_path/ --maxreads 270000000
supernova mkoutput --asmdir supernova_run_output_270mreads/outs/assembly --outprefix meme_270mreads_1_pseudohap --style pseudohap
```

## BUSCO
```
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/dovetail/menidia_menidia_tidy.filter1000.fasta -o busco_dovetail_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/dovetail/menidia_menidia_26Jul2018_hZUSG.fasta -o busco_dovetail_chicago_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/dovetail/menidia_menidia_30Jul2018_eLAsD.fasta -o busco_dovetail_hic_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
python /workdir/anna/busco/scripts/run_BUSCO.py -i /workdir/anna/chromium/meme270m500.fasta -o busco_10x_actino -l /workdir/anna/busco/actinopterygii_odb9 -m genome
```
