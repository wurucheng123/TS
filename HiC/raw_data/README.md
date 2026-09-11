juicer_tools  https://github.com/aidenlab/Juicebox/releases/download/v2.20.00/juicer_tools.2.20.00.jar
3DMax.jar     https://github.com/BDM-Lab/3DMax/releases/tag/v1.0

```bash
hic2cool convert /data/data/Turner/HiCData/ZWJFamilyHiCAfterMapping/outPut/step04_hicpro2juicebox/hicfiles/Z-naive.allValidPairs.hic Z-naive.allValidPairs.cool -p 4
hic2cool convert /data/data/Turner/HiCData/ZWJFamilyHiCAfterMapping/outPut/step04_hicpro2juicebox/hicfiles/F-naive.allValidPairs.hic F-naive.allValidPairs.cool -p 4
hic2cool convert /data/data/Turner/HiCData/ZWJFamilyHiCAfterMapping/outPut/step04_hicpro2juicebox/hicfiles/M-naive.allValidPairs.hic M-naive.allValidPairs.cool -p 4
hic2cool convert /data/data/Turner/HiCData/ZWJFamilyHiCAfterMapping/outPut/step04_hicpro2juicebox/hicfiles/Z-IgD-.allValidPairs.hic Z-IgD-.allValidPairs.cool -p 4
hic2cool convert /data/data/Turner/HiCData/ZWJFamilyHiCAfterMapping/outPut/step04_hicpro2juicebox/hicfiles/F-IgD-.allValidPairs.hic F-IgD-.allValidPairs.cool -p 4
hic2cool convert /data/data/Turner/HiCData/ZWJFamilyHiCAfterMapping/outPut/step04_hicpro2juicebox/hicfiles/M-IgD-.allValidPairs.hic M-IgD-.allValidPairs.cool -p 4
```

```bash
wget -c https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
# gtf2bed4trackc -g Homo_sapiens.GRCh37.75.gtf -o Homo_sapiens.GRCh37.75.bed13 --biotype2bed13
trackc gtf2bed Homo_sapiens.GRCh37.75.gtf -o Homo_sapiens.GRCh37.75.bed13 --biotype2bed13
```