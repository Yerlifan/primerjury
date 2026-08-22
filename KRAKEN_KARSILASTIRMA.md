# The Kraken2 rerun, four methods side by side

Every row is a bin. The bins are piles of reads separated
according to the original Kraken output (`reads_<taxid>.fastq`).

| # | Column | What it means |
|---|---|---|
| 1 | the original Kraken output | The bin's label. The reads were already separated by that claim. |
| 2 | PlusPFP esik=None | The same data at a high confidence threshold. Weak assignments drop out. |
| 3 | PlusPFP esik=0 | A database of broad coverage, with no threshold. |
| 4 | Alignment based identity | Our own measurement (`kimlik_sonuc.csv`). |

An empty cell: no identity reached %20 of the reads, so there is no decision.

| Bin (the source study's Kraken) | Threshold raised | PlusPFP | Alignment | Result |
|---|---|---|---|---|
| Alistipes finegoldii | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Bacteroides ovatus | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Bacteroides thetaiotaomicron | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Blastochloris tepida | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Ca. Cloacimonas acidaminovorans | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Ca. Methanomassiliicoccus intesti… | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Ca. Nitrosocosmicus hydrocola | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Ca. Sulfurimonas baltica | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Cloacibacillus porcorum | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Colletotrichum higginsianum | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Cyrtohymena siliat | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Helotiales askomikot | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanobrevibacter sp. AbM4 | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanocorpusculum labreanum | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanoculleus bourgensis | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanoculleus chikugoensis | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanoculleus receptaculi | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanofollis liminatans | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanosarcina barkeri | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanosarcina hadiensis | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanosarcina mazei | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanosarcina sp. WH1 | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanosarcina thermophila | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Methanothrix soehngenii | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Microascaceae askomikot | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Oceanobacillus sp. | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Petrimonas mucosa | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Petrimonas sulfuriphila | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Podospora pseudocomata | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Podospora pseudopauciseta | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Proteiniphilum propionicum | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Proteiniphilum saccharofermentans | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Pseudodifflugia amip | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Rhizoctonia solani | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Schizosaccharomyces osmophilus | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Schizosaccharomyces pombe | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Sphaerochaeta associata | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Toxoplasma gondii | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Trichoderma atroviride | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Trichoderma breve | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Tulosesus callinus | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Ustilaginoidea virens | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Vampyrellid amip | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |
| Zoopagomycota mantari | _karar yok_ | _karar yok_ | _yok_ | PlusPFP karar vermedi |


## Sonuc

There is NO PlusPFP run. The table is incomplete and cannot be interpreted.
To run it: `bash kraken_tool.sh esik-a`
