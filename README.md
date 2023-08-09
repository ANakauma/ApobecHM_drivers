# Driver APOBEC-derived hotspot mutations in whole genomes 

This repository contains R scripts to identify driver hotspot mutations in whole-genomes of APOBEC-enriched tumors. The original script was developed for whole-genome sequencing data as it assumes that most hotspot mutations (>1 mutations in the same genomic position) are passengers, which is the base for the background mutational rate. The background mutational rate also considers the trinucleotide context, the stability of hairpin loops, the sequence within the loops and the genomic region. These variables are provided or can be estimated with the scripts of this repository.

The scripts could easily be adapted to include other relevant co-variates or for specific cases. 


## Running the code
The user must provide a table with mutations indicating the sample, chr, position, etc. (see example file `data/SNV_data.txt`). This repository contains the following scripts:

- `1.prepareData.R` identifies APOBEC positive samples and calculates the APOBEC fold enrichment score in tri- and tetra-nucleotide context.
- `2.Hairpinloops_stability.R` will predict hairpin loops around all provided genomic positions (SNVs) and calculates the thermodynamic stability. It will also identify twin mutations and didymi.
- `3. driverHotspotMuts.R` identifies driver hotspot mutations modeled as a poisson process.


## Citation

If you use this script in your work, please cite the original study:

Nakauma-González JA, Rijnders M, Noordsij MTW, Martens JWM, van der Veldt AAM, Lolkema MPJ, Boormans JL, van de Werken HJG. Whole-genome mapping of APOBEC mutagenesis in metastatic urothelial carcinoma identifies driver hotspot mutations and a  novel mutational signature. *Under Review*
