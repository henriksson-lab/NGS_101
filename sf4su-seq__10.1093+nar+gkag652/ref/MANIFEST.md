# Source manifest

```sh
python3 tools/get_sources.py --slugs sf4su-seq__10.1093+nar+gkag652
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/sf4su-seq__10.1093+nar+gkag652
```

Originals include PMC13309788 and its Supplementary Methods, the official NEB E7300
manual Version 9.0_07/25, the versioned bioRxiv preprint and Dölken 2008
(10.1261/rna.1136108). The independently read Dölken publisher original is also cached in
the TT-seq source bundle under `High_resolution_gene_expression_profiling_for_simultaneous_kinetic_parameter_analysis_of_RNA_synthesis_and_decay.pdf`.
The preprint v2 is dated 24 October 2025; its DOI year alone does not identify its version.

No original PDFs, extracted text or inspection images are committed. Source reading,
public wash-difference evidence and remaining actual linked targets are recorded in the
protocol review JSON.

The bounded deposited-library follow-up also registers GSE179306 sample-only full SOFT
(`targ=gsm`, excluding platform tables), original ENA XML for SRX11439177 and NCBI full
SRA experiment metadata (`efetch.fcgi`). These are small metadata originals; raw sequencing
files are not fetched. The version-1 original used for comparison is cached in the SNU
source bundle as `2021.07.14.452379v1.full.pdf`.
