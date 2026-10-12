# GenIE-ATAC — source notes

Defining paper: Cooper et al., *Nucleic Acids Research* (2023), DOI
[10.1093/nar/gkad332](https://doi.org/10.1093/nar/gkad332).

🟢 GenIE-ATAC loads Tn5 with MEDS-A alone, performs biotinylated locus-specific linear
PCR and streptavidin capture, then a semi-nested locus/MEDS-A PCR and indexed MiSeq PCR.
The final schematic must keep the locus-specific primer as a parameter rather than imply
one universal sequence.

## Exact transposome and representative target

🟢 MEDS-A is `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`: the canonical S5 primer-entry handle followed by the 19-base Tn5 mosaic end. It anneals only across that mosaic end to 5′-phosphorylated `CTGTCTCTTATACACATCT` MEDS-REV, whose source notation ends in `/intT/`. The handle remains single-stranded in the loaded transposome. Both ends transferred to an accessible fragment are MEDS-A; the protocol deliberately does not use the usual heterologous A/B pair.

🟢 Tn5 transfer leaves the usual 9-base staggered gaps. The schematic now keeps that newly transferred state apart from the filled, amplifiable duplex.

🟢 Supplementary Table 2 supplies a complete concrete example for `rs12269414`: `/5Biosg/CCAGGAGTTGGAACCTGC` for linear extension and `TCGGCATTCCTGCTGAACCGCTCTTCCGATCTCCTCTTGCCTTGGCTGTG` for the nested forward oligo. The first sequence enriches a locus-to-MEDS-A strand by biotin capture; the second contributes a sequencing tail plus a different, nested genomic binding site. The common Tn5A primer is `TCGTCGGCAGCGTCAG`.

🟢 The defining experiment sequenced 2 × 150 bases on MiSeq and merged overlapping reads for amplicons below 295 bp. 🟡 The paper says its second PCR added published MiSeq barcodes but does not print which index-primer pair was used. The terminal indexed shell is therefore inferred, while the source-exact internal boundaries are shown explicitly.
