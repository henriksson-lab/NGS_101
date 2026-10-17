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

🟢 The nested forward oligos of Supplementary Tables S1 and S2 all carry the same 5′ tail
`TCGGCATTCCTGCTGAACCGCTCTTCCGATCT`. That is the Illumina paired-end / Genome-Analyzer-era
adaptor arm belonging to the cited 96-plex barcoding method (ref 31, Kozarewa & Turner
2011); it is *not* the current TruSeq Read 2 arm `{{= illumina.TRUSEQ_READ2 }}`. The
reverse oligos, which belong to the matched gDNA-amplicon branch rather than the ATAC
amplicon, instead carry the TruSeq Read 1 tail `{{= illumina.TRUSEQ_READ1 }}`. The ATAC
amplicon's other end carries no primer-added tail at all: it is the tagmentation-installed
MEDS-A handle, amplified with the Tn5A primer.

🟡 The schematic therefore draws the source-exact tail and MEDS-A end and marks only the
outer flow-cell/index shell inferred. Which amplicon end receives P5 is not printed either,
so the drawn orientation (PE tail on the P5 side, MEDS-A on the P7 side) is a model choice.

🟢 The defining experiment sequenced 2 × 150 bases on MiSeq (Methods, "Sequencing and read
alignment") and merged overlapping reads with FLASH for amplicons below 295 bp. 🔴 Neither
the barcoding primers nor index-read lengths are printed, so index cycles are not shown as
acquired reads.

## Reactions read in the 2023 Methods

🟢 Nuclei from 0.5 × 10⁶ hiPSCs are tagmented 30 min at 37 °C, cleaned on MinElute, then
linearly amplified from the biotinylated locus primer with KAPA HotStart (95 °C 5 min; 50
cycles of 98/60/72 °C), spiked with non-hotstart Kapa for 50 further cycles. Biotinylated
products are captured on Dynabeads MyOne Streptavidin C1 (overnight in 2× BW buffer,
washed), and the beads are used directly as template for the nested PCR with PowerUp SYBR
Green (95 °C 10 min; 35 cycles of 95/57/72 °C), giving a 400–800 bp smear. Barcodes are
added in a second PCR from 1 µl of that product.
