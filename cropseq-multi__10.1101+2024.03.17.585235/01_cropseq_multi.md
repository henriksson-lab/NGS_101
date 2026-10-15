# CROPseq-multi

**Defining paper:** Walton et al., [bioRxiv doi:10.1101/2024.03.17.585235](https://doi.org/10.1101/2024.03.17.585235). The July 2025 v2 preprint and its supplementary oligo workbook were read for this source review.

## Molecular boundary

🟢 A reverse-oriented U6–tRNA–sgRNA array is inserted into the lentiviral 3′ LTR and duplicated during integration. RNases P and Z process the Pol III transcript into multiple guides without cleaving the lentiviral genomic RNA or the Pol II selection transcript.

🟢 Each guide carries a 12-nt internal barcode in the synthetic crRNA–tracrRNA loop, 19 bp from its linked spacer. The iBAR is transcribed in polyadenylated mRNA and in the mature guide; CSMv2 adds a T7 promoter for optical pooled screens.

## Evidence boundary

🟢 Vector sequences and design resources are published with the paper. Guide spacers and iBARs are library members rather than fixed adapters, so the schematic records lengths, linkage and readout routes as structured features.

## Oligo source audit

🟢 The author preprint's workbook `2024.03.17.585235_media-1.xlsx`, sheet `oligos`, supplies exact oligos for in-situ guide detection, 5′ direct-capture RNA readout, library cloning and NGS, including synthesis modifications. For the in-situ branch it lists biotinylated RT primers (for example `oRW893`, `/5Biosg/GTTCAGTGGTAGAATTCTCGCCTG`), phosphorylated padlock probes (for example `oRW963`, `/5Phos/gcattggccgggaatcgCTCCTGTTCGACAGTCAGCGCCATCTCCGACTTATTgaacttgctaagcagtc`) and sequencing-by-synthesis primers. The selected tRNA/guide branch determines which of these oligos is used.

🟢 The same sheet gives 5′ direct-capture primers `oRW1220` and `oRW1221` with terminal phosphorothioate bonds, and the custom Illumina sequencing primers `oRW1101`, `TCTTTCCCTACACGACGCTCTTCCGATCTGTTCGATTCCCGGCCAATGCA`, and `oRW1102`, `TTCAGACGTGTGCTCTTCCGATCTGCCTTATTTCAACTTGCTATGCTGTT`. It prints multiple indexed P5/P7 PCR primers separately. This establishes exact sequence for particular readout routes while the guide spacers and iBARs remain library members.

🟡 The current schematic is role-level across several demonstrated detection modes. An exact final construct/read layout needs a selected mode, vector version, guide/barcode member and indexing primer. The source does not justify combining the in-situ padlock path with the 5′ capture or direct NGS path into one molecule.

🔴 For the 5′ direct-capture mode, the preprint Methods ("Direct-capture single-cell RNA sequencing with 10X Genomics 5’ capture") explicitly follow 10x guide CG000735 with custom primers and amplification conditions. That vendor guide is absent from the local source bundle and remains unreviewed, so the 10x-derived parts of this branch cannot yet be certified from the preprint alone.
