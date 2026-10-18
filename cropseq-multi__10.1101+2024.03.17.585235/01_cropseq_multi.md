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

## 5′ direct-capture branch (read 2026-10-17)

🟢 Methods, "Direct-capture single-cell RNA sequencing with 10X Genomics 5' capture" (v2 preprint lines 842–874, HeLa cells): the workflow follows 10x Genomics CG000735 (Chromium GEM-X Single Cell 5′ Reagent Kits v3 with Feature Barcode technology for CRISPR Screening) with two substitutions. The RT mastermix (RT Reagent E 2001106, Poly-dT RT Primer B 2001110, Reducing Agent B 2000087, RT Enzyme E 2001105/2001146; 27.5 µL) receives 0.5 µL each of 10 µM `oRW1220` and `oRW1221`, because the kit's CRISPR Poly-dT Primer Mix B (2001145) "does not have perfect complementarity to both sgRNAs". After GEM generation, RT and cDNA amplification, SPRIselect size-separates the cDNA; the large fraction goes to the standard gene-expression protocol (10x CG000733) and the small fraction carries the sgRNA products.

🟢 The sgRNA fraction is dialled out in one reaction with three primers — universal forward `oRW1223` plus `oRW1226` (sgRNA 1) and `oRW1227` (sgRNA 2) — in 1X NEBNext Ultra II Q5, 98 °C 1 min, 11 × (98 °C 10 s, 67 °C 10 s, 72 °C 15 s), 72 °C 1 min, then 1X SPRIselect and the commercial CG000735 sample-indexing PCR. Gene-expression and sgRNA libraries are pooled 10:1 and sequenced on an Illumina NovaSeq X with 28 cycles read 1, 10 cycles per index and 90 cycles read 2.

🟡 The graph now declares these named oligos, the SPRIselect split with its gene-expression partner library and the indexing PCR; the kit-internal cDNA amplification and gene-expression library construction remain warned "magic happens" boundaries because CG000735/CG000733 are not reproduced in the preprint.

🟢 Screen-library (not single-cell) readouts are separate: plasmid and genomic-DNA libraries are single-step PCRs with NEBNext Ultra II Q5 + EvaGreen (10 cycles for plasmid, 24 for gDNA), purified on AMPure XP or QIAquick, and sequenced on MiSeq — 300-cycle kit (MS-102-2002) when primers flank both scaffolds to score deletion recombination, otherwise a 150-cycle kit (MS-102-3001) with custom read primers and 75 cycles per read, 8 cycles per index (lines 735–754, 814–823).

🔴 For the 5′ direct-capture mode, the preprint Methods ("Direct-capture single-cell RNA sequencing with 10X Genomics 5’ capture") explicitly follow 10x guide CG000735 with custom primers and amplification conditions. That vendor guide is absent from the local source bundle and remains unreviewed, so the 10x-derived parts of this branch cannot yet be certified from the preprint alone.

## In situ mRNA and T7-IVT branches (read 2026-10-17)

🟢 The v2 preprint's Methods ('In situ amplification - mRNA protocol', lines 876–912; 'In situ amplification - T7 IVT protocol', lines 914–937; 'Sequencing by synthesis', lines 939–953) read both branches optically, not by PCR. The mRNA protocol runs RevertAid H minus RT with 1 µM of each biotinylated RT primer, then streptavidin binding and post-fixation in 3% formaldehyde and 0.1% glutaraldehyde. Padlock gap fill and ligation follow (0.1 µM of each padlock, RNase H, TaqIT, Ampligase), then 16 h Phi29 RCA, SBS-primer hybridization and in situ SBS with MiSeq Nano v2 reagents. The T7-IVT protocol first transcribes in situ with T7 RNA Polymerase Mix (NEB E2040L, usually 0.2X), post-fixes, and then follows the mRNA protocol from RT onwards. The workbook names the RT primers (oRW893 for iBAR 1; oRW887, oRW965 and oRW967 for iBAR 2 by tRNA), the padlocks (oRW963; oRW968, oRW964 and oRW966) and the SBS primers.

🟡 The graph now models these branches with the in situ steps. Earlier it wrongly sent them through an endpoint PCR. Padlock arm landing and the gap-filled span are not assembled, so gap fill, ligation and RCA form one unmodelled step.

## Local validation migration (2026-10-18)

Re-read the defining v2 preprint direct-capture Methods, lines 842–875, and the
v2 supplementary workbook `oligos` rows 29–33. Ran the bundle oligo scraper and
checked its five capture/dial-out primer hits against those rows. Also inspected
CG000735 Rev A introduction and final-library sequence section as a candidate vendor
revision; the paper does not pin that revision, so its N16/N12 bead boundaries are
not promoted into the defining model.

🟢 The direct-capture run is 28/10/10/90, with two 10-base sample indexes. The
old shared endpoint incorrectly imposed one 8-base index. Both direct-capture
terminals now declare the run; the optical RCA terminals retain their separate
in situ SBS readout. The five printed RT/dial-out primers are explicit graph inputs,
including the RT primers' terminal phosphorothioates. The model distinguishes the
Pol II identity transcript, mature guide RNA and CSMv2 DNA used for T7 transcription;
the alternative assay routes are an explicit preparation boundary.

🟡 Exact bead/capture, custom primer landing and kit index-PCR junction geometry
remain unassembled. The source selects a kit identifier without a revision; the
cached Rev A is a candidate, not proof of the study's revision. No sequence-absence
claim or completed full chemistry audit follows from this ownership migration.
