# inDrops-2

**Defining source:** doi:10.1093/nar/gkae1312.

## Sources read (round 2026-10-17-breadth-first-11)

Published article PMC11724362 (JATS `nar_gkae1312_PMC11724362.xml.txt`): Methods lines 105 (encapsulation, RT-lysis mixes, photo-release, RT), 113 (inDrops-2 IVT library), 117 (inDrops-2 TS library), 121–125 (V1/V2 bead synthesis and full-length primers). Corrected supplementary workbook, Table S1 (lines 16–55) and the heads of Tables S2–S5. Supplementary Protocols 1–4 were not re-read in this capped pass. The bioRxiv preprint is superseded and was not used.

## Molecular boundary

🟢 The two inDrops-2 implementations differ from encapsulation onward, so the model splits them before the droplet step (Methods line 105):

- **IVT:** v1 beads (T7 promoter, PE1 Read 1 site, two-part cell barcode, UMI, oligo(dT)); RT-lysis mix with Maxima H minus; UV photo-release; RT in droplets. Then second-strand synthesis (NEBNext Ultra II non-directional module, 2.5 h 16 °C), HiScribe T7 IVT (15 h 37 °C), RNA fragmentation (Ambion, 70 °C 1.5 min), random-primed RT with `PE2-N6` (Maxima H minus), and 10 cycles of KAPA HiFi PCR with `PE1` and `PE2-Ind` primers, 0.6–0.8× dual SPRI (line 113).
- **TS:** v2 beads; the RT-lysis mix itself contains 50 µM TSO, so template switching happens in the droplet RT. Then 14 cycles of KAPA cDNA PCR with cDNA FWD/REV primers, NEB FS DNA Library Prep fragmentation (8 min 37 °C) with end repair/dA-tailing (30 min 65 °C), ligation of the annealed modified adapter (55 nM, 15 min 20 °C), and 14 cycles of NEBNext Ultra II Q5 PCR with the Forward PCR index primer and `PE2-ind` primers (line 117).

🟡 Model correction: the earlier graph applied Tn5 tagmentation to both branches and placed template switching after a shared RT. Neither matches the article; both were replaced.

## Read layout

🟢 IVT: R1 54 cycles, i7 8, R2 ≥35 (NextSeq550/HiSeq2500). TS: R1 28 cycles, i7 8, R2 35–92 (line 117). The Table S1 PE2 indexes are 6 nt, so the 8-cycle i7 read runs 2 cycles past the printed index (🟡 derived).

## Published discrepancies

🟢 **v1 bead primer.** The Methods (line 121) give the full-length v1 primer as `/5Acryd/PC/CGATGACG TAATACGACTCACTATAGGG ATACCACCATGG CTCTTTCCCTACACGACGCTCTTCCGATCT[11 positions] GAGTGATTGCTTGTGACGCCTT[8 positions] NNNNNNNN T19`. Table S1 prints `CGATTGATCAACG [T7] ATACCACCATGG ...CCGATCT 1234567890 GAGTGATTGCTTGTGACGCCTT 12345678 NNNNNN TTTTTTTTTTTTTTTTTTTV`. The Table S2 bc2 primers (e.g. `bc2.1`) carry `AAAAAAAAAAAAAAAAAAA NNNNNNNN [8]`, which matches N8/T19. The Table S2 bc1 windows vary in length (8–11 nt). The model draws role-level barcode/UMI segments and does not choose between these values.

🟢 **v2 bead stub.** The Methods (line 123) give `/5Acryd/CGATGACG(PC)CTACACGACGCTCTTC-3′`. Table S1 lists a "V4 DNA stub" `/5Acryd/CGATGACG/PC/CTACACGACGCTCTTCCGATCT`. Table S4 `P_Bcd1` oligos begin `/5Phos/CGATCT…`. Ligating them to the Methods stub would rebuild `…CTTCCGATCT` (🟡 derived). The full-length v2 primer ends in T18 in the Methods (line 125) and T19 in Table S1. The page shows both discrepancies.

## Evidence boundary

🟢 The two amplification implementations remain explicit branches of one named protocol rather than being blended into one reaction.

## Oligo source audit

🟢 The published article's Methods (`nar_gkae1312_PMC11724362.xml.txt`, “Synthesis of V1 beads” and “Synthesis of V2 beads”) prints the full hydrogel-bead primer architecture for both versions, including the 5′ acrydite and photocleavable spacer. The corrected supplementary workbook `nar_gkae1312_supp_gkae1312_supplemental_files_gkae1312_Supplementary_Tables_Corrected.xlsx`, Table S1, gives the other reaction oligos with their modifications; Tables S2–S5 give the barcode primers and whitelists.

🟢 For the **IVT** branch, Table S1 identifies the v1 bead primer, `PE2-N6` "2nd step RT primer" (`GTG ACT GGA GTT CAG ACG TGT GCT CTT CCG ATC TNN NNN N` as spaced in the table), `PE1` PCR primer and six `PE2-Ind` indexing primers. The article describes a T7-promoter-containing v1 bead stub as `5′-/5Acryd/PC/CGATTGATCAACG TAATACGACTCACTATAGGGATACCATCTACACTCTTTCCCTACACGACGCTCTTCCG-3′`.

🟢 For the **template-switching** branch, Table S1 specifies a stub labelled "V4 DNA stub" `/5Acryd/CGATGACG/PC/CTACACGACGCTCTTCCGATCT`, template-switching oligo `AAGCAGTGGTATCAACGCAGAGTACATrGrGrG`, cDNA forward/reverse primers `CTACACGACGCTCTTCCGATCT` and `AAGCAGTGGTATCAACGCAGAGT`, and adapter-ligation oligos `/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC/3ddC` and `/5AmMC6/GCTCTTCCGATCT`. The table also gives the indexed PCR primers and marks the terminal `*` as a phosphorothioate linkage.

🟡 The barcode windows shown as numbered positions in Table S1 vary by bead. Their exact members are in Tables S2–S5; one arbitrary member must not be substituted for the family. The source therefore resolves the invariant structural and amplification oligos, while exact read placement requires choosing either the v1 IVT or v2 template-switching branch and a barcode/index instance.
