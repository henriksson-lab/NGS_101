# inDrops-2

**Defining source:** doi:10.1093/nar/gkae1312.

## Molecular boundary

🟢 Cells and barcoded hydrogel beads are co-encapsulated; captured 3′ cDNA can follow either IVT linear amplification or template-switching PCR amplification.

## Evidence boundary

🟢 The two amplification implementations remain explicit branches of one named protocol rather than being blended into one reaction.

## Oligo source audit

🟢 The published article's Methods (`nar_gkae1312_PMC11724362.xml.txt`, “Synthesis of V1 beads” and “Synthesis of V2 beads”) prints the full hydrogel-bead primer architecture for both versions, including the 5′ acrydite and photocleavable spacer. The corrected supplementary workbook `nar_gkae1312_supp_gkae1312_supplemental_files_gkae1312_Supplementary_Tables_Corrected.xlsx`, Table S1, gives the other reaction oligos with their modifications; Tables S2–S5 give the barcode primers and whitelists.

🟢 For the **IVT** branch, Table S1 identifies the v1 bead primer, `PE2-N6` second-strand RT primer (`GTG ACT GGA GTT CAG ACG TGT GCT CTT CCG ATC TNN NNN N` as spaced in the table), `PE1` PCR primer and six `PE2-Ind` indexing primers. The article describes a T7-promoter-containing v1 bead stub as `5′-/5Acryd/PC/CGATTGATCAACG TAATACGACTCACTATAGGGATACCATCTACACTCTTTCCCTACACGACGCTCTTCCG-3′`.

🟢 For the **template-switching** branch, Table S1 specifies the v2 bead stub `/5Acryd/CGATGACG/PC/CTACACGACGCTCTTCCGATCT`, template-switching oligo `AAGCAGTGGTATCAACGCAGAGTACATrGrGrG`, cDNA forward/reverse primers `CTACACGACGCTCTTCCGATCT` and `AAGCAGTGGTATCAACGCAGAGT`, and adapter-ligation oligos `/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC/3ddC` and `/5AmMC6/GCTCTTCCGATCT`. The table also gives the indexed PCR primers and marks the terminal `*` as a phosphorothioate linkage.

🟡 The barcode windows shown as numbered positions in Table S1 vary by bead. Their exact members are in Tables S2–S5; one arbitrary member must not be substituted for the family. The source therefore resolves the invariant structural and amplification oligos, while exact read placement requires choosing either the v1 IVT or v2 template-switching branch and a barcode/index instance.
