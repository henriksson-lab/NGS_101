# stLFR — research notes

Defining source: Wang *et al.*, [10.1101/gr.245126.118](https://doi.org/10.1101/gr.245126.118).

## Sources and selected route

🟢 The original final article (PMC6499310, Results/Figure 1 and Methods) and Cheng *et al.* [10.1038/protex.2018.116](https://doi.org/10.1038/protex.2018.116), original October 2018 manuscript, describe two routes. The page selects the single-Transposon1T, 3′-branch-ligation route used for stLFR-1–3. It does not combine this with the two-transposon route used for stLFR-4.

The protocols.io import explicitly warns that its internal references may differ from the original. The original manuscript is attachment 5, the original assembly figure attachment 6, and Supplement_0.xlsx attachment 7. The workbook is also available as bioRxiv supplement `324392_324392-2.xlsx`; these copies are identical.

🟢 Original manuscript steps 61–81 (pp. 9–10) give capture/barcode ligation, Exonuclease I + III adapter digestion, SDS treatment, UDG + APE1 excess-capture-oligo digestion, ET SSB binding, 3′-branch ligation, then PCR1/PCR2. The Exonuclease I treatment of the separate two-transposon route is not added to the UDG/APE1 stage here.

🟢 The original workbook provides the common-bead, RNA bridge, transposon, PCR, branch, BGI readout and alternative Illumina oligos, with their modifications. Each of its three barcode sets contains 1,536 ten-base members. The page displays their architecture; the source workbook remains the authoritative member list. Its alternative Illumina PCR2 contains literal `[i5]`, with no index width assigned here.

🟢 The cited [3′-branch-ligation study](https://doi.org/10.1093/dnares/dsy037), Methods §2.1/§2.3 and Figure 1 caption, describes joining a blunt donor's 5′ phosphate to a recessed, nicked or gapped acceptor's 3′ hydroxyl without requiring complementary ends. Its Supplementary Table S1 contains its own test substrates and adapters. Those experimental oligos are not assigned to stLFR merely because the mechanism is shared.

## Published PCR-programme difference (round 2026-10-17-breadth-first-19)

🟢 Final article Methods, "stLFR with 3′ branched adapter ligation" (PMC6499310 JATS line 212), gives Exonuclease I + III in TA buffer 10 min 37°C, 1% SDS stop/Tn5 removal, UDG + APE1 30 min 37°C, 300 pmol second adapter with 4000 U T4 ligase 2 h, and states that the PCR mix and conditions were the same as the two-transposon process (line 208: PfuCx, 400 pmol each primer, 10 min 72°C then 12 cycles). The article Methods does not mention the SSB binding step that the original manuscript step 74–75 includes; this is treated as an omission in the summary, not a conflict.

🟢 Original manuscript step 80 (p. 9) instead gives 95°C 3 min then 10–12 cycles of 95°C 10 s / 58°C 30 s / 72°C 2 min for the branch route; the 72°C 10 min first step appears only in two-transposon step 58. The page shows this published difference; the model selects neither programme nor a cycle count.

## Published common-primer variants

🟢 Original Supplement_0.pdf Figure 1 prints bead-common top `/52-Bio/AAAAAAAACCCGTAGCCATGTCGTTCTG`: eight A bases followed by `CCCGTAGCCATGTCGTTCTG`. Original Supplement_0.xlsx instead prints `/52-Bio/AAAAAAAAAATGTGAGCCAAGGAGTTG`: ten A bases followed by `TGTGAGCCAAGGAGTTG`, with matching PCR1 `TGTGAGCCAAGGAGTTG`.

The model uses the ordered workbook and its matching PCR1. The original manuscript pp. 5–6 (text lines 194–195) explicitly allows changes to the common-adapter PCR sequence for different sequencing technologies. Thus these are unresolved source variants, not a demonstrated publication error. The sources inspected do not establish which physical preparations used each common region. The difference is also shown on the public schematic.

## Barcode assignment and sequencing

🟢 Barcode assignment is statistical. The final article reports multiple long molecules per barcode and condition-dependent uniqueness; the original protocol likewise describes most beads capturing fewer than ten long molecules and gives its uniqueness estimate for a specified DNA/bead condition. Capture does not physically reject molecules merely because their barcode will be shared.

🟢 The original manuscript's demultiplexing section (text 356–368) specifies paired-end 100-base genomic sequencing and a 42-position barcode/common-adapter string associated with Read 1. The three ten-base codes are separated by two six-base linkers. The workbook separately labels the two-transposon and branch-route sequencing primers.

🟢 The [current author demultiplexer repository](https://github.com/stLFR/stLFR_read_demux), master README retrieved 2026-10-16, describes paired-end inputs with barcode appended to Read 2, plus 42- and 54-position variants; its examples use paired-end 54 and single-end 42RC scripts. This is a versioned processing declaration, not proof that the legacy protocol's physical readout should be relabelled. The page preserves both scopes.

## Model limits and remaining source review

🟡 The graph follows the selected published reactions and renders all declared waste products. Its reported adapter ends do not establish a complete base-resolved capture/branch junction. The physical PCR product and single-stranded DNB circle are not yet constructed; circularization/DNB formation is an explicit warning-producing step.

The final article cites [Drmanac 2010](https://doi.org/10.1126/science.1181498) for circularization. The standard fetch route did not retrieve that original. This is an actual outstanding precursor source, not evidence that circle information was unpublished. The source review also retains the original supplements, platform citations and relevant linked scope for further checking. No missing-information claim has been accepted.
