# M3-seq — research notes

Defining source: McNulty et al., 2023, [doi:10.1038/s41564-023-01462-3](https://doi.org/10.1038/s41564-023-01462-3).

## Sources read

- 🟢 PMC full text and Methods.
- 🟢 Supplementary Figure 1 workflow and Supplementary Table 3 oligo workbook.
- 🟢 The authors' `brwaang55/m3seq_scripts` repository linked by the paper, including the
  demultiplexing read structure.

## Published oligos recovered

🟢 Supplementary Table 3 prints 288 phosphorylated round-one RT primers. Each has the
structure `TruSeq Read 1 (33 nt) – N8 UMI – 13-nt BC1 – N6 random primer`; Barcode 1 uses
BC1 `AAGTGATTAGCAA`.

🟢 The same table prints the full bridge oDS025
`CGTCGTGTAGGGAAAGAGTGTGACGCTGCCGACGA/3ddC/`. Its sequence is the reverse complement
of the 10x s5 arm followed by the first 21 bases of the RT-primer handle, fixing the
splinted nick geometry rather than merely naming a ligation.

🟢 It also prints the S³ second-strand randomer oBW140, the S³-handle PCR primer oDS028
`AAGCAGTGGTATCAACGCAGAGT`, P5 primers oDS027/oBW170 and oBW171, Read-2-only Tn5 oligos oDS029/oDS030, the T7–s7 primer oBW168, and the complete
indexed P7–s7 primer series oDS033–oDS040 and oBW176–oBW183. The page uses oBW176 as a
representative final i7 primer.

## Molecular path

🟢 Round-one random RT installs a 13-base BC1 and 8-base UMI in fixed bacteria. Round-two
10x ATAC droplet ligation installs a 16-base BC2. BC1+BC2 is the cell identity.

🟢 Methods (PMC10522482 JATS, sections “Round-one indexing” to “Final library prep”, xml.txt
lines 154–190), in order: Maxima H Minus in-situ RT with 7.5% PEG; pooling; Chromium scATAC
v1.1 Chip H loading with Ampligase and oDS025; 12 cycles of 98 °C/59 °C droplet ligation;
emulsion break, Lysis-T/proteinase K 55 °C 1 h, Dynabeads MyOne Silane cleanup; RNase H strip
of the template RNA and 1.8× SPRI; Klenow exo− second strand with oBW140 and 1.8× SPRI; Q5
PCR with oDS028 + oBW170 (14 cycles) and two 1.2× SPRI (the “GEX library”); 10 ng tagmented
with home-assembled Read 2-only EZ-Tn5 (oDS029:oDS030), SDS stop, 1.2× SPRI; Q5 PCR with
oBW170 + oBW168 (72 °C 3 min, 9 cycles), 1.2× SPRI; T7 IVT of 100 ng, DNase I, 2× SPRI; rRNA
probe hybridization, thermostable RNase H + RNase H 50 °C 45 min, DNase I removal of probes,
2× SPRI; final RT primed from the P5 end with oBW170 + oBW171 (so only full-length constructs
are copied), 1.2× SPRI; 9-cycle Q5 PCR with oBW170 + a P7 index primer; two 0.8× SPRI; qPCR
QC with oBW170/oBW176 and oBW141/oBW176.

🟡 Fragments carrying Read 2 adapters at both ends, and RNase-H-cleaved rRNA pieces, are not
shown as physically removed: they lack the ends required by the next primers and are left as
unresolved terminals on the page.

🟢 The paper specifies 26-cycle Read 1, 30-cycle i5, 8-cycle i7 and 74-cycle Read 2.
The authors' demultiplexing scripts further identify the 13-base BC1, 8-base UMI and
16-base BC2 windows.

🟢 The commercial 10x ATAC barcode candidates are distributed in the vendor's
[Cell Ranger ATAC software](https://github.com/10XGenomics/cellranger-atac)
(`lib/python/barcodes/737K-cratac-v1.txt`), and the [10x CG000209 v1.1 guide](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000209_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v1.1_UserGuide_RevG.pdf),
Rev G, shows the 16-base barcode between the P5 and Read 1N arms (guide pp. 63–64,
Appendix A). 🔴 Which candidate a particular M3-seq cell received is determined from
that experiment's reads; the paper's oligo table does not assign a fixed BC2 to a cell.
The page therefore shows BC2 as a variable role token while retaining the source-known
adapter arms and bridge geometry. The author repository lists a
`metadata/737K-cratac-v1.reverse_complement.csv` file; it is now fetched
(`gh_737K-cratac-v1.reverse_complement.csv`) but was not inspected in the 2026-10-17 pass, so
its members and orientation remain unchecked.

## Bounded pass 2026-10-17 (round breadth-first-12)

Model correction (repository error, not a publication error): oBW168 is printed as
`TAATACGACTCACTATAGGGTCTCGTGGGCTCGGAG` (Table S3-Standard oligos). The earlier model
wrote the T7 promoter plus `GGG` before s7 `GTCTCGTGGGCTCGG`, one G too many; it now uses
`TAATACGACTCACTATAGG` + s7 + `AG`, and its annealing 3′ end covers s7 plus the first two
mosaic-end bases. The graph now also includes the RNase H strip, silane cleanup, the
oDS028/oBW170 GEX PCR, DNase steps, the SPRI cleanups and the P5-primed final RT, and gives
the printed oligos exact states so the local source guard checks them.

The RT-oligo sheet is named “Table S3-RT oligos” but its title row reads “Table S4”; this is a
label inconsistency only. Not inspected in this pass: the 10x CG000209 Rev A guide cited in
the Methods (the notes above refer to Rev G), the probe-design reference 19 and the SPLiT/
scifi-RNA-seq predecessor methods (refs 9, 16, 42).
