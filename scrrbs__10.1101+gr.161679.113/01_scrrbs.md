# scRRBS — single-cell reduced-representation bisulfite sequencing

Guo H, Zhu P, Wu X, Li X, Wen L, Tang F. *Single-cell methylome landscapes of mouse embryonic stem cells and early embryos analyzed using reduced representation bisulfite sequencing.* Genome Research 23:2126–2135 (2013). [doi:10.1101/gr.161679.113](https://doi.org/10.1101/gr.161679.113), PMC3847781; GEO GSE47343.

## Sources and active scope

🟢 The source-owned graph follows the published defining article, `gr.161679.113_PMC3847781.xml.txt`, Methods “Construction of single-cell RRBS sequencing libraries” lines 103–105 and data processing lines 111–113. The separate sperm MspI/HpaII nested-PCR assay at 107–109 is a validation experiment; its Supplemental Table 5 primers are not assigned to library construction.

The later detailed protocol [Guo et al., Nature Protocols 2015](https://doi.org/10.1038/nprot.2015.039), the cited Gu 2011 bulk predecessor and full publisher supplement remain deferred originals. Their unavailable chemistry is not filled from the upstream scg schematic or a modern vendor adapter design. Prior source-reading records remain in the source-review ledger; this bounded migration does not change the missing-information audit status.

Evidence: 🟢 source statement; 🟡 molecular interpretation or derived relation. An unassigned native role records the model boundary, not confirmed publication absence.

## Source construction route

1. 🟢 Pick one mouse cell into 5 µl lysis buffer (20 mM Tris-EDTA pH 8.0, 20 mM KCl, 0.3% Triton X-100). Add 1 mg/ml Qiagen protease and 60 fg unmethylated lambda DNA. Lyse for 3 h at 50 °C; heat-inactivate for 30 min at 75 °C. Mouse DNA and lambda spike-in are independent polymer populations.
2. 🟢 Digest with 9 U MspI in 18 µl for 3 h at 37 °C; source cut notation `C^CGG`. 🟡 The resulting two-base 5′-CG overhang is a local interpretation of that cut, not an adapter assignment. Both mouse and lambda substrates enter digestion.
3. 🟢 Fill and add an extra A to the 3′ blunted ends with 5 U Klenow exo-minus, 1 mM dATP, 0.1 mM dGTP and 0.1 mM dCTP in 20 µl. 🟡 Newly incorporated ordinary dCTP is distinct from any methylated genomic C; it is not a native methylation measurement.
4. 🟢 Ligate “Illumina standard premethylated indexed adaptors” with 30 U concentrated T4 DNA ligase in 25 µl. Native strands, exact protected positions, junction geometry, index member/width and adapter concentration remain unassigned in this bounded source scope.
5. 🟢 Bisulfite conversion: MethyCode kit as printed (Invitrogen), 150 µl, 98 °C for 10 min then 64 °C for 2.5 h. 🟡 Unmethylated C becomes U when converted, while native genomic 5mC and protected adaptor C retain cytosine identity. Former genomic partner strands therefore need not remain complementary. The two converted originals are represented independently, without reconstructing an unverified adapter duplex.
6. 🟢 On-column desulfonation and Zymo-Spin purification with 10 ng tRNA carrier; elute in 30 µl. Carrier RNA is an independent molecular class. The column yields retained eluate and physically removed flowthrough/washes; class-specific DNA/RNA recovery is unresolved.
7. 🟢 First PCR: 50 µl with 1 U uracil-tolerant PfuTurbo Cx, 95 °C for 2 min then 25 cycles of 95 °C for 20 s, 60 °C for 30 s and 72 °C for 1 min. Exact primers and native sites remain unassigned. Its geometry receives the recovered U-bearing population, not a completed canonical amplicon.
8. 🟢 Second PCR: 50 µl with 1 U Phusion, 98 °C for 2 min then 22 cycles of 98 °C for 10 s, 60 °C for 30 s and 72 °C for 1 min. The inter-round transfer fraction, dilution, cleanup and relationship of the two rounds' primer species remain unresolved. The graph exposes this boundary rather than inventing a physical aliquot or cleanup. Second-PCR geometry receives that native transfer population.
9. 🟢 Recover 200–500-bp DNA after 12% native polyacrylamide TBE electrophoresis. Fragment Analyzer checks size and a standard-curve qPCR assay quantifies final libraries; analytical aliquot volumes are unassigned. The graph represents physical gel retention and exclusion without asserting removal of every co-migrating RNA, original DNA or adapter molecule.

## Modified originals and ordinary PCR daughters

🟡 PCR copying does not recreate template 5mC protection, methylated-adaptor chemistry or uracil in ordinary DNA daughters. Converted-site DNA readout is T/A for originally converted C and C/G for retained cytosine, with strand orientation affecting which base is reported. The original U-bearing and protected methylated polymers can still coexist with daughters; their physical recovery through transfer, column and gel steps remains unresolved. The graph preserves those independent populations rather than relabelling all originals as ordinary amplicons.

The source ligated index belongs to native adapter-bearing polymers; its structured feature has no assigned coordinate, width, member or acquisition role. It is not a free index oligo, inline UMI or imported six-base cell barcode. No canonical TruSeq Y-fork, adapter strands or sequencing primer sites are supplied.

## Sequencing and computational fates

🟢 The defining article specifies paired-end Illumina HiSeq2000 sequencing. Read lengths, index acquisition and native primer sites remain unassigned. Lambda reads provide a bisulfite-conversion reference. Adapter/quality trimming, Bismark alignment to mm9 and the additional lambda reference are computational operations; they do not physically remove molecules from a library.

🟢 Three negative controls omit the cell during picking and are separately carried through library preparation and sequencing to assess contamination. They are described here as independent control preparations, not as physical discard products of the illustrated single-cell lineage. Bulk-cell RRBS and the separate sperm nested-PCR validation are outside this representative single-cell library graph.

## Bounded native-graph migration (2026-10-18, batch 12)

One published original inspected: PMC3847781 JATS Methods lines 103–113, including the separate validation-assay boundary and sequencing/data-processing passages. The mandatory source-bundle oligo scraper produced hits only in the excluded upstream schematic; none was adopted as an original library oligo. A scraper search does not establish absent information.

Removed queued `lib/rrbs.py` source defaults and the canonical adapter/readout detail panels. Local native population constructors now preserve mouse/lambda/tRNA, protected originals and ordinary PCR copies, actual incoming PCR states, source-known physical column/PAGE fates and explicitly unresolved inter-round transfer. Source positive lysis, spike-in, enzymes, reaction volumes, cycles and sequencing platform remain local.

Deferred: full publisher supplement including validation Table 5; Gu 2011 bulk predecessor; full 2015 detailed protocol; matching historical adapter/primer specifications; native ligation/PCR geometry, modified-adaptor positions, run/index cycles, inter-round physical allocation and class-specific recoveries. This capped pass is not a full protocol review or either independent missing-information pass; no publication-wide absence is accepted.
