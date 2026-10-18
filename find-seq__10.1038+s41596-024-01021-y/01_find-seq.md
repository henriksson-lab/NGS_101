# FIND-seq — research notes

## Active native migration, 2026-10-18

Selected demonstrated Smart-seq3 scFIND route. Native 5′acrydite-T6-ACGAGCATCAGCAGCATACGA-T30V is matrix linked; the 5′biotin-AGAGACAGATTGCGCAATG-N8-rGrGrG TSO is a separate physical oligo. Source forward primer ends A*T*G (two phosphorothioate bonds); reverse source typography “…AC*G*−3′” is retained rather than used to infer an unprinted bond. Copied TSO DNA retains UMI identity while shedding the original RNA/biotin chemical annotations.

Actual RNA and genomic DNA are independent populations. Matrix-bound ss cDNA is distinct from free WTA daughters. Source cold wash1, two wash2, two RT-buffer washes, EDTA stop and three Tween washes precede detection-PCR sorting. Diagnostic hydrolysis probes and marker PCR products are not sequencing barcodes or library arms. One sorted droplet is a single-cell allocation; bulk bins, validation aliquots and cycle-calibration reactions are separate experiments. Freeze/thaw/60°C aqueous recovery does not cleave the acrydite linkage. Actual recovered matrix-associated ss cDNA enters WTA, followed by source0.8X AMPure and20µl elution. Original/matrix/RNA/marker populations have unresolved species-specific bead recovery.

Step102 delegates to Nextera XT manufacturer instructions. The exact native guide revision, ATM inventory, run and cleanup/PCR internals are not established by this capped pass; a bounded commercial stage represents that genuine source boundary. A canonical final duplex or acquired index recipe is not substituted. Smart2 alternatives and source enrichment-marker optimization remain outside the selected graph.

This bounded pass read 3 originals and reran the mandatory bundle primer scraper. Original/source identity was checked in the external manifest. Historical review status and the two-reader missing-information gate are unchanged.

## Historical observations

The earlier observations below retain their original review scope; current source-owned declarations above replace earlier canonical endpoint and recipe defaults.


Primary protocol: Shin *et al.*, *Nature Protocols* (2024), [10.1038/s41596-024-01021-y](https://doi.org/10.1038/s41596-024-01021-y) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC11537836/)). Defining study: Clark *et al.*, *Nature* (2023), [10.1038/s41586-022-05613-0](https://doi.org/10.1038/s41586-022-05613-0).

## 1. What it is

Focused interrogation of cells by nucleic-acid detection and sequencing (FIND-seq) physically sorts rare cells using RNA- or DNA-specific TaqMan fluorescence, while preserving their polyadenylated transcriptomes for downstream bulk-bin or single-cell RNA sequencing.

## Sources read

- 🟢 2024 protocol, PMC11537836 JATS text (`s41596-024-01021-y_PMC11537836.xml.txt`): development overview and Experimental design (lines 87–115), Materials → Oligonucleotides (lines 155–183), drop breaking, RT reaction, bead counting and cDNA validation WTA (lines 600–727), bead re-injection, TaqMan thermocycling and droplet sorting (lines 733–861), WTA cycle calibration, WTA of sorted material and DNA library preparation (lines 863–931), Anticipated results (lines 979–987). Bounded pass 2026-10-17.
- The Materials list gives the Smart-seq2 and Smart-seq3 capture, template-switching and amplification oligos with their 5′ and terminal modifications, plus the env and Ψ TaqMan sets.
- 🔴 Not inspected: the Clark *et al.* 2023 Nature paper (earlier RNA-marker implementation), the Smart-seq2/Smart-seq3 predecessor protocols, the Nextera XT kit guide and the publisher DWG device drawings. An earlier version of this note said the Nature Methods had been read; no source record supports that, so treat it as unread.
- The source fetch commands are `python3 tools/get_sources.py --doi 10.1038/s41596-024-01021-y` and `python3 tools/get_sources.py --doi 10.1038/s41586-022-05613-0`.

## Construct-changing and selection steps

1. 🟢 Acrydite-modified oligo(dT) is covalently linked to allyl agarose with APS/TEMED.
2. 🟢 A bubble-triggered microfluidic device co-encapsulates single cells, molten agarose and lysis buffer (1,000 mM LiDS, 2 μg/μl Proteinase K; 55 °C). Cooling at 4 °C gels a genome-capture bead: genomic DNA is physically entrapped and poly(A) RNA hybridizes to the oligo(dT) matrix (line 101).
3. 🟢 The emulsion is broken with 20% PFO and the beads are washed cold in high-salt buffer, then in 5× RT buffer; the washes remove LiDS/Proteinase K, which inhibit cDNA synthesis (lines 103, 600–625).
4. 🟢 RT on the beads: Maxima H Minus, 2 μM TSO, 1 M betaine, 7.5% PEG-8000, 6 mM MgCl₂; 30 min at room temperature then 90 min at 42 °C; stopped with EDTA (lines 635–665). The demonstrated run used the Smart-seq3 TSO (line 981). An aliquot is checked by WTA with 14/16/18 cycles (lines 691–727).
5. 🟢 Beads are re-injected (~70 μm drops) with multiplex TaqMan PCR (TaqPath ProAmp; demonstrated Ψ-FAM and env-VIC IPDA assays) and thermocycled 88 °C 10 min, then 55× (88 °C 30 s, 59 °C 1 min) (lines 733–795).
6. 🟢 Double-positive droplets are sorted either in bulk into one tube (FIND-seq) or singly into PCR strips or plates (scFIND-seq). The marker amplicon is a selection signal, not the transcriptome library.
7. 🟢 Sorted material is frozen at −80 °C, melted at 60 °C for 10 min, and the aqueous phase is amplified with the Smart-seq3 forward/reverse primers and KAPA HiFi (bulk 16–20 cycles, single-droplet 20–26 cycles to start). "Amplified double-stranded WTA products are free in solution, distinct from cDNA bound to agarose" and are purified with 0.8× AMPure XP (lines 863–925).
8. 🟢 Tagmentation and library preparation use the Nextera XT kit "following the instructions provided by the manufacturer" (line 931).

## Published oligos

- 🟢 The Materials list prints both TaqMan probes with an `MGBNHQ` quencher label (env: 5′-VIC-CCTTGGGTTCTTGGGA- MGBNHQ −3′; Ψ: 5′-FAM-TTTTGGCGTACTCACCAGTMGBNHQ-3′). The page reproduces that label as printed.
- 🟢 **Smart-seq2 branch:** 5′-acrydite-TTTTTTAAGCAGTGGTATCAACGCAGAGTAC(T)₃₀VN-3′ capture primer; 5′-AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3′ template-switch oligo; 5′-AAGCAGTGGTATCAACGCAGAGT-3′ Smart PCR primer.
- 🟢 **Smart-seq3 branch:** 5′-biotin-AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG-3′ template-switch oligo; 5′-acrydite-TTTTTTACGAGCATCAGCAGCATACGA(T)₃₀V-3′ capture primer; 5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G-3′ forward primer and 5′-ACGAGCATCAGCAGCATAC*G*-3′ reverse primer. The source marks `*` as a phosphorothioate bond. The eight random bases in the Smart-seq3 TSO supply a UMI.
- 🟢 The same Materials list gives exact env and Ψ TaqMan primer/probe sets for the demonstrated viral-marker selection. They assay the bead before transcriptome library preparation and do not enter that library; other marker choices require their own primers and probes.

## Published naming discrepancy

- 🟢 Materials (lines 167–169) lists the acrydited Smart-seq3 oligo(dT) primer separately from the Smart-seq3 forward primer, which has a Nextera handle and no acrydite or dT. Anticipated results (line 981) says the agarose was "conjugated with the oligo dT primer (Smart-seq3 forward primer)".
- 🟡 The model treats the parenthesis as a naming slip, because only the acrydited oligo(dT) can be linked to allyl agarose and capture poly(A) RNA. The schematic page shows this discrepancy.

## Final library and read layout

- 🟢 The final transcriptome is a Nextera XT library made from Smart-seq3 WTA product. The procedure uses the Smart-seq3 primers; Smart-seq2 oligos are listed in Materials only.
- 🔴 The protocol text gives no sequencer, read lengths or UMI read layout within the inspected sections. The Smart-seq3 forward primer carries the Nextera Read 1 handle plus TSO sequence, so 5′ UMI fragments could be read without tagmentation of that end. This was not modelled because no read layout is stated.
- 🟡 The public page uses the standard Nextera primer geometry for the tagmented internal cDNA fragment, with no acquired cycle counts. Selection-target TaqMan primers do not become part of that library.

## Boundary

FIND-seq here means nucleic-acid cytometry for rare-cell transcriptomics, not the unrelated “full interrogation of nuclease DSBs” acronym.
