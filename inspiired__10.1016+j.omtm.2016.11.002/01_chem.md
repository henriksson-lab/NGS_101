# INSPIIRED: integration-site wet preparation

Sherman et al. (2017 issue; online2016), **INSPIIRED: A Pipeline for Quantitative Analysis of Sites of New DNA Integration in Cellular Genomes**. DOI[10.1016/j.omtm.2016.11.002](https://doi.org/10.1016/j.omtm.2016.11.002). This is the biochemical-preparation paper, not the companion.003 analysis/visualization paper.

This private incomplete source graph follows the published preparation, with a conditional HIV-vector entry and source iSL1 linker inventory example; these are not assigned to one clinical study sample. Quantitative inference from independently sheared junctions is distinct from physical barcode chemistry. No affinity, perfect suppression, purity or clonal-frequency prediction is made.

## Selected source scope

🟢 Published PMC5363316 Methods177–185; published supplementalProtocol `mmc1.pdf.txt`, steps1–10 and TablesS9–12; workbook `mmc2.xls.txt` selected iSL1 row only. Same-article combined `mmc3.pdf` header identified as a duplicate presentation, not an independent chemistry authority. Mandatory bundle scraper164 candidates; TableS9 was also rendered from the original PDFp6 (`/tmp/inspiired19-blocker.png`) because chemical notation matters. Exact source manifests live externally; [retrieval](ref/MANIFEST.md).

Source TableS9 HIV blocker,5′→3′:

`(bC)A(bG)(bG)GACTTGAAA(bG)(bC)GAAA(bG)(bG)(bG)AAA(bC)C-3′AmMO`

Nine positions are bridged nucleic acids. Source explicitly defines(bC) as iBNA-meC and(bG) as iBNA-G. The sugar/base modification identities and3′amino block remain in the typed oligo and returned strand; the blocker is not an ordinary extensible DNA primer. Source says blockers are optional but effective; the conditional selected branch includes one, without asserting perfect suppression. New PCR daughter DNA does not automatically inherit the blocker chemistry.

TableS10 negative linker: `5′Phos-GTCCCTTAAGCGGAG-3′AmMO`. Positive iSL1 contains the printed sample-specific fixed region, anN12 random region and common `CTCCGCTTAAGGGACT`, ending in the source3′dT overhang. The Scene computes the15bp duplex from the reverse complement, retaining the3′T and blocked negative terminus. The randomN12 sequence is not called a UMI by assumption; the source's fixed unique linkers identify sample preparations. Replicate indexes are12nt Golay sequences, with member/body/landing unassembled.

## Source operational route

1. Entry is isolated source gDNA, here conditionally HIV-vector-containing material; uninfected gDNA and no-template controls are prepared independently. Optional pre-sonication re-purification is outside the selected graph, not silently imposed on all inputs.
2. CovarisM220 fragmentation targets800–900bp. Source0.7×AMPure purification precedes historical NEBNextUltraEndRepair/dA-Tailing20°C30min and65°C30min. Main Methods establish5′phosphorylation and3′dA. Source unique annealed linkers undergo overnight16°C Blunt/TA ligase treatment, followed by0.7×AMPure purification.
3. The source says at least four PCR replicates. Actual allocation/scaling is explicitly unresolved: SI lists15µL ligation eluate, then four PCR1 reactions each consuming15µL. The model does not manufacture60µL from15µL or silently change the source values.
4. Each represented replicate undergoes PCR1 with sample-specific linker primer, source HIVLTR primer and optional modified blocker. SI describes5“linear” and20“exponential” cycles; stringent MLV80°C settings are not assigned to the HIV example. Its printed HIV70/67°C alternative is retained with exact stage assignment unassembled.
5. Approximately5µL PCR1 diagnostic gel material is separate from the production remainder.2µL source PCR1 template enters nestedPCR2, which adds a replicate-specific Golay barcode and Illumina arms;5“linear” and15“exponential” cycles are retained as printed. A separate blocker stock aliquot is supplied, while carried PCR1 companions are not assumed gone. Other PCR1 fractions retain unresolved fates.
6. PCR2 gets its separate diagnostic gel aliquot.10µL from each replicate is allocated to the within-sample pool; unused fractions are not discarded by assertion. The pool is0.7×AMPure purified and20µL eluted. PicoGreen analytical material is separate; equal-mass or desired-coverage sample pooling precedes independent ScreenTape/KAPAqPCR analyses. Only production material reaches acquisition.

🟢 The source calls the instrument MiSeq with a300-cycleV2 kit and prints R1=179, I1=12, R2=143 cycles. These values are preserved together; their total and kit label are not silently reconciled into a current validated run. The private preview states this limitation. Full run/primer landing, exact Golay members, all modified oligos and historical kit geometries remain unassembled.

## Ancestry and deferrals

The paper develops integration-site isolation with internal-vector blocking and sonication-defined junction counts; no source recipe is imported from TRIP inversePCR, UDiTaS tagmentation or the companion computational pipeline. Direct predecessor integration-site methods, remaining workbook rows and source kit documents are deferred, not confirmed absent. Source guard limitations for(bC)/(bG) notation are reported separately from original-transcription checks and chemical completeness.

## Source-format transcription exception

Two independent readers inspected the published ProtocolSIp6/TableS9 whole HIV blocker and its modification legend. Eight repeated one-strand blocker inputs (four replicate PCR1/PCR2 pairs) have a narrow source-format exception because the guard does not parse the printed(bC)/(bG) notation. The full original modified-string test remains, and the actual strand preserves every bridged position and3′AmMO. This is not a geometry or absence exemption and no synthetic source text was created. Final guard result:3other strands checked,0failed,8strands exempted,23inputs without exact states across the two models.
