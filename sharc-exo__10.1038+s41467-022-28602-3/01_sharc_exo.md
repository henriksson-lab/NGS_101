# SHARC-exo

Van Damme et al. (2022), *Chemical reversible crosslinking enables measurement of RNA 3D distances and alternative conformations in cells*, DOI https://doi.org/10.1038/s41467-022-28602-3 . SHARC and SHARC-exo/no-exonuclease controls are branches of this original.

## Sources

🟢 Published main PMC8854666 Methods175–223, overview/Fig2, published SI opening reagent characterization/Fig1 synthesis, and reporting-summary identity/date were inspected. The reporting summary is not a new chemistry recipe, and unchecked OCR numbers/plots are not source facts. Mandatory scraper inspected original modification pointers; two large source-data files were skipped by default and remain unreviewed. Published authority supersedes the automatically retrieved preprint pointer; no preprint chemistry is imported.

## RNA crosslinks, trimming and contact readout

CDI-activated dipicolinic acid (DPI) forms reversible ribose2′OH crosslinks, with monoadducts/unmodified molecules and reaction remnants coexisting. In-cell source concentrations0/5/12.5/25mM30minRT are distinct source experiments, not deterministic site/distance assignments. A selected two-chain constructor receives actual RNA sequences, anchor positions, and explicit5′/3′ end chemistry. It accepts caller-selected exonuclease-retention bounds only when they retain the supplied anchors; it does not predict where RNaseR stops. Removed suffix strings identify residue origins, not a claim of actual eliminated product formulas.

Proximity ligation selects explicitly supplied3′OH/5′phosphate partners and creates a new RNA phosphodiester junction while retaining the crosslink. Source reversal removes the chemical bridge while retaining this new RNA sequence/junction. Reagent-hydrolysis products remain explicitly unassembled. A complementary-insert reference is not an operational full firstDNA product: RT-primer-derived bases and exact copying geometry stay unassembled. No fixed five-nucleotide offset, guaranteed trimming resistance, universal intramolecular origin, perfect ligation or complete crosslink reversal is assumed.

## Source operation order

In vivo: DPI labeling precedes source solution removal/PBS washes, GuSCN/needle TNA lysis, PK1mg/ml37°C1h, phenol/isopropanol TNA recovery, TURBODNase37°C20min and RNA recovery. Protein/DNA digestion products persist until those physical recovery steps. RNaseIII/5mMMn2+37°C5min fragmentation is stopped with phenol, followed by ethanol recovery and8%/16%denaturing MOPS DD2D gels; the above-diagonal crosslinked fraction is physically excised/extracted.

RNaseR/ATP5mM45°C2/12/24h conditions and no-RNaseR controls remain separate. The P4P6 in-vitro branch explicitly splits DD2D RNA into two physical fractions and uses37°C2h trimming for one; its paired untrimmed companion is preserved. Trimming products remain until source ethanol recovery. Rnl1 proximity ligation overnightRT precedes boil2min, precipitate removal/ethanol recovery, then1xboratepH11/45°C2h reversal. Do not replace this with the separate model-oligo boratepH10/37°C experiment. Ethanol recovery follows.

Main209 adapter is `/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/`. Terminal3ddC contributes an actual C nucleotide plus3′ block, preserved in the local adapter; it is not erased as a bare annotation. Rnl1 room-temperature3h ligation precedes RecJf/5′Deadenylase37°C1h and ZymoRNAclean. Digestion is not immediate physical removal of free adapters.

SSIV uses a **barcoded RT primer**, original body/identifier arrangement still unassembled, with Mn2+25°C15min/42°C10h/80°C10min. RNA and reaction remnants persist until RNaseH/A/T1 digestion37°C30min and ZymoDNA cleanup. CircLigaseII60°C100min/80°C10min precedes source-cited libraryPCR, then6%nativeTBE175bp-and-above purification. Source-known late steps reach the sequencing endpoint even though exact primer/end/backbone/PCR geometry remains unreviewed. No spacer-free circle or modern kit is invented. Main221 explicitly reports MiSeq for P4P6; this platform is not silently assigned to the in-vivo branch. Acquired cycles/read/index layout remain open.

## Ancestry and limits

SourceMethods181/201/217 explicitly adapts or cites actual TNA extraction, ligation and PCR predecessors1/2/3. Their exact versions/components remain unreviewed; no PARIS/CLIP recipe is imported by endpoint resemblance. PARIS/CLIP integration in7SK validation is not imposed on every SHARC library. No full-publication promotion, two-reader absence audit or accepted missing-information claim.

## Bounded unlinked-page adjudication (2026-10-30)

🟢 Published main175–217 supplies actual DPI, RNaseIII/Mn/buffer, RNaseR/ATP/buffer, proximity-ligation ATP/SuperaseIn/buffer, adapter-ligation DMSO/PEG/DTT/SuperaseIn/buffer, free-adapter-digest buffer/inhibitor and RT/circularization co-components. The model now separates these physical molecular, catalyst and buffer classes rather than hiding co-present nucleotides inside enzyme labels. Aggregate reagent inventories are bookkeeping with unassigned packaging/stock count, not undocumented premixes or duplicate fresh doses.

RNA,1pmol barcoded primer and1ul10mM dNTPs enter the source65C5min/ice preparation before SSIV, Mn buffer, DTT and SuperaseIn. RT consumes this mixture without a new primer/dNTP addition. CircLigaseII, actual1ul50mM MnCl2 and buffer are separate supplies. Source-delegated PCR now has generic positive primer, catalyst, nucleotide and buffer functional classes; no exact reference3 recipe, primer body, nucleotide composition or cycle count is borrowed. Known selected PCR-DNA library kind survives native-gel recovery.

Actual source rApp/physical terminal ddC, chemical-bridge-versus-RNA-junction distinction and P4P6 paired fractions remain intact. In-vivo0/5/12.5/25mM concentrations are distinct source experiments; the generic trace does not assign one actual dose or represent them as a mixture. Full concentration/control producers and quantitative stock/QA allocations remain unresolved scope. MiSeq remains qualified to P4P6; no cycles or in-vivo platform are invented.

The original pending predecessor/barcoded-RT/circle/PCR-source item remains required and uncovered. No graph-generated source coverage, complete library geometry, absence claim or endpoint-only promotion is used. No PAGE_FILE/page builder exists, and none is introduced for this incomplete model.
