# nano-COP: direct sequencing of nascent RNA

Drexler et al., *Molecular Cell* (2020), [doi:10.1016/j.molcel.2019.11.017](https://doi.org/10.1016/j.molcel.2019.11.017).

This is an accepted-manuscript-qualified preparation. The selected poly(I) graph now follows the source through supportive RT, motor-adapter attachment, two RNA cleanups and flow-cell loading of native RNA. Complete occupied molecular geometry, final published procedural applicability and gated supplements remain unresolved. Poly(A) and no-tailing alternatives remain separate source-local prefixes.

## Sources and actual reading

- `j.molcel.2019.11.017_PMC7060811.xml.txt` is an accepted author manuscript, explicitly marked unedited at 269, not the publisher version of record. Inspected introduction/results 44–68; 4sU purification 160; native-RNA library methods 176; read processing 180/192; Fig. 1 legend 688; historical kits in key resources 764. The whole available bundle scraper returned zero oligo hits.
- The manifest records four actual gated supplement originals (one PDF and three workbooks). They were not read or declared irrelevant. The main article explicitly locates custom ONT_oligoA/ONT_oligoB_C10 in Table S4 and per-sample kit/tail assignments in Table S1; those remain blocked.

No direct adapter sequence is declared absent because the incomplete source bundle yielded zero hits. No preprint or modern ONT manual is substituted for the historical source.

## Source preparation and preserved physical material

🟢 Cells receive an 8-min 4sU pulse. Cellular fractionation enriches chromatin-associated RNA. RNA undergoes Biotin-HPDP treatment, chloroform/isopropanol recovery, streptavidin capture, source washes, DTT elution and miRNeasy cleanup with DNase treatment. This selects recently transcribed chromatin-associated RNA; the local model does not turn 4sU into U or assign positions/occupancy to the label.

Released/soluble reaction material remains molecular material until a documented extraction, precipitation or wash separates it. The graph marks those physical rejects explicitly. Source rRNA depletion uses Ribo-Zero Gold or RiboMinus v2 depending on sample; Table S1 is gated. Its internals and rejected-material fate remain a marked unresolved operation rather than an invented kit mechanism.

## Native-RNA branches

| Branch | Tail preparation | Ligation recognition |
|---|---|---|
| Poly(A) | NEB E. coli PAP, 37 °C 1 h for unperturbed samples; source also reports Takara PAP, 37 °C 7.5 min for PlaB/DMSO samples | Provided T10 splinted adapter |
| Poly(I) | Yeast PAP with **ITP**, 37 °C 30 min | Custom C10 splinted adapter from Table S4 |
| No added tail | No enzymatic poly(A)/poly(I) addition | Captures an endogenous-tail-compatible subpopulation; no universal ligatable end is assumed |

The default local graph scopes poly(I); `graph_for_tail` also exposes A and no-tail controls. The current A instance represents the specified unperturbed NEB route; the alternative Takara reagent/timing remains a named source instance, not a separate protocol entry.

🟢 Historical SQK-RNA001/RNA002 libraries use T4 DNA ligase for 3′ adapter ligation and SuperScript III RT. The source explicitly says the resulting cDNA is **not sequenced**: it reduces RNA secondary structure. Motor-adapter ligation and cleanup follow. Native RNA travels through the pore from its 3′ end toward its 5′ end; the source runs FLO-MIN106 flowcells for 48 h. No fixed base-cycle count follows from that duration.

The local `NativeRNA` class preserves **inosine** in the artificial tail, distinguishes C10 from T10, and retains RNA as the sequenced polymer after supportive RT. It does not replace the tail with generic DNA/N bases, assign tail length, or assume the complete gated adapter geometry. RNA/I recognition is a source-specific state declaration, not a new global IUPAC complement rule. The exact supporting RT-primer identity remains unreviewed.

## Lineage and readout limits

Explicit conceptual predecessor: Garalde 2018, [doi:10.1038/nmeth.4577](https://doi.org/10.1038/nmeth.4577), adopted for direct RNA sequencing at manuscript Results 52/56 and reference R27. Scope: native-RNA nanopore readout strategy. Historical RNA001/002 components are not equated with modern RNA004 adapters. POINT-nano's cDNA sequencing and nano-COP's native RNA readout remain chemically distinct.

Nascent RNA 3′ ends support transcription-position/splicing analyses, but alignment/basecalling uncertainty, native polyadenylation and splicing intermediates require the source's computational filtering; none is a physical molecular discard in the graph. The present model does not implement those analyses.

The public scope is the source-selected poly(I) native-RNA preparation, not complete chemical reconstruction of every historical kit or sample. Deferred: VOR comparison, gated TablesS1/S4, exact adapter/RT/motor geometry, rRNA-kit mechanics and remaining source review. No accepted missing-disclosure claim.

## Source-positive sequencing preparation (bounded follow-up, 2026-10-30)

🟢 Associated-manuscript Methods176 explicitly purifies the supportive-RT sample with RNACleanXP beads, ligates a sequencing adapter with preloaded motor protein, performs a second RNACleanXP purification, resuspends in ElutionBuffer, mixes with RNA RunningBuffer and loads a primed FLO-MIN106 flowcell for48h using MINKNOW direct-RNA defaults. These operations are now modeled. Native RNA is the sequencing-library substrate; its copresent supportive cDNA is not a substitute DNA endpoint. No fixed base-cycle count is inferred from run duration.

The selected source500ng RNA contribution and final loading contribution retain separately marked untransferred stocks with unassigned presence/amount/fate. SourceHPDP and DTT are molecular contributions; precipitation, columns and washes provide the actual physical rejects. The two named custom adapter oligos and preloaded motor/adapter are components of their aggregate source inventories, not invented additional free-reagent doses. Stock count/packaging remains unassigned unless the original explicitly specifies a preloaded stock.

The partial C10 recognition view places unknown context on both sides: it does not assign a reviewed whole TableS4 oligo, terminal polarity or exact ligation-ready end. Supportive RT remains an opaque source-positive RNA-to-cDNA operation with exact priming/nucleotide relation unassembled, avoiding a guessed fresh primer. Motor attachment likewise does not select an unreported fresh ligase or modern motor identity. RNA001 versus RNA002 and rRNA-kit alternatives remain unresolved sample assignments from unread TableS1; this is not a whole-version-equivalence claim.

A literal reviewed input/step inventory covers the selected Methods160/176 operations only. No SourceItem was deleted or excluded. The existing alternative A/no-tail prefixes retain their earlier limitations. This follow-up uses the same one accepted-manuscript original and mandatory cached-bundle scraper, with no acquisition or source expansion.
