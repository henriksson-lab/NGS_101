# nano-COP: direct sequencing of nascent RNA

Drexler et al., *Molecular Cell* (2020), [doi:10.1016/j.molcel.2019.11.017](https://doi.org/10.1016/j.molcel.2019.11.017).

This is a capped accepted-manuscript source pass. The graph retains native RNA and its poly(A)/poly(I) branch identity through supportive RT, then stops before the complete historical ONT motor-adapter/readout model. Final published procedural applicability and gated supplements remain unresolved.

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

Proposed metadata: RNA; bulk cells; metabolic labeling plus chromatin enrichment; affinity purification; alternative poly(A)/poly(I) tailing; splinted adapter ligation; supportive RT; **native RNA** nanopore readout; no amplification; published instance. Keep `notes`/`pending_source_review`. Deferred: VOR comparison, gated Tables S1/S4, exact adapter/RT/motor geometry, rRNA-kit mechanics and remaining source review. No accepted missing-disclosure claim.
