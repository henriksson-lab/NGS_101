# PAR-CLIP: photoactivatable ribonucleoside crosslinking

Hafner et al., *Cell* (2010), [doi:10.1016/j.cell.2010.03.009](https://doi.org/10.1016/j.cell.2010.03.009).

This is a bounded first source pass. The local graph captures the 4SU crosslink/recovery front end and deliberately stops before the unresolved library implementation. It is not a sequencing-ready schematic. Source availability and remaining linked sources are recorded in `catalogue/source_reviews/par-clip__10.1016+j.cell.2010.03.009.json`; no missing-disclosure claim is accepted.

## Sources inspected

- Defining article `j.cell.2010.03.009_PMC2861495.xml.txt`: Summary 59–61; Introduction 71–79; Results 101–103 (T→C controls); Discussion 179–189 (mutation interpretation and 6SG); Methods PAR-CLIP 201–205; library-predecessor citation 365–367.
- The defining chemistry supplement `NIHMS195398-supplement-Supplement.doc` and Supplementary Tables S1–S7 are actual linked originals, but the standard fetcher could not retrieve them. Their manifest gives each gated URL. Their contents were not inspected; tables were not judged irrelevant from filenames.
- Directly cited predecessor Hafner et al., *Methods* (2008), [doi:10.1016/j.ymeth.2007.09.009](https://doi.org/10.1016/j.ymeth.2007.09.009), `doi_10.1016+j.ymeth.2007.09.009/j.ymeth.2007.09.009_PMC2847350.html.txt`: overview 1143–1163; platform-specific Solexa oligos 1352–1399; 3′ ligation 1402–1430; 5′ ligation 1439–1463; RT/hydrolysis/PCR 1470–1504. Its explicit alternatives are retained as alternatives; the missing PAR-CLIP supplement is not silently replaced.

The complete PAR-CLIP bundle scraper returned zero sequence hits. The predecessor bundle scraper found 16, and the Solexa-specific hits were checked in the original table. Zero hits in the incomplete defining bundle establish no absence.

## Defining chemistry

🟢 Live cultured cells incorporate 4-thiouridine, followed by 365 nm UV irradiation. Isolated crosslinked RNA yields enriched T→C changes in sequenced cDNA. The original also explores 6-thioguanosine, with G→A changes and lower crosslink efficiency; the local graph scopes the principal 4SU route only.

The mutation is a sequencing signature, not a declaration that irradiation physically converts every RNA U into C. The original observes background changes in non-crosslinked material and context-dependent variation; the graph therefore does not rewrite the RNA sequence or assume complete modification/crosslink occupancy. The source recommends interpreting conversion enrichment across read clusters rather than replacing high read counts with an assertion of binding.

## Original preparation order

1. 🟢 Grow FLAG/HA-tagged RBP-expressing HEK293 cells overnight with 100 µM 4SU; irradiate live cells at 365 nm.
2. 🟢 Lyse with NP40, digest with RNase T1, immunoprecipitate with anti-FLAG/Protein G Dynabeads, then digest again on beads.
3. 🟢 Calf intestinal alkaline phosphatase dephosphorylates the RNA; PNK with radioactive ATP labels it. The compact graph records these operations without claiming nucleotide-resolution labeling occupancy.
4. 🟢 SDS-PAGE, electroelution, proteinase K digestion, acidic phenol/chloroform and ethanol precipitation recover RNA. The graph preserves the biological-RNA provenance and explicitly accounts for IP flow-through.
5. 🟢 The article says the recovered RNA was converted to cDNA “as described” by Hafner 2008 and Solexa sequenced. The 2008 platform-specific library uses preadenylated blocked DNA 3′ adapter, RNA 5′ adapter, RT, RNA hydrolysis and PCR. Applicability of every detailed alternative to each PAR-CLIP instance remains deferred until the defining supplement is read.

## Predecessor Solexa inventory: supporting source, not substituted PAR-CLIP constants

🟢 Hafner 2008's platform-specific table prints:

| Component | Verbatim sequence |
|---|---|
| 3′ adapter | `AppTCGTATGCCGTCTTCTGCTTG-L` |
| 5′ adapter | `rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC` |
| PCR 5′ primer | `AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA` |
| PCR 3′ primer | `CAAGCAGAAGACGGCATACGA` |

Its table defines `L` as a 3′-OH blocking group; Fig. 1 caption identifies an aminolinker. Its protocol permits truncated Rnl2 and the K227Q alternative for the first ligation, then Rnl1 with ATP for the second. Which alternative was used for the defining PAR-CLIP experiments is not assigned here. These historical Solexa PCR adapters are not modern P5/P7 sequences with missing tails filled from a contemporary kit.

## Reported readout, metadata and remaining model work

🟢 The original discusses 32-nt Solexa sequence segments as representative (Discussion 187), and its mutation-control example discusses 40-nt reads (Results 103). These contextual values are not promoted to one universal run length. The current graph exposes no exact sequencing plan, UMI, sample-sheet index or instrument primer.

Proposed source-supported metadata: RNA; bulk cultured cells; metabolic 4SU labeling; 365 nm crosslink/IP; RNase fragmentation; ligation/RT/PCR library as an explicitly cited downstream route; Solexa sequencing; published instance. Keep `pending_source_review` until the complete molecular/readout scope is curated.

Deferred: retrieve the defining supplement and linked tables; compare its actual library choices with the 2008 alternatives; separate source-known retained/rejected gel fractions and RNA-end chemistry; build exact adapter/RT/PCR and sequencing states. No additional fetcher or challenge workaround was attempted.

## Component ancestry candidates

Explicit predecessor: Hafner 2008 (DOI above) → this PAR-CLIP directory, limited to recovered-RNA cDNA library preparation; evidence is original Methods 205 and reference 365–367. The predecessor's own sequence tables were read, but the original supplement remains unavailable and exact variant inheritance is not established. That predecessor does not yet have a registered directory.

Conceptual CLIP precedent is explicit in Introduction 71–73 (Ule 2003 and later CLIP originals); its component scope is UV crosslink/IP/RNA-fragment isolation, not adapters or photoactivatable labeling. Existing eCLIP's original also cites PAR-CLIP as an earlier CLIP improvement, but that comparative citation does not prove direct component inheritance. Similar primers establish no ancestry.
