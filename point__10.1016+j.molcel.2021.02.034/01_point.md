# POINT: polymerase-associated intact nascent transcripts

Sousa-Luís et al., *Molecular Cell* (2021), [doi:10.1016/j.molcel.2021.02.034](https://doi.org/10.1016/j.molcel.2021.02.034).

This bounded first pass reads the defining published preparation and preserves POINT-seq, POINT-5 and POINT-nano as branches of one protocol. The selected POINT-seq graph now follows a source-qualified preparation from chromatin through a reported PCR DNA library and sequencing contribution. Historical kit internals and exact read geometry remain unresolved; POINT-5 and POINT-nano retain distinct partial models.

## Sources and actual inspection

- `j.molcel.2021.02.034_PMC8122139.xml.txt`: published version of record (manifest: CC BY, author-manuscript false). Inspected limitations 323–325; key resources 405–409; STAR Methods 487–503; read processing/orientation 515–521. Reference `bib37` identifies the cited 2016 mNET-seq protocol.
- Published `j.molcel.2021.02.034_supp_mmc1.pdf.txt`: title/identity; Fig. S1 caption 346–373 (Empigen/DNase/Pol-II isolation validation); Fig. S5 caption 1303–1312 (polyadenylated RNA and pre-library >500 nt selection). Original panels were not systematically inspected.
- `mmc2.pdf` opens as the published article and combines article/supplement contents; only its identity and the repeated orientation-primer hit were inspected. It is not evidence for a separate chemistry version. The whole-bundle scraper returned 14 hits across these three files; its principal library-associated hit was checked at JATS 521.

The automatically listed 2020 preprint is superseded and excluded from current chemical evidence. Its failed retrieval is not a needed-source blocker while the published original is available. No source-absence conclusion is made.

## Isolation and molecular scope

🟢 The method begins from the mNET-seq preparation with source-specific changes. Modified NUN/Empigen treatment prepares chromatin; DNase digestion solubilizes chromatin; clarification retains supernatant. Anti-Pol-II antibody beads recover intact polymerase-associated nascent RNA, followed by six washes and two Direct-zol purifications with intervening DNase treatment. The source uses 10⁷ cells for POINT-seq/POINT-5 and 4×10⁷ for POINT-nano.

The graphs retain DNase reaction products as molecular material until clarification/extraction. Their physical pellet, IP flow-through and cleanup rejects are explicit. They do not call DNA digestion products automatic waste, infer nucleotide-resolution RNA ends, or replace the source's modified extraction with an imported mNET-seq constructor.

## Three source branches

| Branch | Published preparation | Acquisition | Local model boundary |
|---|---|---|---|
| POINT-seq | Fragment RNA at 90 °C for 10 min, targeting 150–200 nt; NEBNext Ultra II Directional RNA kit E7760S | NovaSeq 6000 | Source-positive kit RNA-to-PCR-DNA interface is explicit; exact historical kit chemistry/polymers remain unassembled. RNA below ~150 nt is lost during library preparation, not by fragmentation alone |
| POINT-5 | Random-N6/template-switch approach; SMARTer Stranded RNA-seq kit 634836; no initial RNA fragmentation; amplified library selected to 150–800 bp | NovaSeq 6000 | Stops at purified RNA before exact kit oligos/RT/template-switch/PCR states; does not inherit POINT-seq fragmentation |
| POINT-nano | E. coli PAP polyadenylation; SPRISelect 0.6× removes RNA below 500 nt; ONT direct-cDNA kit SQK-DCS109 | PromethION | Artificial poly(A) remains an RNA tail of unassigned length before the unreviewed direct-cDNA kit; no PCR amplification is assigned |

🟢 POINT-nano sequencing is of cDNA. This differs from native-RNA sequencing even though both begin with nascent RNA. The local `NascentRNA` class describes the **pre-kit RNA**, not an asserted RNA sequencing terminal. ATP as the polyadenylation substrate is an explicit biochemical inference; its concentration and exact kit conditions are not invented.

The orientation-analysis sequence `GAAGATAGAGCGACAGGCAAGT` is printed at JATS 521. It is a sequence searched in reads, not authority for filling the remainder of an ONT primer. Source processing also names NBD104/NBD114 barcodes; actual selected barcode members, sequence structures and acquisition cycles remain unassigned.

## Lineage, metadata and deferrals

Explicit component predecessor: Nojima 2016 mNET-seq, [doi:10.1038/nprot.2016.012](https://doi.org/10.1038/nprot.2016.012), cited in Methods 489 and `bib37`. Scope: initial chromatin/nascent-RNA purification strategy **with stated alterations**. This does not establish whole-recipe inheritance from a differently versioned mNET-seq catalogue entry or equal end chemistry.

Proposed metadata: RNA; bulk cells; polymerase affinity capture; one fragmented directional short-read branch, one unfragmented template-switch short-read branch, one tailed long-read **cDNA** branch; published instance. The public schematic selects POINT-seq; alternative branch prefixes remain unassembled. Exact historical guides, variant-specific molecular adapter/RT/PCR states, read cycles, index assignments and remaining source panels are deferred. No confirmed missing-information claim follows from this capped pass.

## Source-qualified POINT-seq completion, 30 October 2026

Published non-manuscript main Methods487–503 and original resource rows405–409 were re-read; mandatory bundle oligo scraping was repeated. A second reader independently checked the same source scope, selected nuclear/chromatin input boundary and kit-to-DNA statement. This is positive preparation evidence, not a full linked-source or absence review. Main493 explicitly calls the E7760S output a PCR library; kit internal enzymes, ligation/RT/PCR sequence orders, index identities and acquired read roles/cycles remain unassigned rather than filled from current vendor material.

The graph starts with the source-prepared nuclear/chromatin pellet, preserving the original model boundary. Initial crude-nuclear isolation delegates to the cited mNET-seq original and is outside the schematic. Actual NUN1/modifiedNUN2 treatment, pellet/PBS recovery, DNase/clarification, tenfold NET2E dilution, antiPolII capture/six washes, firstDirect-zol, fresh intervening DNase and secondDirect-zol are now explicit source-local stages. DNA-digestion products remain molecular material until physical workup.

Actual4150TapeStation RNA size analysis is an analytical allocation distinct from production stock. Source RNA contributes to90°C10min fragmentation, then the explicitly opaque historical E7760S library handoff. Its product is a source-reported PCR DNA-library class, not RNA relabelled as directly sequenced. SmallRNA loss is placed in the kit's physical workup, not automatically at the fragmentation reaction. RNA/library stock classes have unassigned quantities and presence; no nonzero leftover or full-stock transfer is asserted. A selected DNA contribution proceeds to sourceNovaSeq6000 with no invented run cycles or primer geometry.

The constructor's literal source inventory was checked against the actual Methods, rather than generated from graph-node enumeration. No existing required source item was removed, excluded or renamed. POINT-5's random-N6/template-switch kit and POINT-nano's poly(A)/SPRI/direct-cDNA kit remain separate RNA-prefix models; completion of this selected POINT-seq preparation is not completion of those alternatives or of the historical kit's inner molecular relations. No missing sequence or reaction detail is declared absent.
