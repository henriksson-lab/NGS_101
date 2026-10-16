# SNU-Seq

[Gerlevik et al., NAR (2026)](https://doi.org/10.1093/nar/gkag652) enriches RNA containing pulse-incorporated 4sU, adds artificial poly(A) tails and uses a 3′-region RNA library kit. Intentional sonication is omitted. A labelled molecule may include pre-existing RNA tracts; enrichment does not prove that every nucleotide was synthesized during the pulse.

## Original versions

The actual [Lorenz preprint v1](https://www.biorxiv.org/content/10.1101/2021.07.14.452379v1.full.pdf) was recovered separately and its printed **July 14, 2021** header verified. Its methods (current text 1310–1380) report bacterial PAP on 150 ng RNA for 45 min at 37°C, isopropanol recovery, thirteen PCR cycles and **Ion Chef/Ion Proton**. They do not specify Illumina FWD or SR100 in that passage. Its Hep3B pulse conditions and yeast spike-in are version-specific.

The previously cached v2 is printed **October 24, 2025**, despite the 2021 DOI. That version and the 2026 supplement describe Ion Torrent and Illumina FWD options. The 2026 original Supplementary Methods (text 306–429) reports thirteen PCR cycles and single-end 100-nt acquisition on Ion Proton or Illumina NextSeq 2000. The graph selects that later +bPAP route, with dated kit references qualified below. The size-selected branch and no-bPAP/no-4sU/rRNA-depletion controls are separate source variants.

## Capture, release and PAP

The 2026 SNU route inherits its described TT-seq enrichment while omitting sonication: HPDP-biotinylation, cleanup before affinity capture, µMACS streptavidin-column washing, **DTT release**, miRNeasy recovery with on-column DNase I, then bacterial PAP and precipitation. The graph begins at recovered total RNA after the pulse and source spike-in addition.

The unbound/wash fractions are physical separations. The eluted fraction is enriched for pulse-labelled RNA; complete capture, release or separation of all transcript age classes is not established. Remaining column contents and individual RNA-degradation-product fates are explicitly unresolved.

### Source-version wash difference

The printed 2026 Supplementary Methods p. 18 (text line 338) specifies **1 mM NaCl** in the capture wash. The printed Dölken 2008 Methods p. 1970 (text 966–968) specifies **1 M NaCl**. Both original pages were visually inspected. The 2026 procedure explicitly includes modifications, so this is a qualified source-version difference, not a proven typographical error. The schematic displays both values and does not select a concentration.

The TT-seq 2016 supplement delegates purification to Dölken 2008 rather than printing this entire recipe itself. Dölken's original room-temperature capture and RNeasy recovery are not substituted for the 2026 conditions and miRNeasy recovery.

## Two kit-specific declarations

The [Ion Torrent guide 012UG036V0130](https://www.lexogen.com/wp-content/uploads/2018/06/012UG036V0130_QuantSeq-IonTorrent.pdf), printed June 13, 2018, describes:

- P1-side oligo(dT)-primed first-strand synthesis.
- RNA-template degradation with RS1/heat/RS2.
- A-side random-primed second-strand synthesis, introducing an inline barcode.
- Bead cleanup, PCR and post-PCR bead purification.

Its Appendix E prints the finished architecture and a `TCAG–barcode–CGAT–RNA` read start. The barcode family is 9–11 nt; no separate index read is required in this guide. The actual study barcode member and physical loaded sequencing primer remain unassigned.

The [Illumina FWD guide 015UG009V0260](https://www.lexogen.com/wp-content/uploads/2021/08/015UG009V0260_QuantSeq_Illumina_2021-08-16.pdf), revision table dated August 16, 2021, describes oligo(dT) first strand, RNA removal with RS/heat, random second strand using the Read 1 side, bead cleanup, PCR-introduced i7 and post-PCR purification. Appendix G supplies the finished library and complete sequencing-primer strings. The PDF also includes a separately numbered critical-recommendations insert; its presence does not change the main guide's printed revision identity. Optional UMI, UDI and REV/CSP modules are not selected.

**Inferred study assignment:** both displayed endpoint constructors use these source-positive dated guide architectures as conditional references. Exact experimental kit revisions, full ordered kit reagents and selected barcode/index members remain unassigned. This is model data awaiting assignment, not absence of published vendor architectures. Fixed reference bases default to inferred; biological RNA and copied poly(A) lengths stay variable. The two platform branches are joined by an explicit warning-producing comparison boundary, not a claim that one experimental sample was physically split.

## Readout and processing limits

Both chosen guide families read **towards poly(A)** from the random second-strand side. FWD does not use the REV/CSP read direction. A longer read can be needed to reach the native RNA 3′ boundary; the first acquired base does not universally identify it.

The 2026 analysis clips poly(A), masks genomic A-rich regions, filters alignments and computes the **aligned-read 3′ coordinate** using `bedtools genomecov -3`. The 2021 original likewise uses `-bg -3`. Those operations do not prove native-end reach for every read or physically remove library molecules. The computed FWD reference preview stops at the variable biological interval instead of inventing read-through into a fixed tail.

## Linked source scope

The original GSM5240711 record confirms an Ion Torrent QuantSeq SNU sample with yeast spike-in, biotin capture and PAP. The three article-linked GEO series were scoped separately: GSE165251 contains nuclear/cytoplasmic polyadenylated-RNA controls, GSE172053 contains multiple Hep3B assays, and GSE179306 includes SNU, size-selected and TT-seq contexts. Their child sample/experiment bodies are not all inspected. No raw reads or result archives were downloaded.

The actual author repository landing and embedded README describe preprocessing and analysis for the later study; child script bodies, archived revision equivalence, further vendor sources and remaining supplementary contents remain open in the review. No missing-information claim is accepted.
