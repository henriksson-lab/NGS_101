# STARR-seq: source review and implementation queue

Arnold CD et al. (2013).
[Genome-wide quantitative enhancer activity maps identified by STARR-seq](https://doi.org/10.1126/science.1232542).
Science 339:1074–1077.

## Current evidence boundary

A limited author-deposit graph is now implemented below; this is not a completed chemistry review.
The defining article's published title, authors and DOI have been identified.
The source fetcher has not retrieved its full original or supplement in this pass.
A university-hosted published PDF is identified, but the download timed out;
that is a route failure rather than proof of publisher restriction or source absence.

The authors' project page,
[Arnold Science 2013](https://starklab.org/data/arnold_science_2013/), identifies
the same published article and links source-owned vector records for DSCP, hsp70
and human STARR-seq. This establishes named retrieval targets, not interchangeable
constructs. No vector sequence has yet been read or assigned to the model.
The page also links large sequencing datasets; none is needed for this source pass.

## Named follow-up targets

- Defining article, DOI 10.1126/science.1232542: original main text and figures.
- Its published Supporting Online Material: construction methods, primer table,
  readout and source-specific vector identity.
- The author-linked vector record for the modelled experiment. Keep DSCP, hsp70
  and human variants separate; verify each record's source identity before use.
- The author's GEO accession, if required to locate the published preparation
  metadata; raw sequencing data are outside this initial review.

Run `tools/scrape_primers.py` when an article/supplement bundle has been retrieved
and before interpreting its chemistry. An empty bundle is not a negative oligo audit.

## Intended implementation scope

Model reporter-library construction and the RNA-derived sequencing path, with a
separate input-library path where the source supports it. Keep the tested insert's
identity distinct from a separately assigned reporter barcode in barcoded MPRA.
Derive any transcript boundaries, intron/splice changes, amplification junctions,
adapter segments and read acquisition from the selected published source variant.

No actual primer, vector sequence, ligation geometry or run recipe is assigned by
this acquisition-only note. Exact molecular assignments, a finished public schematic and independent
source review remain pending. No confirmed-missing information is claimed.

## Bounded acquisition follow-up, 2026-10-21

Retrieved and read the author project page, its linked DSCP vector GenBank
original (`STARR-seq_DSCP.gbk`), and deposited input-library sample
[GSM1000405](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1000405).
The external manifest identifies each original URL. These are three originals;
this pass does not substitute a later STARR-seq recipe for the published method.

🟢 GSM1000405 is an approximately 500-bp genomic input library for the OSC
experiment. Its deposited extraction protocol specifies shearing/size selection,
NEB E6000L preparation with Illumina PE-400-1001 adapters, PCR addition of
In-Fusion homology arms, recombination into the DSCP/sgGFP STARR-seq vector,
and plasmid recovery. The record names GAIIx acquisition and paired-read mapping.
It does not by itself assign the expression-library primer set or splice geometry.
No complete preparation-kit chemistry is borrowed from its family name.

The DSCP vector has an internal source discrepancy: its LOCUS line declares
4,080 bp, while counting all uppercase and lowercase bases in ORIGIN gives
5,592 bp, matching the `source 1..5592` feature. Some other annotated coordinates
also exceed the LOCUS length. Both values are retained here; the vector has not
been assigned wholesale to the model or silently repaired. The file contains many
primer annotations, including alternative/control primers; their presence alone
does not establish which were used for a particular published library.
The main text/supplement and sequence-level applicability checks remain necessary.

Ran `tools/scrape_primers.py`: it scanned zero eligible files in this bundle.
The GenBank file and deposited plain-text record were therefore read directly;
zero scraper hits are not evidence about primer disclosure.
The publisher candidate supplement route returned HTTP 403; the university
article route again timed out. These are recorded route failures, not absence.
Next targets are the published article/supplement, the expression-library sample
methods, and applicable vector-region/primer assignments. No molecular graph or
confirmed-missing claim has been added from these partial records.

A subsequent openly indexed CiteSeerX candidate copy (the paper's title matched
the search result) returned HTTP404 through the normal validated PDF fetcher.
No candidate document was saved or used as evidence. Its full request URL and
failure remain in the external manifest and the retrieval declaration remains in
`catalogue/extra_sources.tsv`; the original main/supplement requirements remain open.

## Original OSC deposit and limited graph

Read the original GSE40739 title/identity and overall-design fields, and the
complete OSC sample GSM1000403 record. The current sample exports state a
last-update date of 2019-05-15; they are accession-identified author records, not
an archived 2013 revision of the published Methods. It is **OSC biological replicate2**,
not replicate1: sample molecule poly(A) RNA, description amplified STARR-seq
transcript cDNA, instrument Genome Analyzer IIx. The extraction field repeats
the reporter DNA-cloning preparation and does not itself assign the RT/PCR
geometry. The input sample GSM1000405 is the source approximate500-bp input
library used for OSC. Its identity, description and series membership are
checked independently of the transcript sample. The series explicitly separates
Drosophila paired-end STARR-seq from human single-end STARR-seq; no human or
later protocol recipe is merged into this graph.

`tools/starr_seq.py` now declares the source-supported reporter construction,
separate RNA/input readout paths and paired acquisition, and renders native SVG.
The entry remains `notes`: the graph is reviewable locally but is excluded from
the public site and sequence tools. All exact oligo, adapter/index, splice,
transcript and primer-site assignments remain open. The author GenBank is not
imported wholesale or silently repaired; its recorded discrepancy remains above.

Known source operations are explicit: genomic fragmentation/size selection,
kit-named preparation, adapter ligation, homology-arm PCR and In-Fusion. Kit
preparation, exact In-Fusion chemistry and biological propagation are unexpanded
boundaries. The reporter assay and RNA-to-cDNA conversion are also unexpanded;
no particular RT primer, intron removal, transcript end or RNA-selection recipe
is inferred. RNA is shown as a single strand; cDNA strand state is left
unassigned, rather than manufacturing a duplex before it is established.

Input plasmids and transcript cDNA remain separate populations. Source PCR
readouts have unassigned primer reagents and sequencing-library structures.
Paired acquisition is recorded without cycles, sequencing primers or index
architecture. No independent MPRA tag or UMI is added by analogy. The expected
reporter arrangement is inferred, not an exact vector map or a verified genotype.

Seven local checks cover the limited claim set and original accession/population
identities. No exact oligo state is available for transcription checking; the
source guard reports seven unassigned inputs. Chemistry remains unresolved,
not established by a green source-identity or graph-shape check. Reran the
primer scraper after fetching the deposits; it found zero eligible source files,
so the original deposited text was read directly. This is not a negative oligo audit.

The published article/supplement, exact original kit documents and applicable
RNA/readout primer definitions remain required before this becomes a finished
public chemistry page. No confirmed-missing conclusion is accepted.

## Supplied published originals, 2026-10-22 — current model scope

This section supersedes the acquisition-only limitations above where explicitly resolved.
The uploaded published article (`science.1232542.pdf`, own DOI verified) and corrected
SOM (`1232542.arnold.sm.revision1.pdf`) are now available in the external source store.
The SOM cover says published 17 January 2013, corrected 6 February 2013, with a
posting-date correction. The older uploaded `arnold.sm.pdf` is retained as history;
its chemistry has not been assumed equivalent wholesale. Twelve supplementary workbooks
were extracted from the uploaded archive; their individual tables remain unreviewed.

🟢 Read article first-page description and Fig1A caption; corrected SOM printed pp5–8
(screening vector, input construction, RNA recovery, RT and nested PCR) and Mapping p11.
Ran the primer scraper on the full source directory and checked its relevant pointers
against the printed primer strings. Table S1/S12 hits concern luciferase validation;
they are not automatically the pooled STARR-seq readout primers.

🟢 The Drosophila source places candidates downstream of sgGFP in the reporter 3′ UTR,
uses DSCP/mhc16/sgGFP/ccdB context, Covaris S220 fragmentation, a 500–700-bp size class,
PE-400-1001 adapters with E6000L preparation, exact cloning primers and AgeI-HF/SalI-HF
recipient preparation before In-Fusion. The source says eight cloning PCR reactions,
then “the five PCR reactions were pooled” in the same p6 paragraph. The page preserves
both counts and models the common reaction class without selecting a batch count.

🟢 Oligo-dT25 beads select poly(A) RNA before Turbo DNase; they are not the RT primer.
RT uses SuperScript III and `CAAACTCATCAATGTATCTTATCATG`, followed by RNase A.
The reporter first PCR uses `AAGCCACCATGGAAAAGGCCAT` and
`TATCATGTCTGCTCGAAGCGG`; terminal five and three nucleotides respectively have
phosphorothioate protection as stated in the source (the printed base strings have
no bond symbols). The reporter forward spans the mhc16 splice junction. Input uses
`CTAGGCACACCGAAACGACTAAC` and the same reverse sequence. Source modifications of
the different input forward are not assigned by analogy. Second PCR uses NEBNext
Multiplex Oligos E7335L and KAPA HiFi HotStart ReadyMix. Exact selected historical
kit oligos/index values remain unassigned. Mapping specifies 36-mer paired reads for
Drosophila, versus 50-mer single reads for human; these variants are kept separate.

The notes-only graph now represents these reagents, RNA selection/DNase/RT/RNase,
separate first and second PCRs and 36/36 acquisition. Exact full vector/adapter states,
source-derived intron/splice coordinates, primer landing sites and sequencing endpoints
remain unfinished modelling work. The published p5 vector sequence is available and
must be assessed directly; its existence is not a missing-information claim. No whole
GenBank map is silently repaired or substituted. This bounded review is not an independent
full source review. Entry remains notes until those exact molecular assignments are ready.

## Partial DNA geometry and source splice conflict, 2026-10-22

🟢 Inspected the corrected SOM rendered original pp5 and8. The complete2502-base
Drosophila cassette printed on p5 matches the author-linked DSCP GenBank ORIGIN
exactly. This establishes inheritance of that specific cassette, not validity of the
GenBank LOCUS length or every annotation. The input-DNA forward site lies in that
printed cassette. Both15-nt In-Fusion homologies also match it. The short downstream
first-PCR reverse/RT context is contiguous in the same author record and matches
the independently printed p8 primer sequences. Only that short context is used;
the complete GenBank is neither checked in nor imported into the molecular model.

🟡 The intended input-DNA construct now has these source-qualified anchors and the
primer-written adapter contexts, separated by explicit unknown intervals. Its input
PCR sites face inward and its partial amplicon follows the full primer ends. All
source-to-intended-product placements remain inferred, so the molecular check is
conditional; the schematic character positions are not genomic/plasmid coordinates
or a predicted full amplicon length. Whole kit/ligation states remain unresolved.

The original author record has two different intron annotations,399..469 and402..468
(one-based inclusive). Applying the latter67-base deletion to the source-matching
cassette gives `AAGCCACCATGGAAAAGGCCAA` at the proposed forward junction; the
published p8 primer is `AAGCCACCATGGAAAAGGCCAT`. They differ at the3′ terminal
base. Applying the other annotation also does not give a complete exact primer site.
The note and rendered page preserve both source values. Neither intron annotation
is silently chosen as the experimental transcript, and the primer is not repaired.
No inference about mismatch tolerance, experimental failure or confirmed absence
follows. Exact reporter-RNA PCR geometry remains unresolved for this reason.

This pass used two originals (corrected SOM and the already-linked author record),
with a strict primer search against the corrected source. It is a bounded geometry
check, not either independent full source-review pass. Local suite now has19 checks.

## 2026-10-30 public preparation review

The OSC published preparation can be shown with its existing source-qualified workflow and literal reviewed inventory. This bounded pass inspected the published article identity/Fig1A, corrected SOM Methods pp5–8 and Mapping p11 (rendered pp6/8), and the author DSCP GenBank identity/cassette/tail scope. The bundle oligo scraper was rerun. The public page retains the source cloning-PCR count conflict, unresolved reporter splice-primer context, and explicit GenBank LOCUS4080 versus ORIGIN5592 discrepancy. It does not import the full map or assign historical kit adapter/index sequences. This is not a completed independent source review or an accepted missing-information claim.
