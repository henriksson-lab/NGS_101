# RiboMethSeq (Illumina)

**Defining source:** Marchand *et al.*, *Nucleic Acids Research* (2016),
doi:10.1093/nar/gkw547 (PMC5027498) — the Illumina implementation this page draws.
**Original method (not drawn):** Birkedal *et al.*, *Angewandte Chemie* (2015),
doi:10.1002/anie.201408362, an Ion Torrent implementation. Its Wiley Supporting
Information was obtained and read on 2026-10-17; the article itself is still unobtained.

*Catalogue boundary (2026-10-17):* this entry was re-anchored from the Birkedal 2015 DOI
to Marchand 2016 because every drawn oligo, condition and step comes from the 2016
Illumina paper. Birkedal 2015 is a different library chemistry (see below) and must not
be described by this page's construct.

## Evidence used

- 🟡 The defining method's principle — protection of the phosphodiester bond immediately
  3′ of a 2′-O-methylated ribose and inference from depleted fragment endpoints — is
  taken from Marchand *et al.*'s description of Birkedal *et al.* (PMC5027498 JATS
  lines 40, 94); the Birkedal SI Figure S5 legend (lines 487–495) states the same read-end
  logic first-hand.
- 🟢 The Illumina method compares alkaline, MgCl₂ and ZnCl₂ fragmentation and selects
  alkaline hydrolysis for its precision.
- 🟢 It specifies alkaline hydrolysis of 1–250 ng RNA in 50 mM bicarbonate pH 9.2 at
  95 °C with no gel purification (30–200-nt fragments), 3′ dephosphorylation (Antarctic
  Phosphatase, heat-inactivated), 5′ phosphorylation (T4 PNK, 1 mM ATP), RNeasy MinElute
  cleanup, NEBNext Small RNA E7330S library preparation, and HiSeq 1000 SR50 or MiSeq
  paired-end runs (JATS lines 62–78); Supplementary Table S2 gives the MiSeq runs as 2×75.
- 🟢 The 39-nt synthetic 2′-O-methylated RNA `5′-CCUUAAGGCNmNAUAACAUACCAGAUCGCCACCCGCGCUC-3′`
  (JATS line 60) is a quantification substrate, not a library oligo.
- 🟡 Kit adapter and primer members remain unassigned until a version-qualified manufacturer original is checked.

## Oligo source boundary

🟢 The **original** Birkedal Ion Torrent library (Wiley SI, read 2026-10-17; see the section
below) used a kinase-dead tRNA ligase, a 5′-phosphorylated 3′ RNA adaptor and an in vitro
transcribed 2′,3′-cyclic-phosphate 5′ RNA adaptor. Its oligos are recorded below as a
description of the undrawn original only. Do not transfer them, or the later NEBNext Illumina
oligos, between the two chemistries.

🟢 Marchand *et al.*'s **later Illumina implementation** uses direct 3′ and 5′ ligation with the NEBNext Small RNA Library kit after end repair (PMC5027498, “Comparison of fragmentation approaches”). Its supplementary methods also print `5′-GTTGAATTTCTTCACTTTGACATTCAGAGCACTGGGCAGAAATCACATTG-3′` as a 25S-rRNA capture oligo for the validation assay. That capture oligo is not part of the sequenced library construct.

## Construct logic

The informative absence is local: protected bonds produce fewer 5′ and 3′ fragment ends.
End repair must occur after hydrolysis, and direct adapter ligation is essential because
tailing or template-switching would blur the nucleotide-resolution cleavage coordinate.

## Bounded curation round 4 (2026-10-16)

Inspected: Related published Marchand2016 PMC5027498 JATS1–29 (identity),62–78 (fragmentation, end repair, library and reads).

Made explicit that the displayed NEBNext/Illumina recipe is Marchand2016, not a verified reconstruction of the catalogue-defining Birkedal2015 protocol. Named Antarctic-phosphatase3-prime repair followed by T4PNK/ATP5-prime phosphorylation.

Deferred: Defining 2015 article/supplements remain unavailable in this bundle; protocol-boundary/catalogue correction and exact version-specific oligos require a separate pass. This pass does not establish a complete source review or publication-wide absence.


## Bounded curation round 5 (2026-10-17)

Retrieval retried: `get_sources.py` (DOI, both Wiley SI routes → HTTP 403), Europe PMC
(PMID 25417815, no PMC/author manuscript), Unpaywall and OpenAlex (closed, no repository
copy; also the German-edition DOI 10.1002/ange.201408362), Semantic Scholar (closed),
University of Copenhagen research portal (record only, no file). The defining article and
SI remain unobtained; manual URLs are in the source MANIFEST.

Inspected: Marchand2016 JATS lines 40–44, 58–78, 94, 104, 198 and Supplementary Data File008
Tables S2 and S4 (Birkedal2015 vs this work); File009 is figure panels only.

🟡 Marchand Table S4, second-hand about the defining protocol: ≥1000 ng of PAGE-purified
20–40-nt alkaline fragments (or 10 µg input), no end repair, “proprietary direct ligation
to 5′-OH/3′-P (>P) fragments”, Ion Torrent, <1 M reads/sample, Score A/Score C. The main
text (line 104) adds the mutant RNA ligase and a specially prepared 2′,3′-cyclophosphate
5′ adapter. These are shown on the page caveat as Marchand's description, not as checked
facts of the 2015 article.

Corrections: note no longer marks the defining principle 🟢 from an unread article; the graph
now includes the RNeasy MinElute cleanup (with discarded flow-through) after end repair; MiSeq
read layout stated as 2×75 (Table S2).

Deferred: Birkedal2015 article and Wiley SI (manual retrieval); their Ion Torrent adapters,
ligase and primers; NEBNext E7330S manual for kit oligos other than the 3′ adapter. The
catalogue identity question (displayed chemistry belongs to Marchand2016, a distinct
Illumina implementation) is recorded for maintainers rather than changed here.


## Deferred follow-up 02 (2026-10-17)

Retrieval re-checked: Europe PMC (PMID 25417815, subscription only), PMC ID converter (not
in PMC), OpenAlex (closed; DOI, PubMed and University of Copenhagen portal record only). The
Wiley article and SI remain unobtained (manual URLs in the source MANIFEST). New route: the
GEO series cited with the defining article, GSE54655 (`Series_pubmed_id` 25417815; samples
GSM1321380–85; full SOFT now in the source bundle with a MANIFEST row). SRA SRP036188
experiment records repeat the same library text.

🟢 GEO `Series_overall_design`: yeast RNA "was degraded at denaturing conditions into small
fragments … (20-40 nucleotides). The fragments were ligated to RNA oligos using a tRNA
ligase that was mutated to remove its kinase activity, reverse transcribed and the cDNA used
as input for Ion semiconductor sequencing." Samples: Ion Torrent PGM (5) and Proton (1);
mature rRNA gel-purified; data processing removes a 4-nt barcode and the 3′ adaptor
(Cutadapt v1.2, reads <15 nt discarded) and counts 5′ and 3′ read ends. This corroborates,
from a defining-article-linked record, Marchand's second-hand description (20–40-nt
fragments, mutant ligase, Ion Torrent); it supplies no adapter, primer or ligase sequence.
The page caveat now cites it. Nothing of the defining Ion Torrent construct is drawn.

Deferred: Birkedal2015 article Methods and Wiley SI (adapter structure, ligase identity,
hydrolysis conditions, RT/PCR primers); NEBNext E7330S manual; catalogue boundary
(Marchand2016 vs Birkedal2015) for maintainers.


## Deferred follow-up 03 (2026-10-17): Birkedal 2015 Supporting Information

Source: `anie_201408362_sm_miscellaneous_information.pdf(.txt)`, hand-downloaded by the
maintainer (MANIFEST row). Identity checked: printed header "Supporting Information ©
Wiley-VCH 2014", title "Profiling of Ribose Methylations in RNA by High-Throughput
Sequencing", authors Birkedal, Christensen-Dalsgaard, Krogh, Sabarinathan, Gorodkin, Nielsen,
filename anie_201408362 — the SI of the cited predecessor. The article itself is still
missing (manual https://doi.org/10.1002/anie.201408362). Read `.txt` lines 1–175 (Materials
and Methods to Sequencing and read mapping) and 455–495 (Figure S2–S5 legends); the score
formulae and qPCR/primer-extension tables (lines 175–450) were only identified.

🟢 Verbatim original oligos (lines 105–137), undrawn:

- 3′ adaptor CR30 `5′-CGCACUGAGUCGGAG` or CR35 `5′-CGAGCUGAGUCGGAG` (RNA; "Ion Torrent A
  sequence including a 4 nt barcode (#1: CGCA or #2 CGAG)"); 5′-phosphorylated with T4 PNK
  before use.
- 5′ adaptor P1 IVT `5′-GGCCACUACGCCUCCGCUUUCCUCUCUAUGGGCAGUCGGUGAU>P` (5′ triphosphate,
  2′,3′-cyclic phosphate; 43 nt; HDV-ribozyme release). Line 108 prints `…GUCGGTGAU>P` with
  a T; line 137 prints U. 🟡 The IVT template (line 117) encodes this position as T, which is
  U in the transcript, so the T is read as a typographical DNA letter (page caveat shows both).
- IVT-template primers C894, C887, C895 and the AtRNL mutagenesis primers C813–C818 are
  reagent-construction oligos, not library parts.
- 🔴 The RT primer ("short version of the Ion Torrent adaptor A with a 4 nt barcode (#1; CGCA)
  complementary to the 3′ adaptor CR30") is not printed base-by-base.

🟢 Original steps (lines 76–160): 1 µg/µl RNA + equal volume 100 mM Na₂CO₃ pH 9.0, 90 °C
1–10 min; ethanol/KAc precipitation; 10% denaturing PAGE, 20–40-nt band, elution in near-neutral
TEN/PCI to preserve cyclic phosphates. 3′ adaptor ligation: 1 µg fragments, 100 pmol
phosphorylated CR30, AtRNLΔPNK (S701A), 14% PEG 8000, 1 mM ATP, 30 °C 1 h; PCI/ethanol.
One T4 PNK reaction phosphorylates the ligation product 5′ end and removes cyclic phosphates
of unligated fragments; 10% urea PAGE. 5′ adaptor (P1 IVT) ligation under similar conditions;
10% urea PAGE. FastAP removes the junction 2′-phosphates. First-strand cDNA (Invitrogen kit)
with the barcoded short Ion A primer; qPCR quantification (Ion Library Quantitation Kit);
emulsion PCR (60 cycles), Ion Torrent/Proton, 250 or 500 flows. Reads without barcode
discarded; Cutadapt removes the 3′ adaptor (<15 nt discarded); Bowtie 2 to RDN37-1; 3′ and
5′ read ends counted, the 5′ ends shifted one position (Figure S5).

🟡 Marchand Table S4's "no end repair" is consistent with no repair before the first ligation,
but the SI shows T4 PNK between ligations and FastAP afterwards; this is a summary difference,
not a published error.

Page change: the caveat now gives this first-hand SI description instead of Marchand's
second-hand account. The drawn Marchand graph is unchanged.

Deferred: Birkedal2015 article main text (Angewandte; manual DOI); RT primer bases; Schutz et
al. (SI ref. 7) tRNA-ligase cloning predecessor; NEBNext E7330S manual for drawn kit oligos;
catalogue identity (Marchand2016 vs Birkedal2015) for maintainers.

## Local validation migration (2026-10-18)

Re-read defining Marchand 2016 PMC5027498 Methods lines 62–80 and direct-ligation
rationale lines 100–104; original supplementary File008 Table S2 key and rows 1–4.
Ran the source-bundle oligo scraper. The drawn alkaline branch now owns its
30–200-nt fragment population, precipitation, Antarctic phosphatase/inactivation,
T4 PNK with ATP and MinElute cleanup. No gel size selection is inserted before repair.

Direct 3′/5′ adapter ligation, RT and PCR are explicitly scoped to the named E7330S
kit route; exact kit members, end chemistry, revision and primer landing remain
unassigned until the matching manufacturer original is checked. The previously
shown generic 3′ adapter string had no version-qualified local source and is not
promoted into this study's molecular declaration. Birkedal's different Ion Torrent
chemistry remains undrawn, with its existing discrepancy notice preserved.

Table S2 assigns different libraries to HiSeq 1000 SR50 and MiSeq PE2×75. The graph
keeps separate run/library alternatives; it does not assign both settings to every
physical sample or invent an index read length. Source-review checklist roles now
reflect the existing Marchand defining-DOI boundary, while historical passes and
blocked absence-review obligations remain intact. This migration does not establish
full reaction validation or missing-information acceptance.
