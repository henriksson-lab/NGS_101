# SMART-seq sequence review for the animation

31 October 2026. Bounded source/component review; no protocol chemistry was edited.
This is not a complete protocol or missing-information audit. The original SMART-seq
movie and the new SMART-seq2 example must retain separate source/version identities.

## Originals inspected

Source root: `/husky/henriksson/NGS101_papers/sources/`. The bundle primer scraper was
run for `smart-seq2__10.1038+nmeth.2639`, `smart-seq__10.1038+nbt.2282` and
`nextera-xt-dna-library-prep`; findings below were checked in the original text twins.
Upstream schematic/FASTA scraper hits were excluded as evidence.

| Original | Exact inspected location | Identity and scope |
|---|---|---|
| [Picelli 2013 article](https://doi.org/10.1038/nmeth.2639) | `smart-seq2__10.1038+nmeth.2639/nmeth.2639.pdf.txt`, 598–693 | Publisher version of record, DOI/title checked in existing manifest; Online Methods RT, preamplification, tagmentation/enrichment |
| [Picelli 2013 Supplementary Table 2](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.2639/MediaObjects/41592_2013_BFnmeth2639_MOESM301_ESM.xlsx) | `nmeth.2639_41592_2013_BFnmeth2639_MOESM301_ESM.xlsx.txt`, 1–16 | Published TSO orders and modifications, particularly `rGrG+G` row10 |
| [Picelli 2014 protocol](https://doi.org/10.1038/nprot.2014.006) | `nprot.2014.006.pdf.txt`, 325–348 | Published Reagent Setup only; later procedure/back-end not imported into 2013 |
| [Illumina Adapter Sequences, document1000000002694v22](https://support-docs.illumina.com/SHARE/AdapterSequences/1000000002694_22_ilmn-adapter-sequences.pdf) | `nextera-xt-dna-library-prep/illumina_1000000002694_22_ilmn-adapter-sequences.pdf.txt`, 250–283, 13318–13344, 13504–13555, 22460–22477 | Shared canonical-sequence original, September2025; Nextera transposon/PCR orders, index-family names/example rows and explicitly obsolete Epicentre section |
| Ramsköld 2012 cached PMC manuscript | `smart-seq__10.1038+nbt.2282/nbt.2282_PMC3467340.xml.txt`, 109–117 | Existing manifest identifies author-manuscript PMC Cloud material. Selected Methods excerpt read; this pass does not newly establish whole-version equivalence with publisher version. |

The three SMART-seq2 originals supply the scoped protocol reading. The vendor original
was checked separately for the existing shared canonical sequence definitions. The
original SMART-seq excerpt was checked only to understand its distinct branch.
No new source was downloaded; originals remain outside the repository.

## RT primer, TSO and amplification handle

🟢 Picelli2013 Online Methods630–645 prints the anchored oligo(dT) primer, specifies
that V is A/C/G and N is any base, and names SuperScriptII. Lines662–666 print the
ISPCR primer and name KAPA HiFi HotStart ReadyMix for preamplification.

- The known ISPCR handle is `AAGCAGTGGTATCAACGCAGAGT`.
- The oligo(dT) is that handle, then `AC`, then30T, then `VN`.
- The entire5′ handle/linker is physically present before RT and remains unpaired
  when only poly(dT)/VN binds the RNA. It must not appear only after extension.
- A concrete VN member and transcript are illustrative choices, not recovered
  study-specific sequences. The VN anchor should pair to bases preceding poly(A).

🟢 The2013 Supplementary Table2 `rGrG+G` row prints **AC** before two riboguanosines
and one LNA guanosine. The2014 Reagent Setup331–335 prints **ACAT** before the same
modified terminal triplet. Reagent Setup337–347 independently prints the same
oligo(dT) and ISPCR as the2013 article.

The existing module deliberately supplies ACAT while marking its trailingAT inferred
relative to the defining2013 paper. A new movie can explicitly use the **2014 TSO
form**, but should say so. It must not describe ACAT as the unqualified2013 order.
The GGG base identities and the RNA/RNA/LNA sugar modifications are separate facts.
The TSO is not an ordinary all-DNA primer that is being extended: RT copies it as a
template. The RNA–cDNA hybrid should remain present during this copying step.

## Nextera components and final-library anatomy

🟢 The2013 article677–693 names **Nextera DNA Sample Preparation**, Zymo DCC-5,
N7xx/Index1, N5xx/Index2, PPC, and the enrichment programme. Its first incubation is
72°C3min, followed by98°C30s and five PCR cycles. This is not a blanket authorization
to substitute the later NexteraXT procedure. These paragraphs do not print complete
adapter/index-primer sequences or assign a study-specific index pair.

🟢 The shared vendor original272–283 explicitly prints both transferred
`S5+ME`/`S7+ME` orders and both index-PCR orders:

- N5: P5 + `[i5]` + S5.
- N7: P7 + `[i7]` + S7.

The component strings should be imported from `lib/illumina.py` and `lib/nextera.py`,
not retyped in the animation. The vendor's index-family tables map N7xx and N5xx to
**Nextera Index Kit** and give eight-base example rows. This is positive support for
the canonical example architecture. It does not identify which index pair or source-era
revision was used for the2013 study. At this initial checkpoint, the `final_library()`8N windows were treated as
example placeholders. The focused historical-source follow-up below subsequently
verified source-era eight-base widths and component orders; study-specific index-pair
assignment remains distinct from that positive finding.
PPC is a named reagent whose full ordering is not assembled here; do not invent a
further product component for it.

For the defined top strand, show all nine source-model components, preserving order:

1. P5.
2. i5 placeholder.
3. S5.
4. ME.
5. Selected cDNA insert.
6. Reverse-complement ME.
7. Reverse-complement S7.
8. i7 placeholder (orientation must follow the represented strand).
9. Reverse-complement P7.

The opposite strand is computed from those components. For placeholder sequences,
retain identifier identity rather than treating placeholder characters as actual
sample bases. Index interpretation is not inferred from colour.

## Index-PCR continuity: derived geometry

🟡 This is an explicit sequence derivation from the supplied canonical primer orders,
not a new claim about experimental yields or exact historical kinetics.

Let A=S5+ME, B=S7+ME and I=the retained insert. Gap-repaired strands are
`A + I + rc(B)` and `B + rc(I) + rc(A)`. Let F=P5+i5+S5 and R=P7+i7+S7.
The first index-PCR extension produces:

- `F + ME + I + rc(B)` from the reverse template.
- `R + ME + rc(I) + rc(A)` from the forward template.

Each daughter initially carries only its **own newly introduced tail**. The first
extensions do not magically create both P5/i5 and P7/i7 outer tails on one duplex.
After melting and copying a tailed daughter in the following cycle, both outer tails
are present. Show that incorporation, or clearly label a transition as a later-cycle
product while retaining molecule/sequence correspondence.

The final canonical length is computed as
`len(P5)+8+len(S5)+len(ME)+len(I)+len(ME)+len(S7)+8+len(P7)`.
With the present definitions this is136+insert length; primer-binding footprints are
only the3′ S5/S7 regions. The5′ P5/index and P7/index overhangs need not bend away
from a correctly placed primer: a straight oligo can retain a visibly unpaired tail.

## Original SMART-seq remains a distinct example

The inspected cached2012 Methods describes tailed oligo(dT) priming, MMLV terminal
additions, SMARTerIIA template switching and modified **Epicentre Nextera**. It calls
for PB stripping before SPRI and nine cycles of standard Nextera PCR. It does not
justify using SMART-seq2's known TSO/ISPCR ordering as the original commercial kit.
The exact native tail extent and historical kit identity remain modelling questions,
not accepted absent information.

A useful deferred source target is already cached: the vendor original's **obsolete
Epicentre** section22460–22477 prints different transferred adapters, optional-barcode
adapter orders and PCR primers. Those family components can guide a later historical
identity check. They must not be silently swapped into the2012 movie without linking
the exact source recipe/version. This note does not establish that equivalence.


## Focused source-era Nextera follow-up: positive gap closure

A second capped pass inspected three originals: the already cached2013 Supplementary
Table5, a source-era Illumina kit insert, and a source-era Illumina sequence letter.
The latter two were read through the web tool from publicly reachable originals.
No modern revision was substituted for an old document. Whole-source downloads and
registration were requested from the coordinating worker; until cached, their durable
locators are the URLs and printed page numbers below. Attempts to inspect rendered
PDF pages via the web screenshot tool returned cache misses; the searchable PDF text
and explicit table headers/values were inspected, not rendered screenshots.

1. **Picelli2013 Supplementary Table5**, external twin
   `nmeth.2639_41592_2013_BFnmeth2639_MOESM304_ESM.xlsx.txt`19–37, was reread.
   Rows27–29 identify FC-121-1031, FC-121-1011 and FC-121-1003 for the actual
   SMART-seq2+Nextera preparation.
2. **[Using Nextera DNA Sample Preparation Kits and Enabling Dual Indexing,
   part15029964 RevC](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexteradna/15029964_C_Nextera_DNA_Sample_Preparation_Insert.pdf)**,
   printed©2012, both pages read. Page1's kit tables name those same catalogue IDs;
   page2 specifies eight index cycles for each side. This establishes the source-era
   kit workflow, not the acquired cycle configuration of every deposited sample.
3. **[Illumina customer sequence letter, printed May24,2012](https://www.seqanswers.com/filedata/fetch?id=304385)**,
   printed pages1–3 read. The filename's apparentMay23 date is not used as identity;
   the original letter itself printsMay24. Printed page3 explicitly associates the
   transferred-adapter orders with FC-121-1031/1030 and the index-PCR orders with
   FC-121-1011/1012. The base orders match the current shared canonical definitions.
   Its index tables have eight-base entries and distinguish adapter bases from
   sample-sheet bases.

These readings positively establish source-era support for both transferred adapters,
N5/N7 primer architecture and eight-base index widths for the source-named kit family.
The earlier concern that those properties were supported only by a2025 sequence
listing is therefore superseded within this precise component scope. They do not
recover the per-sample index assignment or disclose every NPM/PPC composition.

A concrete source-kit example can replace genericN windows for an exact-copy animation:
N501 has adapteri5 `TAGATCGC`; N701 has adapteri7 `TCGCCTTA`, with sample-sheeti7
`TAAGGCGA`. This is an available illustrative kit pair, not a claim about a particular
study well. Use the adapter-oriented string to construct an oligo. The reverse
complement relationship must be computed for the strand/read being represented.

The letter's individual-sequence publication permission on printed page2 requires this
notice with presentations/publications: **Oligonucleotide sequences ©2007–2012
Illumina, Inc. All rights reserved.** Do not redistribute the entire letter in this
repository. Sequence reproduction and attribution remain distinct from an assertion
that arbitrary synthetic replacements reproduce proprietary reagent modifications.

The full Nextera guide15027987 RevB was sought using its original vendor URL and the
linked LaJolla mirror; both were unavailable through the web tool. A search-visible
mirror transcription identifies that revision and also describes two8-base indices,
but was not promoted to original evidence. The accessible2016v01 replacement was
not substituted. These deferred full-guide targets do not negate the positive
component findings from the actual2012 sequence letter and kit insert.
