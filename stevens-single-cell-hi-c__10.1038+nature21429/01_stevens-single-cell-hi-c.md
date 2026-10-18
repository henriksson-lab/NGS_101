# Stevens single-cell Hi-C

## 1. What it is

A single-cell Hi-C workflow in which nuclei are isolated before their individual restriction, ligation and amplification reactions.

## 2. Sources and evidence

- 🟢 Stevens et al., *Nature* 2017, DOI [10.1038/nature21429](https://doi.org/10.1038/nature21429), defining paper and methods.
- 🟢 The contact digest is **MboI**, followed by biotin-14-dATP fill-in and proximity
  ligation. AluI is a later fragmentation step on captured contact DNA.
- 🟢 The supplement prints the customized 3-bp-tagged Illumina adapters and both PCR
  primers. It also reports Nextera XT dual-index preparation as an alternative route.

## 3. Distinguishing chemistry

The essential distinction from Nagano 2013 is early isolation of individual nuclei. The
contact junction is the biotin-filled MboI junction; AluI must not be drawn there.

## 4. Bounded source pass — 2026-10-17

Read published PMC5385134 JATS identity and introductory scope100–114,
and the published Springer supplementary Methods PDF `MOESM26_ESM.pdf.txt`
Single cell Hi-C protocol192–219, both library variants227–292, population
boundary294–295 and source adapter/primer tables361–410. Manifest identifies
both originals. The initial archived Methods DOCX extraction was inspected
at 1–5; embedded citation fields obscure its text, so the published PDF was
used for chemistry. Ran the bundle oligo scraper, checked used sequences
and modifications against the PDF tables, and excluded encoded citation
strings as oligos.

🟢 The selected route sorts single haploid G1 nuclei, images them, overlays
low-melting-point agarose and performs MboI digestion, Klenow large-fragment
fill-in and T4 ligation in wells. Overnight 65 C crosslink reversal precedes
M280 capture and removal of agarose-containing solution. The graph now
represents the agarose condition and crosslink reversal before capture.
Other nuclei remain outside this lineage with unresolved fate, not a claimed
discard.

🟢 For the standard adapter route, AluI fragmentation is followed by Klenow
exo-minus dA tailing, T4 adapter ligation, Platinum Pfx 25-cycle PCR and 0.7x
AMPure XP recovery. Pooled libraries are double-size-selected at 300–700 bp and
sequenced on MiSeq with 2x75-base reads. These choices now appear in the
graph. The separately described Nextera XT route has distinct KAPA/index-PCR
stages; those settings are not imported into this adapter route.

🟢 The independently printed CAA adapter is
`ACACTCTTTCCCTACACGACGCTCTTCCGATCTCAA*T` paired with
`pTTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG`; p means 5′ phosphate and * is
a phosphorothioate linkage. The drawing now uses its complete 17-base core
and both unpaired arms rather than a generic 12-base adapter stem. The
printed library PCR primers carry terminal phosphorothioates; source-derived
primer graph states preserve them. The historical endpoint helper is reused
only for these independently source-supported arms and CAA tag.

🟡 The old ligation-detail panel showed a PCR-completed endpoint immediately
after ligation. It now leaves forked pre-PCR ends unrepresented and reserves
the completed molecule for the final-library view. This corrects a repository
model error; it does not validate the whole ligation/PCR transition.

Deferred: remaining article/source sections, imaging/selection predecessors
and vendor originals; complete agarose-removal/wash intermediates and exact
forked ligation products; population Hi-C and Nextera XT branches. No source
absence is accepted and no complete chemistry validation is claimed.


## Bounded native migration — 2026-10-18 batch08

Inspected 2 verified originals: Published PMC5385134 JATS identity/introduction100–113; original PMCCloud Supplementary_methods.docx text1115–1121(single-cell/imaging/agarose/contact/capture, standardadapter and NexteraXT branches, pooling/selection/PE75), ordered adapters1128–1165 and PCRprimers1166–1173. Embedded citationfield sections excluded from chemistry reading; these actual procedure paragraphs inspected. Manifest identities checked; mandatory bundle scraper reviewed.

Removed shared ContactWorkflow/endpoint/batchseg. Local MboI5′GATC/fill/contact graph has actual sourceY-adapter strands and unknownfold PCRgeometry on real ligated input; dT/dA junction retained in derived daughter endpoint. Native barcodefeatures/modifications/physical washes/recovery retained. Published NexteraXT alternative now separately constructed/rendered from common captured-contact precondition, not an invented physicalaliquot: AluI/wash→tagment/neutralize/wash→KAPA9cycles→beadremoval1×cleanup→KAPA16cycles→1×cleanup→pool/doublesize selection. No CAA tag or assigned indexwidth imported intoXT; installedi5/i7features appear only afterPCR. SourceMiSeqPE75 retained; indexacquisition unassigned.

Deferred: Complete native forked/XT strand geometry and exact study indexmembers/widths/acquisition; Component-specific original/offproduct/agarose recovery and sourcevariants beyond representedroutes; Imaging/selection predecessors, vendor originals and populationHiC route not freshly reread. This source-capped pass is not a complete protocol review or an acceptance of source absence. Historical audit paragraphs above retain their original scope; this entry states the current model.

Independent implementation review retained source terminal phosphorothioate bonds on the physically incorporated forward/reverse PCR-primer footprints. Original adapter backbone modifications and biotin are not copied into daughter DNA. Native Y-ligated originals retain their own source A*T bond; this is distinct from the C*T bond supplied by a PCR primer.
