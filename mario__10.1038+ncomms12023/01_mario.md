# MARIO — research notes

Defining source: Nguyen *et al.*, [10.1038/ncomms12023](https://doi.org/10.1038/ncomms12023).

## Evidence used

- 🟢 MARIO maps RNA–RNA contacts mediated by a single protein in unperturbed cells.
- 🟢 UV crosslinking preserves RNPs; a biotinylated RNA linker is ligated so a selected chimera has the form RNA1–linker–RNA2.
- 🟢 Streptavidin purification selects linker-containing chimeras and a random barcode supports PCR-deduplication.
- 🟢 The defining paper's Methods section publishes every custom library oligo, rather than merely naming them: the internally biotinylated RNA linker, phosphorothioate-protected linker-complement oligo, activated RT adapter, three RT-primer variants, Cut_oligo, DP5/DP3 and full Illumina PE1/PE2 PCR primers.
- 🟢 The RNA linker is printed as `rCrUrArG/iBiodT/rArGrCrCrCrArUrGrCrArArUrGrCrGrArGrGrA`; the DNA oligo used to remove terminal linker is `T*C*G*C*ATTGCATGGGCTACTAGCAT`.
- 🟢 The pre-adenylated RT adapter is `/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/`. The three 5′-phosphorylated RT primers contain random bases, the source's four-base sample barcode, a BamHI-bearing body, and the adapter-binding end.
- 🟢 The final amplification uses the fully printed Illumina PE PCR Forward Primer 1.0 and Reverse Primer 2.0, followed by paired-end HiSeq sequencing (read length not stated). PE 2.0 carries no i7 index; samples are separated by the inline four-base barcode. The exact PE Read2 sequencing-primer assignment remains unresolved in this bounded model; its acquired biological read is preserved as unplaced, while the compatible TruSeqRead1 site is computed from PE1.0.

## Protocol boundary

MARIO is protein-contact capture, not a psoralen duplex assay; its shared linker is therefore represented separately from PARIS, SPLASH and LIGR-seq.

The oligo table on the schematic preserves the source's internal modifications and all three experimental RT-primer variants. The diagram draws the ES-1 sample; the other variants change identifiers rather than the reaction topology.

- 🟢 **Published discrepancy:** the original Methods section ‘Reverse transcription’ gives the generic RT-primer form `/5Phos/NNXXXXNNNN…` and says the random part has six bases; Supplementary Fig. 1 calls the whole inline barcode 10 nt. The original ‘Oligonucleotides used’ list prints the ES-1 primer `/5Phos/NNAGGTNNN…`, which has five random bases plus a four-base sample index (nine total), while the ES-2 and ES-indirect primers each have six random bases (ten total). The model copies the explicit ES-1 oligo; whether the ordered ES-1 primer actually had five random bases remains uncertain.

## Additional source conflict and linked-source audit

- 🟢 **Published enzyme conflict:** Methods ‘Proximity ligation’ (JATS line 149) says RNA 3′ phosphate is removed using PNK phosphatase activity, but its mixture names T4 PNK (3′-phosphatase minus). The preceding linker-ligation section (line 143) explicitly distinguishes this variant from phosphatase-active PNK. [NEB M0236 product information](https://www.neb.com/en-us/products/m0236-t4-polynucleotide-kinase-3-phosphatase-minus) says the named variant has no 3′ phosphatase activity. The model retains the intended ligation outcome; the exact dephosphorylation reagent remains unresolved.
- The article links its analysis software at [MARIO tools](https://mariotools.ucsd.edu/). The portal identifies v0.3 as the publication implementation; the linked legacy documentation could not be retrieved during this pass. Source coverage is therefore incomplete.
- The RT-primer/circularization procedure cites Huppertz et al., iCLIP, [10.1016/j.ymeth.2013.10.011](https://doi.org/10.1016/j.ymeth.2013.10.011), and terminal-linker removal cites T7-exonuclease and phosphorothioate-protection originals. Those linked method sources still need explicit independent inspection before accepting any missing-information conclusion.

## 2026-10-17 bounded pass (round breadth-first-12)

Read: PMC4931010 JATS Methods lines 116–261 (cross-linking through the oligo list) and
Supplementary Information Fig. 1 legend (`supp…-s1.pdf.txt` lines 3–17) and Fig. 26 legend
(lines 227–236). Scraper hits for every printed oligo were checked against the oligo list.

- 🟡 **Model correction — reaction order.** The previous graph selected linker-containing
  chimeras on streptavidin immediately after linker ligation and then went straight to
  RT/PCR. In the source, RNPs are first biotinylated on protein cysteines (IPB) and
  immobilized on MyOne T1 beads; PNK (3′-phosphatase minus) 5′ phosphorylation, linker
  ligation (T4 Rnl1), 3′ dephosphorylation, linker 5′ phosphorylation and dilute proximity
  ligation follow on beads; terminal linkers are degraded by guard-oligo/T7 exonuclease
  treatment; RNase III shearing, pre-adenylated RT-adapter ligation (Rnl2 tr KQ) and
  SuperScript III RT precede the linker-biotin selection on MyOne C1 beads with RNase H
  elution; CircLigase II circularization, Cut_oligo/BamHI relinearization, DP5/DP3
  pre-amplification and PE 1.0/2.0 PCR complete the library. All are now graph steps.
- 🟡 **Model correction — final library.** The page previously used a generic TruSeq shell
  with an i7 index and the TruSeq Read 2 adapter. It is now computed from the printed
  oligos (ES-1 sample): PE 1.0 (P5 + TruSeq Read 1), the reverse
  complement of the ES-1 inline barcode (NNN·ACCT·NN), RNA1 fragment, linker copy, RNA2
  fragment, the RT-adapter copy and the reverse complement of PE 2.0. The earlier graph labelled DP5/DP3 and PE primer landings exact; the current migration below corrects their actual incoming strand states.
- The two published discrepancies (RT-primer random-base count; phosphatase-minus PNK named
  for 3′ dephosphorylation) are now graph-owned annotations on the page.
- Deferred: iCLIP RT-primer/circularization source (Huppertz 2014, ref. 62), T7-exonuclease
  and phosphorothioate originals (refs 63–64), DSN and GeneRead rRNA-depletion branches,
  ES-indirect EGS/formaldehyde branch, MARIO tools v0.3 documentation, read lengths.


## Native graph migration — 2026-10-18, batch02

Re-read published PMC4931010 JATS Methods116–261 and the published supplement
Figure1 legend1–17, with the mandatory source-bundle scraper. Whole modified-oligo
transcription checks verify all11 printed inventory forms, including internal
biotin, phosphorothioates, /5Phos/, /5rApp/, /3ddC/ and /3InvdT/.
No cited iCLIP, exonuclease-protection, vendor or software original was newly read.

The module now declares its own metadata, native graph, sequencing plan and
reviewed-target list. Protein-bound RNA partners and the later RNA/cDNA hybrid
are independent strands. The complete RNA linker is a molecular input; the RT
adapter's terminal ddC remains on the RNA template through RT and is not copied
as part of the primer's core-binding interval. Proteinase-derived peptides, RNA, terminal-linker/guard and DSN
hydrolysate fates are unresolved where actual cleanup partition is not modelled.

The selected ES-1 route now includes the published DSN preparation: pre-amplify
until80–100ng, skip the separate150–350bp gel route, reanneal and treat with DSN,
then amplify the surviving population with full PE primers. Both PCR geometries
use actual ss or post-DSN precursors. Exact duplex first-cycle geometry is not
claimed. The DP5-derived pre-amplification arm excludes the terminal BamHI-site
G beyond the DP5 binding interval, a repository modelling correction.
The final1.8×SPRI/250–550bp selection is explicit. Other sample/depletion branches
remain separate variants.

A scoped circle-order relation compares the same cDNA parts before and after
neutral ssDNA closure; unknown insert extent makes this conditional. It does not
validate circularization yield, BamHI reopening or the whole reaction route.
Feature identities persist on RT primer, circular/reopened cDNA and the expected
target amplicon. The post-DSN population has unassigned native ends and pairing.
The source-positive paired biological readout is preserved, including unplaced
Read2, without inventing read lengths or an external index acquisition recipe.

Both published discrepancies remain on the public page: the explicit ES-1
five-random-base primer versus generic/supplementary descriptions, and the
phosphatase-minus enzyme named for intended3′dephosphorylation. The model uses
the printed ES-1 primer and labels the intended phosphate-removal outcome inferred.
This pass does not accept any missing-information claim. Deferred targets include
linked predecessor/vendor/software sources, other sample branches, actual
3′phosphatase identity, read lengths and native PCR/restriction/recovery proof.
