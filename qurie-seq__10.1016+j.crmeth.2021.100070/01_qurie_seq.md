# QuRIE-seq — intracellular epitopes and RNA through separate library branches

Rivello et al. (2021), *Single-cell intracellular epitope and transcript detection reveals signal transduction dynamics*, [10.1016/j.crmeth.2021.100070](https://doi.org/10.1016/j.crmeth.2021.100070).

## Bounded source scope

Two published originals were inspected: `j.crmeth.2021.100070_PMC9017125.xml.txt`,
Results overview and STAR Methods antibody labeling, BJAB preparation,
microfluidic barcoding, post-RT recovery and both library-preparation sections
(lines 281–327); and published `j.crmeth.2021.100070_supp_mmc1.pdf.txt`,
Supplementary Table S3 (lines 1073–1125). The rendered Table S3, PDF page 14
(printed 13), was read to confirm sequence characters and columns. The source
bundle oligo scraper was run after retrieval; published hits were verified.
Downloaded superseded preprint material is excluded from active chemistry.
The downloaded mmc2 PDF is a combined article/supplement counterpart, not an
independent chemistry authority in this pass. No absence claim is made.

## Source-defined chemistry

🟢 Antibodies are labeled through DBCO-S-S-NHS ester and **5′ azide** DNA.
Cells are reversibly fixed with DSP/SPDP, quenched/permeabilized and stained.
They are co-encapsulated with inDrop-based primer-loaded beads and RT mix.
DTT reduces cell cross-links; UV releases bead primers. Poly(dT) primes mRNA
and the antibody tags' poly(dA), yielding barcoded products. Following emulsion
break, debris/oil recovery and Exonuclease I treatment, **both** products of the
0.6× AMPure fractionation are retained: long RNA-derived products on beads and
short protein-derived DNA reporters in the supernatant.

🟢 The RNA branch receives second-strand synthesis, T7 IVT, RNA fragmentation,
random-primer RT and PCR. The protein branch receives T7 IVT, a stored post-IVT
backup aliquot, sequence-specific RT and PCR. The source states **37°C for
15 hours** for RNA IVT and **50°C for 15 hours** for protein IVT. These are
preserved as branch-specific source conditions, not reconciled by a shared
recipe. Both branches use qPCR to select amplification cycles, undergo their
own cleanup and yield separately represented sequencing libraries.

## Published oligo discrepancy

🟢 The rendered Table S3 calls PE2-N6 a 38-bp oligo and prints the 32-base arm
`TCGGCATTCCTGCTGAACCGCTCTTCCGATCT` followed by **five** Ns. The main RNA
Methods calls the primer a random hexamer. The printed characters total 37,
not the stated 38. The local model preserves the five Ns as written and leaves
its effective priming geometry unresolved. The preview identifies both source
locations; it does not silently add a sixth N or claim which oligo was ordered.

🟢 PE1 contains canonical P7, a six-base library index and
`CTCTTTCCCTACACGA`. PE2 and the custom sequencing primers are transcribed from
Table S3 in the local module using canonical Illumina components where they
match. Methods additionally names PE2-NNNN-Next1 and PE2-NNNN-BioHash2 for the
protein branch. Their exact sequence/version mapping remains unassigned in this
bounded reading; names alone do not justify importing another protocol's recipe.

## Model scope and open work

The graph follows both RNA and protein branches and records cleanup waste,
qPCR aliquots and the explicitly stored post-IVT protein backup. Antibody-tag
identity, UMI/cell-barcode structure, complete bead-primer/promoter geometry,
full gene templates and per-cycle sequencing layouts remain unassigned.
Sequence-known individual primers do not establish an entire final library.
NextSeq500 output is source-supported; this pass does not invent read lengths
or exact priming sites. Catalogue status remains notes.

Deferred: remainder of article/supplement including full antibody dictionary;
exact bead-primer and protein RT oligo designs; full inDrop protocol version and
vendor source; complete transcription substrates/promoter geometry; source
reconciliation of PE2-N6; read-cycle declarations and final QC details;
independent source review. None is declared confirmed missing.

## Lineage

The published Results and microfluidic Methods explicitly state that the
platform is based on inDrop and cite Klein et al. (2015) and Zilionis et al.
(2017). This supports an inDrop → QuRIE component edge for droplet/bead capture,
subject to version/source qualification; it does not license importing all
inDrop molecular defaults. Gerlach et al. (2019), cited for reversible fixation
and intracellular profiling, is another component-lineage candidate whose
original and precise scope remain deferred. CITE/REAP are contextual comparisons.
