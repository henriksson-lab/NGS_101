# scM&T-seq — research notes

## Current source scope

This note describes the source-qualified native graph migration, not a complete chemistry
or missing-information review. Historical source-review records remain in the catalogue.

- 🟢 Angermueller et al. 2016, [defining article](https://doi.org/10.1038/nmeth.3728),
  published PMC4770512 JATS, Online Methods lines 82–88 and the main-text primer
  contamination observation. This pass did not reread Supplementary Fig. 3e; earlier
  review records describe its first-strand-primer concatemer observation.
- 🟢 Smallwood et al. 2014, [cited scBS-seq predecessor](https://doi.org/10.1038/nmeth.3035),
  PMC4117646 Methods, Single-Cell Library Preparation. Inherited component scope is
  limited by the explicit scM&T changes below.
- 🟢 Macaulay et al. 2015, [cited G&T-seq predecessor](https://doi.org/10.1038/nmeth.3370),
  original publisher article PDF including Online Methods, cached as
  `zenodo895966_nmeth.3370_article.pdf.txt` in the G&T source bundle, printed pages
  13–14 / text lines 315–350. This is the cited 2015 article, not the later detailed
  G&T protocol revision.

The mandatory scM&T bundle oligo scraper ran: 19 hits in five hit files; five large
analysis files were skipped. These counts do not establish absence. This capped pass
inspected three originals. G&T S1, Smallwood's binary supplementary Word file,
period-specific Nextera XT documentation and linked index details remain deferred.

## Physical separation and RNA route

🟢 Cells enter RLT Plus with SUPERase-In. The cited G&T method uses magnetic
streptavidin beads and a 5′-biotin-triethyleneglycol capture primer:
`AAGCAGTGGTATCAACGCAGAGTAC` followed by `T30VN`. RNA and genomic DNA are distinct
physical populations; bead washes are retained with the DNA. Shared cell identity
is well provenance, not an invented covalent barcode.

🟢 The cited 2015 G&T original prints TSO
`AAGCAGTGGTATCAACGCAGAGTACrGrG+G` and ISPCR
`AAGCAGTGGTATCAACGCAGAGT`. Its RT and cDNA amplification reference uses SuperScript II,
KAPA HiFi and 1:1 AMPure recovery. 🟡 The model explicitly qualifies these as inherited
2015 reference components. The exact study-specific TSO assignment and complete copied
end geometry remain inferred; a later ACAT-containing reference does not override the
cited AC form. The public page retains the caveat, inferred caption and inferred linker.

The graph shows first-strand cDNA beside the independent RNA template, followed by a
switched single-strand precursor. That actual precursor enters one-species ISPCR;
a completed duplex is not substituted for its input. Native first-cycle geometry is
unassigned. Recovered cDNA enters the defining article's one-fifth-volume Nextera XT
route. Native transferred strands, gaps, selected index primers and opposite-end
architecture remain explicit boundaries, rather than default completed XT ends.

## DNA route and inherited changes

🟢 The defining article specifies DNA elution into 10 µl water, EZ Methylation Direct
conversion with AMPure XP present, and magnetic-bead purification/desulphonation into
first-strand mix. Physical converted DNA contains uracil before copying; a sequenced
T call is not substituted for the converted molecule.

🟢 Smallwood prints oligo1 handle `CTACACGACGCTCTTCCGATCT` and oligo2 handle
`TGCTGAACCGCTCTTCCGATCT`, each followed by N9 in its own method. Angermueller explicitly
changes **both primers to N6**, omits first-strand biotin capture and adds **0.8× AMPure
between second-strand synthesis and PCR**. Omission of capture does not establish
that the original 5′ biotin primer modification was removed; native 5′ chemistry
remains unassigned.

🟢 The cited scBS recipe supplies repeated Klenow exo-minus first-strand synthesis,
Exonuclease I excess-primer digestion, 0.8× bead recovery, oligo2 copying,
PE1.0/iPCRTag amplification and 0.8× post-PCR recovery. The graph retains originals
and independently copied strand classes. Exonuclease hydrolysis is chemical conversion,
not physical removal. No streptavidin/alkali capture is inherited into the modified
scM&T branch. Native complete copy geometry, selected cycle choice, indexed iPCRTag
strand/member and index width remain unassigned.

The actual oligo1/oligo2-derived mixture enters PCR. The source PE1.0 sequence is
computed from canonical Illumina components; no guessed eight-base iPCRTag scaffold
supplies the opposite end.

## Sequenced products and read layout

🟢 The defining article reports substantial primer contamination among sequenced DNA
and points to Supplementary Fig. 3. Earlier independent reads record a first-primer
concatemer observation there. Its assembly and formation stage remain unresolved.
The graph treats primer-derived non-genomic DNA as a **sequenced coexisting library
class**, including through final bead recovery. Computational failure to align is
not an experimental discard; bead cleanup does not establish complete removal.

🟢 RNA and DNA libraries are sequenced on HiSeq2000 with 125-bp paired-end reads.
The DNA analysis trims six random-priming bases and uses nondirectional alignment.
Index acquisition, cycles, selected members and native sequencing-primer landing
geometry remain unassigned. Unknown physical index widths are not filled from a
canonical kit or inferred from read cycles.

Open follow-up: the deferred originals named above; study-specific RNA oligo choices;
complete first-copy and native XT/PCR geometry; iPCRTag identity and acquired indexes;
primer-derived product assembly and physical class recovery. No publication-wide
absence claim is accepted.
