# SPLASH — research notes

Defining source: Aw *et al.*, [10.1016/j.molcel.2016.04.028](https://doi.org/10.1016/j.molcel.2016.04.028). Detailed protocol: [10.3791/55255](https://doi.org/10.3791/55255).

## Evidence used

- 🟢 SPLASH expands to Sequencing of Psoralen crosslinked, Ligated, and Selected Hybrids.
- 🟢 A biotinylated psoralen crosslinks paired RNA; streptavidin enrichment selects those duplexes before proximity ligation.
- 🟢 The crosslink is reversed after ligation and the joined RNA partners are converted into a sequencing library.
- 🟢 The author protocol prints the activated 3′ adapter as `/5rApp/CTGTAGGCACCATCAAT/3ddC/`.
- 🟢 It prints the 96-nt RT primer as `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC/iSp18/CACTCA/iSp18/TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG`. The 3′ terminal 17 nt are complementary to the adapter.
- 🟢 After RT, the source gel-selects a product 96 bases longer than the ligated RNA, circularizes the single-stranded cDNA, and amplifies with reagents called “Universal PCR Primer” and “Index (X) Primer”. It specifies paired-end 150-base sequencing.
- 🟢 JoVE section 5 treats crosslinked RNA ends with T4 PNK to convert 3′ cyclic phosphate to 3′-OH and then 5′-OH to 5′ phosphate before proximity ligation. After photoreversal, section 7 ligates the 5′-activated, 3′-blocked RNA linker, reverse-transcribes, and circularizes first-strand cDNA with CircLigase II before PCR. The graph now declares these operations in that order.
- 🟢 The JoVE materials list identifies NEB E7335L as its universal and indexed PCR primer kit. The [NEB E7335 manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7335_e7500_-e7710_e7730.pdf) prints the universal and Index 1 primer sequences, including each 3′ phosphorothioate-linked terminal T (Set 1 table). The schematic shows those exact reagents; Index 1 is one representative of the kit's indexed family.
- 🔴 The retrieved author PDF does not state the RT primer's 5′ end group or give a base-resolved circular cDNA junction. The named primer kit alone does not establish an exact circle-to-amplicon or read-primer geometry.

## Protocol boundary

The biotin–streptavidin enrichment is the defining distinction from the 2D-gel selection in PARIS and RNase R enrichment in LIGR-seq.

The two disclosed custom oligos and the vendor-defined PCR primers are rendered on the schematic page. Circle and amplicon geometry remain explicitly unresolved.
