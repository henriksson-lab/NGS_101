# Abseq — droplet linkage of antibody and cell identities

Shahi et al. (2017), *Abseq: Ultrahigh-throughput single cell protein profiling with droplet microfluidic barcoding*, [10.1038/srep44447](https://doi.org/10.1038/srep44447).

## Bounded originals inspected

The defining published article was read in the external source store as
`srep44447_PMC5349531.xml.txt`: Results, Optimization of Abseq molecular biology;
Microfluidic workflow for Abseq; Bioinformatic workflow; Methods, Antibody
conjugation and purification, Cell staining and in-droplet cell lysis,
Droplet SOE-PCR and Library preparation/NGS, and Bioinformatics. Figure 1 was
also inspected in the rendered original PDF, page 2. The source-bundle scraper
was run; sequence hits were checked against the Methods. This pass read one
original, not all cited methods or an exhaustive absence audit.

## Molecular workflow

🟢 Original Abseq labels antibodies with DNA through a 5′ thiol and SM(PEG)6
maleimide/NHS chemistry, then purifies conjugates. Washed stained cells are
encapsulated, lysed with proteinase K, and heated to inactivate the protease.
A separate dilute-barcode emulsion is PCR-amplified to give clonal barcode drops.
One lysate drop, one barcode drop and PCR reagents are merged, mixed and divided
into four portions before SOE-PCR. The overlap extension physically joins a
cell barcode to antibody identity/UMI material. After breaking the emulsion,
column recovery and Exonuclease I digestion of unlinked ssDNA, bulk PCR supplies
sequencing input. The article specifies MiSeq v2 single-read 50-cycle sequencing.
This is not a poly(A)-captured antibody RT branch or a current BD AbSeq kit recipe.

## Source sequences and current representation

🟢 Methods prints P5 `AATGATACGGCGACCACCGA` (the source's shortened P5 primer),
P7 `CAAGCAGAAGACGGCATACGAGAT`, FAM probe 3Taq1
`AGATCGGAAGAGCGTCGTGTAGG`, and a bioinformatics-filter homology sequence
`GCAGTGGTATCAACGCAGAG`. The model imports the canonical Illumina components,
retaining the source's P5 truncation; the probe and filtering motif are different
reported sequences, not silently treated as the same overlap.

🟢 Methods identifies CD3 by `TTATAAC` and CD19 by `TTAATTG`; Results assigns an
eight-random-base UMI to each antibody tag. These are represented as separate
identity-bearing components, with structured feature metadata. Their full strand
positions, terminal linkers and complete input oligos remain unassigned.
The reported substrate/product lengths alone do not reconstruct those molecules.

## Published disagreement

The Results describes random **10-mers** for barcode generation and **10-bp**
cell barcodes. Methods Bioinformatics prints `NNNN NNNN`, eight characters.
The model assigns no exact cell-barcode width. Whether this is abbreviated
notation or a different construct is unresolved; no silent correction is made.
The preview displays both source locations and the model's choice.

## Model limits and next sources

The graph starts with cells and already purified antibody–DNA conjugates. It
tracks washing, protease treatment/inactivation, clonal barcode PCR, droplet
merger, SOE-PCR, recovery, Exonuclease I and bulk PCR. PCR geometry is explicitly
unavailable from this bounded reading; no exact overlap/product relation or
per-cycle read map is claimed. The primer-free and oligo-free graph nodes are not
proof that the original lacks that information. Catalogue status remains notes.

Deferred: source-linked droplet barcoding original (Lan et al., 2016,
10.1038/ncomms11784); cited protease-droplet methods (Eastburn et al.,
10.1093/nar/gku606 and 10.1021/ac402057q; identities/relevance still to verify);
complete barcode/tag design and sequencing primer/cycle mapping; complete source
inventory and independent reading. The fetched PMC supplementary archive contains
article figures; its archive member names were inspected, not accepted as an
independent oligo source. No sequence or condition is confirmed absent.

## Lineage

The article says droplet barcoding was recently described for DNA molecules and
transcriptomes (references 26, 29–30), and describes the protease workflow as
previously shown (references 33–34). These are component-level predecessor
candidates, pending reading of the cited originals. Similarity to CITE-seq/REAP-seq
is not evidence of ancestry. No parent molecular constructor is imported.
