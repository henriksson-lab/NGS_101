# PARIS — research notes

Defining source: Lu *et al.*, [10.1016/j.cell.2016.04.028](https://doi.org/10.1016/j.cell.2016.04.028).

## Evidence used

- 🟢 AMT and long-wave UV reversibly crosslink base-paired RNA regions in living cells.
- 🟢 Partial digestion is followed by two-dimensional gel purification of crosslinked duplexes and proximity ligation.
- 🟢 Short-wave UV reverses the psoralen crosslink; the covalent RNA ligation junction remains for cDNA library preparation.
- 🟢 The authors' detailed PARIS protocol (Lu, Gong and Zhang, 2018, [10.1007/978-1-4939-7213-5_4](https://doi.org/10.1007/978-1-4939-7213-5_4)) prints the preadenylated/3′-biotin-blocked RNA adapter, phosphorylated barcoded RT-primer-1 with two internal iSp18 spacers, P3Tall/P6Tall primers, P3/P6 Solexa primers and P6 custom sequencing primer.
- 🟢 RT-primer-1 defines 24 six-base sample indices and carries 3 W plus 3 N bases before the index and 5 N bases after it for duplicate discrimination.
- 🟡 The detailed protocol says its library-preparation steps were optimized after the original Cell protocol. These oligos document that later PARIS workflow. The original 2016 Extended Experimental Procedures have now been recovered and distinguish two original library routes; neither should silently inherit the later RT-primer design (details below).

## Protocol boundary

SPLASH and LIGR-seq use the same broad psoralen–ligation principle but different duplex-enrichment chemistry, so they have separate pages.

## Original 2016 library routes (source audit, 2026-10-16)

🟢 The publisher's original `mmc1.pdf`, Extended Experimental Procedures,
PARIS Step 3 (extracted lines 92–110), distinguishes **HeLa random-hexamer RT**
from **HEK293/mES adapter ligation followed by adapter-primed RT**.
The HeLa ordering oligo is printed as:

```
/5phos/DDDNNXXXXNNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN
```

D means A/G/T; XXXX is a four-base sample barcode, whose ordering is inherited
from Spitale 2015. The first N6 tract identifies duplicates; terminal N6 primes RNA.
This differs from the later 2018 W3N3/six-base-index/N5 adapter-RT design.
These are source-labelled implementations, not an apparent publication error.

🟢 The HEK293/mES route expressly inherits Flynn 2016 icSHAPE **step 26 onward**
([doi:10.1038/nprot.2016.011](https://doi.org/10.1038/nprot.2016.011)).
Its original reagent list prints a distinct adapter-annealing RT-primer-1 example:

```
/5phos/DDDNNAACCNNNNAGATCGGAAGAGCGTCGTGGA/iSp18/GGATCC/iSp18/TACTGAACCGC
```

AACC is the printed example index. Flynn's steps 58–61 circularize cDNA with
CircLigase II before library PCR. The main graph now follows this 2016 adapter route
and includes adapter ligation and cDNA circularization; exact molecular end/junction
geometry still requires reconstruction. The original supplement states >40-nt inserts
were retained and libraries were sequenced on NextSeq. It prefers direct adapter
ligation for subsequent libraries because of greater RT efficiency.

The six original supplementary workbooks contain sequencing summary, duplex groups,
conservation/comparison coordinates, alternative structures and XIST structural results.
Their sequence-like entries describe observed RNA structures, rather than synthesized
ordering oligos; their scope is recorded separately in the source checklist.

## Historical four-base barcode evidence

🟢 The linked Spitale 2015 original ([doi:10.1038/nature14263](https://doi.org/10.1038/nature14263),
Methods, RT/circularization/library sections) prints a phosphorylated RT primer
with `DDDNN-AACC-NNNN`, two internal iSp18 spacers and the body
`AGATCGGAAGAGCGTCGTGAT`. Its Methods specify a 1×100-base HiSeq2500 run.
The Flynn 2016 example above instead prints `AGATCGGAAGAGCGTCGTGGA`;
these are kept as source-specific versions, without calling the difference a
publication error or replacing the PARIS HeLa primer with either example.

🟢 The article's author-linked icSHAPE repository was inspected at historical
revision `175bd97170287d8a34901c85e4e85eaf38754271` (2015-10-01).
The README and `splitFastq.pl` require caller-supplied barcode/library mappings
and position/length values; `5:4` is the splitting example. The configuration's
13-base leading trim is compatible with Spitale's DDDNN + four-base index + N4
scaffold. It does not establish a trim or ordering for the different PARIS HeLa
N6-containing scaffold.

🟡 AACC is a verified example, while a complete historical barcode ordering has
not been established by these readings. The original PARIS four-base placeholder
can remain unassigned. No conclusion that the full list is absent is accepted:
the historical member mapping has not been assigned to this model. No specific
inaccessible original has been identified as a barcode-order blocker. Supplementary Tables 1/2 describe RBP motifs and
miRNA-seed reactivity results, rather than synthesized library barcodes.

The four original Spitale publisher figure slides were subsequently rendered and
read visually. Figure 1 contains the NAI-N3/click-biotin/RT/capture workflow;
Figures 2–4 contain native RNA reactivity, motif, translation and m6A/Rbfox
results. Their biological motif labels are not library index members. These
figures do not change the distinction between a printed example, an unassigned
model index and a claim of absent published information. Exact circular-junction
and read-entry reconstruction remains work on the model.

## Bounded pass 2026-10-17 (round breadth-first-15)

🟢 Reaction conditions now named in the graph, from the original `mmc1.pdf.txt`
Steps 1–3 (lines 25–110): AMT 0.5 mg/ml with UV365 nm for 30 min on ice; S1 nuclease
in the urea/SDS lysate, proteinase K (50 °C, 30 min) and TRIzol; ShortCut RNase III
(MnCl2, 37 °C, 20 min); 12% native then 20% urea 2D gel with the region above the
diagonal excised; T4 RNA ligase 1 (NEB M0437M) for 16–20 h at room temperature,
stopped by boiling for 2 min; and UV254 nm for 15 min on ice.

🟢 The inherited Flynn 2016 steps 26–61 (PMC4896316 lines 532–800) supply the
adapter-route conditions. T4 RNA ligase 1 (high concentration) with 15% PEG8000 ligates
the preadenylylated 3′-blocked linker `/5rApp/AGATCGGAAGAGCGGTTCAG/` at 25 °C for 3 h.
Urea PAGE then size-selects the product; PARIS says >40-nt inserts were kept. RT uses
SuperScript III, followed by denaturing PAGE to separate cDNA from unextended RT primer.
CircLigase II with MnCl2 runs at 60 °C for 2 h, and qPCR-monitored Phusion HF PCR uses
the P5/P3-Solexa primers. The graph now includes the two gel size selections and their
discarded linker/primer waste.

🟡 Open: Flynn uses a 3′-biotin linker for DMSO samples and a 3′-ddC linker for NAI-N3
samples. It then performs streptavidin capture and RNase/biotin elution in steps 43–50.
The PARIS supplement does not say which linker variant it used. Whether that capture
applies is therefore not modelled, and the linker's 3′ block is left unassigned. The
read length is not given in Step 3; only NextSeq is stated.
