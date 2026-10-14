# Digenome-seq — research notes

Kim et al., *Nature Methods* (2015), [doi:10.1038/nmeth.3284](https://doi.org/10.1038/nmeth.3284).

## Assay target

The sequenced material is purified human genomic DNA cleaved in vitro by Cas9, sequenced as an ordinary whole genome. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a stock TruSeq whole-genome library kit, with no cut-end-specific oligo). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining article's supplementary text/figures (`41592_2015_BFnmeth3284_MOESM426_ESM.pdf`), including the 5′-end plots, Supplementary Notes and Supplementary Tables 1–6. The **main article is paywalled** and is not in PMC; `get_sources.py` records it as a `(manual)` row. Nothing in the fetched supplement mentions Covaris, shearing, ligation, RNase, TruSeq, adapters, indexes or PCR-free preparation — it is tables and figure legends only.
- 🟢 Kim et al., *Genome Research* (2016), [doi:10.1101/gr.199588.115](https://doi.org/10.1101/gr.199588.115) (PMC4772022, open access), the same group's multiplex Digenome-seq paper, read verbatim under "Whole-genome and Digenome sequencing": RNase A at 50 µg/mL, Covaris fragmentation to 400–500 bp, 1 µg of fragmented DNA "ligated with adaptors using TruSeq DNA library prep kit", Illumina HiSeq X Ten, 30–40× depth, and adapter clipping with `AGATCGGAAGAGC`/`GCTCTTCCGATCT`. Registered in `catalogue/extra_sources.tsv`.
- 🟢 Illumina *TruSeq DNA PCR-Free Reference Guide* 1000000039279 v00 (October 2017) and *Illumina Adapter Sequences* 1000000002694 v22, both registered in `catalogue/extra_sources.tsv` and fetchable. The adapter document states for TruSeq single indexes: "Index adapter sequences are six bases as underlined."
- 🟡 Kim, Kang and Kim, *Nature Protocols* (2021), [doi:10.1038/s41596-020-00453-6](https://doi.org/10.1038/s41596-020-00453-6). **Not read.** Nature Protocols is paywalled and the article has no PMCID, so there is no legitimate fetch route; it is deliberately *not* in `extra_sources.tsv`. Claims attributed to it below are marked 🟡 until someone obtains the PDF by hand and drops it into `$CHEM_DATA/sources/digenome-seq__10.1038+nmeth.3284/`.

## Construct logic

1. 🟢 Purified genomic DNA is cleaved in vitro with Cas9–sgRNA. Digenome-seq therefore measures biochemical cleavage of naked DNA, not editing inside a cell.
2. 🟢 RNase A (50 µg/mL) removes guide RNA, protein is removed, and the DNA is repurified with the DNeasy kit. No barcode or affinity handle is attached to a cut end. (*Genome Research* 2016 methods, verbatim.)
3. 🟢 The entire digested genome is fragmented to 400–500 bp with Covaris. The resulting pool contains a small class of fragments bounded by a Cas9 cut and a random shear, among the much larger population bounded only by random shearing. (*Genome Research* 2016 methods, verbatim.)
4. 🟢 1 µg of fragmented DNA is "ligated with adaptors using TruSeq DNA library prep kit" and sequenced on an Illumina HiSeq X Ten at 30–40× depth; the authors' read clipping uses the TruSeq adapter stems `AGATCGGAAGAGC`/`GCTCTTCCGATCT`. (*Genome Research* 2016 methods, verbatim.)
   🟡 *Which* TruSeq DNA kit is not stated, and **no readable source says "PCR-free"**. "TruSeq DNA library prep kit" does not distinguish the PCR-Free product from the amplified TruSeq DNA Sample Prep / TruSeq Nano kits, the fetched 2015 supplement never mentions library preparation at all, and this is a 2016 companion paper rather than the defining one. The page's "TruSeq PCR-free" framing is therefore an inference, not a transcription — see "Unsubstantiated on the page" below.
5. 🟢 After alignment, a cleavage site produces forward- and reverse-strand reads sharing precise 5′ endpoints. Random shearing distributes endpoints, so it does not make the same local pileup.

## What is—and is not—encoded in the molecule

There is no Digenome-specific oligo, UMI, cut-site barcode, or selective adapter. The standard sample index identifies a library, but it does not identify a cut end. The assay signal is the genomic coordinate of the physical insert boundary, interpreted across many ordinary WGS molecules.

The page consequently shows a standard single-index TruSeq construct and declares the genomic insert boundary as the measurement. It does not add a fictional proximal adapter or enrichment step.

## What the targeted-deep-sequencing primers do and do not establish

🟢 Supplementary Tables 4 and 6 ("List of primers") are the validation assay, not the WGS
library. They are 1st/2nd-PCR primer pairs for HBB, VEGF-A and the Digenome-captured
candidate sites, and their tails are verbatim `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` (163
occurrences) and `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` (165). Those independently confirm
`il.TRUSEQ_READ1` and `revcomp(il.INDEX1_PRIMER)+T` as this group's amplicon tails, but
they say nothing about the whole-genome library's adapters, index length or whether that
library was amplified. The model is right to claim no protocol-specific oligo.

## Unsubstantiated on the page (open)

Three claims on `digenome-seq.html` have no backing in any source that can be read here.
They are recorded rather than silently repaired because the Digenome-seq model lives in
`lib/dsb_capture.py`, which this task was instructed not to edit.

1. 🟡 **"PCR-free" is derived, not transcribed — but it is supported.** The authors state
   the absence of amplification directly: "because no PCR amplification steps are involved
   with Digenome-seq prior to WGS, it is much simpler and easier to carry out than other
   methods" (Genome Research 2016 Discussion,
   `PMC4772022_PMC4772022.xml.txt:114`). That is a method-comparison claim in a companion
   paper — its contrast is GUIDE-seq, "where PCR is used to amplify
   oligonucleotide-captured sites" — and it names no kit, so it establishes *that no PCR
   precedes sequencing*, not *which* product was used. The only library sentence available,
   "TruSeq DNA library prep kit", does not distinguish the PCR-Free catalogue number from
   the amplified TruSeq DNA Sample Prep and Nano kits. So `amplification=PCR-free` in
   `catalogue/properties.tsv` stands on the authors' own statement, while the adapter
   *construction* stays inferred: nothing here gives one adapter base, and this page does
   not use the full-length-adaptor builder (`truseq_pcr_free_library`) that models a real
   PCR-Free kit, because that would assert the pre-indexed two-arm adaptor as this
   protocol's chemistry. Reading the 2021 *Nature Protocols* article would settle the
   product; it is paywalled with no legitimate route and is deliberately not registered.
2. 🟡 **Source-backed adapter styling.** `_digenome_library()` calls `truseq_library(…,
   dual_index=False)` with the default `inferred_adapters=False`, so P5, the Read 1 arm,
   the dA junction, the Index 1/Read 2 arm, i7 and P7_rc all render as transcribed from a
   source. What is actually published is a kit *family name* in a companion paper and two
   adapter clipping stems; no source gives this library's adapter oligos, its index format
   or even which TruSeq kit was used. Per AGENTS.md ("uncertainty is a constructor
   property") this should be
   `inferred_adapters=True`, which also supplies the `<inf>` base styling that the page
   currently lacks; the step caption and page caveat have to say INFERRED alongside it.
3. 🟡 **An 8-base i7 index.** `truseq_library` defaults to `i7_nt=8` and `READ_LENGTHS`
   declares an 8-cycle Index 1 read. Illumina's *Adapter Sequences* document
   (1000000002694 v22) states that TruSeq **single** index adapters are six bases
   (AD001–AD027, `ATCACG`, `CGATGT`, …); eight-base indexes arrived with the dual-index
   IDT-for-Illumina TruSeq DNA UD revision, which postdates the 2014–15 HiSeq X Ten work.
   A single-index TruSeq library with an 8-nt i7 is a product Illumina never sold, so the
   drawn geometry is internally inconsistent. The defensible model is `i7_nt=6` with
   `inferred_adapters=True` (and `READ_LENGTHS["Index 1 (i7)"] = 6`) — the same treatment
   DISCOVER-Seq already gets in `lib/editing_outcomes.py` for a historical single-index
   TruSeq kit. No source on disk states the length, so it stays an inference either way.
