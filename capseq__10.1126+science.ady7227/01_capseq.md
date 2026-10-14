# CapSeq — capsule-based single-cell RNA sequencing

> **Evidence marks:** 🟢 stated in a source examined here · 🟡 interpretation or
> derived consequence · 🔴 unresolved. This is a reference note, not a finished
> base-by-base schematic: the barcode oligo tables cited by the methods still need review.

**Defining paper:** Baronas D, Norvaisis S, Zvirblyte J, et al., including Mazutis L.
"High-throughput single-cell omics using semipermeable capsules." *Science* (online
2025; issue 2026), doi:[10.1126/science.ady7227](https://doi.org/10.1126/science.ady7227).
The author's [2025 bioRxiv preprint](https://doi.org/10.1101/2025.03.14.642805)
explicitly names the RNA workflow **CapSeq (Capsule-based RNA-Sequencing)**. 🟢

## Sources examined and access boundary

| Source | What could be checked | Limit |
|---|---|---|
| [Science paper](https://doi.org/10.1126/science.ady7227) and [PubMed record](https://pubmed.ncbi.nlm.nih.gov/41411411/) | identity, publication, authors, abstract | full Science text and supplement were not accessible to the source fetcher |
| [bioRxiv preprint, version 1](https://www.biorxiv.org/content/10.1101/2025.03.14.642805v1.full) and [PMC mirror](https://pmc.ncbi.nlm.nih.gov/articles/PMC11957016/) | full preprint, Supplement 1 methods and Supplement 2 figures, fetched via `tools/get_sources.py` into `$CHEM_DATA/sources/` | methods cite Tables S21–S25, but these oligo tables are absent from the two archived supplementary PDFs |
| [Authors' analysis repository](https://github.com/mazutislab/capseq) | repository scope and link to the STARsolo wrapper | analysis code is not an oligo-ordering source |
| [GEO GSE291189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291189), especially [GSM8830586](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8830586) | submitter's library construction description, two adapter oligos, read processing lengths | does not give Table S25 or the complete upstream R1–R3 oligo designs; no raw reads were examined |
| [WO2026190526A1](https://eureka.patsnap.com/patent/WO2026190526A1), Example 3 and Tables 1–5 | later inventors' disclosure names CapSeq-REV/FWD, TSO, hairpin, BC1–BC3 plates and downstream adapter; its adapter pair and processed layout agree with GEO | patent text extraction merges some table rows and misrenders `µL` as `pL`; each sequence needs checking against the patent PDF and the Science supplement before transcription into a model |

`tools/get_sources.py --slugs capseq__10.1126+science.ady7227` retrieved the
preprint, both supplementary PDFs and PMC XML. The Science publisher text remains
`(manual)` in the source manifest. 🔴

## What is new

🟢 The biological partition is a **semipermeable capsule** with a liquid core. Cells
are isolated, lysed and washed in separate capsules, retaining individual RNA while
diffusible inhibitors can be removed. Barcoding then uses plates for successive
split-and-pool rounds. The capsule is a reaction compartment, not a sequence segment.

🟢 The preprint's Results section gives this order:

1. Isolate one cell per semipermeable capsule; lyse and wash to obtain an individual
   transcriptome in each capsule.
2. Distribute capsules into a 96-well plate (approximately 1000 capsules per well).
   Reverse transcribe poly(A) RNA with a poly(dT) primer carrying barcode round **R1**
   and a **UMI**.
3. Pool and redistribute capsules into a second plate. Ligate barcode round **R2** to
   cDNA.
4. Repeat the split-and-pool ligation for barcode round **R3**.
5. Pool and divide capsules into sublibraries. Add sequencing barcode round **R4** by
   indexing PCR, then sequence and group reads by the R1–R4 combination.
6. 🟢 In the [GSM8830586 submitter protocol](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8830586), amplified cDNA is fragmented and end-repaired/dA-tailed with NEBNext Ultra II FS, ligated to an adapter duplex, then amplified with two indexing primers (the submitter refers to **Table S25**). This is a downstream library-preparation stage after the capsule barcoding. Its order and exact adapter sequence must be represented before a final construct can be drawn.

🟡 The four rounds identify a cell by its route through wells and sublibraries;
the UMI identifies original RNA molecules within that cell. The paper reports over
10 million combinations and an approximately 300,000-cell run at a 1% multiplet
rate. Those capacity figures do not define the base lengths of the barcodes.

🟢 Supplement 1, "CapSeq procedure" (PDF text lines 276–340), confirms the three
barcoding reactions, then RNase H/A/T1 digestion of the RNA:cDNA hybrid, overnight
Hairpin-TSO ligation with T4 DNA ligase, and 15-cycle cDNA PCR *inside capsules*
before capsule dissolution. The later patent and GEO also discuss a no-hairpin
comparison; keep this variant separate until its full oligos are verified. 🔴

## Oligos and molecular geometry still to establish

| Component | Established here | Still required for a molecular model |
|---|---|---|
| R1 reverse-transcription primer | 🟢 poly(dT), first barcode and UMI are present; GEO processing reports BC1 8 nt and UMI 8 nt in the *processed read* | 🔴 full 5′→3′ oligo including modifications, spacer/handle and anchored versus unanchored dT |
| R2 and R3 ligation oligos | 🟢 attached in two successive ligation rounds; GEO processing reports BC2 10 nt and BC3 8 nt in the *processed read* | 🔴 exact sequences and modifications, bridge/splint design, end polarity and sealed junctions |
| Post-amplification library adapter | 🟢 submitter gives two modified oligos verbatim below | 🔴 source-confirmed relationship of this adapter to the earlier barcoded cDNA, including exact junction after end repair |
| Template-switching oligo, hairpin and cDNA PCR primers | 🟢 preprint Supplement 1 prints the TSO verbatim below and names Hairpin-TSO, CapSeq-REV and CapSeq-FWD | 🔴 verify hairpin and PCR primer sequences in the cited Tables S21–S25 or the patent PDF |
| R4 indexing PCR primers | 🟢 sublibrary barcode introduced by PCR; GEO points to Table S25 | 🔴 primer sequences and index position |
| Final library and sequencing | 🟢 paired-end NextSeq 2000 data in GEO; processed Read 1 carries BC3–BC2–BC1–UMI after linker removal | 🔴 complete unprocessed read structure, primer binding sites and exact cycle layout |

🟢 The [GSM8830586 construction protocol](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8830586)
gives this adapter pair in IDT notation, both HPLC purified:

```
forward ligation primer  /5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC/3ddC
reverse ligation primer  /5AmMC6/GCTCTTCCGATCT
```

The submitter says the pair is annealed before ligation to fragmented, repaired,
dA-tailed cDNA. The 5′ phosphate, 3′ dideoxy C and 5′ amino modification are part of
the transcribed oligos, not plain DNA bases. 🟢 The same GEO entry says barcode linker
sequences were removed from Read 1 before the processed layout
`BC3[8 nt]–BC2[10 nt]–BC1[8 nt]–UMI[8 nt]` was supplied to STARsolo. That is an
analysis layout after trimming; it is not an exact unprocessed molecule or run recipe.

🟢 Supplement 1 (PDF text line 282) prints the template-switching oligo as
`/5Biosg/AAGCAGTGGTATCAACGCAGAGTACATrGrGrG`. Its 5′ biotin and three
3′ RNA guanosines are source-declared modifications. The same supplement refers
to the 96-well BC1, BC2 and BC3 oligos in Tables S22–S24 and the Hairpin-TSO in
Table S21; those tables are not present in its two PDF attachments. 🔴

Do not substitute the ligation geometry or barcode sequences of SPLiT-seq, sci-RNA-seq,
or the Atrandi genome barcoding kit. Their similar split-and-pool topology does not
establish the CapSeq oligos. 🟡 The authors' analysis pipeline or deposited FASTQ may
help check read-cycle boundaries, but observed reads alone do not establish an ordering
oligo's modifications, ligation ends or all invariant bases.

The [later CapMux paper](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1846065/full)
describes a **split** UMI (U1/U2), while the older GEO processing description combines
it into one 8-nt UMI after linker removal. 🟡 These can be consistent views of one
construct, but that relationship must be checked against the unprocessed read and
the primary oligo table. The newer GEO sample
[GSM9578506](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9578506)
reports raw Read 1 starts and lengths for BC3, BC2, BC1 and the split UMI; it is a
useful independent check, not a substitute for the synthesis oligos.

## Next source review

Obtain the Science supplementary oligo tables or the patent PDF; read the detailed
oligo tables and record
every oligo verbatim with its source location. Then construct the RT and ligation intermediates with `lib/`
semantics, attach structured features to R1–R4 and the UMI, and derive the sequencing
layout from the final construct. Until that review, this entry stays `notes` and its
title has no public schematic link.

## Name collision

This **CapSeq** is distinct from the repository's **CAP-seq**
(`cap-seq__10.1371+journal.pgen.1001134`), which uses CXXC affinity capture of
unmethylated CpG-rich DNA. The similar names do not imply related chemistry.
