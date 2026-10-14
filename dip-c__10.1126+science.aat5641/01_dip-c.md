# Dip-C

## 1. What it is

A single-cell chromatin-conformation method that combines proximity ligation with multiplex end-tagging amplification (META) to obtain dense diploid contact maps.

## 2. Sources and evidence

- 🟢 Tan et al., *Science* 2018, DOI [10.1126/science.aat5641](https://doi.org/10.1126/science.aat5641), defining paper.
- 🟢 Tan et al., *Nature Protocols* 2021, [PMC8225968](https://pmc.ncbi.nlm.nih.gov/articles/PMC8225968/), detailed Dip-C/META procedure.
- 🟢 The workflow omits biotin fill-in and streptavidin enrichment.
- 🟢 The detailed protocol permits MboI, DpnII, NlaIII, or combinations; the schematic draws its default MboI example rather than claiming one universal Dip-C cut site.
- 🟡 Canonical outer sequencing arms are a kit-level reconstruction.

## 3. Distinguishing chemistry

After contact ligation and single-nucleus isolation, META introduces transposon end tags and amplifies genomic material. That branch is neither PCR-only classical scHi-C nor phi29 MDA-based snHi-C.

🟢 The detailed protocol notes that each read begins with 39 bases of META transposon sequence before usable genomic sequence; the page keeps these regions explicit.

## 4. Oligo evidence still needed

🔴 The accessible defining-paper text and 2021 protocol identify the META read-prefix
length, but do not print the complete transposon oligos, their modifications, or the
outer sequencing shell. The defining paper points to Materials and Methods and Tables
S1–S2 in its supplemental package; those files are listed by PMC but require manual
retrieval in this source workflow. The same paper deposits raw reads as **SRP149125**
and processed data as **GSE117876**. A consensus of read prefixes from the raw reads
could support a labelled inference about sequenced bases, after verifying sample and
run provenance; it cannot establish the ordered oligos or chemical modifications.
An inventor-associated [META patent](https://patents.google.com/patent/WO2018217912A1/en)
(Example VII and Example XI) publishes a twenty-member family of transposon priming-site
sequences and explicitly discusses Dip-C. 🟡 These are candidate source sequences, not
proof of which members or complete transposomes were used for the 2018 libraries; no
patent sequence is promoted into the schematic without that cross-check.
