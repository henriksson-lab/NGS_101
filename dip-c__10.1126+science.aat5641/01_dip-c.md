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
