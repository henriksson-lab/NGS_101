# lentiCRISPR v2 screening

Sanjana, Shalem and Zhang introduced the higher-titre lentiCRISPR v2 all-in-one vector and the GeCKO v2 libraries. Its guide cassette uses the same BsmBI cloning logic as v1, but its cPPT position changes which genomic-readout primer pairs make a short product.

## Sources read

- 🟢 Sanjana, Shalem & Zhang 2014, *Nature Methods*, [doi:10.1038/nmeth.3047](https://doi.org/10.1038/nmeth.3047) — defining v2 vector and libraries.
- 🟢 Addgene plasmid #52961 and the Zhang-lab cloning protocol — vector map and BsmBI guide cloning.
- 🟢 Joung et al. 2017, [doi:10.1038/nprot.2017.016](https://doi.org/10.1038/nprot.2017.016) — compatible vector-independent one-PCR readout.

## Oligo scope

🟢 The Sanjana 2014 Supplementary Methods print the **ArrayF** and **ArrayR** primers
used to amplify the synthesized 74-nt GeCKO v2 guide-oligo pool (the PDF line-wrap
spaces are removed here):

```
ArrayF  TAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG
ArrayR  ACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC
```

Those are **pool-construction primers**, not the PCR used to count guides from screened
cells. The supplement says Illumina adapters were added for library-representation
validation, but does not print that validation primer pair. The page's later Joung readout
is sourced to its own protocol rather than assigned to the original v2 experiment.

## Important distinction

The guide scaffold is not the v1/v2 difference. The cPPT/CTS moves from downstream of the guide cassette in v1 to upstream in v2. Readout primers must therefore be chosen from the actual vector geometry.

## Defining pooled construction (source audit, 2026-10-16)

🟢 Sanjana 2014 Supplementary Methods, “GeCKO library pooled synthesis and cloning”
(lines 286–319 in the extracted publisher supplement), starts with synthesized **74-nt
oligos**, amplifies them with ArrayF/ArrayR into **140-bp inserts**, and inserts those
into a **BsmBI-digested, FastAP-treated vector using Gibson assembly**. The same
paragraph appears in the original preprint (lines 586–618). The page and graph now
show this pooled front end. The former single-guide duplex ligation display described
an individual-guide cloning route and omitted these defining pooled-library operations;
that omission was a repository modelling error.

The primer strands are source-transcribed, but the graph does not yet reconstruct an
individual synthesized 74-nt member and its exact Gibson junctions. Its assembly step
therefore carries an explicit unresolved molecular-transition warning. This warning
records the present model scope and does not assert that the source omits the information.

## Reaction conditions and the paper's own sequencing endpoint (bounded pass, 2026-10-17)

🟢 Re-read of the same supplement paragraphs (publisher MOESM647 text lines 286–309 and
319–332) adds the conditions the graph was missing:

- pool PCR with **Phusion HS Flex (NEB)**, 0.1 µl template, **63 °C annealing, 15 s
  extension, 20 cycles**; the **140-bp amplicon is size-selected on a 2% agarose E-Gel EX**;
- backbone **BsmBI (Fermentas)** digestion, **FastAP at 37 °C for 2 h**, gel-purified on a
  1% E-Gel EX;
- **20 µl Gibson reaction (NEB)** with **10 ng insert and 25 ng vector** for lentiCRISPR v2
  (5 ng insert for lentiGuide-Puro, a vector this page does not draw);
- electroporation into **Lucigen electrocompetent cells**, enough parallel transformations
  plated on 245 × 245 mm **carbenicillin (50 µg/ml)** plates to give **30–200× coverage**,
  colonies scraped and combined, then **Endotoxin-Free Plasmid Maxiprep (Qiagen)**.

🟢 The only sequencing this paper performs on its own libraries is **representation QC**,
not a screen readout: the plasmid library is PCR-amplified for **16 cycles with Phusion
Flash High-Fidelity**, with primers that add Illumina adaptors, **in duplicate and
barcoded** to quantify PCR bias, and sequenced on a **MiSeq, v3 150-cycle kit, 10% PhiX**.
Reads are demultiplexed (FASTX-Toolkit), adaptor-trimmed (cutadapt) and aligned with
Bowtie allowing one mismatch. The graph now ends that branch at a sequencing terminal whose
plan carries only this limitation: 🔴 the supplement prints neither the adapter-adding
primer strands nor the barcodes, so no library construct or primer site is drawn. That
remains an unresolved source question — the GeCKO website and Supplementary Data are still
unchecked for it — not an accepted absence.

The Joung one-PCR library shown at the end of the page is a separately cited compatible
readout for this vector, not part of Sanjana et al.'s experiment; the page now says so.

The original preprint continues with “Supplementary Data” (line 650 onward), including
labelled lentiCRISPR v2 and lentiGuide-Puro sequence material. The local Addgene #52961
full map also remains available in the pooled-vector reference directory; its provenance
and restoration route are recorded in that directory's MANIFEST.md. The historic
`genome-engineering.org/gecko/` link could not be opened in the current check; that
failed route does not establish that the original library data are unpublished.
