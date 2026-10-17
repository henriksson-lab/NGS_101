# CAGE

## 1. What it is

CAGE enriches capped RNA 5′ ends by cap trapping and converts them into short sequencing tags for promoter and transcription-start-site mapping.

## 2. Authoritative source

🟢 Takahashi et al., *5′ end-centered expression profiling using cap-analysis gene expression (CAGE) and next-generation sequencing*, Nature Protocols (2012), DOI `10.1038/nprot.2012.005`.

## 3. Construct-changing operations

🟢 First-strand cDNA remains paired to RNA while the cap is oxidized and biotinylated. RNase digestion and streptavidin capture select cap-complete hybrids before linker installation and tag construction.

## 4. Exact oligos and read geometry

🟢 Table 1 gives RT-N15-EcoP as `AAGGTCTATCAGCAGNNNNNNNNNNNNNNN`, PCR forward as `AATGATACGGCGACCACCGACAGGTTCAGAGTTC`, PCR reverse as `CAAGCAGAAGACGGCATACGA`, and the custom sequencing primer as `CGGCGACCACCGACAGGTTCAGAGTTCTACAG`.

🟢 The 5′ linker contains a three-base sample barcode, `CAGCAG`, and a terminal N6/GN5 ligating tract. The 3′ linker begins with a two-base protrusion. EcoP15I releases a 27-nt cap-derived tag.

🟡 The final duplex is assembled from the overlapping PCR-primer and linker sequences. The custom run primer is located against that same construct rather than assumed to be standard TruSeq.

## Ligation geometry

🟢 For the AGA representative, the N6 upper strand is `CCACCGACAGGTTCAGAGTTCTACAGAGACAGCAGNNNNNN` with a 3′ phosphate; the GN5 alternative replaces the first random position with G. The matching lower is `/5Phos/CTGCTGTCTCTGTAGAACTCTGAACCTGTCGGTGG/3AmMO/`. The fixed 35-base portions form an exact antiparallel duplex. The upper strand's terminal N6/GN5 tail anneals transiently at a captured cDNA end, aligning the lower strand's 5′ phosphate beside the cDNA 3′-OH. This is a source-defined splinted ligation, not an unknown linker end.

🟢 The 3′ upper linker is `NNTCGTATGCCGTCTTCTGCTTG` and the lower is `CAAGCAGAAGACGGCATACGA`. The initial `NN` is therefore a **5′ two-base protrusion** on the upper strand. The earlier schematic placed it at the 3′ end; the oligo table fixes both its position and polarity.

🟢 The defining runs were single-read (36 cycles on GAIIx or 50 cycles on HiSeq 2000). The custom sequencing primer first reads the inline sample barcode and EcoP15I site, then the cap-derived transcript tag.

## Bounded source pass, 2026-10-17

🟢 Read: PMC4094379 author-manuscript JATS text (extracted `.xml.txt`)
Introduction/Experimental design lines 59–89 (cap trapping, SSLLM,
EcoP15I, PCR, sequencing), Materials lines 149–209, Procedure steps 1–61
(lines 400–912), Table 1 lines 1328–1507 and Box 1 lines 1567–1593.
The Springer ESM `MOESM283` is a mapping shell script; `MOESM282` is a
figure image. The PMC supplementary PDF remains a `(manual)` gated download.

🟢 Reaction order and reagents now in the model: PrimeScript RT with
sorbitol/trehalose (steps 1–4); NaIO4 oxidation then long-arm biotin
hydrazide (cap and RNA 3′ diols); RNase ONE digestion of unprotected
single-stranded RNA (step 16); capture on tRNA-coated MPG Streptavidin with
four wash buffers; release of single-stranded cDNA with 50 mM NaOH and
Tris neutralization (steps 23–25); 5′ linker ligation with DNA Ligation
Kit Mighty Mix (16 h); two AMPure purifications to remove 5′ linker, with
barcoded samples pooled (steps 33–34); second strand with the biotinylated
2nd SOL primer and LA Taq (94/42/68 °C, step 35–36); Antarctic phosphatase
then AMPure (37–39); EcoP15I with sinefungin, MgCl2 and 65 °C inactivation
(40–42); 3′ linker ligation with T4 DNA ligase at 16 °C for 16 h (43–44);
a second MPG Streptavidin capture that leaves 3′-linker dimers and
non-biotinylated material in solution (45–50); Phusion test-cycle and bulk
PCR from the beads (51–54); Exonuclease I and QIAquick (55–57).

🟡 Model correction: the earlier final construct placed two extra linker
bases after a full 27-nt tag and a 7-base P5 prefix, giving 98 bp while
labelled 96 bp. The paper states a 96-bp product (steps 53, 59) and a
27-nt tag. The model now uses the 13-base forward-primer tail
(`{{= cage.FWD_TAIL }}`) before the 26-base 5′ linker handle, and places the
3′ upper linker's 5′ NN at tag positions 26–27 opposite the cDNA-derived
two-base overhang; the 25/27 placement is derived from these lengths, not
stated in the paper. PCR primer binding is now computed exactly on the
pre-PCR dual-linker tag.

🟢 Published discrepancy, shown on the page: the Introduction says GN5 is
used "at a ratio of 5:1 with the random N6", while Box 1 mixes N6:GN5 at
1:4 (3 µg + 3 µg vs 12 µg + 12 µg per strand). The model draws both
species and chooses no ratio. Table 1 also labels the ACT lower strand
`5’SOL-lower-AGT` (a duplicate label); the AGA representative drawn here is
unaffected. The Introduction calls AdoMet "an analogue of sinefungin" while
the reagent list adds sinefungin; this wording does not change the model.

Open: the PMC supplementary PDF (gated), the cited predecessor SSLLM and
cap-trapper methods (refs 22, 23, 28) and the Illumina cluster/sequencing
kit documents were not inspected in this capped pass.
