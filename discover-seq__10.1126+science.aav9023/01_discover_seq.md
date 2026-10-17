# DISCOVER-Seq — research notes

Wienert et al., *Science* (2019),
[doi:10.1126/science.aav9023](https://doi.org/10.1126/science.aav9023).

## Assay target

The sequenced material is human chromatin fragments enriched by MRE11 occupancy at nuclease-induced double-strand breaks. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a stock ChIP-seq library kit). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining Science author manuscript and associated description of MRE11-bound
  read endpoints.
- 🟢 Wienert et al., *Nature Protocols* (2020),
  [doi:10.1038/s41596-020-0309-5](https://doi.org/10.1038/s41596-020-0309-5),
  the authoritative stepwise DISCOVER-Seq procedure.
- 🟢 Illumina, *TruSeq ChIP Sample Preparation Guide*, part 15023092 Rev. B, the manual
  explicitly named by the detailed protocol.
- 🟢 van de Venn, Kontarakis & Wienert, *DISCOVER-Seq: MRE11 ChIP Seq* V.2 (30 March 2023),
  [dx.doi.org/10.17504/protocols.io.kqdg349ypl25/v2](https://dx.doi.org/10.17504/protocols.io.kqdg349ypl25/v2),
  the authors' bench version of the same procedure (protocol integer ID 69486).

## Construct logic

1. 🟢 CRISPR nuclease is delivered to living cells or tissue and MRE11 is allowed to bind
   newly made double-strand breaks. The appropriate capture time depends on delivery;
   the detailed protocol suggests roughly 8–16 h for RNP editing.
2. 🟢 Cells are crosslinked in 1% formaldehyde for 15 min at room temperature, quenched with
   125 mM glycine, lysed, and the nuclei are sheared on a Covaris S2 (Nature Protocols
   lines 657–667, 707–729). The same article states the target chromatin size twice:
   “approximately 150–500 bp” in the sonication overview (line 113) and “median size of
   200–500 bp” in Anticipated results (line 994); the authors' protocols.io v2 expects a
   200–500 bp smear (line 292). The page shows this discrepancy and the model leaves the
   window unassigned, because the authors state sonication must be optimized per cell type.
3. 🟢 MRE11-bound chromatin is immunoprecipitated with a human/mouse cross-reactive
   anti-MRE11 antibody (Novus NB100-142, 10 µL per ChIP, or 10 µg ab208020) prebound to
   Protein A Dynabeads (Nature Protocols lines 679–687; protocols.io lines 208–236).
   Crosslinks are reversed at 65 °C for at least 6 h; RNase A (37 °C, 30 min) and
   proteinase K (55 °C, 1 h) are followed by Qiagen MinElute purification and elution in
   41 µL water (lines 765–791). ChIP-qPCR at the intended target is a checkpoint before
   library preparation.
4. 🟢 The detailed procedure delegates the remaining chemistry to the Illumina TruSeq
   ChIP library kit (Nature Protocols lines 799–845). The kit manual specifies end repair,
   3′ adenylation, indexed-adapter ligation, size selection and PCR enrichment.
5. 🟢 The historical kit is single-indexed. Its guide says to record six index bases;
   AR001 is `ATCACG`, with the parenthesized seventh base explicitly excluded from the
   sample-sheet index (TruSeq guide lines 2109–2186).
6. 🟢 Individual libraries are quantified by Qubit, pooled in equal amounts and cleaned of
   adapter dimers with 0.8× AMPure XP beads, repeating the cleanup (or using PippinPrep or
   a gel) if dimers persist; the pooled library is expected to average 300–500 bp
   (Nature Protocols lines 809–841, 996). Libraries are sequenced as 100-bp paired-end
   reads on a HiSeq 4000, NovaSeq or similar at about 30–40 million read pairs
   per ChIP sample (Nature Protocols lines 843–849). The BLENDER analysis uses the pileup
   and orientation of fragment endpoints around a cut rather than an assay-specific UMI.

## Why the end is informative

One end of an enriched fragment can be the cellular nuclease break occupied by MRE11;
the other is normally a random sonication boundary. A real cut therefore produces many
MRE11-enriched fragments whose reads start or end at a common genomic coordinate, while
their distal endpoints vary. The paired-end construct on the public page preserves that
asymmetry without pretending the genomic sequence is fixed.

## Commercial-sequence boundary

The six-base adapter indexes and the final TruSeq sequencing architecture are public.
The discontinued ChIP kit guide does not provide every supplied reagent as a complete
ordering oligo, so the page does not invent such sequences or call the library dual-indexed.
