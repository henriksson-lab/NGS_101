# DNase-seq (Boyle 2008)

Boyle et al. (2008), High-Resolution Mapping and Characterization of Open Chromatin across the Genome. DOI: https://doi.org/10.1016/j.cell.2007.12.014

This entry represents the explicitly described Solexa/Illumina branch in the 2008 paper. The paper also sequences DNase tags on 454 and includes a different DNase-chip array branch; those preparations are not merged into the Solexa graph.

## Sources and scope

🟢 Associated accepted author manuscript, PMC2669738, identity text lines1003–1006; Results initial library paragraph; Experimental Procedures / DNase Sample Preparation (line1303) and High-Throughput Sequencing and Data Analysis (1307–1309). The full article links a 9.2 MB supplemental PDF, whose generic fetch returned the PMC “Preparing to download” gate. That relevant original is uninspected, not absent.

## Verbatim oligos from DNase Sample Preparation

All sequences are written 5′→3′ as in the original; modifications belong to the particular oligo. Their original linker names are retained.

- 🟢 Linker I top: 5′ Bio-ACAGGTTCAGAGTTCTACAGTCCGAC.
- 🟢 Linker I bottom: 5′ P-GTCGGACTGTAGAACTCTGAAC.
- 🟢 Linker II top: 5′ P-TCGTATGCCGTCTTCTGCTTG.
- 🟢 Linker II bottom: 5′ CAAGCAGAAGACGGCATACGANN.
- 🟢 PCR primer: 5′ CAAGCAGAAGACGGCATACGA.
- 🟢 PCR primer: 5′ AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA.
- 🟢 Sequencing primer: 5′ CCACCGACAGGTTCAGAGTTCTACAGTCCGAC.

The terminal NN in linker II is a degenerate ligation end, not an asserted UMI, index or barcode. The sequencing primer differs from the PCR primer: the latter does not include the final C of linker I. Copying that remaining C is a candidate extension relation, not a fully validated final-library assembly.

## Preparation and readout

1. 🟢 Digest intact primary human CD4+ T-cell nuclei with DNase I. The accessible main text does not itself resolve the source-local nuclease titration, end polishing or exact substrate-end distribution; the supplement remains a relevant follow-up. No canonical repair step is inserted merely because ligation needs a compatible end.
2. 🟢 Ligate biotinylated linker I to DNase ends. Its 3′ end contains the stated MmeI recognition site. A model can compute linker oligo pairing, but pairing alone does not establish the actual genomic-end ligation geometry.
3. 🟢 Digest with MmeI, then enrich DNase ends on streptavidin beads. Cleaved unbound genomic DNA stays in the reaction after digestion and is physically separated only at capture/washing. MmeI offsets and exact tag-end geometry remain deferred pending the chemistry supplement/original enzyme conditions.
4. 🟢 Ligate linker II, then amplify with the printed linker-specific primers and PAGE-purify. The partial model does not invent bead-release conditions, cycler settings or an N-window identifier meaning.
5. 🟢 Sequence with the printed linker-I primer; the analysis explicitly reports 20 bp Illumina tags. The 454 ~108 bp data constitute another readout branch, not a second read of each Illumina molecule.
6. 🟢 Uniquely align and collapse reads sharing a 5′ start on the same strand. This is computational duplicate filtering after sequence acquisition; no UMI or molecular pre-sequencing discard is implied.

## Modelling and source boundaries

🟡 Exact oligo strings and computed linker-overlap states are owned locally; the complete PCR/readout topology is unassembled. Adapter input states retain their original biotin/phosphate modifications. A valid primer landing proves a supplied sequence relation, not that an unread supplement establishes all joining substrates, enzyme behavior or amplification conditions. Source caveats remain on the preview.

The author manuscript and printed main-text primers are source-positive. Complete VOR/supplement review, primary preparation details and the 454 branch remain deferred. No missing-information claim has passed the independent two-reader gate.

## Lineage

The introduction cites Crawford et al. (2006), DOI 10.1101/gr.4074106, as earlier DNase-seq using MPSS, and refCrawford2006a for DNase-chip. These establish historical relationships at assay scope; the Solexa MmeI/linker recipe must not be imported from the different MPSS or sonication/array branches without specific original evidence.


## Public-page hold adjudication, 2026-10-30

**Hold: required source scope remains unreviewed.** The actual normal graph still
reports `authoritative source item 'unreviewed_original' is not accounted for`.
Its requiredness and coverage are unchanged. Source-exact oligos and local
pairing checks do not cover the relevant original chemistry supplement.

The only validated local original in this pass is the associated accepted
manuscript, [PMC2669738](https://pmc.ncbi.nlm.nih.gov/articles/PMC2669738/),
external `j.cell.2007.12.014_PMC2669738.html(.txt)`. Its own header at text1003–1006
says “Author manuscript” and separately cites the final Cell2008 publication;
complete accepted/final-version equivalence is not asserted. Experimental
Procedures1303 prints the seven modified oligos and Solexa
linkerI→MmeI→streptavidin→linkerII→PCR/PAGE order; readout1309 reports20bp
Illumina tags. The454 and array branches remain separate. The degenerate
linkerII NN end is not assigned UMI or sample-index meaning.

The actual linked [NIHMS106456 supplemental PDF](https://pmc.ncbi.nlm.nih.gov/articles/instance/2669738/bin/NIHMS106456-supplement-supplemental.pdf)
has only a manual “Preparing-to-download” manifest row and no validated local
body. Nuclei preparation/titration, DNase-end preparation, MmeI end/offset
geometry and intervening workup/ligation conditions remain source-review
obligations. This is an unreviewed/unavailable original, **not** evidence those
details are absent from the publication. No new source was fetched or
substituted in this pass.

Mandatory scraper: `/tmp/unlinked-dnase7-scraper.log`; relevant main-text orders
and modifications were checked. Nine preflight source/guard checks passed;
chemistry7steps/4localpairingrelations/8unresolved/0contradictory. These results
validate only the stated local relations and intentional publication hold, not
a finished preparation. No code or graph gate was changed; no tests repeated
for this documentation-only adjudication.
