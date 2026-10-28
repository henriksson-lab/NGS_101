# DAb-seq — jointly barcoded genomic amplicons and antibody tags

Demaree et al. (2021), *Joint profiling of DNA and proteins in single cells to dissect genotype-phenotype associations in leukemia*, [10.1038/s41467-021-21810-3](https://doi.org/10.1038/s41467-021-21810-3).

## Source scope and revision identity

The published article `s41467-021-21810-3_PMC7952600.xml.txt` was inspected at
Results overview (72–80) and Methods conjugation, staining, microfluidic capture,
library preparation and sequencing/bioinformatics (138–176). The published
supplement `s41467-021-21810-3_supp_41467_2021_21810_MOESM1_ESM.pdf` was read at
Tables 1, 2 and 6, including rendered PDF pages 14, 15 and 19. These tables are
molecular authorities; the separate reviewer-report PDF is not a bench protocol.
The full-bundle oligo scraper was run; published sequence/modification hits were
checked against these originals. Superseded preprint material was excluded.

The automatically retrieved detailed protocol identifies itself as **Version 2,
28 October 2020**, DOI **10.17504/protocols.io.bn4ymgxw**, on its printed cover.
That identity matches the article's cited DOI. Only its cover/identity was
inspected in this capped pass; its chemistry remains a deferred third original.
The Mission Bio protocol inherited by the article is **AML V1**. A current
Tapestri export is not assumed equivalent. The source-data ZIP holds derived
expression/genotype files, not a replacement for source oligo tables.

## Distinct library branches

🟢 Antibodies receive DBCO-PEG5-NHS and 5′-azide DNA tags by copper-free click
chemistry. After staining/washing, Tapestri droplets lyse/proteolyze cells,
then heat-inactivate the protease. UV-released barcoded primers amplify DNA
loci and capture antibody tags in the same compartment. **RUNX1_4's forward
primer** is reused to capture tags; Table 1 explicitly omits its reverse primer.
The model does not treat this as an ordinary unmodified genomic primer panel.

🟢 After emulsion recovery, 0.75× AMPure separates genomic amplicons (beads)
from shorter antibody-tag material (retained supernatant). The latter is
hybridized to a **5′ biotin/3′ C6** probe and captured on streptavidin beads.
Genomic products receive custom P5/Nextera P7 library PCR (10 cycles);
protein reporters receive custom P5/antibody-specific P7 PCR (20 cycles).
Each branch is cleaned independently with 0.7× AMPure.

## Exact source oligos and computed binding

🟢 Table 2 prints the antibody input as:

`/5AzideN/GGCTTGTTGTGATTCGACGATNNNNN[barcode]NNNNNAAGCTTGTTTCTGTGCACTGAG`

The two five-base random windows jointly form the source's ten-base UMI and
flank the eight-base antibody identity (CD3 example `AACGCTTC`). The model
preserves both UMI components with shared identity metadata. It does not invent
cell barcode values or a complete final library.

🟢 Table 1 gives the forward primer
`CTCAGTGCACAGAAACAAG*C*T*T`. Its complementarity to the tag's 3′ site is
computed, not eyeballed. The three phosphorothioate bonds are retained in the
source-written record; they are not converted into a generic terminal modification.

🟢 Methods prints capture probe `/5Biosg/GGCTTGTTGTGATTCGACGA/3C6/`.
The probe has the **same** sequence as the input tag's 5′ prefix; it anneals to
the PCR-generated complementary strand following denaturation, not directly
to the source-written input strand. The local Scene computes this relationship.
The blocked capture probe is distinct from library-PCR primers.

🟢 Table 6 gives P5 indexes 1–16 and antibody P7 indexes 1–4, including
phosphorothioates in the latter. The local module includes source index-1
representatives as examples of the documented primer pools, not assignments
of that pair to any particular patient. Full PCR primer-to-product coordinates
remain unresolved because the entire Tapestri barcoded primer is unassigned.

## Readout and limits

🟢 All runs are dual-indexed and use the custom Read 1 primer
`GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAG`. In the patient branch, DNA is sequenced
2×150 on NovaSeq6000; antibody tags use NextSeq550 Read1=38/Read2=39. Cell
barcodes are parsed from Read1, while antibody identity and UMI are extracted
from Read2. The earlier three-cell control uses MiSeq 2×150 and is kept separate
from the patient-run plan. Exact per-cycle barcode spans and complete terminal
primer geometry remain unassigned; no sequencing platform defaults fill them.

Deferred: detailed protocol V2 chemistry; historical AML V1 guide and complete
bead barcode/primer design; all source primer pools and patient index assignment;
final-library/primer-coordinate reconstruction; complete original/linked-source
inventory and independent review. This is a capped source pass, not accepted
absence evidence or a finished protocol. Catalogue status remains notes.

## Lineage

The published Results explicitly states that DAb-seq modifies Tapestri to support
**Abseq**, its own reference 23 (Shahi et al. 2017). This is a source-qualified
Abseq → DAb-seq edge for antibody-DNA measurement and joint droplet capture,
not inheritance of Abseq's thiol conjugation or clonal randomer barcode recipe:
DAb-seq supplies different click chemistry and a commercial bead-PCR workflow.
Mission Bio AML V1 is another explicit platform/component parent whose historical
original remains deferred. No RNA branch is implied by the multiomic label.
