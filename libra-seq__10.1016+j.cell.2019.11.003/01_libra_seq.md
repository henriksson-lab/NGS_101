# LIBRA-seq — receptor sequences linked to antigen specificity

Setliff et al. (2019), *High-Throughput Mapping of B Cell Receptor Sequences to Antigen Specificity*, [10.1016/j.cell.2019.11.003](https://doi.org/10.1016/j.cell.2019.11.003).

## Bounded sources

The final-version author manuscript of the published article was read as
`j.cell.2019.11.003_PMC7158953.xml.txt`: antigen-barcode design, conjugation,
FACS, single-cell processing and read processing (198–220), plus Figure 1
caption (618). This is not an active preprint. The fetched source bundle was
scraped for oligos and relevant hits checked against these Methods. The six
linked PMC supplement routes were gated; none had a PMC Cloud counterpart.
Their contents are deferred, not confirmed missing. No challenge was bypassed.

## Source chemistry and distinct branches

🟢 Biotinylated antigen proteins receive 5′-amino DNA tags by Solulink
4FB/HyNic conjugation. Excess oligos are removed by FPLC. Tagged antigens are
complexed with streptavidin-PE, used to stain B cells, and antigen-positive cells
are sorted before Chromium B-cell Single Cell V(D)J processing.

🟢 The antigen tag contains a truncated small-RNA sequencing handle, antigen
identity and a sequence that anneals to the bead-delivered template-switch
oligo. It is **not** the poly(dA) reporter used by other protein assays.
Bead-delivered oligos index BCR transcripts and antigen tags in the same cells.
An additive antigen primer is included during cDNA amplification. SPRI separates
antigen reporter products from cellular cDNA; both branches are retained.
BCR preparation proceeds through the cited 10x V(D)J solution; antigen products
receive another 10–12 PCR cycles and 1.6× purification. Both libraries are
sequenced on NovaSeq6000. Detailed kit revision and complete molecular products
remain unassigned, rather than borrowed from a modern kit.

## Exact written antigen oligos

🟢 The common source pattern is
`CCTTGGCACCCGAGAATTCCANNNNNNNNNNNNNCCCATATAAGA*A*A`, with a 5′ amino
modification. The three starred bonds are preserved as phosphorothioates in the
source-written record. The additive primer is `CCTTGGCACCCGAGAATT*C*C`.
The prefix is imported from the reverse complement of canonical Illumina
small-RNA RA3 in the local model, after source verification.

The cell-line/NIAID45 branch uses 13-base antigen identifiers: BG505
`CATGATTGGCTCA`, CZA97 `TGTCCGGCAATAA`, influenza H1 New Caledonia
`GATCGTAATACCA`. The N90 branch instead lists 15-base identifiers. The current
graph represents the 13-base source branch; it does not silently apply its
window width to N90. No 10x cell-barcode length is supplied from memory.

🟢 The authors process reverse-complemented antigen R2 to locate the cell
barcode, then the immediately following ten-base UMI and antigen identifier.
That computational orientation is not automatically an exact physical read map.
Exact acquired cycles and terminal priming sites remain unassigned here.

## Scope and lineage

The graph separates antigen and BCR outputs; receptor heavy/light identities are
linked through cellular barcodes, not through an invented physical VH:VL junction.
Complete 10x V(D)J amplification/library chemistry and exact antigen library
primer geometry remain source-limited handoffs. The public page is a source-qualified partial published instance; it does not reconstruct complete historical kit molecules.

The Methods explicitly follows **CITE-seq** library preparation with more antigen
PCR cycles; this supports a component-scoped CITE-seq → LIBRA edge, not transfer
of CITE antibody or poly(A) capture chemistry. The cited Chromium B-cell V(D)J
solution is another explicit platform parent, with exact historical revision
unresolved. Deferred: supplements, cited CITE protocol version, 10x guide,
full bead/TSO/primer mapping, antigen panel variants and independent review.
No absence claim is accepted.

## Bounded public-page preparation

The same final-version author manuscript Methods 198–220 supports the public
13-base antigen branch. Whole version-of-record equivalence remains unestablished;
no supplements or historical kit guides were added to this pass. 4FB/S-HyNic
are molecular reactants, while FPLC is an apparatus-qualified operation with
unassembled column/mobile-phase inventory. Actual BCA/NanoDrop and conjugate
SDS-PAGE checks consume separate analytical allocations; no aliquot amounts
are invented. PBS/1% BSA and staining-marker contributions remain molecular
inputs. Streptavidin-PE is supplied incrementally from its actual stock at the
source 5:1 antigen:streptavidin molar ratio.

RT and the two PCR stages expose copresent enzyme/nucleotide/buffer functional
classes from their respective actual aggregate supplies. This bookkeeping
does not assert commercial premix packaging, new duplicate doses, exact
formulations or primer bodies. Final library loading allocations and unresolved
stock remainders stay separate for BCR and antigen outputs; no cross-branch
physical pool is inferred. The full workflow remains visible on the page.
