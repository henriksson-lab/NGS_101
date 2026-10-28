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
primer geometry remain source-limited handoffs. Status remains notes.

The Methods explicitly follows **CITE-seq** library preparation with more antigen
PCR cycles; this supports a component-scoped CITE-seq → LIBRA edge, not transfer
of CITE antibody or poly(A) capture chemistry. The cited Chromium B-cell V(D)J
solution is another explicit platform parent, with exact historical revision
unresolved. Deferred: supplements, cited CITE protocol version, 10x guide,
full bead/TSO/primer mapping, antigen panel variants and independent review.
No absence claim is accepted.
