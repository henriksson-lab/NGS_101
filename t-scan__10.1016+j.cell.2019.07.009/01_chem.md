# T-Scan: bounded source-qualified chemistry

Kula et al. (2019), **T-Scan: A Genome-wide Method for the Systematic Discovery of T Cell Epitopes**, Cell. [Defining DOI](https://doi.org/10.1016/j.cell.2019.07.009).

The positive wet-preparation boundary is functional T-cell recognition of antigen-library epitope-discovery cells (EDCs), fluorogenic granzyme-B reporter selection, and antigen-cassette sequencing from recovered cells. This is a cell-selection front end, rather than an ordinary target panel. Premade engineered EDC and effector-T-cell stocks are inputs; their complete construction is unassembled.

The operational original is the version-qualified PMC6939866.1 author manuscript, overview85, NGS Methods289–295 and reference673–677. One original inspected; no superseded preprint used. Mandatory scraper29 candidates/one file; every selected whole primer body and index pattern was checked against the actual Methods. Final published main comparison and the unreviewed supplements remain ordinary scope, not confirmed missing disclosure.

🟢 Sorted EDCs are mixed with 2e5 wild-type HEK293T carrier cells **before** GeneJET gDNA extraction. Carrier DNA is not a UMI. A separate unsorted input-library preparation samples at least100-fold library representation; it is a control, not a physical mixture with the selected DNA.

🟢 The selected custom-primer screen route uses three Q5 hot-start PCRs. PCR1 uses1ug gDNA per100ul reaction; reactions for one screen replicate are pooled and1ul enters50ul PCR2. Then1ul PCR2 enters50ul PCR3. The graph preserves both untransferred stocks and does not insert an undocumented intermediate cleanup. Actual reaction counts and cycle choices remain unassigned. The Methods print the vendor as “Thermo”; this wording is retained rather than silently corrected from present-day product knowledge.

🟢 Source literal primers, written5′→3′ without additional inferred modifications:

| Primer | Verbatim source bases |
|---|---|
| TCell-PCR1-F | `TCCAGTCAGGTGTGATGCTCGGGGATCCAGGAATTCAGTTTGTACAAAAAAGCAGGCTCA` |
| TCell-PCR1-R | `CGAGCTTATCGTCGTCATCCCCACTTTGTACAAGAAAGCTGGGTCA` |
| TCell-PCR2-F | `ACACTCTTTCCCTACACGACTCCAGTCAGGTGTGATGCTC` |
| TCell-PCR2-R | `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCCGAGCTTATCGTCGTCATCC` |
| TCell-PCR3-F | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACTCCAGT` |
| TCell-PCR3-R | `CAAGCAGAAGACGGCATACGAGATxxxxxxxGTGACTGGAGTTCAGACGTGT` |
| T7-Pep2.2_Illumina_SP | `GGTGTGATGCTCGGGGATCCAGGAATTC` |

Here `xxxxxxx` is the source's seven-base sample-index window. Its actual sample-specific value is not assigned. It is a fixed sample label, not a molecule UMI; the model records that structured meaning. Canonical Illumina substrings are imported in the constructor. Full genomic endpoints and primer-site geometry remain unassembled; literal primers alone do not validate all PCR products.

🟢 Indexed sample PCR3 products are physically pooled before gel purification with QIAquick Gel Extraction columns. MiSeq or NextSeq sequencing uses the custom T7-Pep2.2_Illumina_SP primer. Exact run cycles/index-primer placement remain open. The separate human-genome-library eight-stagger-primer/standard-sequencing route is not merged into this custom-primer case; stagger diversity is not relabelled as a UMI.

The main identifies modified Xu2015 sequencing preparation as **Viral immunology. Comprehensive serological profiling of human populations using a synthetic human virome** (10.1126/science.aaa0698, reference673–677). This supports a scoped NGS component-adaptation link; the original itself remains unread. That directly cited preparation and historical GeneJET/QIAquick/Q5 documentation are deferred at depth one, not inherited molecular constructors. Source-selected molecular states, exact read geometry, engineered-stock construction and alternate screen branches remain incomplete. No source absence or full chemical validation is claimed.
