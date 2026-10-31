# hmC-Seal

Defining DOI: https://doi.org/10.1038/nbt.1732

## Sources read and scope

🟢 `nbt.1732_PMC3107705.xml.txt`: Associated author manuscript NIHMS295550 identity, abstract73, Results79–95 and Fig1/2/4. Main delegates Methods to online publication. Selective azido-glucose/cyclooctyne-biotin chemistry and actual mouse genomic deep sequencing; separate HeLaHRMS and streptavidin-stop validations.

🟢 `nbt.1732_41587_2011_BFnbt1732_MOESM6_ESM.pdf.txt`: Published supplement with DOI header; synthesis26–137 and selected Fig5/7 captions180–217. Long gene/data lists not chemically reviewed. Published supplement title differs slightly from final article title, but DOI/authors/figure identity match. Synthesis/controls inspected, not a genomic library recipe.

Capped pass: at most10minutes/3originals, published supplements preferred over associated-manuscript aliases where available. No linked-reference recursion, two-reader absence audit or finished-protocol promotion.

## Source chemistry and model

🟢 The original uses T4 β-glucosyltransferase to transfer 6-azido-glucose from UDP-6-N3-Glu onto pre-existing5hmC, then copper-free click chemistry with cyclooctyne-biotin compound1. The selected chemical state is biotin-N3-5gmC: the genomic cytosine remains a cytosine nucleotide with a covalent substituent. No C→T conversion or guaranteed all-site tagging is implied. The local constructor requires caller-declared5hmC-bearing cytosines and preserves nucleotide identity; it does not label every genomic C.

🟢 Main describes sonicated100–500bp genomic DNA, biotin-affinity enrichment, multiple actual mouse-cerebellum deep-sequencing libraries and input/no-biotin controls (regular UDP-Glu). Source input and no-biotin controls are distinct experiments; this selected graph follows the chemically labeled branch and does not duplicate input DNA into it. Actual capture/recovery reagent, end preparation, adapters, amplification and raw run schedule remain unassembled because the associated main delegates online Methods. Their source-known library/acquisition boundary remains in the graph without fabricated recipe details or molecular disposal. No DTT/disulfide chemistry is borrowed from the later SMRT method.

🟢 The accessible published supplement identifies UDP-azido-glucose synthesis and selected validation controls. Its source title says “intragenic enrichment ... in mouse cerebellum” while the final article says “genome-wide distribution”; DOI10.1038/nbt.1732 and matching authors/figures establish the supplement association. This ordinary title difference is not asserted as a chemical conflict. Much of the large text extraction is data; no claim is made to have inspected every gene row or image-based oligo table. Mandatory scraper0hits is not absence evidence. The actual PMC supplement aliases remain gated, and their equivalence to the published supplement was not assumed from filename.

### Boundaries and readout

Region enrichment from fragment sequencing is not itself a single-base assay. Main's streptavidin-dependent polymerase-stop experiment is a separate proposed locus readout; it is not inserted into the actual genomic NGS library. HeLa monomeric-avidin purification followed by digestion/HPLC/HRMS is also separate and does not establish the exact mouse NGS capture/recovery recipe. No selective cleavage, ideal purification, universal label yield or synthetic oligo sequence is imposed on the genomic branch.

Actual online Methods and historical library/vendor originals remain relevant direct targets, not confirmed missing disclosures. Final publication chemistry applicability remains explicit. Nano-hmC-Seal, cell-free-DNA applications and the later disulfide-SMRT method are not imported as substitutes.

## Actual online Methods target and access history

The relevant final online Methods target is https://www.nature.com/articles/nbt.1732 (associatedmain106 delegates online Methods). The prior discovery browser request returned a tool500 after redirect to an inaccessible `idp.nature.com/authorize` endpoint. This is an access-route failure, not evidence of absent Methods. Current standard retrieval obtained the associated main, whose excerpt omits those online Methods. The available publisher supplement was already read; its gated PMC aliases are recorded below as retrieval history, **not requests to reacquire an already inspected original merely because byte equivalence is unverified**.

## Actual access barriers and deferred originals

- https://pmc.ncbi.nlm.nih.gov/articles/instance/3107705/bin/NIHMS295550-supplement-data.pdf — PMC supplementary file; PMC download gate (not in the PMC Cloud bucket: licence TDM); the publisher's ESM files were fetched, so probably redundant
- https://pmc.ncbi.nlm.nih.gov/articles/instance/3107705/bin/NIHMS295550-supplement-supplement_1.pdf — PMC supplementary file; PMC 'Preparing to download' gate; the publisher's ESM files were fetched, so probably redundant

Manifest wording such as “probably redundant” is a fetcher hint, not an identity/equivalence conclusion. Exact inspected scopes above determine authority; ordinary deferred historical geometry is not an access failure.

## Validation

Protocol-owned `tools/chem.py` has a zero-argument graph and focused state/fate checks. It is an incomplete graph-only implementation, with no public page or completion marker. Source-known late operations/acquisition are retained; unassembled adapters, modifications and fates are not inferred from other families.

## Bounded unlinked-page adjudication (2026-10-30)

🟢 Associated manuscript76/80 identifies β-glucosyltransferase, UDP-6-N3-Glu and cyclooctyne-biotin compound1 as actual chemical contributors. The genomic graph now separates the glycosyl donor as a molecular reagent from the β-GT catalyst, and exposes compound1 as the molecular click reagent rather than a condition. Genomic doses, actual packaging, transferred fractions and full chemical products remain unassembled. The separate 11mer validation's enzyme fractions, reaction times and yield are not assigned to genomic DNA. No actual new reaction or bench condition is supplied by the inventory correction.

The actual online-Methods/library SourceItem remains uncovered. Source multiple biological samples/libraries, input and regular-UDP-Glu controls stay distinct obligations; no one-product array, modern library kit, HeLa recovery procedure or later disulfide-cleavage route is substituted. Typing an opaque boundary as a library would not establish its complete producer. The source pass remains closed on the associated main74–106 and selected published supplement synthesis/validation contexts, with source-version limits retained. No new original or access request, source-absence conclusion, public page or completion marker is introduced.
