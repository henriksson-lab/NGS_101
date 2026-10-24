# Exponential Megapriming PCR (EMP) cloning

Ulrich A, Andersen KR, Schwartz TU. *Exponential Megapriming PCR (EMP) Cloning—Seamless DNA Insertion into Any Target Plasmid without Sequence Constraints.* PLOS ONE **7**, e53360, published31December**2012**. [Defining article](https://doi.org/10.1371/journal.pone.0053360), [PMC3534072](https://pmc.ncbi.nlm.nih.gov/articles/PMC3534072/). The PMID entry's later indexing does not change the defining publication year.

## Sources inspected

Read the original published PDF, Figure1A on p2, p3 **Exponential Megapriming PCR**, **Synthesis of the EMP megaprimer**, **Insertion of the EMP megaprimer**, **In-vitro ligation and transformation of EMP product**, and Discussion pp8–9 inverse-PCR comparison. Inspected the rendered original published TableS1, `pone.0053360.s004.tif`, from the validated PMC supplementary archive. Its table image contains the actual EMP protocol even though it has no extracted-text twin. The archive's article-linked inventory establishes its identity; the table's printed title is **Table S1. Exponential megapriming PCR protocol**.

Ran the primer scraper on the external source bundle: zero candidate sequences were reported. Read the published PDF text and rendered original directly. This result establishes no absence of primer sequences. No exact oligo sequence is transcribed or target sequence designed in this schematic. This bounded pass covers two original documents; it is not a full independent absence review.

## Primer architecture and molecular order

🟢 **F1** matches the insert's 5′ end and is reused during second-stage amplification. **R1** has a 5′ recipient-homology tail, reverse complementary to the recipient region immediately downstream of insertion, followed by the reverse complement of the insert's 3′ end. **R2** binds the recipient region upstream of insertion. Methods describe typical F1 and R2 lengths of20–25nt and R1 of40–50nt; these ranges are not converted into a chosen primer sequence or a single asserted length.

🟡 Neutral payload, landing site B and retained backbone are inferred role tokens. No source test plasmid, functional insert or experimental clone is reconstructed. The roles retain the source-defined order, but do not establish base-pair identity, unique priming or ligation-ready end geometry.

1. 🟢 Phusion PCR with F1 and R1 amplifies the insert and adds the recipient-homology segment B. Original Figure1A shows a **duplex product**, with B present on both strands. The phrase “3′ overhang” in Methods describes the megaprimer's active strand in the insertion mechanism; it does not justify drawing an already purified duplex as permanently single-stranded at that end.
2. 🟢 Purify the first product. Upon denaturation, the active megaprimer strand's 3′ homology binds B downstream of insertion on the circular recipient. R2 binds upstream. Initial extension generates insert–recipient molecules, after which F1 and R2 amplify the linear product exponentially. This is the original two-step branch, with F1 included in stage2; the optional F1-omission branch is not silently substituted.
3. 🟢 Purify the second-PCR product, then T4 PNK supplies 5′ phosphates. T4 DNA ligase circularizes the product. The circle's closure lies between the insert start and upstream recipient backbone; exact end bases and a complete physical junction are not assigned in this role-only instance.
4. 🟢 DpnI acts **after phosphorylation and ligation**, removing remaining parental-plasmid background. Digest fragments are shown as still present until colony recovery; purification is not assumed to prove complete parental-template clearance.
5. 🟡 Recover and screen colony candidates. Colony recovery does not establish genotype or insert function. The final constructed-DNA state is an expected architecture, with no NGS run or sequencing library implied.

## Published discrepancy

**Methods, p3, Insertion of the EMP megaprimer:** specifies **25 cycles** for the second EMP PCR. **TableS1, second-EMP cycle scheme:** specifies **35 cycles**. The lower-left TableS1 reagent block is headed **“2nd RF PCR”**, although the table title, adjoining cycle block, F1/R2 ingredients and megaprimer treatment identify the EMP context. These are values/headings printed in the published originals, not errors introduced by extraction.

The note follows the main Methods'25-cycle second-PCR statement while preserving the TableS1 value; the model's molecular transformations do not depend on cycle count. The intended correction and whether the experiments used one or both counts are unconfirmed. The public schematic displays this discrepancy rather than silently reconciling it.

## Boundaries

No complete recipient or insert map, exact F1/R1/R2 oligo, uniquely located site, ligation end sequence or particular clone genotype is claimed. All generic body/site roles are constructed with inferred styling. One-step EMP, multi-insert EMP, and the RF comparison are separate source branches and remain outside this original two-step schematic. TableS2–S4, supplementary result figures, all study-specific constructs and cited predecessor inverse-PCR methods are deferred, not judged irrelevant or absent. Source review remains in progress with no confirmed missing-information claims.
