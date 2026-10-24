# SLIM — Site-directed, Ligase-Independent Mutagenesis

Chiu J, March PE, Lee R, Tillett D. *Site-directed, Ligase-Independent Mutagenesis (SLIM): a single-tube methodology approaching 100% efficiency in 4 h.* Nucleic Acids Research **32**, e174 (2004). [Defining published article](https://doi.org/10.1093/nar/gnh172), [PMC535700](https://pmc.ncbi.nlm.nih.gov/articles/PMC535700/).

## Source scope

Read the published article's Abstract, Methods **SLIM primer design**, **SLIM PCR amplification**, **SLIM hybridization**, Table1, Results **Design of SLIM primers**, and original Figure1A–C with its caption. The source bundle and URL outcomes are external; see `ref/MANIFEST.md`. Ran `tools/scrape_primers.py` on the bundle: twenty candidate sequences in the article. Inspected their Table1 context; this schematic uses only the insertion tail, not a source target gene or its specific binding-core sequences. The figure was inspected as its original article-linked JPG, including the separate productive 5′- and 3′-overhang classes and the nicked circular products.

This is a bounded, two-original pass (article and its directly linked figure), not a complete absence audit. Original insertion, deletion and substitution examples concern a T7 DNA-polymerase plasmid; no function-specific mutation or target-gene design is implemented here. The depicted insertion architecture uses the source-defined tag tail and a symbolic retained plasmid body.

## Oligo components

🟢 Table1 insertion primers contain the following **18-base 5′ tails**. These are components of longer primers, not complete oligos:

```text
forward insertion tail: CACCATCACCATCACCAT
reverse insertion tail: ATGGTGATGGTGATGGTG
```

The reverse tail is computed as the reverse complement of `slim.INSERT`. Both are independently checked against the published codon194 tailed-primer records. Table1's footnote identifies the underlined 18-base adapters; the source's full oligos include target-specific cores not used in this topology scope. No phosphate or other primer modification is assigned by this model.

🟢 FT and FS share a 3′ template-binding core; RT and RS share the opposite 3′ binding core. FT and RT add complementary tails. 🟡 In the model, the target-specific cores and retained plasmid body are inferred role tokens. No complete plasmid sequence or exact primer coordinates are claimed.

## Construct-changing operations

1. 🟢 Whole-plasmid inverse PCR combines FT, FS, RT and RS. The source uses Taq plus Platinum Pfx DNA polymerases. Four product classes are shown separately: FT+RS, FS+RT, FT+RT and FS+RS. They are separate molecules, not covalently concatenated.
2. 🟢 DpnI digestion follows amplification and acts on parental-template carryover. Digest fragments remain in the tube; no physical cleanup is invented.
3. 🟢 Denaturation and reannealing produce sixteen possible duplex classes, two of which are productive hybrids formed from opposite strands of Products1 and2. The drawing shows one class with exposed 5′ tails and one with exposed 3′ tails. Symbolic body alignment is not counted as established base pairing.
4. 🟢 Complementary exposed tails circularize noncovalently. Original Figure1A labels the products **nicked circular products**; the model does not introduce an in-vitro ligase or depict this stage as a covalently sealed plasmid. Because the exact nick positions are unassigned, this intermediate uses an explicit unsealed-circular-duplex limitation state, retaining the declared components without drawing a sealed reference strand or invented nick coordinates. The host-repaired final circle is a separate state.
5. 🟡 Host repair and recovery yield an expected insertion-construct candidate. Exact host-repair chemistry and a particular clone's genotype are not inferred from colony recovery. The fourteen other potential hybrids retain an unresolved disposition in the graph; they are not claimed to have been removed during hybridization.

The circle contains one retained insertion tail, not two covalently retained copies. No quantitative fraction, assembly efficiency, clone identity or target-protein property follows from this topology model.

## Final construct and readout

Expected repaired circular DNA comprises the retained plasmid body and one insertion. This is a **constructed-DNA endpoint**, with no sample-sheet index, UMI, NGS library or assigned sequencing run. The article's colony PCR and clone sequencing are downstream QC; they are not redefined as NGS preparation.

## Remaining scope

Deletion and substitution branches in original Figure1B–C need their own target-map declarations. Full source target maps, exact binding coordinates, all sixteen duplex structures and downstream QC are unexpanded. The directly cited predecessor Tillett and Neilan1999 and Molecular Cloning manual are deferred; no claim of missing disclosure or complete source review is made. The source checklist remains in progress and asserts no confirmed absence.
