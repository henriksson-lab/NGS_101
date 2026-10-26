# BioBrick / RFC10 — illustrative 3A assembly

[Shetty, Endy and Knight2008](https://doi.org/10.1186/1754-1611-2-5), published article, Results “Assembling BioBrick parts using a new BioBrick vector”, Figure3 and Methods “Assembly of BioBrick parts using the new BioBrick vectors”. The [noncoding BioBrick standard](https://openwetware.org/mediawiki/index.php?title=Biobrick_standard&oldid=118773) supplies the prefix/suffix strings. The separately archived [RFC10](https://hdl.handle.net/1721.1/45138), dated May3,2007, is an identified follow-up document; its original export has not been inspected in this capped pass.

## Scope and inputs

This is an **invented illustrative example**, not one of the paper’s experimentally assembled constructs. It selects the paper’s three-fragment part-assembly branch. Its short neutral bodies and recipient reference have no replicon, functional genes or selection markers; they do not constitute a usable vector. “RFC10” here identifies the noncoding prefix/suffix convention, not a claim of complete standard compliance or validated cloning performance.

🟢 The source branch digests prefix DNA with EcoRI/SpeI, suffix DNA with XbaI/PstI, and destination DNA with EcoRI/PstI. The Methods dephosphorylate the destination using Antarctic phosphatase, purify the preparations, ligate using T4 DNA ligase, and transform cells. The example preserves this dephosphorylation rather than assuming all strands become sealed in vitro.

The standard noncoding flanks are rendered directly from the local declaration:

```
{{= biobrick.PREFIX }}
{{= biobrick.SUFFIX }}
```

## Molecular relationships

🟡 Full duplex reference arcs precede each toy digest. Explicit local recognition-site and cut-coordinate declarations derive retained physical strands. EcoRI and SpeI/XbaI have 5-prime cohesive ends; PstI has a 3-prime cohesive end. The coordinate checker compares those declared intervals with the retained fragments; it does not predict enzyme specificity or cleavage distributions. Complete source parent maps are outside this reference-arc scope.

The selected circular order is prefix, suffix, recipient. Matching ends are paired by Scene; the reference retains the SpeI/XbaI fusion and its two flanking standard bases once, giving the eight-base `TACTAGAG` noncoding scar. External standard flanks remain available in the intended composite. Exact short bodies and the recipient arc are inferred model choices.

Recipient 5-prime hydroxyls are displayed explicitly. The ligation product is an unsealed circular class: at least its corresponding recipient-strand nicks cannot be closed by the stated ligase reaction. Exact nick positions, ligation fractions and repair are unassigned. Host repair/recovery is an unresolved boundary before a conditional intended circular reference. No clone, biological function or quantitative yield is established by that reference.

Unretained parent intervals, incomplete digests and alternative assemblies are unresolved product classes. They are not concatenated into invented excised molecules or declared physically removed. Source preparation cleanup is represented without treating it as a perfect separation. NGS adapters, sample indexes and sequencing read plans are inapplicable to this endpoint.

## Review boundary

The primer scraper was run on the defining article and extracted bundle. Its source-primer hits concern functional source constructs, not these neutral examples; none is adopted. The published article sections above and OpenWetWare noncoding standard text were inspected. The archived RFC10 original, figure graphics, functional vector maps and complete linked-source inventory remain deferred. No missing-information conclusion is accepted; the independent two-reader review remains unfinished.
