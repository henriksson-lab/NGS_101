# in4mer: four-guide array construction and sequencing

Esmaeili Anvar N et al. (2024).
[Efficient gene knockout and genetic interaction screening using the in4mer CRISPR/Cas12a multiplex knockout platform](https://doi.org/10.1038/s41467-024-47795-3).
Nature Communications 15:3577.

## Sources inspected

Bounded passes, 2026-10-21. The source-scoped four-guide graph and schematic are implemented.

- Defining published article, PMC11055879: own DOI verified; Methods “One-component
  CRISPR/enCas12a vector”, “7mer library production”, “In4mer prototype library
  production”, “Inzolia library production”, and both genomic-DNA sequencing
  subsections. Local JATS text lines 197–231 and 251–261.
- Supplementary Data 10: original workbook
  `s41467-024-47795-3_41467_2024_47795_MOESM13_ESM.xlsx`, sheet
  `Supplementary Data 10`, oligonucleotide names and sequences inspected.
  The article's Data 10 is **MOESM13**, not MOESM10.
- Ran `tools/scrape_primers.py` across the retrieved article and supplements with
  two displayed hits per file. Followed the article array declarations and the
  Data 10 primer rows directly; scraper hits are not the evidence authority.

Originals are under
`$CHEM_DATA/sources/in4mer__10.1038+s41467-024-47795-3/`.
This pass neither establishes source completeness nor confirms missing information.

## Distinct published branches

🟢 The article contains a seven-guide test, a four-guide prototype and the final
Inzolia library. The four-guide arrays preserve ordered guide positions separated
by distinct direct-repeat sequences. A guide combination is not an unordered set.

🟢 Prototype and Inzolia use different oligo-pool amplification handles and different
published assembly reagents. Prototype construction names BsmBI-v2 assembly into
pRDA-550. Inzolia names Esp3I and T7 ligase, with pRDA-550 and pRDA-052 vector
branches. A shared final guide-array core is not evidence of identical preparation.

🟢 Vector-cassette orientation is specified in the article's one-component-vector
section. Keep the transfer cassette separate from the unmodelled rest of the vector;
do not infer a full vector map from the oligo-pool sequence.

## Oligos inspected verbatim

The amplification primers in Supplementary Data 10, written 5′→3′:

| Branch | Primer | Sequence |
|---|---|---|
| Prototype / 7mer | Forward | {{= in4mer.PROTOTYPE_F }} |
| Prototype / 7mer | Reverse | {{= in4mer.PROTOTYPE_R }} |
| Inzolia | Forward | {{= in4mer.INZOLIA_F }} |
| Inzolia | Reverse | {{= in4mer.INZOLIA_R }} |

🟢 Prototype custom sequencing forward primer:
`{{= in4mer.CUSTOM_READ1 }}`.
The workbook also labels `ACCGACTCGGTGCCACTTTTTCAAGACCAG` as the reverse
sequencing primer for the **7mer only**. It must not be silently assigned to the
single-end four-guide run.

The last workbook block repeats the prototype/7mer genomic-library heading while
listing vector-specific pRDA-052 and pRDA-550 primers. Its forward primers include
the literal token `[s]`. Preserve that token in any verbatim transcription until its
meaning and experiment assignment are checked against an original definition. Do not
replace it with invented bases or silently drop it from a source-exact oligo.

## Recovery and sequencing

🟢 Both four-guide branches recover array-bearing amplicons from genomic DNA by
one-step indexed PCR, then pool and purify the products.

🟢 Prototype: NextSeq 500, custom primers, read format 151–8–8, single-end.
The 7mer test separately uses 151–8–8–151, paired-end.

🟢 Inzolia: the genomic-PCR subsection names standard eight-base Illumina index
sets D501–D508 and D701–D712 and HiSeq2500 Rapid Run, 200-cycle kit. A kit capacity
is not an exact assignment of cycles to individual reads; do not copy the prototype
151–8–8 recipe into this branch without further evidence.

## Implementation and review queue

The current graph has three alternative construction/readout routes: prototype into
pRDA-550, and Inzolia into pRDA-550 or pRDA-052. Inzolia preparation is repeated in
the view to show the two recipients; it does not describe two sequential assemblies.
The routes preserve their own amplification handles, polymerases, cleanup reagents,
assembly enzymes and ligases. Source carbenicillin selection is a shared real reagent,
without implying equal growth conditions, recovery or genotype fidelity.
The printed pool templates have computed lengths of
{{= in4mer.physical_length(in4mer.synthesis('prototype')) }} nt for the prototype and
{{= in4mer.physical_length(in4mer.synthesis('inzolia')) }} nt for Inzolia, including
their distinct amplification handles and the declared guide windows.

The guide windows have structured positional guide-identity features. The sample
indexes have separate i5/i7 features, added only by the genomic-PCR primers. The
pool-synthesis handles are amplification sites; they are not assigned sample-index
meaning merely because they match Illumina sequence components.

🟡 Selected type-IIS cut coordinates derive the transferred insert from each printed
pool design and give `{{= in4mer.LEFT_END }}`/`{{= in4mer.RIGHT_END }}` ends. A separate
restriction-site implementation checks those terminal sites. Unknown guide windows
leave additional internal sites unresolved; this is not a full digest proof.
Each selected insert core retains all four guide windows and the distinct intervening
repeats. No array deduplication
is performed. The constructor rejects reordered positional identities.

Recipient end geometry is independently unassigned. The ligation view therefore
shows the known insert beside an unresolved opened recipient rather than inventing
compatible receiving ends. Intended product junctions are inferred by default and
marked in preamble, affected captions and bases. No complete cohesive-product
relation is claimed for those assemblies.

The downstream biological assay and genomic recovery form an explicit unmodelled
boundary. It does not establish a molecular transformation, delivery, guide processing,
editing outcome or recovered genomic genotype. The genomic-PCR outputs are intended
target classes; their complete vector/genomic binding coordinates remain unassigned.

🟡 Prototype custom Read1 context is inferred from its reported sequencing primer.
Within that declared context, the {{= in4mer.PROTOTYPE_CYCLES['Read 1'] }}-cycle read
covers all four ordered {{= in4mer.GUIDE_LENGTH }}-base guide windows. Index sites
are computed from the matching canonical PCR arms and use the source's reported
eight-cycle index reads. This geometry does not establish activity or representation.
The source also reports a 30% custom random-library sequencing spike; its identity
and molecular composition are not expanded in the target-library view.

The vector-specific workbook rows are candidates for the matching Inzolia recipients,
not a resolved experiment assignment. This is marked on those primer materials,
derived libraries and the public page. The literal `[s]` is preserved without
assuming a stagger, modification, sequence or length. Local source-dependent checks
compare all four literal candidate rows, separately from their unresolved use.
Potential canonical sequencing sites can be shown, but neither acquisition nor
guide-cycle coordinates are assigned to Inzolia from kit capacity or unknown `[s]`.
No absence claim follows from this unresolved source interpretation.

Follow-up targets:

- Reconcile the repeated workbook heading and `[s]` with an original experiment definition.
- Read original vector maps and source-specific cassette boundaries before claiming
  exact recipient geometry, genome binding sites or full-map digestion.
- Assign Inzolia's actual acquisition recipe from an applicable original; the
  200-cycle kit capacity is insufficient on its own.
- Expand prototype's reported custom random-library spike only once its source
  identity/composition are assigned. No mixture balance is inferred here.
- The 7mer validation branch remains separate and outside this four-guide schematic.
- Supplementary Information original figure captions were inventoried for follow-up;
  no additional vector geometry or `[s]` meaning was taken from them. Full figure
  reading and independent source review remain open.
- Validate declared molecular relations separately from source transcription and
  experimental abundance. No representation or editing-efficiency guarantee follows
  from a correctly assembled example sequence.
