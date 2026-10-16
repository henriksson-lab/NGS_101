# Small-seq sequencing evidence and model scope

## Published versions

The schematic models Hagemann-Jensen et al. 2018, DOI
[10.1038/s41596-018-0049-y](https://doi.org/10.1038/s41596-018-0049-y).
Its Reagent setup prints the modified adapters and PCR primers; Supplementary
Table 1 supplies 192 complete SRX primers with eight-base indexes. Steps 37–38
refer instrument operation to the manufacturer's protocol and recommend a single
read of at least 51 bases.

The original Faridani et al. 2016 paper, DOI
[10.1038/nbt.3701](https://doi.org/10.1038/nbt.3701), uses six-base RPI indexes
and reports 43- or 51-base single reads (Online Methods, library preparation).
Its Supplementary Table 3 distinguishes 33- and 41-base RNA portions after the
ten-base UMI/linker removal. These are distinct published implementations.

## What the model establishes

🟢 The 2018 adapter-ligation section specifies sequencing from the UMI, followed
by CA and the RNA insert. The original 2016 supplementary software README
extracts eight UMI bases and trims two additional bases; its adapter FASTA lists
`TGGAATTCTCGGGTGCCAAGG`, matching the RA3 body.

🟡 The displayed Read 1 entry is therefore a binding-site query spanning RA5
and ending immediately before the UMI. It is not an identified loaded sequencing
oligo. RP1 ends five bases earlier and would read those bases before the UMI.
The canonical TruSeq DNA Read 1 sequence has no exact binding site on this
modeled library. These calculations concern this construct, not instrument
performance with every possible reagent mixture.

🟢 The SR001 indexed PCR primer and its separately reported index read establish
the eight-base index geometry. An index-primer binding calculation does not by
itself document the authors' loaded primer mixture.

## Remaining source work

The exact loaded Read 1 primer is outside the present model's established
source scope. This does not establish that the oligo was unpublished. The
manufacturer instructions referenced by Step 37 still require a scoped source
check before accepting any absence claim. The per-source checklist and two
independent review passes are maintained in
`catalogue/source_reviews/small-seq__10.1038+nbt.3701.json`.

The 2018 protocol was recovered from a pre-existing ignored publisher PDF and
migrated into the source bundle. Calling it inaccessible in earlier audit records
was a repository source-discovery failure, not a publication omission.
