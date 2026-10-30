# m1A-seq

Dominissini et al.(2016), *The dynamic N1-methyladenosine methylome in eukaryotic messenger RNA*. DOI https://doi.org/10.1038/nature16998 .

## Sources and version

🟢 Published SI MOESM118 Notes1–3 pp1–2 inspected first, especially Note2 preparation/control logic and Note3 approximate trough resolution. Associated accepted manuscript PMC4842015 is explicitly identified by provenance72; Methods144/160 and mismatch analysis184 inspected for scoped procedural detail. The main published Methods were not independently reconciled with this manuscript. Its procedures are AM-qualified, not silently promoted to VOR facts. The publisher-mirrored published SI is distinct from the gated NIHMS supplement route. Mandatory full-bundle oligo scraper run, no sequences found in the selected bundle; this is not an absence conclusion.

## Source-owned chemistry

🟢 Published SI Note2 establishes RNA-fragment immunoprecipitation, competitive elution using free m1A mononucleoside, then the alkaline m1A→m6A rearrangement control before RT. This is partial rearrangement and reduces stops/mismatch signatures; RNA m1A is rearranged, not removed or deleted from RNA. The source also describes native readthrough with mismatches. Peak centers/troughs approximate positions; a source trough is not guaranteed to identify an exact single chemical residue.

🟢 Associated AM Methods160: fragmentation70°C15min followed by RNA Clean&Concentrator column; anti-m1A pre-coupled ProteinG Dynabeads3h4°C. Two alternative elutions are stated: ProteinaseK/TRIZOL/ethanol, or free-mononucleoside competitive elution. The local graph selects the latter, retaining residual bead-associated material explicitly; it does not blend the alternatives.

The AM rearranges input/IP fragments in50mMNa2CO3/2mMEDTApH10.4,1h at50°C or60°C alternatives; no single temperature is assigned to every library. Column cleanup follows. Historical NEBNext Ultra RNA Library Prep for Illumina supplies first-strand/library processing, but exact edition/internal strand steps/oligos remain unassembled. AM readout names HiSeq2500/TruSeqSRv3; exact cycle count/index/run-primer geometry remains unreviewed. Whole-mRNA expression and PARS/ribosome analyses in the same article are separate experiments, not extra operations in this library.

🟡 The selected-site local `dimroth(positions)` transforms m1A to m6A at caller-specified residues. It keeps the modified methylated nucleotide, coordinate count and untreated template identity. Complementary firstDNA is derived separately; native m1A readthrough requires an explicit caller-selected copying outcome. It does not encode a universal mutation distribution or complete conversion.

## Boundary and lineage

Distinct from Li2016 m1A-ID-seq demethylation/AMV comparison. Published SI explicitly calls m1A-seq modeled after m6A-seq, supporting an explicit component-strategy lineage to the cited m6A-seq predecessor, not inheritance of all adapters/end chemistry or a later kit edition. That predecessor and historical kit originals remain deferred. Main VOR applicability, exact full library geometry and controls remain open; no absence claim or finished schematic.
