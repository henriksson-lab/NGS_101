# BAsE-Seq

Hong et al., *BAsE-Seq: a method for obtaining long viral haplotypes from short sequence reads*, Genome Biology (2014), [defining published article](https://doi.org/10.1186/s13059-014-0517-9).

## Sources inspected

🟢 Published PMC4269956 JATS: identity/licence, Results68–72/Figure1 overview, Methods346–356, sequencing/processing360–364, and separate Deep-Seq378. Published Additional file2, **last updated2July2014**, overview, oligo tablep3 and all Steps1–9 (text225–707); AppendicesA–C describe barcode carryover, amplification chimeras and dsDNA-circularization optimization. The rendered original tablep3 was inspected to verify underlined arms, the red internal-standard AT insertion, internal biotin-dT and phosphorothioate notation. Additional file1 FigureS1/caption explains barcode assignment versus later error removal; other result tables were not inspected. All three are source-qualified originals, acquired through PMC Cloud with assigned manifest URLs and checksums. The mandatory scraper returned28 hits across four extracted files; the original detailed table, not the scraper parsing of `/iBiodT/`, establishes that physical modified T.

This is a capped source pass, not a complete disclosure/linked-source audit. The historical NexteraXT Sample Preparation Guide is explicitly cited by the detailed recipe but has not been version-resolved/read. Available source facts are distinguished from product geometry not yet assembled. See [retrieval instructions](ref/MANIFEST.md).

## Distinct chemistry and input scope

🟢 Source long HBV DNA molecules receive a random genome barcode. Exonuclease-mediated nested deletion and intramolecular **dsDNA** ligation bring different portions of each amplified genome adjacent to that barcode; NexteraXT then fragments/tags the circular preparation. Short reads are grouped computationally by barcode to assemble viral haplotypes. This differs from independent droplet/bead co-barcoding and does not require a modified long-read sequencer. The same barcode is intentionally shared by amplified descendants: it is a linked-read genome identifier, not a fresh tag for every PCR copy.

🟢 The graph starts at full-length HBV genomes. Upstream plasmid construction, serum extraction and viral-load qPCR are documented in the article but outside this implementation. The separate Deep-Seq comparator uses another primer/library/read route and is not silently merged into BAsE-Seq. An internal-standard assignment primer carries an additional AT immediately before its random window; the standard is separately barcoded before source pooling.

## Source oligos

These are source ordering strings. Asterisks mark phosphorothioate bonds and `/iBiodT/` is a physical modified thymidine, not a nucleotide-free spacer. Canonical flow-cell and Nextera sequence components are imported; no current kit reaction recipe is substituted for the historical source.

| Source name | 5′→3′ ordering string |
|---|---|
| A-SbfI-ID20-HBV_For | {{= base_seq.assignment_primer().top() }} |
| A-SbfI-ID20v2-HBV_For | {{= base_seq.assignment_primer(internal_standard=True).top() }} |
| B-HBV_Rev | {{= base_seq.UNIB + base_seq.HBV_R }} |
| Uni-A-SbfI-Bio | {{= base_seq.UNI_A_BIO.ordering() }} |
| Uni-B | {{= base_seq.UNIB }} |
| NexteraXT-PE1 | {{= base_seq.FINAL_F }} |
| Illumina-PE2_HBV-For.RevComp | {{= base_seq.FINAL_R }} |
| Illumina_PE2 SVS.SeqPrimer | {{= base_seq.R2_PRIMER.ordering() }} |

The local modified-primer type preserves the internal biotin position and each PS bond between physical bases. A polymerase-made complementary copy is not assigned copied biotin or PS chemistry. The entire resulting biotinylated duplex is not yet assembled; the source reagent is an explicit physical input.

## Source operations

1. 🟢 Two Long-PCR cycles assign the20random-base barcode. ExonucleaseI digestion is followed by heat inactivation and1xAMPure recovery. Digestion does not itself remove soluble products or prove complete primer removal. The source warns that carryover assignment primers can assign further barcodes later.
2. 🟢 Two clonal-amplification stages use Uni-A-SbfI-Bio/Uni-B. Separate real-time reactions choose log-linear cycle counts; PCR1 is cleaned before distributing PCR2 among eight wells, pooling and cleaning again. No fixed universal cycle count, chimera-free copying, collision-free barcoding or unbiased recovery is inferred.
3. 🟢 SbfI-HF digestion generates a barcode-proximal4nt3′overhang that inhibits proximal ExonucleaseIII deletion. Source1.8xAMPure cleanup precedes the ExoIII time course. Aliquots at1/3/5/7/9minutes are quenched in EDTA, heat-inactivated, pooled and purified. S1 nuclease then processes exposed single-stranded regions, followed by EDTA quench and another cleanup. The heterogeneous full genomic ends and substrate-specific efficiencies remain unassembled.
4. 🟢 Streptavidin Dynabeads capture biotinylated DNA for15h. Magnet/washes remove the unretained fraction. **Two PBI dissociations at50°C**, each followed by magnetic separation, transfer DNA-containing supernatants onto the same MinElute column. The model does not retain all DNA permanently on the beads, discard the released fraction, or infer perfect biotin-class purity.
5. 🟢 T4DNApolymerase/T4PNK treatment supplies blunt ending/5′phosphorylation, then cleanup. T4DNAligase circularizes at10°C16h using≤100ng in45ml. The intended intramolecular dsDNA circles coexist with uncircularized DNA and possible intermolecular products. Appendices describe abandoning the less effective CircLigaseII ssDNA route; that is not the represented circle chemistry.
6. 🟢 Precipitation/washing/recovery precedes Lambda-exonuclease plus ExonucleaseI digestion. Products remain in the reaction until AMPure cleanup; SpeedVac concentrates the recovered preparation. Circle enrichment is not a complete purity or yield prediction.
7. 🟢 Historical NexteraXT transposomes fragment/tag the preparation. Source NPM PCR uses the custom primers rather than the kit's index primers, beginning with72°C3minutes and using the source14cycle program unless input is low. A0.65xAMPure size-selection step precedes quantification and sequencing. The source historical kit formulation, exact transposome ends, duplicated-target geometry and full final genomic strands are not reconstructed yet.

The graph records physical mixture/workup order, with separate retained/rejected capture and cleanup outputs. Soluble products persist until actual source separations. This does not establish every molecular claim: exact end/product assembly and linked-source coverage keep the graph incomplete and its preview private.

## Sequencing and computational fates

🟢 MiSeq paired-end150 is reported. The custom second-read primer anneals to the HBV forward priming sequence so the read starts at the genome barcode. A Scene computes that local body/arm binding relation; it is not a declaration of the entire final library. Exact physical final primer/index acquisition geometry remains unassembled. No additional sequencer sample index is invented from the presence of a NexteraXT Index Kit in the reagent list.

🟢 Read grouping, duplicate filtering, haplotype assembly and mismatch/error thresholds are computational operations. They do not physically remove molecules during preparation. The20random bases do not guarantee uniquely assigned genomes, and assignment-cycle errors can remain shared among all descendants of one barcode. No source sensitivity, error rate or chimerism estimate is turned into a constructor invariant.

## Remaining work

Full HBV-specific assignment/amplification strands, internal modified-primer propagation, SbfI/ExoIII/S1/end-repair product geometry, circle closure coordinates, historical NexteraXT guide identity/applicability, final library/read/index geometry, uninspected result tables and an independent complete original/linked-material review remain open. The available documents support preparative chemistry reading; they do not establish confirmed missing information or a finished schematic.
