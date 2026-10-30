# XYZeq: spatial RT and PCR-well indexing

Lee et al. (2021), **XYZeq: Spatially resolved single-cell RNA sequencing reveals expression heterogeneity in the tumor microenvironment**. [Published defining article](https://doi.org/10.1126/sciadv.abg4755), [published PMC original](https://pmc.ncbi.nlm.nih.gov/articles/PMC8059935/).

Source scope: one published main original, JATS Results136, Methods210–222/238 and references17/18; its rendered PDF Fig1A/B on page3 was inspected because the extraction does not encode the library diagram. Other main sections, detailed supplements, predecessor originals and historical reagents are not fully audited. Mandatory oligo scraper found six candidates across the PDF/XML twins; the relevant source lines and rendered figure were checked. These are not six different oligos. The fetched ePMC supplementary ZIP contains only main figure images. The exact Supplementary Materials index printed in the main returned HTTP403. No SI filename was guessed, and no missing-information claim is accepted.

🟢 Unlike the separate 10x comparator, the selected XYZeq tissue preparation fixes and cryopreserves tissue, partitions tissue on spatially indexed microwells, performs RT in place, recovers intact cells, then sorts them into PCR wells for a second RT and indexing PCR. Source main136/Fig1B identifies a cellular index made from spatial location and PCR-well identity. The PCR-well i5 is **a cell-barcode component**, rather than an ordinary sample index.

🟢 Main214 supplies the RT primer, without additional modification notation. The main's rendered Fig1B labels the random region as UMI and the adjacent region as spatial barcode. Here the spatial code is a placeholder, not a fabricated whitelist member:

```
{{scene: chem.rt_primer()}}
```

🟢 Main222 supplies the P5–i5–R1 indexing-primer pattern. The model imports canonical fixed bodies; the printed `[i5]` stays a variable PCR-well code. The read's eight-cycle index acquisition does not by itself assign a concrete well code or separately prove its molecular length:

```
{{scene: chem.pcr_index_primer()}}
```

🟢 Selected source preparation, main210–222:

1. Harvested organ spends30minutes in ice-cold DMSO-free freezing medium, then30minutes in ice-cold DSP with10percent FCS. Neutralize in ice-cold20mMTrisHClpH7.5, seal in a cryomold and slowly freeze overnight at−80C. The diagram does not predict crosslink occupancy or reversal products.
2. Place a25µm tissue section on a preloaded, cooled microwell array. The source adds DSP-fixed HEK293T control cells to the RT/digestion mix; actual distribution across wells varies. Seal the chip and warm gradually to42C, then incubate20minutes. Exact intermediate temperature ramps, full cocktail and template-switch chemistry remain unassembled.
3. Recover the chip into50ml1xSSC with25percent FCS, vortex, centrifuge1000rcf10minutes and remove excess liquid. Retain the cells; recovery and carryover are not assumed perfect. Filter, stainDAPI and sort into96-well plates preloaded with5µl second RT mix. DAPI is a stain, not a sequence index. Uncollected sorting material is an unresolved physical fraction.
4. Second RT runs42C for1.5hours. PCR uses2xKapaHotStartReadyMix and the printed P5–i5–R1 primer. The unassembled partner primer, full RT mix and cycle program are explicit gaps in this selected scope; no sci-RNA-seq or SPLiT-seq constructor fills them.
5. Physically pool PCR plate contents into2mltubes and purify cDNA using AMPureXP SPRI beads. The private graph currently aggregates the well classes, so this is an explicitly incomplete physical-pool boundary rather than an exact validated pooling relation.
6. Tagment cDNA, amplify with an Illumina Nextera p7 index, assess BioAnalyzer and Qubit aliquots, then sequence. Historical tagmentation loaded oligos, exact strand ends, final PCR partner and workups remain unassembled.

🟢 Main222 reports NovaSeq or HiSeq4000 acquisition: R1=26, R2=98, I1=8 and I2=8cycles. Fig1B places P7/i7/R2 toward transcript, followed by oligo(dT), spatial code, UMI, R1, PCR-well i5 and P5. That source schematic is not a complete constructed duplex or computed primer-binding validation. The private graph does not invent exact read-cycle feature spans.

Main136 credits sci-RNA-seq and SPLiT-seq for combinatorial indexing. This is a strategy-level development relationship; it does not establish matching RT cocktails, template-switch oligos, ligation steps or all library chemistry. Main238 delegates computational read processing to Cao2017; DecontX, UMI deduplication, barcode assignment and expression filtering are computational operations, not physical removal of DNA from this preparation.

This is an incomplete private molecular graph with known oligo inventory and unresolved endpoints. Needed next: the actual published supplementary methods and relevant chemistry/oligo tables from the printed index, then directly cited original components where needed. No full source review, public schematic or accepted absence claim is asserted.
