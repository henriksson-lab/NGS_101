# Duplex Proximity Sequencing (Pro-Seq)

Duplex Proximity Sequencing (Pro-Seq): A method to improve DNA sequencing accuracy without the cost of molecular barcoding redundancy. PLoS ONE (2018), [doi:10.1371/journal.pone.0204265](https://doi.org/10.1371/journal.pone.0204265).

Source scope: 🟢 published main PMC6168144, Results70–99, Methods194–226 and S2Fig caption284–286. S1 workbook header/table coordinates and S2/S3 workbook headers were administratively inspected, totaling FOUR originals including main. This exceeded the three-original investigation cap; investigation stopped, and this is not a successful strict capped pass. Mandatory scraper reviewed (6 candidates in main PDF/XML). S2 architecture TIFF retrieval was gated; its original was not inspected. No complete source or absence audit.

🟢 DNA Pro-Seq is separate from RNA precision-run-on PRO-seq. The selected route partitions non-denatured dsDNA so both original senses can remain together. Multiplex gene-specific PCR and5′PEG-linked universal primers generate physically linked copies that seed one optical cluster. The target is mixed-cluster signal redundancy, not a tandem DNA concatemer or a barcode-family computational collapse. The paper says only about50% of two-linker complexes represent both senses; duplex classification uses custom index evidence after sequencing. No universal two-original-strand guarantee is modeled. The manuscript's whole-genome PEG-loop/strand-displacing single-extension scheme is a distinct, unmodeled route.

🟢 Plasma is spun2000g10min4°C; QIAampCirculatingNucleicAcidKit extracts DNA with2×100µl0.1×TE10min/20000g3min elutions. Monarch cleanup elutes15µl0.1×TDTE; material not used within24h is frozen−20°C. COG5/ALB qPCR serial dilution/standard curves estimate input genomic copies. The source40µl droplet mix contains Q5HotStart0.02U/µl,dNTP0.2mM,Q5buffer,RDTstabilizer,index1/index2 gene-specificforward primers each25nM,gene-specificreverse50nM,universalreverse400nM and linked universalprimer200nM. RainDropSource/ThunderBolts produces droplets; printed approximate count is `8,000,0000`, retained without silently repairing typography. Occupancy is statistical.

🟢 Universal reverse primer as written5′→3′ (source spacing retained):

```
CAAGCAGAAGACGGCATACGAGAT GCGACGGTTAGACGAACGGTACG
```

🟢 Linked ordering as printed, with the first arm written3′→5′ and a5′–5′ `/Sp9/` connection:

```
3'-AGTCAACGTCGTCTTCTTCCACATCTAGAGCCACCAGCGGCATAGTAA-5'/Sp9/5'-AATGATACGGCGACCACCGAGATCTACACTCCCTCCTATCATGGACAC-3'
```

The first arm must be reversed, not reverse-complemented, when expressing both arms5′→3′. The PEG connection is not a nucleotide/phosphodiester join. Exact gene-specific and custom sequencing primer inventories are unassembled; the source S2 original remains relevant.

🟢 PCR98°C30s;38cycles98°C10s/60°C30s/72°C30s,1°C/s ramp,4°Chold. Emulsions are broken with62.5µl destabilizer;0.8:1AMPure elutes20µl. T7exonuclease6.7U+RecJf41.7U in50µlNEBuffer4 digests unlinked DNA37°C1h/70°C20min.1.6:1AMPure elutes25µl. Cleavage is not physical removal until this cleanup, and neither digestion nor purification is assumed perfect. Standard amplicon controls use unlinked primers and omit digestion; they are not the selected route.

🟢 KAPALibraryQuantKit QC precedes loading:18µllibrary+2µl1NNaOH5minRT/ice,5%denaturedPhiX,2µl1NHCl,HT1to600µl/5.5pM. MiSeq acquisition uses overlappingR1/R2 lengths selected per amplicon. CustomR1 is a1:1 index1/index2primer mix plus standardR1; each index read uses its custom primer plus standardR1; customR2 is3.5µl100µM inwell14. Index2 is sequenced before turnaround and dark cycles omitted by XMLchanges. Intensity files, non-pass-filter clusters and indexFASTQs are retained. Doubly-seeded checks and mixed-base masking are computational, not physical library selection.

Open: S2TIFF, exact whole linked-complex architecture, target/read-primer identities, predecessor/vendor originals and whole-genome route. No absence claims, whole review or source-qualified parent relation established in this bounded pass.
