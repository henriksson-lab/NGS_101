# Record-seq

Schmidt, Cherepkova and Platt(2018), *Transcriptional recording by CRISPR spacer acquisition from RNA*, [published article](https://doi.org/10.1038/s41586-018-0569-1).

## Source scope

🟢 The defining Nature original is available through a university mirror, with its printed title, DOI, authors, Nature562:380–385 and2018 copyright identity checked. Selected Methods1385–1458 and ExtendedDataFigure2/text2168–2198 were read. Published SupplementaryInformation MOESM1 supplies Notes9–12, PrimerDesignNotes1–2 and Tables4–7. Original PDFpages14/29 (printed13/28) were rendered to verify the duplicate reverse-primer schematic and FS_0911/0963/0964. Two originals close this bounded pass; no later author-protocol revision or cited predecessor was substituted.

The mandatory acquired-bundle scraper was run after fetching the university original; its pointers to Table6 DR arms and Table7 adapters/primers were followed in the original. Reporting form, functional-screen workbook and stress-response workbook are cached but their full contents remain unreviewed. No confirmed missing-information or complete independent source audit.

## What is recorded

🟢 FsRT–Cas1/Cas2 enables acquisition into a plasmid CRISPR array. RNA-only splice-junction and RNA-virus evidence supports RNA-derived acquisition, but published SI explicitly does not exclude acquisition from DNA. The model therefore does not label every acquired spacer as RNA-derived, map a fixed RT-before-integration pathway, or predict sequence preference or acquisition efficiency. Source authors discuss semi-quantitative population recording, rare multiple acquisitions and limited temporal information. This is not an assured chronological transcript record, a single-cell lineage barcode or a UMI.

The selected SENECA route uses pFS_0235, a minimal FsArray2 recorder with leader, one direct repeat and FaqI-site context. Recording is induced in source BL21(DE3)Gold,3mlTB/100µMIPTG,37°C300rpm, typically12–14h to saturation. Two millilitres are harvested for standard plasmid miniprep/PicoGreen. Plasmid construction, orthologue screening, alternative host conditions and separate RNA-seq validation remain outside this selected route.

## Source library and molecular selectivity

🟢 Source FS_0963 and FS_0964 are annealed95°C5min then cooled to20°C at0.12°C/s and diluted1:100 inTE. Their actual source strings are represented using imported canonical sequences: FS_0963 is the exact source read2 arm (the canonical sequence without its finalT); FS_0964 is the source4nt `AAAG` overhang plus the complementary core. Scene computes antiparallel pairing rather than hand-placing columns. Source phosphates, complete plasmid cut geometry and ligated junction remain unassembled; a compatible donor duplex alone does not prove ligation.

Forty femtomoles plasmid, FaqI, T7DNA ligase, SAM/ATP/DTT and annealed adapters undergo99cycles37°C3min/20°C3min, then55°C15min. The published figure explains selectivity: digestion/ligation leaves a truncated DR in both parental and expanded arrays; an expanded array also retains an intact DR and new spacer. The parental class remains eligible on only the adapter side, while the expanded class can present both primer sites. The local ArrayClass compares this source regional eligibility only, conditional on the depicted cut context and intact copied repeats. It does not digest arbitrary spacer sequences or assess additional acquired restriction sites. It does not simulate PCR kinetics, perfect enrichment, loss of parental DNA or a guaranteed maximum spacer count.

FirstPCR uses source forward stagger family FS_0968–0974 and universal FS_0911 with NEBNextHighFidelity mix for22cycles98/57/72. Source2–8N staggers provide sequence diversity, not molecular identity. The forward arm/stagger/tail is one typed incomplete primer inventory, not three free oligo reagents; exact joined native template/primer geometry remains unassembled. Source firstPCR pooling is retained with unassigned membership, followed by1:1.6 magnetic-bead cleanup into25µlTE. Both source array classes and other carryover remain possible before actual separation; no perfect bead partition is assumed.

SecondPCR adds P5/BC1 and P7/BC2 sample-index families (8nt each, published Table4), using5µl input in20µl and8cycles98/65/72. Indexed sample pooling precedes QIAquickPCR desalting and2%E-gelEX/QIAquickgel extraction. The211–300bp and280–350bp windows belong to the separate classical readout and are not borrowed into SENECA. Source MiSeq/NextSeq500 uses150cycle kits,1.4–1.6pM loading and5–10%PhiX; a kit cycle count does not establish single-read length, paired-end configuration or index acquisition. All remain explicit read-plan limitations.

## Published ambiguities retained

- Main Methods specifies a10µl adapter ligation but10.38µl ligated input for one firstPCR. No volume, dilution or extra reaction is silently invented to reconcile this.
- Main cites PrimerDesignNote3, while the published SI lists only Notes1–2; Table6 identifies DR primer-binding arms. The crossreference remains unresolved.
- SI Note1 prints the same P5-arm string beneath both forward/P5 and reverse/P7 schematics. The model uses actual Table7 FS_0911/0963/0964 and Note2 secondPCR families; it does not silently replace the printed reverse schematic or impose it on every reverse primer.

Full recording intermediates, exact plasmid/substrate/end chemistry, complete firstPCR primers/products and physical index/read geometry remain unassembled. No source absence or finished schematic follows from this pass.


## Bounded publication follow-up (2026-10-30)

Bounded corrections retain actual2mlculture allocation, PicoGreen analytical allocation and40fmolplasmid contribution with explicit remainingstocks. The source adapter dilution/1ulcontribution also retains untransferredstock. ATP/SAM/DTT are molecularGoldenGate inputs, not only conditions. The explicitly named firstPCR master mix exposes copresent polymerase/nucleotide/othercomponents from its same actual stock. The secondPCR source5ulfrom25ul is allocated separately; its extra supplies/enzyme identity and fresh-versus-carried provenance remain unassembled rather than inheriting a named master mix. Actual first/indexed pools have selected own and distinct other source contributions with membership/amounts unassigned, not one-product relabelling. Gel input, qPCRQC and libraryloading each retain source allocations/remainders; PhiX5–10percent is an actual molecularDNA spike. A fixed literal reviewed recording/adapter/SENECA inventory replaces runtime builder-accumulated coverage. The required unassembled item remains uncovered; published10ul/10.38ul, Note3/Notes1–2 and duplicated schematic ambiguities remain unchanged. No completepage marker, source obligation waiver or public promotion.
