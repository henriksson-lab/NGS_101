# LIGR-seq — research notes

Defining source: Sharma *et al.*, [10.1016/j.molcel.2016.04.030](https://doi.org/10.1016/j.molcel.2016.04.030).

## Evidence used

- 🟢 LIGR-seq means LIGation of interacting RNA followed by high-throughput sequencing.
- 🟢 AMT photocrosslinking retains base-paired RNA duplexes; RNase R is used to deplete uncrosslinked linear RNA.
- 🟢 Proximity ligation joins the two RNA partners before crosslink reversal and cDNA library construction.
- 🟢 Step order (Experimental Procedures, coauthor-hosted PDF text lines 520–545; Supplemental Experimental Procedures lines 798–853; Figure 1A legend lines 249–251; Results lines 212–223): 293T cells + 20 µg/ml AMT, 365 nm for 30 min on ice → TRI-reagent extraction, TURBO DNase I → Ribo-Zero Gold depletion (4 µg) → limited S1 nuclease digestion of 400 ng (30 min room temperature, phenol–chloroform stop) → circRNA ligase (Epicentre) on 50 ng: 2 min 85 °C denaturation, ice, ATP + ligase 1 h 60 °C → RNase R (10 min 37 °C), phenol–chloroform, ethanol → 254 nm 10 min on ice, ethanol precipitation removes AMT → SMARTer stranded library. **Ligation precedes RNase R.**
- 🟢 Library preparation (Supplemental lines 838–853): "SMARTER strand library preparation kit (Clontech)" with modifications; first-strand primer adds a 5′ universal sequence; the RT adds nucleotides that anneal a universal forward primer for second-strand synthesis; PCR with primers to both added sequences, the reverse primer carrying a barcode index; 20 cycles (replicate 1) and 17 cycles (replicate 2); AMPure SPRI at 1:2 sample:beads after RT and PCR; 6% TBE PAGE, products above 200 nt excised and eluted overnight; Bioanalyzer check.
- 🟢 RNase R also resists circular ligation products (Figure S1A legend, lines 695–698), and Aligater detects circular products as concatemers of the two RNAs (Figure S1C, lines 708–710).
- 🟡 Model correction (2026-10-17): the previous schematic enriched with RNase R before ligation and drew a dual-indexed TruSeq shell with a dA ligation junction and 8-nt i5/i7 windows. None of that is supported: the source ligates first, and the cited kit is a SMART template-switching PCR library "without the need for … adapter ligations" (2013 components list). The library is now a role endpoint with the random barcode and kit index as structured features.
- 🟢 The paper describes a first-strand primer that adds a universal sequence, terminal-transferase activity followed by a universal second-strand primer, and indexed PCR. Its analysis explicitly trims a 5′ random barcode from single-end strand-specific reads.
- 🟢 The available 2013 Clontech *Product Components List* for catalog numbers 634836–634839 says the kit includes SMARTer Stranded RNA-Seq components and an Illumina Indexing Primer Set for Illumina-compatible libraries. Its second page points to a separate user manual.
- 🟢 An online transcript of the 2014 Clontech *SMARTer Stranded RNA-Seq Kit User Manual* (revision 050714, pp. 6–8) marks its SMARTer Stranded Oligonucleotide and N6 primer as proprietary sequences. It prints the twelve six-base Illumina Indexing Primer Set barcode values, including `ATCACG` for index 1. The original vendor PDF is unavailable in this source bundle. A 2019 Takara revision confirms the same kit family and lists a distinct HT index set; neither revision identifies which kit generation or index the LIGR-seq study used.
- 🔴 Model-unassigned: full native SMARTer/N6 oligos, the selected historical index primer and internal junctions. This is not a confirmed disclosure-absence finding; linked original vendor scopes remain unresolved.

## Protocol boundary

The page retains LIGR-seq's RNase R selection rather than presenting a generic RNA-duplex family schematic.

The schematic shows the library by role: SMARTer forward-primer end, the paper-reported 5′ random barcode, the two chimera partners, and the indexed reverse-primer end. Exact kit priming junctions, random-barcode length/origin, the selected index and the acquired read length remain unresolved; the scraper's only sequence hits in the article PDF (lines 1068–1090) are RNase-protection and RT-qPCR validation primers, not library oligos.

## Bounded source pass 2026-10-17 (round breadth-first-12)

Inspected: coauthor-hosted article + Supplemental Experimental Procedures PDF (Results lines 210–230, Figure 1A legend 249–251, Discussion 483–498, Experimental Procedures 520–556, Figure S1 legend 690–725, Supplemental Experimental Procedures 798–866); 2013 Clontech components list for 634836–634839 (whole file). Publisher full text remains a manual download (manifest). Not followed: 2014/2019 SMARTer manual originals (both manifest rows HTTP 403), RNase R and psoralen predecessor references, GEO GSE80167 run metadata for read length.


## Native graph migration — 2026-10-18, batch02

This pass re-read the manifest-identified coauthor article plus complete supplement
(`Sterne-Weiler_MolCell_2016.pdf.txt`, Methods520–545, FigureS1A/C690–710,
Supplemental Experimental Procedures798–866). The bundle scraper found validation
oligos, whose source context was checked; they are not assigned as library reagents.
No vendor manual or predecessor was newly inspected.

The module now owns its metadata, molecular states, endpoint and reviewed-target list.
Unknown RNA lengths and paired extent are role tokens; no shared 24-base partner
recipe remains. The diagram distinguishes a representative linear chimera and the
source-reported circular class without inventing a closure sequence. DNase/RNaseR
hydrolysates and RT RNA-template recovery remain molecular or unresolved rather
than automatic physical discards. The PCR uses the actual cleaned cDNA state.
The source reports single-end reads and random-barcode trimming; index acquisition,
read length, native kit oligos and first-cycle copying remain unassigned.

**Published size wording:** Supplemental Library Preparation says the 1:2
sample:beads ratio obtains fragments **<100bp**, while the next paragraph selects
products **above200nt** by PAGE. The `<100` mark was verified in the rendered
original PDF page21; the above200nt text is on page22. The public schematic shows
both locations and values. The model uses the printed bead ratio without claiming
a SPRI size cutoff and retains the explicit final PAGE threshold. Whether the
<100bp phrase describes retained or excluded material remains unresolved.

This is a bounded ownership/native-graph migration, not a complete source review
or acceptance of missing information. Deferred targets remain the source-specific
kit revision/primers, read/index acquisition, RNaseR/psoralen linked originals,
and complete native copying/ligation/recovery proof.
