# HDST — research notes

Vickovic et al., *Nature Methods* (2019),
[doi:10.1038/s41592-019-0548-y](https://doi.org/10.1038/s41592-019-0548-y).

## Sources read

- 🟢 Defining article and Online Methods (PMCID: PMC6765407).
- 🟢 Jemt et al. automated Spatial Transcriptomics library preparation,
  [doi:10.1038/srep37137](https://doi.org/10.1038/srep37137).

## Bead oligo and decoding

🟢 The published precursor is
`/AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-Spatial_barcode_Pool1-Bridge1`.
It contains an amine coupling modification, five uracils for later USER release, a T7
promoter, an Illumina handle, the first 15-base spatial block and its bridge.

🟢 Two more split-and-pool ligations add 14/15-base barcode pools and bridges, producing
65 × 211 × 211 combinations on 2-µm silica beads. Online Methods "Bead production"
(JATS line 121): the precursor's bridge is 7 bp; `/Pho/Bridge2-Spatial_barcode_Pool2-Bridge3`
is ligated to Bridge1 through a complementary 14-bp `Bridge1'Bridge2'` helper with T4 DNA
ligase (NEB, manufacturer's protocol), added 2:1 over the precursor; the product's Bridge3
then accepts `/Pho/Bridge4-Spatial_barcode_Pool3` through a `Bridge3'Bridge4'` helper.
Pool 3 is followed by the UMI and 20 poly(dT)VN. Pool oligos were made on Illumina's
"Big Bird" synthesizer. 🟡 Bridge2 is therefore the last seven printed bases
(`CTAGAGC` of `GACTTGTCTAGAGC`); the helper bases are not printed and the
Bridge3/Bridge4 split is not stated. The model draws only this round-2 splinted donor end. The original Online Methods,
"Bead production" (JATS line 121), prints the ligated bridge sequences as
`GACTTGTCTAGAGC` (Bridge1Bridge2) and `TGATGCCACACTACTC` (Bridge3Bridge4), and
states that a **five-base** UMI and 20-base poly(dT)VN capture region follow the spatial
address. The individual barcode and decoder sequences were an Illumina product and are
explicitly described as proprietary; this does not apply to the two printed bridges.

🟢 Beads are randomly loaded into more than 1.4 million hexagonal 2-µm wells. Fourteen
successive FAM/Cy3/dark decoder-hybridization images identify all three barcode blocks and
therefore map every bead to an x–y well before tissue is applied.

## Capture and library construction

1. 🟢 Place a 10- or 16-µm frozen section on the decoded array, fix, H&E stain and image.
2. 🟢 Permeabilize (exonuclease I buffer 30 min 37 °C, then 0.1× pepsin pH 1 10 min
   37 °C; breast cancer: 14 U collagenase I 20 min before pepsin) and synthesize cDNA
   overnight at 42 °C with SuperScript III (RNaseOUT, 1× FS buffer, 5 mM DTT, 0.5 mM dNTP,
   50 µg/ml actinomycin D, 0.19 mg/ml BSA; JATS line 141). The
   bead-bound poly(dT) primer covalently links spatial address, UMI and transcript cDNA.
3. 🟢 Proteinase-K digest the section (1 h, 56 °C) and release barcode–cDNA molecules with USER for
   2 h at 37 °C at the five-uracil tract.
4. 🟢 In the cited Jemt et al. Spatial Transcriptomics workflow, released cDNA undergoes
   second-strand synthesis, end repair and T7 in vitro transcription (JATS lines 49–53;
   Fig. 1 caption, lines 308–310).
5. 🟢 Jemt et al. then ligate a sequencing adapter to the amplified RNA, synthesize cDNA
   again and add sample indexes by PCR. Jemt's Methods (line 89) delegate the exact
   adapter and primer sequences to the supplementary information of Ståhl et al.
   (2016), doi:10.1126/science.aaf2403. That original supplement remains inaccessible;
   these steps do not establish the exact final HDST adapter junction. The graph keeps
   them as one warned unresolved operation before the indexing PCR.

## Reads

🟢 Libraries were sequenced at 1.08 pM on a NextSeq 500 with v2 chemistry (JATS line 141);
index-read lengths are not stated. Mouse olfactory-bulb libraries used R1 125 and R2 175; breast-cancer libraries used
150 + 150. Read 1 positions 1–15, 31–45 and 61–76 were retained as the three-part spatial
barcode, and positions 77–82 are the six-base UMI. Read 2 is transcript sequence.

🟡 The published Online Methods disagree internally. "Bead production" (JATS line 121)
specifies barcode pools of 15/14/15 bases, ligated bridges of 14 and 16 bases, and a
**five-base UMI**. "Sequence preprocessing" (JATS line 145) extracts Read 1 ranges of
1–15, 31–45, 61–76 and 77–82, which are 15/15/16 barcode bases and a **six-base UMI**.
The model uses the explicitly described oligo synthesis lengths and printed bridge
bases for its physical construct, while the page quotes the published read mask separately.
The coordinate offset is unresolved; neither source value is silently repaired. The
compatible final Illumina shell remains inferred because the cited library construction
does not give an exact HDST-specific end sequence.

🟡 The authors' earlier bioRxiv v2 Methods introduce a third UMI length: the bead
production passage (PDF text lines 220–235) says **seven bases**, whereas its own
preprocessing passage (lines 302–319) also selects positions 77–82, or six bases.
The model follows the final published bead-production value of five bases; neither
version establishes how the reported read coordinates align to that physical oligo.
