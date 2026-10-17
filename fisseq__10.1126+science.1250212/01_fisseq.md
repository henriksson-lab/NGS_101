# FISSEQ

## 1. What it is

FISSEQ converts cellular RNA into immobilized cDNA circles, copies each circle into a
localized rolling-circle amplicon (a rolony), and reads the repeated cDNA directly in a
fixed cell or tissue by fluorescent sequencing-by-ligation.

Defining paper: Lee *et al.*, “Highly multiplexed subcellular RNA sequencing in situ,”
*Science* (2014), doi:[10.1126/science.1250212](https://doi.org/10.1126/science.1250212).

## 2. Sources read

- 🟢 The defining Science paper's open author manuscript, PMCID **PMC4140943**.
- 🟢 Lee *et al.*, “Fluorescent in situ sequencing (FISSEQ) of RNA for gene expression
  profiling in intact cells and tissues,” *Nature Protocols* (2015),
  doi:[10.1038/nprot.2014.191](https://doi.org/10.1038/nprot.2014.191), PMCID
  **PMC4327781**. This is the authors' detailed protocol and the source for exact oligos.
- Bounded pass 2026-10-17 inspected: Science JATS (`science.1250212_PMC4140943.xml.txt`)
  abstract line 76, Results lines 95–160 and Fig. 1–2 legends lines 372–393; Nature Protocols
  text (`PMC4327781.html.txt`) procedure overview and Fig. 1 legend lines 1286–1330, Box 1
  lines 1745–1765, reagent oligos lines 2035–2083, and Steps 2–17 and 22–28 lines 2334–2860.
- 🔴 The 2014 Science supplementary material (figs. S1–S10, methods) remains unread because
  the PMC download is gated (see `catalogue/source_reviews/`).

## 3. Reverse transcription and immobilization

🟢 The random-hexamer RT primer is written
`/5phos/TCTCGGGAACGCTGAAGANNNNNN`: an 18-nt common adapter followed by six random
bases. M-MuLV reverse transcriptase copies RNA in situ while aminoallyl-dUTP introduces
primary amines into the cDNA (2.5 µM RT primer, 40 µM aminoallyl-dUTP; 10 min at 4 °C, then
37 °C for 1–2 h to overnight; Steps 2–3). BS(PEG)9 then cross-links that cDNA to the
cellular protein matrix (1 h, then 1 M Tris quench; Steps 4–5). RNase and RNase H remove
the RNA before circularization (1 h at 37 °C; Step 6).

## 4. Direct cDNA circularization

🟢 CircLigase II acts for one hour at 60 °C with 2.5 mM MnCl2 and 0.5 M betaine (Steps 8–9). It joins the cDNA's free 3′-OH directly to
the RT primer's 5′ phosphate, producing a covalently closed single-stranded cDNA circle.
This is direct intramolecular end joining, not padlock capture and not plasmid cloning.

🟡 `lib/circular.py` represents that topology explicitly. `circularize_ssdna()` refuses
unphosphorylated or 3′-blocked input and records the last-to-first closure. The schematic's
circle and ligation marker are generated from that boundary.

## 5. Rolling-circle amplification

🟢 The RCA primer is `TCTTCAGCGTTCCCGA*G*A`; the final two internucleotide bonds are
phosphorothioates. Its base sequence is exactly the reverse complement of the RT adapter.
After 500 nM primer hybridization for 1 h at 60 °C (Steps 11–12), Phi29 DNA polymerase amplifies overnight at 30 °C with
aminoallyl-dUTP. The amine-bearing rolonies are cross-linked again.

🟡 `rolling_circle()` locates the primer on the circle and derives one complete
complementary repeat. The displayed concatemer repeats that computed unit; neither the
repeat sequence nor its boundary is entered independently.

## 6. Sequencing-by-ligation

🟢 Five 5′-phosphorylated sequencing primers bind the adapter-complement portion of each
RCA repeat. Primer N is `TCTCGGGAACGCTGAAGA`; N-1 through N-4 remove one successive base
from its 5′ end. Seven rounds of fluorescent probe ligation and cleavage are performed
per offset so the five primer starts cover intervening positions. Each round uses T4 DNA
ligase with SOLiD ToP eight-base probes for 45 min, imaging, then Cleave Solutions 1 and
2.1 (Steps 23–28). The first ligation image for recessed primers N-2 to N-4 is excluded.

🟡 Ligation direction (model correction 2026-10-17). The Fig. 1b legend says fluorophore
cleavage leaves "a phosphorylated 5’ end at the ligation complex suitable for additional
ligation cycles". All five primers carry /5phos/ and share the same 3′ end, so the offsets
work only if probes are added at the primer's 5′-phosphate end. The earlier schematic drew
a 3′-extension arrow; it now marks 5′-end ligation. On the drawn rolony, that end faces the
transcript copy, which is consistent with the Science paper's mCherry transcription-start-site
reads (Results lines 137–146). This direction is derived; it is not stated as a strand map.

🟢 Published discrepancy shown on the page: Box 1 lists start positions 1-6-11-16-**20**-26-31
for primer N and 2-7-12-**16**-22-27-32 for primer N-4, while primer N-1 is 0-5-10-15-20-25-30
and the text describes a five-base stride. A strict stride predicts 21 and 17. The model draws
no per-position read map, so it does not choose between them.

🟢 Read length: the Science abstract says "30-base reads from 8,742 genes" (line 76), while
the fibroblast Results report "sequencing reads of 27 bases" for 4,171 genes (line 150);
the Nature Protocols sequencing section is timed "for 30 cycles". These may describe
different datasets. No acquired read length is modelled.

🔴 Not modelled: the Science paper's partition sequencing with pre-extended sequencing
primers (Results lines 118–124; Fig. 2A). The extended primer sequences are not printed in
the inspected text.

This is SOLiD-style imaging chemistry on a rolony, not an Illumina flow-cell library. The
schematic therefore shows these five published sequencing primers at their actual repeated
binding site rather than inventing Read 1, index, or Read 2 primers.

🟢 The detailed protocol (PMC4327781, reagent list) also prints three fluorescent
adapter-specific probes with this same 18-base sequence and 5′ `/56-FAM/`,
`/5TYE563/` or `/5TEX615/` labels, plus four 18S and seven 28S rRNA detection
primers with 5′ biotin. Those are ancillary imaging and control oligos; they do
not become part of the sequenced cDNA circle or its RCA product.
