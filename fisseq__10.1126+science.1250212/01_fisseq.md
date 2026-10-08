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

## 3. Reverse transcription and immobilization

🟢 The random-hexamer RT primer is written
`/5phos/TCTCGGGAACGCTGAAGANNNNNN`: an 18-nt common adapter followed by six random
bases. M-MuLV reverse transcriptase copies RNA in situ while aminoallyl-dUTP introduces
primary amines into the cDNA. BS(PEG)9 then cross-links that cDNA to the cellular protein
matrix. RNase and RNase H remove the RNA before circularization.

## 4. Direct cDNA circularization

🟢 CircLigase II acts for one hour at 60 °C. It joins the cDNA's free 3′-OH directly to
the RT primer's 5′ phosphate, producing a covalently closed single-stranded cDNA circle.
This is direct intramolecular end joining, not padlock capture and not plasmid cloning.

🟡 `lib/circular.py` represents that topology explicitly. `circularize_ssdna()` refuses
unphosphorylated or 3′-blocked input and records the last-to-first closure. The schematic's
circle and ligation marker are generated from that boundary.

## 5. Rolling-circle amplification

🟢 The RCA primer is `TCTTCAGCGTTCCCGA*G*A`; the final two internucleotide bonds are
phosphorothioates. Its base sequence is exactly the reverse complement of the RT adapter.
After hybridization at 60 °C, Phi29 DNA polymerase amplifies overnight at 30 °C with
aminoallyl-dUTP. The amine-bearing rolonies are cross-linked again.

🟡 `rolling_circle()` locates the primer on the circle and derives one complete
complementary repeat. The displayed concatemer repeats that computed unit; neither the
repeat sequence nor its boundary is entered independently.

## 6. Sequencing-by-ligation

🟢 Five 5′-phosphorylated sequencing primers bind the adapter-complement portion of each
RCA repeat. Primer N is `TCTCGGGAACGCTGAAGA`; N-1 through N-4 remove one successive base
from its 5′ end. Seven rounds of fluorescent probe ligation and cleavage are performed
per offset so the five primer starts cover intervening positions.

This is SOLiD-style imaging chemistry on a rolony, not an Illumina flow-cell library. The
schematic therefore shows these five published sequencing primers at their actual repeated
binding site rather than inventing Read 1, index, or Read 2 primers.
