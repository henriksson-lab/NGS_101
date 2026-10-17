"""EMP ITS1f–ITS2 fusion-primer library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Earth Microbiome Project ITS amplicon sequencing"
NOTES = "01_emp-its.html"
SOURCE = '<a href="https://earthmicrobiome.org/protocols-and-standards/its/">Earth Microbiome Project ITS Illumina Amplicon Protocol (EMP.ITSkabir)</a> and its November 2016 oligo workbook.'
SUMMARY = "One fusion-primer PCR amplifies fungal ITS1, installs the flow-cell arms and a twelve-base Golay sample barcode, and creates landing sites for three locus-specific custom sequencing primers."
CAVEAT = ("The published Read 2 primer is a degenerate pool (V and R). The landing diagram uses one valid member, V=G and R=A, so its exact binding can be computed; the ordered degenerate sequence is retained below. The fungal insert length and Golay identity vary. "
          "<b>Published barcode discrepancy:</b> the EMP ITS web page, reverse-primer construct, and its 2018 protocols.io version, p. 4, print <code>NNNNNNNNNN</code> (10 bases). The linked November 2016 primer workbook, Reverse Primer Constructs sheet, prints 960 individual 12-base Golay barcodes, and Smith and Peay 2014, Molecular methods, p. 3, says both the barcode and index read are 12 bases. This model uses 12 bases because the ordering workbook gives the actual oligos; whether the 10-N web/protocol placeholder was intended as shorthand or a different construct is unresolved. "
          "<b>Published amplicon-size discrepancy:</b> the current EMP ITS page states \"Amplicon size: ~250–600 bp\" under Thermocycler conditions but \"Expected band size for ITS1f-ITS2 is ~230 bp\" in Amplification step 3; the 2018 protocols.io version gives ~230 bp in both places. The model draws the fungal insert as a variable-length placeholder and does not choose a size.")

FWD_ARM = "AATGATACGGCGACCACCGAGATCTACAC"
FWD_LINK = "GG"
FWD_LOCUS = "CTTGGTCATTTAGAGGAAGTAA"
REV_ARM = "CAAGCAGAAGACGGCATACGAGAT"
REV_LINK = "CG"
REV_LOCUS = "GCTGCGTTCTTCATCGATGC"

# November 2016 EMP workbook, sheet "Forward and Sequencing Primers" (rows 4-9).
READ1_PRIMER = "TTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC"
READ2_PRIMER_ORDERED = "CGTTCTTCATCGATGCVAGARCCAAGAGATC"
READ2_PRIMER_REPRESENTATIVE = "CGTTCTTCATCGATGCGAGAACCAAGAGATC"  # V=G, R=A
INDEX1_PRIMER = "TCTCGCATCGATGAAGAACGCAGCCG"

# The workbook primers extend beyond the PCR primer's locus-specific 3' ends.
R1_EXTENSION = READ1_PRIMER[len(FWD_LOCUS) - 1:]
R2_EXTENSION = READ2_PRIMER_REPRESENTATIVE[len(REV_LOCUS) - 4:]

FINAL_LIBRARY = Construct([
    seg("P5-side fusion arm", FWD_ARM, "p5"),
    seg("forward linker", FWD_LINK),
    seg("ITS1f site", FWD_LOCUS, "r1"),
    seg("Read 1 locus extension", R1_EXTENSION),
    seg("fungal ITS1 insert", "N" * 34, placeholder=True),
    seg("Read 2 locus extension reverse complement", revcomp(R2_EXTENSION)),
    seg("ITS2 site reverse complement", revcomp(REV_LOCUS), "r2"),
    seg("reverse linker", revcomp(REV_LINK)),
    seg("twelve-base Golay barcode", "B" * 12, "cbc", placeholder=True,
        feature=feature("sample_golay", "sample_index", "whitelist", whitelist="EMP 12-base Golay barcode set")),
    seg("P7-side fusion arm reverse complement", revcomp(REV_ARM), "p7"),
], name="EMP ITS fusion-primer library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "EMP ITS Read 1", READ1_PRIMER,
              "EMP November 2016 workbook, Forward and Sequencing Primers sheet"),
    sp.custom("Index 1 (i7)", "EMP ITS barcode primer", INDEX1_PRIMER,
              "EMP November 2016 workbook, Forward and Sequencing Primers sheet"),
    sp.custom("Read 2", "EMP ITS Read 2 (V=G, R=A member)",
              READ2_PRIMER_REPRESENTATIVE,
              "EMP November 2016 workbook, Forward and Sequencing Primers sheet",
              "One explicit member of the published degenerate primer pool."),
)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2")):
    raise ValueError("EMP ITS final library: " + "; ".join(problems))

FINAL_CAPTION = "EMP.ITSkabir paired-end, single-index fusion-primer product. The twelve-base Golay barcode is read in its own custom index read."
SEQUENCING_INTRO = "The cleaned pool is sequenced with 5–10% PhiX added for complexity (EMP ITS page). Three published custom primers sequence the fungal insert from both ends and read the twelve-base Golay barcode. There is no i5/Index 2 read. The Read 2 oligo shown for landing is one member of the ordered V/R-degenerate pool. The EMP page defers run settings to Caporaso et al. 2012; acquired cycle counts for EMP ITS are not fixed here."
SHOW_READ_CONTENT = True

def sections():
    return [
        ("One PCR installs the complete fusion-primer construct", [
            Row(chunks=[("P5—GG—ITS1f → fungal ITS1 ← ITS2—CG—Golay12—P7", "cbc", False)]),
        ], "Both locus selection and platform-arm installation occur in the same amplification."),
        ("Custom sequencing primers cross into the locus", [
            Row(chunks=[("Read 1: ITS1f suffix + 19 locus bases →", "r1", False)]),
            Row(chunks=[("← Read 2: ITS2 suffix + 15 degenerate locus bases", "r2", False)]),
            Row(chunks=[("Index 1: four locus-extension bases + ITS2/CG → Golay12", "cbc", False)]),
        ], "The longer locus-matching primers raise melting temperature in the low-complexity amplicon. The index primer approaches the barcode through the reverse-primer side."),
    ]


PAGE = "EMP ITS Illumina Amplicon Protocol web page (retrieved 2026-10-17)"
PROTOCOLS_IO = "EMP ITS Illumina Amplicon Protocol, protocols.io doi:10.17504/protocols.io.pa7dihn (Apr 14, 2018)"


def forward_primer() -> Construct:
    return Construct([seg("P5-side fusion arm", FWD_ARM, "p5"), seg("forward linker", FWD_LINK),
                      seg("ITS1f site", FWD_LOCUS, "r1")], name="EMP.ITSkabir forward primer (ITS1f)")


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("community_dna","fungal community DNA (1 µl per 25 µl reaction)",MaterialKind.MOLECULAR,
                    primary=True, source=PAGE)
    forward=graph.input("fusion_forward","P5-ITS1f fusion primer, 0.2 µM final",MaterialKind.OLIGO,
                        state=forward_primer(), source=PAGE)
    reverse=graph.input("fusion_reverse","Golay12-P7-ITS2 fusion primer (one barcode per sample), 0.2 µM final",
                        MaterialKind.OLIGO, source="EMP November 2016 workbook, Reverse Primer Constructs sheet")
    polymerase=graph.input("pcr_polymerase","Platinum Hot Start PCR Master Mix (2x), 0.8x final",
                           MaterialKind.CATALYST, source=PAGE)
    cycling=graph.input("pcr_cycling","94 °C 1 min; 35 × (94 °C 30 s, 52 °C 30 s, 68 °C 30 s); 68 °C 10 min",
                        MaterialKind.CONDITION, source=PAGE)
    amplicon=graph.product("sample_amplicon","one sample's amplicon (three replicate 25 µl PCRs pooled to 75 µl)",
                           MaterialKind.MOLECULAR,state=FINAL_LIBRARY, source=PAGE)
    graph.react("fusion_pcr","Amplify fungal ITS1 in triplicate and pool the replicates of each sample",
                StepKind.PCR,
                ((dna,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST),
                 (cycling,InputRole.CONDITION)),
                (amplicon,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The fusion-primer sequences are exact, but the variable fungal genomic template is not represented as one fixed Construct."),
                source=PAGE)
    others=graph.input("other_samples","other samples' amplicons with distinct Golay barcodes",
                       MaterialKind.MOLECULAR, source=PAGE)
    pooled=graph.product("equimolar_pool","pool of 240 ng per sample (agarose-gel checked, PicoGreen-quantified)",
                         MaterialKind.MOLECULAR,state=FINAL_LIBRARY, source=PAGE)
    graph.react("pool_samples","Combine an equal amount (240 ng) of each sample's amplicon",
                StepKind.POOL,
                ((amplicon,InputRole.TEMPLATE),(others,InputRole.TEMPLATE)),(pooled,),
                source=PAGE)
    column=graph.input("cleanup_kit","MoBio UltraClean PCR Clean-Up Kit",MaterialKind.SELECTION,
                       source=PAGE)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY, source=PAGE)
    waste=graph.product("cleanup_waste","primers, dNTPs and PCR reagents in the clean-up flow-through",
                        MaterialKind.WASTE, source=PAGE)
    graph.react("pool_cleanup","Clean the amplicon pool (optional gel extraction of half the pool not drawn)",
                StepKind.CLEANUP,
                ((pooled,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(library,waste),
                source=PAGE, changes_molecule=False)
    graph.discard(waste,"removed by the clean-up kit")
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "community_dna","fusion_forward","fusion_reverse","pcr_polymerase","pcr_cycling",
        "fusion_pcr","other_samples","pool_samples","cleanup_kit","pool_cleanup",
        "library"),SOURCE)
    return graph.require_publishable()
