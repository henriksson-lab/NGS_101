"""Yeast TeloTag nanopore workflow (Sholes et al. 2022)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import (Construct,Row,Scene,Segment,complement_segments,feature,oligo,
                      feature_rows)
from telomere import TelomereEnd
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

SOURCE=('Sholes et al. 2022, <a href="https://doi.org/10.1101/gr.275868.121">'
        'doi:10.1101/gr.275868.121</a>.')
TITLE="Yeast TeloTag nanopore sequencing"

END=TelomereEnd(repeat="TGGTGT",duplex_copies=3,overhang_copies=2)
# Supplemental Table S5 (the Methods text points to S6, which is the plasmid table).
TELOTAG_CORE="CAGTCTACACATATTCTCTG"
TELOTAG_PRIMERS=(
    TELOTAG_CORE+"T"*35,
    TELOTAG_CORE+"T"*36+"CAC",
    TELOTAG_CORE+"T"*36+"ACA",
    TELOTAG_CORE+"T"*36+"CCA",
    TELOTAG_CORE+"T"*36+"ACC",
)
POLYA="A"*35  # representative minimum paired extent; TdT tail length varies
SAMPLE_INDEX=feature("ont_native_sample_index", "sample_index", "whitelist",
                     whitelist="Oxford Nanopore EXP-NBD104 native barcode set")

def tail_scene():
    top=[Segment("subtelomere","X"*16,placeholder=True),Segment("telomere",END.duplex,"r1"),Segment("native overhang",END.overhang,"r2"),Segment("poly(A) tail",POLYA,"r3")]
    sc=Scene(); sc.strand("G-rich strand",top,label="chromosome")
    sc.anneal("C-rich strand",complement_segments(top[:2]),to="G-rich strand",pair=("telomere'","telomere"),label="chromosome")
    tag=[Segment("TeloTag core",TELOTAG_CORE,"r3"),Segment("oligo(dT)","T"*35,"r2")]
    sc.anneal("TeloTag",tag,to="G-rich strand",pair=("oligo(dT)","poly(A) tail"),label="TeloTag",unpaired=("TeloTag core",))
    sc.footer("** T4 ligase seals the filled TeloTag–chromosome nick","TeloTag","TeloTag core")
    return sc

def nanopore_library() -> Construct:
    return Construct([
        Segment("motor-loaded ONT adapter", "[motor adapter]", "r2", placeholder=True),
        Segment("left ONT native barcode", "[ONT barcode]", "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
        Segment("TeloTag-labelled chromosome", "[tagged chromosome]", placeholder=True),
        Segment("right ONT native barcode", "[ONT barcode]", "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
    ], name="barcoded TeloTag nanopore library")

def nanopore_library_rows():
    lib=nanopore_library(); sc=Scene(); sc.strand("library",list(lib),label="nanopore library")
    sc.junction("library","motor-loaded ONT adapter","left ONT native barcode","adapter assembly")
    prefix=len("nanopore library") + 1 + len("5'- ")
    return [*sc.rows(),*feature_rows(lib,prefix_width=prefix)]

def oligos():
    return tuple(oligo(f"TeloTag mix {i}", (Segment("TeloTag primer", sequence),))
                 for i, sequence in enumerate(TELOTAG_PRIMERS, 1))

def panels():
    return [
        ("Terminal transferase adds poly(A) at the native 3′ chromosome end",
         tail_scene(),
         "The poly(A) tail preserves a terminal landmark; the 35-base paired footprint represents the first of five published TeloTag primers."),
        ("Sulfolobus DNA polymerase IV fills the end and T4 ligase seals the nick",
         [Row(chunks=[("[duplex chromosome][telomere][poly(A):oligo(dT)][TeloTag] **",None,False)])],
         "Fill-in makes a blunt tagged end; ** identifies the sealing ligation."),
        ("Attach native barcodes and the Oxford Nanopore sequencing adapter",
         nanopore_library_rows(),
         "EXP-NBD104 supplies the sample barcode used at both molecule ends; the motor-loaded sequencing adapter is assembled after barcoding."),
    ]


def protocol_graph() -> ProtocolGraph:
    """Protocol-owned graph from the final Methods (PMC8997346 JATS lines 185-203).

    Assembled locally rather than through the queued shared telomere recipe so
    that each reagent, condition and cleanup is the one this article prints.
    """
    g=ProtocolGraph(TITLE, authority=SOURCE)
    gdna=g.input("genomic_dna","2 µg high-molecular-weight yeast genomic DNA with native telomere ends",
                 MaterialKind.BIOLOGICAL,primary=True,
                 state=molecule("chromosome ends",(role("subtelomere"),role("native telomere 3-prime overhang"))))
    tdt=g.input("tdt","terminal transferase (NEB M0315), 20 U",MaterialKind.CATALYST)
    datp=g.input("datp","5 mM dATP in 1x terminal transferase buffer with 1x CoCl2",MaterialKind.CONDITION)
    tailed=g.product("polya_tailed","poly(A)-tailed DNA 3-prime ends (tail length not measured)",
                     MaterialKind.MOLECULAR,state=gdna.state)
    g.react("polya_tail","Terminal transferase adds poly(A) to 3-prime ends, 1 h 37 °C; 10 min 70 °C stop",
            StepKind.TAILING,((gdna,InputRole.TEMPLATE),(tdt,InputRole.CATALYST),(datp,InputRole.CONDITION)),
            (tailed,))
    tag=g.input("terminal_tag","TeloTag primer mix (five Supplemental Table S5 oligos)",
                MaterialKind.OLIGO,state=TELOTAG_PRIMERS)
    pol=g.input("gap_fill_polymerase","Sulfolobus DNA Polymerase IV (NEB M0327S), 4 U, with 2.5 mM dNTPs in ThermoPol buffer",
                MaterialKind.CATALYST)
    filled=g.product("gap_filled","TeloTag annealed to poly(A) and extended to the recessed chromosome 5-prime end",
                     MaterialKind.MOLECULAR,state=tailed.state)
    g.react("anneal_fill","Anneal TeloTag mix and fill in, 1 min 56 °C then 10 min 72 °C",StepKind.EXTENSION,
            ((tailed,InputRole.TEMPLATE),(tag,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(filled,))
    ligase=g.input("nick_ligase","T4 DNA ligase (NEB M0202), 400 U, 1 mM ATP present",MaterialKind.CATALYST)
    sealed=g.product("sealed","nick-sealed TeloTag-labelled chromosome end",MaterialKind.MOLECULAR,state=filled.state)
    g.react("seal_nick","Seal the filled TeloTag–chromosome nick, 20 min 12 °C",StepKind.NICK_SEALING,
            ((filled,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(sealed,))
    beads=g.input("ampure","AMPure XP beads",MaterialKind.SELECTION)
    cleaned=g.product("tag_cleaned","bead-purified tagged gDNA",MaterialKind.MOLECULAR,state=sealed.state)
    unbound=g.product("tag_unbound","excess TeloTag oligos, enzymes and nucleotides in the supernatant",MaterialKind.WASTE)
    g.react("cleanup_tagged","Purify tagged gDNA with AMPure XP beads",StepKind.CLEANUP,
            ((sealed,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(cleaned,unbound),changes_molecule=False)
    g.discard(unbound,"bead supernatant removed")
    sre=g.input("short_read_eliminator","Circulomics Short Read Eliminator XS (SS-100-121-01)",MaterialKind.SELECTION)
    large=g.product("size_selected","tagged gDNA retained above 10 kb",MaterialKind.MOLECULAR,state=cleaned.state)
    short=g.product("short_fragments","DNA below the kit's 10-kb retention threshold",MaterialKind.WASTE)
    g.react("size_select","Size-select tagged gDNA with Short Read Eliminator XS",StepKind.SELECTION,
            ((cleaned,InputRole.TEMPLATE),(sre,InputRole.SELECTION)),(large,short),changes_molecule=False)
    g.discard(short,"removed by the size-selection kit")
    barcodes=g.input("native_barcodes","Oxford Nanopore EXP-NBD104 native barcodes",MaterialKind.ADAPTER,
                     inferred=True,source=SOURCE+" (barcode sequences and kit end-preparation reagents are not printed)")
    barcoded=g.product("barcoded","native-barcoded tagged gDNA (1 µg input)",MaterialKind.MOLECULAR,
                       state=Construct([Segment("left ONT native barcode","[ONT barcode]","cbc",placeholder=True,feature=SAMPLE_INDEX),
                                        Segment("TeloTag-labelled chromosome","[tagged chromosome]",placeholder=True),
                                        Segment("right ONT native barcode","[ONT barcode]","cbc",placeholder=True,feature=SAMPLE_INDEX)],
                                       name="native-barcoded tagged gDNA"))
    g.react("native_barcode","Native barcoding per EXP-NBD104 kit workflow",StepKind.MAGIC,
            ((large,InputRole.TEMPLATE),(barcodes,InputRole.STRUCTURAL)),(barcoded,),
            reason="The article cites EXP-NBD104/SQK-LSK109 without printing the kit's end-preparation, barcode-ligation reagents, barcode sequences or end geometry; the historical kit protocol was not inspected.")
    motor=g.input("motor_adapter","SQK-LSK109 motor-loaded Oxford Nanopore sequencing adapter",MaterialKind.ADAPTER,
                  inferred=True,source=SOURCE+" (the kit- and revision-specific motor-bearing assembly is unavailable)")
    motor_ligase=g.input("motor_ligase","SQK-LSK109 adapter ligation reagents",MaterialKind.CATALYST)
    library=g.product("library","barcoded TeloTag nanopore library",MaterialKind.LIBRARY,state=nanopore_library())
    g.react("ligate_motor","Attach the motor-loaded nanopore sequencing adapter",StepKind.LIGATION,
            ((barcoded,InputRole.TEMPLATE),(motor,InputRole.STRUCTURAL),(motor_ligase,InputRole.CATALYST)),(library,),
            geometry=AdapterLigationGeometry.unavailable(
                "The cited SQK-LSK109 workflow does not disclose the motor-bearing adapter assembly or its end geometry."))
    g.sequence(library,SequencingPlan(acquired_roles=("Sequencing output",),
        limitation="MinION/Flongle R9.4.1 nanopore reads; motor entry is platform-driven and no sequencing primer binds. Only reads with barcodes on both ends pass, and telomere/TeloTag reads are then selected computationally (TideHunter, SeqKit), not by a physical separation."))
    g.declare_reviewed_scope(("genomic_dna","tdt","datp","polya_tail","terminal_tag","gap_fill_polymerase",
                              "anneal_fill","nick_ligase","seal_nick","ampure","cleanup_tagged",
                              "short_read_eliminator","size_select","native_barcodes","native_barcode","motor_adapter","motor_ligase",
                              "ligate_motor","library"),SOURCE)
    return g.require_publishable()
