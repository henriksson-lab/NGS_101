"""sci-Plex nuclear-hashing model.

Source scope (bounded pass 2026-10-17): Srivatsan et al. Science 2020 author
manuscript PMC7289078 main text and figure legends, plus the authors' twelve IDT
hash-oligo order forms.  The Science Supplementary Materials (Materials and
Methods, Tables S1-S8) are gated and were not read, so the version-specific
sci-RNA-seq RT, ligation and PCR oligos are not drawn here.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, oligo, revcomp
from protocol_graph import (Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind, Terminal,
                            AnnotationSlot)

# Author IDT plate-1 order form, well A1 (sciPlex_0), transcribed 5'->3'.
HASH_HANDLE="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
HASH_BARCODE="TTCTCGCATG"
if HASH_HANDLE!=nx.READ2_PRIMER:
    raise AssertionError("the order form's 'Nextera Read2' arm must equal nextera.READ2_PRIMER")
HASH_FEATURE=feature("sciplex_hash","sample_index","whitelist",
                     whitelist="author IDT hash oligo plates 1-12 (1,152 rows)")
CELL_FEATURE=feature("sciplex_cell","cell_barcode","combinatorial",group="cell",
                     part="version-specific sci-RNA-seq index rounds")
UMI_FEATURE=feature("sciplex_umi","umi","random")


def hash_oligo():
    return [Segment("Nextera Read 2 handle",HASH_HANDLE,"r2"),
            Segment("10-base hash",HASH_BARCODE,"cbc",feature=HASH_FEATURE),
            Segment("B as printed in IDT form","B",placeholder=True),
            Segment("poly(A)32","A"*32)]


def hash_library():
    """Hash-derived molecule: oligo-dT-primed copy of the hash, P5 side left.

    The right-hand arm is the reverse complement of the order-form hash oligo; the
    sci-RNA-seq RT primer, cell-index rounds and both outer PCR arms are role tokens
    because the inspected sci-Plex sources do not print them.
    """
    return Construct([
        Segment("P5-side sci-RNA-seq arm","[P5 · i5 · Read 1 site]",placeholder=True),
        Segment("cell index rounds","[sci-RNA-seq cell indexes]",placeholder=True,
                feature=CELL_FEATURE),
        Segment("UMI","[UMI]",placeholder=True,feature=UMI_FEATURE),
        Segment("oligo-dT copy of hash poly(A)","[oligo-dT · poly(T)]",placeholder=True),
        Segment("complement of printed B","[B complement]",placeholder=True),
        Segment("hash barcode (reverse complement)",revcomp(HASH_BARCODE),"cbc",
                feature=HASH_FEATURE),
        Segment("Nextera Read 2 handle (reverse complement)",revcomp(HASH_HANDLE),"r2"),
        Segment("P7-side sci-RNA-seq arm","[i7 · P7]",placeholder=True),
    ],name="hash-derived sci-Plex library molecule")


FINAL_LIBRARY=hash_library()
SEQ_PRIMERS=(sp.NEXTERA["R2"],)
TITLE="sci-Plex — nuclear hashing before combinatorial RNA indexing"
NOTES="01_sci-plex.html"
SOURCE='<a href="https://doi.org/10.1126/science.aax6234">Srivatsan et al., <i>Science</i> 367:45–51 (2020)</a>; hash oligos from the authors\' IDT order forms (github.com/cole-trapnell-lab/sci-plex).'
SUMMARY="Each treated sample is lysed in the presence of a well-specific polyadenylated ssDNA hash oligo and the nuclei are fixed. Pooled nuclei then undergo sci-RNA-seq (two-level in the barnyard pilot, sci-RNA-seq3 in the 188-compound screen), so hash oligos and endogenous mRNAs receive the same combinatorial cell indexes."
CAVEAT="Only the hash oligo is drawn base by base. The sci-RNA-seq RT primer, cell-index rounds, UMI and outer PCR arms differ between the two sci-RNA-seq versions used and are in the Science Supplementary Materials, which were not available to this review; they appear as bracketed role tokens."
FINAL_CAPTION="Hash-derived molecule: the right-hand arm is the reverse complement of the order-form hash oligo, so its Nextera Read 2 handle supplies the Read 2 site without tagmentation. Bracketed regions are version-specific sci-RNA-seq parts that are not drawn. Transcript-derived molecules share the left-hand arm but are not drawn."
SEQUENCING_INTRO="Read 2 is computed from the hash oligo's Nextera Read 2 handle and reads the hash barcode as printed. Read 1 and the index reads come from the version-specific sci-RNA-seq arms and carry the cell index combination and UMI; their primers and cycle counts are not drawn here."
SHOW_READ_CONTENT=False
READ_SPECS=(
    sp.read_spec("Read 1",None,"P5-side sci-RNA-seq arm (cell-index/UMI content); version-specific primer and layout not printed in the inspected sources","undetermined"),
    sp.read_spec("Index 1 (i7)",None,"P7-side sci-RNA-seq arm; primer, index content and cycles not printed in the inspected sources","undetermined"),
    sp.read_spec("Index 2 (i5)",None,"P5-side sci-RNA-seq arm; primer, index content and cycles not printed in the inspected sources","undetermined"),
    sp.read_spec("Read 2",None,"hash-derived molecules: computed Nextera Read 2 entry reads the 10-base hash as printed, then B and poly(A); transcript-derived molecules: cDNA (not drawn)","extract_then_include"),
)


def oligos():
    """Author IDT plate-1 ordering form, well A1 (sciPlex_0)."""
    yield oligo("sciPlex_0 nuclear hash",hash_oligo())


def protocol_graph():
    source=("Srivatsan2020 author manuscript PMC7289078 Results/Fig. 1 and Fig. 3 legends; "
            "hash oligo from author IDT plate-1 order form well A1")
    graph=ProtocolGraph("sci-Plex",authority=source)
    nuclei=graph.input("nuclei","treated cells, one perturbation per well",
                       MaterialKind.BIOLOGICAL,primary=True)
    hash_in=graph.input("hash_oligo","well-specific polyadenylated ssDNA hash (sciPlex_0 shown)",
                        MaterialKind.OLIGO,state=tuple(hash_oligo()))
    pfa=graph.input("paraformaldehyde","paraformaldehyde fixation",MaterialKind.CONDITION)
    hashed=graph.product("hashed_nuclei","hash-labelled fixed nuclei",MaterialKind.BIOLOGICAL)
    unbound=graph.product("unbound_hash","hash oligo not associated with nuclei",
                          MaterialKind.MOLECULAR)
    graph.react("hash_lysis_fixation",
                "Lyse cells in the presence of the well-specific hash oligo, then fix nuclei",
                StepKind.BINDING,
                ((nuclei,InputRole.TEMPLATE),(hash_in,InputRole.STRUCTURAL),
                 (pfa,InputRole.CONDITION)),(hashed,unbound),
                source="PMC7289078 Results ¶3 (nuclear lysis with 96 well-specific hash oligos, paraformaldehyde fixation) and Fig. 1B legend")
    graph.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,
        "The inspected main text does not describe washing or removal of free hash oligo; "
        "the supplementary Materials and Methods were not read."))
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    transcripts=graph.product("transcript_library","transcript-derived sci-RNA-seq molecules",
                              MaterialKind.MOLECULAR)
    graph.react("sci_rna_indexing",
                "Pool hashed nuclei and run sci-RNA-seq (two-level pilot; sci-RNA-seq3 screen)",
                StepKind.MAGIC,((hashed,InputRole.TEMPLATE),),(library,transcripts),
                reason=("The main text states only that pooled hashed nuclei undergo a two-level "
                        "sci-RNA-seq experiment (barnyard pilot) or sci-RNA-seq3 (large screen) and "
                        "that the polyadenylated hash is indexed like mRNA. RT primer, ligation and "
                        "PCR oligos, enzymes and conditions are in the unread Science supplement."),
                source="PMC7289078 Results ¶3 and Fig. 3A legend")
    graph.mark_terminal(Terminal(transcripts.id,Disposition.UNRESOLVED,
        "Sequenced in the same pool, but their version-specific structure (Tn5-derived Read 2 end) "
        "is not drawn from the inspected sci-Plex sources."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=("Read 2",),
        read_uses={"Read 2":"extract"},
        limitation=("Read 1 and index-read primers, cycle counts and instrument are not printed in "
                    "the inspected main text; only the Read 2 site carried by the hash is computed.")))
    graph.annotate(Highlight("library","Model correction (2026-10-17)",
        "An earlier model drew a TruSeq Read 2 arm after a dA junction. The author hash oligo "
        "instead carries the Nextera Read 2 handle (s7 + mosaic end), so the hash-derived "
        "molecule ends in that handle; a TruSeq Read 2 arm is not supported.",
        slot=AnnotationSlot.SEQUENCING))
    graph.declare_reviewed_scope((
        "nuclei","hash_oligo","paraformaldehyde","hash_lysis_fixation","unbound_hash",
        "sci_rna_indexing","transcript_library","library"),source)
    return graph.require_publishable()
