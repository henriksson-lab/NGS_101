"""Sequence model for ASTAR-seq.

Provenance of every base drawn here:

* the three C1 oligos (`C1_HANDLE`, `C1_DT_SPACER`, `C1_DT`, `C1_TSO`, `C1_PCR`)
  are transcribed from the authors' own oligo table, printed identically in the
  defining preprint (bioRxiv 829960, Supplementary Table 5, `media-5.xlsx`) and in
  the journal version (Genome Research 30:1027, Supplemental Table 5);
* every other concrete base is a canonical constant imported from `illumina` or
  `nextera` (P5, P7, s5/s7, mosaic end), never retyped here;
* the i5 and i7 windows are placeholders.

Nothing in this module comes from the scg_lib_structs drawn page: that page is
inspiration only, no copy of it was fetched for this protocol, and agreement with
it is not evidence.  The RNA library uses the same Nextera back end as the ATAC
library after the separated cDNA is tagmented.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
from illumina import P5, P7
import nextera as nx
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "ASTAR-seq — chromatin accessibility and transcriptome from one cell"
SOURCE = 'Defining source: <a href="https://doi.org/10.1101/829960">Xing et al. (2019)</a>.'


C1_HANDLE = "GGCGACAACACCGATTGATCA"
C1_DT_SPACER = "CG"
C1_DT_LEN = 31
C1_DT = C1_HANDLE + C1_DT_SPACER + "T" * C1_DT_LEN
C1_TSO = C1_HANDLE + "GGG"       # ordered as an all-RNA oligo
C1_PCR = C1_HANDLE

# Modifications as printed in Supplementary Table 5.
C1_DT_MOD = "/5BiosG/"
C1_PCR_MOD = "/5BiosG/"
C1_TSO_MOD = "all bases are ribonucleotides"

ATAC_QPCR_S5 = nx.ADAPTOR_S5
ATAC_QPCR_S7 = nx.ADAPTOR_S7
INDEX_LEN = 8


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def c1_dt_oligo() -> list[Segment]:
    """C1-P2-T31, exactly as Supplementary Table 5 prints it (5' /5BiosG/)."""
    return [_seg("C1 handle", C1_HANDLE, "tso"), _seg("CG", C1_DT_SPACER),
            _seg("poly(T)", "T" * C1_DT_LEN)]


def c1_tso_oligo() -> list[Segment]:
    """C1-P2-RNA-TSO, as printed; every base is ordered as a ribonucleotide."""
    return [_seg("C1 handle", C1_HANDLE, "tso"), _seg("GGG", "GGG", "tso")]


def c1_pcr_oligo() -> list[Segment]:
    """C1-P2-PCR-2, as printed (5' /5BiosG/): the bare handle."""
    return [_seg("C1 handle", C1_PCR, "tso")]


def c1_cdna() -> Construct:
    """Representative full-length cDNA after template switching and C1 PCR.

    The transcript body is a placeholder; both terminal handles derive from the
    shared sequence on the oligo-dT and TSO.
    """
    return Construct([
        _seg("C1 handle, oligo-dT end", C1_HANDLE, "tso"),
        _seg("CG", C1_DT_SPACER),
        _seg("poly(T)", "T" * C1_DT_LEN),
        _seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("untemplated CCC", "CCC", "tso"),
        _seg("C1 handle, TSO end", revcomp(C1_HANDLE), "tso"),
    ], name="ASTAR full-length cDNA")


def nextera_library(insert_name: str = "insert") -> Construct:
    """Dual-indexed Nextera library used for either ASTAR modality."""
    lib = Construct([
        _seg("P5", P5, "p5"),
        _seg("i5", "N" * INDEX_LEN, "cbc", placeholder=True,
             feature=feature("sample_i5", "sample_index", "unknown")),
        _seg("s5", nx.S5, "s5"),
        _seg("ME", nx.ME, "me"),
        _seg(insert_name, "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("ME'", nx.ME_RC, "me"),
        _seg("s7'", nx.S7_RC, "s7"),
        _seg("i7'", "N" * INDEX_LEN, "cbc", placeholder=True,
             feature=feature("sample_i7", "sample_index", "unknown")),
        _seg("P7'", revcomp(P7), "p7"),
    ], name=f"ASTAR {insert_name} library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError(f"invalid ASTAR {insert_name} library: " + "; ".join(problems))
    return lib


def atac_library() -> Construct:
    return nextera_library("accessible genomic DNA")


def rna_library() -> Construct:
    return nextera_library("cDNA")


SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cell=graph.input("single_cell","one cell in a Fluidigm C1 chamber",
                     MaterialKind.BIOLOGICAL,primary=True)
    transposome=graph.input("atac_transposome","s5/s7-loaded ATAC transposomes",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged_dna=graph.product("tagged_dna","tagmented accessible genomic DNA",
                             MaterialKind.MOLECULAR)
    cellular_rna=graph.product("cellular_rna","RNA retained in the same chamber",
                               MaterialKind.MOLECULAR)
    tag_waste=graph.product("atac_waste","non-library tagmentation products",
                            MaterialKind.WASTE)
    graph.react("atac_tagment","Tagment accessible chromatin in the intact single-cell reaction",
                StepKind.TAGMENTATION,
                ((cell,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_dna,cellular_rna,tag_waste))
    graph.discard(tag_waste,"non-amplifiable strand-transfer products")
    dt=graph.input("biotin_dt","biotinylated C1 oligo-dT primer",MaterialKind.OLIGO,
                   state=tuple(c1_dt_oligo()))
    rt=graph.input("reverse_transcriptase","SuperScript IV reverse transcriptase",
                   MaterialKind.CATALYST)
    first=graph.product("first_strand","biotinylated first-strand cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Reverse-transcribe polyadenylated RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cellular_rna,InputRole.TEMPLATE),(dt,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    tso=graph.input("tso","all-RNA C1 template-switch oligo",MaterialKind.OLIGO,
                    state=tuple(c1_tso_oligo()))
    switched=graph.product("switched_cdna","template-switched full-length cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("template_switch","Install the second C1 handle",StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    c1_primer=graph.input("c1_primer","biotinylated C1 PCR primer",MaterialKind.OLIGO,
                          state=tuple(c1_pcr_oligo()))
    pre_pol=graph.input("preamp_polymerase","on-chip cDNA PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","biotinylated full-length cDNA",
                            MaterialKind.MOLECULAR,state=c1_cdna())
    graph.react("cdna_preamp","Amplify cDNA for five on-chip cycles",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(c1_primer,InputRole.STRUCTURAL),
                 (pre_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.symmetric(
                    c1_cdna(),PCRPrimer("C1-P2-PCR-2",C1_PCR),
                    statement="The one published C1 handle primer binds the shared handle at both cDNA ends."))
    mixed=graph.product("mixed_modalities","mixed tagmented DNA and biotinylated cDNA",
                        MaterialKind.MOLECULAR)
    graph.react("combine_modalities","Retain both products in the same chamber output",
                StepKind.POOL,
                ((tagged_dna,InputRole.TEMPLATE),(amplified,InputRole.TEMPLATE)),(mixed,),
                changes_molecule=False)
    beads=graph.input("streptavidin","streptavidin beads",MaterialKind.SELECTION)
    bead_cdna=graph.product("bead_cdna","streptavidin-bound biotinylated cDNA",
                            MaterialKind.MOLECULAR)
    atac_fraction=graph.product("atac_fraction","unbound tagmented genomic DNA",
                                MaterialKind.MOLECULAR)
    graph.react("separate","Separate RNA and ATAC products with streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((mixed,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (bead_cdna,atac_fraction))
    gap_pol=graph.input("gap_fill_polymerase","ATAC gap-fill polymerase",MaterialKind.CATALYST)
    filled=graph.product("filled_atac","gap-filled ATAC DNA",MaterialKind.MOLECULAR)
    graph.react("atac_gap_fill","Fill the nine-base transposition gaps",StepKind.EXTENSION,
                ((atac_fraction,InputRole.TEMPLATE),(gap_pol,InputRole.CATALYST)),(filled,))
    atac_i5=graph.input("atac_i5","per-cell v2_Ad1.N i5 primer",MaterialKind.OLIGO)
    atac_i7=graph.input("atac_i7","per-cell v2_Ad2.N i7 primer",MaterialKind.OLIGO)
    index_pol=graph.input("index_polymerase","Nextera index PCR polymerase",MaterialKind.CATALYST)
    atac_final=atac_library()
    atac_lib=graph.product("atac_library",atac_final.name,MaterialKind.LIBRARY,state=atac_final)
    graph.react("atac_index_pcr","Index the ATAC fraction",StepKind.PCR,
                ((filled,InputRole.TEMPLATE),(atac_i5,InputRole.STRUCTURAL),
                 (atac_i7,InputRole.STRUCTURAL),(index_pol,InputRole.CATALYST)),(atac_lib,),
                geometry=PCRGeometry.unavailable(
                    atac_final,reason="The heterogeneous gap-filled ATAC precursor is not reconstructed base-for-base."))
    graph.sequence(atac_lib,SequencingPlan(library=atac_final,primers=SEQ_PRIMERS,
        read_lengths={"Read 1":50,"Index 1 (i7)":8,"Index 2 (i5)":8,"Read 2":50}))
    heat=graph.input("heat_release","heat strand-release condition",MaterialKind.CONDITION)
    released=graph.product("released_cdna","released non-biotinylated cDNA strand",
                           MaterialKind.MOLECULAR,state=c1_cdna())
    bound=graph.product("bound_cdna_waste","bead-bound biotinylated cDNA strands",
                        MaterialKind.WASTE)
    graph.react("release_cdna","Heat-release the non-biotinylated cDNA strand",
                StepKind.STRAND_SEPARATION,
                ((bead_cdna,InputRole.TEMPLATE),(heat,InputRole.CONDITION)),(released,bound))
    graph.discard(bound,"biotinylated complementary strands remain bead-bound")
    rna_transposome=graph.input("rna_transposome","Nextera XT transposome",MaterialKind.ADAPTER)
    rna_tagged=graph.product("rna_tagged","tagmented full-length cDNA",MaterialKind.MOLECULAR)
    rna_waste=graph.product("rna_tagment_waste","non-amplifiable cDNA tagmentation products",
                            MaterialKind.WASTE)
    graph.react("rna_tagment","Tagment the separated cDNA",StepKind.TAGMENTATION,
                ((released,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(rna_tagged,rna_waste))
    graph.discard(rna_waste,"homotypic and otherwise non-amplifiable Tn5 products")
    rna_i5=graph.input("rna_i5","RNA-library Nextera i5 primer",MaterialKind.OLIGO)
    rna_i7=graph.input("rna_i7","RNA-library Nextera i7 primer",MaterialKind.OLIGO)
    rna_final=rna_library()
    rna_lib=graph.product("rna_library",rna_final.name,MaterialKind.LIBRARY,state=rna_final)
    graph.react("rna_index_pcr","Complete and index the RNA library",StepKind.PCR,
                ((rna_tagged,InputRole.TEMPLATE),(rna_i5,InputRole.STRUCTURAL),
                 (rna_i7,InputRole.STRUCTURAL),(index_pol,InputRole.CATALYST)),(rna_lib,),
                geometry=PCRGeometry.unavailable(
                    rna_final,reason="The heterogeneous Tn5 precursor is not reconstructed base-for-base."))
    graph.sequence(rna_lib,SequencingPlan(library=rna_final,primers=SEQ_PRIMERS,
        read_lengths={"Read 1":101,"Index 1 (i7)":8,"Index 2 (i5)":8,"Read 2":101}))
    graph.declare_reviewed_scope(("single_cell","atac_transposome","tn5","atac_tagment",
      "biotin_dt","reverse_transcriptase","reverse_transcribe","tso","template_switch",
      "c1_primer","preamp_polymerase","cdna_preamp","combine_modalities","streptavidin",
      "separate","gap_fill_polymerase","atac_gap_fill","atac_i5","atac_i7",
      "index_polymerase","atac_index_pcr","atac_library","heat_release","release_cdna",
      "rna_transposome","rna_tagment","rna_i5","rna_i7","rna_index_pcr","rna_library"),SOURCE)
    return graph.require_publishable()
