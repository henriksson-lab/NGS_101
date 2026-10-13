"""Original Schmitt et al. Duplex Sequencing library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import duplex_family_rows, seg, truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Row, Scene, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceCoverage, SourceItem, StepKind)

TITLE = "Duplex Sequencing — complementary-strand consensus"
NOTES = "01_duplex-sequencing.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1073/pnas.1208715109">Schmitt et al., <i>PNAS</i> (2012)</a>.'
SUMMARY = "Complementary double-stranded random tags identify both strands of each original DNA duplex. PCR descendants first form two single-strand consensuses, which are paired by their swapped end tags to form a duplex consensus."
CAVEAT = "α and β denote the two independently random 12-nt end tags on one original fragment. They are molecule identifiers, not fixed sequences; complementary strands carry the reciprocal αβ/βα relationship by construction."
TAG_ALPHA=feature("duplex_tag_alpha","umi","random",group="duplex_tag",part="alpha")
TAG_BETA=feature("duplex_tag_beta","umi","random",group="duplex_tag",part="beta")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("duplex tag α", "N" * 12, "umi", placeholder=True, feature=TAG_ALPHA),
    seg("genomic insert", "X" * 38, placeholder=True),
    seg("duplex tag β", "N" * 12, "umi", placeholder=True, feature=TAG_BETA)],
    "Duplex Sequencing library", dual_index=False)
FINAL_CAPTION = "Each end contributes a 12-nt duplex tag, giving a 24-nt αβ identity for one strand family. The complementary family is recognized by the reciprocal βα tag relationship."
SEQUENCING_INTRO = "Paired-end sequencing observes the two end tags and the intervening insert. Consensus construction is part of the molecular design: strand families are linked only when their tags are complementary and swapped."
SHOW_READ_CONTENT = True

def sections():
    adapter = [seg("asymmetric PCR arm", "X" * 12, placeholder=True),
               seg("12-nt random tag", "N" * 12, "umi", placeholder=True,
                   feature=feature("duplex_tag_adapter", "umi", "random")),
               seg("3-prime dA", "A")]
    return [
        ("Synthesize a double-stranded random-tag adapter",
         Scene.duplex(adapter, label="Duplex Tag adapter").rows(),
         "A single-stranded N12 region is copied to make a complementary duplex tag. Extended polymerase incubation completes the adapter’s 3′ dA."),
        ("Ligate tagged adapters to both ends of T-tailed DNA",
         [Row(chunks=[("adapter—N12(α)—A:T ** genomic duplex ** T:A—N12(β)—adapter", "umi", False)])],
         "The original implementation T-tailed the sheared insert and A-tailed the adapters. ** marks the sealed adapter–insert boundaries."),
        ("Build paired strand families and require duplex agreement",
         duplex_family_rows(),
         "Reads sharing αβ form one SSCS; reads sharing the reciprocal βα tags form its complementary SSCS. Only agreeing bases survive into the DCS."),
    ]


def tag_adapter_geometry():
    adapter=[seg("adapter arm","X"*12,placeholder=True),
             seg("12-nt random tag","N"*12,"umi",placeholder=True,
                 feature=feature("duplex_tag_adapter_geometry","umi","random")),
             seg("adapter 3-prime dA","A")]
    adapter_scene=Scene.duplex(adapter,label="Duplex Tag adapter")
    adapter_scene.ligation_end("top","adapter 3-prime dA",kind="cohesive",
                               polarity="3-prime",owner="adapter",text="adapter 3-prime dA")
    substrate=Scene()
    substrate.strand("T-tailed insert",[seg("genomic duplex end","X"*18,placeholder=True),
                                         seg("insert 3-prime dT","T")],label="T-tailed DNA")
    substrate.overhang("T-tailed insert","insert 3-prime dT",polarity="3-prime",
                       owner="substrate",text="insert 3-prime dT")
    return AdapterLigationGeometry.known(
        "cohesive",adapter_scene,substrate,
        statement="The original method ligates 3-prime-dA tag adapters to 3-prime-dT-tailed genomic fragments.")


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("dna","sheared genomic DNA duplexes",MaterialKind.MOLECULAR,
                    primary=True)
    tailing_enzyme=graph.input("tailing_enzyme","DNA 3-prime dT-tailing enzyme",
                               MaterialKind.CATALYST)
    t_tailed=graph.product("t_tailed","3-prime-dT-tailed genomic duplexes",
                           MaterialKind.MOLECULAR)
    graph.react("tail_insert","Add 3-prime dT to sheared genomic DNA",StepKind.TAILING,
                ((dna,InputRole.TEMPLATE),(tailing_enzyme,InputRole.CATALYST)),
                (t_tailed,))
    tag_oligo=graph.input("tag_oligo","single-stranded N12 Duplex Tag adapter precursor",
                          MaterialKind.ADAPTER)
    adapter_polymerase=graph.input("adapter_polymerase","Duplex Tag adapter-copy polymerase",
                                   MaterialKind.CATALYST)
    tag_adapter=graph.product("tag_adapter","double-stranded N12 dA-tailed Duplex Tag adapter",
                              MaterialKind.ADAPTER,state=tuple(tag_adapter_geometry().rows()))
    graph.react("copy_adapter","Copy the random tag into a complementary duplex adapter",
                StepKind.EXTENSION,
                ((tag_oligo,InputRole.TEMPLATE),(adapter_polymerase,InputRole.CATALYST)),
                (tag_adapter,))
    ligase=graph.input("ligase","Duplex Tag adapter ligase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","duplex-tagged genomic molecules",
                         MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("ligate_tags","Ligate complementary tags to both genomic-fragment ends",
                StepKind.LIGATION,
                ((t_tailed,InputRole.TEMPLATE),(tag_adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(tagged,),geometry=tag_adapter_geometry())
    forward=graph.input("pcr_forward","Duplex Sequencing PCR forward primer",
                        MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","Duplex Sequencing PCR reverse primer",
                        MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify duplex-tagged strand families",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The source establishes the PCR arms, while sample-index identity varies by library."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "dna","tailing_enzyme","tail_insert","tag_oligo","adapter_polymerase",
        "copy_adapter","ligase","ligate_tags","pcr_forward","pcr_reverse",
        "pcr_polymerase","library_pcr","library"),SOURCE)
    graph.add_source_item(SourceItem(
        "duplex_consensus","SSCS pairing and duplex-consensus computation",SOURCE))
    graph.cover_source_item(SourceCoverage(
        "duplex_consensus","excluded",
        reason="Consensus calling occurs after sequencing and does not transform a library-preparation molecule."))
    return graph.require_publishable()
