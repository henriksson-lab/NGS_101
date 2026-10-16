"""Rohland et al. 2015 partial-UDG double-stranded ancient-DNA library."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from base_conversion import PartialUDGProduct
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Double-stranded ancient DNA with partial UDG"
NOTES = "01_partial-udg-ds-adna.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1098/rstb.2013.0624">Rohland et al., <i>Philosophical Transactions B</i> (2015)</a>.'
SUMMARY = "USER removes internal uracils before end repair, while terminal damage is retained as an authenticity signal; two molecular barcodes are ligated before indexed PCR."
CAVEAT = ("The article establishes the two 7-nt molecular barcodes and completed dual-index Illumina geometry, while the original supplementary oligo table remains deferred in this pass. The paired inline barcodes identify a DNA extract, not individual molecules or UMIs. "
          "Adapter bases are therefore shown as an inferred standard TruSeq-compatible endpoint.")
P5_BARCODE = "ATCGATT"  # representative source-printed barcode in Table 2
P7_BARCODE = "GACTTAT"  # paired representative source-printed barcode in Table 2


def library():
    def adapter(name, bases, tag=None, **kw):
        return seg(name, bases, tag, inferred=True, **kw)
    lib = Construct([
        adapter("P5", il.P5, "p5"),
        adapter("7-nt i5", "J"*7, "cbc", placeholder=True,
                feature=feature("sample_i5", "sample_index", "unknown")),
        adapter("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        seg("P5-side 7-nt molecular barcode", P5_BARCODE, "cbc",
            feature=feature("molecular_barcode_p5", "sample_index", "whitelist", whitelist="published molecular barcode set", group="extract_id", part="P5 side")),
        seg("partial-UDG ancient DNA insert", "X"*36, placeholder=True),
        seg("P7-side 7-nt molecular barcode", P7_BARCODE, "cbc",
            feature=feature("molecular_barcode_p7", "sample_index", "whitelist", whitelist="published molecular barcode set", group="extract_id", part="P7 side")),
        adapter("Read 2 site start", "A"),
        adapter("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        adapter("7-nt i7 reverse complement", "I"*7, "cbc", placeholder=True,
                feature=feature("sample_i7", "sample_index", "unknown")),
        adapter("P7 reverse complement", il.P7_RC, "p7"),
    ], name="partial-UDG ancient-DNA library")
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
    problems=sp.verify(lib,primers)
    if problems: raise ValueError("invalid partial-UDG library: "+"; ".join(problems))
    return lib,primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "INFERRED — completed dual-index Illumina endpoint around the two source-printed 7-nt molecular barcodes."
SEQUENCING_INTRO = "Read 1 and Read 2 each begin with a molecular barcode; i5 and i7 distinguish experiments after indexing PCR."


def partial_user_rows(product):
    damaged = [
        seg("left terminal U", "U", "w1"),
        seg("left surviving interior", "ACG"),
        seg("internal U cleavage site", "U", "w1"),
        seg("right surviving interior", "TCG"),
        seg("right terminal U", "U", "w1"),
    ]
    before = Scene()
    before.strand("damaged strand", damaged)
    before.mark("damaged strand", "internal U cleavage site", "USER cleavage")

    survivor = Scene()
    survivor.strand("left cleavage fragment", [seg("potentially retained left terminal U", product.left_terminal, "w1"), seg("left surviving interior", "ACG")])
    survivor.strand("right cleavage fragment", [seg("right surviving interior", "TCG"), seg("potentially retained right terminal U", product.right_terminal, "w1")])
    survivor.footer("Separate cleavage fragments: no covalent rejoining across the excised internal U is asserted.", "right cleavage fragment")
    return [*before.rows(), *survivor.rows()]



def sections():
    survivor=PartialUDGProduct("U","ACGTCG","U")
    lig=Scene.duplex(list(FINAL_LIBRARY),label="barcoded library")
    lig.junction("top","P5-side 7-nt molecular barcode","partial-UDG ancient DNA insert","adapter ligation")
    lig.junction("top","partial-UDG ancient DNA insert","P7-side 7-nt molecular barcode","adapter ligation")
    return [
        ("Partial USER treatment", partial_user_rows(survivor),
         "UDG and Endonuclease VIII act before T4 polymerase and kinase; terminal uracils remain inefficiently removed."),
        ("End repair and barcode-adapter ligation", AdapterLigationGeometry.unavailable("Rohland et al. specify blunt end repair and two seven-base molecular barcodes, while the original supplementary adapter oligo table remains deferred; the adapter-facing ends cannot be verified.").rows()+lig.rows(),
         "T4 polymerase/PNK blunt the surviving fragments; distinct P5- and P7-side 7-mers are ligated. ** marks both ligations."),
        ("INFERRED — fill in and indexed PCR", [Row(chunks=[("short barcoded library → fill-in → i5/i7 PCR → sequencing library",None,False)])],
         "The paper explicitly describes this order; unavailable supplementary oligo bases are not presented as source-transcribed."),
    ]

def workflow():
    s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])))
    g=AdapterLigationGeometry.unavailable("Rohland et al. specify blunt end repair and two seven-base molecular barcodes, while the original supplementary adapter oligo table remains deferred; the adapter-facing ends cannot be verified.")
    w.ligate_adapter(s[1][0],g,s[1][1],note=s[1][2]);w.react(s[2][0],s[2][1],note=s[2][2]);return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("ancient_dna","double-stranded ancient DNA",
                    MaterialKind.MOLECULAR,primary=True)
    udg=graph.input("udg","limited uracil-DNA glycosylase",MaterialKind.CATALYST)
    endoviii=graph.input("endonuclease_viii","Endonuclease VIII",MaterialKind.CATALYST)
    treated=graph.product("treated","partial-UDG-treated terminal-damage-retaining DNA",
                          MaterialKind.MOLECULAR,
                          state=tuple(partial_user_rows(PartialUDGProduct("U","ACGTCG","U"))))
    cleaved=graph.product("cleaved_internal","other internally uracil-cleaved fragments",
                          MaterialKind.MOLECULAR)
    graph.react("partial_user","Remove internal uracils while retaining terminal damage",
                StepKind.DIGESTION,
                ((dna,InputRole.TEMPLATE),(udg,InputRole.CATALYST),
                 (endoviii,InputRole.CATALYST)),(treated,cleaved))
    graph.mark_terminal(Terminal(cleaved.id, Disposition.UNRESOLVED, "Internal cleavage creates separate fragments; their library retention depends on subsequent end repair and recovery, not cleavage alone."))
    ugi=graph.input("ugi","UDG inhibitor (UGI), 30 min at 37°C",MaterialKind.CONDITION)
    inhibited=graph.product("inhibited","UDG-inhibited treated DNA",MaterialKind.MOLECULAR,state=treated.state)
    graph.react("inhibit_udg","Inhibit UDG with UGI before end repair",StepKind.QUENCHING,
                ((treated,InputRole.TEMPLATE),(ugi,InputRole.CONDITION)),(inhibited,),changes_molecule=False)
    repair_enzymes=graph.input("repair_enzymes","T4 DNA polymerase and polynucleotide kinase",
                               MaterialKind.CATALYST)
    repaired=graph.product("repaired","blunt end-repaired partial-UDG DNA",
                           MaterialKind.MOLECULAR,state=treated.state)
    graph.react("end_repair","Blunt surviving fragments and phosphorylate their ends",
                StepKind.END_REPAIR,
                ((inhibited,InputRole.TEMPLATE),(repair_enzymes,InputRole.CATALYST)),
                (repaired,))
    adapters=graph.input("barcoded_adapters","paired seven-base molecular-barcode adapters",
                         MaterialKind.ADAPTER,inferred=True,
                         source=SOURCE+" (supplementary adapter-end bases unavailable)")
    ligase=graph.input("adapter_ligase","blunt-end adapter ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated","dual molecular-barcode ancient-DNA library precursor",
                          MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (completed compatible shell inferred)")
    geometry=AdapterLigationGeometry.unavailable(
        "The paper specifies blunt end repair and two seven-base molecular barcodes, while the original supplementary adapter oligo table remains deferred.")
    graph.react("ligate_adapters","Ligate the two molecular-barcode adapters",
                StepKind.LIGATION,
                ((repaired,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=geometry)
    fill_primer=graph.input("fill_primer","adapter fill-in primer",MaterialKind.OLIGO,
                            inferred=True,source=SOURCE+" (supplementary sequence unavailable)")
    fill_polymerase=graph.input("fill_polymerase","Bst polymerase, large fragment",
                                MaterialKind.CATALYST)
    filled=graph.product("filled","filled-in short-adapter library",
                         MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,
                         source=ligated.source)
    graph.react("fill_in","Fill in the ligated short adapters",StepKind.EXTENSION,
                ((ligated,InputRole.TEMPLATE),(fill_primer,InputRole.STRUCTURAL),
                 (fill_polymerase,InputRole.CATALYST)),(filled,))
    prehyb_forward=graph.input("prehyb_forward","PreHyb PCR forward primer",
                               MaterialKind.OLIGO,inferred=True,
                               source=SOURCE+" (supplementary sequence unavailable)")
    prehyb_reverse=graph.input("prehyb_reverse","PreHyb PCR reverse primer",
                               MaterialKind.OLIGO,inferred=True,
                               source=SOURCE+" (supplementary sequence unavailable)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    prehyb=graph.product("prehyb","PreHyb-amplified short-adapter library",
                         MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True,
                         source=ligated.source)
    graph.react("prehyb_pcr","Amplify the short-adapter library before capture: 30 cycles with Pfu Turbo Cx",
                StepKind.PCR,
                ((filled,InputRole.TEMPLATE),(prehyb_forward,InputRole.STRUCTURAL),
                 (prehyb_reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (prehyb,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The accessible paper names PreHyb PCR with primer assignment still unresolved."))
    index_forward=graph.input("index_forward","i5 indexing PCR primer",MaterialKind.OLIGO,
                              inferred=True,source=SOURCE+" (selected index varies)")
    index_reverse=graph.input("index_reverse","i7 indexing PCR primer",MaterialKind.OLIGO,
                              inferred=True,source=SOURCE+" (selected index varies)")
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,source=ligated.source)
    graph.react("index_pcr","Complete i5 and i7 arms by indexing PCR",StepKind.PCR,
                ((prehyb,InputRole.TEMPLATE),(index_forward,InputRole.STRUCTURAL),
                 (index_reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Index identities vary by experiment and the supplementary adapter table remains deferred."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "ancient_dna","udg","endonuclease_viii","partial_user","repair_enzymes",
        "ugi","inhibit_udg","end_repair","barcoded_adapters","adapter_ligase","ligate_adapters",
        "fill_primer","fill_polymerase","fill_in","prehyb_forward","prehyb_reverse",
        "pcr_polymerase","prehyb_pcr","index_forward","index_reverse","index_pcr",
        "library"),SOURCE)
    return graph.require_publishable()
