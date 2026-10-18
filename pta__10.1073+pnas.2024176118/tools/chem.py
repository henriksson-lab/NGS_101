"""Source-owned PTA populations and native, unresolved KAPA/IDT library ends."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment
from pcr import PCRGeometry
from protocol_graph import (InputRole as I, MaterialKind as M,
                            ProtocolGraph, SequencingPlan, StepKind as K)

TITLE = "PTA — primary template-directed whole-genome amplification"
NOTES = "01_pta.html"
SOURCE = ('<a href="https://doi.org/10.1073/pnas.2024176118">Gonzalez-Pena et al., '
          '<i>PNAS</i> (2021)</a>, Results, Methods and Supplementary Figure S1.')
AUTHORITY = "https://doi.org/10.1073/pnas.2024176118; PMC8214697 JATS Results 83–85, Methods 172–177; published SI Fig S1"
SUMMARY = ("Random-primed phi29 synthesis incorporates exonuclease-resistant alpha-thio-ddNTPs. "
           "The entire purified PTA reaction enters KAPA HyperPlus without fragmentation, "
           "followed by IDT unique-dual-index adapter ligation, 15 PCR cycles and gel recovery.")
CAVEAT = ("Native random-primer chemistry, KAPA intermediate ends and IDT UDI strands remain "
          "unassigned in this bounded source scope. Terminated originals retain their no-3′-OH chemistry; "
          "ordinary PCR copies do not inherit a template phosphorothioate modification. "
          "Component recovery and kit processing of blocked ends are unresolved. No canonical adapter "
          "bases, index widths or read mask are substituted.")
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "NovaSeq 6000 sequencing of the recovered 400–700-bp library; native primer sites and run/index mask remain unassigned."
READ_SPECS = (sp.read_spec("Sequencing output", None, "Native PTA genomic library; exact read/index acquisition unassigned", "include"),)


def role(name):
    return Segment(name, "[" + name + "]", placeholder=True)


def amplification_product():
    """One product class; not a claim that every recovered molecule is terminated."""
    return Construct([
        role("random-primer-derived end; native modification pattern unassigned"),
        role("primary-template-derived DNA"),
        Segment("irreversible terminator", "[alpha-thio-ddN]", placeholder=True,
                note="alpha-thio-ddN: no 3′-OH; exonuclease-resistant phosphorothioate bond"),
    ], name="terminated PTA daughter class")


def populations(stage):
    """Independent polymers; role tokens deliberately supply no invented end bases."""
    sc = Scene()
    if stage == "genome":
        sc.strand("genome", [role("denatured genomic DNA")], label="genomic template population")
    elif stage in {"amplified", "purified"}:
        sc.strand("original", [role("primary genomic templates; recovery unassigned")], label="independent original DNA")
        sc.strand("terminated", list(amplification_product()), label="terminated daughter class", mod3="no 3′-OH")
        sc.strand("other", [role("other PTA DNA and residual primer classes; identities and recovery unassigned")], label="independent possible carried populations")
    elif stage == "prepared":
        sc.strand("prepared", [role("KAPA-prepared PTA DNA; strand/end chemistry and blocked-end processing unassigned")], label="native kit substrate population")
        sc.strand("carried", [role("possible carried genomic DNA, terminated products and primers; class allocation unassigned")], label="independent possible carried originals")
    elif stage == "ligated":
        sc.strand("ligated", [role("native UDI-ligated PTA DNA; adapter strands, sealed junctions and index widths unassigned")], label="native ligation-product population")
        sc.strand("carried", [role("possible unligated/terminated originals and free adapters; physical recovery unassigned")], label="independent possible carried populations")
    elif stage in {"pcr", "selected"}:
        sc.strand("copies", [role("ordinary DNA PCR daughters; native primer-derived ends and UDI member unassigned")], label="independent PCR daughter population")
        sc.strand("carried", [role("possible carried original DNA, alpha-thio-ddN products, adapters and primers; recovery unassigned")], label="independent possible original populations")
        sc.note("copies", "DNA copying does not recreate a template phosphorothioate bond. Gel selection does not establish component-specific removal within the retained window.")
    else:
        raise ValueError(stage)
    return sc


def protocol_graph():
    g = ProtocolGraph("PTA", authority=AUTHORITY)
    genome = g.input("genome", "denatured single-cell genomic DNA after source lysis/neutralization", M.BIOLOGICAL, state=populations("genome"), primary=True)
    primer = g.input("random_primer", "exonuclease-resistant random primer; native bases/length/modification pattern unassigned", M.OLIGO)
    phi29 = g.input("phi29", "REPLI-g amplification mix / phi29", M.CATALYST)
    terminators = g.input("terminators", "equal alpha-thio-ddATP/ddCTP/ddGTP/ddTTP mixture", M.MOLECULAR)
    amplified = g.product("pta_product", "PTA reaction: genomic originals and daughter populations", M.MOLECULAR, state=populations("amplified"))
    g.react("pta_amplification", "Random-prime and displace DNA at 30 °C for 8 h; heat-stop 65 °C for 3 min", K.EXTENSION,
            ((genome,I.TEMPLATE),(primer,I.STRUCTURAL),(phi29,I.CATALYST),(terminators,I.STRUCTURAL)),(amplified,))
    beads = g.input("beads", "AMPure XP, 2:1 beads:reaction", M.SELECTION)
    purified = g.product("purified", "purified entire PTA reaction; class recovery unassigned", M.MOLECULAR, state=populations("purified"))
    bead_waste = g.product("bead_waste", "physically removed bead-cleanup fraction; composition unassigned", M.WASTE)
    g.react("cleanup", "Purify PTA reaction with AMPure XP", K.CLEANUP, ((amplified,I.TEMPLATE),(beads,I.SELECTION)),(purified,bead_waste),changes_molecule=False)
    g.discard(bead_waste,"Source bead-cleanup removal; no specific terminated-product depletion is inferred.")
    kit = g.input("kapa", "KAPA HyperPlus preparation without fragmentation; historical internal steps unassigned", M.CONDITION)
    prepared = g.product("prepared", "native KAPA-prepared population", M.MOLECULAR,state=populations("prepared"))
    g.react("native_preparation", "Native KAPA preparation without fragmentation; end chemistry unresolved", K.MAGIC,
            ((purified,I.TEMPLATE),(kit,I.CONDITION)),(prepared,),
            reason="Methods177 names KAPA HyperPlus without fragmentation. Native kit intermediates and conversion/recovery of alpha-thio-ddN-blocked ends are not established in the inspected article/SI; an extendable dA-tailed duplex is not assumed.")
    adapters = g.input("udi_adapters", "IDT UDI adapters, 2.5 µM; native strands/member unassigned", M.ADAPTER)
    ligase = g.input("ligase", "native KAPA ligation mix", M.CATALYST)
    ligated = g.product("ligated", "native UDI-ligated reaction", M.MOLECULAR,state=populations("ligated"))
    g.react("adapter_ligation", "Ligate source IDT unique-dual-index adapters", K.LIGATION,
            ((prepared,I.TEMPLATE),(adapters,I.STRUCTURAL),(ligase,I.CATALYST)),(ligated,),
            geometry=AdapterLigationGeometry.unavailable("Methods177 identifies IDT UDI ligation but native physical strands, ends, junctions and blocked-end substrate allocation remain unassigned in this bounded reading."))
    forward = g.input("pcr_forward", "native library PCR forward primer; ordered strand unassigned", M.OLIGO)
    reverse = g.input("pcr_reverse", "native library PCR reverse primer; ordered strand unassigned", M.OLIGO)
    polymerase = g.input("pcr_polymerase", "native library PCR mix", M.CATALYST)
    amplified_library = g.product("library", "15-cycle library PCR population", M.MOLECULAR,state=populations("pcr"))
    g.react("library_pcr", "Amplify native ligation reaction for 15 cycles", K.PCR,
            ((ligated,I.TEMPLATE),(forward,I.STRUCTURAL),(reverse,I.STRUCTURAL),(polymerase,I.CATALYST)),(amplified_library,),
            geometry=PCRGeometry.unavailable(ligated.state,reason="Actual incoming native ligation population; complete primers and first-cycle junctions unassigned. A completed canonical amplicon is not used as PCR substrate."))
    gel = g.input("gel", "1% agarose E-Gel / ZymoClean Gel DNA Recovery", M.SELECTION)
    selected = g.product("selected", "400–700-bp recovered library fraction", M.LIBRARY,state=populations("selected"))
    rejected = g.product("gel_waste", "physically excluded gel fraction", M.WASTE)
    g.react("size_select", "Excise and recover 400–700-bp library fraction", K.SELECTION,
            ((amplified_library,I.TEMPLATE),(gel,I.SELECTION)),(selected,rejected),changes_molecule=False)
    g.discard(rejected,"Physically outside recovered gel slice; no computational read filtering is represented as molecular disposal.")
    g.sequence(selected, SequencingPlan(acquisition="unspecified", limitation=SEQUENCING_ENDING))
    g.declare_reviewed_scope(("genome","random_primer","phi29","terminators","pta_amplification",
                             "beads","cleanup","kapa","native_preparation","udi_adapters","ligase",
                             "adapter_ligation","pcr_forward","pcr_reverse","pcr_polymerase",
                             "library_pcr","gel","size_select","selected"),AUTHORITY)
    return g.require_publishable()
