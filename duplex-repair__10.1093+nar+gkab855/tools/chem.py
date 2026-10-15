"""Protocol-owned source, oligos, workflow and reaction graph."""
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from endprep import ThreePrimeOverhangAdapter, ligation_product_scene, react_dA_tail_then_ligate, repair_and_dA_tail
from protocol_model import closed_circle, molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from accuracy_protocols import _workflow


def workflow() -> Workflow:
    damaged = molecule("damaged or nicked DNA duplex", (role("damaged duplex with ragged ends"),))
    excised = molecule("lesion-excised duplex", (role("short gaps"), role("trimmed overhangs")))
    blunt = molecule("blunt repaired duplex", (role("restricted fill-in"), role("remaining nicks")))
    sealed = molecule("nick-sealed duplex", (role("original duplex interior"),))
    ata = molecule("dA-tailed duplex", (role("original duplex interior"), role("3′ dA", length_bp=1)))
    adapted_construct = Construct((
        role("custom duplex UMI adapter A", tag="umi",
             molecular_feature=feature("duplex_umi_a", "umi", "random")),
        role("repaired insert"),
        role("custom duplex UMI adapter B", tag="umi",
             molecular_feature=feature("duplex_umi_b", "umi", "random"))),
        name="duplex-UMI library")
    adapted = ligation_product_scene(
        adapted_construct, "repaired insert", adapter_group="duplex UMI adapter A",
        right_adapter_group="duplex UMI adapter B")
    return _workflow(damaged, (
        ("Excise damaged bases and remove single-stranded overhangs", excised,
         "EndoIV, Fpg, UDG, T4 PDG and EndoVIII act with ExoVII."),
        ("Blunt ends and restrict fill-in", blunt,
         "T4 PNK and T4 DNA polymerase fill short gaps and residual short 5′ overhangs."),
        ("Seal nicks with HiFi Taq DNA ligase", sealed,
         "Nick sealing precedes dA-tailing so polymerase does not start from open nicks."),
        ("Restricted dA-tail with dATP as the only nucleotide", ata,
         "Klenow fragment exo− and Taq add the ligation-compatible 3′ dA."),
        ("Ligate custom dual-index duplex-UMI adapters and PCR-amplify", adapted,
         "Complementary UMI families recover both original strands."),
    ))


FINAL_LIBRARY = Construct([
    role("custom left duplex-UMI adapter", tag="umi", molecular_feature=feature("duplex_umi_a", "umi", "random")),
    role("repaired insert"),
    role("custom right duplex-UMI adapter", tag="umi", molecular_feature=feature("duplex_umi_b", "umi", "random")),
], name="Duplex-Repair sequencing library")
SEQ_PRIMERS = ()
TITLE = 'Duplex-Repair — damage-aware end repair for duplex sequencing'
NOTES = '01_duplex_repair.html'
SOURCE = '<a href="https://doi.org/10.1093/nar/gkab855">Xiong et al., <i>Nucleic Acids Research</i> (2022)</a>.'
SUMMARY = 'Lesion excision, controlled blunting, nick sealing and dATP-only tailing are separated to minimize copying of one-strand damage into an apparent duplex mutation.'
CAVEAT = 'The paper identifies custom dual-index duplex-UMI adapters but does not publish their sequences; adapter geometry and primer landing positions therefore remain unspecified.'
FINAL_CAPTION = 'Published roles of the custom duplex-UMI library; exact adapter bases are unavailable.'
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The authors sequenced custom duplex-UMI libraries on HiSeq 2500 or HiSeq X in paired-end mode.'
SEQUENCING_UNAVAILABLE = 'The custom dual-index duplex-UMI adapter and primer oligos are not given in the accessible paper or supplement.'
def oligos():
    return ()

def protocol_graph():
    final = FINAL_LIBRARY
    source = "Xiong et al., Nucleic Acids Research (2022); doi:10.1093/nar/gkab855"
    graph = ProtocolGraph("Duplex-Repair", authority=source)
    damaged = graph.input(
        "damaged_dna", "damaged or nicked DNA duplex",
        MaterialKind.BIOLOGICAL,
        state=molecule("damaged or nicked DNA duplex", (
            role("damaged duplex with ragged ends"),)), primary=True)
    lesion_enzymes = graph.input(
        "lesion_enzymes", "EndoIV/Fpg/UDG/T4-PDG/EndoVIII/ExoVII repair cocktail",
        MaterialKind.CATALYST)
    excised = graph.product(
        "excised_dna", "lesion-excised duplex with trimmed overhangs",
        MaterialKind.MOLECULAR,
        state=molecule("lesion-excised duplex", (role("short gaps"), role("trimmed overhangs"))))
    graph.react(
        "lesion_excision", "Excise damaged bases and remove single-stranded overhangs",
        StepKind.CHEMICAL_CONVERSION,
        ((damaged, InputRole.TEMPLATE), (lesion_enzymes, InputRole.CATALYST)),
        (excised,))
    blunting = graph.input("blunting_enzymes", "T4 PNK and T4 DNA polymerase",
                           MaterialKind.CATALYST)
    blunt = graph.product(
        "blunt_dna", "blunt repaired duplex with restricted fill-in",
        MaterialKind.MOLECULAR,
        state=molecule("blunt repaired duplex", (
            role("restricted fill-in"), role("remaining nicks"))))
    graph.react("blunt_ends", "Blunt ends and restrict fill-in", StepKind.END_REPAIR,
                ((excised, InputRole.TEMPLATE), (blunting, InputRole.CATALYST)),
                (blunt,))
    nick_ligase = graph.input("nick_ligase", "HiFi Taq DNA ligase",
                              MaterialKind.CATALYST)
    sealed = graph.product(
        "sealed_dna", "nick-sealed repaired duplex", MaterialKind.MOLECULAR,
        state=molecule("nick-sealed duplex", (role("original duplex interior"),)))
    graph.react("seal_nicks", "Seal nicks before tailing", StepKind.NICK_SEALING,
                ((blunt, InputRole.TEMPLATE), (nick_ligase, InputRole.CATALYST)),
                (sealed,))
    tailing = graph.input(
        "restricted_tailing_enzyme", "Klenow fragment exo-minus and Taq polymerase",
        MaterialKind.CATALYST)
    datp = graph.input("datp_only", "dATP as the only tailing nucleotide",
                       MaterialKind.MOLECULAR)
    tailed = graph.product(
        "tailed_dna", "restricted dA-tailed repaired duplex",
        MaterialKind.MOLECULAR,
        state=molecule("dA-tailed duplex", (
            role("original duplex interior"), role("3-prime dA", length_bp=1))))
    graph.react("restricted_a_tail", "Restricted dA-tail with dATP only",
                StepKind.TAILING,
                ((sealed, InputRole.TEMPLATE), (datp, InputRole.STRUCTURAL),
                 (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("duplex_umi_adapters", "custom dual-index duplex-UMI adapters",
                           MaterialKind.ADAPTER)
    adapter_ligase = graph.input("adapter_ligase", "DNA adapter ligase",
                                 MaterialKind.CATALYST)
    adapted = graph.product("adapted_dna", "duplex-UMI adapter-ligated DNA",
                            MaterialKind.MOLECULAR, state=final)
    graph.react("adapter_ligation", "Ligate custom duplex-UMI adapters",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper identifies custom duplex-UMI adapters but does not publish their ligating-end sequences."))
    forward = graph.input("pcr_forward", "custom duplex-library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "custom duplex-library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "Duplex-Repair sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "PCR-amplify the duplex-UMI library", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The accessible source does not disclose the custom duplex-library primer sequences."))
    graph.sequence(library, SequencingPlan(
        limitation="Paired-end HiSeq 2500 or HiSeq X sequencing is reported, but custom adapter and sequencing-primer bases are unavailable."))
    graph.declare_reviewed_scope((
        "damaged_dna", "lesion_enzymes", "lesion_excision", "blunting_enzymes",
        "blunt_ends", "nick_ligase", "seal_nicks", "restricted_tailing_enzyme",
        "datp_only", "restricted_a_tail", "duplex_umi_adapters", "adapter_ligase",
        "adapter_ligation", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()
