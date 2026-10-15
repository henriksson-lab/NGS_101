"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('inDrops-2', '10.1093/nar/gkae1312', 'Juzėnas et al., Nucleic Acids Research (2025)', 'Open-source droplet capture supports either IVT linear amplification or template-switching PCR amplification of bead-barcoded 3′ cDNA.', 'The barcode read contains cell barcode and UMI; the insert read reports the transcript 3′ end.', '01_indrops_2.html')
PAYLOAD = (CELL, UMI, role("3-prime transcript insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ""
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted."

def workflow():
    return chain(molecule("cell and barcoded hydrogel bead in a droplet", (role("poly(A) RNA"), CELL, UMI)), (
        ("Lyse the cell and reverse-transcribe RNA on the bead oligo",
         molecule("barcoded first-strand cDNA", (CELL, UMI, role("3-prime transcript cDNA"))),
         "Cell barcode and UMI are copied into every captured molecule."),
        ("Choose linear IVT amplification or template-switching PCR amplification",
         molecule("amplified inDrops-2 cDNA", (CELL, UMI, role("transcript cDNA"), role("IVT or TSO handle"))),
         "The defining protocol validates two distinct amplification branches without treating them as identical chemistry."),
        ("Fragment or tagment the amplified product and complete Illumina adapters",
         molecule("inDrops-2 library", (CELL, UMI, role("3-prime transcript insert"), SAMPLE)),
         "Both branches retain the same barcode and UMI readout."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    droplet = graph.input(
        "droplet", "cell and barcoded hydrogel bead in a droplet",
        MaterialKind.BIOLOGICAL,
        state=molecule("cell and barcoded hydrogel bead", (
            role("poly(A) RNA"), CELL, UMI)), primary=True)
    bead_primer = graph.input("bead_primer", "barcoded bead oligo(dT) RT primer",
                              MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase",
                     MaterialKind.CATALYST)
    first = graph.product(
        "barcoded_cdna", "barcoded first-strand cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("barcoded first-strand cDNA", (
            CELL, UMI, role("3-prime transcript cDNA"), role("amplification handle"))))
    graph.react("reverse_transcription", "Lyse the cell and reverse-transcribe on the bead oligo",
                StepKind.REVERSE_TRANSCRIPTION,
                ((droplet, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    ivt_input = graph.product("ivt_input", "cDNA assigned to linear IVT amplification",
                              MaterialKind.MOLECULAR, state=first.state)
    tso_input = graph.product("tso_input", "cDNA assigned to template-switch PCR amplification",
                              MaterialKind.MOLECULAR, state=first.state)
    graph.react("choose_branch", "Retain the two validated amplification alternatives",
                StepKind.SPLIT, ((first, InputRole.TEMPLATE),), (ivt_input, tso_input),
                changes_molecule=False)
    reviewed = ["droplet", "bead_primer", "reverse_transcriptase",
                "reverse_transcription", "choose_branch"]

    t7 = graph.input("ivt_t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    amplified_rna = graph.product(
        "ivt_rna", "linearly amplified inDrops-2 RNA",
        MaterialKind.MOLECULAR,
        state=molecule("IVT-amplified RNA", (CELL, UMI, role("transcript RNA copy"))))
    graph.react("ivt", "Linearly amplify the T7-bearing cDNA by IVT",
                StepKind.TRANSCRIPTION,
                ((ivt_input, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (amplified_rna,))
    ivt_rt_primer = graph.input("ivt_rt_primer", "IVT-product reverse-transcription primer",
                                MaterialKind.OLIGO)
    ivt_rt = graph.input("ivt_reverse_transcriptase", "IVT-product reverse transcriptase",
                         MaterialKind.CATALYST)
    ivt_cdna = graph.product("ivt_cdna", "IVT-amplified inDrops-2 cDNA",
                             MaterialKind.MOLECULAR, state=first.state)
    graph.react("ivt_back_rt", "Convert amplified RNA back to cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((amplified_rna, InputRole.TEMPLATE), (ivt_rt_primer, InputRole.STRUCTURAL),
                 (ivt_rt, InputRole.CATALYST)), (ivt_cdna,))
    reviewed.extend(("ivt_t7_polymerase", "ivt", "ivt_rt_primer",
                     "ivt_reverse_transcriptase", "ivt_back_rt"))

    tso = graph.input("tso", "inDrops-2 template-switch oligo", MaterialKind.OLIGO)
    tso_rt = graph.input("tso_reverse_transcriptase", "template-switching reverse transcriptase",
                         MaterialKind.CATALYST)
    switched = graph.product("tso_cdna", "template-switched inDrops-2 cDNA",
                             MaterialKind.MOLECULAR, state=first.state)
    graph.react("template_switch", "Install the template-switch PCR handle",
                StepKind.TEMPLATE_SWITCHING,
                ((tso_input, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (tso_rt, InputRole.CATALYST)), (switched,))
    amp_forward = graph.input("tso_amp_forward", "template-switch amplification forward primer",
                              MaterialKind.OLIGO)
    amp_reverse = graph.input("tso_amp_reverse", "template-switch amplification reverse primer",
                              MaterialKind.OLIGO)
    amp_polymerase = graph.input("tso_amp_polymerase", "cDNA amplification polymerase",
                                 MaterialKind.CATALYST)
    tso_cdna = graph.product("tso_amplified_cdna", "PCR-amplified inDrops-2 cDNA",
                             MaterialKind.MOLECULAR, state=switched.state)
    graph.react("tso_amplification", "PCR-amplify the template-switched branch",
                StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (amp_forward, InputRole.STRUCTURAL),
                 (amp_reverse, InputRole.STRUCTURAL), (amp_polymerase, InputRole.CATALYST)),
                (tso_cdna,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The open-source implementation permits primer/index choices; no one immutable amplification pair is universal."))
    reviewed.extend(("tso", "tso_reverse_transcriptase", "template_switch",
                     "tso_amp_forward", "tso_amp_reverse", "tso_amp_polymerase",
                     "tso_amplification"))

    def finish(prefix: str, template, label: str):
        transposome = graph.input(prefix + "_transposome", label + " loaded library transposome",
                                  MaterialKind.ADAPTER)
        tn5 = graph.input(prefix + "_tn5", label + " Tn5 transposase",
                          MaterialKind.CATALYST)
        tagmented = graph.product(prefix + "_tagmented", label + " tagmented cDNA",
                                  MaterialKind.MOLECULAR, state=template.state)
        losses = graph.product(prefix + "_losses", label + " non-library Tn5 products",
                               MaterialKind.WASTE)
        graph.react(prefix + "_tagment", "Tagment the " + label + " product",
                    StepKind.TAGMENTATION,
                    ((template, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                     (tn5, InputRole.CATALYST)), (tagmented, losses))
        graph.discard(losses, "does not retain an amplifiable barcode-bearing insert")
        reviewed.extend((prefix + "_transposome", prefix + "_tn5", prefix + "_tagment"))
        _finish_pcr(graph, prefix, tagmented, final, label + ": " + d.endpoint,
                    reviewed, label=label)

    finish("ivt_library", ivt_cdna, "inDrops-2 IVT branch")
    finish("tso_library", tso_cdna, "inDrops-2 TSO-PCR branch")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
