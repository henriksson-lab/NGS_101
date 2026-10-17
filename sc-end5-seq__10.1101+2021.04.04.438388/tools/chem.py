"""Protocol-owned source, molecular workflow and reaction graph."""
from chemdraw import Construct, Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('sc-end5-seq', '10.1101/2021.04.04.438388', 'Moody et al., bioRxiv (2021)', 'Oligo(dT) or random priming followed by 10x template switching retains nucleotide-resolution RNA 5′ ends with cell barcodes and UMIs.', "The insert's genomic 5′ boundary reports a candidate transcription start site; cell barcode and UMI support single-cell counting.", '01_sc_end5_seq.html')
PAYLOAD = (CELL, UMI, role("RNA 5-prime insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
RT_HANDLE = "AAGCAGTGGTATCAACGCAGAGTAC"
DT_RT_PRIMER = RT_HANDLE + "T" * 30 + "VN"
RANDOM_RT_PRIMER = RT_HANDLE + "N" * 6
# Preprint Methods (v2 PDF text line ~1810): the 13-nt TS-oligo sequence aligned
# against genomic sequence to remove strand invaders.  It is printed as read DNA,
# not as the bead oligo with its terminal ribonucleotides or its 5' portion.
TS_OLIGO_ANALYSIS_3PRIME = "TTTCTTATATGGG"
# Preprint Methods line ~1421-1431: cDNA PCR cycles per variant.
CDNA_PCR_CYCLES = {"sc-end5-dT": 12, "sc-end5-rand": 14}
CAVEAT = ("The preprint prints both sc-end5 RT primers. Its nearby ACATrGrGrG "
          "template-switch oligo belongs to the sc-end3-dT control, not the sc-end5 "
          "gel bead; for sc-end5 analysis the preprint prints only the 13-nt TS-oligo "
          "3′ sequence TTTCTTATATGGG (as DNA). Library steps after cDNA PCR are "
          "delegated to the kit guide and drawn as an unresolved boundary. A later 10x CG000086 Rev K guide gives same-kit bead architecture, "
          "but the experiment's exact kit revision and chosen i7 well are unresolved. "
          "The reported run used 75-base paired reads on a HiSeq 2500.")
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The original kit revision and selected i7 well are not identified by the preprint. The later guide's primer placement is a comparator, not proof of the experiment's exact run."

def workflow():
    return chain(molecule("single-cell RNA in a 10x droplet", (role("poly(A) and non-poly(A) RNA"),)), (
        ("Prime first-strand synthesis with oligo(dT) or random hexamers",
         molecule("nascent first strand", (role("transcript cDNA"),)),
         "sc-end5-dT and sc-end5-rand are explicit priming variants of the named method."),
        ("Template-switch at the RNA 5′ end onto the gel-bead oligo",
         molecule("5-prime-end cDNA", (CELL, UMI, role("template-switch junction"), role("transcript cDNA"))),
         "Barcode and UMI become associated with a nucleotide-resolution transcript 5′ end."),
        ("Amplify and make the sequencing library",
         molecule("sc-end5-seq library", (CELL, UMI, role("RNA 5-prime insert"), SAMPLE)),
         "The genomic end of the insert reports a candidate transcription start site."),
    ))


def oligos():
    return (
        "sc-end5-dT RT primer, 5′-" + RT_HANDLE + "-T(30)-VN-3′ (preprint Methods)",
        "sc-end5-rand RT primer, 5′-" + RANDOM_RT_PRIMER + "-3′ (17 µM; preprint Methods)",
        "TS-oligo 3′ sequence used for strand-invader filtering, 5′-…" + TS_OLIGO_ANALYSIS_3PRIME
        + "-3′ (printed as DNA read sequence; full bead oligo and ribonucleotides not printed)",
    )


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "single_cell_rna", "single-cell poly(A) and non-poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("single-cell RNA in a 10x droplet", (
            role("poly(A) and non-poly(A) RNA"),)), primary=True)
    dt_rna = graph.product("dt_input", "RNA assigned to sc-end5-dT",
                           MaterialKind.MOLECULAR, state=rna.state)
    random_rna = graph.product("random_input", "RNA assigned to sc-end5-rand",
                               MaterialKind.MOLECULAR, state=rna.state)
    graph.react("split_priming_variants", "Separate oligo(dT) and random-hexamer variants",
                StepKind.SPLIT, ((rna, InputRole.TEMPLATE),), (dt_rna, random_rna),
                changes_molecule=False)
    reviewed = ["single_cell_rna", "split_priming_variants"]

    def branch(prefix: str, template, primer_label: str, primer_bases: str, label: str):
        primer = graph.input(prefix + "_primer", primer_label, MaterialKind.OLIGO,
                             state=Construct((Segment(primer_label + " handle", RT_HANDLE),
                                              Segment(primer_label + " variable tail",
                                                      primer_bases[len(RT_HANDLE):],
                                                      placeholder=True))))
        rt = graph.input(prefix + "_rt_enzyme", "10x template-switching reverse transcriptase",
                         MaterialKind.CATALYST)
        first = graph.product(prefix + "_first_strand", label + " nascent first strand",
                              MaterialKind.MOLECULAR,
                              state=molecule("nascent first strand", (
                                  role("transcript cDNA"),)))
        graph.react(prefix + "_rt", "Prime first-strand synthesis with " + primer_label
                    + " in droplets (53 °C, 45 min)",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (first,))
        bead_tso = graph.input(prefix + "_bead_tso", "cell-barcode/UMI gel-bead template-switch oligo",
                               MaterialKind.OLIGO)
        switched = graph.product(
            prefix + "_switched", label + " 5-prime-end cDNA",
            MaterialKind.MOLECULAR,
            state=molecule("5-prime-end cDNA", (
                CELL, UMI, role("template-switch junction"), role("transcript cDNA"))))
        graph.react(prefix + "_template_switch", "Template-switch at the RNA 5-prime end",
                    StepKind.TEMPLATE_SWITCHING,
                    ((first, InputRole.TEMPLATE), (bead_tso, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (switched,))
        reviewed.extend((prefix + "_primer", prefix + "_rt_enzyme", prefix + "_rt",
                         prefix + "_bead_tso", prefix + "_template_switch"))
        cycles = CDNA_PCR_CYCLES[label]
        # The preprint names one "cDNA primer mix"; splitting it into a primer for
        # each cDNA end is an inference required by exponential PCR, not a source fact.
        cdna_primers = graph.input(prefix + "_cdna_primers", "kit cDNA primer mix: bead-oligo-end primer (inferred component; bases not printed)",
                                   MaterialKind.OLIGO, inferred=True)
        cdna_primers_rt = graph.input(prefix + "_cdna_primers_rt_end", "kit cDNA primer mix: RT-primer-end primer (inferred component; bases not printed)",
                                      MaterialKind.OLIGO, inferred=True)
        cdna_polymerase = graph.input(prefix + "_cdna_polymerase", "kit cDNA amplification polymerase",
                                      MaterialKind.CATALYST)
        amplified = graph.product(
            prefix + "_amplified_cdna", label + " amplified 5-prime-end cDNA",
            MaterialKind.MOLECULAR,
            state=molecule("amplified 5-prime-end cDNA", (
                CELL, UMI, role("template-switch junction"), role("transcript cDNA"))))
        graph.react(prefix + "_cdna_pcr", f"Amplify cDNA with the kit cDNA primer mix ({cycles} cycles)",
                    StepKind.PCR,
                    ((switched, InputRole.TEMPLATE), (cdna_primers, InputRole.STRUCTURAL),
                     (cdna_primers_rt, InputRole.STRUCTURAL), (cdna_polymerase, InputRole.CATALYST)), (amplified,),
                    geometry=PCRGeometry.unavailable(
                        switched.state,
                        reason="The preprint names the kit cDNA primer mix but prints neither its bases nor its binding sites."))
        prepared = graph.product(
            prefix + "_adapter_ready", label + " kit-processed library precursor",
            MaterialKind.MOLECULAR,
            state=molecule("kit-processed library precursor", (
                CELL, UMI, role("RNA 5-prime insert"))))
        graph.react(prefix + "_kit_steps", "Standard kit library steps after cDNA amplification",
                    StepKind.MAGIC, ((amplified, InputRole.TEMPLATE),), (prepared,),
                    reason="Preprint Methods say only 'followed by the standard steps according to "
                           "manufacturer's instructions' for the initial 5-prime kit; the guide revision "
                           "used is not identified, so fragmentation, adapter attachment and cleanups are not drawn.")
        reviewed.extend((prefix + "_cdna_primers", prefix + "_cdna_primers_rt_end", prefix + "_cdna_polymerase",
                         prefix + "_cdna_pcr", prefix + "_kit_steps"))
        _finish_pcr(graph, prefix, prepared, final, label + ": " + d.endpoint,
                    reviewed, label=label)

    branch("dt", dt_rna, "oligo(dT) RT primer", DT_RT_PRIMER, "sc-end5-dT")
    branch("random", random_rna, "random-hexamer RT primer", RANDOM_RT_PRIMER, "sc-end5-rand")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
