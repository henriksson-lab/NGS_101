"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
import seqprimers as sp
from rna_end_readout import (StopReadout, BaseConversionReadout, SelectiveCleavageReadout,
                             _finish_pcr, _apply_linear_step)


TITLE = 'GLORI — absolute m6A quantification by A-to-I conversion'
NOTES = '01_glori.html'
SOURCE = 'Liu et al., Nature Biotechnology (2023), <a href="https://doi.org/10.1038/s41587-022-01487-9">doi:10.1038/s41587-022-01487-9</a>; detailed protocol: Shen et al., Nature Protocols (2024), <a href="https://doi.org/10.1038/s41596-023-00937-1">doi:10.1038/s41596-023-00937-1</a>.'
SUMMARY = 'Glyoxal-assisted nitrite treatment converts unmethylated adenosine to inosine while m6A remains adenosine; A retention therefore reports absolute methylation stoichiometry.'
CAVEAT = "This page follows the modified-eCLIP library branch of the Liu 2023 Methods ('RNA protection, deamination and deprotection'; 'Library preparation and sequencing'), which print the 3′ linker, RT primer and UMI 5′ adapter with their 5′ adenylation/phosphate and 3′ C3-spacer blocks. The library PCR primers are named only as 'P5 and Illumina-indexed P7 primers' (2023) or NEBNext Universal and Index primers (2024 Nature Protocols Step 92); their bases are not printed in either text, so the outer P5/P7 arms and the sample-index window stay inferred. The authors also validated NEBNext Small RNA and SMARTer Pico v3 kit alternatives, not drawn."
# Role order derived 2026-10-17 from the 2024 Supplementary Table 1 oligos: the RNA
# adaptor is the reverse complement of the TruSeq Read 1 primer 3' end (the RT primer
# is that 3' end), and the DNA adaptor arm is the reverse complement of the TruSeq
# Read 2 primer 3' end.  Hence the cDNA strand reads Read1 site - insert - N10 UMI -
# Read2 site.  P5 on the Read 1 side is the standard Illumina convention, inferred:
# the library PCR primers are not printed in the inspected sources.
FINAL_LIBRARY = Construct((role('P5', inferred=True), role('Read 1 site from RNA adaptor / RT primer', tag='r1'), role('A-to-I-converted cDNA (antisense to RNA)'), role('10-nt UMI', tag='umi', length_bp=10, molecular_feature=feature('umi', 'umi', 'random')), role('Read 2 site from 3-prime DNA adaptor', tag='r2'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7', inferred=True)), name='GLORI modified-eCLIP library')
FINAL_CAPTION = 'Modified-eCLIP branch, role order derived from the published adaptors. The 10-nt UMI lies between the cDNA 3-prime end and the Read 2 site; converted A is read as G in the RNA sense while m6A is retained as A.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'Liu 2023 Methods: the eCLIP-branch libraries were sequenced on Illumina HiSeq X Ten as paired-end 2 × 150 base reads; the 2024 Nature Protocols Step 121 gives the same paired-end 150-bp setting. Only read R2 was analysed, because R2 is sense to the original RNA; the 10-bp UMI at the 5′ end of R2 is used for duplicate removal and then trimmed. This matches the model, where Read 2 is primed from the UMI adapter arm. GEO GSE210563 additionally lists NovaSeq 6000 runs, and the deposited SRA runs (PRJNA866163) are single-end 150-base spots, consistent with only R2 being deposited. The sample-index read length is not stated.'
SEQUENCING_UNAVAILABLE = "The 2023 Methods print the 3′ linker, RT primer and UMI adapter. The library PCR primers ('P5 and Illumina-indexed P7', 2023; NEBNext Universal/Index primers, 2024) and the sequencing primers are not printed in the inspected texts, so primer landing on the final library is not reconstructed."

def oligos():
    return ("Liu 2023 Methods, 3′ linker — 5′rAPP-AGATCGGAAGAGCGTCGTG-3SpC3 (2024 Supplementary Table 1 prints rAppAGATCGGAAGAGCGTCGTG; 2024 Reagents state the 5′ adenylation and 3′ C3 chain)",
            'Liu 2023 Methods, reverse transcription primer — ACACGACGCTCTTCCGATCT',
            "Liu 2023 Methods, 5′ adapter (UMI) — 5′Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3SpC3 (2024 '3′ DNA adaptor'; Reagents state the 5′ phosphate and 3′ C3 spacer)")

def workflow() -> Workflow:
    readout = BaseConversionReadout("glyoxal plus nitrite", "adenosine", "inosine", "G", "m6A")
    return chain(molecule("fragmented poly(A) RNA", (role("RNA containing A and m6A"),)), (
        ("Protect guanosines with glyoxal",
         molecule("glyoxal-protected RNA", (role("RNA"), role("protected G"), role("A and m6A exposed"))),
         "G protection suppresses unwanted guanosine deamination during the nitrite reaction."),
        ("Deaminate unmethylated A with nitrite",
         molecule("GLORI-converted RNA", (role(f"{readout.converted_base} converted to {readout.product}"), role(f"{readout.protected_base} retained as A"))),
         "Liu 2023 Methods: 750 mM NaNO2, 40 mM MES pH 6.0 and glyoxal for 8 h at 16 °C (an earlier repository note's extra 4 °C overnight incubation is not in the Methods or 2024 Step 23)."),
        ("Heat-deprotect guanosines and repair RNA ends",
         molecule("deprotected converted RNA", (role("A-to-I converted RNA with intact m6A"), role("ligatable termini"))),
         "Liu 2023 Methods: deprotection buffer (triethylammonium acetate/formamide) at 95 °C for 10 min; T4 PNK then removes 3′ phosphates before linker ligation."),
        ("Ligate the modified-eCLIP 3-prime RNA adapter",
         molecule("GLORI adapter-ligated RNA", (role("converted RNA"),
                  Segment("published 3-prime RNA adapter", "AGATCGGAAGAGCGTCGTG", "r1"))),
         "Liu 2023 Methods: 5′-adenylated, 3′-C3-blocked linker ligated by T4 RNA ligase 2 truncated KQ; excess linker removed with 5′ deadenylase and RecJf."),
        ("Reverse-transcribe with the published RT primer",
         molecule("GLORI cDNA", (Segment("published RT primer", "ACACGACGCTCTTCCGATCT", "r1"),
                  role("converted cDNA"))),
         "Liu 2023 Methods: SuperScript III with this primer; excess primer digested by exonuclease I."),
        ("Ligate the ten-base-UMI DNA adapter to cDNA",
         molecule("UMI-tagged GLORI cDNA", (Segment("published RT primer", "ACACGACGCTCTTCCGATCT", "r1"),
                  role("converted cDNA"), role("10-nt UMI"),
                  Segment("published 3-prime cDNA adapter arm", "AGATCGGAAGAGCACACGTCTG", "r2"))),
         "Liu 2023 Methods: 5′-phosphorylated, 3′-C3-blocked adapter with ten degenerate UMI bases, ligated to the cDNA 3′ end on silane beads with T4 RNA ligase 1 (high concentration)."),
        ("INFERRED — amplify with P5 and indexed P7 primers (bases not printed) and sequence",
         molecule("GLORI library", tuple(FINAL_LIBRARY)),
         f"Unmodified A is read as {readout.read_as} in the RNA sense; m6A remains A, enabling site-level stoichiometry. The 2023 Methods name P5 and Illumina-indexed P7 primers without bases; 2 × 150 paired-end HiSeq X Ten, R2 (RNA-sense, UMI first) analysed."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

M23 = "Liu 2023 Methods"
LINKER = Construct((Segment("3′ linker (5′ rApp, 3′ SpC3)", "AGATCGGAAGAGCGTCGTG", "r1"),), name="GLORI 3′ linker")
RT_PRIMER = Construct((Segment("RT primer", "ACACGACGCTCTTCCGATCT", "r1"),), name="GLORI RT primer")
UMI_ADAPTER = Construct((Segment("10-nt UMI", "NNNNNNNNNN", "umi", feature=feature("umi", "umi", "random")),
                         Segment("adapter arm (3′ SpC3)", "AGATCGGAAGAGCACACGTCTG", "r2")),
                        name="GLORI 5′ adapter (5′ phosphate)")

# (identifier, label, kind, keyword facts for _apply_linear_step)
STEPS = (
    ("fragment", "Fragment mRNA", StepKind.FRAGMENTATION,
     dict(reagent_label=f"94 °C for 4 min in NEB fragmentation buffer (E6150S) ({M23})")),
    ("protect_g", "Protect guanosines with glyoxal", StepKind.CHEMICAL_CONVERSION,
     dict(reagent_label=f"200 ng RNA with 6 µl 8.8 M glyoxal and 20 µl DMSO in 40 µl, 50 °C 30 min; "
                        f"then 10 µl saturated H3BO3, 50 °C 30 min ({M23})")),
    ("deaminate", "Deaminate unmethylated adenosine with nitrite", StepKind.CHEMICAL_CONVERSION,
     dict(reagent_label=f"750 mM NaNO2, 40 mM MES pH 6.0, 10 µl 8.8 M glyoxal; 8 h at 16 °C; "
                        f"ethanol precipitation ({M23})")),
    ("deprotect", "Heat-deprotect guanosines", StepKind.THERMAL_TREATMENT,
     dict(reagent_label=f"deprotection buffer (triethylammonium acetate pH 8.6, formamide), 95 °C "
                        f"10 min; ethanol precipitation and RNA Clean & Concentrator ({M23})")),
    ("repair", "Remove RNA 3′ phosphates with T4 PNK", StepKind.DEPHOSPHORYLATION,
     dict(catalyst_label=f"T4 PNK (NEB M0201S), 3′-end dephosphorylation; ethanol precipitation ({M23})")),
    ("ligate", "Ligate the pre-adenylated 3′ linker", StepKind.LIGATION,
     dict(reagent_label=f"3′ linker 5′rAPP-AGATCGGAAGAGCGTCGTG-3SpC3 ({M23})", reagent_state=LINKER,
          catalyst_label=f"T4 RNA ligase 2, truncated KQ (NEB M0373S) ({M23})",
          geometry=AdapterLigationGeometry.unavailable(
              "Single-stranded ligation of a 5′-adenylated, 3′-C3-blocked DNA linker to RNA 3′-OH "
              "ends of random fragments; the RNA-side bases are not defined, so no base junction is drawn."))),
)
STEPS_AFTER_RT = (
    ("ligate_cdna", "Ligate the UMI adapter to the cDNA 3′ end on silane beads", StepKind.LIGATION,
     dict(reagent_label=f"5′Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3SpC3 ({M23})", reagent_state=UMI_ADAPTER,
          catalyst_label=f"T4 RNA ligase 1, high concentration (NEB M0437M), on MyOne Silane beads ({M23})",
          geometry=AdapterLigationGeometry.unavailable(
              "Single-stranded ligation of the 5′-phosphorylated UMI adapter to random cDNA 3′ ends; "
              "the cDNA-side bases are not defined."))),
)


def _remove_excess(graph, current, identifier, label, catalyst_label, waste_label, reviewed):
    enzyme = graph.input(identifier + "_enzyme", catalyst_label, MaterialKind.CATALYST)
    kept = graph.product(identifier + "_product", label + " retained material", MaterialKind.MOLECULAR,
                         state=current.state)
    waste = graph.product(identifier + "_excluded", waste_label, MaterialKind.WASTE)
    graph.react(identifier, label, StepKind.DIGESTION,
                ((current, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (kept, waste),
                changes_molecule=False)
    graph.discard(waste, "digested")
    reviewed.extend((identifier + "_enzyme", identifier))
    return kept


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    start = graph.input("input", "oligo(dT)-purified mRNA with m6A spike-in probes", MaterialKind.BIOLOGICAL,
                        state=molecule("source RNA", (role("input RNA"),)), primary=True)
    reviewed = ["input"]
    current = start
    for identifier, label, kind, facts in STEPS:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed, **facts)
    current = _remove_excess(graph, current, "linker_cleanup", "Deadenylate and digest excess linker",
                             f"5′ deadenylase (NEB M0331S) and RecJf ({M23}); ethanol precipitation",
                             "deadenylated, digested excess 3′ linker", reviewed)
    current = _apply_linear_step(graph, current, "rt", "Reverse-transcribe converted RNA",
                                 StepKind.REVERSE_TRANSCRIPTION, reviewed,
                                 reagent_label=f"RT primer ACACGACGCTCTTCCGATCT, 2 µM × 1 µl ({M23})",
                                 reagent_state=RT_PRIMER,
                                 catalyst_label=f"SuperScript III (Invitrogen 18080093) in FS buffer ({M23})")
    current = _remove_excess(graph, current, "primer_cleanup", "Digest excess RT primer",
                             f"exonuclease I (NEB M0293S); silane-bead purification ({M23})",
                             "digested excess RT primer", reviewed)
    for identifier, label, kind, facts in STEPS_AFTER_RT:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed, **facts)
    forward = graph.input("library_forward", "P5 primer (2023 Methods; bases not printed) — 2024 Step 92: "
                          "NEBNext Universal PCR Primer for Illumina", MaterialKind.OLIGO, inferred=True)
    reverse = graph.input("library_reverse", "Illumina-indexed P7 primer (2023 Methods; bases not printed) — "
                          "2024 Step 92: NEBNext Index Primer (E7335S)", MaterialKind.OLIGO, inferred=True)
    polymerase = graph.input("library_polymerase", "PCR polymerase not named in the 2023 Methods; 2024 Step 92–93: "
                             "NEBNext Ultra II Q5 Master Mix, 10–12 cycles (input-dependent)", MaterialKind.CATALYST)
    amplified = graph.product("library_amplified", "amplified GLORI library", MaterialKind.MOLECULAR,
                              state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify with P5 and Illumina-indexed P7 primers", StepKind.PCR,
                ((current, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The 2023 Methods name P5 and Illumina-indexed P7 primers and the "
                                          "2024 protocol names NEBNext primers, but neither prints their bases."))
    beads = graph.input("library_beads", f"1.8× AMPure XP beads ({M23})", MaterialKind.SELECTION)
    cleaned = graph.product("library_cleaned", "bead-purified GLORI library", MaterialKind.MOLECULAR,
                            state=FINAL_LIBRARY)
    unbound = graph.product("library_bead_excluded", "primers and short products not retained by 1.8× beads",
                            MaterialKind.WASTE)
    graph.react("library_cleanup", "Purify the PCR product with 1.8× XP beads", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (cleaned, unbound),
                changes_molecule=False)
    graph.discard(unbound, "not retained by the beads")
    gel = graph.input("library_gel", f"8% TBE gel, 165–250 bp excised ({M23}; 2024: 170–250 bp)",
                      MaterialKind.SELECTION)
    library = graph.product("library_library", "GLORI sequencing library", MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    off = graph.product("library_gel_excluded", "fragments outside 165–250 bp, including adapter dimers",
                        MaterialKind.WASTE)
    graph.react("library_size_select", "Select 165–250 bp on an 8% TBE gel", StepKind.SELECTION,
                ((cleaned, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, off),
                changes_molecule=False)
    graph.discard(off, "outside the excised gel region")
    graph.sequence(library, SequencingPlan(limitation=SEQUENCING_ENDING))
    reviewed.extend(("library_forward", "library_reverse", "library_polymerase", "library_pcr",
                     "library_amplified", "library_beads", "library_cleanup", "library_gel",
                     "library_size_select", "library_library"))
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
