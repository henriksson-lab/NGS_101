from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, FourSUReadout, InputRole, MaterialKind, ProtocolGraph, StepKind, finish_cdna, molecule, role
DEFINITION = Definition(
    "TimeLapse-seq — nucleoside-recoding RNA dynamics",
    "10.1038/nmeth.4582",
    "Schofield et al., Nature Methods (2018)",
    "Oxidative nucleophilic substitution converts 4sU in newly synthesized RNA into a cytidine-like base, producing T-to-C mismatches in ordinary RNA-seq reads.",
    "Standard RNA-seq reads identify newly synthesized molecules from distributed T-to-C conversions.",
    "01_timelapse_seq.html",
    caveat=(
        "Published pulse discrepancy: the original Supplementary Protocol (Springer MOESM3, step 2) "
        "says 100 mM final 4sU for K562 cells, whereas the original Supplementary Note (PMC "
        "supplement-1, 'Important parameters') says the four-hour K562 pulse succeeded at 100 µM "
        "and reports a successful range of 10 µM–1 mM. This page follows the latter reported "
        "experiment and leaves the graph's pulse concentration unspecified; the intended "
        "value of the Protocol's 100 mM entry remains uncertain. The named SMARTer Pico "
        "kit's 2016 guide establishes a ligation-free N6/random-priming, LNA-template-switch, "
        "dual-index PCR workflow, but the study's exact kit revision and index choices are "
        "unverified and the guide labels its RT/TSO oligo mixes proprietary. "
        "The TT-TimeLapse K562 experiment additionally uses Duffy et al. 2015 "
        "MTS-biotin/streptavidin enrichment before recoding; this graph shows the recoding "
        "path shared with the non-enriched experiments."
    ),
)

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=True, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The paper reports paired-end 75 bp for targeted Actb/Gapdh, paired-end 100 bp for transcriptome-wide MEF, and paired-end 150 bp for K562, but does not expose one shared sequencing-primer geometry; kit adapters are shown by role.')))

READOUT = FourSUReadout('oxidation plus 2,2,2-trifluoroethylamine nucleophilic substitution', 'U-to-C mismatch', False)

def four_su_graph(d: Definition, final: Construct) -> ProtocolGraph:
    cfg = READOUT
    source = f'{d.citation}; doi:{d.doi}'
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input('cells', 'living cells during the 4sU pulse', MaterialKind.BIOLOGICAL, state=molecule('living cells', (role('RNA synthesis during pulse'),)), primary=True)
    four_su = graph.input('four_su', '4-thiouridine pulse-label nucleotide', MaterialKind.MOLECULAR)
    polymerase = graph.input('cellular_rna_polymerase', 'cellular RNA polymerases', MaterialKind.CATALYST)
    labelled = graph.product('labelled_rna', '4sU pulse-labelled RNA', MaterialKind.MOLECULAR, state=cfg.labelled())
    graph.react('pulse_label', 'Pulse-label newly synthesized RNA with 4sU', StepKind.TRANSCRIPTION, ((cells, InputRole.TEMPLATE), (four_su, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (labelled,))
    reviewed = ['cells', 'four_su', 'cellular_rna_polymerase', 'pulse_label']
    chemistry = graph.input('recoding_chemistry', cfg.chemistry, MaterialKind.CONDITION)
    selected = graph.product('recoded_rna', 'chemically recoded 4sU RNA', MaterialKind.MOLECULAR, state=cfg.product())
    graph.react('recode', 'Chemically recode 4sU before reverse transcription', StepKind.CHEMICAL_CONVERSION, ((labelled, InputRole.TEMPLATE), (chemistry, InputRole.CONDITION)), (selected,))
    reviewed.extend(('recoding_chemistry', 'recode'))
    finish_cdna(graph, 'library', selected, final, d.endpoint, reviewed,
                rt_label='SMARTer Stranded Total RNA-Seq Pico kit RT primer')
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return four_su_graph(d, final)
