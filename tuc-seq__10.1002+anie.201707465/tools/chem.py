from nascent_protocols import make_spec

from nascent_protocols import Construct, Definition, FourSUReadout, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, molecule, role
DEFINITION = Definition("TUC-seq — thiouridine-to-cytidine conversion sequencing", "10.1002/anie.201707465", "Riml et al., Angewandte Chemie (2017)", "Osmium-mediated chemistry converts incorporated 4sU into native cytidine, encoding newly transcribed RNA as U-to-C sequence changes.", "Targeted reads of selected transcripts distinguish labelled molecules through U-to-C conversions; the original library-end chemistry has not been recovered.", "01_tuc_seq.html", "The 2017 supplement and full methods remain inaccessible. The later same-author TUC-seq DUAL supplement describes targeted RT and amplicon sequencing, but cannot establish the exact 2017 library construction.")

FINAL_LIBRARY = Construct((role("unresolved left sequencing end"), role("target-derived cDNA"),
                           role("unresolved right sequencing end")), name=DEFINITION.title + " library")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The accessible 2017 abstract establishes targeted sequencing but does not specify the library adapters, index, read primer, or read layout. The original supporting information is inaccessible.')))

READOUT = FourSUReadout('osmium-tetroxide-mediated 4sU-to-cytidine conversion', 'U-to-C mismatch', False)

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
    library = graph.product('library', 'target-derived sequencing library; assembly unresolved',
                            MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react('library_preparation', 'Prepare target-derived sequencing library (details unresolved)',
                StepKind.MAGIC, ((selected, InputRole.TEMPLATE),), (library,),
                reason='The 2017 full methods and both publisher supplements are inaccessible. The later TUC-seq DUAL supplement describes random-hexamer RT, gene-specific PCR and vendor library preparation, but does not establish the original 2017 assembly or adapter ends.')
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    reviewed.extend(('library_preparation', 'library'))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return four_su_graph(d, final)
