from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, FourSUReadout, InputRole, MaterialKind, ProtocolGraph, StepKind, finish_cdna, molecule, role
DEFINITION = Definition("SNU-Seq — single-nucleotide-resolution 4sU sequencing", "10.1093/nar/gkag652", "Gerlevik et al., Nucleic Acids Research (2026)", "Biotin-selected intact 4sU RNA is polyadenylated with bacterial poly(A) polymerase and captured by a 3′ RNA-seq workflow to preserve native endpoint information.", "The 3′-anchored read reports the native end of an intact pulse-labelled RNA molecule.", "01_snu_seq.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=False, platform="platform-specific")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The paper names Ion Torrent and Illumina FWD QuantSeq branches without assigning a kit revision or barcode to each sample. The official 2018 Ion and 2021 FWD guides print their distinct finished library shells and sequencing entries; this combined schematic keeps the platform-specific ends by role.')))
CAVEAT = ('Lexogen QuantSeq Ion Torrent and Illumina FWD are distinct library branches. '
          'Official <a href="https://www.lexogen.com/wp-content/uploads/2018/06/012UG036V0130_QuantSeq-IonTorrent.pdf">Ion guide 012UG036V0130</a> '
          '(Appendix E) and <a href="https://www.lexogen.com/wp-content/uploads/2021/08/015UG009V0260_QuantSeq_Illumina_2021-08-16.pdf">Illumina FWD guide 015UG009V0260</a> (Appendix G) '
          'print their finished library architecture. The original paper and GEO do not '
          'assign an exact guide revision or barcode well to every SNU-Seq library; '
          'the final panel represents their shared assay insert and leaves kit-specific '
          'ends symbolic.')

READOUT = FourSUReadout('thiol biotinylation, streptavidin purification and bacterial poly(A) polymerase', 'native RNA ends', True)

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
    biotin = graph.input('thiol_biotin', '4sU thiol-biotinylation reagent', MaterialKind.MOLECULAR)
    condition = graph.input('biotinylation_condition', 'thiol-biotinylation condition', MaterialKind.CONDITION)
    biotinylated = graph.product('biotinylated_rna', 'biotinylated 4sU RNA', MaterialKind.MOLECULAR, state=labelled.state)
    graph.react('biotinylate', 'Thiol-biotinylate incorporated 4sU', StepKind.CHEMICAL_CONVERSION, ((labelled, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL), (condition, InputRole.CONDITION)), (biotinylated,))
    streptavidin = graph.input('streptavidin', 'streptavidin RNA-capture reagent', MaterialKind.SELECTION)
    enriched = graph.product('enriched_rna', 'streptavidin-enriched 4sU RNA', MaterialKind.MOLECULAR, state=cfg.product())
    unlabelled = graph.product('unlabelled_rna', 'unlabelled RNA background', MaterialKind.WASTE)
    graph.react('enrich', 'Affinity-purify biotinylated 4sU RNA', StepKind.AFFINITY_CAPTURE, ((biotinylated, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, 'does not carry the pulse-installed 4sU affinity label')
    reviewed.extend(('thiol_biotin', 'biotinylation_condition', 'biotinylate', 'streptavidin', 'enrich'))
    pap = graph.input('poly_a_polymerase', 'bacterial poly(A) polymerase', MaterialKind.CATALYST)
    atp = graph.input('tailing_atp', 'ATP for RNA polyadenylation', MaterialKind.MOLECULAR)
    selected = graph.product('selected_rna', 'polyadenylated intact 4sU RNA', MaterialKind.MOLECULAR,
                             state=molecule('polyadenylated intact 4sU RNA',
                                            (role('native intact 4sU-labelled RNA'),
                                             role('bacterial poly(A) tail'))))
    graph.react('polyadenylate', 'Polyadenylate intact enriched RNA', StepKind.TAILING, ((enriched, InputRole.TEMPLATE), (atp, InputRole.STRUCTURAL), (pap, InputRole.CATALYST)), (selected,))
    reviewed.extend(('poly_a_polymerase', 'tailing_atp', 'polyadenylate'))
    finish_cdna(graph, 'library', selected, final, d.endpoint, reviewed,
                rt_label='Lexogen QuantSeq 3-prime oligo(dT) capture primer')
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    return four_su_graph(d, final)
