"""Source-owned original TUC-seq boundary; no inherited 4sU population recipe."""
from chemdraw import Construct, Scene, Segment
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
import seqprimers as sp

TITLE='TUC-seq — original thiouridine-to-cytidine conversion'
NOTES='01_tuc_seq.html'
SOURCE='<a href="https://doi.org/10.1002/anie.201707465">Riml et al., Angewandte Chemie (2017)</a>; original publisher/PubMed record. Later TUC-seq DUAL and protocol chapters are version-scoped comparisons.'
SUMMARY='The original TUC-seq record describes osmium-mediated conversion of 4sU into native cytidine, followed by detection in selected transcripts. This native-C product differs from the substituted cytidine analogue in TimeLapse chemistry.'
CAVEAT=('The original 2017 full methods and named supporting-information originals have not been read in this audit: current standard publisher routes return access failures. The model therefore declares the source-positive recoding boundary while leaving sample preparation, library assembly and run primers unassigned. '
        'A later TUC-seq DUAL amplicon workflow is a separate implementation, not a substitute for the original methods. '
        'INFERRED — RNA classes distinguished by pulse-label incorporation and a possible unconverted labelled class express population possibilities, not physical fractions or quantitative conversion yields. The classes are not concatenated into one strand; transcript flanks may predate pulse incorporation.')
FINAL_LIBRARY=Construct((role('original sequencing end unassigned'),role('selected-transcript-derived sequence; length unassigned'),role('other sequencing end and indexing configuration unassigned')),
                        name='original TUC-seq sequencing architecture unassigned')
FINAL_CAPTION='Unassigned original 2017 sequencing architecture. No particular adapter family, index width, RT/PCR route or loaded run primer is inferred from later comparisons.'
SEQ_PRIMERS=()
SEQUENCING_ENDING='The original publication record reports detection of pulse-labelled mRNA from seven selected mammalian genes and endogenous tRNA 4sU. These are distinct experimental applications; no universal acquired run is assigned here.'
SEQUENCING_UNAVAILABLE='The original study-specific sequencing architecture and acquired read layout are not declared. Actual full methods and supporting information remain open source obligations.'
READ_SPECS=(sp.read_spec('Read 1','unassigned','selected-transcript conversion readout; exact original run/read primer unassigned','undetermined'),)

def rna_pool(converted=False):
    sc=Scene()
    sc.strand('old',[role('RNA without a pulse-installed label')],label='RNA without a pulse-installed label')
    mark=Segment('native cytidine formed from 4sU','C','r1') if converted else Segment('pulse-labelled 4sU','4sU',placeholder=True,length_bp=1)
    sc.strand('new',[role('RNA flank before the marked site; synthesis time unassigned'),mark,role('RNA continuation after the marked site; synthesis time unassigned')],
        label='RNA with a representative pulse-derived converted site' if converted else 'RNA with a representative pulse-incorporated 4sU site')
    sc.note('new','Representative site only; conversion completeness and fraction sizes are not asserted.')
    if converted:
        sc.strand('possible_unconverted',[role('possible pulse-labelled RNA with unconverted 4sU',inferred=True)],label='INFERRED — possible unconverted labelled RNA')
    return sc

def protocol_graph():
    graph=ProtocolGraph(TITLE,authority='Original 2017 Riml publisher/PubMed abstract; full methods unassigned')
    cells=graph.input('cells','living cells during metabolic pulse',MaterialKind.BIOLOGICAL,state=molecule('cells',(role('cellular RNA populations'),)),primary=True)
    label=graph.input('four_su','4-thiouridine metabolic label',MaterialKind.MOLECULAR)
    polymerase=graph.input('cellular_rna_polymerase','cellular RNA polymerases',MaterialKind.CATALYST)
    labelled=graph.product('labelled_rna','RNA classes distinguished by pulse-label incorporation',MaterialKind.MOLECULAR,state=rna_pool())
    graph.react('pulse_label','Metabolically label newly synthesized RNA',StepKind.TRANSCRIPTION,
        ((cells,InputRole.TEMPLATE),(label,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(labelled,))
    ready=graph.product('rna_for_conversion','RNA prepared for original conversion; processing unassigned',MaterialKind.MOLECULAR,state=labelled.state,inferred=True)
    graph.react('original_rna_preparation','INFERRED — original RNA handling before conversion unassigned',StepKind.MAGIC,
        ((labelled,InputRole.TEMPLATE),),(ready,),reason='The original 2017 main/SI scientific methods have not been recovered/read; exact extraction and handling order are not borrowed from later DUAL.')
    chemistry=graph.input('recoding_chemistry','source osmium-tetroxide-mediated 4sU-to-native-C chemistry; exact original mixture unassigned',MaterialKind.CONDITION)
    recoded=graph.product('recoded_rna','RNA with representative native cytidine conversion',MaterialKind.MOLECULAR,state=rna_pool(True))
    graph.react('recode','Convert 4sU to native cytidine (original source boundary)',StepKind.CHEMICAL_CONVERSION,
        ((ready,InputRole.TEMPLATE),(chemistry,InputRole.CONDITION)),(recoded,))
    library=graph.product('library','selected-transcript sequencing material; original assembly unassigned',MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    graph.react('library_preparation','INFERRED — original downstream preparation and library assembly unassigned',StepKind.MAGIC,
        ((recoded,InputRole.TEMPLATE),),(library,),reason='The exact original 2017 downstream steps, physical adapter states, indexing and runout are unassigned while its full methods/SI remain uninspected. Later DUAL RT/amplicon/Eurofins PE150 choices are not inherited.')
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(('cells','four_su','cellular_rna_polymerase','pulse_label',
        'original_rna_preparation','recoding_chemistry','recode','library_preparation','library'),graph.authority)
    return graph.require_publishable()
