"""Exact CG000239 Rev A pre-fragmentation cDNA boundary used by Spatial VDJ.

Separate from the Rev F sequencing route. Printed amplified duplex is a source
output, not proof of native first-cycle PCR geometry or the slide release strand.
"""
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import molecule,role
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Disposition,Terminal
from protocol_handoff import export_material
from pcr import PCRGeometry
import illumina
import rt

PROTOCOL='10x-visium-spatial-gene-expression'
VERSION='CG000239 Rev A'
SOURCE_URL='https://assets.ctfassets.net/an68im79xiti/3pyXucRaiKWcscXy3cmRHL/a1ba41c77cbf60366202805ead8f64d7/CG000239_VisiumSpatialGeneExpression_UserGuide_Rev_A.pdf'
SOURCE='10x Genomics CG000239 Rev A, Appendix printed p75 and Steps4.2–4.3 printed pp49–50; exact source cited by Engblom et al.2023'
FORWARD=illumina.TRUSEQ_READ1[-22:]
REVERSE=rt.SMART_HANDLE[:22]
TSO_DNA=rt.SMART_HANDLE+'ACAT'  # Rev A PN-3000228 DNA portion, before rGrGrG
SPATIAL=feature('visium_rev_a.spatial','cell_barcode','whitelist',whitelist='10x-visium-spot-barcode-set',note='Spatial capture spot, not individual cell')
UMI=feature('visium_rev_a.umi','umi','random')
LIMIT='Source-printed amplified duplex pattern; actual transcript, spot/UMI values and recovery are variable. Native first-cycle PCR and earlier slide-release strand chemistry remain unvalidated. No fragmentation, ligation adapter, P5/P7 or sample index is present at this handoff.'


def full_cdna():
    return Construct([Segment('partial Read 1',FORWARD,'r1'),
        Segment('spatial barcode','N'*16,'cbc',placeholder=True,length_bp=16,feature=SPATIAL),
        Segment('UMI','N'*12,'umi',placeholder=True,length_bp=12,feature=UMI),
        Segment('poly(dT)30','T'*30),Segment('anchoring V','V',placeholder=True,length_bp=1,bottom='b'),Segment('anchoring N','N',placeholder=True,length_bp=1,bottom='n'),
        Segment('transcript cDNA','X'*24,placeholder=True,note='Unknown full transcript; glyph width is not physical length'),
        Segment('RT-added CCC','CCC'),Segment('copied TSO DNA',revcomp(TSO_DNA),'tso')],name='CG000239 Rev A full-length amplified spatial cDNA before fragmentation')


def cdna_primer_scene():
    sc=Scene();sc.strand('PN2000089_forward',[Segment('Partial Read 1',FORWARD,'r1')]);sc.strand('PN2000089_reverse',[Segment('Partial TSO',REVERSE,'tso')]);return sc


def intermediate_graph():
    g=ProtocolGraph('Visium CG000239 Rev A full-length cDNA handoff',authority=SOURCE,required_output=Disposition.STORED)
    recovered=g.input('recovered_cdna','Rev A Step3.2 recovered cDNA, before qPCR-guided bulk amplification',K.MOLECULAR,primary=True,
        state=molecule('Recovered cDNA strand identity unresolved',(role('actual released spatial cDNA; strand identity not assigned'),)),
        chemical_state={'stage':'recovered_cdna','release_strand':'unassigned'},source=SOURCE+' Earlier slide prefix is outside this exported preparation boundary.')
    f=g.input('cdna_forward','cDNA Primers PN-2000089 forward',K.OLIGO,state=Construct([Segment('Partial Read 1',FORWARD)]),source=SOURCE)
    r=g.input('cdna_reverse','cDNA Primers PN-2000089 reverse',K.OLIGO,state=Construct([Segment('Partial TSO',REVERSE)]),source=SOURCE)
    pol=g.input('amp_mix','Rev A Amp Mix PN-2000047/2000103',K.CATALYST,source=SOURCE)
    programme=g.input('programme','Step4.2 Cq-guided cycles;98C3min,98C15s/63C20s/72C1min;72C1min final',K.CONDITION,source=SOURCE)
    amplified=g.product('amplified_cdna','Source-printed full-length spatial cDNA duplex',K.MOLECULAR,state=Scene.duplex(list(full_cdna())),source=SOURCE)
    g.react('cdna_amplification','Amplify full-length spatial cDNA before fragmentation',S.PCR,((recovered,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST),(programme,R.CONDITION)),(amplified,),
        geometry=PCRGeometry.unavailable(recovered.state,reason='Exact source primer orders and amplified duplex are printed. Native released strand identity and its first-cycle primer geometry remain unassigned; the output is not substituted for that incoming substrate.'))
    beads=g.input('beads','Step4.3 SPRIselect0.6X, two80%ethanol washes; transfer40ulEB',K.SELECTION,source=SOURCE)
    product=g.product('full_length_cdna','Purified full-length spatially barcoded cDNA before fragmentation',K.MOLECULAR,state=amplified.state,
        chemical_state={'stage':'amplified_full_length_cdna','polymer':'DNA','strandedness':'duplex_pattern','fragmented':False,'spatial_barcode_nt':16,'umi_nt':12,'sample_index_present':False,'cdna_primers_part_number':'PN-2000089','validation_scope':LIMIT},source=SOURCE)
    waste=g.product('cleanup_waste','Actual nonretained cleanup fractions; molecular composition unassigned',K.WASTE)
    g.react('cdna_cleanup','Purify amplified cDNA, retaining full-length molecules',S.CLEANUP,((amplified,R.TEMPLATE),(beads,R.SELECTION)),(product,waste),changes_molecule=False);g.discard(waste,'Physical bead/supernatant workup fractions; perfect specificity not asserted.')
    g.mark_terminal(Terminal(product.id,Disposition.STORED,'Reusable Rev A intermediate at Step4.3m, before Step5.1 fragmentation.'))
    g.declare_reviewed_scope(('recovered_cdna','cdna_forward','cdna_reverse','amp_mix','programme','cdna_amplification','beads','cdna_cleanup','full_length_cdna'),SOURCE)
    g.require_publishable()
    export_material(g,'full_length_cdna_rev_a',product,protocol=PROTOCOL,version=VERSION,page='visium.html',source=SOURCE,limitation=LIMIT)
    return g


def full_length_export():return intermediate_graph().exports['full_length_cdna_rev_a']
