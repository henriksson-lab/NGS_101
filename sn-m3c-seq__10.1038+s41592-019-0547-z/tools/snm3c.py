"""Lee2019 sn-m3C-seq: owned contact front end and source-qualified2018 native finish."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, feature, complement_segments, revcomp
from restriction import DPNII
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, SequencingPlan, Terminal, Disposition
import seqprimers as sp
TITLE='sn-m3C-seq — contacts and DNA methylation from one nucleus'
NOTES='01_sn-m3c-seq.html'
SOURCE='Lee et al., Nature Methods2019, <a href="https://doi.org/10.1038/s41592-019-0547-z">doi:10.1038/s41592-019-0547-z</a>; published PMC6765423 SupplementalMethods Hi-C/3C, FANS, library preparation and data processing.'
CITATION=SOURCE
SUMMARY='Fixed chromatin is DpnII-digested and proximity-ligated without biotin fill-in. Sorted nuclei enter bisulfite conversion and the explicitly inherited Luo2018 snmC-seq2 native finishing procedure.'
CAVEAT='The inherited2018 scope is random-priming/Adaptase/indexing/recovery, applied here to contact DNA. Native index/random-primer members and Adaptase end geometry remain unassigned. PE150 and computational25-base5′/3-base3′ trims do not establish molecular adapter lengths or exact primer landing. The original contact motif can become physical U/C at bisulfite-converted cytosines; its unconverted drawing is not the pre-PCR substrate.'
WORKFLOW=SimpleNamespace(amplification='bisulfite-PCR',marking='none')
RANDOM_PRIMER_TRIM_NT=25
ADAPTASE_TAIL_TRIM_NT=3
RP_HANDLE='TTCCCTACACGACGCTCTTCCGATCT'
CELL_FEATURE=feature('snm3c_inline','cell_barcode','unknown',group='cell_identity',note='2018-family inline cell identifier; selected native member unassigned.')
I5_FEATURE=feature('snm3c_i5','sample_index','unknown',group='cell_identity',note='Dual library-pool index; member/width allocation unassigned for this study.')
I7_FEATURE=feature('snm3c_i7','sample_index','unknown',group='cell_identity',note='Dual library-pool index; member/width allocation unassigned for this study.')
def random_parts():
 return [Segment('inherited2018 RP-H handle',RP_HANDLE),role('inline cell barcode; member unassigned',molecular_feature=CELL_FEATURE),Segment('inherited2018 H9 tract','H'*9,placeholder=True)]
def contact_product():
 return Construct([role('contact locus A'),Segment('unfilled cohesive contact junction',DPNII.site),role('contact locus B')],name='unconverted DpnII contact-product example')
FINAL_LIBRARY=Construct([role('native P5/indexed end; geometry/member unassigned',molecular_feature=I5_FEATURE),*random_parts()[1:],role('contact-derived DNA copy; converted-parent signal encoded in bases'),role('native Adaptase-dependent span and P7/index end; geometry/member unassigned',molecular_feature=I7_FEATURE)],name='sn-m3C native source-limited library')
SEQ_PRIMERS=()
FINAL_CAPTION='Native contact-library endpoint roles; complete opposite-end geometry remains unassigned.'
SEQUENCING_INTRO='PE150: MiSeq/HiSeq4000 for cultured-cell samples and HiSeq4000/NovaSeq6000 for human cortex. Both reads are computationally trimmed25bases at5′ and3bases at3′; this is not a molecular tail-length claim.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1',150,'Contact/methylation sequence; source processing trims first25 and last3 bases','extract_then_include'),sp.read_spec('Read 2',150,'Contact/methylation sequence; source processing trims first25 and last3 bases','include'),sp.read_spec('Index 1 (i7)','unassigned','Native pool-index acquisition unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Native pool-index acquisition unassigned','undetermined'))
def final_library():return FINAL_LIBRARY

def role_scene(text,label):
 sc=Scene();sc.strand('DNA',[role(text)],label=label);return sc

def copied_state(*,adapted=False):
 sc=Scene();sc.strand('original',[role('original converted contact DNA; physical U at converted C; protected methyl-C retained')],label='independent original template; later fate unassigned')
 parts=[*random_parts(),role('new contact-DNA copy; methylation represented by copied base sequence, not retained chemical5mC')]
 if adapted:parts.append(role('native Adaptase-added3′ end; extent/sequence unassigned'))
 sc.strand('copy',parts,label='single-strand Adaptase/PCR precursor' if adapted else 'new copy; original/copy duplex occupancy unassigned before denaturation',mod5='SpC3')
 return sc

def post_pcr_state():
 sc=Scene.duplex(list(FINAL_LIBRARY),label='expected amplified library population',unpaired=tuple(p.name+"'" for p in FINAL_LIBRARY if p.is_role_token()))
 sc.strand('residual0',[role('possible co-carried original converted contact DNA; U/C chemistry/integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 sc.strand('residual1',[role('possible co-carried5′C3 RP-H/Adaptase contact-copy templates; integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 return sc

POST_PCR_STATE=post_pcr_state()

def cut_scene():
 sc=Scene()
 left=[role('left chromosomal flank; sequence/length unassigned')]
 left_bottom=[Segment('left fragment5′ cohesive end',revcomp(DPNII.overhang)),*complement_segments(left)]
 sc.strand('left top',left,label='left DpnII fragment',mod3='OH')
 sc.anneal('left bottom',left_bottom,to='left top',pair=(left[0].name+"'",left[0].name),unpaired=(left[0].name+"'",),label='left DpnII fragment',mod5='p')
 sc.overhang('left bottom','left fragment5′ cohesive end',polarity='5-prime',owner='substrate',text='DpnII5′ GATC')
 right=[Segment('right fragment5′ cohesive end',DPNII.overhang),role('right chromosomal flank; sequence/length unassigned')]
 right_bottom=complement_segments(right[1:])
 sc.strand('right top',right,label='right DpnII fragment',mod5='p')
 sc.anneal('right bottom',right_bottom,to='right top',pair=(right[1].name+"'",right[1].name),unpaired=(right[1].name+"'",),label='right DpnII fragment',mod3='OH')
 sc.overhang('right top','right fragment5′ cohesive end',polarity='5-prime',owner='substrate',text='DpnII5′ GATC')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 sample=g.input('sample','formaldehyde-fixed chromatin sample',MaterialKind.BIOLOGICAL,state=role_scene('chromosomal DNA in fixed chromatin','fixed chromatin'),primary=True)
 dpnii=g.input('dpnii','DpnII',MaterialKind.CATALYST)
 cut=g.product('cut','DpnII-cut chromatin',MaterialKind.MOLECULAR,state=cut_scene())
 g.react('digest','DpnII digestion37°C overnight',StepKind.DIGESTION,((sample,InputRole.TEMPLATE),(dpnii,InputRole.CATALYST)),(cut,))
 heat=g.input('heat','62°C10min DpnII inactivation',MaterialKind.CONDITION)
 inactive=g.product('inactive','enzyme-inactivated DpnII-cut chromatin',MaterialKind.MOLECULAR,state=cut.state,chemical_state=cut.state)
 g.react('inactivate','Heat-inactivate DpnII; no biotin fill-in',StepKind.THERMAL_TREATMENT,((cut,InputRole.TEMPLATE),(heat,InputRole.CONDITION)),(inactive,))
 ligase=g.input('ligase','T4 DNA ligase',MaterialKind.CATALYST)
 contact=contact_product();contact_scene=Scene.duplex(list(contact),label='unconverted expected cohesive contact',unpaired=tuple(s.name+"'" for s in contact if s.is_role_token()))
 contacts=g.product('contacts','unselected in situ contact products',MaterialKind.MOLECULAR,state=contact_scene)
 g.react('ligate_contacts','Contact ligation4h room temperature; no biotin pull-down',StepKind.PROXIMITY_LIGATION,((inactive,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(contacts,))
 qc_amount=g.input('qc_amount','10% sample for conventional-library shallow MiSeq QC',MaterialKind.CONDITION)
 remainder=g.product('remaining_contacts','contact-bearing sample for nucleus sorting',MaterialKind.MOLECULAR,state=contacts.state)
 qc=g.product('qc_fraction','10% conventional-library QC branch',MaterialKind.QC)
 g.react('take_qc','Take source10% conventional-library QC aliquot',StepKind.SAMPLING,((contacts,InputRole.TEMPLATE),(qc_amount,InputRole.CONDITION)),(remainder,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,reason='Conventional shallow MiSeq library QC is outside the selected bisulfite branch.'))
 gates=g.input('sorting','Hoechst/singlet/DNA-content FANS gates;384-well deposition',MaterialKind.CONDITION)
 nucleus=g.product('nucleus','single contact-bearing nucleus',MaterialKind.BIOLOGICAL,state=contacts.state)
 others=g.product('other_nuclei','nuclei outside selected single-nucleus lineage',MaterialKind.BIOLOGICAL)
 g.react('sort_nucleus','Sort individual intact Hoechst-labelled nuclei',StepKind.SAMPLING,((remainder,InputRole.TEMPLATE),(gates,InputRole.CONDITION)),(nucleus,others));g.mark_terminal(Terminal(others.id,Disposition.UNRESOLVED,reason='Other nuclei are outside this selected lineage; no physical discard inference.'))
 protease=g.input('protease','sorted-nucleus proteinase treatment; detailed input conditions unassigned',MaterialKind.CATALYST)
 released=g.product('released','released contact DNA',MaterialKind.MOLECULAR,state=role_scene('unconverted contact DNA; native fragments and junctions','released DNA'))
 g.react('release_dna','Digest sorted-nucleus proteins before inherited bisulfite workflow',StepKind.DEPROTEINIZATION,((nucleus,InputRole.TEMPLATE),(protease,InputRole.CATALYST)),(released,))
 bs=g.input('bisulfite','inherited2018 bisulfite conversion reagents',MaterialKind.CONDITION)
 converted=g.product('converted','physical converted contact DNA strands',MaterialKind.MOLECULAR,state=role_scene('converted contact DNA: original unmethylated C→U; protected methyl-C retained; contact junction C/U unassigned','original converted strands'))
 g.react('convert','Bisulfite convert contact DNA before library copying',StepKind.CHEMICAL_CONVERSION,((released,InputRole.TEMPLATE),(bs,InputRole.CONDITION)),(converted,))
 binding=g.input('binding','inherited2018 binding/desulfonation/wash recovery',MaterialKind.SELECTION)
 recovered=g.product('recovered','recovered converted contact DNA',MaterialKind.MOLECULAR,state=converted.state)
 bs_waste=g.product('bs_waste','conversion-recovery non-retained fractions',MaterialKind.WASTE)
 g.react('recover_bs','Bind/desulfonate/wash and recover converted contact DNA',StepKind.CLEANUP,((converted,InputRole.TEMPLATE),(binding,InputRole.SELECTION)),(recovered,bs_waste));g.discard(bs_waste,'Source-method physical recovery fraction; composition unassigned.')
 primer=g.input('random_primer','inherited2018 blocked RP-H indexed random primer; native member unassigned',MaterialKind.OLIGO)
 denature=g.input('denature','95°C3min then ice2min',MaterialKind.CONDITION)
 annealed=g.product('annealed','denatured converted contact templates with RP-H',MaterialKind.MOLECULAR,state=recovered.state)
 g.react('denature_originals','Inherited2018 converted-DNA denaturation',StepKind.THERMAL_TREATMENT,((recovered,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(denature,InputRole.CONDITION)),(annealed,))
 klenow=g.input('klenow','inherited2018 Klenow exo-minus random extension',MaterialKind.CATALYST)
 copied=g.product('copied','indexed contact-DNA copies alongside original templates',MaterialKind.MOLECULAR,state=copied_state())
 g.react('random_prime','Inherited random extension4°C5min,25°C5min,37°C60min',StepKind.EXTENSION,((annealed,InputRole.TEMPLATE),(klenow,InputRole.CATALYST)),(copied,))
 exo=g.input('exo','ExonucleaseI',MaterialKind.CATALYST);rsap=g.input('rsap','rSAP',MaterialKind.CATALYST)
 treated=g.product('treated','treated expected contact-DNA copies',MaterialKind.MOLECULAR,state=copied.state)
 chemical=g.product('chemical_products','ExoI/rSAP chemical products; physical partition unassigned',MaterialKind.MOLECULAR)
 g.react('exo_rsap','Inherited ExoI/rSAP37°C30min',StepKind.DIGESTION,((copied,InputRole.TEMPLATE),(exo,InputRole.CATALYST),(rsap,InputRole.CATALYST)),(treated,chemical));g.mark_terminal(Terminal(chemical.id,Disposition.UNRESOLVED,reason='Hydrolysis/inactivation is distinct from physical cleanup.'))
 other_reactions=g.input('other_reactions','seven independently inline-indexed reactions for inherited eight-way pooling',MaterialKind.MOLECULAR)
 pooled=g.product('pooled','pooled contact libraries; selected nucleus lineage shown',MaterialKind.MOLECULAR,state=treated.state)
 g.react('pool_eight','Inherited2018 eight-way inline-indexed reaction pooling',StepKind.POOL,((treated,InputRole.TEMPLATE),(other_reactions,InputRole.TEMPLATE)),(pooled,))
 spri=g.input('spri','SPRI beads',MaterialKind.SELECTION)
 clean=g.product('clean','recovered native copies in10µlEB',MaterialKind.MOLECULAR,state=pooled.state);waste=g.product('priming_waste','0.8×SPRI non-retained fraction',MaterialKind.WASTE)
 g.react('clean_primed','Inherited0.8×SPRI before Adaptase',StepKind.CLEANUP,((pooled,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(clean,waste));g.discard(waste,'Physical recovery fraction; molecular depletion unassigned.')
 adapt_heat=g.input('adapt_heat','95°C3min then ice2min',MaterialKind.CONDITION)
 ss=g.product('ss_copies','denatured contact-DNA copies for Adaptase',MaterialKind.MOLECULAR,state=clean.state)
 g.react('denature_copies','Denature recovered copies before Adaptase',StepKind.THERMAL_TREATMENT,((clean,InputRole.TEMPLATE),(adapt_heat,InputRole.CONDITION)),(ss,))
 adapter=g.input('adapter','native Adaptase3′ completion end; bases unassigned',MaterialKind.ADAPTER);adaptase=g.input('adaptase','Adaptase module',MaterialKind.CATALYST)
 adapted=g.product('adapted','native single-strand Adaptase-completed copies',MaterialKind.MOLECULAR,state=copied_state(adapted=True))
 g.react('adaptase_tag','Inherited Adaptase37°C30min,95°C2min',StepKind.ADAPTER_ATTACHMENT,((ss,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(adaptase,InputRole.CATALYST)),(adapted,))
 f=g.input('forward','native2018-family P5 indexed PCR primer; member unassigned',MaterialKind.OLIGO);r=g.input('reverse','native2018-family P7 indexed PCR primer; member unassigned',MaterialKind.OLIGO);kapa=g.input('kapa','KAPA HiFi HotStart library PCR',MaterialKind.CATALYST)
 amp=g.product('amplified','indexed contact-library pool',MaterialKind.LIBRARY,state=POST_PCR_STATE)
 g.react('library_pcr','Inherited2018 native library PCR15cycles',StepKind.PCR,((adapted,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(kapa,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(adapted.state,reason='Actual single-strand native Adaptase contact copies; selected primers and native adapter-site geometry unassigned.'))
 clean1=g.product('clean1','first recovered pool',MaterialKind.LIBRARY,state=POST_PCR_STATE);w1=g.product('waste1','first final cleanup non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup1','Inherited final0.8×SPRI cleanup1',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(clean1,w1));g.discard(w1,'Physical cleanup fraction.')
 clean2=g.product('clean2','second recovered pool',MaterialKind.LIBRARY,state=POST_PCR_STATE);w2=g.product('waste2','second final cleanup non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup2','Inherited final0.8×SPRI cleanup2',StepKind.CLEANUP,((clean1,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(clean2,w2));g.discard(w2,'Physical cleanup fraction.')
 library=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=POST_PCR_STATE);w3=g.product('waste3','third final cleanup non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup3','Inherited final0.8×SPRI cleanup3',StepKind.CLEANUP,((clean2,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(library,w3));g.discard(w3,'Physical cleanup fraction.')
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),read_lengths={'Read 1':150,'Read 2':150},limitation='PE150 source-stated; exact native primer landing and index acquisition unassigned. Computational25/3 trimming does not establish chemical adapter lengths.'))
 g.declare_reviewed_scope(('sample','dpnii','digest','cut','heat','inactivate','inactive','ligase','ligate_contacts','contacts','qc_amount','take_qc','remaining_contacts','sorting','sort_nucleus','nucleus','protease','release_dna','released','bisulfite','convert','converted','binding','recover_bs','recovered','random_primer','denature','denature_originals','annealed','klenow','random_prime','copied','exo','rsap','exo_rsap','treated','other_reactions','pool_eight','pooled','spri','clean_primed','clean','adapt_heat','denature_copies','ss_copies','adapter','adaptase','adaptase_tag','adapted','forward','reverse','kapa','library_pcr','amplified','cleanup1','clean1','cleanup2','clean2','cleanup3','library'),SOURCE)
 return g.require_publishable()
