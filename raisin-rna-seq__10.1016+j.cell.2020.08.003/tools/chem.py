"""RAISIN2020 plate route: native particle/RNA populations and actual PCR inputs."""
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp
from protocol_model import role,molecule
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition
from types import SimpleNamespace
import seqprimers as sp
RAISIN_CST=(("CHAPS","0.49% w/v"),("NaCl","146 mM"),("CaCl2","1 mM"),("MgCl2","21 mM"),("Tris pH 8.0","10 mM"))
RAISIN_RT_PRIMER='AAGCAGTGGTATCAACGCAGAGTACT(30)VN'
RAISIN_RT_HANDLE=RAISIN_RT_PRIMER.removesuffix('T(30)VN')
RAISIN_LOCKED_TSO='AAGCAGTGGTATCAACGCAGAGTACrGrG+G'
TSO_HANDLE=RAISIN_LOCKED_TSO.removesuffix('rGrG+G')
RAISIN_PCR_PRIMER='AAGCAGTGGTATCAACGCAGAGT'
if not all(x.startswith(RAISIN_PCR_PRIMER) for x in (RAISIN_RT_HANDLE,TSO_HANDLE)):raise ValueError('Source SMART amplification primer must match both native handles')
TITLE='RAISIN RNA-seq — intact nuclei with ribosome-bound RNA'
NOTES='01_raisin_rna_seq.html'
SOURCE='Drokhlyansky et al., Cell2020, <a href="https://doi.org/10.1016/j.cell.2020.08.003">doi:10.1016/j.cell.2020.08.003</a>; PMC8358727 Methods1032–1076/KeyResources3032–3040; expressly inherited Shekhar2016 plateMethods290.'
SUMMARY='CST preserves intact nuclei with associated ribosomes/RNA; selected plate SMART-seq2 route uses source oligo-dT/lockedTSO,21cycle preamplification and NexteraXT with custom barcode adapters.'
CAVEAT='The separate GemCode droplet RAISIN implementation is not merged. Particle components and nuclear/ribosome-associated RNA are independent populations. Physical custom XT barcode/adapter species and native tagmentation/PCR geometry remain unassigned. The named75-cycle NextSeq kit does not assign acquired read/index allocation; no generic8nt indices or fixed transcript lengths are inserted.'
SEQ_PRIMERS=()
READ_LENGTHS=None
SEQUENCING_INTRO='Source plate libraries pool384 uniquely barcoded wells and use NextSeq500 FC-404-2005 (75-cycle high-outputv2 kit); exact biological/index acquisition remains unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE='Actual custom XT adapters and native sequencing-primer geometry are not reconstructed.'
READ_SPECS=tuple(sp.read_spec(x,None,'Acquisition allocation unassigned for source75-cycle kit; no predecessor run mask substituted','undetermined') for x in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))
FINAL_LIBRARY=Construct([role('custom XT left end; native index placement/width unassigned',molecular_feature=feature('raisin_index_left','sample_index','unknown')),role('tagmented transcript-derived fragment; extent unassigned'),role('custom XT right end; native index placement/width unassigned',molecular_feature=feature('raisin_index_right','sample_index','unknown'))],name='RAISIN source-limited custom XT endpoint')
FINAL_CAPTION='Expected source-kit endpoint; full-length cDNA and native tagged substrate are distinct earlier molecules.'
def rt_parts():return [Segment('SMART RT handle',RAISIN_RT_HANDLE),Segment('oligo(dT)30','T'*30),Segment('anchoredVN','VN',placeholder=True,length_bp=2)]
def cdna_parts(*,switched=True):
 p=[*rt_parts(),role('RNA-derived first-strand cDNA interval')]
 if switched:p += [Segment('RT terminal C annealing tract','CCC',note='model-derived representative3C/3G template-switch relation; exact native junction heterogeneity unassigned',inferred=True),Segment('copied locked TSO handle',revcomp(TSO_HANDLE),inferred=True)]
 return p

def cdna_scene(*,switched=True):
 sc=Scene();sc.strand('cDNA',cdna_parts(switched=switched),label='INFERRED — representative switched ssDNA' if switched else 'actual first-strand target',mod3='OH')
 if switched:sc.note('cDNA','Source locked3rG/rG/+G template is not covalently attached; copied handle is a derived representative, native terminal-C variation unassigned.')
 return sc

def particle_scene():
 sc=Scene()
 for k,label in (('nucleus','intact nucleus'),('ribosomes','nucleus-associated roughER/ribosome population'),('RNA','associated poly(A) RNA population')):sc.strand(k,[role(label)],label=label)
 return sc

def rna_scene():
 sc=Scene();sc.strand('nuclearRNA',[role('nuclear RNA templates')],label='independent nuclear RNA pool');sc.strand('exteriorRNA',[role('ribosome-associated poly(A) RNA templates')],label='independent associated RNA pool');return sc

def tagged_scene():
 sc=Scene()
 for side in ('A','B'):sc.strand('tag'+side,[role('source NexteraXT transferred arm'+side+'; physical strand/gap unassigned'),role('tagmented cDNA-derived strand interval'+side)],label='independent tagged-strand class; native pairing unassigned')
 return sc

def oligos():return (oligo('Source SMART oligo-dT RT primer',rt_parts()),'5′-'+RAISIN_LOCKED_TSO+'-3′ — source locked SMART TSO; rG/rG/+G as printed',oligo('Source single SMART cDNA PCR species',[Segment('SMART PCR primer',RAISIN_PCR_PRIMER)]))

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 tissue=g.input('tissue','frozen tissue',MaterialKind.BIOLOGICAL,state=molecule('frozen tissue',(role('frozen tissue sample'),)),primary=True)
 cst=g.input('cst','CST extraction: '+', '.join(k+' '+v for k,v in RAISIN_CST),MaterialKind.CONDITION)
 particles=g.product('particles','RAISIN particles: intact nuclei plus associated ribosomes/RNA',MaterialKind.MOLECULAR,state=particle_scene())
 g.react('extract','Extract RAISIN particles using source CST',StepKind.EXTRACTION,((tissue,InputRole.TEMPLATE),(cst,InputRole.CONDITION)),(particles,))
 rule=g.input('sort_rule','Ruby-stained single RAISIN into each well;96well plate',MaterialKind.SELECTION)
 selected=g.product('single_raisin','one sorted RAISIN per well',MaterialKind.MOLECULAR,state=particles.state)
 other=g.product('other_particles','particles outside selected single-RAISIN well',MaterialKind.MOLECULAR)
 g.react('sort','Sort one RAISIN per source well',StepKind.SAMPLING,((particles,InputRole.TEMPLATE),(rule,InputRole.SELECTION)),(selected,other),changes_molecule=False);g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Not selected into this well; collection/disposal of other particles is unassigned.'))
 lysis=g.input('lysis','5µLTCL+1%βME;centrifuge800g3min,freeze−80°C',MaterialKind.CONDITION)
 lysate=g.product('lysate','lysed nuclear and associated RNA populations',MaterialKind.MOLECULAR,state=rna_scene())
 g.react('lyse','Lyse the sorted particle before RNA purification',StepKind.EXTRACTION,((selected,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(lysate,))
 beads=g.input('rnaclean','RNACleanXP; inherited Shekhar plate route2.2×,no final elution',MaterialKind.SELECTION)
 rna=g.product('purified_rna','bead-purified RNA; intact particle no longer the RT substrate',MaterialKind.MOLECULAR,state=rna_scene())
 washes=g.product('rnaclean_washes','excluded lysate/supernatant and bead washes',MaterialKind.WASTE)
 g.react('rna_cleanup','Purify RNA from lysed RAISIN wells',StepKind.CLEANUP,((lysate,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,washes));g.discard(washes,'Source bead/supernatant removal, not retention of intact nuclei as RT substrate.')
 primer=g.input('rt_primer','source SMART anchored oligo-dT',MaterialKind.OLIGO,state=Construct(rt_parts(),name='source RT primer'))
 rt=g.input('reverse_transcriptase','MaximaHMinus EP0753; selected source replacement',MaterialKind.CATALYST)
 first=g.product('first_strand','first-strand RAISIN cDNA target',MaterialKind.MOLECULAR,state=cdna_scene(switched=False))
 residual=g.product('rt_rna','RNA templates after copying',MaterialKind.MOLECULAR,state=rna.state)
 g.react('reverse_transcribe','Oligo-dT reverse-transcribe purified RNA',StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Reverse transcription does not consume RNA. Template recovery/removal after this target-copying boundary is not reconstructed.'))
 tso=g.input('locked_tso','source locked SMART TSO; rGrG+G physical3′end',MaterialKind.OLIGO,state=RAISIN_LOCKED_TSO)
 switched=g.product('switched','representative template-switched ssDNA target',MaterialKind.MOLECULAR,state=cdna_scene(),inferred=True)
 remaining_tso=g.product('remaining_tso','TSO molecules after copying',MaterialKind.MOLECULAR,state=RAISIN_LOCKED_TSO)
 g.react('template_switch','INFERRED — template-switch copy within same RT incubation',StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,remaining_tso));g.mark_terminal(Terminal(remaining_tso.id,Disposition.UNRESOLVED,'TSO is a separate template reagent; copying does not join the oligo covalently or establish physical disposal.'))
 f=g.input('preamp_forward','source SMART cDNA PCR primer species',MaterialKind.OLIGO,state=RAISIN_PCR_PRIMER)
 r=g.input('preamp_reverse','same source SMART PCR species at opposing handle',MaterialKind.OLIGO,state=RAISIN_PCR_PRIMER)
 kapa=g.input('preamp_polymerase','KAPAHiFiHotStartReadyMix;21cycles',MaterialKind.CATALYST)
 con=Construct(cdna_parts(),name='representative amplified full-length SMART cDNA')
 ampstate=Scene.duplex(list(con),unpaired=tuple(s.name+"'" for s in con if s.is_role_token()))
 amplified=g.product('amplified','21cycle amplified duplex cDNA',MaterialKind.MOLECULAR,state=ampstate,inferred=True)
 g.react('preamp','INFERRED representative product — amplify actual switched ssDNA21cycles',StepKind.PCR,((switched,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(kapa,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(switched.state,symmetric_species=True,reason='Actual single-strand switched cDNA precedes first-copy phase. Source single primer and derived representative handles do not establish exact native switch junction/copy geometry.'))
 ampure=g.input('ampure','source two consecutiveAMPureXP cDNA purifications',MaterialKind.SELECTION)
 clean=g.product('purified_cdna','twice purified amplified duplex cDNA',MaterialKind.MOLECULAR,state=amplified.state,inferred=True)
 cw=g.product('cdna_wash','AMPure supernatant and wash fractions',MaterialKind.WASTE)
 g.react('cdna_cleanup','Purify amplified cDNA twice',StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(clean,cw));g.discard(cw,'Physically excluded source bead supernatant/washes.')
 xt=g.input('transposome','NexteraXT FC-131-1096; inherited0.075ngcDNA/¼reactions;custombarcodeadapters',MaterialKind.ADAPTER)
 tn5=g.input('tn5','source NexteraXT transposase',MaterialKind.CATALYST)
 tagged=g.product('tagged','native tagged cDNA substrate before indexPCR',MaterialKind.MOLECULAR,state=tagged_scene())
 nonamp=g.product('nonamplifying','candidate same-end/non-library tagged molecules',MaterialKind.MOLECULAR)
 g.react('tagment','Tagment amplified cDNA using inherited plate route',StepKind.TAGMENTATION,((clean,InputRole.TEMPLATE),(xt,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,nonamp));g.mark_terminal(Terminal(nonamp.id,Disposition.UNRESOLVED,'Failure to form an amplifiable endpoint does not physically remove a molecule; recovery is unassigned.'))
 lf=g.input('library_forward','custom XT index-adapter primer; physical sequence/width unassigned',MaterialKind.OLIGO)
 lr=g.input('library_reverse','opposing custom XT index-adapter primer;physicalsequence/widthunassigned',MaterialKind.OLIGO)
 pol=g.input('library_polymerase','NexteraXT PCR system;¼source predecessor reaction',MaterialKind.CATALYST)
 lib=g.product('library','RAISIN indexed XT library',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 g.react('library_pcr','Index-amplify actual tagged cDNA targets',StepKind.PCR,((tagged,InputRole.TEMPLATE),(lf,InputRole.STRUCTURAL),(lr,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(tagged.state,reason='Actual native tagged strands, not final indexed duplex. Custom XT primer members/widths and native gap/first-copy geometry unassigned.'))
 g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_ENDING+' '+SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(('tissue','cst','extract','sort_rule','sort','lysis','lyse','rnaclean','rna_cleanup','rt_primer','reverse_transcriptase','reverse_transcribe','locked_tso','template_switch','preamp_forward','preamp_reverse','preamp_polymerase','preamp','ampure','cdna_cleanup','transposome','tn5','tagment','library_forward','library_reverse','library_polymerase','library_pcr','library'),SOURCE)
 return g.require_publishable()
SPEC=SimpleNamespace(**{k:v for k,v in globals().copy().items() if k.isupper()},oligos=oligos,protocol_graph=protocol_graph)
