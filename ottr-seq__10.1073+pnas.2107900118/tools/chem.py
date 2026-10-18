"""Published OTTR universal-adapter branch; native ordered relay and ss precursor."""
from chemdraw import Construct,Scene,Segment,feature
from protocol_model import molecule,role
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition
import seqprimers as sp
TITLE='OTTR-seq — Ordered Two-Template Relay'
NOTES='01_ottr_seq.html'
SOURCE='Upton et al., PNAS(2021), <a href="https://doi.org/10.1073/pnas.2107900118">doi:10.1073/pnas.2107900118</a>; published PMC8594491 Fig3/Methods186,194.'
SUMMARY='BoMoC adds one nonextendable ddR to RNA, then orders two template jumps: a +1Y DNA primer copies RNA, and its added dG recruits the terminating3′C adapter template.'
CAVEAT='Selected universal-adapter/Q5 branch. Full-length PCR-free adapters are a separate source alternative. Published SI TableS1 is still unverified; superseded preprint oligos are excluded. Source-known overhangs/blocks are retained, while complete primer/AT sequence and block identities remain unassigned. The relay phases below occur within one RT incubation; no intervening separation is implied.'
SEQ_PRIMERS=()
READ_LENGTHS=None
SEQUENCING_INTRO='MiniSeq75-cycle high-output kit; read/index allocation is unassigned, not inferred from the kit name.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE='Native sequencing-primer sites and acquired read/index mask have not been reconstructed.'
READ_SPECS=tuple(sp.read_spec(name,None,'Acquisition allocation unassigned; source specifies75-cycle kit only','undetermined') for name in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))

def copy_parts(*,relay=False):
 p=[role('universal DNA primer; complement of R2'),Segment('primer terminalY','Y',placeholder=True,length_bp=1,note='source +1T/+1C mixture; not identifier'),role('input-RNA-derived cDNA excluding matched added ddR')]
 if relay:p += [Segment('single non-templated dG relay base','G'),role('copy of terminating adapter beyond matched3′C; complement of R1')]
 return p

def copy_scene(*,relay=False):
 sc=Scene();sc.strand('cDNA',copy_parts(relay=relay),label='actual ssDNA relay product' if relay else 'first-jump ssDNA intermediate',mod3='OH');return sc

def capped_rna():
 sc=Scene()
 for b in ('A','G'):
  sc.strand('dd'+b,[role('input RNA'),Segment('single dd'+b,b)],label='dd'+b+'-capped RNA population',mod3='dideoxy; no3′OH')
 return sc

def starter_scene():
 sc=Scene();sc.strand('DNA',[role('universal DNA primer; complement of R2'),Segment('DNA +1Y','Y',placeholder=True,length_bp=1)],label='source DNA primer',mod3='OH')
 sc.strand('RNA',[role('RNA partner of starter duplex; exact sequence unassigned')],label='independent primer-duplex RNA strand',mod3='blocked; identity unassigned')
 sc.note('DNA','Source describes RNA/DNA duplex with1nt DNA3′Y overhang. Exact paired core and blocking chemistry await published TableS1; no fictitious pairing columns are assigned.')
 return sc

FINAL_LIBRARY=Construct([role('primer-side P7-derived end; sequence/orientation unassigned'),role('i7; member/width unassigned',molecular_feature=feature('ottr_i7','sample_index','unknown')),*copy_parts(relay=True),role('i5 complement; member/width unassigned',molecular_feature=feature('ottr_i5','sample_index','unknown')),role('AT-side P5-derived end; sequence/orientation unassigned')],name='OTTR universal-adapter source-limited PCR endpoint')
FINAL_CAPTION='Primer-oriented endpoint; actual PCR input is ss relay cDNA, not this expected amplified product.'
def oligos():return ()

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 rna=g.input('rna','10ng source small RNA; miRXplore, cellular orEV pools',MaterialKind.BIOLOGICAL,state=molecule('native RNA',(role('input RNA'),)),primary=True)
 bomoc=g.input('tailing_bomoc','BoMoC0.7µM,2mMMnCl2,ddATP250µM ±ddGTP;30°C1.5–2h',MaterialKind.CATALYST)
 capped=g.product('capped_rna','single-ddR capped RNA populations',MaterialKind.MOLECULAR,state=capped_rna())
 g.react('cap_end','Add one blocking ddA/ddG to each input3′end',StepKind.TAILING,((rna,InputRole.TEMPLATE),(bomoc,InputRole.CATALYST)),(capped,))
 phosph=g.input('sap','65°C5min then MgCl2/SAP0.5U;37°C15min,EGTA,65°C5min',MaterialKind.CATALYST)
 toggled=g.product('ddntp_depleted','capped RNA after free-ddNTP dephosphorylation',MaterialKind.MOLECULAR,state=capped.state)
 free=g.product('dephosphorylated_ddntp','dephosphorylated free ddNTP products',MaterialKind.MOLECULAR)
 g.react('deplete_ddntp','Inactivate free ddNTPs and toggle Mn/Mg conditions',StepKind.DEPHOSPHORYLATION,((capped,InputRole.TEMPLATE),(phosph,InputRole.CATALYST)),(toggled,free));g.mark_terminal(Terminal(free.id,Disposition.UNRESOLVED,'Chemical conversion does not physically remove nucleotide products.'))
 starter=g.input('starter','90nM RNA/DNA primer duplex,+1T/+1C',MaterialKind.OLIGO,state=starter_scene())
 rt=g.input('relay_bomoc','BoMoC0.5µM;200µMdGTP,40µMdTTP/dCTP,2µMdATP,150µMDAP;Mg/PEG conditions',MaterialKind.CATALYST)
 first=g.product('first_copy','first-jump cDNA intermediate',MaterialKind.MOLECULAR,state=copy_scene())
 residual=Scene();residual.strand('inputRNA',[role('ddR-capped input RNA templates')],label='original RNA remains',mod3='ddR');residual.strand('starterRNA',[role('blocked RNA starter partner')],label='independent starterRNA',mod3='blocked')
 rnas=g.product('rt_rnas','RNA templates/primer partner after copying',MaterialKind.MOLECULAR,state=residual)
 g.react('first_jump','RT phase1: +1Y engages ddR and copies inputRNA end to end',StepKind.TEMPLATE_SWITCHING,((toggled,InputRole.TEMPLATE),(starter,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,rnas))
 ats=Scene();ats.strand('AT',[role('universal R1 adapter template'),Segment('AT terminalC','C')],label='terminating adapter template',mod5='jump-blocking group; identity unassigned',mod3='blocked; identity unassigned')
 at=g.input('terminating_at','180nM terminating3′C adapter template',MaterialKind.OLIGO,state=ats)
 dual=g.product('dual_adapter_cdna','ssDNA after ordered second template jump',MaterialKind.MOLECULAR,state=copy_scene(relay=True))
 at_residual=g.product('at_residual','adapter templates after copying',MaterialKind.MOLECULAR,state=ats)
 g.react('second_jump','RT phase2: add dG, engage AT3′C and copy to its5′block',StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(at,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(dual,at_residual));g.mark_terminal(Terminal(at_residual.id,Disposition.UNRESOLVED,'Template copying does not consume AT; nucleic-acid type/block-specific recovery is not reconstructed.'))
 rnase=g.input('rnase','Methods186/194 processed-as-above RNaseA0.5µg/µL and thermostableRNaseH0.5U',MaterialKind.CATALYST)
 digested=g.product('digested_cdna','ss relay cDNA target after RNA digestion',MaterialKind.MOLECULAR,state=dual.state)
 debris=g.product('rna_debris','RNA hydrolysis products',MaterialKind.MOLECULAR)
 g.react('digest_rna','Digest input/starterRNA; retain the cDNA target',StepKind.DIGESTION,((dual,InputRole.TEMPLATE),(rnas,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(digested,debris));g.mark_terminal(Terminal(debris.id,Disposition.UNRESOLVED,'RNA fragments are molecular products; their partition through subsequent extraction/precipitation remains unassigned.'))
 extraction=g.input('rt_recovery','PCI extraction/ethanol precipitation;resuspend10µL',MaterialKind.SELECTION)
 clean=g.product('purified_cdna','recovered ssDNA relay targets',MaterialKind.MOLECULAR,state=dual.state)
 waste=g.product('rt_recovery_waste','physically separated extraction/precipitation fractions',MaterialKind.WASTE)
 g.react('recover_cdna','Recover RT products before universal-adapter PCR',StepKind.CLEANUP,((digested,InputRole.TEMPLATE),(extraction,InputRole.SELECTION)),(clean,waste));g.discard(waste,'Source PCI phase separation and precipitation-supernatant removal; detailed partition unassigned.')
 aliquot=g.product('pcr_input','5µL of10µL purified RT product',MaterialKind.MOLECULAR,state=clean.state)
 rest=g.product('remaining_cdna','RT product not used in this PCR',MaterialKind.MOLECULAR,state=clean.state)
 sample=g.input('aliquot_rule','5µL of recovered10µL RT product',MaterialKind.CONDITION)
 g.react('aliquot','Use half of the recovered RT product for PCR',StepKind.SAMPLING,((clean,InputRole.TEMPLATE),(sample,InputRole.CONDITION)),(aliquot,rest),changes_molecule=False);g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Not used in the source PCR aliquot; physical disposal is unassigned.'))
 f=g.input('pcr_forward','P7/i7 indexing primer; physical species unassigned',MaterialKind.OLIGO)
 r=g.input('pcr_reverse','P5/i5 indexing primer; physical species unassigned',MaterialKind.OLIGO)
 q5=g.input('q5','Q5 high-fidelity polymerase;4–8cycles',MaterialKind.CATALYST)
 amp=g.product('amplified_library','universal-adapter indexed PCR product',MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 g.react('library_pcr','Amplify actual ss universal-adapter relay cDNA',StepKind.PCR,((aliquot,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(aliquot.state,reason='Actual ss relay cDNA enters Q5; complete source primer/AT sequences and first-copy geometry await published TableS1 verification.'))
 select=g.input('final_recovery','Column cleanup,7.5%ureaPAGE;37°C overnight diffusion,PCI/ethanol',MaterialKind.SELECTION)
 lib=g.product('library','purified universal-adapter OTTR library',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 off=g.product('final_waste','column/PAGE excluded fractions including adapter dimer',MaterialKind.WASTE)
 g.react('purify_library','Remove adapter dimers and recover the amplified library',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(select,InputRole.SELECTION)),(lib,off));g.discard(off,'Physically excluded column/gel and extraction fractions.')
 g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_ENDING+' '+SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(('rna','tailing_bomoc','cap_end','sap','deplete_ddntp','starter','relay_bomoc','first_jump','terminating_at','second_jump','rnase','digest_rna','rt_recovery','recover_cdna','aliquot_rule','aliquot','pcr_forward','pcr_reverse','q5','library_pcr','final_recovery','purify_library','library'),SOURCE)
 return g.require_publishable()
