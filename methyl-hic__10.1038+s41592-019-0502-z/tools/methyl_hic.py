"""Li2019 bulk Methyl-HiC native protocol; study adapters remain source-bounded."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,complement_segments,feature
from restriction import DPNII,Digest
from protocol_model import role
from endprep import EndPreparedDuplex,dA_tailed_scene
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as R,StepKind as K,Terminal,Disposition,SequencingPlan
import seqprimers as sp
TITLE='Methyl-HiC — bulk chromatin contacts and DNA methylation'
NOTES='01_methyl-hic.html'
SOURCE='<a href="https://doi.org/10.1038/s41592-019-0502-z">Li et al., Nature Methods2019</a>, Methods Methyl-HiC and published MOESM1 bulk protocol. Single-cell Adaptase and WGBS controls are distinct routes.'
SUMMARY='DpnII overhang fill-in, proximity ligation, recovery and sonication precede streptavidin capture. Native adapter ligation comes before physical bisulfite conversion and uracil-tolerant PCR.'
CAVEAT='Native indexed adapter strands, methylation/protection and PCR-primer sequences remain unassigned; no canonical TruSeq shell or8nt index is substituted. The main article names DpnII and published MOESM1 Biotin Fill-in step1 likewise says inactivateDpnII; the older MboI wording was not this published original and is not an active source conflict. Bulk conversion recovery remains a real source boundary: main Methods say purified convertedDNA, whereas MOESM1 pools column eluate and beads before PCR. The graph retains mixed original/copy/lambda classes with exact component recovery unassigned. Source conversion does not turn the physical pre-PCR U substrate directly into T, and PCR does not copy biotin/methylation modifications.'
DIGEST=Digest(DPNII,'X'*16,'X'*16)
FILLED=DIGEST.fill_in(biotin_base='A')
INDEX=feature('methyl_hic_index','sample_index','unknown',note='native adapter index member/count/width unassigned')
SEQ_PRIMERS=()
BULK_READ_NT=100
def contact_product():
 return Construct([role('locus A; native methylation/extent unassigned'),Segment('DpnII filled contact junction',FILLED.junction),role('locus B; native methylation/extent unassigned')],name='bulk DpnII contact reference')
def contact_scene():
 parts=list(contact_product());sc=Scene();sc.strand('top',parts,label='original genomic/copy fill strand; biotin-dA at junction positions'+str(FILLED.junction_biotin_top))
 bottom=complement_segments(parts);sc.anneal('bottom',bottom,to='top',pair=(bottom[0].name,parts[0].name),unpaired=tuple(p.name for p in bottom if p.is_role_token()),label='opposite original strand; drawn biotin-dA positions'+str(FILLED.junction_biotin_bottom));return sc

def restriction_scene():
 sc=Scene()
 for key in ('left','right'):
  body=role(key+' fragment body; native methylation retained')
  sc.strand(key,[Segment(key+'5′GATC','GATC'),body],label='DpnII cohesive end',mod5='p')
  opp=complement_segments([body]);sc.anneal(key+' opposite',opp,to=key,pair=(opp[0].name,body.name),unpaired=(opp[0].name,),label='recessed3′OH')
  sc.overhang(key,key+'5′GATC',polarity='5-prime',owner='substrate',text='DpnII5′GATC')
 return sc

def adapted_scene():
 sc=Scene();sc.strand('adapted contact',[role('native adapter end A; protection/geometry unassigned'),*contact_product(),role('native indexed adapter end B',molecular_feature=INDEX)],label='actual adapter-ligated original contact strand; opposing native assembly unassigned')
 sc.strand('other originals',[role('opposite original strands/offtarget adapted products; native state/partition unassigned')],label='independent source mixture');return sc

def converted_scene():
 sc=Scene()
 for key in ('contact','lambda control'):
  sc.strand(key,[role('native converted adapter end; protection/retention unassigned'),role(key+' strand; native order/extent unassigned'),Segment('representative unmethylated cytosine after bisulfite','U',note='physical deoxyuridine'),role('protected methylcytosine remains5mC; presence/site unassigned' if key=='contact' else 'unmethylated lambda remaining converted sequence'),role('opposite native adapter end',molecular_feature=INDEX)],label='actual single-stranded converted '+key+' precursor; U is not yet amplified T')
 sc.strand('conversion residuals',[role('possible original/incompletely converted/hydrolyzed or damaged molecules; recovery unassigned')],label='independent residual original/background population');return sc

def final_library():
 return Construct([role('native primer-derived end A'),role('contact sequence; original methylation readout C/T, native junction conversion classes'),role('native sample-index-bearing end B',molecular_feature=INDEX)],name='Methyl-HiC native amplified endpoint')
FINAL_LIBRARY=final_library()
FINAL_CAPTION='Native source-limited daughter library; C/T is amplified readout, not the physical U-bearing substrate. No generic index width assigned.'
SEQUENCING_INTRO='Main Methods: paired-end HiSeq2500/4000, at least50% balanced libraries/PhiX for converted preparations. Supplement bulk sonication says for PE2×100; source-specific index acquisition remains unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1','unassigned','paired-end; bulk preparation PE2×100 reference does not establish acquired cycles','include'),sp.read_spec('Read 2','unassigned','paired-end opposite converted contact arm; acquired cycles unassigned','include'),sp.read_spec('Index 1 (i7)','unassigned','Native sample-index acquisition unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Native index acquisition unassigned','undetermined'))
def daughter_scene():
 sc=Scene.duplex(list(FINAL_LIBRARY),label='normal daughter DNA; C/T readout',unpaired=tuple(p.name+"'" for p in FINAL_LIBRARY))
 sc.strand('lambda daughter',[role('normalDNA lambda control copies; same native adapter class')],label='independent control daughter population')
 sc.strand('residual templates',[role('possible U/methylated/biotinylated original templates and offproducts; physical recovery unresolved')],label='independent original lineage; modifications not copied');return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
 def out(i,l,state=None,k=M.MOLECULAR):return g.product(i,l,k,state=state)
 def op(i,l,kind,x,reagents,state,geometry=None):
  y=out(i+'_out',l,state);g.react(i,l,kind,((x,R.TEMPLATE),*reagents),(y,),geometry=geometry);return y
 def wash(i,label,x,state=None):
  select=inp(i+'_selection',label,M.SELECTION);y=out(i+'_out',label+' retained fraction',state or x.state,x.kind);w=out(i+'_waste','physical removed fraction',k=M.WASTE)
  g.react(i,label,K.CLEANUP,((x,R.TEMPLATE),(select,R.SELECTION)),(y,w));g.discard(w,'Physical source recovery fraction; class-specific purity unassigned.');return y
 sample=g.input('sample','fixed/permeabilized bulk nuclei; main1%formaldehyde10min and MOESM1 nuclei isolation/SDS',M.BIOLOGICAL,state=Scene.duplex([role('chromosomal DNA; native methylation retained')],label='fixed genomic DNA',unpaired=("chromosomal DNA; native methylation retained'",)),primary=True)
 enz=inp('dpnii','DpnII100U overnight37°C',M.CATALYST);x=op('digest','Digest fixed nuclei',K.DIGESTION,sample,((enz,R.CATALYST),),restriction_scene())
 heat=inp('heat','published MOESM1 DpnII inactivation65°C20min then cool');x=op('heat_inactivate','Heat-inactivate DpnII',K.THERMAL_TREATMENT,x,((heat,R.CONDITION),),x.state)
 fill=inp('klenow_fill','Klenow plus biotin-14-dATP/dCTP/dGTP/dTTP;37°C90min',M.CATALYST)
 x=op('fill','Fill5′GATC to blunt biotin ends',K.EXTENSION,x,((fill,R.CATALYST),),Scene.duplex([role('independent filled blunt restriction fragments; source GATC/biotin-dA extension')],label='filled fragment reference',unpaired=("independent filled blunt restriction fragments; source GATC/biotin-dA extension'",)))
 lig=inp('contact_ligase','T4 DNA ligase RT4h',M.CATALYST);x=op('contact_ligation','In-situ proximity ligation',K.PROXIMITY_LIGATION,x,((lig,R.CATALYST),),contact_scene())
 pk=inp('pk','ProteinaseK/SDS55°C30min;NaCl68°Covernight',M.CATALYST);x=op('reverse_crosslinks','Reverse crosslinks',K.DEPROTEINIZATION,x,((pk,R.CATALYST),),x.state)
 x=wash('precipitate','Split/pool bulkDNA ethanol/NaAc precipitation and70%ethanol washes',x)
 son=inp('sonication','Covaris~400bp; bulk PE2×100 reference');x=op('sonicate','Shear contact DNA',K.FRAGMENTATION,x,((son,R.CONDITION),),x.state)
 x=wash('spri','1×SPRI cleanup',x)
 beads=inp('beads','MyOneT1 streptavidin RT15min',M.SELECTION);bound=out('captured','biotin-enriched bead fraction',x.state);background=out('background','retained optional methylation background-control supernatant',k=M.MOLECULAR)
 g.react('capture','Capture biotin DNA; retain optional background',K.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,background));g.mark_terminal(Terminal(background.id,Disposition.STORED,'Source optionally ethanol-precipitates capture supernatant for separate control; downstream route outside scope.'))
 x=wash('capture_wash','Two55°C TWB washes and ligase-buffer exchange',bound)
 repair=inp('end_repair_mix','T4PNK/T4DNApolymerase/Klenow+dNTPs RT30min',M.CATALYST);x=op('end_repair','Repair ends/remove biotin from unligated ends',K.END_REPAIR,x,((repair,R.CATALYST),),x.state)
 x=wash('repair_wash','Source TWB washes and buffer2 exchange',x)
 tail=inp('tail_mix','Klenow exo−/dATP37°C30min',M.CATALYST);tailed=dA_tailed_scene(EndPreparedDuplex(contact_product(),five_prime_phosphate=True,three_prime_overhang='A'),label='dA-tailed original contact fragments')
 x=op('da_tail','Add3′dA ends',K.TAILING,x,((tail,R.CATALYST),),tailed);x=wash('tail_wash','TWB washes and QuickLigation buffer exchange',x)
 split=inp('control_split','optional10/40 split to conventionalHiC versus bisulfiteHiC; different adapterindexes');selected=out('bisulfite_aliquot','retained bisulfite aliquot',x.state);ctrl=out('normal_hic_control','optional normalHiC aliquot',x.state)
 g.react('control_partition','Optional physical10/40 control split',K.PARTITION,((x,R.TEMPLATE),(split,R.CONDITION)),(selected,ctrl));g.mark_terminal(Terminal(ctrl.id,Disposition.STORED,'Optional normalHiC control branch outside selected bulk methylation library; not destroyed.'))
 ad=inp('adapter','native indexed Illumina adapter; exact strands/protection unassigned',M.ADAPTER);al=inp('adapter_ligase','DNAQuickLigase RT15min',M.CATALYST)
 x=op('adapter_ligation','Ligate native adapters before conversion',K.LIGATION,selected,((ad,R.STRUCTURAL),(al,R.CATALYST)),adapted_scene(),AdapterLigationGeometry.unavailable('Selected bulk adapter duplex/protection not assigned; WGBS-control TruSeq names do not establish this adapter.'))
 x=wash('adapter_wash','TWB washes/Tris resuspension20µl',x)
 lam=inp('lambda','0.5%same-adapter-ligated unmethylated lambda; independent molecules',M.MOLECULAR);bs=inp('bisulfite','EZDNA Methylation-Gold kit native conversion',M.CONDITION)
 x=op('convert','Bisulfite-convert on beads with lambda spike',K.CHEMICAL_CONVERSION,x,((lam,R.SPIKE_IN),(bs,R.CONDITION)),converted_scene())
 recovery=inp('conversion_recovery_conditions','source column eluate/beads pooling boundary; exact component retention unassigned')
 y=out('purified_converted','actual recovered converted mixture',converted_scene());g.react('conversion_recovery','Recover convertedDNA for PCR',K.MAGIC,((x,R.TEMPLATE),(recovery,R.CONDITION)),(y,),reason='Main Methods purifiedDNA vs supplement column-plus-beads instruction; full transfer/retention unresolved, no physical source class discarded.');x=y
 f=inp('forward','native TruSeq-named PCR forward primer; sequence unassigned',M.OLIGO);r=inp('reverse','native TruSeq-named PCR reverse primer; member/extent unassigned',M.OLIGO);pol=inp('pcr_polymerase','KAPA HiFi HotStart Uracil+; supplement starts10cycles',M.CATALYST)
 amp=out('amplified','normal daughterDNA and possible original/control carryover',daughter_scene(),M.LIBRARY)
 g.react('library_pcr','Uracil-tolerant PCR on actual converted SS substrates',K.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(x.state,reason='Actual converted U-bearing contact/lambda strands; exact adapter protection and PCR primer boundaries unassigned.'))
 x=wash('size_selection','AMPure double0.6–0.15 selection mainly400–700bp',amp)
 g.sequence(x,SequencingPlan(library=FINAL_LIBRARY,primers=(),acquired_roles=('Read 1','Read 2'),limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(('sample','dpnii','digest','heat','heat_inactivate','klenow_fill','fill','contact_ligase','contact_ligation','pk','reverse_crosslinks','precipitate','precipitate_selection','sonication','sonicate','spri','spri_selection','beads','capture','capture_wash','capture_wash_selection','end_repair_mix','end_repair','repair_wash','repair_wash_selection','tail_mix','da_tail','tail_wash','tail_wash_selection','control_split','control_partition','adapter','adapter_ligase','adapter_ligation','adapter_wash','adapter_wash_selection','lambda','bisulfite','convert','conversion_recovery_conditions','conversion_recovery','forward','reverse','pcr_polymerase','library_pcr','size_selection','size_selection_selection','size_selection_out'),SOURCE)
 return g.require_publishable()
