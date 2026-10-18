"""Original2010 PARS: source-owned V1/S1 molecular branches and SOLiD boundary."""
from chemdraw import Construct,Scene,Segment,feature,oligo
from protocol_model import molecule,role
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition
import seqprimers as sp
TITLE='PARS — parallel analysis of RNA structure'
NOTES='01_pars.html'
SOURCE='Kertesz et al., Nature(2010), <a href="https://doi.org/10.1038/nature09322">doi:10.1038/nature09322</a>; published SupplementaryMethods90–143; SOLiD SmallRNAExpressionKit guide4399434 RevC(04/2009), AppendixA.'
SUMMARY='Separate V1/S1 digests mark structure-sensitive5′phosphates. Subsequent alkaline fragmentation and selective adapter ligation retain cleavage-defined read starts in independent SOLiD libraries.'
CAVEAT='This graph follows the defining2010 supplementary methods, not the later2013 protocol. Exact MixA ligating species/ends and native RT/PCR geometry remain unassigned. Known2009 kit common P1 primer is retained; a selected indexed reverse-primer member is not invented. Generated35–50bp sequence lengths are distinct from35bp mapping truncation and do not define color-space cycle chemistry.'
P1='CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT'
SEQ_PRIMERS=()
READ_LENGTHS=None
SEQUENCING_INTRO='Both V1 and S1 libraries undergo source SOLiD emulsionPCR/bead deposition. SupplementaryMethods reports35–50bp sequences generated; analysis truncates to35bp. Index acquisition/color-space run-cycle allocation is unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE='Native MixA sequencing entry/primer sites and instrument cycle chemistry are not reconstructed.'
READ_SPECS=(sp.read_spec('Sequencing output','35–50bp generated;35bp used for mapping','nuclease-cleavage-defined RNA start; source2009 reverse-primer barcode6nt, selected member/acquisition unassigned; color-space cycle geometry unassigned','include'),)

def endpoint(branch):
 return Construct([Segment('source common SOLiD P1 PCR end',P1),role(branch+' cleavage-defined DNA; insert geometry unassigned'),role('SOLiD P2 primer-derived end; six-base barcode location/member unassigned',molecular_feature=feature('pars_'+branch.lower()+'_index','sample_index','unknown'))],name=branch+' PARS source-limited SOLiD endpoint')
FINAL_LIBRARY=endpoint('V1')
FINAL_CAPTION='Representative V1 endpoint; independent S1 branch has its own terminal plan. Exact selected P2/adapter geometry is unassigned.'
def oligos():return (oligo('Common5′PCR primer —2009 kit AppendixA, all10sets',[Segment('common P1',P1)]),)
def rna_scene(label,*,five=False,three=False,mod5='',mod3=''):
 p=[]
 if five:p.append(role('MixA5′adapter; native sequence/geometry unassigned'))
 p.append(role(label))
 if three:p.append(role('MixA3′adapter; native sequence/geometry unassigned'))
 sc=Scene();sc.strand('RNA',p,label=label,mod5=mod5,mod3=mod3)
 if five:sc.junction('RNA',p[0].name,p[1].name,'RNA ligation')
 if three:sc.junction('RNA',p[-2].name,p[-1].name,'RNA ligation')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 initial=g.input('input','twice poly(A)-selected yeast RNA plus separate folded-RNA controls',MaterialKind.BIOLOGICAL,state=rna_scene('poly(A)+ transcript pool; independent RNA molecules'),primary=True)
 folding=g.input('folding','90°C,ice,slowRT in10mMTris/100mMKCl/10mMMgCl2',MaterialKind.CONDITION)
 folded=g.product('folded','in-vitro refolded RNA pool',MaterialKind.MOLECULAR,state=initial.state)
 g.react('refold','Refold RNA before structure-selective probing',StepKind.FOLDING,((initial,InputRole.TEMPLATE),(folding,InputRole.CONDITION)),(folded,))
 v1=g.product('v1_aliquot','independent RNaseV1 aliquot',MaterialKind.MOLECULAR,state=folded.state)
 s1=g.product('s1_aliquot','independent S1 aliquot',MaterialKind.MOLECULAR,state=folded.state)
 g.react('split','Split into matched independent V1/S1 aliquots',StepKind.SPLIT,((folded,InputRole.TEMPLATE),),(v1,s1),changes_molecule=False)
 for prefix,template,label,amount in (('v1',v1,'V1','0.01U RNaseV1'),('s1',s1,'S1','1000U S1 nuclease')):
  enzyme=g.input(prefix+'_nuclease',amount+';100µL,RT15min',MaterialKind.CATALYST)
  cut=g.product(prefix+'_cleaved',label+' cleaved RNA pool',MaterialKind.MOLECULAR,state=rna_scene(label+' nuclease-created cleavage fragment',mod5='P'))
  g.react(prefix+'_digest','Probe '+label+'-sensitive RNA positions',StepKind.DIGESTION,((template,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(cut,))
  buffer=g.input(prefix+'_alkaline','alkaline hydrolysis buffer;95°C3min',MaterialKind.CONDITION)
  fragstate=rna_scene(label+' cleavage5′P-retaining fragment;3′P/cyclicP heterogeneous',mod5='P',mod3='P/cyclicP; heterogeneous')
  fragstate.strand('randomOH',[role('random-fragmentation-only5′OH products')],label='independent non-cleavage-selected population',mod5='OH',mod3='P/cyclicP; heterogeneous')
  fragmented=g.product(prefix+'_fragmented',label+' alkaline-fragmented RNA populations',MaterialKind.MOLECULAR,state=fragstate)
  g.react(prefix+'_fragment','Alkaline-fragment probed RNA before adapter ligation',StepKind.FRAGMENTATION,((cut,InputRole.TEMPLATE),(buffer,InputRole.CONDITION)),(fragmented,))
  page=g.input(prefix+'_rna_page','6%denaturingPAGE;50–200nt excision,freeze/thaw,4°C overnightdiffusion,ethanol',MaterialKind.SELECTION)
  selected=g.product(prefix+'_size_selected',label+'50–200nt RNA fraction',MaterialKind.MOLECULAR,state=fragmented.state)
  off=g.product(prefix+'_offsize','RNA outside50–200nt selected window',MaterialKind.WASTE)
  g.react(prefix+'_size_select','Recover50–200nt fragments before selective5′ligation',StepKind.SELECTION,((fragmented,InputRole.TEMPLATE),(page,InputRole.SELECTION)),(selected,off));g.discard(off,'Not included in excised source RNA gel fraction.')
  mix=g.input(prefix+'_mixa','SOLiD SmallRNAExpressionKit AdapterMixA; selected species/ends unassigned',MaterialKind.ADAPTER)
  ligase=g.input(prefix+'_ligase','kit ligation enzyme mix;16°Covernight',MaterialKind.CATALYST)
  first=g.product(prefix+'_five_ligated',label+'5′adapter-ligated cleavage-target RNA',MaterialKind.MOLECULAR,state=rna_scene(label+' cleavage-defined RNA',five=True,mod3='P/cyclicP; heterogeneous'))
  other=g.product(prefix+'_unligated','unligated RNA including5′OH fragmentation products',MaterialKind.MOLECULAR)
  g.react(prefix+'_ligate5','Select nuclease-created5′P by first MixA ligation',StepKind.LIGATION,((selected,InputRole.TEMPLATE),(mix,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(first,other),geometry=AdapterLigationGeometry.unavailable('Original supplementary method specifies5′P-selective MixA ligation; native adapter mix sequences/ends and pairing geometry are unassigned.'))
  g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Failure to ligate is not physical disposal; subsequent phosphatase/second-ligation/recovery of these classes is not reconstructed.'))
  phosph=g.input(prefix+'_phosphatase','AntarcticPhosphatase;37°C1h,65°C7min',MaterialKind.CATALYST)
  dephos=g.product(prefix+'_dephosphorylated',label+'5′adapted RNA with3′OH',MaterialKind.MOLECULAR,state=rna_scene(label+' cleavage-defined RNA',five=True,mod3='OH'))
  g.react(prefix+'_dephosphorylate','Remove terminal phosphate after selective5′ligation',StepKind.DEPHOSPHORYLATION,((first,InputRole.TEMPLATE),(phosph,InputRole.CATALYST)),(dephos,))
  more=g.input(prefix+'_second_mixa','re-added source AdapterMixA plusligation system;16°C6h',MaterialKind.ADAPTER)
  dual=g.product(prefix+'_adapted',label+'dual-adapter RNA targets',MaterialKind.MOLECULAR,state=rna_scene(label+' cleavage-defined RNA',five=True,three=True))
  g.react(prefix+'_ligate3','Re-add MixA for second3′adapter ligation',StepKind.LIGATION,((dephos,InputRole.TEMPLATE),(more,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(dual,),geometry=AdapterLigationGeometry.unavailable('Source reports second MixA addition after phosphatase; actual ligating adapter species/end geometry remains unassigned, not borrowed from later SOLiD protocols.'))
  rtprimer=g.input(prefix+'_rt_primer','source kit RT primer; native sequence/landing unassigned',MaterialKind.OLIGO)
  rt=g.input(prefix+'_array_script','ArrayScript reverse transcriptase',MaterialKind.CATALYST)
  ss=molecule(label+' actual ss first-strand cDNA',(role('kit RT primer'),role('adapter-derived copy beyond RT site; extent unassigned'),role(label+' cleavage-defined cDNA'),role('copied5′adapter')))
  cdna=g.product(prefix+'_cdna',label+'first-strand cDNA target',MaterialKind.MOLECULAR,state=ss)
  rnas=g.product(prefix+'_rt_rna',label+'retained adapted RNA templates',MaterialKind.MOLECULAR,state=dual.state)
  g.react(prefix+'_rt','Reverse-transcribe '+label+' adaptedRNA',StepKind.REVERSE_TRANSCRIPTION,((dual,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,rnas))
  rnase=g.input(prefix+'_rnaseh','RNaseH; original supplementary method RNAremoval',MaterialKind.CATALYST)
  target=g.product(prefix+'_pcr_input',label+'ss cDNA target after RNaseH',MaterialKind.MOLECULAR,state=ss)
  debris=g.product(prefix+'_rna_debris','RNaseH RNA hydrolysis products',MaterialKind.MOLECULAR)
  g.react(prefix+'_remove_rna','Hydrolyse templateRNA with RNaseH',StepKind.DIGESTION,((cdna,InputRole.TEMPLATE),(rnas,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(target,debris));g.mark_terminal(Terminal(debris.id,Disposition.UNRESOLVED,'RNaseH digestion does not physically discard hydrolysates; later product-specific recovery is unassigned.'))
  f=g.input(prefix+'_pcr_forward','common SOLiD P1 PCR primer;2009RevC guide',MaterialKind.OLIGO,state=P1)
  r=g.input(prefix+'_pcr_reverse','SOLiD reverse-primer set member unassigned;kit6nt barcode',MaterialKind.OLIGO)
  pol=g.input(prefix+'_polymerase','source kit PCR;18–20rounds',MaterialKind.CATALYST)
  final=endpoint(label)
  lib=g.product(prefix+'_library',label+' amplified SOLiD library',MaterialKind.LIBRARY,state=final)
  g.react(prefix+'_pcr','Amplify actual '+label+' first-strand target18–20cycles',StepKind.PCR,((target,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(target.state,reason='Actual postRNaseH ss cDNA enters PCR. Native adapter/RT/reverse-primer species and first-copy geometry are not reconstructed; the expected final endpoint is not the template.'))
  g.sequence(lib,SequencingPlan(library=final,acquired_roles=('Read 1',),limitation=label+' branch. '+SEQUENCING_ENDING+' '+SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(('input','folding','refold','split','v1_nuclease','v1_digest','v1_alkaline','v1_fragment','v1_rna_page','v1_size_select','v1_mixa','v1_ligase','v1_ligate5','v1_phosphatase','v1_dephosphorylate','v1_second_mixa','v1_ligate3','v1_rt_primer','v1_array_script','v1_rt','v1_rnaseh','v1_remove_rna','v1_pcr_forward','v1_pcr_reverse','v1_polymerase','v1_pcr','v1_library','s1_nuclease','s1_digest','s1_alkaline','s1_fragment','s1_rna_page','s1_size_select','s1_mixa','s1_ligase','s1_ligate5','s1_phosphatase','s1_dephosphorylate','s1_second_mixa','s1_ligate3','s1_rt_primer','s1_array_script','s1_rt','s1_rnaseh','s1_remove_rna','s1_pcr_forward','s1_pcr_reverse','s1_polymerase','s1_pcr','s1_library'),SOURCE)
 return g.require_publishable()
