"""Vendor-revision-owned native molecular declaration."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, MaterialKind as M, InputRole as I, StepKind as K, SequencingPlan, Disposition, Terminal, SourceDiscrepancy, SourceClaim, AnnotationSlot
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from endprep import EndPreparedDuplex,dA_tailed_scene
import seqprimers as sp
SEQ_PRIMERS=()
ACQUIRED_ROLES=()
SEQUENCING_INTRO=''
SEQUENCING_UNAVAILABLE='Native complete primer/adapter strands and selected run are unassigned in this bounded source pass.'
SEQUENCING_ENDING='Illumina-compatible source kit endpoint; complete native primer geometry is unassigned.'
CAVEAT='Unassigned native adapter/primer regions are role tokens. Trim motifs constrain read-through sequence; they do not determine complete kit strands, ends or index placement. Recommendations are distinguished from an acquired run.'
SHOW_READ_CONTENT=False
I5=feature('sample_i5','sample_index','unknown')
I7=feature('sample_i7','sample_index','unknown')
def ss(label,parts):return molecule(label,tuple(parts))
def endpoint(label,payload,left,right):return Construct((left,*payload,right),name=label)
def recover(g,key,template,label,state=None,kind=M.MOLECULAR):
 reagent=g.input(key+'_selection',label,M.SELECTION)
 kept=g.product(key+'_retained','recovered '+template.label,kind,state=template.state if state is None else state)
 waste=g.product(key+'_waste','physically removed supernatant/washes or excluded bead fraction',M.WASTE)
 g.react(key,label,K.CLEANUP,((template,I.TEMPLATE),(reagent,I.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Physical non-retained fraction; complete species-specific depletion is not inferred.')
 return kept

def pcr(g,key,template,state,label,forward,reverse,polymerase,symmetric=False):
 product=g.product(key+'_product',label+' product',M.MOLECULAR,state=state)
 residual=g.product(key+'_residual','original-template and other co-carried molecular classes; recovery unassigned',M.MOLECULAR,state=template.state)
 inputs=[(template,I.TEMPLATE),(forward,I.STRUCTURAL),(polymerase,I.CATALYST)]
 if reverse is not None:inputs.insert(2,(reverse,I.STRUCTURAL))
 g.react(key,label,K.PCR,inputs,(product,residual),geometry=PCRGeometry.unavailable(template.state,reason='Actual incoming native precursor is shown; first-cycle/end geometry remains unvalidated. Completed library is a separate output.',symmetric_species=symmetric))
 g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'PCR does not chemically erase original templates or physically purify co-carried polymers; later recovery is unassigned.'))
 return product

def sequence(g,current,endpoint):
 lib=g.product('library','source-qualified sequencing library',M.LIBRARY,state=endpoint)
 g.react('retain_library','Retain recovered library for sequencing',K.GENERIC,((current,I.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(library=endpoint,limitation=SEQUENCING_UNAVAILABLE))
 return lib

TITLE='NEBNext Enzymatic Methyl-seq v2'
NOTES='01_nebnext-em-seq-v2.html'
SOURCE='<a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale8015.pdf?rev=224af38f17a145749b296dab1c27cea0">NEB E8015 manual v2.0_7/26</a>; LV UDI manual v2.0_6/26.'
SUMMARY='Selected ≤10ng input branch: separate sample/lambda/pUC19 DNA, post-ligation carrier addition, TET2/T4-BGT protection, physical C-to-U deamination and Q5U copying.'
TRIM_READ1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
TRIM_READ2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
LEFT=role('native LV UDI end A; physical i5 N8 member/position unassigned',molecular_feature=I5)
RIGHT=role('native LV UDI end B; physical i7 N8 member/position unassigned',molecular_feature=I7)
FINAL_LIBRARY=endpoint('native EM-seq v2 library',(role('sample/control DNA copied from converted originals; U readout T/A and protected sites C/G'),),LEFT,RIGHT)
FINAL_CAPTION='Native EM-seq endpoint. Physical eight-base LV indexes are source-known; native outer strands, selected well and acquired cycles remain unassigned.'
READ_SPECS=(sp.read_spec('Sequencing output',None,'Converted DNA: protected5mC/5hmC reports C, unmodified genomic C reports T after copying; read/index acquisition unassigned','include'),)
NAMES=('sample DNA','unmethylated lambda control','CpG-methylated pUC19 control')
def population(stage,*,adapted=False,ss_dna=False,carrier=False):
 sc=Scene()
 for n,name in enumerate(NAMES):
  description=stage
  if name=='unmethylated lambda control' and '5caC/5gmC' in stage:
   description=('APOBEC C-to-U target in unmethylated control; no protected-site assignment; efficiency unassigned' if 'U at' in stage else 'unmethylatedC after protection reaction; no protected-site assignment')
  for strand in range(2):
   parts=[role(name+' '+description+' '+str(strand+1))]
   if adapted:parts=[role('native EM adaptor boundary '+name+' '+str(strand)),*parts,role('opposite native EM adaptor boundary '+name+' '+str(strand))]
   sc.strand(str(n)+'_'+str(strand),parts,label=('separate physical ssDNA: ' if ss_dna else 'strand class of dsDNA: ')+name)
 if carrier:sc.strand('carrier',[role('independent carrier DNA; native protection/deamination/recovery unassigned')],label='carrier added after adapter ligation; no invented adapter ends')
 return sc

def dA_state():
 rows=[]
 for name in NAMES:
  prep=EndPreparedDuplex(Construct([role(name+' prepared insert')]),five_prime_phosphate=False,three_prime_overhang='A',overhang_inferred=False)
  rows.extend(dA_tailed_scene(prep,label=name+' source3-primeA;5-prime chemistry unassigned').rows())
 return rows

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 dna=g.input('dna','selected0.1–10ng sample dsDNA',M.BIOLOGICAL,state=ss('sample DNA',(role('sample dsDNA'),)),primary=True)
 lam=g.input('lambda','unmethylated lambda control',M.MOLECULAR)
 puc=g.input('puc','CpG-methylated pUC19 control',M.MOLECULAR)
 shear=g.input('shear','UltraShear or Covaris fragmentation to approximately350bp',M.CONDITION)
 frag=g.product('fragmented','fragmented sample and independent controls',M.MOLECULAR,state=population('fragmented DNA'))
 g.react('fragment','Combine controls and fragment §1.1',K.FRAGMENTATION,((dna,I.TEMPLATE),(lam,I.SPIKE_IN),(puc,I.SPIKE_IN),(shear,I.CONDITION)),(frag,))
 prep=g.input('endprep','UltraII EndPrep enzymes:20C15min/65C15min',M.CATALYST)
 repaired=g.product('repaired','end-repaired three DNA populations',M.MOLECULAR,state=population('end-repaired DNA'))
 g.react('repair','Repair fragment ends §1.2',K.END_REPAIR,((frag,I.TEMPLATE),(prep,I.CATALYST)),(repaired,))
 tailed=g.product('tailed','opposed3-prime dA DNA ends',M.MOLECULAR,state=dA_state())
 g.react('tail','Add3-primeA in same endprep reaction',K.TAILING,((repaired,I.TEMPLATE),(prep,I.CATALYST)),(tailed,))
 adapter=g.input('adapter','native EM-seq adaptor E7165AVIAL; strands unassigned',M.ADAPTER)
 ligase=g.input('ligase','UltraII ligation mix/enhancer:20C15min',M.CATALYST)
 adapted=g.product('adapted','native adaptor-ligated DNA',M.MOLECULAR,state=population('adapter-ligated insert',adapted=True))
 g.react('ligate','Ligate native EM adaptor §1.3',K.LIGATION,((tailed,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable('Source3-primeA substrate known; native adapter strands/terminal activation and cytosine protection are unassigned in this pass.'))
 ligated=recover(g,'ligated_cleanup',adapted,'§1.4 1X93µl beads; selected≤10ng branch27µl elution')
 carrier=g.input('carrier','§1.4.9B add1µl NEBNext CarrierDNA after recovery',M.MOLECULAR)
 mixed=g.product('carrier_added','three adapted DNA populations plus independent carrier',M.MOLECULAR,state=population('unconverted adapted DNA',adapted=True,carrier=True))
 g.react('add_carrier','Add carrier for selected≤10ng branch',K.GENERIC,((ligated,I.TEMPLATE),(carrier,I.SPIKE_IN)),(mixed,),changes_molecule=False)
 protect=g.input('protect','TET2/T4-BGT/UDP-glucose/FeII:37C1h; source low-input dilution',M.CATALYST)
 protected=g.product('protected','physical5caC/5gmC protected DNA',M.MOLECULAR,state=population('protected5caC/5gmC plus unmodifiedC',adapted=True,carrier=True))
 g.react('protect_dna','TET2 oxidation and T4-BGT glucosylation §1.5',K.CHEMICAL_CONVERSION,((mixed,I.TEMPLATE),(protect,I.CATALYST)),(protected,))
 stop=g.input('stop','StopReagent37C30min',M.CONDITION)
 stopped=g.product('stopped','stopped protection reaction',M.MOLECULAR,state=protected.state)
 g.react('stop_oxidation','Stop oxidation',K.QUENCHING,((protected,I.TEMPLATE),(stop,I.CONDITION)),(stopped,),changes_molecule=False)
 cleaned=recover(g,'protected_cleanup',stopped,'§1.6 1X50µl beads;16µl protected DNA')
 denat=g.input('denaturant','selected formamide branch:4µl/85C10min then ice; NaOH alternative separate',M.CONDITION)
 singles_state=population('physical protected original strand; unmodifiedC',adapted=True,ss_dna=True)
 top_state=Scene();bottom_state=Scene()
 for key,strand in singles_state.strands.items():
  target=top_state if key.endswith('_0') else bottom_state
  target.strand(key,strand.segs,label=strand.label)
 singles=g.product('single_strands_top','protected original top-strand classes',M.MOLECULAR,state=top_state)
 bottoms=g.product('single_strands_bottom','protected original bottom-strand classes',M.MOLECULAR,state=bottom_state)
 carrier_original=g.product('carrier_denatured','independent native carrier population after denaturation',M.MOLECULAR,state=ss('carrier',(role('carrierDNA; strands/chemical state unassigned'),)))
 g.react('denature','Denature protected DNA §1.7',K.STRAND_SEPARATION,((cleaned,I.TEMPLATE),(denat,I.CONDITION)),(singles,bottoms,carrier_original))
 enzyme=g.input('apobec','APOBEC and deamination reagents §1.8; no cleanup before PCR',M.CATALYST)
 converted=g.product('converted','physical C-to-U original strands',M.MOLECULAR,state=population('U at formerly unmodifiedC;5caC/5gmC retained',adapted=True,ss_dna=True,carrier=True))
 g.react('deaminate','Deaminate unmodifiedC to physicalU',K.CHEMICAL_CONVERSION,((singles,I.TEMPLATE),(bottoms,I.TEMPLATE),(carrier_original,I.TEMPLATE),(enzyme,I.CATALYST)),(converted,))
 f=g.input('forward','LV UDI forward; source physicalN8 index, selectedwell/strand unassigned',M.OLIGO)
 r=g.input('reverse','LV UDI reverse; source physicalN8 index, selectedwell/strand unassigned',M.OLIGO)
 pol=g.input('q5u','Q5U source input-dependent4–14cycles §1.9',M.CATALYST)
 amp=pcr(g,'pcr',converted,FINAL_LIBRARY,'Q5U copying of converted strands',f,r,pol)
 clean=recover(g,'pcr_cleanup',amp,'§1.10 0.8X72µl bead cleanup')
 sequence(g,clean,FINAL_LIBRARY)
 g.annotate(SourceDiscrepancy('pcr','LV index kit identifiers in E8015',(
 SourceClaim('E8015v2.0_7/26 required materials/overview','E3390 and E3392'),SourceClaim('same manual Fig1caption','E3900 and E3902'),SourceClaim('LVUDI manualv2.0_6/26 title/indexsets','E3390 and E3392')),'Use the identified E3390/E3392 kits; selected well remains unassigned.','Required-materials list agrees with dedicated LVUDI original.','No vendor correction confirmed; older revisions not compared.',slot=AnnotationSlot.INTRO))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('dna','lambda','puc','shear','fragment','endprep','repair','tail','adapter','ligase','ligate','ligated_cleanup_selection','ligated_cleanup','carrier','add_carrier','protect','protect_dna','stop','stop_oxidation','protected_cleanup_selection','protected_cleanup','denaturant','denature','apobec','deaminate','forward','reverse','q5u','pcr','pcr_cleanup_selection','pcr_cleanup','retain_library','library')
