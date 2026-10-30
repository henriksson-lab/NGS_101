"""Protocol-owned, source-qualified incomplete preparation graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry,reported_adapter_end

def unknown(label):
 s=Scene();s.strand('preparation',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete formulation/geometry not assembled'):
 role=R.CATALYST if kind in {S.DIGESTION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
 st=g.input(key+'_stock',reagent,K.CATALYST if role==R.CATALYST else K.OLIGO if role==R.STRUCTURAL else K.CONDITION);y=g.product(key+'_product',label+'; all substrates/products/companion species and reaction remnants retained until stated workup',K.MOLECULAR)
 kw={'reason':'Known source operation; exact molecular states, ends and geometry not yet assembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
 g.react(key,label,kind,((x,R.TEMPLATE),(st,role)),(y,),**kw);return y

def cleanup(g,key,x,label):
 st=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source recovered preparation; recovery/carryover not presumed perfect',K.MOLECULAR);w=g.product(key+'_waste','Actual separated physical workup fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(st,R.CONDITION)),(y,w));g.discard(w,'Source physical separation, not computational filtering.');return y

def split(g,key,x,labels):
 ys=[g.product(key+'_'+str(i),lab,K.MOLECULAR) for i,lab in enumerate(labels)];g.react(key,key,S.SPLIT,((x,R.TEMPLATE),),ys);return ys

def assay(g,key,x,label):
 y,a=split(g,key,x,('Retained preparation after source assay allocation; assay amount unassigned',label));g.mark_terminal(Terminal(a.id,Disposition.QC,'Separate analytical material; no implicit return to preparation.'));return y

def pcr(g,key,x,label,primers):
 ps=[g.input(key+'_p'+str(i),lab,K.OLIGO,state=state) for i,(lab,state) in enumerate(primers)];e=g.input(key+'_enzyme',label+' source catalyst',K.CATALYST);y=g.product(key+'_product','New DNA copies plus original DNA/free primers/reagents',K.MOLECULAR);geom=PCRGeometry.unavailable(unknown('Whole operational PCR geometry not assembled'),reason='Source primer identities do not by themselves establish all inward sites/products.',symmetric_species=len(ps)==1);g.react(key,label,S.PCR,((x,R.TEMPLATE),*((p,R.STRUCTURAL) for p in ps),(e,R.CATALYST)),(y,),geometry=geom);return y

def oligo_scene(name,seq,mod5='',mod3=''):
 s=Scene();s.strand(name,[Segment(name,seq)],mod5=mod5,mod3=mod3);return s

def finish(g,x,limitation):
 lib=g.product('library','Source recovered sequencing preparation; complete molecular state unassembled',K.LIBRARY);g.react('library_identity','Source recovered sequencing library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=limitation));g.add_source_item(SourceItem('pending','Full source-specific substrate/adapter/end/read geometry and historical kit originals','Selected original scope only; ordinary unresolved work is not confirmed source absence.'));return g

from dataclasses import dataclass
from illumina import P7
SOURCE='Published PNAS2020 article primaryHTML Design of probes250–283/ExperimentalProtocol285–304/Computational310–320. Local generic/Cloud acquisition unavailable; primarybrowseroriginal read. SI Appendix unavailable, not confirmed absent. Not conventional HCR-chain or RNA-copy sequencing.'
INIT5='GGAGGGCAGCAAACGG'
INIT3='GAAGAGTCTTCCTTTACG'
PBS='CTCGACCGTTAGCAAAGCTC'
B1H1='CGTAAAGGAAGACTCTTCCCGTTTGCTGCCCTCCTCGCATTCTTTCTTGAGGAGGGCAGCAAACGGGAAGAG'
READOUT='CTTACGGATGTTGCACCAGCAAGAAAGAATGCGA'
BEAD_SUFFIX='GTTGGCACCAGGCTTACGGATGTTGCACCAGC'
READ1='GACACATGGGCGGAGCTTTGCTAACGGTCGAG'
INDEX1='GCTGGTGCAACATCCGTAAGCCTGGTGCCAAC'
INDEX2='CTCGACCGTTAGCAAAGCTCCGCCCATGTGTC'
INDEXING_NOTICE='Published sequencing paragraph says Index1=8/Index2=12, but custom bead has12J positions and computational demultiplexing calls Index1 bead barcode and Index2 well/sample barcode. ActualSI FigS2 andDatasetS3 unresolved. Keep both printed statements; no silently repaired mask.'
@dataclass(frozen=True)
class SelectedProbePair:
 target_reference:str
 inferred:bool=True
 def __post_init__(self):
  if len(self.target_reference)!=52 or set(self.target_reference)-set('ACGT'):raise ValueError('Explicit caller-selected52nt target DNA-alphabet reference required')
 def probes(self):
  left,right=self.target_reference[:25],self.target_reference[27:]
  five=Scene();five.strand('selected 5prime initiator',[Segment('source initiator',INIT5,inferred=self.inferred),Segment('source spacer','AA',inferred=self.inferred),Segment('right25nt reference homology',revcomp(right),inferred=self.inferred),Segment('10nt source eventUMI','N'*10,placeholder=True,inferred=self.inferred,feature=feature('hypr_umi','umi','random')),Segment('source primer binding site',PBS,inferred=self.inferred)],mod5='5Phos')
  three=Scene();three.strand('selected3primeinitiator',[Segment('left25nt reference homology',revcomp(left),inferred=self.inferred),Segment('source spacer','TA',inferred=self.inferred),Segment('source initiator',INIT3,inferred=self.inferred)])
  return five,three

def bead_stock():
 s=Scene();s.strand('source custom68mer bead oligo',[Segment('source P7',P7),Segment('12position bead identifier','N'*12,placeholder=True,feature=feature('hypr_bead_barcode','cell_barcode','unknown',note='Source writes12Js; multiple bead barcodes may be coencapsulated and merged computationally.')),Segment('source bead primer tail',BEAD_SUFFIX)],mod5='bead-linker-PC-linker');return s

def hybrid(g,key,x,name,state,label):
 pr=g.input(key+'_probe',name,K.OLIGO,state=state);y=g.product(key+'_product','Source target/probe complexes plus unbound/free reagent and all biological companion classes',K.MOLECULAR);g.react(key,label,S.HYBRIDIZATION,((x,R.TEMPLATE),(pr,R.STRUCTURAL)),(y,));return y

def protocol_graph():
 g=ProtocolGraph('HyPR-seq — source reporter ligation and droplet readout',authority=SOURCE,representation='published_instance')
 x=g.input('cells','Source dissociated cells containing chosen RNA targets; actual target probe panel/biologicalspecimen unassigned',K.BIOLOGICAL,primary=True)
 x=op(g,'fix',x,'Source5Mcells example:4%formaldehyde/PBS0.1%Tween20,1Mcells/ml,1hRT rocking')
 x=cleanup(g,'pbstrinse',x,'Source2x5mlPBS0.2%Tween wash;postfixation850g5minRT')
 x=op(g,'ethanol',x,'Source70%icecoldEtOH5ml,4C10min permeabilization')
 x=cleanup(g,'post_ethanol',x,'Source2x5mlPBST wash and transfer to2ml round-bottom tubes')
 x=op(g,'prehybridize',x,'Source500ul5xSSC/30%formamide/0.1%Tween37C5min')
 x=hybrid(g,'initiators',x,'Source two25nt RNA-homology initiators20nM each;5prime5Phos/10ntUMI/PBS; actual wholeDatasetS1 target probes unassembled',unknown('Actual whole source target-specific two-probe panel'), 'Source500ul pooled initiator probes37C overnight')
 x=cleanup(g,'probe_washes',x,'Source4x500ulhybridizationbuffer washes10min37C;500ul5xSSCT equilibration5minRT')
 x=hybrid(g,'b1h1',x,'Source MolecularInstrumentsB1H1 hairpin; separate95C90s then30minRT stock cooling',oligo_scene('B1H1',B1H1),'Source200ul75nM cooledB1H1,37C1h')
 x=cleanup(g,'b1h1_wash',x,'Source2x200ul5xSSCT washes after first hairpin; no HCR-chain primer default')
 x=hybrid(g,'readout',x,'Source adaptedB1H2 readout oligo, separately snapcooled95C90s→30minRT',oligo_scene('readout',READOUT),'Source200ul75nM cooledreadout37C1h')
 x=cleanup(g,'readout_wash',x,'Source3x500ul5xSSCT wash and200ul1xT4ligasebuffer equilibration')
 x=op(g,'ligate_reporter',x,'Source200ul1xT4ligasebuffer/1:100T4DNAligase,1hRT; readout already in input complex, no second invented addition',S.NICK_SEALING,'Source T4 DNA ligase; complete source hairpin/splint/junction geometry unassembled')
 x=cleanup(g,'ligated_cell_wash',x,'Source3x500ulPBST washes and20um filtration; retained cells carry ligatedDNA reporters, actualcellloss not predicted')
 x=assay(g,'count_cells',x,'Separate source hemocytometer/single-cell-suspension check allocation')
 used,rest=split(g,'cell_allocation',x,('Source2000cells for25ul emulsionPCR premix','Unused source ligated-cell stock; amount/fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No implicit whole-cell-stock consumption.'))
 beads=g.input('beads','Source20000custom68mer photocleavable primer beads per25ul premix',K.OLIGO,state=bead_stock());idx=g.input('indexing_primer','Source500nM indexing primer from actualDatasetS3; whole body/samplevalue unreviewed',K.OLIGO,state=unknown('Source DatasetS3 indexing primer'))
 mix=g.input('evagreen','Source12.5ul2xEvaGreenSupermix plus buffers/water to25ul premix; catalyst and droplets oil components not fully assembled',K.CONDITION);premix=g.product('pcr_premix','Source25ul cells+beads+indexingprimer+PCRreagents physical mixture',K.MOLECULAR);g.react('premix','Source co-loaded cell/bead/PCR mixture before emulsification',S.GENERIC,((used,R.TEMPLATE),(beads,R.STRUCTURAL),(idx,R.STRUCTURAL),(mix,R.CONDITION)),(premix,),changes_molecule=False)
 used,rest=split(g,'cartridge_load',premix,('Source20ul25ul premix loaded to cartridge','Source untransferred5ul premix; subsequentfate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Loaded aliquot, not entire25ul source reaction.'))
 x=op(g,'droplet_uv_pcr',used,'SourceQX200with70uloil;UV365nm6.5J/cm2 lamp5min releases bead primers;~50ul droplets transferred/sealed;94C30s;30x94C5s/64C30s/72C30s;72C5min;printed50%ramp2C/min retained')
 # Primers/catalyst were already co-loaded and remain inside the input mixture.
 # Keep this known compound operation opaque until full operational PCR geometry is assembled;
 # adding the external stocks again merely to satisfy PCR ports would duplicate material.
 # Selected source-allowed one-well breakage case; two-to-four well pools are a separate option.
 pool=x
 x=cleanup(g,'break_emulsions',pool,'Source40ul97%PFO,vortex5s/1000g1min; recover TOPaqueous layer, separate oil/otherphase waste')
 x=cleanup(g,'spri',x,'Source1.8xSPRI cleanup of aqueous amplicons')
 x=assay(g,'gel_qc',x,'Separate source gelQC around206bp and primer-dimercheck; not preparative gel size extraction')
 x=assay(g,'qubit_qc',x,'Separate source Qubitquantification aliquot')
 return finish(g,x,'SourceMiSeq6pMorNextSeq5501.8pM;Read1=35bp customREAD1,customINDEX1,both;NextSeq also customINDEX2. '+INDEXING_NOTICE+' ReporterUMIs label hybridization/ligation events; many probes can bind sameRNA, not oneUMI-per-nativeRNA guarantee. Bead merging is computational, not chemical barcode fusion. UnselectedactualRNAvariant/modification is not directly sequenced.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();p=SelectedProbePair('ACGT'*13);five,three=p.probes();c('Caller-selected reference defaults inferred',p.inferred,True);c('Selected 5prime probe uses downstream target25',five.strands['selected 5prime initiator'].segs[2].top,revcomp(p.target_reference[27:]));c('Cell and PCR stock remainders retained',all(k in g.terminals for k in ('cell_allocation_1','cartridge_load_1')),True);c('PCR mixture does not duplicate external primer additions',g.steps['droplet_uv_pcr'].inputs[0].material,'cartridge_load_0');c('Gel is separate QC, not invented preparative selection','gel_qc_1' in g.terminals,True);c('Exact index/read layout remains unassembled',bool(g.publication_problems()),True)
 src=Source(sources_root()/'hypr-seq__10.1073+pnas.2010738117'/'pnas.2010738117_PMC7776864.xml.txt','Publishedmainwholehairpin/readout/customsequencingprimers')
 if have(c,src,label='Locallycachedpublishedmain'):
  t=''.join(src.path.read_text().replace('\u200b','').split());c('Whole source B1H1/readout/readprimers verified',all(x in t for x in (B1H1,READOUT,READ1,INDEX1,INDEX2)),True)
 c.report()
if __name__=='__main__':selftest()
