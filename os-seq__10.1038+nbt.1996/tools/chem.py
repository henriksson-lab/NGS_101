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
from illumina import TRUSEQ_READ1
COMMON='CGAGATCT'+TRUSEQ_READ1
READ1=TRUSEQ_READ1[4:]
READ2='CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
SOURCE='PMC4336783 author manuscriptMethods68–98;publishedSI MOESM13ESM primerTablep2/moleculareventspp3–5 controls exactoligos. Finalmainnotcompared.'
@dataclass(frozen=True)
class SourceAdapter:
 index:str|None=None
 def __post_init__(self):
  if self.index is not None and self.index not in {'TGCTAA','AGGTCA','AACCTG'}:raise ValueError('Only source-knownselected6baseindexes;no syntheticdefault')
 @property
 def duplex_top(self):return COMMON+self.index if self.index else COMMON[:-1]
 @property
 def top(self):return self.duplex_top+'T'
 @property
 def bottom(self):return revcomp(self.duplex_top)
 @property
 def ps_after(self):return len(self.top)-2
 def scene(self):
  sc=Scene()
  if self.index:
   index_feature=feature('os_seq_'+self.index.lower(),'sample_index','fixed')
   top=[Segment('pairedtop',COMMON),Segment('source sample index',self.index,feature=index_feature),Segment('source3T','T')]
   bottom=[Segment('sample index complement',revcomp(self.index),feature=index_feature),Segment('pairedbottom',revcomp(COMMON))]
  else:top=[Segment('pairedtop',self.duplex_top),Segment('source3T','T')];bottom=[Segment('pairedbottom',self.bottom)]
  sc.strand('top',top,mod3='PS bond before terminal T')
  sc.anneal('bottom',bottom,to='top',pair=('pairedbottom','pairedtop'),mod5='P');sc.overhang('top','source3T',polarity='3-prime',text='Source 3primeT cohesive end');return sc

def protocol_graph():
 g=ProtocolGraph('OS-Seq — originalunindexedcolumn-probe capture case',authority=SOURCE,representation='published_instance')
 x=g.input('genomic','Source1ug genomicDNA; exactsampleunassigned; optionalindexedcasesseparate',K.MOLECULAR,primary=True)
 x=op(g,'fragment',x,'CovarisE210R source500bpmeanfragments',S.FRAGMENTATION)
 x=op(g,'end_repair',x,'50ulKlenow/T4DNApolymerase/T4PNK,dNTP,T4ligaseATPbuffer25C45min',S.END_REPAIR)
 x=op(g,'tail',x,'80ulTaq3.2U,dATP100uM,Mg1.5mM,72C15min',S.TAILING)
 x=cleanup(g,'prelig_cleanup',x,'SourceFermentasPCRpurificationbeforeligation')
 inv=SourceAdapter();ad=g.input('adapter','Sourceunindexedadapterduplex,topterminalPSbond/bottom5phosphate',K.ADAPTER,state=inv.scene(),chemical_state=inv);e=g.input('ligase','2000UT4DNAligase/ATPbuffer25C1h;sourceadapterannealing85→20Cthen4Cseparatestockconditioning',K.CATALYST);y=g.product('ligated','Sourceadapter-endedDNAandfreeadapter/otherligationclasses/remnants;wholejunctionunassembled',K.MOLECULAR);g.react('ligation','Sourcesingle-adapterligation',S.LIGATION,((x,R.TEMPLATE),(ad,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=reported_adapter_end('OSadapter',kind='cohesive',polarity='3-prime',terminal_sequence='T',substrate_scenes=(inv.scene(),),statement='PublishedSITableactualterminalPS-T/bottom5P andmain3A-tailedDNA;fullselectedgenomicjunctionunassembled.'))
 x=cleanup(g,'postlig_cleanup',y,'SourceFermentasPCRpurification');x=pcr(g,'library_pcr',x,'PhusionHotStart,1uM SINGLEsourceAd_topprimer,20–22cycles98C10s/65C30s/72C30s;selectedexactcyclesunassigned',(('SourceAd_top_FC_capture_A_tailPCRprimer',oligo_scene('Ad_top',inv.top,mod3='PSbondbeforeterminalT')),))
 x=cleanup(g,'postpcr_cleanup',x,'SourceFermentasPCRpurification');x=assay(g,'library_quant',x,'Sourceconcentrationassay;method/amountunassigned')
 lawn=g.input('flowcell','IndependentsourceIlluminaGAIIxflowcellC/Dsurfaceprimers;C2deoxyU/D8oxoG/surfacechemistryunassembled',K.MOLECULAR,primary=True);probes=g.input('capture_oligos','Sourcecolumn-synthesizedOSSeq366target40meroligopool;actualmembersunreviewed;notUMIs',K.OLIGO);st=g.input('flowcell_conditions','Sourceoligostock95C5min;100nM30ul/lane96→40C18min;Illuminav4extension/wash/denature',K.CONDITION);modified=g.product('modified_lawn','ExtendedimmobilizedDprimer-probes plusremainingC/Dlawn;removedtemplateoligosnotincorporatedasDNA',K.MOLECULAR);removed=g.product('probe_templates_removed','Physicallywashedsourcecaptureoligotemplates/inactivecomplementarystrands',K.WASTE);g.react('modify_lawn','Source Dprimer extension fromhybridizedcaptureoligo then templatewash',S.MAGIC,((lawn,R.TEMPLATE),(probes,R.STRUCTURAL),(st,R.CONDITION)),(modified,removed),reason='Actualsourceanneal/extension/denature sequence;surfaceoligo8oxoG/deoxyU andcompletegeometrynotassembled.');g.discard(removed,'Sourceoriginalcapturetemplateandsecondstrandwashedoffflowcell;newcovalentDextensionretained.')
 library_stock=x;x,rest=split(g,'capture_input',library_stock,('Source22.5ul40–56.6ng/ullibraryallocation;95C5min thenbuffer30ultotal','Unusedlibraryamount/fateunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Noimplicitentirestockconsumption.'))
 mix=g.product('flowcell_mixture','Sourcegenomiclibrary plusmodifiedprimerlawn',K.MOLECULAR);g.react('load_capture','Source30ul/lanelibraryandcaptureflowcell',S.POOL,((x,R.TEMPLATE),(modified,R.TEMPLATE)),(mix,));x=op(g,'target_capture',mix,'Source65C20hhybridization andv4extension;originalgenomictemplatesretaineduntilactualwash',S.MAGIC,'SourceIlluminav4extensionreagents;Dprimer-probesalreadycovalentlypresentininputlawn,nofreeprimeradded')
 x=cleanup(g,'template_wash',x,'SourcepostextensiondenatureandwashremovesoriginalgenomictemplateDNA')
 x=op(g,'bridge',x,'Source96→40C18minhybridizelibrary12basetailtoC,extendbothC/targetends;twosurfacecopies',S.EXTENSION,'SourceIlluminav4DNApolymerase;fullC/Dsurfacechemistryunassembled')
 x=op(g,'cluster_processing',x,'Sourcestandardv4bridgeamplification/clusterprocessing;actualCdeoxyU/D8oxoGhandlingnotassembled')
 return finish(g,x,'GAIIxOSSeq366PE40/OSSeq11kPE60 areseparateknownsourcecases;selectedcolumn366casePE40. Read1genomicDNA,Read2syntheticprimer-probe,notassumedbothindependentgenomicmoleculeends. Optionalindexes6basesnotmodelUMIs;historicalC/Dmodifiedsurfaceandfullflowcell/readgeometryunassembled.')
def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check()
 for index in (None,'AACCTG','TGCTAA','AGGTCA'):
  a=SourceAdapter(index);c('Sourceadapterpairingbuilds '+str(index),bool(a.scene()),True);c('PSbondretainsactualterminalT '+str(index),a.top[a.ps_after+1:],'T')
 g=protocol_graph();c('SingleprimerPCRdeclaredassingleinput',sum(i.role==R.STRUCTURAL for i in g.steps['library_pcr'].inputs),1);c('Probe-templateworkupisphysical','probe_templates_removed' in g.terminals,True);c('Unassembledhistoricalsurfacepreventsfinishedmarker',bool(g.publication_problems()),True)
 src=Source(sources_root()/'os-seq__10.1038+nbt.1996'/'nbt.1996_41587_2011_BFnbt1996_MOESM13_ESM.pdf.txt','PublishedSIprimerTablep2;renderoriginalmodifications')
 if have(c,src,label='PublishedSI'):
  t=''.join(src.path.read_text().split());c('Wholemodifiedunindexedtopandbottom',COMMON[:-1]+'xT' in t and 'p-'+SourceAdapter().bottom in t,True);c('Wholeindexedtopandphosphatebottom',all((COMMON+i+'xT') in t and ('p-'+SourceAdapter(i).bottom) in t for i in ('AACCTG','TGCTAA','AGGTCA')),True)
 c.report()
if __name__=='__main__':selftest()
