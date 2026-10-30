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
 role=R.CATALYST if kind in {S.DIGESTION,S.DEPROTEINIZATION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
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
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2,TRUSEQ_SMALL_RNA_RA3
SOURCE='Published PMC13000175.1 main115/285–321;10Feb2026epubmetadata.18Mar2026publisher-version identity/body comparison unverified (actualNatureURL identityredirect). Downloadedpublished facts source-qualified, no corrected-body equivalence claim. SeqScope2024reference15, not imported2021constructor.'
ANCHOR_WRITTEN='/5Acryd/TT+TT+TT+TT+TT+TT+TT+T/3ATTO590N/'
@dataclass(frozen=True)
class SourceAnchor:
 written:str=ANCHOR_WRITTEN
 @property
 def lna_positions(self):
  seq=self.written.split('/')[2];positions=[];i=0
  for token in seq:
   if token=='+':positions.append(i)
   else:i+=1
  return tuple(positions)
 def scene(self):
  s=Scene();s.strand('source Anchor590ST',[Segment('source LNA-dT' if i in self.lna_positions else 'source dT','T',note='Source locked-nucleic-acid thymidine' if i in self.lna_positions else '') for i in range(15)],mod5='Acryd',mod3='ATTO590N');return s

def custom_antibody_tag():
 s=Scene();s.strand('source oYo custom tag',[Segment('sourceTR2body',TRUSEQ_READ2[11:]),Segment('source10ntAbbarcode','N'*10,placeholder=True,feature=feature('seqscopex_antibody','feature_barcode','unknown')),Segment('source10ntUMI N9B','N'*9+'B',placeholder=True,feature=feature('seqscopex_umi','umi','degenerate')),Segment('sourcepolyA32','A'*32,note='Source terminaltwo phosphorothioate bonds, ...A*A')]);return s

AT_F_WRITTEN='TCTTTCCCTACACGACGC*T*C'
AT_R_WRITTEN='CCTTGGCACCCGAGAATTC*C*A'
@dataclass(frozen=True)
class SourcePSPrimer:
 written:str
 @property
 def phosphorothioate_positions(self):
  positions=[];i=0
  for token in self.written:
   if token=='*':positions.append(i)
   else:i+=1
  return tuple(positions)
 def scene(self,name):
  s=Scene();s.strand(name,[Segment(name,self.written.replace('*',''),note='Source phosphorothioate ordering: '+self.written)]);return s
def psprimer(name,written):return SourcePSPrimer(written).scene(name)
def index_primer(name):
 if name=='P5TR1':left=P5;right=TRUSEQ_READ1[:27];written='AATGATACGGCGACCACCGAGATCTACAC [10-base index barcode] ACACTCTTTCCCTACACGACGCTCTT*C'
 else:left=P7;right=TRUSEQ_READ2[:12]+revcomp(TRUSEQ_SMALL_RNA_RA3);written='CAAGCAGAAGACGGCATACGAGAT [10-base index barcode] GTGACTGGAGTTCCTTGGCACCCGAGAATTCC*A'
 s=Scene();s.strand(name,[Segment('sourceflowcellhandle',left),Segment('10ntfixedsampleindex','N'*10,placeholder=True,feature=feature('seqscopex_'+name.lower(),'sample_index','fixed')),Segment('sourceprimerbody',right,note='Source terminal phosphorothioate bond: '+written)]);return s

def index_source_ordering(name):
 segs=index_primer(name).strands[name].segs
 return segs[0].top+'[10-baseindexbarcode]'+segs[2].top[:-1]+'*'+segs[2].top[-1]

def graph_for_assay(assay):
 if assay not in {'rna','totalseq_a'}:raise ValueError('Actual selected RNA orTotalSeq-A branchrequired;customoYo not silently substituted')
 g=ProtocolGraph('Seq-Scope-X — source '+assay+' expansion preparation',authority=SOURCE,representation='published_instance')
 x=g.input('section','Source10–16um frozen tissue section; actualspecimen unassigned; RNA and protein branches separate',K.BIOLOGICAL,primary=True)
 if assay=='totalseq_a':
  x=op(g,'antibody_stain',x,'Sourcemethanol−20C20min,washes/block1h;TotalSeqA Human/Mousecocktail1:3overnight4C;3x37C10minwashes;4%formaldehyde10minpostfixation. ExactADToligos/vendorconjugatesunassembled')
 else:
  x=op(g,'fix',x,'Source4%formaldehyde10minRT;3xPBSwashes;icecoldmethanol−20C20–30min;3x2xSSCwashes')
 anchor=g.input('anchor','SourceAnchor590ST2mM printedstock;15mer LNAoligodT with5Acryd/3ATTO590N',K.OLIGO,state=SourceAnchor().scene());reag=g.input('anchor_conditions','Source2xSSC/30%formamide/yeasttRNA1mgml/dextransulfate10%/SUPERaseIn1%;37C48h after formamidepreincubation',K.CONDITION);y=g.product('anchored_section','Source anchor/polyA target complexes plus unboundRNA/probes/reagents and biologicalmatrix; no perfect capture prediction',K.MOLECULAR);g.react('anchor_hybridize','SourcepolyA anchor hybridization;LNApattern andacrdyite retained',S.HYBRIDIZATION,((x,R.TEMPLATE),(anchor,R.STRUCTURAL),(reag,R.CONDITION)),(y,));x=cleanup(g,'anchor_wash',y,'Source3x2xSSCwash beforegelembedding;unbound fractionsphysicallyremoved')
 x=op(g,'embed',x,'Source2MNaCl/12%sodiumacrylate/3%19:1acryl-bis/60mMTrismonomer;0.2%APS+0.2%TEMEDinitiation,no4hydroxyTEMPOdelay;1–2hRT polymerization;acrydite probesincorporated')
 x=op(g,'clear_tissue',x,'Source3–4.5ml2xSSC/2%SDS/0.5%Triton/1%ProteinaseK37C24h',S.DEPROTEINIZATION,'Source ProteinaseK/SDS clearing; target nucleic acids not asserted nuclease-digested')
 x=cleanup(g,'clear_wash',x,'Source digestionbufferremoved/1xSSCgelwash;RNA/protease/waste recovery notperfect')
 x=op(g,'expand',x,'Source0.1xSSC15minexchange/15minexchange/~30minplateau;DAPI final30min;notpurewater;expansionfactor/measuredgeometryunassigned')
 chip=g.input('chip','SourceindependentlypreparedNovaSeq6000S4 SeqScope2024ref15chip withHDMI/oligodT30/TR1;completearrayoligo/sequencedcoordinate mapunassembled',K.OLIGO,state=unknown('Exact2024sourceSeqScopearraymoleculargeometry notborrowedfrom2021construct'))
 assembly=g.product('gel_chip','Source trimmed~7x7mm gel/SeqScopechip sandwich including all gel/array populations',K.MOLECULAR);g.react('sandwich','Source expandedgeloverlaidontospatialarray',S.HYBRIDIZATION,((x,R.TEMPLATE),(chip,R.STRUCTURAL)),(assembly,))
 x=op(g,'partner_switch',assembly,'Sourcehumidity-maintained45C30min:releasefrom15merLNAanchor andcaptureby30merarraydT;sourceTm~39versus>50C,notathermodynamic-yieldprediction')
 kept=g.product('captured_chip','Source retrievedchip withcapturedRNA orADTDNA andcarryover;fullmolecularstatesunassembled',K.MOLECULAR);gel=g.product('removed_gel','Physicalremovedgelwithanchors/residualuntransferredtargets/reagents;actualsubsequentfateunassigned',K.MOLECULAR);g.react('remove_gel','Sourcecarefulphysicalgelremovalaftercapture',S.SPLIT,((x,R.TEMPLATE),),(kept,gel));g.mark_terminal(Terminal(gel.id,Disposition.UNRESOLVED,'Removedgelisnotdeclaredperfectlyemptyorautomaticallydiscarded.'));x=kept
 if assay=='rna':
  x=op(g,'cited_library',x,'Source2024SeqScopeprotocol15 startingRT andfollowingsteps; OMITS tissueProteinaseK digestionbecausenorealtissueonchip. ActualsourceRT/primers/cleanup/rununassembled,notborrowedfrom2021method')
  return finish(g,x,'SourceRNA spatialcDNA library;exact2024array/primer/run/indexgeometry remainsunassembled. Expansion/heat-mediatedpartner-switch isactualprepboundary,notonlyscaledarrayresolution. SpatialcoordinateHDMI andmoleculeidentifiers needactual2024sourceconstruct;noassumedbarcodelengths.')
 x=op(g,'adt_copy',x,'SourceRTenzymeDNA-polymeraseactivity extendsADTDNAprimerfrompolyAtail usingCHIPDNAtemplate encodingoligodT/HDMI/TR1;notRNA reverse-transcriptionorordinarybeadprimerdirection',S.EXTENSION,'SourcealreadydeclaredSeqScopeRTenzyme/buffer;actualhistoricalidentity/geometryunassembled')
 x=cleanup(g,'adt_elute',x,'Source2x30ul0.1NNaOH5minelutions(total60ul),neutralize30ul3MKOAcpH5.5;OligoCleanConcentrator40uleluent. Sourcearrayretainedphysicalfractionnotperfecttargetrelease')
 x=pcr(g,'at_pcr',x,'Source200ulATPCR:100ul2xKAPA+1uleachprinted100mMPSprimers+40ulprimarylibrary+58ulwater;95C3min;15x95C30s/60C1min/72C1min;72C2min',(('SourceRPEPCRF',psprimer('F',AT_F_WRITTEN)),('SourceSR2RPEPCRR',psprimer('R',AT_R_WRITTEN))))
 x=cleanup(g,'at_clean',x,'Source1.8xAMPure,30ulelution')
 used,rest=split(g,'index_aliquot',x,('Source3ulATPCRproduct toindexPCR','Untransferred27ulATPCRstock;subsequentfateunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actual3of30ultransferdoesnotconsumewholeATstock.'))
 x=pcr(g,'index_pcr',used,'Source40ulindexPCR:20ul2xKAPA/4uleachprinted10mMP5TR1andP7SR1/3ulproduct/9ulwater;source6–8cyclesbasedonconcentration,exactchosenbranchunassigned',(('SourceP5TR1PSprimer',index_primer('P5TR1')),('SourceP7SR1PSprimer',index_primer('P7SR1'))))
 x=cleanup(g,'index_clean',x,'Source1xAMPure30uleluate andagarose250bpband/ZymocleanGelDNARecovery')
 return finish(g,x,'SourceTotalSeqAproteinreadoutIlluminaNovaSeqX;dual10ntsampleindexprimerssourcepatterns,actualindexmembers/cyclemaskunassigned. Nativeprotein isreadthroughADTDNA;customoYoTR2oligos/primersremainseparatevariant,notinterchangedwithTotalSeqASR2. Printed100mM/10mMprimerstocksnotnormalizedbyguess. Primarycopydirectionandfullreadgeometryunassembled.')
def protocol_graph():return graph_for_assay('rna')
def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();r=protocol_graph();p=graph_for_assay('totalseq_a');c('Removedgelphysicalfateretained','removed_gel' in r.terminals,True);c('ATstockremainderretained','index_aliquot_1' in p.terminals,True);c('ProteinbranchusesDNAextension',p.steps['adt_copy'].kind,S.EXTENSION);c('No2021SeqScopeconstructimportorcompleteclaim',bool(r.publication_problems()),True);c('CustomUMIretainsconstrainedterminalB',custom_antibody_tag().strands['source oYo custom tag'].segs[2].top,'N'*9+'B')
 src=Source(sources_root()/'seq-scope-expanded__10.1038+s41467-026-69346-8'/'s41467-026-69346-8_PMC13000175.xml.txt','Publishedmainhybridization/TotalSeqAprimers;versionqualificationexplicit')
 if have(c,src,label='PublishedPMCmain'):
  t=''.join(src.path.read_text().split());c('Wholeanchororderingwithmodificationssourceexact',ANCHOR_WRITTEN in t,True);c('BothATPCRPSprimerbodiessourceexact',all(s in t for s in (AT_F_WRITTEN,AT_R_WRITTEN)),True);c('WholeP5TR1andP7SR1sourceorderallfixedbasesindexwindowandPS',all(index_source_ordering(name) in t for name in ('P5TR1','P7SR1')),True)
 c.report()
if __name__=='__main__':selftest()
