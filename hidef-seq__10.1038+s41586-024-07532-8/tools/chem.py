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
from hairpin import Hairpin

@dataclass(frozen=True)
class SelectedHairpin(Hairpin):
 inferred:bool=True

SOURCE='PMC11216816 author manuscriptMethods316–350/402–406 andpublishedSI Note1;selectednative restrictionfragment/nonAtailedcase,notmodernSMRTdefault.'
ADAPTERS={
'1001':('ACGCACTCTGATATGTGATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGATCACATATCAGAGTGCGT','CACATATCAGAGTGCG'),
'1002':('ACTCACAGTCTGTGTGTATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGATACACACAGACTGTGAGT','ACACACAGACTGTGAG'),
'1003':('ACTCTCACGAGATGTGTATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGATACACATCTCGTGAGAGT','ACACATCTCGTGAGAG'),
'1008':('ACGCAGCGCTCGACTGTATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGATACAGTCGAGCGCTGCGT','ACAGTCGAGCGCTGCG')}
def source_hairpin(name='1001'):
 seq,_=ADAPTERS[name];n=max(i for i in range(1,len(seq)//2) if seq[:i]==revcomp(seq[-i:]));return SelectedHairpin('SourcebcAd'+name,seq,n)
def hairpin_scene(name='1001'):
 seq,bc=ADAPTERS[name];feat=feature('hidef_bc'+name,'sample_index','fixed');spans=sorted([(seq.index(bc),bc),(seq.index(revcomp(bc)),revcomp(bc))]);segments=[];i=0
 for start,bases in spans:
  if start>i:segments.append(Segment('adapterbody'+str(i),seq[i:start]))
  segments.append(Segment('samplebarcode'+str(start),bases,feature=feat));i=start+len(bases)
 if i<len(seq):segments.append(Segment('adapterbody'+str(i),seq[i:]))
 s=Scene();s.strand('SourcebcAd'+name,segments,mod5='P');return s

def graph_for_mode(mode):
 if mode not in {'blunt','a_tailed'}:raise ValueError('Sourceblunt/nonAtailed orsourceA-tailed mode required')
 g=ProtocolGraph('HiDEF-seq — '+mode+'sourcecase',authority=SOURCE,representation='published_instance');x=g.input('genomic','SourcenativeextractedDNA500–3000ng;noPCRorcompleteDNAdamagerepairassumed',K.MOLECULAR,primary=True);x=assay(g,'input_qc',x,'SourceNanoDrop/Qubitinputdiagnosticallocations')
 x=op(g,'hpy166ii',x,'Hpy166II0.14U/ul,70ulCutSmart,37C20min;sourcebluntfragmentation',S.DIGESTION,'SourceHpy166II;selectedcutsnotpredictedwithoutcallerreference')
 x=op(g,'dilute1',x,'SourceDNA concentrationcappedat10ng/ulusingpriorQubit;dilutionconditional',S.GENERIC)
 x=cleanup(g,'selection1',x,'Source75percentdilutedAMPurePB,0.8xratio,two80percentEtOHwashes,22ulTrispH8');x=assay(g,'postfrag_qubit',x,'SourcepostfragmentationQubitassay')
 x=op(g,'nick_ligation',x,'SourceEcoliDNAligase15U,NAD+,30ulrCutSmart16C30min,lidoff',S.NICK_SEALING,'SourceNADdependentEcoliDNAligase;exactnicks/ligationextentunknown')
 x=op(g,'dilute2',x,'Sourceconcentrationcappedat10ng/ulusingpreviousQubit;dilutionconditional',S.GENERIC);x=cleanup(g,'selection2',x,'Source0.75xdilutedAMPurePB,two80percentEtOHwashes,22ulTrispH8');x=assay(g,'postnick_qubit',x,'SourcepostnickcleanupQubitassay')
 used,rest=split(g,'klenow_input',x,('Source21ulDNAaliquotto30ulKlenowreaction','ResidualpostnickDNAeluate;exactremainingamountafterQCunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source22ulelutiondoesnotallenter21ulKlenowinput.'))
 x=op(g,'klenow',used,'SourceKlenowexo-minus7.5U,NEBuffer4,37C30min;'+('ddCTP/ddGTP/ddTTP,NOdATP' if mode=='blunt' else 'dATPplusddCTP/ddGTP/ddTTP'),S.EXTENSION,'SourcenuclotidesandKlenow;damage/mismatch/nativeoriginalDNAandallpossibleproductsretained;noadversarialrepairprediction')
 x=op(g,'dilute3',x,'Sourceconcentrationcappedat10ng/ulusingpriorQubit;dilutionconditional',S.GENERIC);x=cleanup(g,'selection3',x,'Source0.75xdilutedAMPurePB,22ulTrispH8;adjustto30ulNEBuffer4/NFW')
 ad=g.input('hairpin_adapter','SourcebcAd1001blunt5phosphateadapterstock;separate95C3min/RT30minannealingbeforecontactDNA' if mode=='blunt' else 'SourcePacBioBarcodedOverhangAdapter;exactwholeoligo/endinventoryunreviewed',K.ADAPTER,state=hairpin_scene() if mode=='blunt' else None,chemical_state=source_hairpin() if mode=='blunt' else None);e=g.input('ligation_stock','SourceSMRTbellExpressTemplatePrep2.0:2.5uladapter+15ulligationmix+0.5uladditive+0.5ulenhancer,20C60min,lidoff',K.CATALYST);y=g.product('ligated','Hairpin-ligated/linear/nickedDNA/freehairpinsandligationremnants;notassumedallSMRTbells',K.MOLECULAR);g.react('adapter_ligation','Sourcehairpinadapterligation',S.LIGATION,((x,R.TEMPLATE),(ad,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=reported_adapter_end('HiDEFbluntadapter',kind='blunt',statement='Publishedsourceblunthairpinand5phosphateinventory;fullselectedDNAnicks/ligationgeometryunassembled.') if mode=='blunt' else AdapterLigationGeometry('unavailable','unavailable',(unknown('Wholehistoricaloverhangadapter/endgeometryunreviewed'),),'Sourcebarcodedoverhangadapterknown;exactgeometrynotassembled.'))
 x=op(g,'nuclease_cleanup',y,'Immediatelyadd4ulSMRTbellEnzymeCleanUp1.0A/B/C/Dmix,37C60min;physicaldigestionofunsealed/nickedclasses',S.DIGESTION,'SourcePacBioEnzymesA2ul/B0.5ul/C0.5ul/D1ul;actualenzymeidentities/endproductsunreviewed')
 x=cleanup(g,'selection4',x,'Source1.2xdilutedAMPurePB,24ulTrispH8elution')
 if mode=='blunt':x=cleanup(g,'blunt_dimer_cleanup',x,'Additionalsource1.1xdilutedAMPurePBremovesresidualadapterdimers')
 x=assay(g,'final_qc',x,'SourceQubit/HighSensitivityD5000ScreenTapeassays;source<5percentmassbelow1kbandabove4kbQCcriterion,notautomaticphysicalcutoff')
 other=g.input('other_run_libraries','Independentsourcebarcodedlibrariesfortherun',K.LIBRARY,primary=True);pool=g.product('pooled','SourceQubitnormalizedrunpool≥50ngin≤15ul;concentrationmethodunassigned',K.MOLECULAR);g.react('pooling','Sourcephysicalrunmultiplexing',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
 if mode=='blunt':pool=op(g,'conditional_pool_cleanup',pool,'Sourceconditional1.3xdilutedAMPurePBifresidualdimersafterpooling;whetherperformedunassigned')
 return finish(g,pool,'SequelII/IIe,8MSMRTCell,Binding2.0/Sequencing2.0,primerv4,1hbinding/diffusionload125–160pM/2hpreextension/30hmovie. NativeDNAnoPCR;A-tailedandbluntbranchesnotmerged. ActualUdamagepreservedchemicallyunlessexplicitlytreated;C>TcallnotliteralU→Tconversion. Ancillarystockanneal95Cseparatefrom≤37ClibraryDNAprocessing.')
def protocol_graph():return graph_for_mode('blunt')
def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();c('Derived fold candidate is inferred by default',source_hairpin().inferred,True)
 for name,(seq,bc) in ADAPTERS.items():
  h=source_hairpin(name);c('Callerindependentexactselfpair '+name,h.left_stem,revcomp(h.right_stem));c('Wholehairpinsequenceconserved '+name,h.left_stem+h.loop+h.right_stem,seq)
 g=protocol_graph();a=graph_for_mode('a_tailed');c('NonAtailedspecificextraAMPurebranch',('blunt_dimer_cleanup' in g.steps,'blunt_dimer_cleanup' in a.steps),(True,False));c('SourceleftoverDNAexplicit','klenow_input_1' in g.terminals,True);c('Noamplificationinvented',all(st.kind!=S.PCR for st in g.steps.values()),True);c('CompleteSMRTbellgeometrynotclaimed',bool(g.publication_problems()),True)
 src=Source(sources_root()/'hidef-seq__10.1038+s41586-024-07532-8'/'s41586-024-07532-8_PMC11216816.xml.txt','Associatedmainhairpinsequences/phosphates')
 if have(c,src,label='SourceAM'):
  t=src.path.read_text();c('Allfourwholephosphorylatedsourcehairpinsverified',all('/5’Phos/'+seq in t for seq,_ in ADAPTERS.values()),True)
 c.report()
if __name__=='__main__':selftest()
