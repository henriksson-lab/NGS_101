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

from illumina import P5,P7,TRUSEQ_READ2,TRUSEQ_SMALL_RNA_RA5
SOURCE='Published2014PMC4179956main68–82/154–166 and actualpublishedSI FigS1–S4/oligopatterns;twooriginals. Selected optimized5primeBiotin TSO8, not all TSO blocking controls. RNA/DNA are branches ofone method.'
DT_BODY=TRUSEQ_READ2[13:]
TSO_BODY=TRUSEQ_SMALL_RNA_RA5
PCR_F=P7+'CGTGAT'+TRUSEQ_READ2
PCR_R=P5+TSO_BODY[:21]
CUSTOM_SEQ=TSO_BODY+'GGG'
def dt_primer():
 s=Scene();s.strand('one-base anchoredRTprimer',[Segment('source adapter',DT_BODY),Segment('source dT30','T'*30),Segment('sourceVanchor','V',placeholder=True)]);return s

def blocked_tso():
 s=Scene();s.strand('source TSO8',[Segment('source DNA adapter body',TSO_BODY),Segment('source RNA guanosines','GGG',note='Exact source3prime rGrGrG, not three DNA dGs')],mod5='Biotin');return s

def indexed_forward():
 s=Scene();s.strand('sourcePCRforward',[Segment('sourceP7',P7),Segment('source6ntsampleindex','CGTGAT',feature=feature('cats_sample_index','sample_index','fixed')),Segment('sourceR2body',TRUSEQ_READ2)]);return s

def graph_for_template(template,rt_aliquot=1):
 if template not in {'RNA','DNA'}:raise ValueError('Explicit source RNA/DNA branch required')
 if rt_aliquot not in {1,10}:raise ValueError('Source permits1ul or10ul RT aliquot; actual chosen route required')
 g=ProtocolGraph('CATS — source '+template+' selected optimized route',authority=SOURCE,representation='published_instance')
 x=g.input('template','Source short '+template+' input; selected1ng synthetic control branch, exact molecule preparation unassembled; native plasma/longfragment branches separate',K.MOLECULAR,primary=True)
 x=op(g,'tail',x,('Source E.coliPAP/1xPAPbuffer/RNaseinhibitor10U/0.1mMATP37C10min;65C20min inactivation' if template=='RNA' else 'SourceTdT/1xTdTbuffer/0.1mMdATP37C30min;70C10min inactivation'),S.TAILING,('SourceE.colipolyApolymerase' if template=='RNA' else 'Source terminaldeoxynucleotidetransferase'))
 primer=g.input('rt_primer','Source0.1uM anchoredpoly-dT primer for1nginput;5pguses0.001uM separately',K.OLIGO,state=dt_primer());enzyme=g.input('smartscribe','Source single SmartScribeRT addition1ul withFirstStrandbuffer,dNTP/DTT/RNaseinhibitor; subsequentlyretained duringTSOincubation',K.CATALYST);y=g.product('first_copy','Source first-copy DNA/inputtemplate/unextendedprimers/reagents; co-retainedSmartScribe tracked as separate same-tube component; variable terminal bases predominantlyC, not inventedexactCCC',K.MOLECULAR)
 retained=g.product('retained_smartscribe','Same-tube SmartScribe component retained afterfirst15min;notphysicalfractionornewenzymeaddition',K.CATALYST)
 g.react('first_synthesis','Source1ul tailedtemplate+2.5ul1xRTbuffercontaining20%DMSO+1ulprimer;72C2min/42C1min BEFOREenzyme master mix;42C15min',S.REVERSE_TRANSCRIPTION if template=='RNA' else S.EXTENSION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,retained,))
 tso=g.input('tso','Source later1ul10uM5primeBiotinTSO8; three3primeRNAguanosines; no5primeblockedbiotintemplate added',K.OLIGO,state=blocked_tso());x=g.product('ts_product','Source template-switched/unswitchedfirst-copy DNA plus originaltemplates, TSO/RTprimers, retainedenzyme and minoroffproducts; no perfecttail/TSyield guarantee',K.MOLECULAR)
 g.react('switch','Source addTSO AFTERfirst15min;another42C15min;70C10min inactivation',S.TEMPLATE_SWITCHING,((y,R.TEMPLATE),(tso,R.STRUCTURAL),(retained,R.CATALYST)),(x,))
 used,rest=split(g,'rt_transfer',x,('Source'+str(rt_aliquot)+'ulRTproduct to100ulPCR','Untransferred RTstock;actualremainingvolume/fateunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No implicit wholeRTconsumption; source1/10ul alternatives remain explicit.'))
 x=pcr(g,'pcr',used,'Source100ul2xQiagenTaq mastermix with250nM eachPCRprimer;actual sample-cycleprogramme not assigned in selectedmain',(('Source indexedforward',indexed_forward()),('Source reverse',oligo_scene('PCRreverse',PCR_R))))
 x=cleanup(g,'column',x,'SourceQIAquickPCRPurification column cleanup; small-empty/concatemer products not declared absent')
 x=cleanup(g,'gel',x,'Source preparative4%agarose/PureLinkGelExtraction fornext-generationsequencing;actualsizewindowunassembled')
 x=assay(g,'bioanalyzer',x,'Separate sourceAgilentHighSensitivityDNAQC aliquot')
 x=op(g,'denature_load',x,'Source5nMlibrary/0.2NNaOH5min then11pM forMiSeqloading;exact aliquot/workup not assembled')
 return finish(g,x,'SourceMiSeqReagentKit50cycles run77cycles; customDNAsequencingprimer '+CUSTOM_SEQ+' fromSI FigS4, not an arbitraryrepair ofstandardprimer. Exactindexcycles/readmask/junctiongeometry unassembled. FigS4 dualTSO8ACA methodseparate, not selected. Nativeplasma/Mg-fragmentedRNA PNKenhancement andDNA95C5min denaturation branchesnotautomaticallyappliedtosyntheticcontrol. NoUMI assigned; CGTGATfixedsampleindex distinctfromnontemplatedC/Gswitchingbases.')
def protocol_graph():return graph_for_template('RNA')
def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();r=protocol_graph();d=graph_for_template('DNA',10);c('RNAcopies useRT, DNAcopy usesextension', (r.steps['first_synthesis'].kind,d.steps['first_synthesis'].kind),(S.REVERSE_TRANSCRIPTION,S.EXTENSION));c('TSO followsfirstcopy',r.steps['switch'].inputs[0].material,'first_copy');c('SmartScribeexternalstockaddedonceandretainedforTS',sum(i.material=='smartscribe' for step in r.steps.values() for i in step.inputs)==1 and r.steps['switch'].inputs[2].material=='retained_smartscribe',True);c('SmallRTaliquotretainsstock','rt_transfer_1' in r.terminals,True);c('TSO8 retains5primeBiotin',blocked_tso().strands['source TSO8'].mod5,'Biotin');c('CompletePCR/readgeometryunassembled',bool(r.publication_problems()),True)
 src=Source(sources_root()/'cats__10.4161+rna.29304'/'rna.29304_PMC4179956.xml.txt','Publishedmain162wholePCRprimerbodies')
 if have(c,src,label='Publishedmain'):
  t=''.join(src.path.read_text().split());c('BothwholePCRprimerssourceexact',all(x in t for x in (PCR_F,PCR_R)),True)
 si=Source(sources_root()/'cats__10.4161+rna.29304'/'rna.29304_supp_rna-11-817-s01.pdf.txt','PublishedSI TSO8/FigS4')
 if have(c,si,label='PublishedSI'):
  t=''.join(si.path.read_text().split());c('WholeBiotinTSOwith3RNAguanosinessourceexact','Biotin-'+TSO_BODY+'rGrGrG' in t,True);c('Customprimerliteralverified',CUSTOM_SEQ.lower() in t.lower(),True)
 c.report()
if __name__=='__main__':selftest()
