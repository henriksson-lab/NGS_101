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

from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
SOURCE='Version-qualified PMC6939866 author manuscriptMethods289–295/sourceoverview85;premadeengineeredstocksunassembled. Selectedcustom-primer screen route,not8staggerhumanlibraryroute.'
PCR1F='TCCAGTCAGGTGTGATGCTCGGGGATCCAGGAATTCAGTTTGTACAAAAAAGCAGGCTCA'
PCR1R='CGAGCTTATCGTCGTCATCCCCACTTTGTACAAGAAAGCTGGGTCA'
PCR2F=TRUSEQ_READ1[:20]+'TCCAGTCAGGTGTGATGCTC'
PCR2R=TRUSEQ_READ2[:-1]+'CGAGCTTATCGTCGTCATCC'
PCR3F=P5+TRUSEQ_READ1[4:20]+'TCCAGT'
SEQ_PRIMER='GGTGTGATGCTCGGGGATCCAGGAATTC'
SOURCE_OLIGO_EXCEPTIONS={'pcr3_p1':'The original prints seven literal lowercase x positions between the complete P7 and R2-prefix bodies. The generic scanner splits that notation; the local source check verifies the complete P7+xxxxxxx+R2-prefix string and exact seven-position count. This exemption does not supply a concrete sample index.'}
def indexed_reverse():
 s=Scene();s.strand('PCR3reverse',[Segment('sourceP7',P7),Segment('source7baseSAMPLEindex','N'*7,placeholder=True,feature=feature('tscan_sample_index','sample_index','fixed')),Segment('sourceR2prefix',TRUSEQ_READ2[:21])]);return s

def protocol_graph():
 g=ProtocolGraph('T-Scan — source sorted-cell antigen-cassette NGS preparation',authority=SOURCE,representation='published_instance')
 edc=g.input('epitope_discovery_cells','Sourcepremadeantigen-libraryEDCswithfluorogenicGzBreporter;upstreamengineeredstockconstructionunassembled',K.BIOLOGICAL,primary=True);tcells=g.input('effector_cells','IndependentsourceTcellpreparation;actualTCR/sampleselectionandupstreamconstructionunassembled',K.BIOLOGICAL,primary=True);culture=g.product('coculture_product','SourceEDC/Tcellmixtureandallunreactive/reactivecellstates;notperfectrecognitionpredicted',K.MOLECULAR);g.react('coculture','SourceTScanfunctionalcellcoculture;conditionsnotassignedinthisselectedNGSscope',S.POOL,((edc,R.TEMPLATE),(tcells,R.TEMPLATE)),(culture,))
 x=cleanup(g,'facs',culture,'SourceGzBreporterpositiveEDCsselectedbyFACS;uncollectedcellfractionsphysical,gatepuritynotpredicted')
 carrier=g.input('carrier_cells','Source2e5wildtypeHEK293Tcarriercells;genomicDNAisnotmoleculeUMI',K.BIOLOGICAL,primary=True);mix=g.product('carrier_mix','SourceFACSrecoverypluscarriercells/carryover',K.MOLECULAR);g.react('carrier_pool','SourcecarriercelladditionbeforegDNAextraction',S.POOL,((x,R.TEMPLATE),(carrier,R.TEMPLATE)),(mix,))
 kit=g.input('genejet_stock','SourceGeneJETgDNApurificationkit;historicalinternalchemistryunreviewed',K.CONDITION);x=g.product('gdna_product','SourcepurifiedgDNAfromsortedcellspluscarrier/carryover;recoveryunquantified',K.MOLECULAR);w=g.product('genejet_removed','Physicalkitremovedfractions;compositionunassembled',K.WASTE);g.react('gdna','SourceGeneJETgDNAisolation',S.MAGIC,((mix,R.TEMPLATE),(kit,R.CONDITION)),(x,w),reason='KnownsourcegDNApurificationwithphysicalseparation;exacthistoricalkitlysis/washchemistryunassembled.');g.discard(w,'Actualkitpurificationfractions;noexactspeciespurityprediction.')
 x=pcr(g,'pcr1',x,'SourcecollectionofQ5hotstartPCR1reactions:1uggDNAper100ul;reactioncount/cyclesunassigned;printedvendorThermo retained',(('SourceTCellPCR1F',oligo_scene('PCR1F',PCR1F)),('SourceTCellPCR1R',oligo_scene('PCR1R',PCR1R))))
 pooled=g.product('pcr1_pool','SourcephysicallypooledPCR1reactionsfromonescreenreplicate;actualreactioncountunassembled',K.MOLECULAR);g.react('pool_pcr1','ActualsourcePCR1replicate-reactionpool',S.GENERIC,((x,R.TEMPLATE),),(pooled,),changes_molecule=False)
 used,rest=split(g,'pcr2_transfer',pooled,('Source1ulPCR1poolto50ulPCR2','UntransferredPCR1pool;volume/fateunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'NoimplicitwholePCR1poolconsumption.'))
 x=pcr(g,'pcr2',used,'Source50ulQ5PCR2;cycleprogrammeunassigned',(('SourceTCellPCR2F',oligo_scene('PCR2F',PCR2F)),('SourceTCellPCR2R',oligo_scene('PCR2R',PCR2R))))
 used,rest=split(g,'pcr3_transfer',x,('Source1ulPCR2to50ulPCR3','UntransferredPCR2stock;fateunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'NoinventedPCR2cleanuporwhole-stockconsumption.'))
 x=pcr(g,'pcr3',used,'Source50ulQ5PCR3addssequencingadaptersand7ntsampleindex;actualindex/cyclesunassigned',(('SourceTCellPCR3F',oligo_scene('PCR3F',PCR3F)),('SourceTCellPCR3Rsampleindexpattern',indexed_reverse())))
 other=g.input('other_samples','IndependentsourceindexedPCR3samplepreparations;actualpoolmembersunassigned',K.MOLECULAR,primary=True);pooled=g.product('sample_pool','Actualsourcephysicalsamplepool',K.MOLECULAR);g.react('pool_samples','SourcePCR3samplespooledbeforegel',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,));x=cleanup(g,'gel',pooled,'Sourcegel-extractionwithQIAquickGelExtractioncolumns;sizewindowunassigned')
 return finish(g,x,'SourceMiSeqorNextSeqwithcustomT7-Pep2.2_Illumina_SP;exactcycle/run/indexprimerbindinggeometryunassembled.7ntPCRsampleindexisnotUMI. Unsorted≥100xrepresentationinputcontrolisseparatepreparation,notinventedmolecularmixwithselectedsample. Human-genome8staggerprimerrouteandstandardsequencingprimerremainseparatecase;staggernotUMI. Readenrichmentcomputational.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();c('CarrierentersbeforegDNA',g.steps['gdna'].inputs[0].material,'carrier_mix');c('BothsmallPCRtransfersretainstock',all(k in g.terminals for k in ('pcr2_transfer_1','pcr3_transfer_1')),True);c('IndexedreversefeaturemarksSAMPLE',indexed_reverse().strands['PCR3reverse'].segs[1].feature.role,'sample_index');c('NoassembledwholePCRorreadclaim',bool(g.publication_problems()),True)
 src=Source(sources_root()/'t-scan__10.1016+j.cell.2019.07.009'/'j.cell.2019.07.009_PMC6939866.xml.txt','AssociatedmainMethods293customprimerroute')
 if have(c,src,label='SourceAM'):
  t=''.join(src.path.read_text().split());c('WholefiveexactPCRprimerbodiesverified',all(s in t for s in (PCR1F,PCR1R,PCR2F,PCR2R,PCR3F)),True);c('Wholeunknown7ntsourcePCR3patternpreserved',P7+'xxxxxxx'+TRUSEQ_READ2[:21] in t,True);c('Wholecustomsequencingprimerverified',SEQ_PRIMER in t,True)
 c.report()
if __name__=='__main__':selftest()
