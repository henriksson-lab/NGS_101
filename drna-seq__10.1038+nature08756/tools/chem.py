"""Published2010dRNA-seq selected ML mock/TEX pair; private incomplete molecular graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from chemdraw import Scene,Segment,feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
BARCODES={'minus':'ACTA','plus':'AGGT'}
SOURCE='Published2010mainidentity/abstractsubscriptionpreview; actualpublishedSI206Methods237–282 andSI208Table2originalp6. Threeoriginals; sourceinvestigationclosed. Corepreparationsourcepositive, fullmainMethods/vendorrecipe notread.'

def unknown(label):
 s=Scene();s.strand('unassembled',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def indexed_primer(arm):
 if arm not in BARCODES:raise ValueError('ExplicitMLminus/plusarmrequired')
 s=Scene();s.strand('source_forward',[Segment('source_5_context','[unassembled5primecontext]',placeholder=True),Segment('source_barcode',BARCODES[arm],feature=feature('drna_'+arm,'sample_index','fixed')),Segment('source_binding','[unassembledbindingarm]',placeholder=True)]);return s

def source_rna_end():
 s=Scene();s.strand('rna',[Segment('source_RNA','[sourceRNAsequenceunassembled]',placeholder=True),Segment('source_polyA','[variablepolyAtail]',placeholder=True)],mod5='P');s.note('rna','Source5P ligation-compatibleRNAclass only; other5endclassesandfailedtailingremainco-present, exactadapterend/junctionunassembled.');return s

def protocol_graph():
 g=ProtocolGraph('dRNA-seq selected2010 MLminus/pluspair',authority=SOURCE,representation='published_instance')
 x=g.input('total_rna','SourceextractedtotalRNAforoneMLmock/TEXpair; upstreamspecimenculture/extractionnotmodeled',K.MOLECULAR,primary=True)
 dn=g.input('dnase_enzyme','SourceDNaseIremovesresidualgenomicDNA;fullconditionsunassigned',K.CATALYST);y=g.product('dnase_product','SourceDNase-treatedRNAplusDNAcleavageproducts/reagents; no perfectDNAremovalclaim',K.MOLECULAR);g.react('dnase','SourceDNaseIbeforeequalRNAallocation',S.DIGESTION,((x,R.TEMPLATE),(dn,R.CATALYST)),(y,))
 arms={a:g.product(a+'_rna','SourceequalRNAamountselectedforML'+a+'arm; actualmassunassigned',K.MOLECULAR)for a in BARCODES};rest=g.product('rna_remainder','UntransferredRNAmaterial;amountunassigned',K.MOLECULAR)
 g.react('pair_allocation','SourceequalamountsforTEXandbuffercontrol',S.SPLIT,((y,R.TEMPLATE),),(*arms.values(),rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No whole extractedRNAstockconsumptionor perfectrecoveryasserted.'))
 for arm,start in arms.items():
  def key(k):return arm+'_'+k.lower()
  def op(k,label,kind,reagent,product,role=R.CATALYST,extra=(),geometry=None):
   nonlocal start
   c=g.input(key(k+'_stock'),reagent,K.CATALYST if role==R.CATALYST else K.CONDITION);y=g.product(key(k+'_product'),product,K.MOLECULAR);kw={'geometry':geometry}if geometry is not None else {};g.react(key(k),label,kind,((start,R.TEMPLATE),(c,role),*extra),(y,),**kw);start=y
  def clean(k,label):
   nonlocal start
   c=g.input(key(k+'_condition'),label,K.CONDITION);kept=g.product(key(k+'_kept'),'SourcephysicallyrecoveredRNApreparation;yield/carryoverunassigned',K.MOLECULAR);w=g.product(key(k+'_waste'),'Sourcephysicalorganic/precipitationworkupfractions',K.WASTE);g.react(key(k),label,S.CLEANUP,((start,R.TEMPLATE),(c,R.CONDITION)),(kept,w));g.discard(w,'Physicalworkup,notcomputationalfilteringorcompletechemicalselection.');start=kept
  if arm=='plus':op('tex','SourceTEX1Uper1ugRNA30C60min',S.DIGESTION,'SourceTerminator5primephosphate-dependentexonucleaseEpicentreTER51020','RNApopulationafterTEX;primary/resistantRNAplusresidualprocessedRNAanddegradationproductsremainuntilphysicalworkup. No absoluteselectivityorallrRNAdepletion.')
  else:op('mock','Sourcebuffer-alone30C60minmatchedmock',S.THERMAL_TREATMENT,'Sourcebufferalone','MatchedmockRNAmixture;notTEX-treated',role=R.CONDITION)
  clean('organic1','Sourcephenol/chloroform/isoamylalcohol25:24:1')
  clean('precipitate1','Sourceovernight2.5volumesethanol/0.1MsodiumacetatepH6.5precipitation;wholeformulationliteral,notaddedfragmentation')
  op('tap','SourceTAP1U37C1h BOTHarms afterTEX/mockworkup',S.CHEMICAL_CONVERSION,'SourceTobaccoAcidPyrophosphataseEpicentreT19100','Source5monophosphateRNAclassesfor5linkerligation,pre-existingmonoP/otherends/incompleteconversionremain. Not5P→5OHdephosphorylation.')
  clean('organic2','Sourceagainphenol/chloroform/isoamylalcohol25:24:1')
  clean('precipitate2','Sourceagainovernightethanol/acetateprecipitationasabove')
  op('polyA','Sourcepoly(A)polymeraseRNA3tailing;exactATP/conditions/taildistributionunassembled',S.TAILING,'Sourcepoly(A)polymerase;cofactor/nucleotideinventoryunassembled','SourcepolyadenylatedanduntailedRNAclasses;noexacttail-lengthorperfectactivityclaim')
  adapter=g.input(key('rna_adapter'),'SourceRNA5adapter;wholeordering/modificationsunassembled,historicalvertisrecipeunread',K.ADAPTER)
  geom=AdapterLigationGeometry('unavailable','unavailable',(source_rna_end(),),'SourceRNA5P substrateknown afterTAP/polyAtailing; actualRNAadapterendgeometry/bodyandsealedjunctionnotassembled. This is not a publicationabsenceclaim.')
  op('ligate','SourceRNAadapterligationto5PRNA AFTERpolyAtailing',S.LIGATION,'Sourcevendorligationcatalystidentity/wholeconditionsunassembled','Source5adapter-ligatedandunligatedRNAmix;no assumedperfectselection',extra=((adapter,R.STRUCTURAL),),geometry=geom)
  rt=g.input(key('rt_primer'),'Sourceoligo(dT)-adapterRTprimer; wholeordering/tailbindinglengthunassembled',K.OLIGO)
  op('rt','SourceMMLV-RNaseHminus42C20min,ramp55C,55C5min',S.REVERSE_TRANSCRIPTION,'SourceM-MLV-RNaseHminusreverse transcriptase','SourcecDNA/RNA/freeprimerreactionmixture; exactstrandretentionandRTjunctionunassembled',extra=((rt,R.STRUCTURAL),))
  f=g.input(key('pcr_fwd'),'SourcePCRprimercomponentwithknownML'+arm+'4ntbarcode,wholeorderingunassembled',K.OLIGO,state=indexed_primer(arm));r=g.input(key('pcr_rev'),'SourceotherPCRprimer,wholeorderingunassembled',K.OLIGO)
  pgeom=PCRGeometry.unavailable(unknown('SourcewholecDNA/primerlanding/amplicongeometryunassembled'),reason='Published4ntsampletagdoesnotestablishwholeoperationalprimer/templategeometry.')
  op('pcr','SourcehighfidelityPCRadds4ntbarcodeto5cDNAend,target20–30ngul;exactcyclecountunassigned',S.PCR,'SourcehighfidelityDNApolymeraseidentityunassembled','SourcebarcodedcDNAampliconsplusoriginalRNA/cDNA/freeprimer/reagents; wholepurification/offproductsunassembled',extra=((f,R.STRUCTURAL),(r,R.STRUCTURAL)),geometry=pgeom)
  lib=g.product(key('library'),'SourceML'+arm+'454sequencinglibrary;exactends/readgeometryunassembled',K.LIBRARY);g.react(key('library_identity'),'Sourceprepared454libraryidentity',S.GENERIC,((start,R.TEMPLATE),),(lib,),changes_molecule=False)
  g.sequence(lib,SequencingPlan(limitation='SourceRoche454FLX,ML'+arm+'barcode'+BARCODES[arm]+'in5endlinker; variablepyrosequencingreadlength,exactfulllinker/primer/rungeometryunassembled. NoIlluminaindexcyclesinvented. SeparatefragmentedSolexaGAII76SEoperonlibraryNOTthisTEXpairedgraph.'))
 g.add_source_item(SourceItem('pending','FullpublishedmainanddirectlycitedhistoricalvertismicroRNArecipe/oligos','SourcepositivepublishedSIcorepreparationonly,notfulllinked-source/absenceaudit.'));return g

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph()
 c('Bothactualsourcearmssequenced',sum(t.disposition==Disposition.SEQUENCING for t in g.terminals.values()),2)
 c('TAPfollowsTEXandphysicalworkup',g.steps['plus_tap'].inputs[0].material,'plus_precipitate1_kept')
 c('MockisnotTEXdigestion',g.steps['minus_mock'].kind,S.THERMAL_TREATMENT)
 c('Noadditional454RNAfragmentation',any(s.kind==S.FRAGMENTATION for s in g.steps.values()),False)
 c('PolyAprecedes5adapterligation',g.steps['plus_ligate'].inputs[0].material,'plus_polya_product')
 c('Whole454pairremainsprivatepartial',bool(g.publication_problems()),True)
 src=Source(sources_root()/'drna-seq__10.1038+nature08756'/'nature08756_41586_2010_BFnature08756_MOESM208_ESM.pdf.txt','PublishedTable2originalp6')
 if have(c,src,label='PublishedTable2'):
  import re
  t=src.path.read_text();c('IndependentactualMLbarcodes',all(re.search(r'ML\s*'+('-'if a=='minus'else r'\+')+r'\s*'+b,t)for a,b in BARCODES.items()),True)
 src=Source(sources_root()/'drna-seq__10.1038+nature08756'/'nature08756_41586_2010_BFnature08756_MOESM206_ESM.pdf.txt','PublishedSI206processedRNA/libraryparagraphs237–282')
 if have(c,src,label='PublishedMethods'):
  t=''.join(src.path.read_text().split());c('SourceTAPmakesmonophosphatebeforelibrary','togenerate5’-mono-phosphatesforlinkerligation' in t,True)
  c('Source454omitsRNAfractionation','omittingtheRNAsizefractionationstep'in t,True)
 c.report()
if __name__=='__main__':selftest()
