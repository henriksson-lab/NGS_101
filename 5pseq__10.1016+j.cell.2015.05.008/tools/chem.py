from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

AUTHORITY='One already selected associatedAMNIHMS688997/PMC4461875 main51/135–141/Fig1Acaption385. Original-source pass closed by targeted owner before implementationhandoff; no final-equivalence or HT5PSeq/PARE recipe inheritance.'

def native_rna_scene():
    s=Scene();s.strand('rna',[Segment('Source native5P RNA subset within totalRNA mixture','[RNA body unassembled]',placeholder=True)],mod5='P');return s

def source_oligo_scene():
    s=Scene();s.strand('source_oligo',[Segment('Source DNA/RNA molecular-barcode oligo; whole order/sugar partition unassembled','[DNA/RNA oligo unassembled]',placeholder=True)]);return s

def random_primer_scene():
    s=Scene();s.strand('primer',[Segment('Source random-hexamer priming preparation; full modifications unassembled','[random-hexamer primer]',placeholder=True)]);return s

def protocol_graph():
    g=ProtocolGraph('5PSeq nativeRNA5P fragment/randomRT/capture preparation',authority=AUTHORITY,representation='published_instance')
    start=g.input('source','Actual50–100ugDNA-freeTOTALRNA mixture; native5P subset displayed, other5endclasses/populationsunassembled',K.MOLECULAR,primary=True,state=native_rna_scene(),chemical_state={'source_native_five_prime':'P on ligatable nativeRNAclass, notalltotalRNA','input_history':'DNAfree totalRNA already prepared; extraction/DNase history unassembled'})
    b=B.__new__(B);b.g=g;b.x=start
    oligo=g.input('barcode_oligo','ActualsourceDNA/RNA molecular-barcode ligation oligo; full ordering/length/sugar/biotinlocationunassembled',K.OLIGO,state=source_oligo_scene(),chemical_state={'molecular_barcode':feature('fivepseq.molecule','umi','random',note='SourceFig1molecularbarcode/main141randombarcodecollapse; length/coordinate/readspanunassembled'),'native_rna_end':'Source5P substrate','adapter_acceptor_3prime':'ConditionalOH requiredbylocalRNA5P ligation interpretation, not directlyassignedsourcewholeoligo','biotin_location':'Unassigned; later streptavidincapture doesnotaloneidentify carrierstrand'})
    geometry=AdapterLigationGeometry('unavailable','unavailable',(native_rna_scene(),source_oligo_scene()),'Source5P RNA subset andDNA/RNAbarcodedoligo positive. Exactadapterstrand/body/ends/sugars/junctionunassembled;3OHacceptor isconditionalchemistryinterpretationonly, notcompletegeometry.')
    b.op('native_ligation','Source direct native5primeP ligation, NO CIP/TAP inproduction',S.MAGIC,((oligo,R.STRUCTURAL),b.stock('ligation_preparation','Actual source RNA5P ligation preparation; ligase/cosubstrates/buffer/dosesunassembled',K.MOLECULAR,R.STRUCTURAL)),geometry=geometry,reason='Source direct native5P ligation positive; completeoccupiedRNA/adapterends/sugarpartition/junctions/offproductclassesunassembled. No PARE/HT5PSeq adapter or vendor recipe inherited.',product='Source ligated/unligated native5P and otherendRNA/oligo classes; no perfectligation or automaticexclusionofcappedRNA claim')
    enriched=g.product('polya_enriched','Actual selected sourcepolyA-enrichedRNA preparation afternativeendligation',K.MOLECULAR);other=g.product('polya_other','Actual otherRNA/enrichmentworkupfractions, fatesunassigned',K.MOLECULAR)
    g.react('poly_a_enrich','Source primarypolyAenrichment branch; rRNAdepletionalternative notmerged',S.MAGIC,((b.x,R.TEMPLATE),b.stock('enrichment_preparation','Actual source polyAenrichment preparation; capturearm/wholeworkup formulationunassembled',K.MOLECULAR,R.STRUCTURAL)),(enriched,other),reason='Positive sourcepolyAenrichment afterRNA5Pligation. Exactcapturetechnology/reagents/endclasses/fractionsunassembled; no invented oligo-dT manufacturer recipe.');g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Actualnonselectedenrichmentfraction notimplicitlysequenced.'));b.x=enriched
    b.op('fragment','Source RNA fragmentation AFTERendligation/polyAenrichment',S.MAGIC,(b.stock('fragment_preparation','ActualsourceRNAfragmentationpreparation, reagent/time/dose/endchemistryunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Source RNA fragmentation positive, no borrowedchemicalhydrolysis/RNaseI or exactfragmentlength/enddistribution.')
    primer=g.input('random_hexamers','Actualsourceexternalrandom-hexamerRTpriming preparation; NOTadditionalmoleculeUMI',K.OLIGO,state=random_primer_scene(),chemical_state={'priming_class':'random hexamer','identifier_meaning':'none; sourceUMI belongs to native5ligationoligo, nottheseprimers'})
    b.op('random_rt','Source random-hexamer reverse transcription ofRNAfragments',S.MAGIC,((primer,R.STRUCTURAL),b.stock('rt_preparation','Actual source RT reaction preparation; enzyme/nucleotides/buffer/doseunassembled',K.MOLECULAR,R.STRUCTURAL)),reason='Positive randomhexamerRT; completeoccupiedRNA/primerproducts/unextendedclassesunassembled. No assignment of biotin to RTprimer merelyfromdownstreamcapture.',product='Source cDNA/RNA/ligatedoligo preparations; intendedcDNAclasses retain5prime-ligatedoligoinformation, exactcovalentbiotincarrierunassigned')
    bound=g.product('streptavidin_bound','Actualsource retainedcDNA preparation containing5prime-ligatedoligo; capturecompetencepositive but exactbiotincarrier/occupancyunassembled',K.MOLECULAR);flow=g.product('streptavidin_other','Actualnonretainedcapturefraction, fateunassigned',K.MOLECULAR)
    g.react('streptavidin_capture','Source streptavidincapture ofcDNAcontainingnative5prime-ligatedoligo',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),b.stock('streptavidin_preparation','Actualsource streptavidinselectionpreparation; support/formulation/doseunassembled',K.SELECTION,R.SELECTION)),(bound,flow));g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Source nonretainedcapturefraction notreturnedtoproduction byassumption.'));b.x=bound
    b.op('library_finish','Source Illumina libraryproduction fromcapturedpreparation',S.MAGIC,(b.stock('finish_preparation','Actualsource libraryfinishingreactionstock; PCR/primer/adapters/index/NT/internalfractionationunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Positive sourceIlluminasequencing endpoint; no kit, readmode, PCRprogramme or wholeadapter recipe invented fromplatformname.')
    b.split('loading','ActualsourceIlluminasequencingloadingallocation, amountunassigned','Selectednative5PSeqloadingaliquot','Untransferredfinishednative5PSeqlibrarystock, amount/fateunassigned')
    lib=g.product('native_library','Source native5PSeqsequencinglibrary, fullmolecules/readgeometryunassembled',K.LIBRARY);g.react('library_load','Source Illumina acquisition; actualloadingallocationunassembled',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='Native5P/molecularbarcode ligation followedbyRNAfragmentation/randomRT/capture; exactprimerlanding/readmode/cycles/index/UMIlength/readspanunassembled. Main141barcodecollapse and17ntstopoffset/3ntperiodicity arecomputational, notnewchemicaltreatments or28ntRPFproduction.'))
    for name,cap_only in (('cap_control',True),('composite_control',False)):
        b.x=g.input(name+'_rna','Actualseparate preparedRNAfor'+('5caponly'if cap_only else'5capplusnative5P')+' sourcecontrol; commonsamplestock/allocationrelationshipunassembled',K.MOLECULAR,primary=True)
        if cap_only:b.op(name+'_cip','SourceCIP treatment removesexisting5primeP beforeTAP',S.MAGIC,(b.stock(name+'_cip_stock','ActualsourceCIPenzyme; fullrecipe/doseunassembled',K.CATALYST,R.CATALYST),),reason='Separate source5capcontrol; exactoccupiedRNAendclasses/productyieldunassembled, notCIPadditiontonativeproduction.')
        b.op(name+'_tap','SourceTAP treatment '+('afterCIP for5capcontrol'if cap_only else'withoutCIP forcomposite5capplusnative5Pcontrol'),S.MAGIC,(b.stock(name+'_tap_stock','ActualsourceTAPenzyme; fullrecipe/doseunassembled',K.CATALYST,R.CATALYST),),reason='Actual sourcecontrol branch; completecap/5Pproducts/internalcleanup/libraryworkupunassembled. No controlRNA silentlymergedintonativeproduction.')
        g.mark_terminal(Terminal(b.x.id,Disposition.UNRESOLVED,'Separate sourcecontrolpreparation; actualdownstream5endligation/libraryfinishing notinstantiatedhere.'))
    g.add_source_item(SourceItem('identifier_and_fate_scope','AssociatedAMFig1Acaption385DNA/RNAoligowithmolecularbarcode andmain141randombarcodecollapse establishmoleculeUMImeaning only.','Wholebarcodebody/length/sugar/readspan andbiotinlocationunassembled. ActualSI1–7rolesgeneric/unverified; sourceDOCX6isnotassertedtobeExtendedProcedures. Noabsenceclaim.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('NativeinputretainsactualtotalRNAmixture',g.steps['native_ligation'].inputs[0].material,'source');c('Productionhasnocapenzymes',all('control'not in x.material for x in g.steps['native_ligation'].inputs),True);c('PolyAafterdirectnativeendligation',g.steps['poly_a_enrich'].inputs[0].material,'native_ligation_product');c('RNAfragmentationafterenrichment',g.steps['fragment'].inputs[0].material,'polya_enriched');c('RandomRTusesactualfragmentedRNA',g.steps['random_rt'].inputs[0].material,'fragment_product');c('CaptureusesRTpreparation',g.steps['streptavidin_capture'].inputs[0].material,'random_rt_product');c('Libraryusesretainedcapturefraction',g.steps['library_finish'].inputs[0].material,'streptavidin_bound');c('CapcontrolTAPfollowsCIP',g.steps['cap_control_tap'].inputs[0].material,'cap_control_cip_product');c('CompositeTAPbypassesCIP',g.steps['composite_control_tap'].inputs[0].material,'composite_control_rna');c('Sequencingusesselectedloadingaliquot',g.steps['library_load'].inputs[0].material,'loading_selected');c('Onlynativebranchsequencinginstantiated',set(g.sequencing),{'native_library'});c('Partialnotpublicready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
