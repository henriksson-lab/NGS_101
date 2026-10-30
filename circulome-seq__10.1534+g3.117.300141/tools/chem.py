"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    sc=Scene();sc.strand('preparation',[Segment('population','['+label+']',placeholder=True)]);return sc

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')



def qc(g,key,x,label):
    keep=g.product(key+'_retained','Sourceproductionremainder; analyticalaliquotamountnotinvented',K.MOLECULAR);test=g.product(key+'_assay','Separateanalyticalaliquotforthe sourceassay',K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,test));g.mark_terminal(Terminal(test.id,Disposition.QC,'Analyticalmaterialdoesnotbecomethepreparativesequencingtemplate.'));return keep

def unknown_pcr(g,key,x,label,reagent):
    f=g.input(key+'_f','Sourceforwardindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);r=g.input(key+'_r','Sourcereverseindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);cat=g.input(key+'_enzyme',reagent,K.CATALYST);y=g.product(key+'_product','SourcePCRcopiesplusoriginaltemplates/unextendedprimersandreactionremnants; exactbackbonesunassembled',K.MOLECULAR)
    g.react(key,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('Sourcefulladapter/primer/endgeometryunassembled'),reason='KnownPCRoperationdoesnotestablishcompletefirstcopyorfinalstrandgeometry.'));return y

SOURCE='PublishedPMC5633380 Methods90–106/Fig1; originalPPTXSupplementFigureS1/S2slidecaptions; DOI10.1534/g3.117.300141. AnalyticalnativeeccDNA, noRCA.'
def selected_single_insertion_filled_insert(circle,cut):
    # Caller-selectedsingle-eventclassafterpolymerasegapfill, notanassayprediction.
    # This isinsert-only: fullTn5adapter/junction/readprimergeometryunassembled.
    if not circle or set(circle)-set('ACGT') or len(circle)<9 or type(cut)is not int or not 0<=cut<len(circle):raise ValueError('explicitcircularDNAandcutwith9bpwindowrequired')
    rotated=circle[cut:]+circle[:cut];return rotated+rotated[:9]

def graph_for_enrichment(enrichment='exo',attenuated=True):
    if enrichment not in {'exo','gradient','both'} or type(attenuated)is not bool:raise ValueError('sourceenrichmentandtagmentationbranchrequired')
    g=ProtocolGraph('Circulome-Seq — '+enrichment+('/attenuated' if attenuated else '/standardlowinput'),authority=SOURCE,representation='published_instance')
    x=g.input('genomic','SourceextractedtotalgenomicDNAwithnativecircles,linearDNAandothermaterial; controlDNApreparationseparate',K.MOLECULAR,primary=True)
    if enrichment in {'exo','both'}:
        x=op(g,'exov',x,'SourceC.elegansExoV400U/µgDNAover3days;heat70°C30minafterrounds',S.DIGESTION,reagent='SourceExoV/ATP/DTT/ampicillin; human5day/500–1000Uµgbranchseparate',catalyst=True,product='Nativecircle-enrichedpopulationplusdigestedlinearDNAresidues/possiblysurvivinglinearDNA andenzymes; digestionnotphysicalcleanup')
    if enrichment in {'gradient','both'}:
        condition=g.input('gradient_condition','SourceCsCl1.55g/ml/ethidiumgradient500237g2.5h; separateplasmidreference tube',K.SELECTION);bottom=g.product('bottom_fraction','Sourcecollectedapproximateclosedcirclefraction; noabsolutecomposition/purityclaim',K.MOLECULAR);top=g.product('other_fractions','Actualothergradientfractionsincludinglinearandnicked/relaxedcircles; notallcirclescollected',K.MOLECULAR)
        g.react('gradient','Sourcebiophysicaltopologyfractionation',S.CLEANUP,((x,R.TEMPLATE),(condition,R.SELECTION)),(bottom,top));g.mark_terminal(Terminal(top.id,Disposition.UNRESOLVED,'Otherfractionsincludesrelaxed/nickedcircles; notforcedintopurelinearDNAwaste.'));x=bottom
        x=cleanup(g,'butanol',x,'SourceCsCl/TEsaturated1butanolethidiumextraction')
        x=cleanup(g,'dialyze',x,'SourceTEdialysis2days4°C')
    x=op(g,'entry_boundary',x,'SourceeccDNApreparedfortagmentation; exactExoVbranchtransfer/cleanupunreviewed',reagent='Do notinventphysicalremovalofExoVresiduesorRCAstep',product='SourceenrichedDNApreparationwithunassembledcarryoverfates')
    tn5=g.input('nextera','SourceNexteraXTtagmentasecomplex; originaladapter/oligolengthsnotreviewed',K.MOLECULAR);enzyme=g.input('tagmentase','SourceTn5 catalyticcomponent',K.CATALYST);y=g.product('tagged','SourceTn5-fragmentedDNA/adaptersplusunreactedcircles,linearDNA andsourcecarryover; targetduplicationnotcompleteuntilgapfill',K.MOLECULAR)
    g.react('tagment','Source0.5µlXT/'+('37°C3minwithoutshakingattenuated' if attenuated else '37°C30minlowinput')+'; notsource1ng/1.5µlgenomiccontrolrecipe',S.TAGMENTATION,((x,R.TEMPLATE),(tn5,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,))
    x=op(g,'gapfill',y,'SourceFig1BfirstPCRpolymerasecyclefillsTn5gaps',S.EXTENSION,reagent='SourcehistoricalXTpolymerase; fullpriminggeometryunreviewed',catalyst=True,product='Conditionalsingle-insertioncircleinsertwithduplicated9bpboundaryplusmultipleinsertions/othertemplatesandreagents; nativecircleisnotidenticaltolinearcopy')
    x=unknown_pcr(g,'library_pcr',x,'SourcePCRcyclechoicebasedongels;10–12sufficientreported, notuniversalfixedcount','SourcehistoricalNexteraXTlibraryPCRenzymes/primersunreviewed')
    x=qc(g,'analytical_gel',x,'Sourcelibraryassessmentaliquotandretainedsequencingpreparation')
    g.sequence(x,SequencingPlan(limitation='Actualpaired-endreadmapping/9bpduplicationanalysisestablishedbyFig1andMethods; preciseplatform/run/adapter/readprimergeometryunreviewed. Computational9bp-circleevidenceisnotaperfectselectionorphysicaldiscardrule; sourceuntreatedgenomiccontrolsseparate.'))
    g.add_source_item(SourceItem('pending','HistoricalNexteraXToligo/primer/kitoriginalandexovcleanupgeometry','PublishedsourceS1/S2supplementcaptionsread; otherdataZIPinventoriesnotchemicaldisclosure. NoRCAorCircleSeq/CIDERfamilyrecipeimported.'))
    return g

def protocol_graph():return graph_for_enrichment()
def selftest():
    from checks import Check
    from chemdraw import revcomp
    c=Check();circle='ACGTGCAATTGCCGTA';product=selected_single_insertion_filled_insert(circle,3)
    c('selectedone-insertionfilledinsertretainscircleonceplus9bp',len(product),len(circle)+9)
    c('duplicationconnectsbothinsertboundaries',product[:9],product[-9:])
    c('oppositeorientedreadstartiscomplementaryboundary',revcomp(product[-9:]),revcomp(product[:9]))
    try:selected_single_insertion_filled_insert('ACGT',0);bad=False
    except ValueError:bad=True
    c('single-insertionclassrejectsinadequateexplicitcircle',bad,True)
    g=protocol_graph();c('ExoVdebrisnotimplicitlydiscardedbeforeentry',g.steps['entry_boundary'].inputs[0].material,'exov_product')
    c('gradientmethodkeepsotherphysicalfractions',graph_for_enrichment('gradient').terminals['other_fractions'].disposition,Disposition.UNRESOLVED)
    c('sourcealternativesremainbranchesnotnewentries',len({graph_for_enrichment(k).name for k in ('exo','gradient','both')}),3)
    c('sourcegeometrylimitsremainprivate',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
