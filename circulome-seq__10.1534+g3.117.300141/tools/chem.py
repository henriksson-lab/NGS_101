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
        added=g.input('exov_supplies','ActualsourceExoV-regimen additions; per-roundstockallocations/countunassembled, no continuousactivityafterheatasserted',K.MOLECULAR)
        enzyme=g.product('exov_enzyme','ActualsourceExoV catalystcontributions overregimen',K.CATALYST);atp=g.product('exov_atp','ActualsourceATP cofactor/substrate contribution, amountunassigned',K.MOLECULAR);other=g.product('exov_other','ActualsourceDTT/100ugmlampicillin molecularcomponents',K.MOLECULAR);cond=g.product('exov_conditions','Source400Uug3days/bufferpH7.9;70C30minafterrounds; humanprotocolseparate',K.CONDITION)
        g.react('exov_parts','Components ofactualsourceadditiongroups, notextrafreshdoses',S.GENERIC,((added,R.TEMPLATE),),(enzyme,atp,other,cond),changes_molecule=False)
        digested=g.product('exov_product','Nativecircle-enrichedDNA plusdigestedDNAresidues/possiblysurvivinglinearDNA/sourceinactiveenzymes; no physicalpurification',K.MOLECULAR)
        g.react('exov','SourceExoVregimen andbetween-roundheat, no perfectcircleenrichment',S.MAGIC,((x,R.TEMPLATE),(enzyme,R.CATALYST),(atp,R.STRUCTURAL),(other,R.STRUCTURAL),(cond,R.CONDITION)),(digested,),reason='Actualmulti-roundenzymaticoperation; roundcount/aliquots/activevsinactiveenzymeinventories unassembled; heatnotcleanup.');x=digested
    if enrichment in {'gradient','both'}:
        condition=g.input('gradient_condition','SourceCsCl1.55g/ml/ethidiumgradient500237g2.5h; separateplasmidreference tube',K.SELECTION);bottom=g.product('bottom_fraction','Sourcecollectedapproximateclosedcirclefraction; noabsolutecomposition/purityclaim',K.MOLECULAR);top=g.product('other_fractions','Actualothergradientfractionsincludinglinearandnicked/relaxedcircles; notallcirclescollected',K.MOLECULAR)
        g.react('gradient','Sourcebiophysicaltopologyfractionation',S.CLEANUP,((x,R.TEMPLATE),(condition,R.SELECTION)),(bottom,top));g.mark_terminal(Terminal(top.id,Disposition.UNRESOLVED,'Otherfractionsincludesrelaxed/nickedcircles; notforcedintopurelinearDNAwaste.'));x=bottom
        x=cleanup(g,'butanol',x,'SourceCsCl/TEsaturated1butanolethidiumextraction')
        x=cleanup(g,'dialyze',x,'SourceTEdialysis2days4°C')
    x=op(g,'entry_boundary',x,'SourceeccDNApreparedfortagmentation; exactExoVbranchtransfer/cleanupunreviewed',reagent='Do notinventphysicalremovalofExoVresiduesorRCAstep',product='SourceenrichedDNApreparationwithunassembledcarryoverfates')
    tn5=g.input('nextera','Actualselected0.5ulsourceNexteraXTloadedtagmentase stock; no separatefreeTn5addition',K.MOLECULAR)
    enzyme=g.product('tagmentase','CopresentTn5 catalystcomponent ofsame loadedstock',K.CATALYST);adapters=g.product('tagment_adapters','CopresentsourceXTadaptorclasses ofsame loadedstock; wholebodies/lengthsunreviewed',K.ADAPTER);other=g.product('tagment_other','Otheractualloaded-stockcomponents; occupancy/freecompanionstatesunassembled',K.MOLECULAR)
    g.react('tagment_parts','Sameactualloaded-stockfunctionalcomponents, notadditionalfresh enzyme/adapterdoses',S.GENERIC,((tn5,R.TEMPLATE),),(enzyme,adapters,other),changes_molecule=False)
    conditions=g.input('tagment_conditions','Source37C3minwithoutshaking attenuated' if attenuated else 'Source37C30min low-input; genomic1ng/1.5ulcontrolnotinherited',K.CONDITION)
    y=g.product('tagged','Tn5-fragmentedDNA plusunreactednativecircles/linearDNA/adaptors/sourcecompanions;9bpduplicationnotcompleteuntilcopying',K.MOLECULAR)
    g.react('tagment','Actualsource0.5ulXT loadedstock tagmentation',S.TAGMENTATION,((x,R.TEMPLATE),(adapters,R.STRUCTURAL),(enzyme,R.CATALYST),(other,R.STRUCTURAL),(conditions,R.CONDITION)),(y,))
    supplies=g.input('pcr_supplies','AggregateactualsourceXT-PCRpreparation supplies; originalkit/polymer/primer/dosesunreviewed, stockcount/packagingunassigned',K.MOLECULAR)
    pol=g.product('pcr_polymerase','ActualsourcePCRpolymerasefunction, identityunreviewed',K.CATALYST);nt=g.product('pcr_nt','ActualsourcePCRnucleotides; dose/identityconventionunassembled',K.MOLECULAR);f=g.product('pcr_f','ActualsourceforwardPCRprimerclass, wholebodyunreviewed',K.OLIGO);r=g.product('pcr_r','ActualsourcereversePCRprimerclass, wholebodyunreviewed',K.OLIGO);buf=g.product('pcr_conditions','ActualsourcePCRconditions, originalrunprogramunreviewed',K.CONDITION)
    g.react('pcr_parts','ActualonePCRpreparationfunctionalcomponents, notnewextra stocks',S.GENERIC,((supplies,R.TEMPLATE),),(pol,nt,f,r,buf),changes_molecule=False)
    filled=g.product('gapfill_product','DNAcontext aftersourcefirstPCRcyclegapfill; same-tube catalyst/unusedprimers/NTfunctions representedseparately; wholebackboneunassembled',K.MOLECULAR)
    retained_pol=g.product('retained_pcr_polymerase','SamePCRenzyme carriedfromgapfillcycle, notfreshaddition',K.CATALYST);retained_nt=g.product('retained_pcr_nt','SamePCRunincorporatednucleotidefunctions, amountunassembled',K.MOLECULAR);retained_f=g.product('retained_pcr_f','Sameactualforwardprimer component',K.OLIGO);retained_r=g.product('retained_pcr_r','Sameactualreverseprimercomponent',K.OLIGO)
    g.react('gapfill','SourcefirstPCRcycle fillsTn5gaps; conditional9bpduplicationclass',S.EXTENSION,((y,R.TEMPLATE),(pol,R.CATALYST),(nt,R.STRUCTURAL),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(buf,R.CONDITION)),(filled,retained_pol,retained_nt,retained_f,retained_r))
    amplified=g.product('library_pcr_product','SourceamplifiedDNA plusoriginaltemplates/reagents; nativecircle/gappedintermediate/linearcopiedlibrary aredistinct',K.MOLECULAR)
    g.react('library_pcr','SamePCRprogram continuation; sourcegelchoosescycles/10–12reportedsufficient, not10–12extracyclesaftergapfill',S.PCR,((filled,R.TEMPLATE),(retained_f,R.STRUCTURAL),(retained_r,R.STRUCTURAL),(retained_pol,R.CATALYST),(retained_nt,R.STRUCTURAL)),(amplified,),geometry=PCRGeometry.unavailable(scene('HistoricalXTfulladapter/primer/endgeometryunassembled'),reason='SourceonePCRprogram containsinitialgapfill; thermal/cycleconvention/wholeproductunassembled.'));x=amplified
    x=qc(g,'analytical_gel',x,'Sourcelibraryassessmentaliquotandretainedsequencingpreparation')
    loaded=g.product('loading_selected','ActualselectedsourceDNA-librarysequencingcontribution; amountunassigned',K.MOLECULAR);rest=g.product('loading_remainder','UntransferredsourceDNAstock; fateunassigned',K.MOLECULAR)
    g.react('loading','ActualsourceDNA-libraryallocation afteranalyticalQA',S.SPLIT,((x,R.TEMPLATE),),(loaded,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actualuntransferredstock.'))
    g.sequence(loaded,SequencingPlan(limitation='Actualpaired-endreadmapping/9bpduplicationanalysisestablishedbyFig1andMethods; preciseplatform/run/adapter/readprimergeometryunreviewed. Computational9bp-circleevidenceisnotaperfectselectionorphysicaldiscardrule; sourceuntreatedgenomiccontrolsseparate.'))
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
    c('Tagmentationusesproducedloaded-stockenzyme',g.steps['tagment'].inputs[2].material,'tagmentase')
    c('LoadedTn5notaseparatefreshinput','tagmentase' in g.steps['tagment_parts'].outputs,True)
    c('ExoVhasactualATPcontributor',g.steps['exov'].inputs[2].material,'exov_atp')
    c('PCRcontinueswithretainedsameenzyme',g.steps['library_pcr'].inputs[3].material,'retained_pcr_polymerase')
    c('PCRstocknotaddedagainaftergapfill',g.steps['gapfill'].inputs[1].material,'pcr_polymerase')
    c('FinalDNAstockremainder retained',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('sourcegeometrylimitsremainprivate',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
