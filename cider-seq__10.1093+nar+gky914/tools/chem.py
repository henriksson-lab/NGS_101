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


SOURCE='PublishedPMC6344846mainMethods67–95,Fig1; publishedsupplementTable1/2andFig8; DOI10.1093/nar/gky914. AnalyticalextractedDNArouteonly, no biologicalengineering.'

def graph_for_fraction(size_selected=False):
    if type(size_selected)is not bool:raise ValueError('sourceSSorNSSfractionrequired')
    g=ProtocolGraph('CIDER-Seq — '+('SS' if size_selected else 'NSS'),authority=SOURCE,representation='published_instance')
    x=g.input('extracted','Analyticalsource extractedtotalDNA; circular/linearDNAmixture, exactsequencesunassigned',K.MOLECULAR,primary=True)
    x=cleanup(g,'initial_size',x,'SourceBluePippin0.75%gelcassette0.8–5kbwindow')
    primer=g.input('random_primer','SourceexoresistantrandomprimerSO181; exactlength/PSinventoryunreviewed',K.OLIGO);phi=g.input('phi29','SourcePhi29/dNTP/pyrophosphatase30°C18h;65°C2minstop',K.CATALYST);y=g.product('rca_mix','SourceRCAamplified/branchedDNA plusoriginalcircles,linearDNA,primersandenzymeremnants; PCR-freeisnotunamplified',K.MOLECULAR)
    g.react('rca','Sourcerandomcircleamplification',S.EXTENSION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(phi,R.CATALYST)),(y,))
    x=cleanup(g,'rca_clean',y,'Sourcesodiumacetate/ethanolprecipitation')
    x=op(g,'debranch',x,'SourcePhi29WITHOUTaddedprimer30°C2h;65°C2min',S.EXTENSION,reagent='sourcePhi29debranchingenzyme',catalyst=True,product='SourcepostPhi29DNAmixture; exactbranchtopologyandremnantsunassembled')
    x=cleanup(g,'debranch_clean',x,'Sourcesodiumacetate/ethanolprecipitation')
    x=op(g,'s1',x,'SourceS1nuclease37°C30min;EDTA/70°C10minstop',S.DIGESTION,reagent='sourceS1nuclease',catalyst=True,product='Sourcebranch-releasedDNA plusdigestednucleotideresidues/unreleasedDNA/reagents; digestionnotphysicalremoval')
    x=cleanup(g,'s1_clean',x,'Sourcesodiumacetate/ethanolprecipitation')
    x=op(g,'repair',x,'SourceT4DNApolymerase/EcoliDNApolI25°C1h;75°C20min',S.END_REPAIR,reagent='sourceT4DNApolymerase/EcoliDNApolymeraseI/dNTPs',catalyst=True)
    x=op(g,'dephos',x,'Sourcealkalinephosphatase37°C10min;75°C5min',S.DEPHOSPHORYLATION,reagent='sourcealkalinephosphatase',catalyst=True,product='Source dephosphorylatedrepairedDNA plusenzyme/reagent/repairremnants')
    x=cleanup(g,'repair_clean',x,'SourceKAPAPureBeads1.5xcleanup/Qubit')
    x=op(g,'library_end_repair',x,'SourcePacBioT4DNApolymerase/PNKendrepair',S.END_REPAIR,reagent='sourcehistoricalSMRTBellBarcodedAdapterCompletePrepKit96enzymes',catalyst=True)
    a=g.input('barcode_adapter','SourcePacBiobarcodedbluntadapter; exactsamplebarcode/hairpinbodyunassigned',K.ADAPTER);cat=g.input('ligase','Sourcehistoricalkitbluntligationenzymes',K.CATALYST);y=g.product('barcoded','Actualsourcebarcodedlibraryclass plusunligatedDNA/adapter/reagents; nohandmadehairpin',K.MOLECULAR)
    g.react('ligate','Sourcebluntbarcodedadapterligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    other=g.input('other_libraries','Other8sourcebarcodedDNApreparations; sourceincludesunrelatedprojects,noassignedsequences',K.MOLECULAR);pooled=g.product('pooled','Sourcepoolof9barcodedpreparationswithpossibleunligatedcompanionmaterial',K.MOLECULAR)
    g.react('pool','Sourcepool9barcodedpreparations',S.POOL,((y,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,))
    x=op(g,'exonucleases',pooled,'Sourcehistoricalkitexonucleasetreatment',S.DIGESTION,reagent='Sourceexonucleaseinventory/conditionsunreviewed',catalyst=True,product='SourceSMRTbell-enrichedpreparationplusdigestionproducts/reagents; noinstantcleanuppresumed')
    fraction=g.product('selected_fraction','Source'+('>3kbSS' if size_selected else 'NSS')+'analyticallibraryaliquot; exactperlibrarypooling/routingtimingunassembled',K.MOLECULAR);otherfraction=g.product('companion_fraction','Separateactual'+('NSS' if size_selected else 'SS')+'comparisonfraction; noassumeddiscard',K.MOLECULAR)
    condition=g.input('routing','Main93barcoded/pool/exonuclease/aliquotrouting boundary; exactper-libraryversuspoolselectiontimingunreviewed',K.CONDITION)
    g.react('fraction_boundary','SourceSS/NSSphysicalaliquotallocation; exactroutingunassembled',S.MAGIC,((x,R.TEMPLATE),(condition,R.CONDITION)),(fraction,otherfraction),reason='KnownactualSS/NSSpairedfractions; sourceper-libraryselectionversuspooltimingnotforcedintoanunsupportedprecisesplit.');g.mark_terminal(Terminal(otherfraction.id,Disposition.UNRESOLVED,'Actualcomparisonlibraryreportedseparately; notduplicatedinselectedacquisition.'));x=fraction
    if size_selected:x=cleanup(g,'post_size',x,'SourceBluePippin>3kbselection')
    retained=g.product('production_remainder','SourceproductionlibraryretainedforP6binding; exactQCaliquotamountunassigned',K.MOLECULAR);diagnostic=g.product('qc_aliquot','SourceBioanalyzer12kbDNAChip/Qubitanalyticalaliquotmaterial; notpreparativesubstrate',K.QC)
    g.react('qc_allocate','AllocateanalyticalQCmaterialseparatelyfromproductionlibrary',S.SPLIT,((x,R.TEMPLATE),),(retained,diagnostic));g.mark_terminal(Terminal(diagnostic.id,Disposition.QC,'Sourceanalyticalqualityinspection/quantificationonly; noBioanalyzer/Qubitreactionmaterialfedintopolymerasebinding.'));x=retained
    poly=g.input('p6_polymerase','SourceP6DNApolymerasebindingkit2.0complexpartner; exactprimergeometryunreviewed',K.MOLECULAR);ready=g.product('ready','Sourcepolymerase-boundSMRTbellpreparation and bindingmixture companions unassembled',K.MOLECULAR)
    g.react('bind','SourceP6DNApolymerasebindingkit2.0complexpreparation',S.BINDING,((x,R.TEMPLATE),(poly,R.STRUCTURAL)),(ready,));x=ready
    g.sequence(x,SequencingPlan(limitation='ActualPacBioRS2SMRTCellv3.0/MagBeadload,360minmovie. DeConcat is computational sequence reduction; itdoesnotphysicallyshortenRCAproductsornativecircles. Historicaladapter/barcode/primergeometryunassembled.'))
    g.add_source_item(SourceItem('pending','HistoricalPacBioSMRTBellBarcodedAdapterCompletePrepKit96/P6originals','Main93exacthistoricalkitidentities; sourcePSrandomprimer,barcodes/hairpinsandexonuclease/recoverygeometryunreviewed. Published2018originalnotmergedwith2020detailedprotocolorits2022correction.'))
    return g

def protocol_graph():return graph_for_fraction()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('nativeinputisnotaliasedasPCRorRCAcopies',g.steps['rca'].outputs,('rca_mix',))
    c('debranchingretainsnoaddednewprimer',len(g.steps['debranch'].inputs),2)
    c('S1digestionproductsremainuntilphysicalprecipitation',g.steps['s1_clean'].inputs[0].material,'s1_product')
    c('sourcePNKendrepairfollowspriorphosphatase',g.steps['library_end_repair'].inputs[0].material,'repair_clean_kept')
    c('NSSdoesnotreceivepostenrichmentselection','post_size' in g.steps,False)
    c('SSreceivessourcepostselection','post_size' in graph_for_fraction(True).steps,True)
    c('analyticalQCdoesnotbecomepolymerasebindingsubstrate',(g.steps['bind'].inputs[0].material,g.terminals['qc_aliquot'].disposition),('production_remainder',Disposition.QC))
    c('pairedcomparisonfractionisretained',g.terminals['companion_fraction'].disposition,Disposition.UNRESOLVED)
    c('geometrylimitsremainunpublishable',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
