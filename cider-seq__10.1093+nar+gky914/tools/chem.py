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

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False,extra_inputs=()):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION))+tuple(extra_inputs),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')


SOURCE='PublishedPMC6344846mainMethods67–95,Fig1; publishedsupplementTable1/2andFig8; DOI10.1093/nar/gky914. AnalyticalextractedDNArouteonly, no biologicalengineering.'

def select(g,key,x,quantity):
    selected=g.product(key+'_selected','Source selected '+quantity,K.MOLECULAR,state=x.state)
    rest=g.product(key+'_remainder','Untransferred stock after '+quantity+' allocation; amount unassigned',K.MOLECULAR,state=x.state)
    g.react(key,'Source actual selected contribution: '+quantity,S.SPLIT,((x,R.TEMPLATE),),(selected,rest),changes_molecule=False)
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred stock; no whole-stock consumption or invented recovery.'))
    return selected

def graph_for_fraction(size_selected=False):
    if type(size_selected)is not bool:raise ValueError('sourceSSorNSSfractionrequired')
    g=ProtocolGraph('CIDER-Seq — '+('SS' if size_selected else 'NSS'),authority=SOURCE,representation='published_instance')
    x=g.input('extracted','Analyticalsource extractedtotalDNA; circular/linearDNAmixture, exactsequencesunassigned',K.MOLECULAR,primary=True)
    x=select(g,'initial_allocate',x,'5ug totalnucleicacid for initialBluePippin')
    x=cleanup(g,'initial_size',x,'SourceBluePippin0.75%gelcassette0.8–5kbwindow')
    x=select(g,'rca_allocate',x,'5ul size-selectedtemplate in20ul RCA reaction')
    nt=g.input('rca_nt','Actual source1mM dNTP contribution to RCA',K.MOLECULAR)
    pp=g.input('rca_ppase','Actual source inorganicpyrophosphatase0.02U contribution',K.CATALYST)
    primer=g.input('random_primer','SourceexoresistantrandomprimerSO181; exactlength/PSinventoryunreviewed',K.OLIGO);phi=g.input('phi29','Actual sourcePhi29 10U contribution;30C18h/65C2min programme',K.CATALYST);y=g.product('rca_mix','SourceRCAamplified/branchedDNA plusoriginalcircles,linearDNA,primersandenzymeremnants; PCR-freeisnotunamplified',K.MOLECULAR)
    g.react('rca','Sourcerandomcircleamplification',S.EXTENSION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(phi,R.CATALYST),(pp,R.CATALYST),(nt,R.STRUCTURAL)),(y,))
    x=cleanup(g,'rca_clean',y,'Sourcesodiumacetate/ethanolprecipitation')
    x=select(g,'debranch_allocate',x,'10ug amplifiedDNA for debranching')
    x=op(g,'debranch',x,'SourcePhi29WITHOUTaddedprimer30°C2h;65°C2min',S.EXTENSION,reagent='sourcePhi29debranchingenzyme',catalyst=True,product='SourcepostPhi29DNAmixture; exactbranchtopologyandremnantsunassembled')
    x=cleanup(g,'debranch_clean',x,'Sourcesodiumacetate/ethanolprecipitation')
    x=op(g,'s1',x,'SourceS1nuclease37°C30min;EDTA/70°C10minstop',S.DIGESTION,reagent='sourceS1nuclease',catalyst=True,product='Sourcebranch-releasedDNA plusdigestednucleotideresidues/unreleasedDNA/reagents; digestionnotphysicalremoval')
    x=cleanup(g,'s1_clean',x,'Sourcesodiumacetate/ethanolprecipitation')
    ecoli=g.input('repair_ecoli','Actual source EcoliDNApolymeraseI10U contribution',K.CATALYST)
    repairnt=g.input('repair_nt','Actual source1mMdNTP contribution to50ulrepair',K.MOLECULAR)
    x=op(g,'repair',x,'SourceT4DNApolymerase/EcoliDNApolI25°C1h;75°C20min',S.END_REPAIR,reagent='Actual sourceT4DNApolymerase3U contribution',catalyst=True,extra_inputs=((ecoli,R.CATALYST),(repairnt,R.STRUCTURAL)))
    x=op(g,'dephos',x,'Sourcealkalinephosphatase37°C10min;75°C5min',S.DEPHOSPHORYLATION,reagent='sourcealkalinephosphatase',catalyst=True,product='Source dephosphorylatedrepairedDNA plusenzyme/reagent/repairremnants')
    x=cleanup(g,'repair_clean',x,'SourceKAPAPureBeads1.5xcleanup; separate analyticalQubit allocation follows')
    repairq=g.product('repair_qc','Source repairedDNA Qubit analytical fraction; volume unassigned',K.QC)
    repairrest=g.product('repair_production','RepairedDNA stock remaining after analyticalQubit fraction',K.MOLECULAR,state=x.state)
    g.react('repair_qc_allocate','Separate actual repairedDNA analyticalQC contribution',S.SPLIT,((x,R.TEMPLATE),),(repairq,repairrest),changes_molecule=False)
    g.mark_terminal(Terminal(repairq.id,Disposition.QC,'Actual sourceQubit analytical material, never returned to production.'));x=repairrest
    x=select(g,'library_allocate',x,'approximately200ng persample for barcodedSMRTlibrary')
    x=op(g,'library_end_repair',x,'SourcePacBioT4DNApolymerase/PNKendrepair',S.END_REPAIR,reagent='sourcehistoricalSMRTBellBarcodedAdapterCompletePrepKit96enzymes',catalyst=True)
    a=g.input('barcode_adapter','SourcePacBiobarcodedbluntadapter; exactsamplebarcode/hairpinbodyunassigned',K.ADAPTER);cat=g.input('ligase','Sourcehistoricalkitbluntligationenzymes',K.CATALYST);y=g.product('barcoded','Actualsourcebarcodedlibraryclass plusunligatedDNA/adapter/reagents; nohandmadehairpin',K.MOLECULAR)
    g.react('ligate','Sourcebluntbarcodedadapterligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    other=g.input('other_libraries','Other8sourcebarcodedDNApreparations; sourceincludesunrelatedprojects,noassignedsequences',K.MOLECULAR);pooled=g.product('pooled','Sourcepoolof9barcodedpreparationswithpossibleunligatedcompanionmaterial',K.MOLECULAR)
    own=select(g,'own_pool_allocate',y,'selected contribution of represented barcodedlibrary; amount unassigned')
    others=select(g,'other_pool_allocate',other,'aggregate selected contributions of other8sourcebarcodedpreparations; individual amounts/identities unassigned')
    g.react('pool','Sourcepool9barcodedpreparations',S.POOL,((own,R.TEMPLATE),(others,R.TEMPLATE)),(pooled,))
    x=op(g,'exonucleases',pooled,'Sourcehistoricalkitexonucleasetreatment',S.DIGESTION,reagent='Sourceexonucleaseinventory/conditionsunreviewed',catalyst=True,product='SourceSMRTbell-enrichedpreparationplusdigestionproducts/reagents; noinstantcleanuppresumed')
    fraction=g.product('selected_fraction','Source'+('>3kbSS' if size_selected else 'NSS')+'analyticallibraryaliquot; exactperlibrarypooling/routingtimingunassembled',K.MOLECULAR);otherfraction=g.product('companion_fraction','Separateactual'+('NSS' if size_selected else 'SS')+'comparisonfraction; noassumeddiscard',K.MOLECULAR)
    condition=g.input('routing','Main93barcoded/pool/exonuclease/aliquotrouting boundary; exactper-libraryversuspoolselectiontimingunreviewed',K.CONDITION)
    g.react('fraction_boundary','SourceSS/NSSphysicalaliquotallocation; exactroutingunassembled',S.MAGIC,((x,R.TEMPLATE),(condition,R.CONDITION)),(fraction,otherfraction),reason='KnownactualSS/NSSpairedfractions; sourceper-libraryselectionversuspooltimingnotforcedintoanunsupportedprecisesplit.');g.mark_terminal(Terminal(otherfraction.id,Disposition.UNRESOLVED,'Actualcomparisonlibraryreportedseparately; notduplicatedinselectedacquisition.'));x=fraction
    if size_selected:x=cleanup(g,'post_size',x,'SourceBluePippin>3kbselection')
    retained=g.product('production_remainder','SourceproductionlibraryretainedforP6binding; exactQCaliquotamountunassigned',K.MOLECULAR);diagnostic=g.product('qc_aliquot','SourceBioanalyzer12kbDNAChip/Qubitanalyticalaliquotmaterial; notpreparativesubstrate',K.QC)
    g.react('qc_allocate','AllocateanalyticalQCmaterialseparatelyfromproductionlibrary',S.SPLIT,((x,R.TEMPLATE),),(retained,diagnostic));g.mark_terminal(Terminal(diagnostic.id,Disposition.QC,'Sourceanalyticalqualityinspection/quantificationonly; noBioanalyzer/Qubitreactionmaterialfedintopolymerasebinding.'));x=retained
    poly=g.input('p6_polymerase','SourceP6DNApolymerasebindingkit2.0complexpartner; exactprimergeometryunreviewed',K.MOLECULAR);ready=g.product('ready','Sourcepolymerase-boundSMRTbellpreparation and bindingmixture companions unassembled',K.MOLECULAR)
    g.react('bind','SourceP6DNApolymerasebindingkit2.0complexpreparation',S.BINDING,((x,R.TEMPLATE),(poly,R.STRUCTURAL)),(ready,));x=ready
    x=select(g,'load_allocate',x,'actual loading contribution from P6-boundlibrary; amount unassigned')
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
    c('Repaired stock reaches analytical allocation after phosphatase and cleanup',g.steps['repair_qc_allocate'].inputs[0].material,'repair_clean_kept')
    c('NSSdoesnotreceivepostenrichmentselection','post_size' in g.steps,False)
    c('SSreceivessourcepostselection','post_size' in graph_for_fraction(True).steps,True)
    c('analyticalQCdoesnotbecomepolymerasebindingsubstrate',(g.steps['bind'].inputs[0].material,g.terminals['qc_aliquot'].disposition),('production_remainder',Disposition.QC))
    c('pairedcomparisonfractionisretained',g.terminals['companion_fraction'].disposition,Disposition.UNRESOLVED)
    c('RCA consumes actual selected source template volume',g.steps['rca'].inputs[0].material,'rca_allocate_selected')
    c('Debranching consumes selected amplifiedDNA mass',g.steps['debranch'].inputs[0].material,'debranch_allocate_selected')
    c('Pool consumes selected contributions rather than all barcoded stocks',tuple(i.material for i in g.steps['pool'].inputs),('own_pool_allocate_selected','other_pool_allocate_selected'))
    c('Original kit scope remains uncovered','pending' in g.source_coverage,False)
    c('geometrylimitsremainunpublishable',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
