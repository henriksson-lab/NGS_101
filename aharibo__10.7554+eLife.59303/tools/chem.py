"""Protocol-owned selected chemical classes; incomplete source-qualified graph."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name,text):
    sc=Scene();sc.strand(name,[Segment(name,text,placeholder=True)]);return sc

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered material; no perfect purity or quantitative partition assumed',K.MOLECULAR)
    waste=g.product(key+'_waste','Actually separated cleanup fractions with unquantified molecular carryover',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,waste));g.discard(waste,'Source physical separation, not perfect molecular classification.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source mixture',catalyst=False,product=None,state=None,chemical_state=None):
    molecular=kind==S.GENERIC or key=='rnase_stop'
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else (K.MOLECULAR if molecular else K.CONDITION))
    y=g.product(key+'_product',product or label+'; all physical reaction companions retained until separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    kwargs={}
    if kind==S.MAGIC:kwargs['reason']='Source operation established; exact internal chemistry/geometry not assembled from unreviewed historical originals.'
    if kind==S.PCR:kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exact source template/primer endpoint geometry unassembled, no ordinary complete library asserted.')
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-known ligation; full substrate/end geometry unassembled]'),),'Unreviewed geometry; not a claim of unavailable disclosure.')
    inputs=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else (R.STRUCTURAL if molecular else R.CONDITION))]
    if key=='rnase_stop':
        env=g.input(key+'_environment','Source 10min on ice',K.CONDITION);inputs.append((env,R.CONDITION))
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
    if kind==S.LIGATION:
        adapter=g.input(key+'_adapter','Source ligating molecular adapter; exact historical oligo/end inventory unassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter geometry unassembled]'));inputs.append((adapter,R.STRUCTURAL))
    if kind==S.PCR:
        for arm in ('f','r'):
            primer=g.input(key+'_primer_'+arm,'Source PCR primer family; exact oligo/endpoint unassembled',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL))
    if kind==S.GENERIC:kwargs['changes_molecule']=False
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained production DNA library',K.LIBRARY)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
from dataclasses import replace
def supplies(g,key,label,items):
    stock=g.input(key+'_supplies',label+'; physical source inventory, stock packaging/count unassembled',K.MOLECULAR)
    components=[g.product(key+'_'+name,text,kind) for name,text,kind in items]
    g.react(key+'_inventory','Bookkeeping of actual source supplies; no extra fresh doses or premix factory asserted',S.GENERIC,((stock,R.STRUCTURAL),),components,changes_molecule=False)
    return components

SOURCE='Published eLife59303v2 mainMethods335–365/405–407/453–465; primarypublisherXML/PDFidentity2021VOR checked. DOI10.7554/eLife.59303. Appendedreviewhistorynotoperativerecipe.'
@dataclass(frozen=True)
class NascentComplex:
    rna_id:str
    peptide_id:str
    azide_sites:tuple
    clicked_sites:tuple=()
    def __post_init__(self):
        if not self.rna_id or not self.peptide_id or len(set(self.azide_sites))!=len(self.azide_sites) or len(set(self.clicked_sites))!=len(self.clicked_sites) or set(self.azide_sites)&set(self.clicked_sites):raise ValueError('ExplicitindependentRNA/peptideidentities anddistinctreactivesitesrequired')
        if not all(isinstance(x,str) and x for x in self.azide_sites+self.clicked_sites):raise ValueError('Named source-selectedAHApeptidesitesrequired')
    def click_selected(self,sites):
        if not sites or len(set(sites))!=len(sites) or not set(sites)<=set(self.azide_sites):raise ValueError('CallerselectedunconsumedpeptideAHAazidesrequired')
        return replace(self,azide_sites=tuple(x for x in self.azide_sites if x not in sites),clicked_sites=self.clicked_sites+tuple(sites))
    @property
    def rna_biotin_sites(self):return ()
    def proteolyze_selected(self):return {'rna_id':self.rna_id,'peptide_remnants_with_biotin':self.clicked_sites,'other_protein_products':'unassembled','rna_direct_biotin':()}

def graph_for_scope(footprints=False,aha_present=True):
    if type(footprints)is not bool or type(aha_present)is not bool:raise ValueError('Explicitsourcepreparation/controlchoice required')
    g=ProtocolGraph('AHARIBO — '+('RPF' if footprints else 'RNAseq'),authority=SOURCE,representation='published_instance')
    x=g.input('cells','SourceAHARIBOmESC/early-neuronRNAexperiment orselectedRPFcellpreparation; actualsampleunassigned',K.MOLECULAR,primary=True)
    x=op(g,'methionine_deplete',x,'Source40minmethionine-depletionmedium',S.GENERIC,reagent='Sourcemethionine-freeculturemedium/leucine aslocallyspecified')
    if aha_present:
        aha=g.input('aha','Source250uML-azidohomoalanine, incorporatedintonascentprotein ratherthanRNA',K.MOLECULAR);medium=g.input('pulse_conditions','Source30minAHA pulse/cultureconditions',K.CONDITION);p=g.product('labeled_cells','NascentproteinAHA-containingcellpopulation plusunlabeledproteins/RNA andcompanions; incorporation/yield/siteoccupancynotpredicted',K.MOLECULAR)
        g.react('aha_pulse','SourceAHAmetabolicpulse',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(aha,R.STRUCTURAL),(medium,R.CONDITION)),(p,));x=p
    else:x=op(g,'no_aha',x,'SourceAHAminuscontrol; nocovalentAHAincorporationasserted',S.GENERIC,reagent='Sourcecontrolmedium')
    if footprints:
        x=op(g,'rpf_preparation',x,'SourceRPFMethods453CYHpretreatment, rapidiceharvest andRL001-1lysis; relationtocommonsBlocktimingnotunified',reagent='SourceRPFlysisRL001-1/CYH; exactschedulingrelativecommonmetabolicpulseunassembled',product='SourceRPFcelllysate beforeclarifiedinput; exactkitlysischemistryunreviewed, no commonRNAseqwash/DNaseI/proteaseinhibitorrecipeimposed')
        x=clean(g,'rpf_clarify',x,'SourceRPFclarifiedlysaterecovery; exactclarificationparametersunreviewed')
    else:
        x=op(g,'sblock',x,'Sourceanisomycin-containingproprietarysBlock10min thenice',S.GENERIC,reagent='SourcesBlockRM8; nocovalentRNAmodificationasserted')
        x=clean(g,'cell_wash',x,'SourcecoldPBS/sBlockwash; recovercellpreparation')
        salts,detergents,sblock,dnase,ribolock,inhibitors= supplies(g,'lysis','Actual source hypotonic lysis supplies',[
            ('buffer','Source 0.01M NaCl/0.01M MgCl2/0.01M Tris-HCl',K.CONDITION),
            ('detergents','Source 1% Triton-X100 and 1% sodium deoxycholate',K.MOLECULAR),
            ('sblock','Source 1x proprietary anisomycin-containing sBlock',K.MOLECULAR),
            ('dnase','Source 5U/ml DNase I',K.CATALYST),
            ('ribolock','Source 200U/ml RiboLock RNase inhibitor',K.MOLECULAR),
            ('inhibitors','Source 1x protease inhibitor cocktail',K.MOLECULAR)])
        lysed=g.product('lyse_product','Hypotonic lysate with source enzyme/inhibitor/detergent companions',K.MOLECULAR)
        g.react('lyse','Source hypotonic cell lysis with actual molecular contributors',S.EXTRACTION,((x,R.TEMPLATE),(salts,R.CONDITION),(detergents,R.STRUCTURAL),(sblock,R.STRUCTURAL),(dnase,R.CATALYST),(ribolock,R.STRUCTURAL),(inhibitors,R.STRUCTURAL)),(lysed,))
        x=lysed
        x=clean(g,'clarify',x,'Sourceclarifiedlysaterecovery; RNAseqmain18000g4C5min removesnuclei/debris')
    if footprints:
        x=op(g,'rnase',x,'SourceRNaseI1.3U45minRT withsBlock',S.DIGESTION,reagent='SourceRNaseIAM2295/WbufferRL001-4',product='RPF-containinglysateplusdigestedRNAproducts/ribosomes/peptides; nopurefootprintclassasserted')
        x=op(g,'rnase_stop',x,'SourceSuperaseIn10minonice',S.QUENCHING,reagent='SourceRNaseinhibitorAM2696')
    ligand=g.input('dbco_biotin','Source50uMdibenzocyclooctyne-PEG4-biotin760749; solubleRNA/RPFcaptureligand,notMag-DBCObeadsubstitution',K.MOLECULAR);conditions=g.input('click_conditions','SourceWbuffer100ul/SuperaseIn,1h4C9rpm',K.CONDITION);clicked=g.product('clicked_mix','Sourcepeptideazide/DBCOcoupledclasses carrybiotin; RNAassociationremainsnoncovalent. Otherligandreactivity/background/freeDBCO/reactionmaterialretained',K.MOLECULAR)
    inhibitor=g.input('click_rnase_inhibitor','Source 40 U Superase-In in RNAseq capture reaction; RPF inhibitor amount unassigned' if not footprints else 'Source Superase-In inhibitor contribution carried from RNase stop; no new 40 U dose asserted',K.MOLECULAR) if not footprints else None
    click_inputs=[(x,R.TEMPLATE),(ligand,R.STRUCTURAL),(conditions,R.CONDITION)]
    if inhibitor is not None:click_inputs.append((inhibitor,R.STRUCTURAL))
    g.react('click','Source solubleDBCO-biotin coupling beforemagneticcapture',S.CHEMICAL_CONVERSION,click_inputs,(clicked,))
    beads=g.input('capture_beads','Source50uleMagSi-cN018-eMS-001',K.SELECTION);bound=g.product('bound','Bead-associatedbiotinpeptide/ribosome/RNAcomplexpopulation, backgroundnotzero',K.MOLECULAR);flow=g.product('capture_flow','Actualmagnetseparatedsourceflowthrough',K.WASTE)
    g.react('capture','Source30min4Ccapture thenmagnetonice',S.AFFINITY_CAPTURE,((clicked,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.discard(flow,'Sourcesupernatantsdiscarded; noidealchemicalpartitionassumed.')
    x=clean(g,'wash',bound,'Sourcetwo500ulPBS/Triton/sBlock/RiboLockwashes, beadsresuspendedWbufferandtransferred')
    sds,pk= supplies(g,'proteolysis','Actual source proteolysis additions',[
        ('sds','Source 20 ul 10% SDS added to 200 ul W-buffer suspension',K.MOLECULAR),
        ('pk','Source 5 ul Proteinase K Qiagen19131',K.CATALYST)])
    proteolysis_env=g.input('proteolysis_conditions','Source 37C 75min',K.CONDITION)
    digested=g.product('proteolyze_product','RNA and proteolysis/peptide-biotin remnants, beads and reagents retained until phase separation; no biotin transfer to RNA',K.MOLECULAR)
    g.react('proteolyze','Source SDS/Proteinase K treatment before phenol extraction',S.DEPROTEINIZATION,((x,R.TEMPLATE),(sds,R.STRUCTURAL),(pk,R.CATALYST),(proteolysis_env,R.CONDITION)),(digested,))
    x=digested
    phenol=g.input('phenol_reagents','Source equal volume acid phenol/chloroform/isoamyl alcohol pH4.5',K.MOLECULAR)
    mixed=g.product('phenol_mixture','Actual phase-extraction mixture before centrifugation/aqueous transfer',K.MOLECULAR)
    g.react('phenol_mix','Add actual organic extraction reagents before phase separation',S.GENERIC,((x,R.TEMPLATE),(phenol,R.STRUCTURAL)),(mixed,),changes_molecule=False)
    x=mixed
    x=clean(g,'phenol',x,'Sourceacidphenol:chloroform:isoamylalcoholpH4.5,14000g5min aqueousphasetransfer')
    precip=g.input('precipitation_reagents','Source 500 ul isopropanol and 2 ul GlycoBlue',K.MOLECULAR)
    precipitated=g.product('precipitate_product','Precipitated RNA plus solvent/carrier companions before centrifugation',K.MOLECULAR)
    env=g.input('precipitation_conditions','Source 3min room temperature then -80C overnight',K.CONDITION)
    g.react('precipitate','Source isopropanol/GlycoBlue precipitation',S.GENERIC,((x,R.TEMPLATE),(precip,R.STRUCTURAL),(env,R.CONDITION)),(precipitated,),changes_molecule=False)
    precipitated=clean(g,'precipitation_spin',precipitated,'Source 14000g 30min; remove supernatant before ethanol wash')
    ethanol=g.input('ethanol_wash','Source 500 ul 70% ethanol pellet wash',K.MOLECULAR)
    wash_mix=g.product('ethanol_wash_mixture','Actual precipitation/wash material before final recovered pellet',K.MOLECULAR)
    g.react('ethanol_add','Source ethanol wash inventory; physical centrifugation boundary follows',S.GENERIC,((precipitated,R.TEMPLATE),(ethanol,R.STRUCTURAL)),(wash_mix,),changes_molecule=False)
    x=wash_mix
    x=clean(g,'pellet',x,'Source ethanol wash 14000g 10min; air-dry pellet and resuspend 10ul water')
    if footprints:
        x=clean(g,'rpf_size',x,'SourceTBE-ureagel25–35ntRNAexcision/recovery')
        x=op(g,'library',x,'SourceRiboLacekitmodule2RL001_mod2 librarypreparation',reagent='Historicalmodule2internaloligos/ends/RT/PCRunreviewed; not2018wholeRiboLacerecipeimport')
        acquisition='SourceHiSeq2500single-end50bp; exactprimer/index/kitendgeometryunassembled.'
    else:
        selected=g.product('library_rna_5ng','Actual 5 ng RNA selected from recovered 10 ul RNA stock',K.MOLECULAR)
        stock=g.product('rna_stock_remainder','Untransferred recovered RNA stock if any; amount and later fate unassembled',K.MOLECULAR)
        g.react('library_allocate','Source 10 ul RNA recovery then 5 ng library-input allocation; no invented concentration',S.SPLIT,((x,R.TEMPLATE),),(selected,stock))
        g.mark_terminal(Terminal(stock.id,Disposition.UNRESOLVED,'Untransferred stock amount/fate unresolved; optional RNA assays not asserted universal.'))
        x=selected
        x=op(g,'library',x,'Source5ngRNA SMARTSeqStranded634443 librarypreparation',reagent='HistoricalTakara634443 internaloligos/RT/PCR/readgeometryunreviewed; notSMARTseqfamilydefault')
        x=qc_split(g,'library_qc',x,'SourceBioAnalyzerDNA1000/QubitdsDNAHS analyticalaliquots,retainedproductionlibrary')
        acquisition='SourceHiSeq2500; exactrunlength/primer/index/kitgeometryunreviewed.'
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=acquisition+' PeptidemetaboliccaptureidentifiesassociatedRNA/RPFs,notRNAcovalentbiotinmark; separateproteomicsonMag-DBCObeadsnotsequencinglibrary.'))
    g.add_source_item(SourceItem('pending','ActualhistoricalSMARTSeqStranded634443/RiboLacemodule2/DBCO/sBlockvendororiginals andfullassociatedsourcefigures/data','Sourceknownbranches throughacquisition, exactkit/modifiedmolecule/primer/readgeometry unreviewed; noacceptedabsence.'))
    return g

def protocol_graph():return graph_for_scope()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=NascentComplex('rna','peptide',('aha1','aha2'));b=a.click_selected(('aha1',))
    c('Clickconsumesonlyselectedpeptideazide',b.azide_sites,('aha2',))
    c('RNAidentityunchangedandnotbiotinylated',(b.rna_id,b.rna_biotin_sites),('rna',()))
    c('Proteolysisretainsbiotinonpeptideremnants',b.proteolyze_selected()['peptide_remnants_with_biotin'],('aha1',))
    try:b.click_selected(('aha1',));bad=False
    except ValueError:bad=True
    c('Cannotreactsameconsumedazidetwice',bad,True)
    g=protocol_graph();r=graph_for_scope(True)
    c('RNAseqdoesnotborrowRPFdigestion','rnase' not in g.steps and 'rnase' in r.steps,True)
    c('RPFpreparationbypassescommonRNAseqwashandduplicatehypotoniclysis','cell_wash' not in r.steps and 'lyse' not in r.steps and r.steps['rnase'].inputs[0].material=='rpf_clarify_kept',True)
    c('RNAseqcaptureswithsolubleligandbeforebeads',g.steps['capture'].inputs[0].material,'clicked_mix')
    c('Actual phenol reagent addition precedes phase separation',g.steps['phenol_mix'].inputs[0].material=='proteolyze_product' and g.steps['phenol'].inputs[0].material=='phenol_mixture',True)
    c('AHAminuscontroldoesnotinventAHAinput','aha' not in graph_for_scope(False,False).materials,True)
    p=Source(sources_root()/'aharibo__10.7554+eLife.59303'/'eLife.59303_v2.xml.txt','get_sources --slugs aharibo__10.7554+eLife.59303')
    if have(c,p,label='publishedsource'):
        t=p.path.read_text();c('SourceexplicitdistinctRNAandproteomicscapture','Mag-DBCO beads' in t and 'dibenzocyclooctyne-PEG4-biotin' in t,True)
    c('Source Proteinase K is a catalyst consumed before phase extraction',g.materials['proteolysis_pk'].kind==K.CATALYST and 'proteolysis_pk' in tuple(z.material for z in g.steps['proteolyze'].inputs),True)
    c('Actual5ngallocationprecedeslibraryfactory',g.steps['library'].inputs[0].material,'library_rna_5ng')
    c('RNA stock remainder is unresolved not library input',g.terminals['rna_stock_remainder'].disposition,Disposition.UNRESOLVED)
    c('Ethanol wash follows actual precipitation centrifugation',g.steps['ethanol_add'].inputs[0].material,'precipitation_spin_kept')
    c('Retained analytical-QC remainder is a DNA library',g.materials['library_qc_retained'].kind,K.LIBRARY)
    c.report()
if __name__=='__main__':selftest()
