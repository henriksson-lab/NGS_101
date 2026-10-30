"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

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

def qc(g,key,x,label):
    keep=g.product(key+'_retained','Sourceproductionremainder; analyticalaliquotamountnotinvented',K.MOLECULAR);test=g.product(key+'_assay','Separateanalyticalaliquotforthe sourceassay',K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,test));g.mark_terminal(Terminal(test.id,Disposition.QC,'Analyticalmaterialdoesnotbecomethepreparativesequencingtemplate.'));return keep

SOURCE='Published2015chemistrySI originaltitle/authors verified:source_pdf_49222db245fceaa530e3.pdf.txt Methods20–32,62–114,164–181/FigS1A223–227; mainpublisherpreview only. Not2021same-acronym5hmC assay.'
@dataclass(frozen=True)
class SelectedCytosine:
    target:str
    protected:bool=False
    inferred:bool=True
    def __post_init__(self):
        if self.target not in {'5fC','5caC'}:raise ValueError('Selected2015sourcefCorcaCclassrequired')
    def install_protection(self):
        if self.protected:raise ValueError('Protectionalreadyinstalled')
        return replace(self,protected=True)
    @property
    def installed_group(self):return ('O-ethylhydroxylamine-derived protection' if self.target=='5fC' else 'EDC/NHS-activated aminomethylbenzylazide conjugate') if self.protected else None
    def selected_bs_base(self):return 'C' if self.protected else 'U'
    def selected_reference_call(self):return 'C' if self.selected_bs_base()=='C' else 'T'
    def selected_first_copy(self):return 'G' if self.selected_bs_base()=='C' else 'A'

def make_lib(g,key,x,label):
    x=op(g,key+'_kit',x,label,reagent='SourceEpiGnomeMethyl-SeqEpicentre-Illumina; historicaladapter/primer/polymerase internalsunreviewed',product='SourcephysicalDNAsequencingpreparation plus unresolved templates/kit remnants; no borrowed modernkitgeometry')
    x=qc(g,key+'_qc',x,'SourceBioanalyzerDNA1000 separateanalyticalaliquot;productionremainder retained')
    lib=g.product(key+'_library','SourceIlluminacompatible physicalDNA library; fullgeometry unassembled',K.LIBRARY);g.react(key+'_identity','Sourcefinalpreparationidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);return lib

def graph_for_target(target='5fC'):
    if target not in {'5fC','5caC'}:raise ValueError('Actual2015sourcefCorcaCchannelrequired')
    g=ProtocolGraph('DIP-CAB-seq2015 — '+target+' selectedchannel',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','Source100ugsonicatedTdgfl/florTdg−/−mESCgenomicDNA average400bp; biologicalgenotype notassigned',K.MOLECULAR,primary=True)
    # Four parallel physical IP preparations, subsequently combined as printed.
    parts=[g.product('ip_'+str(i),'SourceoneoffourparallelIPallocations;amountnotinvented',K.MOLECULAR) for i in range(4)];g.react('allocate_ip','SourcefourparallelIP',S.SPLIT,((x,R.TEMPLATE),),parts)
    recovered=[]
    for i,x in enumerate(parts):
        k='ip'+str(i);ab=g.input(k+'_antibody','Sourceanti5fC1ulActiveMotif' if target=='5fC' else 'Sourceanti5caC0.5ul antiserumgiftYiZhang',K.MOLECULAR);mix=g.product(k+'_complex','AntibodyDNAcomplexes plus freeDNA/antibody/other companions',K.MOLECULAR);g.react(k+'_bind','SourceIPbuffer4C4h',S.BINDING,((x,R.TEMPLATE),(ab,R.STRUCTURAL)),(mix,));beads=g.input(k+'_beads','Source40ulPREWASHEDproteinGDynabeads;priorwashmechanicsunassembled',K.SELECTION)
        bound=g.product(k+'_bound','Sourcebead-associatedantibodyDNA with nonspecificcarryover',K.MOLECULAR);flow=g.product(k+'_flow','SourcephysicalunboundIPfraction',K.WASTE);g.react(k+'_capture','SourceproteinG4Covernightcapture',S.AFFINITY_CAPTURE,((mix,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.discard(flow,'Sourceenrichedfractionselected,noperfectantibodyspecificityassumed.')
        y=cleanup(g,k+'_wash',bound,'Source3x1mlIPbuffer/3xIPbuffer+150mMNaCl actualwashes')
        y=op(g,k+'_pk',y,'SourcePK55C2h',S.DEPROTEINIZATION,reagent='SourceProteinaseK',catalyst=True,product='IPDNA plus proteolytic fragments/beads/remnants until extraction')
        recovered.append(cleanup(g,k+'_organic',y,'Sourcephenol/chloroform extraction and ethanol precipitation'))
    combined=g.product('combined','Combined fourIPrecoveries,unquantifiednativecytosines and molecularcarryover',K.MOLECULAR);g.react('combine','SourcecombinefourIPpreparations',S.POOL,tuple((x,R.TEMPLATE) for x in recovered),(combined,))
    control=g.product('bs_input','Sourcehalf pooledIPfortraditionalBS',K.MOLECULAR);cab=g.product('cab_input','Sourceotherhalf pooledIPforCAB',K.MOLECULAR);g.react('divide','SourceFigS1Adividetwohalves',S.SPLIT,((combined,R.TEMPLATE),),(control,cab))
    if target=='5fC':cab=op(g,'protect',cab,'Source10mMO-ethylhydroxylamine/100mMMESpH5/37C4h',S.CHEMICAL_CONVERSION,reagent='SourceO-ethylhydroxylamine',product='Selected5fCprotectionadductclasses plus unprotectedDNA/othernativebases/reactioncompanions;noquantitativeyield claim')
    else:
        cab=op(g,'activate',cab,'Source20mMNHS/2mMEDC/75mMMESpH5/37C0.5h',S.CHEMICAL_CONVERSION,reagent='SourceEDC/NHS',product='Selectedactivated5caCclasses plus unreactedDNA/EDC/NHS/reactioncompanions')
        cab=cleanup(g,'exchange',cab,'Sourcebufferexchange100mMphosphatepH7.5/150mMNaCl; exchangedfractionsphysical')
        cab=op(g,'protect',cab,'Source10mM(4-aminomethyl)benzylazide37C1h',S.CHEMICAL_CONVERSION,reagent='Source(4-aminomethyl)benzylazide',product='Selected5caCprotectionconjugateclasses plus unprotectedDNA/nativeotherbases and reagents; no click-biotin capture invented')
    cab=cleanup(g,'protect_cleanup',cab,'SourceZymoDNAClean&Concentrationcleanup')
    control=op(g,'traditional_bs',control,'SourceMethylCode traditionalbisulfite treatment',S.CHEMICAL_CONVERSION,reagent='SourceMethylCodekit; historicalinternals/cyclecountunreviewed',product='SelectedunprotectedC/fC/caC→U DNAclasses plus retained/unconvertednativebases/degradationremnants; exactkitfatesunassembled')
    cab=op(g,'cab_bs',cab,'SourceMethylCodeCABbisulfite;'+('sourceprogramruntwicefor5fC' if target=='5fC' else 'source5caCprogram'),S.CHEMICAL_CONVERSION,reagent='SourceMethylCodekit; no numericalconversion/yieldprediction',product='ProtectedtargetC-retaining classes plus unprotected/deaminatedU-DNA/degradationremnants; chemicaladductnotrawreportedTcall')
    # Source50ng allocation after each bisulfite treatment; no discarded remainder invented.
    libraries=[]
    for key,x in [('bs',control),('cab',cab)]:
        used=g.product(key+'_50ng','Source50ngtreatedDNAintolibrarykit',K.MOLECULAR);rest=g.product(key+'_remainder','IndependentunusedtreatedDNA;amount/fateunassigned',K.MOLECULAR);g.react(key+'_allocate','Source50nglibraryallocation',S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No assumedstorageordiscard.'))
        libraries.append(make_lib(g,key,used,'SourceEpiGnomeadaptertaggingandamplificationtoappropriateconcentration'))
    inp=g.input('matched_input','Sourcegenotype-matchedNONENRICHEDgenomicDNAcontrol, independentstock;notIPflowthrough',K.MOLECULAR);libraries.append(make_lib(g,'input',inp,'Sourceparallelmatchednonenrichedinputlibrary;exactcontrolchemical/libraryrecipeunassembled'))
    others=g.input('other_libraries','Otherphysicalsourcepooledlibraries/biologicalreplicates;membersnotassigned',K.LIBRARY);pool=g.product('pooled_library','SourcepooledphysicalDNA libraries; no covalentmerging',K.LIBRARY);g.react('pool','Sourcepoolsequencinglibraries',S.POOL,tuple((x,R.TEMPLATE) for x in libraries+[others]),(pool,))
    g.sequence(pool,SequencingPlan(limitation='SourceHiScanrapidtwo-laneflowcellSE100/TruSeqSBSv3-HS. Protected-versusBSunconvertedrates/Fishercomparison are computationalcalls,notperfecttargetidentity/chemicalyield. RNA or2021deaminationchemistrynotinherited;exactprimers/index/rungeometryunassembled.'))
    g.add_source_item(SourceItem('pending','Publishedmainfullbody/historicalEpiGnomeMethylCodekitgeometry; directlycitedCAB-seqJACS2013original','MainPDFreturnedHTML,EuPMCXML500,PMCcloudnousablefile. PublishedSIpositive fullselectedrecipe authority; linkedCABcompoundpreparationoriginalunread,notabsent.'))
    return g

def protocol_graph():return graph_for_target()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check()
    for target in ('5fC','5caC'):
        x=SelectedCytosine(target);p=x.install_protection();c(target+'physicalUversusreferenceT',(x.selected_bs_base(),x.selected_reference_call(),x.selected_first_copy()),('U','T','A'));c(target+'protectedclassCbutcopyG',(p.selected_bs_base(),p.selected_reference_call(),p.selected_first_copy()),('C','C','G'))
    g=protocol_graph();h=graph_for_target('5caC');c('ActualfourparallelIPbeforepool',len(g.steps['combine'].inputs),4)
    c('CABandBSareindependenthalves',tuple(g.steps['divide'].outputs),('bs_input','cab_input'))
    c('caCactivationprecedestagafterphysicalexchange',h.steps['protect'].inputs[0].material,'exchange_kept')
    c('Protectiondebrisremovedonlyphysicalcleanup',g.steps['protect_cleanup'].inputs[0].material,'protect_product')
    c('InputcontrolisnotIPflowthrough',g.steps['input_kit'].inputs[0].material,'matched_input')
    c('UnusedtreatedDNAretainedunresolved',g.terminals['cab_remainder'].disposition,Disposition.UNRESOLVED)
    c('Kitgeometrypreventsfinishedmarker',bool(g.publication_problems()),True)
    src=Source(sources_root()/'dip-cab-seq-5fc-5cac__10.1038+cr.2015.5'/'source_pdf_49222db245fceaa530e3.pdf.txt','get_sources --slugs dip-cab-seq-5fc-5cac__10.1038+cr.2015.5')
    if have(c,src,label='publishedchemistrySI'):
        t=src.path.read_text();c('SourceactualNGSlibrarynotSangervalidation','EpiGnome Methyl-Seq kit' in t and 'Single-read 100-bp' in t and 'divided into two halves' in t,True)
    c.report()
if __name__=='__main__':selftest()
