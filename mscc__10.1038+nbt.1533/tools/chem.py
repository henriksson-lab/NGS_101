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

from chemdraw import revcomp
from adapter_ligation import reported_adapter_end
SOURCE='AssociatedAMNIHMS103254/PMC3566772 Methods151–167;publishedMOESM19TableS5–8;officialcorrigendumDOI10.1038/nbt0509-485bparagraph642(affiliationonly). Mainfinalversionapplicabilityunreviewed;noBSPPchemistryinherited.'
@dataclass(frozen=True)
class NativeDNAEnd:
    sequence:str
    five_prime:str='P'
    methyl_c:tuple=()
    inferred:bool=True
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT') or self.five_prime not in {'P','OH'} or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='C' for i in self.methyl_c):raise ValueError('ExplicitDNA/endclassandactualnativecytosinemarksrequired')
    def selected_dephosphorylate(self):
        if self.five_prime!='P':raise ValueError('Selectedphosphatasesubstraterequires5P')
        return replace(self,five_prime='OH')
    def selected_unmodified_copy(self):return revcomp(self.sequence)
@dataclass(frozen=True)
class SelectedMmeIProducts:
    top_parts:tuple
    bottom_parts:tuple
    inferred:bool=True
    def rejoin(self):return (''.join(self.top_parts),''.join(self.bottom_parts))
def selected_mmei_cut(top,bottom,top_cut,bottom_cut):
    if not top or set(top)-set('ACGT') or bottom!=revcomp(top) or any(type(i)is not int or not 0<i<len(top) for i in (top_cut,bottom_cut)):raise ValueError('Explicitduplexandbothcaller-selected5to3strandcutcoordinatesrequired;noMmeIoffsetdefaults')
    return SelectedMmeIProducts((top[:top_cut],top[top_cut:]),(bottom[:bottom_cut],bottom[bottom_cut:]))
def source_adapter(which):
    if which not in {'a','b'}:raise ValueError('SourceAorBadapterrequired')
    return reported_adapter_end('SourceMSCCadapter'+which.upper(),kind='cohesive',polarity='5-prime' if which=='a' else '3-prime',terminal_sequence='CG' if which=='a' else 'NN',statement='Sourceendsequence/polarityknown;completecustomadapterbody,modifications,junctiongeometryunreviewed.')
def graph_for_channel(channel='hpaii'):
    if channel not in {'hpaii','mspi','inverse'}:raise ValueError('ThreeexplicitpublishedMSCCchannelsonly')
    g=ProtocolGraph('MSCC — '+channel+'sourcechannel',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','Source2ugPGP1lymphocytegenomicDNA;unboundednativecytosines/marksunassigned',K.MOLECULAR,primary=True)
    if channel=='inverse':
        x=op(g,'hpa_initial',x,'SourceHpaII37C2h/65C20min',S.DIGESTION,reagent='Source20UHpaII/NEBuffer1',catalyst=True,product='SelectedHpaIIcutsplusuncut/methylatedDNA,enzymesandotherfragments;nospecificity/yieldprediction')
        x=op(g,'mask',x,'SourceAntarcticphosphatase37C1h/65C15min',S.DEPHOSPHORYLATION,reagent='Source10UAntarcticphosphatase/buffer',catalyst=True,product='SelectedHpaII-generated5OHendsplus5P/unprocessedDNA/phosphateandenzymes;notalltemplatesassertedblocked')
        x=cleanup(g,'mask_cleanup',x,'Sourcephenol/chloroformfollowedbyethanolprecipitation')
    x=op(g,'primary_digest',x,'Source'+('HpaII20U/NEBuffer1' if channel=='hpaii' else 'MspI40U/NEBuffer2')+'37C2h/65C20min',S.DIGESTION,reagent='Sourceend-specificrestrictionenzyme;nativeCsequenceisnotbisulfiteconverted',catalyst=True,product='SourcecutgenomicDNAanduncutDNAplusenzymes/remnants;inversechannel5OHmaskedpriorendsremainphysical')
    for which,key,label in [('a','ligate_a','SourceadapterA/T4ligation16C4h/65C15min')]:
        geom=source_adapter(which);adapter=g.input('adapter_'+which,'SourcecustomadapterAwith5CGendandMmeIrecognitionsite;fullbody/modificationsunreviewed',K.OLIGO,state=geom.scenes[0]);cat=g.input(key+'_ligase','SourceT4ligase/ATP',K.CATALYST);y=g.product(key+'_product','AdapterA-genomicDNAclassesplusfreeAadapter/otherDNA/remnants;fullligationjunctionunassembled',K.MOLECULAR)
        g.react(key,label,S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=geom);x=y
    x=cleanup(g,'ethanol_a',x,'SourceethanolprecipitationafteradapterAligation')
    x=op(g,'bst',x,'SourceBstfragment/dNTP50C20min/85C20min',S.EXTENSION,reagent='Source8UBstfragment/Thermopol/dNTP;exactboundpriminggeometryunreviewed',catalyst=True,product='Selectedadapter/genomiccopiesplusoriginalnativeDNA/marks,unextendedadaptersandreagents;copyingnotnativeDNAconversion')
    x=cleanup(g,'ethanol_bst',x,'SourceethanolprecipitationafterBstextension')
    x=op(g,'mmei',x,'SourceMmeI/SAM37C2h/80C20min',S.DIGESTION,reagent='Source2UMmeI/50uMSAM/NEBuffer4',catalyst=True,product='18–19bpgenomic-tag-bearingselectedDNAplusbothdistalMmeIdaughters,uncutDNAandreagents;noinstantdistaldiscard')
    geom=source_adapter('b');adapter=g.input('adapter_b','SourcecustomadapterBwith3NNcohesiveend;notUMI/fullbodyunreviewed',K.OLIGO,state=geom.scenes[0]);cat=g.input('b_ligase','SourceT4ligase/ATP',K.CATALYST);y=g.product('ligate_b_product','SourceA/tag/BlibraryclassesplusdistalDNA/freeadapters/reagentsuntilPAGE',K.MOLECULAR)
    g.react('ligate_b','SourceadapterB16C4h/65C15min',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=geom)
    x=cleanup(g,'page',y,'Source6%non-denaturingTBE-PAGE~140bpbandpurification')
    if channel=='hpaii':
        selected=g.product('pcr_input','Sourceapproximately80%gel-purifiedDNAforPCR',K.MOLECULAR);rest=g.product('gel_remainder','Sourceotherapproximately20%gel-purifiedDNA, laterfateunassigned',K.MOLECULAR);g.react('allocate_pcr','Source80%PCRallocationandindependentremainder',S.SPLIT,((x,R.TEMPLATE),),(selected,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Notinferredstoredordiscarded.'))
        x=unknown_pcr(g,'pcr',selected,'SourceiProof8cycles98C10s/67C15s/72C15s;98C30sinitial/72C5minfinal','SourceiProof/HFbuffer/dNTP/customA/Bprimersunreviewed');x=cleanup(g,'qiaquick',x,'SourceQIAquickPCRcleanup')
    lib=g.product('library','SourceMSCCphysicalsequencinglibrary;'+('PCR-amplified' if channel=='hpaii' else 'noPCRaftergelpurification')+';fulladapter/readgeometryunassembled',K.LIBRARY);g.react('library_identity','Sourcefinalpreparationidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='ActualIlluminaGenomeAnalyzerMSCC/HpaII,MspIcontrolandinverselibraries verifiedTableS7. 18–19bpgenomictagsnotrawreadcyclecount. Restriction-tagcounts/methylationestimationcomputational;nativeCnotphysicallyrewrittenasT,noperfectenzymespecificityorlinearcalibrationguaranteed.'))
    g.add_source_item(SourceItem('pending','Finalpublishedmainapplicability/fullcustomadapter/primer/end/readgeometry','Threeoriginalcapclosed:AMmain,publishedSI,officialaffiliationcorrigendumread.61MBMSCCdatatablewasnotsequence-scanned;S6headeridentifieddatacolumns,nosequencedisclosureabsenceclaimed. BSPPmethodinsamepaperexcludedfromMSCCrecipe.'))
    return g
def protocol_graph():return graph_for_channel()
def selftest():
    from checks import Check
    c=Check();a=NativeDNAEnd('ACCGT','P',(1,));b=a.selected_dephosphorylate();c('PhosphataseremovesendPwithoutDNAorcytosinemarkloss',(b.sequence,b.methyl_c,b.five_prime),('ACCGT',(1,),'OH'))
    c('newunmodifiedDNAcopyiscomplementarynotnativeconversion',a.selected_unmodified_copy(),'ACGGT')
    t='ACGTGCTAACGTT';products=selected_mmei_cut(t,revcomp(t),5,9);c('selectedMmeIcutretainsbothwholestrands',products.rejoin(),(t,revcomp(t)))
    c('selectedcutreferenceinferred',products.inferred,True)
    try:selected_mmei_cut(t,t,5,9);bad=False
    except ValueError:bad=True
    c('invalidduplexreferencecannotcreateproduct',bad,True)
    g=protocol_graph();m=graph_for_channel('mspi');i=graph_for_channel('inverse');c('MspIcontrolhasnoaftergelPCR','pcr' not in m.steps,True)
    c('inversebranchmaskingprecedesMspI',i.steps['primary_digest'].inputs[0].material,'mask_cleanup_kept')
    c('bothMmeIdaughtersretaintoPAGE',g.steps['ligate_b'].inputs[0].material,'mmei_product')
    c('unallocatedPCRmassretained',g.terminals['gel_remainder'].disposition,Disposition.UNRESOLVED)
    c('unassembledbody/readgeometrypreventsfinishedmarker',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
