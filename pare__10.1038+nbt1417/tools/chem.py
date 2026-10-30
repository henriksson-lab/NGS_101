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
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    y=g.product(key+'_product',product or label+'; all physical reaction companions retained until separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    kwargs={}
    if kind==S.MAGIC:kwargs['reason']='Source operation established; exact internal chemistry/geometry not assembled from unreviewed historical originals.'
    if kind==S.PCR:kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exact source template/primer endpoint geometry unassembled, no ordinary complete library asserted.')
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-known ligation; full substrate/end geometry unassembled]'),),'Unreviewed geometry; not a claim of unavailable disclosure.')
    inputs=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)]
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
    if kind==S.LIGATION:
        adapter=g.input(key+'_adapter','Source ligating molecular adapter; exact historical oligo/end inventory unassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter geometry unassembled]'));inputs.append((adapter,R.STRUCTURAL))
    if kind==S.PCR:
        for arm in ('f','r'):
            primer=g.input(key+'_primer_'+arm,'Source PCR primer family; exact oligo/endpoint unassembled',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL))
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
SOURCE='PublishedNatureSupplementaryMOESM9_ESM.pdf originalFig1p1rendered +caption65–72; Fig2selectedsignatureorientationoverview. DefiningmainDOI10.1038/nbt1417 unavailable; nofullprotocolreview.'
@dataclass(frozen=True)
class ExplicitDuplexCut:
    top:str
    top_cut:int
    bottom_cut:int
    def __post_init__(self):
        if not self.top or set(self.top)-set('ACGT') or any(type(i)!=int or not 0<i<len(self.top) for i in (self.top_cut,self.bottom_cut)):raise ValueError('ExplicitDNAandcallerselectedinternalcutcoordinatesrequired; noMmeIdefaults')
        if not self.verify_conservation():raise ValueError('BothDNAdaughterstrandsmustberetained')
    def products(self):
        return ((self.top[:self.top_cut],revcomp(self.top[:self.bottom_cut])),(self.top[self.top_cut:],revcomp(self.top[self.bottom_cut:])))
    def verify_conservation(self):
        left,right=self.products()
        return left[0]+right[0]==self.top and right[1]+left[1]==revcomp(self.top)

@dataclass(frozen=True)
class NativeEndRNA:
    sequence:str
    five_end:str='P'
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU') or self.five_end not in {'P','OH','cap'}:raise ValueError('ExplicitRNAandselected5endclassrequired')
    def ligate_source_five_adapter(self,adapter_rna):
        if self.five_end!='P' or not adapter_rna or set(adapter_rna)-set('ACGU'):raise ValueError('SelectedFig1ligationrequiresRNA5PandsourceRNAadapter3OH; noend-repairinferred')
        return adapter_rna+self.sequence

def protocol_graph():
    g=ProtocolGraph('PARE — sourceSupplementaryFig1end-taglibrary',authority=SOURCE,representation='published_instance')
    x=g.input('rna','SourceRNApreparation, polyadenylateddecayproducts included; exactinputworkupunreviewed',K.MOLECULAR,primary=True)
    kept=g.product('poly_a_rna','SourceisolatedpolyadenylatedRNApopulation withselectednative5P3OHclasses; methodandpurityunassembled',K.MOLECULAR);other=g.product('poly_a_other','Otherphysicalisolationfractions, exactsourcefateunreviewed',K.MOLECULAR)
    g.react('poly_a_isolate','SourceFig1step1polyadenylatedRNAisolation',S.MAGIC,((x,R.TEMPLATE),),(kept,other),reason='PublishedsourceestablishespolyAisolation; actualselectionreagents/fractiongeometryunreviewed, nobeadorcaptureassumption.');g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Sourceotherisolationmaterial; nostatedstorage/discardinselectedoriginal.'))
    a=g.input('five_rna_adapter','SourceFig1RNAadapterforRNA5end; bases/endgeometryunassembled',K.ADAPTER,state=scene('five_rna_adapter','[sourceRNA5primeadapterbodyandligationendsunassembled]'));cat=g.input('rna_ligase','Source5adapterligationcatalyticcomponent, exactidentityunreviewed',K.CATALYST);y=g.product('rna_ligated','Source5adapter-RNAclasses plusunligatedRNA/freeadapter/reagents retained',K.MOLECULAR)
    g.react('ligate_five','SourceFig1step2RNA5adapterligation',S.LIGATION,((kept,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(scene('ends','[nativeRNA5P/sourceadapter3OHpopulation; fulljunctionunassembled]'),),'SourceFigure1endclasses; noexactoligoorligationgeometryassembled.'))
    x=op(g,'rt',y,'SourceFig1step3oligo(dT)/adapterfirststrandRT',reagent='Sourceoligo(dT)/adapterprimerandRTcomponentsunreviewed',product='cDNAplusoriginaladapterRNA/unextendedprimers/remnants; sourceVT21printedprimerinventory, completeorientation/bodyunassembled')
    x=op(g,'short_pcr',x,'SourceFig1step4shortPCRgeneratesdsDNAtemplates',S.PCR,reagent='SourcepreMmeIPCRprimers/polymeraseunreviewed',product='Sourceadapter-bearingdsDNAcopies plusoriginalRNA/cDNA/reagentcompanions; fullcopygeometryunassembled')
    x=op(g,'mme_i',x,'SourceFig1step5MmeIdigestion',S.DIGESTION,reagent='SourceMmeI, exactrecognitionposition/cutcoordinatesunreviewed',product='SourceMmeIreactioncontainsboth5tag-bearingDNAanddistalpolyA-endDNAproducts plusuncutDNA/othercompanions; nophysicaldistalfragmentremovalasserted')
    x=op(g,'dna_adapter',x,'SourceFig1step6degenerate-overhangDNAadapterligation',S.LIGATION,reagent='Source3DNAadapterwithdegenerateoverhang; fullbases/polarityunreviewed',product='Adapter-ligatedDNAclasses andbothMmeIdaughterclasses/unligatedadapters/reagents remaininsourcereaction')
    x=op(g,'final_pcr',x,'SourceFig1step7taglibraryPCR',S.PCR,reagent='SourcefinalPCRprimer/polymeraseinventoryunreviewed',product='Sourceamplifiedtaglibraryclasses plusoriginaltemplates/distalDNAdaughters/reactioncompanions; noperfectamplificationselectionassumed')
    x=clean(g,'page',x,'SourcefinalPAGEpurificationoftaglibrary; exactband/quantitativepartitionunreviewed')
    x=qc_split(g,'qc',x,'SourceFig1QCbranchwithretainedSBSproductionpreparation; exactassay/aliquotunreviewed')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='PublishedSupplementaryFig1SBSsequencing; exactinstrument/run/index/primergeometryunreviewed. Mappedend-tagsreportRNAcleavagepositions, notphysicalmiRNA-targetpairedmolecules.'))
    g.add_source_item(SourceItem('pending','Definingpublishedmainandactualdetailedlibrary/oligoreagentsources','MainDOIroutenoopenfulltext. PublishedSIactualworkflowread; sourceadapterends/cutcoordinates/fullRT/PCR/readgeometryunassembled, noabsenceclaim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();cut=ExplicitDuplexCut('ACGTACGT',6,4);left,right=cut.products()
    c('ExplicitstaggeredcutretainsbothDNAproducts',cut.verify_conservation(),True)
    c('Bottomstrandconservationusesphysical5to3order',right[1]+left[1],revcomp(cut.top))
    c('Native5PligationpreservesRNAsequence',NativeEndRNA('ACGU').ligate_source_five_adapter('GGA'),'GGAACGU')
    try:NativeEndRNA('ACGU','OH').ligate_source_five_adapter('GGA');bad=False
    except ValueError:bad=True
    c('CannotinventRNA5endphosphorylation',bad,True)
    g=protocol_graph();c('MmeIdistalproductsnotdiscardedbeforeligation','bothMmeIdaughterclasses' in g.materials['dna_adapter_product'].label,True)
    c('QCassaydoesnotfeedSBSproduction','qc_retained' in g.materials and g.materials['qc_assay'].kind==K.QC,True)
    p=Source(sources_root()/'pare__10.1038+nbt1417'/'nbt1417_41587_2008_BFnbt1417_MOESM9_ESM.pdf.txt','get_sources --slugs pare__10.1038+nbt1417; Fig1renderinspected')
    if have(c,p,label='publishedsupplement'):
        t=p.path.read_text();c('ActualcaptionestablishesMmeIthenoverhangadapterligation','digested with MmeI' in t and 'overhang region' in t,True)
    c.report()
if __name__=='__main__':selftest()
