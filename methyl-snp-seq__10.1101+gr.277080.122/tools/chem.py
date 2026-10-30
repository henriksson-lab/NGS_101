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
SOURCE='Yan et al.2022 DOI10.1101/gr.277080.122 published PMC9808626 mainMethods148–164; SupplementalMaterials protocolsteps1–7, originalSIp8 andmainp2088 rendered.'
NAME_BOUNDARY='Mainp2088 printsEndoVI(M0304); SIprotocolp8 printsEndoIV(M0304). SelectedSI name retained; manufacturer original unreviewed. Main adapterE7536 versusSI E7165/E7140 not silently reconciled.'
HAIRPIN='CCACGACGACGACGACGAGCGTTAGGCTCGTCGTCGTCGTCGUGGT'
@dataclass(frozen=True)
class UracilHairpin:
    sequence:str=HAIRPIN
    phosphorylated:bool=True
    def __post_init__(self):
        if self.sequence.count('U')!=1 or set(self.sequence)-set('ACGTU') or not self.phosphorylated:raise ValueError('selected source5phosphate single-uracilDNAhairpin required')
    def scene(self):
        sc=Scene();sc.strand('hairpin',[Segment('hairpin',self.sequence,note='DNA oligo uracil; not evidence ofRNA sugar. Intramolecular stem/loop and ligation junction unassembled.')],mod5='phosphate');return sc
@dataclass(frozen=True)
class OriginalCopyPair:
    original:str
    biological_methyl_c:tuple=()
    def __post_init__(self):
        if not self.original or set(self.original)-set('ACGT'):raise ValueError('explicit selectedoriginal5to3DNA required')
        if len(set(self.biological_methyl_c))!=len(self.biological_methyl_c) or any(type(i)!=int or not 0<=i<len(self.original) or self.original[i]!='C' for i in self.biological_methyl_c):raise ValueError('selected biological5mC positions must beC')
    @property
    def synthetic_copy(self):return revcomp(self.original)
    @property
    def synthetic_methyl_c(self):return tuple(i for i,b in enumerate(self.synthetic_copy) if b=='C')
    def converted_original(self):return ''.join('U' if b=='C' and i not in self.biological_methyl_c else b for i,b in enumerate(self.original))
    def converted_copy(self):return self.synthetic_copy
    def complementary_dna_descendant(self):return revcomp(self.converted_original().replace('U','T'))
    def scene(self):
        sc=Scene();sc.strand('selected_original',[Segment('original',self.original,note='Original biological5mC coordinates:'+str(self.biological_methyl_c))]);sc.strand('selected_protected_copy',[Segment('synthetic',self.synthetic_copy,note='Every copiedC uses5-mdCTP, not original biological methylation')]);sc.strand('companions',[Segment('remnants','[displaced DNA products, uncleaved/other products, enzymes and source reaction companions retained until source cleanup]',placeholder=True)]);return sc

def graph_for_scope(targeted=False):
    g=ProtocolGraph('Methyl-SNP-seq'+(' — exome branch' if targeted else ' — human WGS'),authority=SOURCE,representation='published_instance')
    x=g.input('genomic_dna','SourceNA12878 human genomicDNA; selectedWGS250bp/exome400–500bp sonication',K.MOLECULAR,primary=True)
    spike=g.input('lambda','Sourceunmethylatedlambda DNA conversioncontrol',K.MOLECULAR);mixed=g.product('spiked_dna','Original human+lambdaDNA mixture',K.MOLECULAR);g.react('spike','Sourceunmethylatedlambda spike-in',S.GENERIC,((x,R.TEMPLATE),(spike,R.SPIKE_IN)),(mixed,))
    x=op(g,'sonicate',mixed,'SourceCovarisS2 sonication',S.FRAGMENTATION,reagent='SourceCovaris physicalsonication')
    x=op(g,'end_repair',x,'SourceNEBUltraIIE7546 endrepair/dA20C30min65C30min',reagent='HistoricalE7546 module; exactendchemistry unassembled')
    x=clean(g,'end_cleanup',x,'SourceZymoOligoClean80ntprotocol')
    h=UracilHairpin();hp=g.input('hairpin','Source5phosphate uracilDNAhairpin, full written sequence retained',K.ADAPTER,state=h.scene(),chemical_state=h);cat=g.input('hairpin_ligase','SourceNEBM0367ligation mix20C1h',K.CATALYST);y=g.product('dumbbell_mix','Two-hairpin closed dumbbellclass plus incompletely ligatedDNA/freehairpins/remnants; closuregeometry unassembled',K.MOLECULAR)
    g.react('ligate_hairpin','Sourcehairpinligation atboth fragmentends',S.LIGATION,((x,R.TEMPLATE),(hp,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(h.scene(),),'Sourcehairpin sequence/closure stated; exactend geometry not assembled.'))
    x=clean(g,'ligation_cleanup',y,'Source0.9xSPRI cleanup')
    x=op(g,'exonucleases',x,'SourceExoIII+ExoVII37C1h digest incompletely protected species',S.DIGESTION,reagent='SourceExoIIIM0206 andExoVIIM0379',product='Closed selectedDNAclass plus digestedDNA products and any incompletely digested material; no perfect enrichment asserted')
    x=clean(g,'exo_cleanup',x,'SourceZymo80ntcleanup afterexonuclease')
    x=op(g,'nick',x,'SourceSI UDG/EndoIV37C15min; printednameboundary retained',S.DIGESTION,reagent='SourceUDGM0280/printedSI EndoIVM0304; '+NAME_BOUNDARY,product='Uracil excision/incision-containing DNA classes, base/cleavageproducts and enzymes retained')
    x=op(g,'nick_translate',x,'SourceDNAPolI/5-mdCTP20C30min added directly without purification',S.EXTENSION,reagent='SourceDNAPolIM0209/dATP/dGTP/dTTP/5-mdCTP',product='Translation/replacement reaction mixture; selected original/protectedcopy class shown with displacedDNA products. Opposingnick encounter creates source half-dumbbell classes; fullcut geometry unresolved')
    x=op(g,'a_tail',x,'SourceTaq/dATP65C30min added directly without purification',S.TAILING,reagent='SourceTaqM0267/dATP')
    x=clean(g,'tail_cleanup',x,'SourceZymo80ntcleanup afterTaq')
    x=op(g,'y_adapter',x,'SourcemethylatedYadapterligation20C1h',S.LIGATION,reagent='SourceSI methylatedadapterE7165/E7140 andM0367; mainE7536 identityunreconciled',product='Y-adapter/one-hairpin linked selected DNA plus remaining source reaction classes; exact adapterjunction unassembled')
    x=clean(g,'y_cleanup',x,'Source0.8xSPRI cleanup')
    x=op(g,'bisulfite',x,'SourceFastBisulfite95C25min, modifiedhistoricalAbcamkit',reagent='Source mainab117127 / SIprintedab11727 kit identifiers retained, vendorversionunreviewed',product='Selected original unprotectedC→U, biological5mC retained; every syntheticcopyC protectedby5mC. Other chemistry/kit recovery unassembled')
    x=op(g,'pcr',x,'SourceQ5UM0597 PCR; SI8cycles',S.PCR,reagent='SourceQ5U andunassembled indexprimer/template geometry',product='PCR descendants plus original convertedtemplates/primers/remnants; C→U physicaltemplate andTreportedDNA are distinct')
    if targeted:
        x=op(g,'target_capture',x,'Source200–300ngindexedlibrary pool thenTwistHumanCoreExome102025 enrichment',reagent='HistoricalTwistcapturekit; actual molecularselection/primergeometry unreviewed')
        x=op(g,'target_pcr',x,'SourceNEBNextQ5M0544+E7603 quantprimer mix enrichmentPCR',S.PCR,reagent='SourceQ5M0544; distinct from initialQ5U')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='SourcehumanNovaSeq6000 PE100 (bacterialcomparison NextSeq550PE150); exactlinkedread/hairpin/cyclegeometry unassembled. Cycle100trimming andjointC/Tvariant/methylation calls are computational. '+NAME_BOUNDARY))
    g.add_source_item(SourceItem('pending','HistoricalNEB adapter/endprep/UDG-EndoIV originals, AbcamFastBisulfiteguide andoptionalTwistcaptureguide','Sourceidentified vendor reagents not independentlyread within cap; full hairpin/Yprimer/read geometry unresolved; no absence claim.'))
    return g

def protocol_graph():return graph_for_scope()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=OriginalCopyPair('ACGCCATG',(3,))
    c('originalphysicalCconversion preserves selectedbiological5mC',p.converted_original(),'AUGCUATG')
    c('syntheticcopyCprotected independent oforiginalmethylpattern',p.converted_copy(),revcomp('ACGCCATG'))
    c('changingoriginalbiologicalmark doesnotchangepreservedsequencecopy',OriginalCopyPair(p.original,()).synthetic_copy,p.synthetic_copy)
    c('copiedDNAreadout containsT whilephysicaloriginalcontainsU','U' in p.converted_original() and 'U' not in p.complementary_dna_descendant(),True)
    g=protocol_graph();c('no cleanup insertedbetweenincisionPolIandTaq',tuple(g.steps).index('a_tail')-tuple(g.steps).index('nick'),2)
    c('targetedbranch addsdistinctpostcapturePCR','target_pcr' in graph_for_scope(True).steps and 'target_pcr' not in g.steps,True)
    q=Source(sources_root()/'methyl-snp-seq__10.1101+gr.277080.122'/'gr.277080.122_PMC9808626.xml.txt','get_sources --slugs methyl-snp-seq__10.1101+gr.277080.122')
    if have(c,q,label='publishedhairpin'):
        c('publishedfullhairpinuracilsequence',HAIRPIN in q.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
