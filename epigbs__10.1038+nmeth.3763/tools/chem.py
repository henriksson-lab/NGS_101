"""Published epiGBS supplement: methylated adapters and nick translation.

The defining main Methods are not acquired. This graph is deliberately private
and incomplete; no GBS constructor or final library recipe is inherited.
"""
from dataclasses import dataclass, replace
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
import illumina as il
from chemdraw import Scene, Segment, revcomp, feature, Construct
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R,
                            StepKind as S, SequencingPlan, SourceItem)
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

SOURCE = 'Published DOI10.1038/nmeth.3763 SI MOESM158, original FigS1p1/p2, S3p4 and S9p11; main Methods unavailable.'
COMMON = 'C' + il.PE_READ2_ARM
# FigS1a A1 and B1 are selected inventory examples, not assigned study samples.
BARCODE = 'AACT'
CODE_A = feature('epigbs.sample.a', 'sample_index', 'fixed', note='FigS1a A1 component example, no sample assignment.')
CODE_B = feature('epigbs.sample.b', 'sample_index', 'fixed', note='FigS1a B1 component example, no sample assignment.')
SOURCE_OLIGO_EXCEPTIONS = {
    'adapter_b': 'Published MOESM158 original FigS1p1 is image-only. Coordinator and independent reader rendered the original, checked all six complete ordering/PCR rows, and verified B top = AACT + AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG; B bottom = CTCGGCATTCCTGCTGAACCGCTCTTCCGATCT + AGTT + TGCA. Only these two failed text matches are exempt; A/PCR strands remain checked. All adapter C residues are5mC and5primeOH retained; A1/B1 are inventory examples, not an assigned sample pair. No synthetic extracted-text twin, absence claim or geometry exemption.'
}

@dataclass(frozen=True)
class SelectedAdapterRepair:
    """Conditional source-selected class, not repair efficiency or an enzyme simulator."""
    ends: tuple[str, str]
    genomic_marks: tuple = ()
    repaired: bool = False

    def __post_init__(self):
        if len(self.ends) != 2 or any(x not in {'A', 'B'} for x in self.ends):
            raise ValueError('Both adapter identities are required')

    def nick_translate_selected(self):
        if self.repaired:
            raise ValueError('Select an unrepaired precursor, not an already repaired class')
        return replace(self, repaired=True)

    @property
    def adapter_cytosine_class(self):
        return '5mC: original oligos and source selected replacement nucleotides'

    @property
    def genomic_methylation(self):
        return self.genomic_marks

def classes():
    # Both orientations and homotypic ends persist; no perfect AB selection.
    return tuple(SelectedAdapterRepair(x) for x in (('A','A'),('A','B'),('B','A'),('B','B')))

def adapter(which):
    sc = Scene()
    if which == 'A':
        sc.strand('top', [Segment('Read1 body',il.TRUSEQ_READ1), Segment('A1 sample index',BARCODE,feature=CODE_A), Segment('PstI cohesive end','TGCA')], mod5='OH')
        sc.anneal('bottom', [Segment('A1 index complement',revcomp(BARCODE),feature=CODE_A), Segment('Read1 complementary body',revcomp(il.TRUSEQ_READ1))], to='top', pair=('Read1 complementary body','Read1 body'), mod5='OH')
        sc.ligation_end('top','PstI cohesive end',kind='cohesive',polarity='3-prime',owner='adapter')
    elif which == 'B':
        sc.strand('top',[Segment('B1 sample index',BARCODE,feature=CODE_B),Segment('common complementary body',revcomp(COMMON))],mod5='OH')
        sc.anneal('bottom',[Segment('common body',COMMON),Segment('B1 index complement',revcomp(BARCODE),feature=CODE_B),Segment('PstI cohesive end','TGCA')],to='top',pair=('common body','common complementary body'),mod5='OH')
        sc.ligation_end('bottom','PstI cohesive end',kind='cohesive',polarity='3-prime',owner='adapter')
    else:
        raise ValueError('Explicit source adapter A or B required')
    return sc

@dataclass(frozen=True)
class MethylatedAdapterInventory:
    name: str
    cytosines: str = 'Every C in both source ordering oligos is 5-methyl-C, including barcode/cohesive-end C'
    ends: str = '5primeOH; source oligos unphosphorylated'

def unresolved_scene(label):
    sc=Scene();sc.strand('population',[Segment('unassembled','['+label+']',placeholder=True)]);return sc

def opaque(g,key,x,label,kind=S.MAGIC):
    catalytic = kind == S.DIGESTION
    stock=g.input(key+'_condition','Source-defined operation; exact conditions in unreviewed main Methods',K.CATALYST if catalytic else K.CONDITION)
    out=g.product(key+'_product',label+'; original molecules, incomplete products and reagent companions retained until an actual workup',K.MOLECULAR)
    kw={'reason':'Opaque source gap: intervening operations, if any, are unreviewed; no cleanup, aliquot allocation or intervening reaction is asserted.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalytic else R.CONDITION)),(out,),**kw)
    return out

def protocol_graph():
    g=ProtocolGraph('epiGBS — source PstI adapter/nick-repair scope',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','Genomic DNA; actual sample/genomic methylation states unassigned, not synthetic adapter methylation',K.MOLECULAR,primary=True)
    x=opaque(g,'digest',x,'Source PstI restriction-fragment entry; Csp6I branch not a separate protocol',S.DIGESTION)
    a=g.input('adapter_a','FigS1 methylated unphosphorylated A adapter; selected A1 inventory example',K.OLIGO,state=adapter('A'),chemical_state=MethylatedAdapterInventory('A'))
    b=g.input('adapter_b','FigS1 methylated unphosphorylated B adapter; selected B1 inventory example',K.OLIGO,state=adapter('B'),chemical_state=MethylatedAdapterInventory('B'))
    lig=g.input('ligation_reagent','Source ligation; detailed reagent/conditions unreviewed',K.CATALYST)
    y=g.product('nicked','Source nicked precursor cohort: AA/AB/BA/BB, free adapters, genomic material and other reaction companions',K.MOLECULAR,chemical_state=classes())
    geometry=AdapterLigationGeometry.known('cohesive',adapter('A'),adapter('B'),unresolved_scene('Full genomic substrates and sealed-side junction products unassembled'),statement='FigS1/S9 establish adapter3primeTGCA and unphosphorylated5primeOH; known adapter end does not assert fully ligated covalent DNA.')
    g.react('ligate','Source ligation leaves a nick at each adapter5primeOH/genomic3primeOH boundary',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(b,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=geometry)
    pol=g.input('pol_i','Source DNA polymerase I and source 5mC nucleotide replacement pool; exact formulation unreviewed',K.CATALYST)
    z=g.product('repaired','Source selected nick-translated adapter classes plus unrepaired/partial DNA, released nucleotides and enzymes; no perfect repair or cleanup assumed',K.MOLECULAR,chemical_state=tuple(c.nick_translate_selected() for c in classes()))
    g.react('nick_translate','FigS9: source selected5prime excision/3prime replacement through adapter end; methylated replacement C',S.EXTENSION,((y,R.TEMPLATE),(pol,R.CATALYST)),(z,))
    x=opaque(g,'pre_conversion_workup',z,'Unreviewed continuation before conversion; physical workup/allocation, if any, unassigned')
    x=opaque(g,'bisulfite',x,'Source bisulfite-treated genomic-fragment preparation; adapter marks and native genomic states remain distinct',S.CHEMICAL_CONVERSION)
    x=opaque(g,'pre_pcr_workup',x,'Unreviewed continuation after conversion; recovery/allocation, if any, unassigned')
    primers=[]
    for key,body in [('primer1',il.NEBNEXT_UNIVERSAL_PRIMER),('primer2',il.PE_PCR_PRIMER_2)]:
        primers.append(g.input(key,'FigS1b historical Illumina PE-PCR oligo; no modern index inferred',K.OLIGO,state=Construct([Segment('source PCR primer',body)])))
    enzyme=g.input('pcr_reagent','Source PCR; polymerase/cycles in unreviewed main',K.CATALYST)
    amplified=g.product('amplified','New PCR descendants, original converted molecules, unextended primers and other reaction companions; no faithful inheritance of synthetic methyl marks asserted',K.MOLECULAR)
    g.react('pcr','Source PE-PCR amplification; full copy/end geometry unresolved',S.PCR,((x,R.TEMPLATE),*((p,R.STRUCTURAL) for p in primers),(enzyme,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(unresolved_scene('Full conversion/PCR geometry unassembled'),reason='PCR copies are not assumed to inherit original5mC chemistry.'))
    x=opaque(g,'finish',amplified,'Unreviewed continuation to acquisition; cleanup/QC/pooling, if any, unassigned')
    g.sequence(x,SequencingPlan(limitation='PublishedSI establishes actual Illumina paired-end/inline barcodes, not complete run chemistry or cycle count. Native strand identity is inferred computationally; read merging/concatenating with10Ns does not create a physical DNA molecule.'))
    g.add_source_item(SourceItem('pending','Defining published full Methods and directly cited GBS component originals','No GBS constructor imported; source SI explicitly attributes adapter design, with child B barcode/methylation/nick repair retained locally. All detailed allocation, formulation, full covalent products, QC and run originals remain unreviewed, not confirmed absent.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();old=SelectedAdapterRepair(('A','B'),((2,'native5mC'),));new=old.nick_translate_selected()
    c('repair does not rewrite native genomic methylation',new.genomic_methylation,old.genomic_methylation)
    c('precursor persists separately from selected descendant',(old.repaired,new.repaired),(False,True))
    c('all ordered adapter classes retained',{x.ends for x in classes()},{('A','A'),('A','B'),('B','A'),('B','B')})
    c.raises('invalid adapter identity rejected',lambda:SelectedAdapterRepair(('A','C')))
    c.raises('already repaired class cannot masquerade as precursor',new.nick_translate_selected)
    for k in ('A','B'):adapter(k) # Scene validates both barcode and fixed-body pairing.
    g=protocol_graph()
    c('repair precedes conversion',g.steps['pre_conversion_workup'].inputs[0].material,'repaired')
    c('post-conversion original population remains PCR template',g.steps['pcr'].inputs[0].material,'pre_pcr_workup_product')
    c('unreviewed geometry prevents finished promotion',bool(g.publication_problems()),True)
    src=Source(sources_root()/'epigbs__10.1038+nmeth.3763'/'nmeth.3763_41592_2016_BFnmeth3763_MOESM158_ESM.pdf.txt','PublishedSI FigS9p11 original rendered; source chemical operation, not absence evidence')
    if have(c,src,label='published SI'):
        t=src.path.read_text();c('source explicitly distinguishes nick repair from bisulfite',all(x in t for x in ['non-phosphorylated adapters','DNA polymerase I','bisulfite']),True)
    c.report()

if __name__=='__main__':selftest()
