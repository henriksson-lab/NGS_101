"""CPT-seq: protein-held contiguity, two indexing tiers, source-limited PCR."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp,workflow_from_protocol_graph,workflow_panel,panel
from page import head,info
from pcr import PCRGeometry
from protocol_model import role
from protocol_graph import (ProtocolGraph,InputRole as R,MaterialKind as K,StepKind as S,
    SequencingPlan,SourceItem,Terminal,Disposition)

PAGE_FILE='cpt-seq.html'
SOURCE='Amini2014 DOI10.1038/ng.3119; associated accepted manuscript Methods116–130 and publisher corrected SupplementaryFigures3–4/Table4 (20May2015)'
T5=('TATAGCCT','ATAGAGGC','CCTATCCT','GGCTCTGA','AGGCGAAG','TAATCTTA','CAGGACGT','GTACTGAC')
T7=('CGAGTAAT','TCTCCGGA','AATGAGCG','GGAATCTC','TTCTGAAT','ACGAATTC','AGCTTCAG','GCGCATTA','CATAGCCG','TTCGCGGA','GCGCGAGA','CTATCGCT')
CONNECTOR_A=nx.S5+'TCCACGC'
CONNECTOR_B=nx.S7+'CTGTCCCTGTCC'
EXTENDED_A='GCGATCGAGGACGGC'+nx.ME
EXTENDED_B='CACCGTCTCCGCCTC'+nx.ME
KNOWN_PRIMERS=(sp.custom('Read 1','Extended ME A',EXTENDED_A,SOURCE),
               sp.custom('Read 2','Extended ME B',EXTENDED_B,SOURCE),
               sp.custom('Index 1 (i7)','Extended ME B reverse complement',revcomp(EXTENDED_B),SOURCE))


def identity(part):
    return feature('cpt_partition_'+part.lower(),'linked_read_barcode','combinatorial',group='cpt_partition',part=part,
                   note='Four-part virtual/physical compartment identity; not a UMI or a guarantee of one molecule per partition')


def transposon(side,index):
    choices=T5 if side=='A' else T7 if side=='B' else ()
    if index not in choices:raise ValueError('Choose a source Table4 transposon member on the specified side')
    connector,extension=(CONNECTOR_A,EXTENDED_A) if side=='A' else (CONNECTOR_B,EXTENDED_B)
    return Construct((Segment('universal connector '+side,connector),
                      Segment('transposon '+side+' index',index,'cbc',feature=identity('T5' if side=='A' else 'T7')),
                      Segment('extended mosaic end '+side,extension)),name='Source transposon '+side)


def adapters():
    sc=Scene()
    for side,index in [('A',T5[0]),('B',T7[0])]:
        con=transposon(side,index)
        prefix=list(con)[:-1]+[Segment('extension flank '+side,con.segments[-1].top[:-len(nx.ME)]),Segment('ME '+side,nx.ME)]
        sc.strand(side,prefix,label='corrected Table4 '+side+' member1')
        sc.anneal('complement '+side,[Segment('pMENTS '+side,nx.ME_RC)],to=side,
                  pair=('pMENTS '+side,'ME '+side),mod5='5Phos',label='non-transferred ME oligo')
    return sc


@dataclass(frozen=True)
class ContiguityState:
    """Main-branch ordering, not a universal physical constraint or yield model.

    The source's pre-dilution SDS control is chemically possible but outside this
    represented post-dilution branch. No intact covalent backbone is implied.
    """
    associated: bool
    stage: str

    def __post_init__(self):
        if self.stage not in {'transposed','quenched','pooled','physical_well','released'}:
            raise ValueError('Unknown represented CPT stage')
        if type(self.associated) is not bool or self.associated != (self.stage!='released'):
            raise ValueError('Association must agree with the represented CPT stage')

    def advance(self,stage):
        permitted={'transposed':'quenched','quenched':'pooled','pooled':'physical_well'}
        if permitted.get(self.stage)!=stage:
            raise ValueError('CPT main branch requires quench,then pool,then physical allocation')
        return ContiguityState(True,stage)

    def dissociate(self):
        if self.stage!='physical_well':
            raise ValueError('Represented main branch requires physical allocation before SDS;pre-dilution control is outside this route')
        return ContiguityState(False,'released')

    def scene(self):
        sc=Scene()
        for key in ['first','second']:
            sc.strand(key,[role('unassigned indexed genomic fragment '+key)],label='representative fragment')
        sc.note('first',('Source-declared association through bound transposase; not a covalent join between the drawn fragments.'
                         if self.associated else 'Released fragments; association no longer preserves physical linkage.'))
        sc.note('second',self.stage+'; fragment count, cut positions and retention probabilities unassigned')
        return sc


def final_reference():
    """Conditional Figure3 part order, with no claimed exact PCR-primer geometry."""
    left=transposon('A',T5[0]);right=transposon('B',T7[0])
    parts=[Segment('P5',il.P5,'p5',inferred=True),
           Segment('PCR i5','N'*8,'cbc',placeholder=True,feature=identity('P5'),inferred=True)]
    parts.extend(Segment(s.name,s.top,s.tag,feature=s.feature,inferred=True) for s in left)
    parts.append(role('genomic insert of unassigned sequence and length',inferred=True))
    parts.extend((Segment('extended mosaic end B reverse complement',revcomp(EXTENDED_B),inferred=True),
                  Segment('transposon B index reverse complement',revcomp(T7[0]),'cbc',feature=identity('T7'),inferred=True),
                  Segment('connector B reverse complement',revcomp(CONNECTOR_B),inferred=True),
                  Segment('PCR i7 reverse complement','N'*8,'cbc',placeholder=True,feature=identity('P7'),inferred=True),
                  Segment('P7 reverse complement',il.P7_RC,'p7',inferred=True)))
    return Construct(parts,name='INFERRED — Figure3 A1/B1 assembly; physical-well index pair unassigned')


def partial_pcr_primer(side):
    if side not in {'i5','i7'}:raise ValueError('Unknown source PCR-primer side')
    prefix,part=(il.P5,'P5') if side=='i5' else (il.P7,'P7')
    return Construct((Segment('Figure3 flowcell arm',prefix,inferred=True),
                      Segment('unassigned physical-well index','N'*8,placeholder=True,feature=identity(part)),
                      role('unassembled historical three-prime annealing arm')),
                     name='Partial Figure3 '+side+' PCR-primer pattern')


def protocol_graph():
    g=ProtocolGraph('CPT-seq — representative transposition/physical-well route',authority=SOURCE)
    dna=g.input('hmw_dna','High-molecular-weight human DNA; one of96 transposition reactions',K.BIOLOGICAL,primary=True)
    stock=g.input('adapter_stock','Representative corrected Table4 A1/B1 oligos annealed to pMENTS',K.ADAPTER,state=adapters())
    protein=g.input('tn5_protein','EZ-Tn5 protein; source adapter loading37C1h',K.MOLECULAR)
    enzyme=g.input('tn5_activity','Catalytic Tn5 activity in the source-loaded complexes',K.CATALYST)
    loaded=g.product('loaded','Loaded transposome stock; exact dimer composition unassigned',K.MOLECULAR,state=adapters())
    g.react('load','Anneal adapters and load indexed transposome stocks',S.BINDING,((stock,R.TEMPLATE),(protein,R.STRUCTURAL)),(loaded,))
    condition=g.input('tagment_condition','1ng DNA/20uL TD buffer;2.5pmol complex;55C10min',K.CONDITION)
    complex_state=ContiguityState(True,'transposed')
    complexes=g.product('complexes','Transposed DNA retaining protein-held contiguity',K.MOLECULAR,state=complex_state.scene(),chemical_state=complex_state)
    g.react('transpose','Transpose before pooling or dilution; bound transposase remains',S.TAGMENTATION,
            ((dna,R.TEMPLATE),(loaded,R.STRUCTURAL),(enzyme,R.CATALYST),(condition,R.CONDITION)),(complexes,))
    edta=g.input('edta','Equal volume40mM EDTA pH8;37C15min; no SDS yet',K.CONDITION)
    quenched_state=complex_state.advance('quenched')
    stopped=g.product('stopped','EDTA-quenched associated complexes',K.MOLECULAR,state=quenched_state.scene(),chemical_state=quenched_state)
    g.react('quench','Quench transposition without assuming protein dissociation',S.QUENCHING,((complexes,R.TEMPLATE),(edta,R.CONDITION)),(stopped,))
    others=g.input('other_transposition_wells','Corresponding products of other95 source transposition reactions',K.MOLECULAR)
    pooled_state=quenched_state.advance('pooled')
    pooled=g.product('pooled','Pool of96 transposition-index pairs; contiguity preserved operationally',K.MOLECULAR,state=pooled_state.scene(),chemical_state=pooled_state)
    g.react('pool_transposition','Pool source transposition wells before dilution',S.POOL,((stopped,R.TEMPLATE),(others,R.TEMPLATE)),(pooled,))
    allocation=g.input('allocation','25pg/uL pool diluted to1pg/uL TE;10pg per physical well;96wells',K.CONDITION)
    allocated_state=pooled_state.advance('physical_well')
    aliquot=g.product('aliquot','Representative physical-well aliquot',K.MOLECULAR,state=allocated_state.scene(),chemical_state=allocated_state)
    reserve=g.product('other_aliquots','Other physical-well aliquots; handled by the same source procedure',K.MOLECULAR,state=allocated_state.scene())
    g.react('allocate_wells','Dilute and allocate while transposase association remains',S.SPLIT,((pooled,R.TEMPLATE),(allocation,R.CONDITION)),(aliquot,reserve))
    g.mark_terminal(Terminal(reserve.id,Disposition.STORED,'Other wells are outside this representative molecular route; their amplified libraries enter the final pool as declared input.'))
    blocker=g.input('i7_blocker','Nextera i7 primer already in physical wells; also serves as surface blocker',K.OLIGO,state=partial_pcr_primer('i7'))
    blocking=g.input('surface_blocking','Preloaded source i7 primer and BSA provide nonspecific surface blocking; no covalent incorporation at SDS step',K.CONDITION)
    sds=g.input('sds','2uL1%SDS per source well;55C15min',K.CONDITION)
    free=allocated_state.dissociate()
    released=g.product('released','Released adapter-tagged fragments within one physical well',K.MOLECULAR,state=free.scene(),chemical_state=free)
    g.react('dissociate','Remove transposase association only after compartmentalization',S.DEPROTEINIZATION,
            ((aliquot,R.TEMPLATE),(blocking,R.CONDITION),(sds,R.CONDITION)),(released,))
    i5=g.input('pcr_i5','Historical FC-121-1012 i5 primer member; exact polymer unassembled',K.OLIGO,state=partial_pcr_primer('i5'))
    npm=g.input('npm','Nextera PCR Master Mix FC-121-1031;72C3min initial extension;20cycles',K.CATALYST)
    amplified=g.product('amplified','INFERRED — representative A/B PCR library',K.MOLECULAR,state=final_reference(),inferred=True)
    other_products=g.product('other_products','Possible other end/amplification classes; source distribution unassigned',K.MOLECULAR)
    g.react('pcr','INFERRED — add physical-compartment PCR index pair',S.PCR,
            ((released,R.TEMPLATE),(i5,R.STRUCTURAL),(blocker,R.STRUCTURAL),(npm,R.CATALYST)),(amplified,other_products),
            geometry=PCRGeometry.unavailable(free.scene(),reason='Historical FC-121-1012 primer polymers and exact annealing boundaries are unassembled. Figure3 supplies part order, not a complete validated circular/linear PCR product. Same-end class behavior is not assumed.'))
    g.mark_terminal(Terminal(other_products.id,Disposition.UNRESOLVED,'No absolute same-primer amplification prohibition or physical deletion is asserted.'))
    other_libs=g.input('other_libraries','Corresponding indexed PCR products from other physical wells',K.MOLECULAR)
    all_libs=g.product('all_libraries','Pooled physical-well PCR products',K.MOLECULAR,state=final_reference(),inferred=True)
    g.react('pool_pcr','Pool50uL from each source PCR reaction',S.POOL,((amplified,R.TEMPLATE),(other_libs,R.TEMPLATE)),(all_libs,))
    column=g.input('column','Zymo Clean&Concentrator500;two washes;2mL elution',K.SELECTION)
    library=g.product('library','Purified pooled CPT library; conditional assembly only',K.LIBRARY,state=final_reference(),inferred=True)
    waste=g.product('column_waste','Actual column flowthrough and washes',K.WASTE)
    g.react('cleanup','Column-purify pooled libraries and assess insert sizes',S.CLEANUP,((all_libs,R.TEMPLATE),(column,R.SELECTION)),(library,waste))
    g.discard(waste,'Physically removed purification fractions; no quantitative enrichment or genotype purity claim')
    g.sequence(library,SequencingPlan(library=final_reference(),primers=KNOWN_PRIMERS,acquisition='unspecified',limitation='Associated manuscript specifies HiSeq2000 PE51; corrected publisher Figure4 reports51/16/16/51 with bipartite index parts. Full physical-cycle mapping across intervening connectors/custom recipe is unassembled; no ordinary contiguous16-cycle read is fabricated.'))
    g.declare_reviewed_scope(('hmw_dna','adapter_stock','tn5_protein','tn5_activity','loaded','load','tagment_condition','complexes','transpose','edta','stopped','quench','other_transposition_wells','pooled','pool_transposition','allocation','aliquot','other_aliquots','allocate_wells','i7_blocker','surface_blocking','sds','released','dissociate','pcr_i5','npm','amplified','other_products','pcr','other_libraries','all_libraries','pool_pcr','column','library','column_waste','cleanup'),SOURCE)
    g.add_source_item(SourceItem('historical_pcr_and_cycle_map','Historical PCR primers and custom physical-cycle recipe',SOURCE+'; explicit unresolved model/source-version boundary'))
    return g


def oligos():
    rows=[]
    for side,choices in [('A',T5),('B',T7)]:
        rows.extend(oligo('transposon '+side+str(i+1),list(transposon(side,x))) for i,x in enumerate(choices))
    return (*rows,oligo('pMENTS',[Segment('ME complement',nx.ME_RC)],mods='5Phos'),
            oligo('Read1 extended ME A',[Segment('source run primer',EXTENDED_A)]),
            oligo('Read2 extended ME B',[Segment('source run primer',EXTENDED_B)]),
            oligo('Index1 primer',[Segment('source run primer',revcomp(EXTENDED_B))]))


def render_page():
    g=protocol_graph();sc=Scene();sc.strand('reference',list(final_reference()),label='conditional Figure3 part order')
    return '\n'.join((head('CPT-seq — bounded source model'),'<main class="wrap"><h1>CPT-seq</h1>',
        '<p class="research-notes"><a href="01_cpt-seq.html">Research notes</a></p>',
        info('INFERRED — A1/B1 example from corrected Table4, with unassigned physical-well PCR indexes. Protein-held contiguity is not an intact covalent DNA backbone. Historical PCR and custom bipartite-index read mapping remain unresolved; this private partial page has no complete-protocol marker.'),
        info('Publisher correction: the first five Table4 rows were originally missing; the corrected20May2015 workbook is used. Figure3 appears to print connectorA ending TCCACCC, whereas that workbook prints TCCACGC. The model uses the workbook; the graphic-letter difference remains unadjudicated. No source value is silently repaired.'),
        *oligos(),workflow_panel(workflow_from_protocol_graph(g,'library')),
        panel(sc.rows(),caption='INFERRED — four component codes separated by the source adapters/connectors; no contiguous index-cycle map assumed.'),
        sp.section(final_reference(),KNOWN_PRIMERS,acquisition_unspecified=True,
                   intro='Known source run-primer sites on the conditional reference; complete physical cycle acquisition remains unassembled.'),
        '<p>Corrected Figure4:Read1 cycles1–51;Index1 T7 then PCRi7,cycles52–67;Index2 PCRi5 then T5,cycles68–83;Read2 cycles84–134. How the custom recipe handles physical intervening connectors is unresolved.</p></main>'))


def selftest():
    import re
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();source=Source(sources_root()/HERE.parent.name/'ng.3119_41588_2014_BFng3119_MOESM2_ESM.xlsx.txt','corrected publisher Table4 with20May2015 notice')
    if have(c,source,label='CPT corrected source oligos'):
        text=re.sub(r'\s+','',source.path.read_text())
        for side,choices in [('A',T5),('B',T7)]:
            for i,x in enumerate(choices):c('source transferred '+side+str(i+1),transposon(side,x).top() in text,True)
        for name,seq in [('R1',EXTENDED_A),('R2',EXTENDED_B),('I1',revcomp(EXTENDED_B)),('ME','5Phos/'+nx.ME_RC)]:c('source '+name,seq in text,True)
    g=protocol_graph()
    c('SDS receives a physical-well aliquot, not the pre-dilution pool',g.steps['dissociate'].inputs[0].material,'aliquot')
    c('Surface-blocking i7 contributes as a primer only at PCR',g._consumers['i7_blocker'],['pcr'])
    try:ContiguityState(True,'transposed').dissociate()
    except ValueError:rejected=True
    else:rejected=False
    c('Main-branch construction rejects SDS before physical allocation',rejected,True)
    c('Quench preserves association and SDS explicitly removes it',
      (g.materials['stopped'].chemical_state.associated,g.materials['released'].chemical_state.associated),(True,False))
    c('Four partition components retain one group and distinct roles',
      sorted(s.feature.part for s in final_reference() if s.feature),['P5','P7','T5','T7'])
    c('No unvalidated historical primer/cycle recipe becomes a finished graph',bool(g.publication_problems()),True)
    c('Known source run primers land on the conditional source reference',
      sp.verify(final_reference(),KNOWN_PRIMERS,required_roles=()),[])
    c.report()
