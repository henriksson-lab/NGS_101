"""LINNAEUS source-qualified conditional Methods routes and published SI primer inventory."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct,Scene,Segment,feature,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SequencingPlan
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='LINNAEUS — transcriptome and RFP-scar recovery'
SOURCE='Conditional accepted-manuscript PMC5942543 Methods82–90,plus separately identified publishedSI MOESM11 Table1 originalPDFpage41;publishedmainHTML is subscriptionpreview,not finalMethods verification'
LIMIT='Private incomplete source-version-qualified model. Accepted-manuscript operational order and published SI primer inventory are identified separately;whole final-version equivalence is unestablished. Historical Chromium guide,full recorder/primer/product strands and exact read-primer geometry remain unassembled.'
TABLE='nbt.4124_41587_2018_BFnbt4124_MOESM11_ESM.pdf.txt'
RFP1='CGGCACGCTGATCTACAAGG'
RFP2=il.TRUSEQ_READ2+'GAGTTCAAGACCATCTACATGGCC'
REVERSE=il.TRUSEQ_READ1[-22:]
BULK_HANDLE='CCTTGGCACCCGAGAATTCCA'
BULK_REVERSE='GTTCAGAGTTCTACAGTCCGACGATC'+'GTCCTCGTTGTGGGAGGTG'
LINEAGE=feature('linnaeus.rfp_scar','lineage_barcode','unknown',note='Inherited RFP scars;multiple independently integrated reporter loci per cell,not a unique current-cell barcode.')
CELL=feature('linnaeus.chromium_cell','cell_barcode','unknown')
UMI=feature('linnaeus.chromium_umi','umi','unknown')

@dataclass(frozen=True)
class RecorderContext:
    independent_integrations:str='16–32 source ZebrabowM RFP integrations'
    guide_target:str='One source RFP target per reporter'
    def rows(self):
        sc=Scene();sc.strand('scar',[Segment('source RFP scar','[source/sample-dependent RFP indel allele]',placeholder=True,feature=LINEAGE)])
        sc.note('scar',self.independent_integrations+';creation/editing upstream of recovery,not an indel/cell-uniqueness prediction.');return sc.rows()

def bulk_forward(which):
    barcode={'a':'CTGGTCGG','b':'AGGCTCTA'}[which]
    return Construct([Segment('source bulkPCR handle',BULK_HANDLE),Segment('source inline sample barcode',barcode,'idx',feature=feature('linnaeus.bulk_sample','sample_index','fixed')),
        Segment('RFP forward arm',RFP1)],'Published bulk forward'+('1' if which=='a' else '2'))

def bead_inventory():
    sc=Scene()
    for key,f in (('cell',CELL),('umi',UMI)):
        sc.strand(key,[Segment('source '+key+' component','[width/placement unassembled]',placeholder=True,feature=f)],label='Known component class,not a reconstructed whole bead primer')
    sc.note('cell','Source polyT capture and adapter regions also stated;historical sequence/order/lengths unassembled.');return sc

def variant_graph(variant):
    if variant not in ('single_cell','bulk_dna','bulk_rna'):raise ValueError('Choose a source-defined recovery branch')
    g=ProtocolGraph(TITLE+' / '+variant,authority=SOURCE)
    def inp(key,label,kind=K.CATALYST,state=None):return g.input(key,label,kind,state=state)
    def op(key,x,label,kind,reagents=(),geometry=None):
        y=g.product(key+'_product',label+' reaction preparation;mixture not presumed pure',K.MOLECULAR);kw={} if geometry is None else {'geometry':geometry}
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(reagents),(y,),**kw);return y
    def clean(key,x,label):
        sel=inp(key+'_selection',label,K.SELECTION);y=g.product(key+'_kept','Source retained preparation',K.MOLECULAR);w=g.product(key+'_waste','Physically removed workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(y,w));g.discard(w,'Actual source workup,not computational scar filtering or complete purity assertion');return y
    def pcr(key,x,f,r,label):
        enzyme=inp(key+'_polymerase','Source PCR reagent/program unassembled;not inherited from later CARLIN/scGESTALT')
        return op(key,x,label,S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),PCRGeometry.unavailable(Scene(),reason='Named source primers/handles known where displayed;full RFP alleles,kit-derived opposite strand and complete product ends unassembled.'))
    if variant=='single_cell':
        cells=g.input('cells','Source scarred ZebrabowM single-cell suspension',K.BIOLOGICAL,primary=True,state=RecorderContext())
        bead=inp('bead','Historical ChromiumPN120233 bead capture components;exact guide/version unreviewed',K.OLIGO,bead_inventory())
        droplets=op('encapsulate',cells,'Source single-cell/barcoded-bead droplet encapsulation',S.GENERIC,((bead,R.STRUCTURAL),))
        lysis=inp('lysis_reagents','Historical source Chromium lysis chemistry unassembled',K.CONDITION)
        lysed=op('cell_lysis',droplets,'Source droplet cell lysis;soluble/insoluble species remain co-present',S.GENERIC,((lysis,R.CONDITION),))
        enzyme=inp('rt_enzyme','Historical source Chromium RT reagents unassembled')
        cdna=op('capture_rt',lysed,'PolyT reverse transcription introduces cell identifier/UMI',S.REVERSE_TRANSCRIPTION,((bead,R.STRUCTURAL),(enzyme,R.CATALYST)))
        pooled=op('pool_cdna',cdna,'Source pooling of droplet-derived material;exact emulsion workup not reviewed',S.GENERIC)
        amp=inp('wta_reagents','Source first-amplification reagents/primers unassembled')
        amplified=op('initial_amplification',pooled,'Source first cDNA amplification',S.GENERIC,((amp,R.CATALYST),))
        gene=g.product('transcriptome_half','Half amplified whole-cDNA mixture for fragmented transcriptome library',K.MOLECULAR)
        scar=g.product('scar_half','Half amplified whole-cDNA mixture retained UNfragmented for targeted RFP',K.MOLECULAR)
        g.react('library_split','Source half-and-half material allocation;not modality purification',S.SPLIT,((amplified,R.TEMPLATE),),(gene,scar))
        frag=inp('gene_fragmentation','Historical source transcriptome fragmentation reagents unassembled')
        gene=op('gene_fragment',gene,'Fragment transcriptome branch only',S.FRAGMENTATION,((frag,R.CATALYST),))
        kit=inp('gene_library_kit','Historical Chromium remaining library components unreviewed',K.CONDITION)
        lib=g.product('gene_library_product','Transcriptome library from historical source kit;polymer/end geometry unassembled',K.MOLECULAR)
        g.react('gene_library','Unassembled historical Chromium downstream library operations',S.MAGIC,((gene,R.TEMPLATE),(kit,R.CONDITION)),(lib,),reason='Defining source refers to manufacturer protocol;exact historical original/version not read.No current kit or source-selected shared recipe substituted.')
        g.sequence(lib,SequencingPlan(acquisition='unspecified',limitation='Source reports NextSeq5002×75 andHiSeq25002×100;mapping those runs/primer lengths to each sample and full library geometry unassembled.'))
        f1=inp('rfp_f1','PublishedSI RFP forward1 inventory',K.OLIGO,Construct([Segment('RFP forward1',RFP1)]))
        f2=inp('rfp_f2','PublishedSI nested RFP forward2 inventory',K.OLIGO,Construct([Segment('nested RFP forward2',RFP2)]))
        reverse=inp('scar_reverse','PublishedSI source reverse handle',K.OLIGO,Construct([Segment('reverse',REVERSE)]))
        scar=pcr('scar_pcr1',scar,f1,reverse,'Source first RFP-targeted PCR;unfragmented substrate')
        scar=pcr('scar_pcr2',scar,f2,reverse,'Source second nested RFP PCR;no invented intermediate cleanup')
        g.sequence(scar,SequencingPlan(acquisition='unspecified',limitation='Source scar reads carry same cell identifier/UMI as transcript reads;NextSeq5002×75 orHiSeq25002×100 reported,full product/primer geometry unassembled. Computational truncation to75nt is not physical DNA trimming.'))
    else:
        rv=inp('bulk_reverse','PublishedSI bulk scar reverse primer',K.OLIGO,Construct([Segment('bulk reverse',BULK_REVERSE)]))
        products=[]
        for key in ('a','b'):
            # Two symbolic samples expose actual physical pooling;the named source
            # primer members do not assign any particular experimental animal.
            sample=g.input(key+'_sample','Independent source sample '+key+';source allele identities unassigned',K.BIOLOGICAL if variant=='bulk_dna' else K.MOLECULAR,primary=True,state=RecorderContext())
            if variant=='bulk_dna':
                alkali=inp(key+'_alkali','50mMNaOH95C20min',K.CONDITION);x=op(key+'_lysis',sample,'Source heated alkaline DNA preparation;debris/reaction species retained',S.EXTRACTION,((alkali,R.CONDITION),))
                tris=inp(key+'_tris','One-tenthvolume1MTris-HClpH8.4',K.CONDITION);x=op(key+'_neutralize',x,'Neutralize without physical purification',S.QUENCHING,((tris,R.CONDITION),))
                used=g.product(key+'_used','20uL source neutralized preparation for PCR',K.MOLECULAR);unused=g.product(key+'_unused','Untransferred DNA preparation',K.MOLECULAR)
                g.react(key+'_aliquot','Take20uL PCR aliquot',S.SPLIT,((x,R.TEMPLATE),),(used,unused));g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Remainder not presumed disposed.'));x=used
            else:
                trizol=inp(key+'_trizol','Source TRIzol isolation;vendor detailed phase chemistry unreviewed',K.CONDITION);x=op(key+'_rna_isolation',sample,'Source TRIzol RNA preparation',S.EXTRACTION,((trizol,R.CONDITION),))
                x=clean(key+'_rna_recovery',x,'Isopropanol RNA precipitation;two75%ethanol pellet washes;air dry')
                dt=inp(key+'_polyt','Source polyT RT primer;length/full sequence unassembled',K.OLIGO)
                rt=inp(key+'_ssii','SuperScriptII/RNaseOUT/dNTPs/sourcebuffer42C2h')
                x=op(key+'_rt',x,'PolyT reverse transcription for bulk RNA scar detection',S.REVERSE_TRANSCRIPTION,((dt,R.STRUCTURAL),(rt,R.CATALYST)))
            f=inp(key+'_bulk_forward','PublishedSI bulk primer member example;experimental sample assignment unassigned',K.OLIGO,bulk_forward(key))
            products.append(pcr(key+'_bulk_pcr',x,f,rv,'Source barcoded RFP scar amplification'))
        pool=g.product('bulk_pool','Pooled sample-barcoded firstPCR products',K.MOLECULAR)
        g.react('bulk_pooling','Pool firstPCR products before adapterPCR',S.POOL,tuple((x,R.TEMPLATE) for x in products),(pool,))
        pool=clean('bulk_ampure',pool,'AMPure magnetic-bead cleanup')
        f=inp('adapter_pcr_forward','Source Illumina-adapter PCR forward primer;full identity not located in this source pass',K.OLIGO)
        r=inp('adapter_pcr_reverse','Source Illumina-adapter PCR reverse primer;full identity not located in this source pass',K.OLIGO)
        out=pcr('bulk_adapter_pcr',pool,f,r,'SecondPCR adds Illumina sequencing adapters')
        g.sequence(out,SequencingPlan(acquisition='unspecified',limitation='Source bulk scar assay usesNGS;exact final purification/read/index/primer geometry unassembled.No single-cell UMI or cell barcode inferred for this bulk route.'))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Historical kit/final-version/whole strand and read geometry',LIMIT));return g

def protocol_graph():return variant_graph('single_cell')

def render_page():
    ends={'single_cell':('scar_pcr2_product','gene_library_product'),'bulk_dna':('bulk_adapter_pcr_product',),'bulk_rna':('bulk_adapter_pcr_product',)}
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_linnaeus.html">Research notes</a></p>',
        ''.join('<h2>'+v+'</h2>'+''.join(workflow_panel(workflow_from_protocol_graph(variant_graph(v),e)) for e in ends[v]) for v in ends),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();bulk=variant_graph('bulk_dna');rna=variant_graph('bulk_rna')
    c('Allbranches private',any(x.publishable for x in (g,bulk,rna)),False);c('No finished marker','data-protocol-complete=' in render_page(),False)
    c('Scar branch is not the fragmented gene branch',g.steps['scar_pcr1'].inputs[0].material,'scar_half')
    c('Nested scar PCRs carry previous reaction',g.steps['scar_pcr2'].inputs[0].material,'scar_pcr1_product')
    c('Bulk pooling precedes cleanup',bulk.steps['bulk_ampure'].inputs[0].material,'bulk_pool')
    c('RNA source uses polyT SSII RT before RFP PCR',rna.steps['a_bulk_pcr'].inputs[0].material,'a_rt_product')
    p=sources_root()/Path(__file__).resolve().parents[1].name/TABLE
    if p.exists():
        source=''.join(p.read_text().split())
        for name,seq in (('RFP1',RFP1),('RFP2',RFP2),('reverse',REVERSE),('bulkReverse',BULK_REVERSE),('bulkF1',bulk_forward('a').top()),('bulkF2',bulk_forward('b').top())):c('PublishedTable1transcription '+name,seq in source,True)
    else:c.skip('PublishedTable1primer transcriptions','External original unavailable')
    c.report()
