"""Source-local bounded partial model; no inherited kit recipe or finished library claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Construct,Scene,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'

def unknown(label):return Construct([Segment(label,'[unassembled]',placeholder=True)])
def cleanup(g,key,x,label,procedure,state=None):
    sel=g.input(key+'_selection',procedure,K.SELECTION)
    kept=g.product(key+'_kept',label,K.MOLECULAR,state=state)
    waste=g.product(key+'_waste','Physically separated '+key+' rejects',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Source-described physical separation; purity and quantitative recovery are not assumed')
    return kept

def aliquot(g,key,x,label):
    chosen=g.product(key+'_chosen',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Untransferred '+key+' material',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(chosen,rest))
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source aliquot remainder fate is not specified in the inspected section; no disposal inferred'))
    return chosen

def limit(g,x):
    g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,LIMIT))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Remaining source library/input chemistry and full strand/read assembly',LIMIT))
    return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="'+NOTE+'">Research notes</a></p>',
        info('Private incomplete source model. '+LIMIT),detail(),
        ''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ENDPOINTS),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Unassembled source and endpoint remain unpublished',g.publishable,False)
    c('Private preview has no completion marker','data-protocol-complete=' in render_page(),False)
    c('Unreviewed geometry does not assert source absence','[adapter overhang not reported]' in render_page(),False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/SOURCE_FILE
    if p.exists():
        source=''.join(p.read_text().split()).upper()
        for seq in TRANSCRIBED:c('Original transcription '+seq,seq.upper() in source,True)
    else:c.skip('Source transcription','External original unavailable')
    source_checks(c,g)
    c.report()
import illumina as il
from adapter_ligation import AdapterLigationGeometry
def unreviewed_ligation_geometry(reason):
    sc=Scene()
    sc.strand('unassembled end',[Segment('source geometry not assembled','[unreviewed molecular geometry]',placeholder=True)],label='Local reading/assembly boundary')
    sc.note('unassembled end',reason+' This is not a claim that the original source fails to report an end.')
    return AdapterLigationGeometry('unavailable','unavailable',(sc,),reason)

TITLE='CRISP-seq — unfragmented UGI aRNA enrichment branch'
NOTE='01_crisp_seq.html'
SOURCE='Published Cell167:1883–1896.e1–e7 DOI10.1016/j.cell.2016.11.039; institutionalCELL2016.pdf STARMethods pp e2–e4'
LIMIT='Source-specific UGI library differs from post-WTA Perturb-seq dial-out. Historical MARS-seq front end, RNA end preparation/ligation enzyme and final PCR/read geometry remain unassembled. Earlier cell7/RMT8/UGI8 roles do not assign the new forward-primer N5 window a UMI function.'
SOURCE_FILE='CELL2016.pdf.txt'
ADAPTER='ATGATCAAGCGACCACCGAG'
RT=revcomp(ADAPTER)
FWD_BODY='CTACACGACGCTCTTCCGATCT'
FWD_BINDING='TCCCCGCGTCGACGGATC'
REV=il.P7+il.TRUSEQ_READ2
TRANSCRIBED=(ADAPTER,RT,FWD_BODY+'N'*5+'XXXX'+FWD_BINDING,REV)
ENDPOINTS=('ugi_pcr','fragmented')

def adapter_reagent():
    sc=Scene();sc.strand('adapter',[Segment('UGI ligation adaptor',ADAPTER)],mod5='phosphate',mod3='C3 nonnucleotide spacer')
    return sc

def adapter_scene():
    sc=adapter_reagent()
    sc.anneal('rt',[Segment('source secondRT primer',RT)],to='adapter',pair=('source secondRT primer','UGI ligation adaptor'))
    return sc

def forward():return Construct([Segment('source partialRead1 handle',FWD_BODY),Segment('random5 window','N'*5,placeholder=True,note='Source random bases; function unassigned, not declared a new UMI'),
    Segment('plate barcode','N'*4,'idx',placeholder=True,feature=feature('crisp.plate','sample_index','unknown')),
    Segment('UGI-enrichment binding arm',FWD_BINDING)])

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    arna=g.input('arna','Source pooled MARS-seq T7-amplified RNA after DNaseI and1.2×SPRI; cell/UMI identity introduced upstream',K.MOLECULAR,primary=True)
    gene=g.product('gene_fraction','Half of source aRNA for MARS-seq transcriptome library',K.MOLECULAR)
    ugi=g.product('ugi_fraction','10–12% of source aRNA for UGI enrichment, UNfragmented',K.MOLECULAR)
    rest=g.product('arna_remainder','Unallocated source aRNA',K.MOLECULAR)
    g.react('aliquot_branches','Take separate aRNA fractions for two physically distinct libraries',S.SPLIT,((arna,R.TEMPLATE),),(gene,ugi,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source proportions leave unassigned remainder; not presumed discarded'))
    condition=g.input('fragmentation','Historical MARS-seq RNA-fragmentation reaction; exact conditions unassembled',K.CONDITION)
    fragmented=g.product('fragmented','Fragmented gene-expression aRNA',K.MOLECULAR)
    g.react('gene_fragmentation','Fragment transcriptome branch only',S.FRAGMENTATION,((gene,R.TEMPLATE),(condition,R.CONDITION)),(fragmented,))
    g.mark_terminal(Terminal(fragmented.id,Disposition.UNRESOLVED,'Remaining MARS-seq ligation/RT/PCR/run geometry unassembled; fragmentation is not applied to UGI branch'))
    adapter=g.input('adapter','Source5primephosphate/3primeC3 UGI ligation adaptor',K.OLIGO,state=adapter_reagent())
    ligase=g.input('ligase','Source MARS-seq ligation reagent identity/conditions unreviewed',K.CATALYST)
    joined=g.product('ligated','UGI aRNA adapter-ligation reaction with free/off-product fractions',K.MOLECULAR)
    g.react('ligation','Ligate source UGI adaptor to unfragmented aRNA',S.LIGATION,((ugi,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(joined,),geometry=unreviewed_ligation_geometry('Source donor5phosphate and3C3 spacer known; original aRNA end chemistry and historical MARS ligation recipe not inspected.'))
    clean=cleanup(g,'ligation_cleanup',joined,'Cleaned UGI ligation material','Same physical ligation cleanup as source MARS-seq; detailed conditions deferred')
    rtprimer=g.input('rtprimer','Printed secondRT UGI primer',K.OLIGO,state=Construct([Segment('secondRT UGI primer',RT)]))
    enzyme=g.input('reverse_transcriptase','Source MARS-seq reverse transcription reagent, exact member unreviewed',K.CATALYST)
    cdna=g.product('ugi_cdna','UGI cDNA plus residual RNA/reaction products',K.MOLECULAR)
    g.react('second_rt','Reverse transcribe UGI material with source-specific primer',S.REVERSE_TRANSCRIPTION,((clean,R.TEMPLATE),(rtprimer,R.STRUCTURAL),(enzyme,R.CATALYST)),(cdna,))
    f=g.input('forward','Printed barcoded forward primer with random5 and plate4 fields',K.OLIGO,state=forward())
    r=g.input('reverse','Printed P7_Rd2 reverse primer',K.OLIGO,state=Construct([Segment('P7_Rd2 primer',REV)]))
    poly=g.input('pcr_mix','Source MARS-seq PCR reagents/program, exact recipe unreviewed',K.CATALYST)
    pcr=g.product('ugi_pcr','Intermediate UGI PCR material after10cycles and plate-barcode introduction',K.MOLECULAR)
    g.react('plate_pcr','Source intermediate10cycle enrichment/pool-barcodePCR',S.PCR,((cdna,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),(pcr,),geometry=PCRGeometry.unavailable(unknown('UGIcDNA'),reason='Primer sequences known, but full original MARS-derived opposite strand/end assembly and source intermediate workup remain unassigned.'))
    return limit(g,pcr)

def detail():return panel(adapter_scene().rows(),caption='Source adaptor and computed complementary secondRT primer; C3 is not an extra cytidine')
def source_checks(c,g):
    c('UGI PCR branch is not fragmented before ligation',g.steps['ligation'].inputs[0].material,'ugi_fraction')
    c('Gene-expression branch is physically separate',g.steps['gene_fragmentation'].inputs[0].material,'gene_fraction')
    c('RT primer is computed full adaptor complement',RT,'CTCGGTGGTCGCTTGATCAT')
    c('Consumed ligation reagent has one strand with source end modifications',tuple((x.mod5,x.mod3) for x in g.materials['adapter'].state.strands.values()),(('phosphate','C3 nonnucleotide spacer'),))
    c('Plate barcode structured independently of source random5 window',forward().get('plate barcode').feature.role,'sample_index')

REVIEWED=('arna', 'gene_fraction', 'ugi_fraction', 'arna_remainder', 'fragmentation', 'fragmented', 'adapter', 'ligase', 'ligated', 'ligation_cleanup_selection', 'ligation_cleanup_kept', 'ligation_cleanup_waste', 'rtprimer', 'reverse_transcriptase', 'ugi_cdna', 'forward', 'reverse', 'pcr_mix', 'ugi_pcr', 'aliquot_branches', 'gene_fragmentation', 'ligation', 'ligation_cleanup', 'second_rt', 'plate_pcr')
