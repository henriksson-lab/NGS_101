from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,K.CATALYST);inputs.append((z,R.CATALYST))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g

TITLE='DNA Typewriter — genomic and single-cell recorder recovery'
SOURCE='Published PMC9352581 Methods200–240 and MOESM4 workbook SupplementaryTables4–5'
LIMIT='Private partial operational models. Recorder manufacture/editing and recovery are distinct. Actual edited allele sequences,full library duplexes, historical CG000205 RevD details,UMI primer identity and PacBio kit geometry remain unassembled. No perfect edit order,editing efficiency or unique lineage-per-cell assumption.'
FWD='GCGTCAGATGTGTATAAGAGACAGGGCTTTATATATCTTGTGGAAAGGACGAAACACC'
REV='TTCAGACGTGTGCTCTTCCGATCTTTGTCTCAAGATCTAGTTACGCCAAGC'
def recorder():
    sc=Scene();sc.strand('recorder',[Segment('Source tape allele','[source/sample-dependent tandem DNA Tape and inserted symbols]',placeholder=True,feature=feature('typewriter.tape','lineage_barcode','unknown'))]);sc.note('recorder','Prime editing installs source symbols+key;actual alleles not inferred from design. Multiple integrated tapes may occur.');return sc

def capture_inventory():
    sc=Scene()
    for key,n,role in (('cell',16,'cell_barcode'),('umi',12,'umi')):
        sc.strand(key,[Segment('source '+key+' component','N'*n,placeholder=True,feature=feature('typewriter.sc_'+key,role,'unknown'))],label='Component inventory only;whole primer/order unassembled')
    return sc

def variant_graph(v):
    if v not in ('bulk','bulk_umi','pacbio','single_cell'):raise ValueError(v)
    g=ProtocolGraph(TITLE+' / '+v,authority=SOURCE);op,clean=operations(g)
    cells=g.input('recorded_cells','Source cells after recorder creation/editing;PE2 and enhanced PE configurations separate upstream experiments',K.BIOLOGICAL,primary=True,state=recorder())
    if v=='single_cell':
        x=clean('cell_prepare',cells,'Dissociate,pellet200RCF5min,resuspend0.04%BSA/PBS')
        x=op('partition_capture',x,'Historical CG000205RevD single-cell/bead partitions and capture;do not borrow gRNA sequence recipe','Historical Chromium components')
        x=op('rt',x,'Historical RT labels transcripts with16ntCBC and12ntUMI;whole primer unassembled','Historical RT chemistry',state=capture_inventory())
        x=op('amplify',x,'Historical cDNA recovery/amplification throughStep2.3;intermediate operational details unreviewed','Historical CG000205RevD reagents')
        a=g.product('fraction_a','Step2.3A largercDNA eluate;containsTape AND gene expression species',K.MOLECULAR)
        b=g.product('fraction_b','Step2.3B smallerfeature eluate;source also detectsTape',K.MOLECULAR)
        r=g.input('selection','Historical size selection exact geometry unreviewed',K.SELECTION)
        g.react('separate','Recover BOTHStep2.3A and2.3B eluates',S.CLEANUP,((x,R.TEMPLATE),(r,R.SELECTION)),(a,b))
        target=[]
        for key,y in (('a',a),('b',b)):
            t=g.product(key+'_target','Tape-directed preparation from source fraction'+key,K.MOLECULAR)
            gene=g.product(key+'_other','Other recovered cDNA incl expression-capable material;not presumed discarded',K.MOLECULAR)
            g.react(key+'_allocation','Source target components recovered;exact allocation/workup unassembled',S.SPLIT,((y,R.TEMPLATE),),(t,gene))
            from protocol_graph import Terminal,Disposition
            g.mark_terminal(Terminal(gene.id,Disposition.UNRESOLVED,'Source also obtains expression library;allocation between these fractions and that library not fully specified in inspected main.'))
            target.append(t)
        pooled=g.product('combined_tapes','CombinedTape components from BOTH recovered eluates',K.MOLECULAR)
        g.react('combine_tapes','Combine recovered TargetBC–Tape components',S.POOL,tuple((y,R.TEMPLATE) for y in target),(pooled,))
        x=op('target_library',pooled,'Paired-end target-library continuation;local whole polymer/end geometry unassembled','Historical source target-library components')
        end(g,x,'Source paired-end recovery of16ntCBC,12ntUMI,TargetBC andTape;gene-expression library also obtained. CBC/UMI/TargetBC/scar functions distinct.No perfect size separation inferred.')
    else:
        x=op('lysis',cells,'Cell lysis/protease preparation cited toprime-editing predecessor;detailed linked procedure unreviewed','Source cited lysis/protease')
        from protocol_graph import Terminal,Disposition
        used=g.product('lysate_aliquot','1.5uL source cell lysate allocated to library reaction',K.MOLECULAR)
        rest=g.product('lysate_remainder','Unallocated lysate',K.MOLECULAR)
        g.react('take_aliquot','Take1.5uL lysate aliquot',S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Remainder not assigned a library or disposal.'));x=used
        if v=='bulk_umi':
            u=g.input('umi_primer','Single source forward primer;10ntUMI;full oligo selection unassembled',K.OLIGO,state=Construct([Segment('UMI component','N'*10,placeholder=True,feature=feature('typewriter.bulk_umi','umi','random'))]))
            out=g.product('linear_product','UMI-extended strands plus original DNA/free reagents',K.MOLECULAR)
            enz=g.input('linear_enzyme','KAPA Robust;source two repeatedanneal/extension blocks separated bydenaturation',K.CATALYST)
            g.react('linear_extension','Single-primer UMI extension;do not treat as exponential PCR',S.EXTENSION,((x,R.TEMPLATE),(u,R.STRUCTURAL),(enz,R.CATALYST)),(out,));x=clean('linear_cleanup',out,'1.5xAMPure recovery')
        f=g.input('forward','Published TAPE_amp_fwd inventory;source variants use other specifically named primers',K.OLIGO,state=Construct([Segment('TAPE_amp_fwd',FWD)]) if v=='bulk' else None)
        r=g.input('reverse','Published TAPE_amp_rev inventory',K.OLIGO,state=Construct([Segment('TAPE_amp_rev',REV)]) if v=='bulk' else None)
        y=g.product('target_pcr_product','TargetedPCR products plus retained reaction reagents',K.MOLECULAR)
        g.react('target_pcr','KAPA Robust targetPCR25–28cycles;bulkUMI uses upstreamUMI handle;PacBio uses source indexedlong primers,not displayed coreprimer pair',S.GENERIC,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL)),(y,))
        if v=='pacbio':
            x=clean('long_cleanup',y,'0.6xAMPure;TapeStationQC')
            x=op('pacbio_library',x,'Source Sequel amplicon submission;detailed SMRTbell preparation unreviewed','Source facility library components')
            end(g,x,'Source PacificBiosciencesSequel CCS. Sample-indexed long-amplicon primer identities and kit polymers unassembled;do not invent a circle geometry.')
        else:
            x=clean('target_cleanup',y,'6%TBEgelassessment then1xAMPure recovery')
            x=op('index_pcr',x,'5–10cycle dual-sample-index/flowcell-adapterPCR','KAPA Robust plus source indexPCR oligos unassembled')
            x=clean('index_cleanup',x,'AMPure ratio unspecified;TapeStationQC')
            x=op('denature',x,'Library denaturation beforeIllumina sequencing','Denaturation conditions unassembled')
            end(g,x,'Source NextSeq/MiSeq;Fig1single-end,Fig2/3paired-end. Sampleindexes distinctfrom10ntbulkUMI andinheritedTape;fullstrand/readgeometryunassembled.')
    return incomplete(g)
def protocol_graph():return variant_graph('bulk')
def render_page():return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_typewriter.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),'denature_product'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();gs=[variant_graph(v) for v in ('bulk','bulk_umi','pacbio','single_cell')]
    c('Allroutes private',any(g.publishable for g in gs),False)
    c('UMI linear extension precedes cleanup',gs[1].steps['linear_cleanup'].inputs[0].material,'linear_product')
    c('PacBio doesnot receiveIlluminaindexPCR','index_pcr' in gs[2].steps,False)
    c('Bothsizefractions contribute',tuple(i.material for i in gs[3].steps['combine_tapes'].inputs),('a_target','b_target'))
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41586-022-04922-8_supp_41586_2022_4922_MOESM4_ESM.xlsx.txt'
    if p.exists():
        t=p.read_text().upper();c('Publishedforwardprimer',FWD in t,True);c('Publishedreverseprimer',REV in t,True)
    else:c.skip('Primertranscription','External workbook absent')
    c.report()
