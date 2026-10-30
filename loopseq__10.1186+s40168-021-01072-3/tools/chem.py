"""Historical LoopSeq SSC v1.7 operational route; proprietary chemistry unassembled."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,workflow_from_protocol_graph,workflow_panel,panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from adapter_ligation import AdapterLigationGeometry
from page import head,info
PAGE_FILE='chem.html'
TITLE='LoopSeq — historical SSC v1.7 operational library model'
SOURCE='Published DOI10.1186/s40168-021-01072-3 PMC8179091 Methods; directly linked LoopSeq16S/16S18S MicrobiomeSSC24Plex manual printedVersion1.7(July2019),SectionsI–IV'
MANUAL='gh_LoopSeq_16S_16S18S_Microbiome_SSC_Manual_v1p7.pdf.txt'
LIMIT='Private partial historical v1.7 model. Kit polymers, barcode widths, mix/enzyme identities, full strand/end geometry and final library/read layout remain unassembled. Barcode distribution is the source-declared operation, not a computed transposase, circularization or nicking mechanism. The current vendor protocol is not substituted.'

def unknown_scene():
    sc=Scene();sc.strand('unknown',[Segment('source polymers','[unassembled kit and molecular geometry]',placeholder=True)])
    sc.note('unknown','Unassembled scope, not a source-disclosure absence claim.');return sc

def identity_components():
    sc=Scene()
    sc.strand('molecule',[Segment('parent molecule UMI','[UMI width unassembled]',placeholder=True,feature=feature('loopseq.parent','umi','random'))])
    sc.strand('sample',[Segment('sample identity','[sample tag width unassembled]',placeholder=True,feature=feature('loopseq.sample','sample_index','unknown'))])
    sc.note('molecule','Separate feature classes established in published Methods; not a sequence or junction reconstruction.')
    return sc

def cleanup(g,key,x,label):
    sel=g.input(key+'_selection',label,K.SELECTION)
    out=g.product(key+'_kept','Eluted DNA after '+label,K.MOLECULAR)
    waste=g.product(key+'_waste','Removed unbound/wash fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(out,waste));g.discard(waste,'Discard unbound supernatant/washes at bead separation; retained DNA eluted into fresh tube')
    return out

def transform(g,key,x,label,reagents,kind=S.GENERIC):
    catalyst=g.input(key+'_reagents',reagents,K.CATALYST)
    y=g.product(key+'_product',label+' reaction including residual reagents/products',K.MOLECULAR)
    g.react(key,label,kind,((x,R.TEMPLATE),(catalyst,R.CATALYST)),(y,));return y

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    products=[]
    for key in ('a','b'):
        x=g.input(key+'_dna','Independent genomic DNA sample '+key+'; symbolic class from source24well workflow',K.MOLECULAR,primary=True)
        mix=g.input(key+'_barcoding','Source predispensed Barcoding Adapters; sequences and mix constituents unassembled',K.CONDITION)
        y=g.product(key+'_assigned','Sample/molecule-barcoded target material in '+key+'; other reaction species retained',K.MOLECULAR)
        # The kit's source-named PCR is intentionally a coarse operation here:
        # neither two primers nor a single symmetric primer species is assigned.
        g.react(key+'_assignment','Section1.2 Barcode Assignment PCR,2cycles; internal kit reactants unassembled',S.GENERIC,((x,R.TEMPLATE),(mix,R.CONDITION)),(y,))
        z=transform(g,key+'_inactivation',y,'Section1.3 PCR inactivation','Source Inactivation Enzyme SSC; chemical mechanism unassembled')
        clean=cleanup(g,key+'_cleanup',z,'Section1.4 SPRIselect cleanup')
        q=g.product(key+'_qc','Separate aliquot for destructive calibration qPCR',K.QC)
        rest=g.product(key+'_remaining','Untransferred cleaned sample '+key,K.MOLECULAR)
        g.react(key+'_qc_aliquot','Take calibration qPCR aliquot',S.SPLIT,((clean,R.TEMPLATE),),(q,rest))
        g.mark_terminal(Terminal(q.id,Disposition.QC,'Dedicated analytical aliquot; calibration amplicons do not enter library pool.'))
        chosen=g.product(key+'_calibrated','Source-calibrated input aliquot for barcode amplification',K.MOLECULAR)
        unused=g.product(key+'_unused','Sample not transferred after calibration',K.MOLECULAR)
        g.react(key+'_normalization','Dilute and choose molecule-count-based aliquot',S.SPLIT,((rest,R.TEMPLATE),),(chosen,unused))
        g.mark_terminal(Terminal(unused.id,Disposition.UNRESOLVED,'Untransferred sample fate not assigned; no disposal inferred.'))
        amp=g.input(key+'_ampmix','Source Amplification Mix SSC; exact oligos/polymerase unassembled',K.CATALYST)
        amplified=g.product(key+'_amplified','Barcoded sample '+key+' after22cycle PCR',K.MOLECULAR)
        g.react(key+'_amplification','Section3.1 Barcode Amplification PCR,22cycles; internal kit reactants unassembled',S.GENERIC,((chosen,R.TEMPLATE),(amp,R.CATALYST)),(amplified,))
        aliquot=g.product(key+'_pool_aliquot','5uL from amplified sample '+key,K.MOLECULAR)
        remainder=g.product(key+'_amp_remainder','Unpooled amplified sample '+key,K.MOLECULAR)
        g.react(key+'_take_pool_aliquot','Take5uL for physical pool',S.SPLIT,((amplified,R.TEMPLATE),),(aliquot,remainder))
        g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Amplification remainder not presumed discarded.'))
        products.append(aliquot)
    pool=g.product('pool','Physical pool of distinct sample-tagged amplification aliquots',K.MOLECULAR)
    g.react('sample_pool','Section3.2:combine5uL from each source24well sample;two symbolic classes shown',S.POOL,tuple((x,R.TEMPLATE) for x in products),(pool,))
    clean=cleanup(g,'pooled_cleanup',pool,'Section3.2 SPRIselect cleanup')
    distributed=transform(g,'distribution',clean,'Section3.3 intramolecular barcode distribution20C15min then75C5min','Distribution Mix and Distribution Enzyme; exact chemistry unassembled')
    activated=transform(g,'activation',distributed,'Section3.4 activation20C2h','Activation Mix and Activation Enzyme; molecular mechanism unassembled')
    neutralized=transform(g,'neutralization',activated,'Section3.4 neutralization37C15min','Neutralization Enzyme; molecular mechanism unassembled')
    out=cleanup(g,'postactivation_cleanup',neutralized,'Section3.5 SPRIselect cleanup')
    fragmented=transform(g,'fragment_endprep',out,'Section4.1 combined fragmentation,end repair,A-tailing32C5min65C30min','Fragmentation Mix and Fragmentation Enzyme; separate polymer end states/enzymatic identities unassembled',S.FRAGMENTATION)
    ligase=g.input('ligation_enzyme','Source Ligation Enzyme10uL',K.CATALYST)
    mix=g.input('ligation_mix','Source Ligation Mix40uL; ligating adapter polymers and other constituents unassembled',K.ADAPTER,state=unknown_scene())
    ligated=g.product('ligated','Section4.2 ligation reaction including residual fragments/reagents',K.MOLECULAR)
    g.react('ligation','Section4.2 Adapter Ligation20C15min',S.LIGATION,((fragmented,R.TEMPLATE),(mix,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),
        geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Historical kit adapter polymer and complete ligation-ready end geometry not assembled;source A-tailing is not evidence of an exact adapter T overhang.'))
    clean=cleanup(g,'postligation_cleanup',ligated,'Section4.3 SPRIselect cleanup')
    indexmix=g.input('index_mix','Source Index Master Mix25uL and5uL of one Index Primer P1,P2,P3 orP4; internal priming species unassembled',K.CONDITION)
    indexed=g.product('indexed','Section4.4 indexPCR reaction after12cycles',K.MOLECULAR)
    g.react('index_pcr','Section4.4 Library Index PCR,12cycles; kit molecular priming geometry unassembled',S.GENERIC,((clean,R.TEMPLATE),(indexmix,R.CONDITION)),(indexed,))
    final=cleanup(g,'postindex_cleanup',indexed,'Section4.5 SPRIselect cleanup')
    g.mark_terminal(Terminal(final.id,Disposition.UNRESOLVED,LIMIT+' Source library operation order modeled throughSection4.5;exact sequenced molecules and run geometry unassembled.'))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Remaining original/kit chemistry and complete sequenced library',LIMIT))
    return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',info(LIMIT),'<p><a href="01_loopseq.html">Research notes</a></p>',
        workflow_panel(workflow_from_protocol_graph(g,'postindex_cleanup_kept')),panel(identity_components().rows(),caption='Source-defined identifier meanings only; no inferred sequence joins'),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Private incomplete model',g.publishable,False)
    c('No false completion marker','data-protocol-complete=' in render_page(),False)
    c('Sample pool contains only post-amplification aliquots',[i.material for i in g.steps['sample_pool'].inputs],['a_pool_aliquot','b_pool_aliquot'])
    c('Distribution starts after physical pool cleanup',g.steps['distribution'].inputs[0].material,'pooled_cleanup_kept')
    c('Neutralization follows activation',g.steps['neutralization'].inputs[0].material,'activation_product')
    c('Ligation follows combined fragmentation/end preparation',g.steps['ligation'].inputs[0].material,'fragment_endprep_product')
    c('IndexPCR follows physical ligation cleanup',g.steps['index_pcr'].inputs[0].material,'postligation_cleanup_kept')
    c('Unassembled geometry does not assert undisclosed ends','[adapter overhang not reported]' in render_page(),False)
    c('Calibration aliquot stays outside pool',g.terminals['a_qc'].disposition,Disposition.QC)
    p=sources_root()/Path(__file__).resolve().parents[1].name/MANUAL
    if p.exists():
        text=' '.join(p.read_text().split())
        for part in ('Version 1.7','July 2019','3.3. Barcode Distribution','3.4. Barcode Activation and Neutralization'):
            c('Historical manual identity/operation '+part,part in text,True)
    else:c.skip('Historical manual identity','External original unavailable')
    c.report()
