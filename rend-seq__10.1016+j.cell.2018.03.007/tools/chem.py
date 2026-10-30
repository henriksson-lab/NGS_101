"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=()):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION, S.DEPHOSPHORYLATION, S.PHOSPHORYLATION}
    stock = g.input(key + '_reagent', reagent, K.CATALYST if catalyst else K.CONDITION)
    inputs = [(x, R.TEMPLATE), (stock, R.CATALYST if catalyst else R.CONDITION)]
    inputs.extend((a, R.STRUCTURAL) for a in structural)
    kwargs = {}
    if kind == S.MAGIC:
        kwargs['reason'] = 'Source operation established; exact molecular transition/kit internals unassembled. No accepted absence claim.'
    if kind == S.PCR:
        kwargs['geometry'] = PCRGeometry.unavailable(Scene(), reason='Exact template/primer endpoint geometry unassembled; selected oligos do not establish whole amplicon.')
        if len(structural) < 2:
            for arm in ('f', 'r'):
                inputs.append((g.input(key + '_' + arm, 'Source PCR primer; exact identity/landing unassembled', K.OLIGO, state=scene(key + '_' + arm, '[source primer]')), R.STRUCTURAL))
    if kind == S.REVERSE_TRANSCRIPTION and not structural:
        inputs.append((g.input(key + '_primer', 'Source RT primer; geometry unassembled', K.OLIGO, state=scene(key + '_primer', '[source RT primer]')), R.STRUCTURAL))
    if kind == S.LIGATION:
        kwargs['geometry'] = AdapterLigationGeometry('unavailable', 'unavailable', (scene(key + '_ends', '[source ligation; complete substrate/end geometry unassembled]'),), 'This incomplete model does not claim publication absence.')
        if not structural:
            inputs.append((g.input(key + '_adapter', 'Source ligating adapter; full geometry unassembled', K.ADAPTER, state=scene(key + '_adapter', '[source adapter]')), R.STRUCTURAL))
    y = g.product(key + '_product', product or label + '; physical companions retained until separation', K.MOLECULAR)
    g.react(key, label, kind, inputs, (y,), **kwargs); return y

def clean(g, key, x, label):
    condition = g.input(key + '_condition', label, K.CONDITION)
    kept = g.product(key + '_kept', 'Source physically recovered fraction; carryover/yield unquantified', K.MOLECULAR)
    waste = g.product(key + '_waste', 'Actual separated cleanup fractions; no perfect classification assumed', K.WASTE)
    g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), (condition, R.CONDITION)), (kept, waste))
    g.discard(waste, 'Actual physical separation rather than computational removal.'); return kept

def qc(g, key, x, label):
    kept = g.product(key + '_kept', 'Source retained preparation', K.MOLECULAR)
    assay = g.product(key + '_aliquot', 'Separate analytical aliquot; no return to library assumed', K.QC)
    g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (kept, assay))
    g.mark_terminal(Terminal(assay.id, Disposition.QC, 'Source analytical aliquot.')); return kept

def aliquot(g,key,x,label):
    take=g.product(key+'_used',label,K.MOLECULAR);rest=g.product(key+'_rest','Physical unused remainder; no downstream assignment',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(take,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Retained outside selected path; not computationally discarded.'));return take

def finish(g,x,readout):
    lib=g.product('library','Source physical sequencing preparation; complete molecular endpoints unassembled',K.LIBRARY)
    g.react('interface','Sequencing interface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation=readout+'; exact primer landing/cycle geometry unassembled.'))
    g.add_source_item(SourceItem('pending','Complete oligos, endpoints and linked chemistry originals remain unreviewed','Positive selected preparation only; not an absence audit.'));return g

from circular import CircularizationGeometry
SOURCE='Published final author-hosted Cell173(2018) source_pdf_cb028449b9993071a147.pdf Methods1549–1568,1966–1989; AMPMC5978003 selected sections consistent; S4 unreviewed.'
def protocol_graph_for(exonuclease=False):
    g=ProtocolGraph('Rend-seq'+(' Terminator control'if exonuclease else' baseline'),authority=SOURCE,representation='published_instance')
    x=g.input('rna','Source bacterial RNA:5ml culture+5mlcoldmethanol→pellet/freezing→RNAeasy/on-columnDNase; reaction scoped after extraction',K.MOLECULAR,primary=True)
    if exonuclease:
        x=op(g,'terminator',x,'7.5ugtotalRNA17ulwater+2ul10xbufferA+1ulTerminator;30C60min',S.DIGESTION,reagent='Terminator5primeP-dependentexonuclease',product='Exonuclease-treated RNA plus digestedcomponents; no perfect elimination assumed')
        x=clean(g,'exo_clean',x,'Source isopropanol precipitation; no rRNAremoval for these samples')
    else:
        x=op(g,'rrna',x,'MICROBExpressrRNAdepletion BEFORE fragmentation',reagent='SourceMICROBExpresskit; actual separation internals unassembled',product='Source physically rRNA-depletedRNA with residualRNA/carryover; not computational rRNAremoval')
        x=clean(g,'rrna_clean',x,'Isopropanol precipitation;40ul10mMTrispH7')
    x=op(g,'preheat',x,'RNA95C2min thenice>=1min',S.THERMAL_TREATMENT,reagent='Source heating/cooling',product='Preheated/cooled physical RNA mixture; no cleavage assigned beforefragmentationbuffer')
    x=op(g,'fragment',x,'4.4ul10xRNAfragmentationbufferonice→95C25sec→5ulstopbuffer→ice',reagent='Sourcezinc-mediated sparsefragmentation;25secbaseline only',product='Mixtureof original-end and internallycleaved fragments; no computational extraction of nativeendfragments')
    x=clean(g,'precip',x,'Source isopropanol precipitation')
    x=clean(g,'short',x,'15percentTBE-ureaPAGE;physicallyexcise15–45ntRNA')
    x=op(g,'pnk',x,'T4PNK20U;37C1h;source3primedephosphorylation',S.DEPHOSPHORYLATION,reagent='SourceT4PNK; no ATP/5primephosphorylation inferred')
    adapter=g.input('linker1','SourceTableS4linker1:5primeadenylated,3primeblockedDNA;exactsequenceunreviewed',K.ADAPTER,state=scene('linker1','[source linker1]'))
    x=op(g,'ligate',x,'3pmolRNA+100pmol linker1;t4RNA ligase2truncatedK277Q;37C2.5h;25percentPEG8000',S.LIGATION,reagent='Source truncatedT4RNA ligase2K277Q',structural=(adapter,),product='3prime-adaptedRNAplusunligatedRNA/freeadapter/reagents; exactjunctionunassembled')
    x=clean(g,'lig_clean',x,'10percentTBE-ureaPAGEsizeexcision')
    primer=g.input('ocj485','SourceS4RTprimerocj485;exactsequenceunreviewed',K.OLIGO,state=scene('ocj485','[source RT primer]'))
    x=op(g,'rt',x,'SuperScriptIII50C45min withocj485',S.REVERSE_TRANSCRIPTION,reagent='SourceSuperScriptIII',structural=(primer,))
    x=clean(g,'rt_clean',x,'10percentTBE-ureaPAGEcDNAsizeexcision; RNA/remnants removed onlytoextentofphysicalseparation')
    cat=g.input('circligase','100UCircLigaseinitial+100Uafter1h',K.CATALYST);y=g.product('circle_product','Circularizedss-cDNA plus uncircularizedtemplates/reagents; fulljunctionunassembled',K.MOLECULAR)
    g.react('circle','CircLigase60C2h;additional100Uat1h',S.CIRCULARIZATION,((x,R.TEMPLATE),(cat,R.CATALYST)),(y,),geometry=CircularizationGeometry.unavailable('Source cDNA circle; full substrate/end sequence geometry unassembled.'));x=y
    x=clean(g,'circle_clean',x,'Source inherited isopropanol cDNA precipitation beforePCR; postcircle rRNAcapture omitted forRend')
    x=op(g,'pcr',x,'Phusion/o231/indexingprimers;6–10PCRcycles',S.PCR,reagent='SourcePhusion;S4o231andindexingprimersunreviewed',product='SourceamplifiedDNAlibrarymixtureplusprimerdimer/carryover; circularRNA not the substrate')
    x=clean(g,'pcr_clean',x,'8percentTBPAGEsizeexcision')
    return finish(g,x,'Published Rend HiSeq2000 orNextSeq500; runlength/primergeometry unassembled')
def protocol_graph():return protocol_graph_for()
def selftest():
    from checks import Check
    c=Check();a=protocol_graph();b=protocol_graph_for(True)
    c('Baseline rRNA processing precedes fragmentation',a.steps['preheat'].inputs[0].material,'rrna_clean_kept')
    c('Exonuclease route omits rRNAremoval','rrna'in b.steps,False)
    c('Exo route actually cleaned',b.steps['preheat'].inputs[0].material,'exo_clean_kept')
    c('15–45ntphysical selection precedesPNK',a.steps['pnk'].inputs[0].material,'short_kept')
    c('Ligated RNA physicallypurifiedbeforeRT',a.steps['rt'].inputs[0].material,'lig_clean_kept')
    c('cDNAcircle followsRT-PAGE',a.steps['circle'].inputs[0].material,'rt_clean_kept')
    c('PCR follows inherited precipitation',a.steps['pcr'].inputs[0].material,'circle_clean_kept')
    c('No postcircle rRNA inheritance','post_circle_rrna'in a.steps,False)
    c('No inherited V.natriegens SMARTerfinish','smarter'in a.steps,False)
    c('Sourcegaps prevent publicendpoint',bool(a.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
