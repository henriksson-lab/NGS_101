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

SOURCE='Published PMC4782308 Methods150–190 and Supplementary Note A opening; selected E.coli route only.'
def capture(g,key,x):
    b=g.input(key+'_beads','30ul prewashed hydrophilic streptavidin beads; prewash2x400ulB/lowsalt then2x400ulA/highsalt',K.SELECTION)
    kept=g.product(key+'_bound','DTB-cap-associated RNA/beads with unquantified carryover',K.MOLECULAR);w=g.product(key+'_unbound','Physical nonretained fraction',K.WASTE)
    g.react(key,'Streptavidin capture20min roomtemperature',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(b,R.SELECTION)),(kept,w));g.discard(w,'Actual separation; no perfect cap classification.')
    kept=clean(g,key+'_wash',kept,'Postbinding physical washes:2x200ulA/highsalt then2x200ulB/lowsalt; retain beads')
    return clean(g,key+'_elute',kept,'30ul bufferB+1mMbiotin competitive elution20min roomtemperature; retain eluate and remove beads')
def protocol_graph_for(enriched=True):
    g=ProtocolGraph('Cappable-seq enriched'if enriched else'Cappable-seq independent unenriched control',authority=SOURCE,representation='published_instance')
    x=g.input('rna','3ug E.coli totalRNA after FastRNABlue/DNaseI/Megaclear; includes diverse5primeendclasses',K.MOLECULAR,primary=True)
    x=qc(g,'rin',x,'Separate Bioanalyzer integrity assay; sourceRIN9')
    x=op(g,'cap',x,'VCE50U/DTBGTP0.5mM/SAM0.1mM;50ul;37C30min',reagent='Source3prime-DTB-GTP/VCE/SAM',product='DTB-tagged caps on reactive5prime-di/triphosphorylatedRNAs plus otherendclasses/unreactedstock; capmethylation/exactgeometryunassembled')
    x=clean(g,'cap_clean',x,'Zymo RNA CleanConcentrator5 >=200nt,4washes;100ul lowTE eluate')
    x=op(g,'fragment',x,'100ulRNA+2.5ul10xPNKbuffer;94C5min',S.THERMAL_TREATMENT,reagent='Source thermal fragmentation buffer; not nuclease stock',product='Fragmented RNA with tagged original5primeends and newlygeneratedRNAends plus companions; fullenddistributionunassembled')
    x=clean(g,'fragment_clean',x,'180ul AMPure+420ulEtOH;2x80percentwash;100ul eluate')
    x=aliquot(g,'pnk_input',x,'75ul of100ul fragmented RNA for PNK;25ul remainder')
    x=op(g,'pnk',x,'75ulRNA+8.2ul10xbuffer+4ul ATP-free T4PNK;15min; temperatureunassigned',S.DEPHOSPHORYLATION,reagent='ATP-free T4PNK prepared by dialysis; no ATP-driven5primekinase assigned',product='Source3prime-phosphatase-treatedRNA mixture; unreacted/carryover retained')
    if enriched:
        x=aliquot(g,'capture_input',x,'50ul of PNK mixture to firstcapture; remainder unassigned')
        x=capture(g,'capture1',x)
        x=clean(g,'capture1_clean',x,'AMPure1.8x plus ethanol1.5x resultingmix;2x80percentwash;60ul eluate')
        x=aliquot(g,'capture2_input',x,'30ul of60ul firstcapture eluate tosecondcapture;30ul remainder')
        x=capture(g,'capture2',x)
        x=clean(g,'capture2_clean',x,'AMPure recovery;30ul lowTE eluate')
    else:
        x=op(g,'control_workup',x,'Independent control omits streptavidin capture; intermediate physical allocation/workup not reconstructed',reagent='Published control sameconditions except streptavidincapture omitted; no invented enriched-route aliquots',product='Unenriched physical control RNA for source downstream RppH; workup volumes unassembled')
    x=op(g,'rpph',x,'RppH15U/3.3ul10xThermopol+30ulRNA;37C60min;0.5ul0.5MEDTA then94C2min',S.DIGESTION,reagent='Source RppH; control reaction volume not independently reconstructed',product='RppH-treatedRNA including5prime-monophosphorylated classes suitablefor5primeRNAadapterligation; chemicalproducts/endgeometryunassembled')
    x=clean(g,'rpph_clean',x,'AMPure recovery to20ul')
    x=op(g,'kit',x,'Historical NEBNext SmallRNA library kit finishing',reagent='Historical source kit original and adapters unreviewed',product='Source DNA library precursor mixture; kitligation/RT/cleanupinternals unassembled')
    x=op(g,'pcr',x,'Source15PCRcycles',S.PCR,reagent='Source kit PCR reagents; exact primers unassembled')
    return finish(g,x,'Published E.coli MiSeqV3 single-end100nt; microbiomeGAII route not modeled')
def protocol_graph():return protocol_graph_for()
def selftest():
    from checks import Check
    c=Check();a=protocol_graph();b=protocol_graph_for(False)
    c('PNK actual75ul aliquot',a.steps['pnk'].inputs[0].material,'pnk_input_used')
    c('First capture actual50ul aliquot',a.steps['capture1'].inputs[0].material,'capture_input_used')
    c('Second capture actual halfeluate',a.steps['capture2'].inputs[0].material,'capture2_input_used')
    c('Capturewashes precede elution',a.steps['capture1_elute'].inputs[0].material,'capture1_wash_kept')
    c('RppH followssecondelution recovery',a.steps['rpph'].inputs[0].material,'capture2_clean_kept')
    c('Control avoids capture',any(k.startswith('capture')for k in b.steps),False)
    c('Control avoids enrichedaliquotinheritance','capture_input'in b.steps,False)
    c('UnusedRNA notidealizedaway',a.terminals['pnk_input_rest'].disposition,Disposition.UNRESOLVED)
    c('Unbound fraction physicallydiscarded',a.terminals['capture1_unbound'].disposition,Disposition.DISCARD)
    c('Kit gaps prevent public endpoint',bool(a.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
