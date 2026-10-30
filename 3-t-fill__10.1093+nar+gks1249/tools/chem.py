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

SOURCE='PublishedPMC3597643 Methods46–59/Fig1; correctionPMC3695521 changesculture sugar1to2percent. SupplementTable1 reviewed; detailedMethods/TableS3unread within3originalcap.'
def protocol_graph():
    g=ProtocolGraph('3primeT-fill — selected published main path',authority=SOURCE,representation='published_instance')
    x=g.input('rna','10ugDNaseI-treatedtotalyeastRNAwithsourceIVTspikecontrols',K.MOLECULAR,primary=True)
    x=op(g,'fragment',x,'RNA80C5minMg-acetatefragmentation',S.FRAGMENTATION,reagent='Tris-acetate/KOAc/MgOAcsourcebuffer')
    x=clean(g,'rna_clean',x,'1.5xAMPure;12.8ulelution')
    rtin=g.product('rt_input','11.2ulof12.8uleluate',K.MOLECULAR);rem=g.product('rna_remainder','Remaining1.6uleluate;fateunassigned',K.MOLECULAR)
    g.react('rna_aliquot','ActualRTinputfraction',S.SPLIT,((x,R.TEMPLATE),),(rtin,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Unassignedremainder'))
    primer=g.input('rt_primer','BiotinylatedP5_dT16VN;wholeTableS3sequenceunreviewed',K.OLIGO,state=scene('rt_primer','[source biotinylated P5-dT16VN]'))
    x=op(g,'rt',rtin,'SSII42C50min/70C15min;actinomycinD/RNasin',S.REVERSE_TRANSCRIPTION,reagent='SSII/dNTP/DTT/actinomycinD/RNasin',structural=(primer,))
    x=clean(g,'rt_clean',x,'1.5xAMPure;40ulelution')
    x=op(g,'second_strand',x,'RNaseH/DNApolymeraseI16C2.5h',S.EXTENSION,reagent='RNaseH/DNApolymeraseI/dNTPbuffer',product='BiotinylatedduplexDNAclasses plusRNAfragments/unextendedsubstrates/reagents')
    x=clean(g,'second_clean',x,'0.8xAMPure;20ulelution')
    beads=g.input('beads','SourceM280streptavidin20ul',K.SELECTION);b=g.product('bound','Biotin-associatedDNA/beadfractionwithcarryover',K.MOLECULAR);w=g.product('unbound','Actualphysicalunboundfraction',K.WASTE)
    g.react('capture','M280capture15min25C',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(b,w));g.discard(w,'Actualunboundseparation')
    x=clean(g,'capture_wash',b,'TwoB&Wwashes/EBwash')
    x=op(g,'end_repair',x,'NEBNextend-repair20C30min',S.END_REPAIR,reagent='HistoricalNEBNextDNAprepMasterMixSet1')
    x=clean(g,'end_wash',x,'TwoB&Wwashes/EBwash')
    x=op(g,'a_tail',x,'Klenowexo-minus/dATP37C30min',S.TAILING,reagent='Klenowexo-minus/NEBuffer2/dATP')
    x=clean(g,'tail_wash',x,'TwoB&Wwashes/EBwash')
    x=op(g,'ligate',x,'P7_T1_Mpxlinker/T4DNAligase20C15min',S.LIGATION,reagent='T4DNAligase/Quickligationbuffer')
    x=clean(g,'ligation_wash',x,'FourB&Wwashes/EBwash;50ulresuspension')
    prod=g.product('pcr_beads','24ulof50ulbeadresuspension',K.MOLECULAR);left=g.product('bead_remainder','26ulbeads;fateunassigned',K.MOLECULAR)
    g.react('pcr_aliquot','Actual24of50ulPCRinput',S.SPLIT,((x,R.TEMPLATE),),(prod,left));g.mark_terminal(Terminal(left.id,Disposition.UNRESOLVED,'Unassignedremainder'))
    x=op(g,'pcr',prod,'Phusion/PE1.0/PE2.0 eighteenPCRcycles',S.PCR,reagent='PhusionMasterMix2x/sourcePE1.0/PE2.0')
    x=clean(g,'pcr_clean',x,'1.8xAMPurefollowedbyeGel2percentSizeSelect300bp')
    x=op(g,'cluster',x,'LoadfinallibrariesontocBot;clusterpreparationinternalsunassembled',reagent='SourcecBotflowcell/clusterreagents',product='Surface-boundclustertemplatesplusphysicalcompanions;no liquidlibraryidentitysubstitution')
    x=op(g,'t_fill',x,'Replaceprimingbuffer;genomicV2primerannealanddT-onlyextension60C5min',S.EXTENSION,reagent='SourceunlabeleddTTP/nonhotstartPhusion/TaqE/HFbuffer',structural=(g.input('seq_primer','GenomicDNASequencingprimerV2;historicalwholesequenceunassembled',K.OLIGO,state=scene('seq_primer','[source genomic V2]')),),product='OnclusterprimerextendedthroughpolyAstretch;remainingtemplates/reagentsretained;stop/slippageefficiencynotpredicted')
    lib=g.product('library','PhysicalT-filledsurface-boundclustermaterial',K.LIBRARY);g.react('library_identity','Physicalacquisitioninterface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='PublishedHiSeq2000;barcodeviaoppositepaired-endread. SourceTableS3andrungeometryunreviewed;noexactread-cycleclaim. Mapping/polyA/internalprimingfilters occurafteracquisition.'))
    g.add_source_item(SourceItem('pending','DetailedsupplementTableS3/historicalkit/readgeometry','Sourceunreviewedwithin3originalcap;noabsenceclaim.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('RTusesactual11.2ulaliquot',g.steps['rt'].inputs[0].material,'rt_input')
    c('PCRusesactual24ulaliquot',g.steps['pcr'].inputs[0].material,'pcr_beads')
    c('Unusedbeadsremainunresolved',g.terminals['bead_remainder'].disposition,Disposition.UNRESOLVED)
    c('Tfillusesclusteredsubstrate',g.steps['t_fill'].inputs[0].material,'cluster_product')
    c('Actualcaptureunboundseparated',g.terminals['unbound'].disposition,Disposition.DISCARD)
    c('Incompleteendpointnotpublishable',bool(g.publication_problems()),True)
    src=Source(sources_root()/'3-t-fill__10.1093+nar+gks1249'/'PMC3695521_PMC3695521.xml.txt','get_sources --slugs 3-t-fill__10.1093+nar+gks1249')
    if have(c,src,label='formalcorrection'):c('Correctionaffectsculturesugar','2% glucose and galactose' in src.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
