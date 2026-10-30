"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment, revcomp
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def oligo(g, key, name, bases, **mods):
    parts = []
    for i, block in enumerate(bases.split('[BtndT]')):
        if i: parts.append(Segment('Internal biotin-dT', 'T', note='Source [BtndT]; physical modification not inherited by copies'))
        if block: parts.append(Segment(key + str(i), block, placeholder=bool(set(block)-set('ACGTU')), note='IUPAC synthesis pattern retained; full landing unassembled'))
    sc = Scene(); sc.strand(key, parts, **mods)
    return g.input(key, name, K.OLIGO, state=sc)

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=()):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION, S.DEPHOSPHORYLATION}
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

def finish(g, x, label, limitation):
    lib = g.product('library', label, K.LIBRARY)
    g.react('library_identity', 'Physical library interface; complete molecular state unassembled', S.GENERIC, ((x, R.TEMPLATE),), (lib,), changes_molecule=False)
    g.sequence(lib, SequencingPlan(limitation=limitation))
    g.add_source_item(SourceItem('pending', 'Historical exact kit, substrate and sequencing geometry; fixed A/B inline-control feature spans unassembled', 'Private partial graph; no full source audit or accepted absence.'))
    return g

SOURCE = 'Associated AM PMC3705217 Methods217–290; published SI MOESM98 overview23446–23545/Table2 p39; AM SI selected comparison only. Main final-version applicability unreviewed.'
CAP_DNA = 'CACTCT'
CAP_RNA = 'GAGCAAUACC'
RT_A = 'TATAGCGGCCGCCTGGAGTTTTTTTTTTTTTTTTVN'
RT_B = 'TATAGCGGCCGCCTCCACTTTTTTTTTTTTTTTTVN'
BIO_A = 'TATAGCGGCCGCTA[BtndT]CACTCTGAGCAATACC'
BIO_B = 'TATAGCGGCCGCTA[BtndT]ATCACTCTGAGCAATACC'
AMP_A = 'CATGTATAGCGGCCGCCTGGAG'
AMP_B = 'ACATGTATAGCGGCCGCCTCCAC'

@dataclass(frozen=True)
class CappingOligo:
    dna: str = CAP_DNA
    rna: str = CAP_RNA
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or not self.rna or set(self.rna)-set('ACGU'):
            raise ValueError('Ordered explicit DNA and RNA blocks required')
    def scene(self):
        sc = Scene(); sc.strand('capping', [Segment('DNA block', self.dna), Segment('RNA block', self.rna)]); return sc
    def dna_copy(self):
        return revcomp((self.dna + self.rna).replace('U', 'T'))
    def flatten(self):
        raise ValueError('Mixed substrate cannot be represented as ordinary DNA')

def branch(g, key, x):
    rt = oligo(g, 'rt_' + key, 'Published Table2 3cDNANotI_' + key + '; tube-specific inline-control feature spans unassembled', RT_A if key == 'a' else RT_B)
    x = op(g, 'rt_' + key + '_copy', x, 'Associated AM oligo-dT SSIII RT42/50/55C', S.REVERSE_TRANSCRIPTION, reagent='SuperScriptIII/dNTPs/DTT/RNase inhibitor', structural=(rt,), product='First-strand DNA plus original RNA/capping DNA/RNA oligo and unextended primers')
    x = op(g, 'rnase_' + key, x, 'RNase cocktail and RNaseH37C30min', S.DIGESTION, reagent='Source RNase cocktail/RNaseH', product='cDNA plus RNA hydrolysates and DNA oligo/primer/reagent remnants')
    x = clean(g, 'cdna_' + key, x, 'Source AMPure XP recovery')
    bi = oligo(g, 'bio_' + key, 'Published Table2 internally biotinylated 5BioNotI_' + key, BIO_A if key == 'a' else BIO_B)
    amp = oligo(g, 'amp_' + key, 'Published Table2 3AmpNotI_' + key, AMP_A if key == 'a' else AMP_B)
    x = op(g, 'pcr_' + key, x, 'Source Phusion10-cycle full-length cDNA PCR', S.PCR, reagent='Phusion HF mixture; paired tube-specific barcode primers', structural=(bi, amp), product='Barcoded biotinylated dsDNA PCR products plus original templates/primers/reagents')
    return clean(g, 'pcr_clean_' + key, x, 'Source AMPure XP recovery before aliquot pooling')

def protocol_graph():
    g = ProtocolGraph('TIF-seq — selected original2013 path, AM main qualified', authority=SOURCE, representation='published_instance')
    x = g.input('total_rna', 'Source DNA-free total yeast RNA plus capped/polyadenylated spike-ins', K.MOLECULAR, primary=True)
    x = op(g, 'sap', x, 'SAP dephosphorylates non-capped RNA37C30min', S.DEPHOSPHORYLATION, reagent='Shrimp alkaline phosphatase/RNase inhibitor', product='RNA classes with exposed noncapped ends dephosphorylated; capped RNA retained; no purity claim')
    x = clean(g, 'sap_clean', x, 'Double phenol extraction/ethanol precipitation')
    x = op(g, 'tap', x, 'TAP decapping37C1h', S.END_REPAIR, reagent='Tobacco acid pyrophosphatase', product='Formerly capped RNA5prime phosphate classes plus original RNA/cleavage/reagent remnants')
    x = clean(g, 'tap_clean', x, 'Phenol:chloroform extraction/ethanol precipitation')
    cap = g.input('cap_oligo', 'Published Table2 dCdAdCdTdCdTrGrArGrCrArArUrArCrC', K.ADAPTER, state=CappingOligo().scene(), chemical_state=CappingOligo())
    x = op(g, 'cap_ligation', x, 'Source T4RNL1 mixed DNA/RNA oligo ligation16C overnight', S.LIGATION, reagent='T4RNA ligase1/DMSO/RNase inhibitor', structural=(cap,), product='Selected mixed-oligo-tagged RNA classes plus unligated RNA/oligo/reagents')
    x = clean(g, 'rna_column', x, 'RNeasy column purification')
    x = qc(g, 'rna_integrity', x, 'Source Bioanalyzer integrity assay')
    a = g.product('aliquot_a', 'Separate RNA aliquot for A barcode branch', K.MOLECULAR)
    b = g.product('aliquot_b', 'Separate RNA aliquot for B barcode branch', K.MOLECULAR)
    g.react('split_barcodes', 'Physical RNA split before separate RT and PCR', S.SPLIT, ((x,R.TEMPLATE),),(a,b))
    a,b = branch(g,'a',a),branch(g,'b',b)
    a=op(g,'noti_a',a,'PublishedSI NotIcohesiveendgeneration in A preparation beforepooling',S.DIGESTION,reagent='SourceNotI; exactconditions/workupnotassignedacrossmainversions',product='NotIdigestedAclassesplusfreeendpieces/uncutDNA/reagents; intermediateworkupunassembled')
    b=op(g,'noti_b',b,'PublishedSI NotIcohesiveendgeneration in B preparation beforepooling',S.DIGESTION,reagent='SourceNotI; exactconditions/workupnotassignedacrossmainversions',product='NotIdigestedBclassesplusfreeendpieces/uncutDNA/reagents; intermediateworkupunassembled')
    x = g.product('pooled', 'Pooled separately barcoded A-A and B-B cDNA preparations', K.MOLECULAR)
    g.react('pool', 'PublishedSI pools sticky-ended preparations afterNotI; AMmainpool-firstorderingnotadopted', S.POOL, ((a,R.TEMPLATE),(b,R.TEMPLATE)),(x,))
    x = op(g,'circle',x,'Source dilute intramolecular T4DNA ligation16C16h',S.LIGATION,reagent='T4DNA ligase,600ul reaction',product='Circularized cDNA classes plus linear substrates/intermolecular chimeras/end pieces/reagents')
    x = op(g,'exonuclease',x,'ExonucleaseI/III37C20min and EDTA/heat inactivation',S.DIGESTION,reagent='ExonucleasesI/III',product='Circular DNA-enriched reaction plus degraded linear DNA fragments and residual uncleaved classes')
    x = clean(g,'circle_clean',x,'Source phenol:chloroform/ethanol precipitation')
    x = op(g,'sonicate',x,'CovarisS220 sonication',S.FRAGMENTATION,reagent='Source sonication conditions',product='Junction-bearing biotin fragments plus other physical DNA fragments; no computational filter at this step')
    x = clean(g,'fragment_clean',x,'Source AMPure recovery')
    beads = g.input('beads','StreptavidinM280 beads',K.SELECTION)
    bound = g.product('bound','Biotin-associated DNA/bead fraction with unquantified carryover',K.MOLECULAR)
    unbound = g.product('unbound','Physical unbound capture fraction',K.WASTE)
    g.react('capture','M280biotin capture30min room temperature',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,unbound));g.discard(unbound,'Actual unbound separation, not perfect junction classification')
    x = clean(g,'wash',bound,'Source bead washing')
    x = op(g,'nebnext',x,'Source NEBNext/Illumina DNA-seq end preparation and forked-adapter installation',reagent='HistoricalNEBNextMasterMixes/barcodedforkedadapters; detailed substeps unassembled',product='Source bead-associated adapted library classes and unresolved reaction companions')
    x = op(g,'final_pcr',x,'Source Phusion20-cycle library enrichment',S.PCR,reagent='HistoricalPhusion/Illuminalibraryprimers')
    x = clean(g,'size',x,'Source eGel2percent SizeSelect300bp recovery')
    return finish(g,x,'Physical TIF-seq junction-fragment DNA library; full oligo landing unassembled','Associated AM and published SI state HiSeq2000 paired-end105bp; exact primer/index geometry unassembled. Chimera filtering follows acquisition and never physically removes library molecules.')

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check(); a=CappingOligo(); c('Mixed capping copy preserves reverse block order',a.dna_copy(),revcomp(CAP_RNA.replace('U','T'))+revcomp(CAP_DNA))
    try:a.flatten(); rejected=False
    except ValueError:rejected=True
    c('Mixed capping substrate rejects ordinary DNA collapse',rejected,True)
    g=protocol_graph();c('Barcodes introduced after physical split',g.steps['rt_a_copy'].inputs[0].material,'aliquot_a')
    c('Second barcode branch uses independent RNA',g.steps['rt_b_copy'].inputs[0].material,'aliquot_b')
    c('Exonuclease does not directly discard linear fragments',g.steps['exonuclease'].outputs,('exonuclease_product',))
    c('Biotin capture explicitly separates physical unbound fraction',g.terminals['unbound'].disposition,Disposition.DISCARD)
    c('PublishedSIpoolingusesdigestedApreparation',g.steps['pool'].inputs[0].material,'noti_a_product')
    c('PublishedSIpoolingusesdigestedBpreparation',g.steps['pool'].inputs[1].material,'noti_b_product')
    c('Partial source graph cannot publish',bool(g.publication_problems()),True)
    src=Source(sources_root()/'tif-seq__10.1038+nature12121'/'nature12121_41586_2013_BFnature12121_MOESM98_ESM.pdf.txt','get_sources --slugs tif-seq__10.1038+nature12121')
    if have(c,src,label='published Table2'):
        c('Published mixed-sugar capping transcription','dCdAdCdTdCdTrGrArGrCrArArUrArCrC' in src.path.read_text(),True)
    c.report()
if __name__=='__main__':selftest()
