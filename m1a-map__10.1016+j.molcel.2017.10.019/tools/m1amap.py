"""m1A-MAP source-owned input/IP/AlkB branches and library operations.

Associated manuscript scope only. No mutation frequency, modification identity
at a native site, or final full strand/read-cycle geometry is predicted.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp, complement, panel, workflow_from_protocol_graph, workflow_panel
from page import head, info
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, SourceItem

PAGE_FILE = 'm1a-map.html'
SOURCE = ('Li et al.2017 DOI10.1016/j.molcel.2017.10.019; associated accepted '
          'manuscript PMC5722686 Methods178–180, processing266–272 and Figure1')
# Exact local printed sequences match these canonical sequence parts; no eCLIP
# recipe or reaction defaults are imported from the cited predecessor.
RIL19 = revcomp(il.TRUSEQ_READ1[-19:])
RT_PRIMER = il.TRUSEQ_READ1[-20:]
CDNA_CORE = revcomp(il.TRUSEQ_READ2[-22:])
PCR_F = il.TRUSEQ_P5_FULL
PCR_R_TAIL = il.TRUSEQ_READ2[:-1]
READ_LENGTH = 150


@dataclass(frozen=True)
class SelectedRTOutcome:
    """Caller-selected reference-oriented base outcome, never a rate predictor.

    RNA's modification state does not change merely because TGIRT copies it.
    None means arrest/no incorporated base at this selected site. A control RNA
    may also miscopy; no perfect-copy rule is imposed on the minus-m1A class.
    """
    rna_status: str
    reference_base: str | None

    def __post_init__(self):
        if self.rna_status not in {'A', 'm1A'} or self.reference_base not in {None, 'A', 'C', 'G', 'T'}:
            raise ValueError('Selected event needs an explicit RNA class and base/arrest outcome')

    @property
    def physical_cdna_base(self):
        return None if self.reference_base is None else complement(self.reference_base)


def cdna_adapter(group):
    return Construct([Segment('random molecular barcode', 'N' * 10, 'umi', placeholder=True,
                             feature=feature('m1amap.' + group + '.umi', 'umi', 'degenerate')),
                      Segment('source cDNA-adapter core', CDNA_CORE)])


def reverse_primer(group):
    return Construct([Segment('P7', il.P7),
                      Segment('sample index', 'N' * 6, 'idx', placeholder=True,
                              feature=feature('m1amap.' + group + '.index', 'sample_index', 'unknown')),
                      Segment('source PCR reverse suffix', PCR_R_TAIL)])


def adapter_scene(which, group='example'):
    sc = Scene()
    if which == 'rna':
        sc.strand('adapter', [Segment('RiL19 source body', RIL19)], mod5='rAPP', mod3='C3 spacer',
                  label='source RNA-linker ordering string; C3 is not an extra cytidine')
    elif which == 'cdna':
        sc.strand('adapter', list(cdna_adapter(group)), mod5='phosphate', mod3='C3 spacer')
    else:
        raise ValueError('Source adapters are rna or cdna')
    return sc


def rt_binding_scene():
    # This is the printed body/primer binding relation, not full RNA-end geometry.
    sc = Scene()
    sc.strand('adapter', [Segment('RiL19 source body', RIL19)], mod3='C3 spacer',
              label='adapter body at RNA 3′ end; upstream junction unassembled')
    sc.anneal('primer', [Segment('RT primer 5′ overhang', RT_PRIMER[:-19]),
                         Segment('RT primer paired arm', RT_PRIMER[-19:])], to='adapter',
              pair=('RT primer paired arm', 'RiL19 source body'), unpaired=('RT primer 5′ overhang',))
    return sc


def unassembled_ligation(which, group):
    sc = adapter_scene(which, group)
    sc.note('adapter', 'Printed reagent ends retained; complete substrate/junction geometry and historical Methods S1 remain unassembled.')
    return AdapterLigationGeometry('unavailable', 'unavailable', (sc,),
                                  'Source reagent known; exact RNA/cDNA substrate end geometry unreviewed, not asserted unreported.')


def protocol_graph():
    g = ProtocolGraph('m1A-MAP — source input/IP/demethylase library branches', authority=SOURCE,
                      representation='published_instance')
    reviewed = []

    def inp(key, label, kind, **kwargs):
        value = g.input(key, label, kind, **kwargs); reviewed.append(key); return value

    def react(key, label, kind, inputs, outputs, **kwargs):
        value = g.react(key, label, kind, inputs, outputs, **kwargs); reviewed.append(key); return value

    def cleanup(key, material, label, procedure):
        condition = inp(key + '_procedure', procedure, K.SELECTION)
        kept = g.product(key, label, K.MOLECULAR)
        reject = g.product(key + '_waste', 'Physical purification rejects; exact composition/recovery unquantified', K.WASTE)
        react(key + '_cleanup', procedure, S.CLEANUP, ((material, R.TEMPLATE), (condition, R.SELECTION)), (kept, reject))
        g.discard(reject, 'Source-described physical separation; no perfect molecular-class removal asserted')
        return kept

    rna = inp('polya_rna', 'Source poly(A)+ RNA; upstream extraction/selection outside this graph', K.MOLECULAR, primary=True)
    magnesium = inp('fragmentation', 'Source magnesium RNA fragmentation to approximately150nt; main-MAP temperature/program unassigned', K.CONDITION)
    fragments = g.product('fragments', 'RNA-fragment population including unmodified/modified RNA; yields and ends unassigned', K.MOLECULAR)
    react('fragment', 'Fragment source40µg poly(A)+ RNA', S.FRAGMENTATION,
          ((rna, R.TEMPLATE), (magnesium, R.CONDITION)), (fragments,))
    input_rna = g.product('input_rna', 'Fragmented-RNA input aliquot, without IP', K.MOLECULAR)
    for_ip = g.product('for_ip', 'Source RNA allocated to anti-m1A enrichment; aliquot sizes unassigned', K.MOLECULAR)
    react('allocate_input', 'Represent source input and IP preparations separately', S.SPLIT, ((fragments, R.TEMPLATE),), (input_rna, for_ip))
    antibody = inp('antibody', 'Anti-m1A immunoprecipitation; cited Li2016 detailed procedure unreviewed', K.SELECTION)
    enriched = g.product('enriched', 'IP-enriched RNA, not a pure m1A-only population', K.MOLECULAR)
    unbound = g.product('ip_reject', 'Physically unretained IP fraction; composition unknown', K.WASTE)
    react('ip', 'Source m1A immunoprecipitation', S.AFFINITY_CAPTURE,
          ((for_ip, R.TEMPLATE), (antibody, R.SELECTION)), (enriched, unbound)); g.discard(unbound, 'Physical IP partition, not perfect modification selection')
    minus = g.product('minus_rna', 'Source IP RNA without demethylase treatment', K.MOLECULAR)
    plus = g.product('plus_rna', 'IP RNA allocated to AlkB treatment; source10ng', K.MOLECULAR)
    react('allocate_demethylase', 'Represent minus/plus demethylase preparations', S.SPLIT, ((enriched, R.TEMPLATE),), (minus, plus))
    alkb = inp('alkb', 'Source ΔN11 AlkB0.4pmol;20µL/MESpH6.5/FeII/alphaKG/ascorbate/RNase inhibitor,37°C2h', K.CATALYST)
    treated = g.product('treated', 'AlkB-treated RNA plus reaction products; residual/intermediate fractions unassigned', K.MOLECULAR)
    react('demethylate', 'Source main-MAP AlkB treatment; no complete-erasure assumption', S.CHEMICAL_CONVERSION,
          ((plus, R.TEMPLATE), (alkb, R.CATALYST)), (treated,))
    edta = inp('edta', 'Source5mM EDTA quench', K.CONDITION)
    stopped = g.product('stopped', 'Quenched AlkB reaction mixture', K.MOLECULAR)
    react('quench_alkb', 'Stop demethylase reaction', S.QUENCHING, ((treated, R.TEMPLATE), (edta, R.CONDITION)), (stopped,))
    treated = cleanup('demethylated_clean', stopped, 'Recovered plus-demethylase RNA', 'Phenol/chloroform extraction followed by ethanol precipitation')

    for group, template in [('input', input_rna), ('minus', minus), ('plus', treated)]:
        def key(name): return group + '_' + name
        pnk = inp(key('pnk'), 'Source PNK3′dephosphorylation37°C1h; heat-inactivate65°C20min', K.CATALYST)
        repaired = g.product(key('repaired'), 'Source end-treated RNA population; precise end chemistry unassembled', K.MOLECULAR)
        react(key('dephosphorylate'), 'Source RNA3′end dephosphorylation', S.DEPHOSPHORYLATION,
              ((template, R.TEMPLATE), (pnk, R.CATALYST)), (repaired,))
        repaired = cleanup(key('repaired_clean'), repaired, 'Ethanol-recovered RNA', 'Source ethanol precipitation after PNK heat inactivation')
        linker = inp(key('rna_linker'), 'RiL19:5′rAPP source body3′C3 spacer', K.ADAPTER, state=adapter_scene('rna', group))
        ligase = inp(key('ligase2'), 'T4 RNA ligase2 truncatedKQ,25°C2h', K.CATALYST)
        ligated = g.product(key('ligated_rna'), 'Adapter-ligated RNA plus excess adapter and reaction-derived material', K.MOLECULAR)
        react(key('rna_ligation'), 'Source3′RNA-linker ligation', S.LIGATION,
              ((repaired, R.TEMPLATE), (linker, R.STRUCTURAL), (ligase, R.CATALYST)), (ligated,), geometry=unassembled_ligation('rna', group))
        dead = inp(key('deadenylase'), '5′deadenylase30°C1h', K.CATALYST)
        deadened = g.product(key('deadened'), 'Reaction after adapter deadenylation; material remains in tube', K.MOLECULAR)
        react(key('deadenylate'), 'First source excess-adapter processing enzyme', S.CHEMICAL_CONVERSION,
              ((ligated, R.TEMPLATE), (dead, R.CATALYST)), (deadened,))
        recj = inp(key('recj'), 'RecJf37°C1h, then70°C20min heat inactivation', K.CATALYST)
        digested = g.product(key('adapter_digested'), 'RNA population and adapter-digestion products before precipitation', K.MOLECULAR)
        react(key('digest_adapter'), 'Second source excess-adapter processing enzyme', S.DIGESTION,
              ((deadened, R.TEMPLATE), (recj, R.CATALYST)), (digested,))
        recovered = cleanup(key('rna_clean'), digested, 'Source ethanol-recovered RNA', 'Ethanol precipitation after adapter digestion/inactivation')
        primer = inp(key('rt_primer'), 'Printed RT primer', K.OLIGO, state=Construct([Segment('RT primer', RT_PRIMER)]))
        heat = inp(key('rt_heat'), 'RNA plus RT primer80°C2min; immediately ice', K.CONDITION)
        heated_rna = g.product(key('heated_rna'), 'Heat-treated RNA in the same primer-containing mixture', K.MOLECULAR)
        heated_primer = g.product(key('heated_primer'), 'Co-incubated RT primer, not physically isolated', K.OLIGO, state=Construct([Segment('RT primer', RT_PRIMER)]))
        react(key('rt_denature'), 'Heat RNA with RT primer before adding TGIRT', S.THERMAL_TREATMENT,
              ((recovered, R.TEMPLATE), (primer, R.STRUCTURAL), (heat, R.CONDITION)), (heated_rna, heated_primer))
        rt = inp(key('tgirt'), 'TGIRT,57°C2h in source RT mixture; no efficiency/mismatch distribution predicted', K.CATALYST)
        cdna_mix = g.product(key('rt_mix'), 'cDNA products plus RNA, unextended primers and possible arrested products', K.MOLECULAR)
        react(key('reverse_transcribe'), 'TGIRT copying; RNA modification and DNA base outcomes remain distinct', S.REVERSE_TRANSCRIPTION,
              ((heated_rna, R.TEMPLATE), (heated_primer, R.STRUCTURAL), (rt, R.CATALYST)), (cdna_mix,))
        exoi = inp(key('exoi'), 'ExonucleaseI,37°C30min', K.CATALYST)
        after_exo = g.product(key('after_exo'), 'Reaction after excess-primer digestion, including RNA and soluble products', K.MOLECULAR)
        react(key('digest_primer'), 'Source RT-primer processing', S.DIGESTION,
              ((cdna_mix, R.TEMPLATE), (exoi, R.CATALYST)), (after_exo,))
        cdna = cleanup(key('cdna_clean'), after_exo, 'Silane-recovered cDNA population; source RNA/carryover fate unquantified', 'Source silane-bead cDNA purification')
        adapter = inp(key('cdna_adapter'), '5′phosphorylated N10/core adapter,3′C3 spacer', K.ADAPTER, state=adapter_scene('cdna', group))
        ligase1 = inp(key('ligase1'), 'High-concentration T4 RNA ligase1,25°C overnight', K.CATALYST)
        adapted = g.product(key('adapted'), 'Source cDNA-adapter ligation products and residual reaction material', K.MOLECULAR)
        react(key('cdna_ligation'), 'Source cDNA3′adapter ligation', S.LIGATION,
              ((cdna, R.TEMPLATE), (adapter, R.STRUCTURAL), (ligase1, R.CATALYST)), (adapted,), geometry=unassembled_ligation('cdna', group))
        adapted = cleanup(key('adapted_clean'), adapted, 'Silane-recovered adapted cDNA', 'Source silane purification after cDNA ligation')
        forward = inp(key('pcr_f'), 'Printed unindexed forward PCR primer', K.OLIGO, state=Construct([Segment('PCR forward', PCR_F)]))
        reverse = inp(key('pcr_r'), 'Printed reverse PCR primer with unknown6base sample index', K.OLIGO, state=reverse_primer(group))
        polymerase = inp(key('pcr_mix'), 'Source PCR amplification; full program/polymerase in uninspected MethodsS1', K.CONDITION)
        amplified = g.product(key('amplified'), 'Library PCR population; full genomic strand/product geometry unassembled', K.MOLECULAR)
        react(key('pcr'), 'Amplify source library with printed primers', S.PCR,
              ((adapted, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (polymerase, R.CONDITION)), (amplified,),
              geometry=PCRGeometry.unavailable(Scene(), reason='Full RNA/cDNA strand assembly and historical detailed MethodsS1 remain unreviewed; known binding subrelations are not complete PCR validation.'))
        library = cleanup(key('library'), amplified, '8%TBE-gel recovered library', 'Source8%TBE-gel purification')
        g.sequence(library, SequencingPlan(acquired_roles=('Read 1', 'Read 2'), read_lengths={'Read 1':READ_LENGTH, 'Read 2':READ_LENGTH},
                                         limitation='Associated manuscript reports HiSeqX10 PE150. Methods268–272 computationally trim the first10bases of Read2 as the random barcode. Exact physical run-primer landing and index acquisition remain unassembled.'))
    g.declare_reviewed_scope(reviewed, SOURCE)
    g.add_source_item(SourceItem('detailed_source_and_geometry', 'Published MethodsS1, prior IP original and complete RNA/cDNA/PCR/read geometry',
                                'Actual cited originals/source-version applicability remain unresolved; no confirmed absence claimed.'))
    return g


def render_page():
    g = protocol_graph()
    return ('<!doctype html><html lang="en">' + head('m1A-MAP — source-limited model') +
            '<main class="wrap"><h1>m1A-MAP</h1><p><a href="01_m1a-map.html">Research notes</a></p>' +
            info('Private incomplete associated-manuscript model. Input/IP-minus/IP-plus are separate; IP is not perfectly pure and AlkB is not assumed to erase every modification. '
                 'RNA modification, reference-oriented mutation call and physical cDNA base are different quantities. Full MethodsS1 and final strand/read geometry remain unassembled.') +
            panel(rt_binding_scene().rows(), caption='Computed printed linker/RT-primer binding; one unpaired5′primer base, not a repaired source sequence.') +
            ''.join('<h2>' + name + '</h2>' + workflow_panel(workflow_from_protocol_graph(g, name + '_library'))
                    for name in ('input', 'minus', 'plus')) + '</main></html>')


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c = Check()
    source = Source(sources_root() / HERE.parent.name / 'j.molcel.2017.10.019_PMC5722686.xml.txt', 'Associated manuscript Methods180')
    if have(c, source, label='m1A-MAP printed library oligos'):
        text = ''.join(source.path.read_text().split())
        for label, seq in [('RNA linker',RIL19),('RT primer',RT_PRIMER),('cDNA adapter','N'*10+CDNA_CORE),('forward PCR',PCR_F),('reverse PCR',il.P7+'X'*6+PCR_R_TAIL)]:
            c('Original transcription '+label, seq in text, True)
    c('C3 spacer is not counted as an extra adapter nucleotide',len(RIL19),19)
    c('Printed RT5prime overhang retained while paired arm complements linker',
      (RT_PRIMER[:-19],RT_PRIMER[-19:]),('A',revcomp(RIL19)))
    c('Computed partial RT binding scene builds',len(rt_binding_scene().strands),2)
    mismatch = SelectedRTOutcome('m1A','G')
    c('Reference-oriented G outcome means physical cDNA C without erasing RNA m1A',
      (mismatch.physical_cdna_base,mismatch.rna_status),('C','m1A'))
    c('Arrest does not manufacture a DNA nucleotide',SelectedRTOutcome('m1A',None).physical_cdna_base,None)
    c('A control RNA is not constrained to perfect copying',SelectedRTOutcome('A','C').physical_cdna_base,'G')
    c('Library branches and unknown source keep publication blocked',protocol_graph().publishable,False)
    c.report()


if __name__ == '__main__':
    selftest()
