"""TruSeq PCR-Free v00 workflow; explicitly selected v22 CD adapter reference."""
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from endprep import EndPreparedDuplex, dA_tailed_scene
from protocol_graph import Disposition, InputRole as R, MaterialKind as M, ProtocolGraph, SequencingPlan, StepKind as K, Terminal

TITLE = 'TruSeq DNA PCR-Free'
NOTES = '01_truseq_dna_pcr_free.html'
SOURCE = ('Illumina <a href="https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-dna-pcr-free-workflow/truseq-dna-pcr-free-workflow-reference-1000000039279-00.pdf">TruSeq DNA PCR-Free Reference Guide 1000000039279 v00 (October 2017)</a>, '
          'Repair Ends and Select Library Size through Check Libraries; '
          '<a href="https://support-docs.illumina.com/SHARE/AdapterSequences/1000000002694_22_ilmn-adapter-sequences.pdf">Illumina Adapter Sequences 1000000002694 v22</a>, TruSeq DNA and RNA CD Indexes, pp.79–80.')
SUMMARY = 'Covaris shearing, end repair, two-sided size selection, dA-tailing and indexed adapter ligation; no library enrichment PCR.'
CAVEAT = ('The workflow is v00; the separately selected CD adapter sequence reference is v22. '
          'This reconstruction does not identify a 2017 kit lot or a sample well. '
          'The eight-base physical index windows do not specify acquired index cycles. '
          'The sealed product retains two forks; flow-cell copying is outside this preparation graph.')
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = ('The guide quantifies cluster-capable molecules by qPCR. Sequencing primer landing '
                     'on the forked substrate, flow-cell copying and the acquired read recipe are not reconstructed here.')

# v22 prints these two arms. Canonical components supply their fixed bases.
# Role windows preserve physical length without selecting a well or its orientation on a sample sheet.
def arms(prefix=''):
    u = [Segment(prefix+'P5', il.P5),
         Segment(prefix+'i5', 'J'*8, placeholder=True, feature=feature('truseq_cd_i5', 'sample_index', 'combinatorial')),
         Segment(prefix+'read1 fork', il.TRUSEQ_READ1[:-1])]
    i = [Segment(prefix+'read2 arm', il.INDEX1_PRIMER),
         Segment(prefix+'i7', 'I'*8, placeholder=True, feature=feature('truseq_cd_i7', 'sample_index', 'combinatorial')),
         Segment(prefix+'P7 reverse complement', il.P7_RC)]
    lengths = [n for n in range(1, min(len(u[-1]), len(i[0]))+1)
               if revcomp(i[0].top[:n]) == u[-1].top[-n:]]
    n = max(lengths)  # computed antiparallel terminal overlap, not an inherited stem default
    u[-1:] = [Segment(prefix+'read1 unpaired arm', u[-1].top[:-n]),
              Segment(prefix+'universal stem', u[-1].top[-n:]), Segment(prefix+'3T', 'T')]
    i[:1] = [Segment(prefix+'indexed stem', i[0].top[:n]),
             Segment(prefix+'read2 unpaired arm', i[0].top[n:])]
    return u, i


def adapter_scene():
    u, i = arms()
    sc = Scene()
    sc.strand('universal', u, label='v22 CD i5 arm')
    sc.anneal('indexed', i, to='universal', pair=('indexed stem', 'universal stem'),
              unpaired=tuple(s.name for s in i[1:]), label='v22 CD i7 arm')
    sc.overhang('universal', '3T', polarity='3-prime', owner='adapter', text='3′ dT')
    return sc


def insert_scene(label):
    return Scene.duplex([Segment('insert', 'N'*36, placeholder=True)], label=label)


def forked_product():
    u, i = arms('top ')
    bu, bi = arms('bottom ')
    body = Segment('genomic insert', 'N'*36, placeholder=True)
    other = complement_segments([body])[0]
    top = [*u, body, Segment('top insert 3A', 'A'), *i]
    bottom = [*bu, other, Segment('bottom insert 3A', 'A'), *bi]
    sc = Scene(); sc.strand('top', top, label='sealed PCR-free strand')
    # Only the two adapter stems, insert and dA/dT pairs are duplex. The platform tails remain forks.
    unpaired = tuple(s.name for s in bottom if s not in (bu[-2], bu[-1], other, bottom[len(bu)+1], bi[0]))
    sc.anneal('bottom', bottom, to='top', pair=(other.name, body.name),
              unpaired=unpaired, label='opposed sealed strand')
    for strand, uni, idx, bname, aname in (
            ('top', u, i, body.name, 'top insert 3A'),
            ('bottom', bu, bi, other.name, 'bottom insert 3A')):
        sc.junction(strand, uni[-1].name, bname, 'ligation')
        sc.junction(strand, aname, idx[0].name, 'ligation')
        sc.component(strand, uni[0].name, uni[-1].name, 'CD adapter')
        sc.component(strand, idx[0].name, idx[-1].name, 'CD adapter')
    return sc


def oligos():
    u, i = arms()
    return ['v22 CD i5 arm, 5′–'+''.join(s.top for s in u).replace('J'*8, '[i5]')+'–3′',
            'v22 CD i7 arm, 5′–'+''.join(s.top for s in i).replace('I'*8, '[i7]')+'–3′']


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input('dna', 'high-molecular-weight genomic DNA', M.BIOLOGICAL,
                  state=insert_scene('genomic DNA; drawn length schematic'), primary=True)
    def change(identifier, label, kind, before, state, reagents, final=False):
        ports = [(before, R.TEMPLATE)]
        for rid, rlabel, rkind, role in reagents:
            ports.append((g.input(rid, rlabel, rkind), role))
        out = g.product(identifier+'_product', label, M.LIBRARY if final else M.MOLECULAR, state=state)
        g.react(identifier, label, kind, ports, (out,))
        return out
    frag = change('shear', 'Covaris-sheared DNA', K.FRAGMENTATION, dna, insert_scene('sheared DNA'),
                  [('covaris', 'Covaris shearing; 350/550-bp workflow choices', M.CONDITION, R.CONDITION)])
    blunt = change('repair', 'ERP2/ERP3 blunt end repair', K.END_REPAIR, frag, insert_scene('blunt repaired DNA'),
                   [('erp', 'End Repair Mix 2 or 3', M.CATALYST, R.CATALYST)])
    def recover(identifier, label, before, state, beadlabel, waste_label, final=False):
        bead = g.input(identifier+'_beads', beadlabel, M.SELECTION)
        kept = g.product(identifier+'_kept', label, M.LIBRARY if final else M.MOLECULAR, state=state)
        waste = g.product(identifier+'_waste', waste_label, M.WASTE)
        g.react(identifier, label, K.SELECTION, ((before,R.TEMPLATE),(bead,R.SELECTION)), (kept,waste))
        g.discard(waste, waste_label)
        return kept
    smaller = recover('upper_cutoff', 'Transfer diluted-SPB supernatant', blunt, blunt.state,
                       'diluted SPB', 'oversized fragments remain on discarded beads')
    selected = recover('lower_cutoff', 'Retain DNA bound to undiluted SPB; wash and elute', smaller, blunt.state,
                        'undiluted SPB', 'small fragments in discarded supernatant and wash fractions')
    prep = EndPreparedDuplex(Construct([Segment('insert','N'*36,placeholder=True)]),
                             five_prime_phosphate=True, three_prime_overhang='A')
    tailed = change('tail', 'Add complementary 3′ dA ends', K.TAILING, selected, dA_tailed_scene(prep),
                    [('atl', 'ATL/ATL2 with CTA', M.CATALYST, R.CATALYST)])
    adapter = g.input('adapter', 'v22 CD reference adapter; sample well unassigned', M.ADAPTER, state=adapter_scene())
    ligase = g.input('ligase', 'LIG2 with CTL', M.CATALYST)
    sealed = g.product('sealed', 'double-ended sealed product retaining forked arms', M.MOLECULAR, state=forked_product())
    other_products = g.product('nonclustering_products', 'DNA lacking two productive adapter ends; source-noted nonclustering fraction', M.MOLECULAR)
    geometry = AdapterLigationGeometry.known('cohesive', dA_tailed_scene(prep), adapter.state,
                   statement='v00: adapter 3′ dT pairs with insert 3′ dA; v22 CD arms explicitly selected.')
    g.react('ligate', 'Ligate at 30 °C for 10 minutes', K.LIGATION,
             ((tailed,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)), (sealed,other_products), geometry=geometry)
    g.mark_terminal(Terminal(other_products.id, Disposition.UNRESOLVED,
        'v00 Check Libraries notes nonclustering molecules lacking adapters at both ends. Their recovery through the two SPB cleanups is unassigned; qPCR distinguishes quantities analytically and does not physically remove them.'))
    stopped = change('stop_ligation', 'Stop ligation with STL', K.QUENCHING, sealed, sealed.state,
                      [('stl','Stop Ligation Buffer',M.CONDITION,R.CONDITION)])
    one = recover('postligation_round1', 'First SPB cleanup, ethanol wash and RSB elution', stopped, sealed.state,
                   '42.5 µl SPB, then 52.5 µl RSB', 'first cleanup discarded supernatant and wash fractions')
    two = recover('postligation_round2', 'Second SPB cleanup, ethanol wash and RSB elution', one, sealed.state,
                   '50 µl SPB, then 22.5 µl RSB', 'second cleanup discarded supernatant and wash fractions', final=True)
    g.sequence(two, SequencingPlan(limitation=SEQUENCING_ENDING, acquisition='unspecified'))
    g.declare_reviewed_scope(('repair','upper_cutoff','lower_cutoff','tail','adapter','ligate',
                              'stop_ligation','postligation_round1','postligation_round2', 'dna','covaris','shear',
                              'erp','upper_cutoff_beads','lower_cutoff_beads','atl','ligase','stl',
                              'postligation_round1_beads','postligation_round2_beads','postligation_round2_kept'), SOURCE)
    return g.require_publishable()
