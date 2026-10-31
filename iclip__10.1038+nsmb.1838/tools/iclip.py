"""Source-scoped partial iCLIP model; intentionally not a publishable library."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from pcr import PCRGeometry
from protocol_views import WorkflowSection,assemble_page
from page import head,info
from chemdraw import Construct, Scene, Segment, revcomp, feature, panel
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import (ProtocolGraph, MaterialKind, InputRole, StepKind,
                            Terminal, Disposition, MolecularRelation, SequencingPlan)
from circular import circularize_ssdna, CircularizationGeometry
from illumina import TRUSEQ_P5_FULL, PE_PCR_PRIMER_2

SOURCE = 'König et al. 2010 doi:10.1038/nsmb.1838, publisher Supplementary Methods pp1–4'
L3_WRITTEN = '5′-phosphate-UGAGAUCGGAAGAGCGGTTCAG-3′-Puromycin'
L3 = 'UGAGAUCGGAAGAGCGGTTCAG'
RT_PRIMERS = ('NNNCAAGATCGGAAGAGCGTCGTGGATCCTGAACCGCTC',
              'NNNGAAGATCGGAAGAGCGTCGTGGATCCTGAACCGC',
              'NNNTGAGATCGGAAGAGCGTCGTGGATCCTGAACCGCTC')
LINEARIZATION_PRIMER = 'GTTCAGGATCCACGACGCTCTTCAAAA'
PCR_PRIMERS = (PE_PCR_PRIMER_2, TRUSEQ_P5_FULL)
BOUNDARY = ('This capped model ends at the original replicate-1 cDNA circle. '
            'BamHI reopening, source-specific PCR reconstruction, variant discrepancies '
            'and exact acquired sequencing-primer geometry remain deferred.')


def linear_cdna():
    # The terminal nine bases pair L3 positions 12–20 (1-based); L3's terminal AG remains beyond the paired interval.
    complement = revcomp(L3.replace('U', 'T'))
    assert RT_PRIMERS[0][-9:] == complement[2:11]
    return Construct([
        Segment('random barcode', 'NNN', 'umi', placeholder=True,
                feature=feature('iclip_umi', 'umi', 'random')),
        Segment('experiment CA', 'CA', 'cbc',
                feature=feature('iclip_replicate', 'sample_index', 'fixed')),
        Segment('inverted PCR handles and BamHI site', RT_PRIMERS[0][5:-9]),
        Segment('RT landing suffix', RT_PRIMERS[0][-9:]),
        Segment('copied L3 upstream', revcomp(L3[:11].replace('U', 'T')),inferred=True),
        Segment('RNA copy to crosslink stop', '[variable cDNA]', placeholder=True,inferred=True),
    ], name='replicate-1 single-stranded cDNA; 5′ phosphate retained')



def l3_scene():
    sc = Scene()
    sc.strand('RNA adapter', [Segment('L3', L3)], label='L3',
              mod5='phosphate', mod3='Puromycin')
    sc.ligation_end('RNA adapter', 'L3', kind='single-strand',
                    polarity='5-prime', owner='adapter')
    return sc


def protocol_graph():
    """Selected replicate-1 post-IP preparation; unknown products stay explicit."""
    K,R,S=MaterialKind,InputRole,StepKind
    g=ProtocolGraph('iCLIP — replicate-1 source preparation',authority=SOURCE,representation='published_instance')
    def state(label):return Construct([Segment(label,'[unassembled occupied molecular state]',placeholder=True,inferred=True)],name=label)
    def condition(i,label):return g.input(i,label,K.CONDITION)
    def clean(i,x,label,carrier=None):
        c=condition(i+'_conditions',label)
        y=g.product(i+'_retained','Retained '+label,K.MOLECULAR,state=x.state)
        w=g.product(i+'_removed','Physically separated recovery rejects; composition unassigned',K.WASTE)
        ports=[(x,R.TEMPLATE),(c,R.CONDITION)]
        if carrier:ports.append((carrier,R.STRUCTURAL))
        g.react(i,label,S.CLEANUP,ports,(y,w));g.discard(w,'Source physical fractionation; no exact recovery/purity implied.')
        return y
    def split(i,x,label,kind=K.MOLECULAR):
        y=g.product(i+'_selected',label,kind,state=x.state)
        z=g.product(i+'_remaining','Untransferred/other actual source fraction; amount/future use unassembled',x.kind,state=x.state)
        g.react(i,label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.mark_terminal(Terminal(z.id,Disposition.UNRESOLVED,'Outside selected source contribution; not silently discarded.'))
        return y
    rna=g.input('bound_rna','Actual post-IP RNase-I-trimmed hnRNP-C/RNA on Dynabeads; upstream UV/IP outside selected graph',K.MOLECULAR,primary=True,state=state('protein-bound source RNA'))
    sap=g.input('sap','Shrimp alkaline phosphatase0.1µl;37C10min',K.CATALYST)
    oh=g.product('rna_3oh','Protein-associated RNA3OH target class',K.MOLECULAR,state=rna.state)
    g.react('dephosphorylate','Remove RNA3prime phosphate',S.DEPHOSPHORYLATION,((rna,R.TEMPLATE),(sap,R.CATALYST)),(oh,))
    oh=clean('preligation_washes',oh,'2high-salt/PNK/RNA-ligase-buffer bead washes')
    l3=g.input('l3',L3_WRITTEN,K.ADAPTER,state=l3_scene(),chemical_state={'5prime':'phosphate','3prime':'Puromycin','sugars':'mixed U/T notation retained; complete sugar partition unassigned'})
    ligase=g.input('rna_ligase','Source RNA ligase0.25µl;15µlL3mix,16Covernight',K.CATALYST)
    ligated=g.product('ligated_rna','Protein-bound RNA–L3 with adapter3Puromycin retained',K.MOLECULAR,state=state('RNA–L3 occupied product unassembled'))
    sc=l3_scene();sc.strand('RNA substrate',[Segment('source RNA','[RNA]',placeholder=True,inferred=True)],label='RNA3OH class',mod3='OH');sc.ligation_end('RNA substrate','source RNA',kind='single-strand',polarity='3-prime',owner='substrate')
    g.react('ligate_l3','Ligate L3 5P to RNA3OH',S.LIGATION,((oh,R.TEMPLATE),(l3,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.known('single-strand',sc,statement='Published L3 donor5P/3Puromycin and source RNA3OH class; full RNA sequence and sealed junction unassembled.'))
    release=condition('bead_release_buffer','Source5µlNuPAGEbuffer,70C5min; collect magnetic eluate')
    released=g.product('released_complex','Actual released protein–RNA–L3 complexes',K.MOLECULAR,state=ligated.state)
    beads=g.product('bead_release_residual','Residual original bead-associated material; recovery unquantified',K.MOLECULAR)
    g.react('release_complex','Collect heat/sample-buffer eluate from beads',S.EXTRACTION,((ligated,R.TEMPLATE),(release,R.CONDITION)),(released,beads));g.mark_terminal(Terminal(beads.id,Disposition.UNRESOLVED,'Source collects eluate; residual bead material unquantified.'))
    membrane=clean('gel_membrane',released,'NuPAGE→nitrocellulose/film→actual two selected membrane bands')
    selected=split('membrane_fraction',membrane,'One actual selected membrane size fraction; other selected fraction retained separately')
    pk=g.input('proteinase_k','Source200µl2mg/mlPK55C30min then130µlPK/7Murea repeat; exact enzyme carryover unassigned',K.CATALYST)
    recovered=g.product('recovered_rna','Peptide-bearing RNA–L3 after two source PK treatments',K.MOLECULAR,state=state('peptide-bearing RNA–L3'))
    g.react('recover_rna','ProteinaseK then PK/urea on selected membrane fragment',S.DEPROTEINIZATION,((selected,R.TEMPLATE),(pk,R.CATALYST)),(recovered,))
    aqueous=clean('phenol_extract',recovered,'Source phenol/CHCl3 phase extraction; aqueous fraction recovered')
    aq450=split('aqueous450',aqueous,'Source450µl aqueous transfer before recentrifugation')
    aq400=split('aqueous400',aq450,'Source400µl supernatant transfer after recentrifugation')
    blue=g.input('rna_glycoblue','Source0.5µlGlycoBlue molecular carrier',K.MOLECULAR)
    recovered=clean('rna_precipitation',aq400,'Source acetate/ethanol precipitation and80percent wash;12µlresuspension',blue)
    primer=g.input('rt_primer','Replicate-1 5P RTprimer;1µl2pmol/µl',K.OLIGO,state=Construct(list(linear_cdna())[:4]),chemical_state={'5prime':'phosphate'})
    nt=g.input('rt_nucleotides','Published literal1 ml10mM dNTP mix; unit applicability unresolved, no silentµlrepair',K.MOLECULAR)
    annealed=g.product('annealed_rna','RNA template component in annealed same-tube mixture',K.MOLECULAR,state=recovered.state)
    retained=g.product('retained_rt_primer','Same-tube retained RTprimer component, not new free addition',K.OLIGO,state=primer.state)
    retained_nt=g.product('retained_rt_nucleotides','Same-tube nucleotide component; no second addition',K.MOLECULAR)
    g.react('anneal_rt','RNA+RTprimer+dNTP70C5min then42C',S.HYBRIDIZATION,((recovered,R.TEMPLATE),(primer,R.STRUCTURAL),(nt,R.STRUCTURAL)),(annealed,retained,retained_nt))
    rt=g.input('superscript_iii','Source RTmix:5µl5xbuffer+1µlDTT+0.5µlSSIII+0.5µlRNasin; sourcecalls6µl although components sum7µl',K.CATALYST)
    cdna=g.product('rt_hybrid','INFERRED — source primer-derived cDNA plus RNA template/RT mixture',K.MOLECULAR,state=linear_cdna(),inferred=True)
    g.react('reverse_transcribe','INFERRED — SSIII42/50/55C; occupied adapter landing/copy derived, not established product yield',S.REVERSE_TRANSCRIPTION,((annealed,R.TEMPLATE),(retained,R.STRUCTURAL),(retained_nt,R.STRUCTURAL),(rt,R.CATALYST)),(cdna,))
    rnase=g.input('rnase_a','Source RNaseA1µl after95C2min;37C20min',K.CATALYST)
    digested=g.product('rnase_mixture','cDNA plus RNA digestion products; fragments not yet removed',K.MOLECULAR,state=linear_cdna(),inferred=True)
    g.react('remove_rna','Heat/RNaseA molecular digestion before physical recovery',S.DIGESTION,((cdna,R.TEMPLATE),(rnase,R.CATALYST)),(digested,))
    blue2=g.input('cdna_glycoblue','Source0.5µlGlycoBlue for cDNA precipitation',K.MOLECULAR)
    recovered=clean('cdna_precipitation',digested,'TE/acetate/ethanol precipitation;6µlresuspension',blue2)
    gel=clean('cdna_size_gel',recovered,'6percentTBEurea gel;100–175nt and175–350nt selected classes')
    sized=split('cdna_size_fraction',gel,'One actual selected cDNA size class; other selected class unexpanded')
    glycogen=g.input('gel_glycogen','Source0.5µlglycogen recovery carrier',K.MOLECULAR)
    cdna=clean('gel_recovery',sized,'Published literal400 mlTE gel crush/SpinX/ethanol; unit unresolved;12µlresuspension',glycogen)
    lig2=g.input('circligase','Source0.75µlCircLigaseII+MnCl2;60C1h',K.CATALYST)
    closed=circularize_ssdna(linear_cdna(),five_prime_phosphate=True)
    circle=g.product('cdna_circle','INFERRED — replicate-1 cDNA circle; complete primer retained',K.MOLECULAR,state=tuple(closed.rows()),chemical_state=closed,inferred=True)
    g.react('circularize','INFERRED — intramolecular source cDNA closure',S.CIRCULARIZATION,((cdna,R.TEMPLATE),(lig2,R.CATALYST)),(circle,),geometry=CircularizationGeometry.known('single-strand intramolecular ligation',statement='Published5Pprimer/CircLigaseII; unassigned captured RNA and efficiency.'),relations=(MolecularRelation.CIRCLE_CLOSURE,))
    opener=g.input('linearization_primer','Source BamHI-site annealing oligo;1µl10µM',K.OLIGO,state=Construct([Segment('linearization primer',LINEARIZATION_PRIMER)]))
    exposed=g.product('exposed_circle','Circle with annealed BamHI-site oligo; full occupied duplex geometry unassembled',K.MOLECULAR,state=circle.state)
    g.react('anneal_bamhi','95C then70steps reducing1C per cycle; expose source BamHI site',S.HYBRIDIZATION,((circle,R.TEMPLATE),(opener,R.STRUCTURAL)),(exposed,))
    bamhi=g.input('bamhi','Source3µlFastDigestBamHI37C30min',K.CATALYST)
    opened=g.product('opened_cdna','Source BamHI-linearized cDNA class; exact whole duplex/cut product unassembled',K.MOLECULAR,state=state('source opened cDNA'))
    g.react('linearize','Source BamHI reopening; no exact full-product claim',S.DIGESTION,((exposed,R.TEMPLATE),(bamhi,R.CATALYST)),(opened,))
    blue3=g.input('opened_glycoblue','Source0.5µlGlycoBlue molecular carrier',K.MOLECULAR)
    opened=clean('opened_precipitation',opened,'Source acetate/ethanol recovery;9µlresuspension',blue3)
    supply=g.input('pcr_primer_mix','Actual source0.3µlpaired-end primer mix;10µM each',K.MOLECULAR)
    f=g.product('pcr_forward','Known source PCR primerP5 component',K.OLIGO,state=Construct([Segment('P5PCR',PCR_PRIMERS[1])]))
    r=g.product('pcr_reverse','Known source PCR primerP7 component',K.OLIGO,state=Construct([Segment('P7PCR',PCR_PRIMERS[0])]))
    g.react('pcr_primer_components','Inventory of one actual primer-mix addition',S.MAGIC,((supply,R.TEMPLATE),),(f,r),changes_molecule=False,reason='Source explicitly supplies one two-primer mix; components not new separate doses.')
    mix=g.input('immomix','Source10µl2xImmomix prepared reaction stock',K.MOLECULAR)
    enz=g.product('pcr_enzyme','Copresent Immomix polymerizing catalyst class',K.CATALYST)
    nts=g.product('pcr_nucleotides','Copresent Immomix nucleotide class; composition unassigned',K.MOLECULAR)
    buf=g.product('pcr_buffer','Copresent Immomix buffer class',K.CONDITION)
    g.react('immomix_components','Inventory decomposition of one actual Immomix addition',S.MAGIC,((mix,R.TEMPLATE),),(enz,nts,buf),changes_molecule=False,reason='Actual source prepared mix; functional classes only, no invented identity or dose.')
    amp=g.product('pcr_product','Source PCR library class; whole occupied amplicon unassembled',K.MOLECULAR,state=state('actual source PCR product'))
    g.react('pcr','Source35cycles PCR; known primers, full binding/product reference unassembled',S.PCR,((opened,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(nts,R.STRUCTURAL),(enz,R.CATALYST),(buf,R.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(opened.state,reason='Printed primer polymers preserved; whole occupied RNA-derived BamHI-opened template and exact primer sites unassembled, not claimed undisclosed.'))
    library=clean('g25_desalt',amp,'SourceMicrospinG25 desalting')
    loaded=split('loading',library,'Selected actual sequencing contribution; pooled-replicate mixing/volumes unassigned',K.LIBRARY)
    g.sequence(loaded,SequencingPlan(limitation='SourceGA2 one lane/three replicates54nt reads. Calling PCR primers paired-end does not establish acquiredRead2; exact loaded primer/second-read acquisition/indexmask unassigned.',read_lengths={'Read 1':54},acquired_roles=('Read 1',)))
    g.declare_reviewed_scope(SOURCE_INVENTORY,SOURCE)
    return g.require_publishable()

SOURCE_INVENTORY=('bound_rna','sap','dephosphorylate','preligation_washes_conditions','preligation_washes','l3','rna_ligase','ligate_l3','bead_release_buffer','release_complex','gel_membrane_conditions','gel_membrane','membrane_fraction','proteinase_k','recover_rna','phenol_extract_conditions','phenol_extract','aqueous450','aqueous400','rna_glycoblue','rna_precipitation_conditions','rna_precipitation','rt_primer','rt_nucleotides','anneal_rt','superscript_iii','reverse_transcribe','rnase_a','remove_rna','cdna_glycoblue','cdna_precipitation_conditions','cdna_precipitation','cdna_size_gel_conditions','cdna_size_gel','cdna_size_fraction','gel_glycogen','gel_recovery_conditions','gel_recovery','circligase','circularize','linearization_primer','anneal_bamhi','bamhi','linearize','opened_glycoblue','opened_precipitation_conditions','opened_precipitation','pcr_primer_mix','pcr_primer_components','immomix','immomix_components','pcr','g25_desalt_conditions','g25_desalt','loading','loading_selected')


def render_page():
    g=protocol_graph()
    return assemble_page(g,(head('iCLIP'),'<div class="wrap" data-protocol-complete="true"><h1>iCLIP</h1>',info('König et al. 2010, <a href="https://doi.org/10.1038/nsmb.1838">Nature Structural &amp; Molecular Biology</a>. Publisher Supplementary Methods governs this selected replicate-1 post-IP preparation.'),
    '<p>Recover crosslink-truncated cDNA by circularizing a phosphorylated RT-primer-derived molecule, reopening its internal BamHI site and amplifying with source Illumina primers. Upstream HeLa UV/RNase-I/hnRNP-C capture is represented by its actual post-IP material boundary.</p>',
    info('INFERRED — the occupied RNA interval and source-derived cDNA reference are unassigned. Exact whole BamHI products, PCR binding sites and acquired sequencing-primer geometry remain unassembled; this is a source-supported preparation, not a complete base-by-base experimental construct.'),
    '<h2>Source oligos</h2>',panel(l3_scene().rows(),caption='Source L3: 5′phosphate, mixed U/T notation and 3′Puromycin retained. Complete sugar partition is unassigned.'),
    '<p>Replicate-1 RT primer: <seq>/5phos/'+RT_PRIMERS[0]+'</seq><br>Linearization primer: <seq>'+LINEARIZATION_PRIMER+'</seq><br>PCR P7: <seq>'+PCR_PRIMERS[0]+'</seq><br>PCR P5: <seq>'+PCR_PRIMERS[1]+'</seq></p>',
    '<h2>Read setup and identifiers</h2><p>The source reports Illumina GA2 sequencing with 54-nt reads. It uses a three-base random barcode and a two-base experiment identifier introduced in RT; this graph selects replicate 1 (CA). The source calls its PCR primers paired-end, which does not alone establish acquired Read 2. No separate index-read mask or exact sequencing-primer landing is invented.</p>',
    info('Printed source values retained: the RT paragraph says 1 ml of 10 mM dNTPs, the gel paragraph says 400 ml TE, and a mixture named 6 µl lists components totaling 7 µl. These are unresolved source quantities/conventions, not silently repaired executable instructions. Replicate-2 GA/shorter primer and the mapping identifier set TG/TC/CA remain unreconciled; no changed primer is substituted.'),
    '<details><summary>Full selected workflow and unresolved molecular checks</summary>',WorkflowSection('Post-IP RNA → source cDNA circle → BamHI → library'),'</details>',
    '<p><a href="01_iclip.html">Source notes and limitations</a></p></div>'))

def selftest():
    from checks import Check
    c = Check()
    # Regression for the exact, non-terminal RT landing interval in the printed L3.
    expected = revcomp(L3.replace('U', 'T'))
    c('printed RT1 lands exactly within L3',
      RT_PRIMERS[0][-9:] == expected[2:11], True)
    c('source-supported selected preparation passes normal gate', protocol_graph().publication_problems(), [])
    g=protocol_graph()
    ports=lambda i:tuple(p.material for p in g.steps[i].inputs)
    c('BamHI consumes actual annealed circle','exposed_circle' in ports('linearize'),True)
    c('PCR consumes recovered actual opened product','opened_precipitation_retained' in ports('pcr'),True)
    c('RT has same-tube retained primer','retained_rt_primer' in ports('reverse_transcribe'),True)
    c('RT does not re-add external primer','rt_primer' in ports('reverse_transcribe'),False)
    c('RNA fragments persist until source physical cleanup',g.materials['rnase_mixture'].kind,MaterialKind.MOLECULAR)
    c('source nucleotide contribution retained','retained_rt_nucleotides' in ports('reverse_transcribe'),True)
    c('source labelled Immomix is one actual stock',ports('immomix_components'),('immomix',))
    c('other selected membrane fraction unresolved',g.terminals['membrane_fraction_remaining'].disposition,Disposition.UNRESOLVED)
    c('loading remainder unresolved',g.terminals['loading_remaining'].disposition,Disposition.UNRESOLVED)
    c('no exact whole PCR reference invented',g.steps['pcr'].geometry.availability,'unavailable')
    from checks import Source,have
    from storage import sources_root
    import re
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'nsmb.1838_41594_2010_BFnsmb1838_MOESM21_ESM.pdf.txt','PublishedSupplementaryMethods15–17')
    if have(c,source,label='iCLIP publisher source oligos'):
        text=re.sub(r'\s+','',source.path.read_text())
        for name,seq in (('L3',L3),('RT1',RT_PRIMERS[0]),('linearization',LINEARIZATION_PRIMER),('PCRP7',PCR_PRIMERS[0]),('PCRP5',PCR_PRIMERS[1])):
            c('whole source '+name,seq in text,True)
    uncovered=protocol_graph();uncovered.input('unreviewed_extra','unreviewed source input',MaterialKind.CONDITION)
    c('new input not silently covered by curator inventory',any('source-inventory' in x and 'unreviewed_extra' in x for x in uncovered.publication_problems()),True)
    c.report()

if __name__ == '__main__':
    selftest()
