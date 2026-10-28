"""hydro-tRNAseq size enrichment, hydrolysis and source-specific end treatment."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Segment, feature
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE='Gogakos et al.2017 doi:10.1016/j.celrep.2017.07.029, accepted-manuscript Methods122–124'
PREDECESSOR='Hafner et al.2012 doi:10.1016/j.ymeth.2012.07.030, accepted-manuscript adapter table92–148'
BOUNDARY=('Barcoded small-RNA library is explicitly inherited from Hafner2012; actual adapter '
          'member, exact fragment-selection stage and acquired read geometry remain deferred. '
          'Gated defining supplements and VOR verification remain open. CIP treatment does '
          'not establish complete resolution of every possible cyclic/2′/3′ phosphate species.')

def barcoded_adapter_inventory():
    """Unassigned member of the cited twenty-adapter family, not fabricated Ad01 use."""
    return Construct([Segment('unassigned source pentamer barcode','N'*5,'cbc',placeholder=True,
                       feature=feature('hydro_sample','sample_index','whitelist',
                                       whitelist='Hafner2012 twenty adapters; actual member unassigned')),
                      Segment('common3′adapterDNA','TCGTATGCCGTCTTCTGCTTG')],
                     name='inherited-family inventory only;5′rApp/3′aminohexyl; actual member unassigned')

def graph_for_replicate(replicate):
    if replicate not in (1,2,3,4): raise ValueError("source replicates are1–4")
    duration=10 if replicate==1 else 60
    g=ProtocolGraph('hydro-tRNAseq — size enrichment and hydrolysis',authority=SOURCE,
                    representation='published_instance')
    total=g.input('total_rna','HEK293 total RNA after TRIzol extraction',MaterialKind.MOLECULAR,primary=True)
    gel=g.input('urea_page','12%urea-PAGE;60–100nt window',MaterialKind.SELECTION)
    enriched=g.product('enriched_rna','60–100nt RNA including mature and precursor tRNAs',MaterialKind.MOLECULAR)
    rejected=g.product('outside_window','RNA outside recovered60–100nt window',MaterialKind.WASTE)
    g.react('size_enrich','Resolve20µg totalRNA; recover60–100nt RNA',StepKind.SELECTION,
            ((total,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(enriched,rejected))
    g.discard(rejected,'Physical PAGE fraction outside the collected window')
    carbonate=g.input('carbonate','10mMNa2CO3+10mMNaHCO3,pH9.7,65°C',MaterialKind.CONDITION)
    mixture=g.product('hydrolysis_mixture','all hydrolysis fragments and any deacylation-derived molecular products',MaterialKind.MOLECULAR)
    g.react('hydrolyze',f'Limited alkaline hydrolysis:{duration}min, source replicate{replicate}',
            StepKind.FRAGMENTATION,((enriched,InputRole.TEMPLATE),(carbonate,InputRole.CONDITION)),(mixture,))
    cip=g.input('cip','calf intestinal phosphatase37°C1h',MaterialKind.CATALYST)
    treated=g.product('cip_treated','CIP-treated fragment mixture; exact residual cyclic/2′/3′phosphate distribution unassigned',MaterialKind.MOLECULAR)
    g.react('dephosphorylate','Source CIP dephosphorylation of hydrolysis products',
            StepKind.DEPHOSPHORYLATION,((mixture,InputRole.TEMPLATE),(cip,InputRole.CATALYST)),(treated,))
    pnk=g.input('pnk','T4PNK37°C1h',MaterialKind.CATALYST)
    atp=g.input('atp','1mMATP in source PNK buffer',MaterialKind.CONDITION)
    repaired=g.product('rephosphorylated','PNK-treated RNA-fragment mixture; no undeclared physical debris removal',MaterialKind.MOLECULAR)
    g.react('rephosphorylate','Rephosphorylate with T4PNK and ATP',StepKind.PHOSPHORYLATION,
            ((treated,InputRole.TEMPLATE),(pnk,InputRole.CATALYST),(atp,InputRole.CONDITION)),(repaired,))
    g.mark_terminal(Terminal(repaired.id,Disposition.UNRESOLVED,BOUNDARY))
    g.declare_reviewed_scope(('total_rna','urea_page','size_enrich','carbonate','hydrolyze',
                             'cip','dephosphorylate','pnk','atp','rephosphorylate'),SOURCE)
    return g

def protocol_graph():
    """Default source instance: replicate1,10min alkaline hydrolysis."""
    return graph_for_replicate(1)

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('molecular hydrolysis products remain before physical separation',
      g.materials['hydrolysis_mixture'].kind,MaterialKind.MOLECULAR)
    c('actual sample identifier remains unassigned',
      list(barcoded_adapter_inventory())[0].feature.role,'sample_index')
    c('partial graph does not assert complete sequencing library',g.publishable,False)
    c.report()

if __name__=='__main__': selftest()
