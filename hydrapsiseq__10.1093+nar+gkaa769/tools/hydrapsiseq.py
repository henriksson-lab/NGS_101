"""HydraPsiSeq source-scoped cleavage and end-repair chemistry, before the kit."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE='Marchand et al. 2020 doi:10.1093/nar/gkaa769, published Methods and Supplementary Figure S1'
BOUNDARY=('Historical NEBNext Small RNA Library kit identity/revision, its exact oligos and '
          'downstream molecular/readout geometry remain unreviewed; no modern kit is substituted.')

@dataclass(frozen=True)
class CleavageEvent:
    """One source-described U scission event, not a deterministic population predictor."""
    upstream_three_prime: str='phosphate'
    downstream_five_prime: str='phosphate'
    removed_residue: str='U'
    debris: str='released residue-derived molecular products present'
    def __post_init__(self):
        if self.removed_residue!='U':
            raise ValueError('only a selected U-cleavage event is modelled; Ψ is relatively protected')
        if self.upstream_three_prime not in {'phosphate','OH'} or self.downstream_five_prime!='phosphate':
            raise ValueError('event must retain the source-specified downstream 5-prime phosphate')
    def dephosphorylate(self):
        return replace(self,upstream_three_prime='OH')
    def scene(self):
        sc=Scene()
        sc.strand('upstream',[Segment('RNA through N−1','[upstream RNA]',placeholder=True)],
                  label='N−1 fragment',mod3=self.upstream_three_prime)
        sc.strand('downstream',[Segment('RNA from N+1','[downstream RNA]',placeholder=True)],
                  label='N+1 fragment',mod5=self.downstream_five_prime)
        sc.note('downstream','The excised U is absent from both RNA fragments; '+self.debris+'.')
        return sc

def protocol_graph():
    g=ProtocolGraph('HydraPsiSeq — hydrazine/aniline cleavage and end repair',authority=SOURCE,
                    representation='published_instance')
    rna=g.input('rna','RNA population: U sites cleave stochastically; Ψ relatively protected',
                MaterialKind.MOLECULAR,primary=True)
    hydrazine=g.input('hydrazine','hydrazine, final50%, on ice30–60min',MaterialKind.CONDITION)
    abasic=g.product('abasic_mixture','RNA abasic sites plus released base-derived molecular products',MaterialKind.MOLECULAR)
    g.react('hydrazinolysis','Hydrazine ring opening/abasic-site formation at susceptible U sites',
            StepKind.CHEMICAL_CONVERSION,((rna,InputRole.TEMPLATE),(hydrazine,InputRole.CONDITION)),(abasic,))
    precip=g.product('hydrazine_pellet','RNA-containing ethanol precipitate; exact debris partition unquantified',MaterialKind.MOLECULAR)
    waste=g.product('hydrazine_supernatant','ethanol supernatant and washes',MaterialKind.WASTE)
    g.react('precipitate_hydrazine','NaOAc/Glycoblue ethanol precipitation and washes',StepKind.CLEANUP,
            ((abasic,InputRole.TEMPLATE),),(precip,waste))
    g.discard(waste,'Physical supernatant/wash separation; no exact composition or quantitative recovery claimed')
    aniline=g.input('aniline','1M aniline,pH4.5,60°C15min in dark',MaterialKind.CONDITION)
    event=CleavageEvent(debris='released residue-derived products; prior hydrazine debris partition unquantified')
    cleaved=g.product('cleavage_mixture','N−1 RNA3′p + N+1 RNA5′p + released nonpolymeric products',
                      MaterialKind.MOLECULAR,state=tuple(event.scene().rows()),chemical_state=event)
    g.react('aniline_scission','Aniline decomposes abasic sites by β/δ elimination',
            StepKind.CHEMICAL_CONVERSION,((precip,InputRole.TEMPLATE),(aniline,InputRole.CONDITION)),(cleaved,))
    cleaned=replace(event,debris='partition of released products after precipitation unresolved')
    fragments=g.product('recovered_fragments','ethanol-recovered RNA-fragment population',
                        MaterialKind.MOLECULAR,state=tuple(cleaned.scene().rows()),chemical_state=cleaned)
    rejected=g.product('aniline_supernatant','aniline-stage precipitation supernatant/washes',MaterialKind.WASTE)
    g.react('precipitate_aniline','Ethanol precipitation after aniline; retain RNA fragments',StepKind.CLEANUP,
            ((cleaved,InputRole.TEMPLATE),),(fragments,rejected))
    g.discard(rejected,'Physical supernatant/wash fraction; released products were molecular material before separation')
    pnk=g.input('pnk','T4 PNK3′phosphatase, ATP-free,37°C6h',MaterialKind.CATALYST)
    repaired=cleaned.dephosphorylate()
    ends=g.product('repaired_fragments','RNA fragment3′OH; downstream5′p retained',MaterialKind.MOLECULAR,
                   state=tuple(repaired.scene().rows()),chemical_state=repaired)
    g.react('dephosphorylate','ATP-free T4PNK3′dephosphorylation, then heat inactivation',
            StepKind.DEPHOSPHORYLATION,((fragments,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),(ends,))
    extract=g.product('kit_input','phenol/chloroform-extracted and ethanol-recovered kit input RNA',MaterialKind.MOLECULAR,
                      state=ends.state,chemical_state=repaired)
    extraction_waste=g.product('extraction_waste','extraction phases and precipitation supernatant',MaterialKind.WASTE)
    g.react('recover_ends','Phenol/chloroform and ethanol recovery after PNK',StepKind.CLEANUP,
            ((ends,InputRole.TEMPLATE),),(extract,extraction_waste))
    g.discard(extraction_waste,'Physical extraction and precipitation reject')
    g.mark_terminal(Terminal(extract.id,Disposition.UNRESOLVED,BOUNDARY))
    g.declare_reviewed_scope(('rna','hydrazine','hydrazinolysis','precipitate_hydrazine',
                             'aniline','aniline_scission','precipitate_aniline','pnk',
                             'dephosphorylate','recover_ends'),SOURCE)
    return g

def selftest():
    from checks import Check
    c=Check();event=CleavageEvent();repaired=event.dephosphorylate()
    c('ATP-free repair preserves5p while removing3p',
      (repaired.upstream_three_prime,repaired.downstream_five_prime),('OH','phosphate'))
    try: CleavageEvent(removed_residue='Ψ')
    except ValueError: rejected=True
    else: rejected=False
    c('U event is not silently reused as a Ψ cleavage prediction',rejected,True)
    c('unreviewed kit prevents sequencing publication',protocol_graph().publishable,False)
    c.report()

if __name__=='__main__': selftest()
