#!/usr/bin/env python3
"""Original library-order and nonphysical-fate regression checks."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import chem as M
from checks import Check, Source, have
from protocol_graph import Disposition, MaterialKind
c=Check()
root=data_root()/"sources"/H.parent.name
g=M.protocol_graph()
source_order=("5' polyphosphatase","5' terminator exonuclease","3' small RNA linker","purification on a gel",
              "alkaline phosphatase","tobacco acid pyrophosphatase","5' small RNA linker","reverse transcription","library amplification")
model_order=('treat_polyphosphate','digest_uncapped','ligate3','purify_linked3','dephosphorylate','decap','ligate5','reverse_transcribe','library_pcr')
for prefix,accession in (('five','GSM463298'),('three','GSM463301')):
    source=Source(root/('geo-'+accession+'-self-full.txt'),'Registered original GEO sample; tools/get_sources.py --slugs '+H.parent.name)
    if have(c,source,label=accession+' original extraction method'):
        text=source.path.read_text()
        protocol=next(x.split(' = ',1)[1] for x in text.splitlines() if x.startswith('!Sample_extract_protocol_ch1'))
        positions=[protocol.index(token) for token in source_order]
        c(accession+' original enzyme/linker/gel order',positions,sorted(positions))
        ids=list(g.steps)
        positions=[ids.index(prefix+'_'+name) for name in model_order]
        c(accession+' model retains all source processing stages in order',positions,sorted(positions))
        c(accession+' original RNA interval retained',f'{M.RNA_RANGE[0]}-{M.RNA_RANGE[1]}nt' in protocol)
        processing=next(x for x in text.splitlines() if x.startswith('!Sample_data_processing'))
        c(accession+' 26 nt is explicitly a trimming rule','trimmed to 26nt' in processing)
        c(accession+' does not turn trimmed 26 nt into acquired cycle count',M.READ_LENGTHS.get('Read 1') is None)
    degraded=prefix+'_digested_uncapped'
    c(accession+' nuclease products have no asserted immediate physical discard',g.terminals[degraded].disposition,Disposition.UNRESOLVED)
    c(accession+' nuclease products remain a molecular product class',g.materials[degraded].kind,MaterialKind.MOLECULAR)
c.report()
