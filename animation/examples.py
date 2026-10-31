"""Explicit animation adapters, retaining protocol-owned data and uncertainty.

The reference fixture is illustrative sequence, never a published study molecule.
Protocol adapters choose individual shots; graph step names alone cannot infer events.
"""
from .model import Event, Movie, Shot, Strand, complement



def reference():
    template='AUGCAGUACCGUAGCUAGCUACGAUCGGAUCC'
    daughter=complement(template[::-1]); primer=8
    return Movie('Polymerase extension',[
        Shot('copy','Copy the template',[
            Strand('template',template,(1270,440),(330,440),'rna','RNA template','RNA'),
            Strand('daughter',daughter,(330,490),(1270,490),'primer','DNA',initially_visible=primer),
        ],[Event('extend','daughter','Reverse transcriptase',template='template',
                  template_start=len(template)-primer-1,primer_length=primer,exact=True,
                  note='dNTPs · primer 3′ OH')],
        'The primer provides a free three-prime hydroxyl. Reverse transcriptase follows the RNA template and extends the DNA in the five-prime to three-prime direction.',
        scope='Illustrative sequence · exact complementary copy',source='Animation reference fixture; synthetic sequence, not a published specimen.')
    ],provenance={'representation':'illustrative_example','sequence_origin':'Synthetic renderer/geometry fixture; no protocol source claim.'}).validate()


EXAMPLES={'reference':reference}

# Imported only for explicit CLI examples; chemistry modules remain protocol-owned.
from . import protocols
EXAMPLES.update({'crispr-mip':protocols.crispr_mip,'boloramis':protocols.boloramis,
                 'atac':protocols.atac,'smallseq':protocols.smallseq,
                 'gibson':protocols.gibson,'hic':protocols.hic,'dropseq':protocols.dropseq,
                 'scrrbs':protocols.scrrbs,'celseq':protocols.celseq})
EXAMPLES['smartseq']=protocols.smartseq

# Public adapter API matches the CLI.
smartseq=protocols.smartseq
