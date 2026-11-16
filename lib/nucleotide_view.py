"""Checked linear sequence maps of selected physical simulator strands.

Self-pairing topology remains in the snapshot; a linear map does not depict
its conformation. Use topology.FoldLayout for the corresponding folded view.
"""
from dataclasses import dataclass
from chemdraw import Scene, Segment


@dataclass
class NucleotideState:
    snapshot: object
    strands: tuple  # (label, tuple of physical nucleotide identities)
    note: str = ''

    def rows(self):
        sc=Scene();placed={}
        for j,(label,ids) in enumerate(self.strands):
            ids=tuple(ids)
            if not ids or len(set(ids))!=len(ids):
                raise ValueError('sequence map requires one occurrence per physical nucleotide')
            for a,b in zip(ids,ids[1:]):
                if self.snapshot.next.get(a)!=b:
                    raise ValueError('sequence map cannot join disconnected backbones')
            parts=[Segment(i,self.snapshot.atoms[i].base,
                           'r1' if 'adapter' in self.snapshot.atoms[i].features else None,
                           inferred='synthetic' in self.snapshot.atoms[i].features) for i in ids]
            name='chain'+str(j)
            end=dict(mod5=self.snapshot.five.get(ids[0],''),mod3=self.snapshot.three.get(ids[-1],''))
            anchor=next(((other,i,self.snapshot.paired[i]) for other,prior in placed.items()
                         for i in ids if self.snapshot.paired.get(i) in prior),None)
            if anchor:
                other,i,mate=anchor
                unpaired=tuple(i for i in ids if self.snapshot.paired.get(i) not in placed[other])
                sc.anneal(name,parts,to=other,pair=(i,mate),label=label,unpaired=unpaired,**end)
            else:
                sc.strand(name,parts,label=label,**end)
            placed[name]=set(ids)
        if self.note:sc.note(name,self.note)
        return sc.rows()
