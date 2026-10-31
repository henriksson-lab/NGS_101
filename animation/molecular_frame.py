"""Replay generated reaction deltas at a timestamp, independent of coordinates.

Connectivity, base births/conversion and released nucleotide pools are molecular
state. The renderer is a projection of this state, never its chemical authority.
"""
from dataclasses import dataclass
from .model import progress


@dataclass
class MolecularFrame:
    atoms: dict
    backbone: set
    duplex: set
    released: set

    def validate(self):
        nxt={};prev={};mates={}
        for a,b in self.backbone:
            if a in nxt or b in prev or a not in self.atoms or b not in self.atoms or a in self.released or b in self.released:raise ValueError('invalid frame backbone')
            nxt[a]=b;prev[b]=a
        for a,b in self.duplex:
            if a in mates or b in mates or a not in self.atoms or b not in self.atoms or a in self.released or b in self.released:raise ValueError('invalid frame duplex')
            mates[a]=b;mates[b]=a
        return self


def molecular_frame(shot,seconds,wpm=125):
    frame=MolecularFrame({a['id']:a for a in shot.initial_atom_records},set(map(tuple,shot.initial_backbone)),{tuple(sorted(p))for p in shot.chemical_pairs},set(shot.initial_released))
    seen=set()
    for e in sorted(shot.events,key=lambda e:e.fraction_start):
        if not e.molecular_operation or seconds<shot.establish+shot.action*e.fraction_start:continue
        op=e.molecular_operation;p=progress(shot,e,seconds,wpm)
        # One Tn5 event changes both target backbones atomically. Multiple visual
        # cut/junction actors referencing it must not perform the reaction twice.
        key=op['operation_id'] if 'operation_id'in op else ('event',shot.events.index(e))
        if key in seen:continue
        seen.add(key)
        adds=op.get('added',[]);count=min(len(adds),int(len(adds)*p+1e-8))
        details=op.get('details',{})
        if details.get('incorporations'):
            for atom,copy in zip(adds[:count],details['incorporations'][:count]):
                frame.atoms[atom['id']]=atom
                for bond in op['backbone_added']:
                    if atom['id']in bond and all(i in frame.atoms for i in bond):frame.backbone.add(tuple(bond))
                old=copy['displaced']
                if old:frame.duplex.discard(tuple(sorted((old,copy['template']))))
                frame.duplex.add(tuple(sorted((atom['id'],copy['template']))))
        elif details.get('removed'):
            removed=details['removed'][:int(len(details['removed'])*p+1e-8)]
            for id in removed:
                frame.released.add(id)
                frame.backbone={edge for edge in frame.backbone if id not in edge}
                frame.duplex={edge for edge in frame.duplex if id not in edge}
        elif adds:
            for atom in adds[:count]:frame.atoms[atom['id']]=atom
            for bond in op.get('backbone_added',[]):
                if all(id in frame.atoms for id in bond):frame.backbone.add(tuple(bond))
        else:
            # Binding/unbinding and conversion can be shown progressively while
            # each individual bond/nucleotide has an unambiguous discrete state.
            if e.kind in {'pair','unpair'}:
                for pair in op.get('pairing_removed',[])[:int(len(op.get('pairing_removed',[]))*p+1e-8)]:frame.duplex.discard(tuple(sorted(pair)))
                for pair in op.get('pairing_added',[])[:int(len(op.get('pairing_added',[]))*p+1e-8)]:frame.duplex.add(tuple(sorted(pair)))
            elif op.get('updated'):
                for atom in op['updated'][:int(len(op['updated'])*p+1e-8)]:frame.atoms[atom['id']]=atom
            elif p>=.5:
                frame.backbone.difference_update(map(tuple,op.get('backbone_removed',[])))
                frame.backbone.update(map(tuple,op.get('backbone_added',[])))
                frame.duplex.difference_update(tuple(sorted(pair))for pair in op.get('pairing_removed',[]))
                frame.duplex.update(tuple(sorted(pair))for pair in op.get('pairing_added',[]))
    return frame.validate()
