"""Compare saved transport events across software environments, without retuning."""
from __future__ import annotations

import numpy as np


def relocate_macro(body, output, configuration):
    replacements={'/gagg/run/outputFile':str(output),'/gagg/sipm/configFile':str(configuration)}
    found={k:0 for k in replacements};lines=[]
    for line in body.splitlines():
        command=line.strip().split(maxsplit=1)[0] if line.strip() else ''
        if command in replacements:
            found[command]+=1;line=command+' '+replacements[command]
        lines.append(line)
    if any(n!=1 for n in found.values()):raise ValueError('Exactly one output and configuration command required')
    return '\n'.join(lines)+'\n'


def macro_semantics(body):
    return '\n'.join(line.strip() for line in body.splitlines() if line.strip() and
                     not line.lstrip().startswith('#') and line.strip().split()[0] not in
                     ('/gagg/run/outputFile','/gagg/sipm/configFile'))


def photon_records(arrays):
    keys=('time_ps','x_mm','y_mm','energy_eV')
    lengths=[len(arrays[k]) for k in keys]
    if len(set(lengths))!=1:raise ValueError('Photon column lengths differ')
    records=np.empty(lengths[0],dtype=[('t','<i8'),('x','<u4'),('y','<u4'),('e','<u4')])
    records['t']=arrays['time_ps']
    for source,target in [('x_mm','x'),('y_mm','y'),('energy_eV','e')]:
        records[target]=np.asarray(arrays[source],np.float32).view(np.uint32)
    return records


def compare_photons(old,new):
    a,b=photon_records(old),photon_records(new)
    shape=len(a)==len(b)
    ordered=shape and np.array_equal(a,b)
    multiset=shape and np.array_equal(np.sort(a,order=['t','x','y','e']),np.sort(b,order=['t','x','y','e']))
    return dict(old_count=len(a),new_count=len(b),same_count=shape,ordered_equal=bool(ordered),
                sorted_exact_equal=bool(multiset),relative_count_change_percent=100*(len(b)/len(a)-1))


def event_statistics(photons, edep, source_time_ps=5000000):
    active=(abs(photons['x_mm'])<=3)&(abs(photons['y_mm'])<=3)
    times=(np.asarray(photons['time_ps'],float)-source_time_ps)/1000
    return dict(edep_MeV=float(edep),photons=len(times),active_photons=int(active.sum()),
        first_active_ns=float(times[active].min()),mean_arrival_ns=float(times.mean()),
        late300_fraction=float(np.mean(times>300)),late1200_fraction=float(np.mean(times>1200)))
