"""Fixed plans and shared-response comparisons, without fitting new parameters."""
import numpy as np

PHASES = (0., .5, 1., 1.5)
ENERGIES = (1.625, 2.5, 5.486, 8.)
CONDITIONS = (('high_weak',8.), ('high_weak',40.), ('high_weak',200.),
              ('high_weak',1200.), ('weak_high',40.))


def make_plans():
    native = []
    for role, count in [('pilot',2), ('formal',32)]:
        for group, energy in enumerate(ENERGIES):
            for block in range(count):
                i = len(native)
                native.append(dict(id=i,role=role,group=group,block=block,source_energy_MeV=energy,
                                   seed=740000000+i,sipm_seed=740100000+i,readout_seed=740200000+i))
    single = [dict(id=r['id']*4+p,kind='single',parents=[r['id']],times_ns=[128.+phase],
                   phase_ns=phase,group=r['group'],block=r['block'],role=r['role'],
                   seed_base=740300000+r['id']*100,repeats=8)
              for r in native for p,phase in enumerate(PHASES)]
    pairs=[]
    for block in range(32):
        weak,high=8+block,104+block
        for condition,gap in CONDITIONS:
            for phase in PHASES:
                pairs.append(dict(id=len(pairs),kind='pair',block=block,condition=condition,gap_ns=gap,
                    phase_ns=phase,parents=[high,weak] if condition=='high_weak' else [weak,high],
                    times_ns=[128.+phase,128.+phase+gap],seed_base=740500000+block*1000,repeats=8))
    return native,single,pairs


def anchored_pair_predictions(old, corrected, times, energies):
    """Compute the expensive shared term once; require identical isolated Q(E).

    The inherited mean correction changes only isolated voltage. This algebraic
    shortcut is additionally compared to direct corrected-model components on
    actual pilot inputs before any formal shared-response data are produced.
    """
    solo_old,q_old=old.isolated(times,energies)
    solo_new,q_new=corrected.isolated(times,energies)
    if not np.array_equal(q_old,q_new):
        raise ValueError('Shared-term reuse requires an unchanged isolated charge map')
    joint_old=old.components(times,energies)
    joint_new=joint_old+(solo_new-solo_old)
    return dict(old_solo=solo_old,corrected_solo=solo_new,old_joint=joint_old,corrected_joint=joint_new)


def paired_mean_ci(values,index):
    values=np.asarray(values,float)
    index=np.asarray(index)
    if values.ndim!=1 or len(values)<2 or index.ndim!=2 or index.shape[1]!=len(values):
        raise ValueError('Bootstrap rows must resample whole mother-pair observations')
    if not np.isfinite(values).all():raise ValueError('Finite block statistics required')
    return np.percentile(values[index].mean(1),[2.5,97.5]).tolist()
