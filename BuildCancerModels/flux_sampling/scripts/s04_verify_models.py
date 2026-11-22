"""s04_verify_models.py -- Reviewer-2 response, checkpoint before sampling.

Confirms that the 20 sprintcore models fetched from the shared Drive folder are
the same models the archived FVA sheets were produced from. fsranalysis_ccle_sprint.m
indexes the FVA columns by findRxnIDs(model, ...), which is only valid if the
xlsx rows are in model.rxns order -- so that ordering is checked too, not just
the reaction sets.
"""
import os
import sys

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xlsxio

ROOT = '/home/pavan/scratch/SPECTRA_review_2'
FVA = '/home/pavan/scratch/spectra/BuildCancerModels/fva'
ALL = ['ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295',
       'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961',
       'ACH-003150', 'ACH-003143', 'ACH-003155', 'ACH-000064', 'ACH-002466',
       'ACH-003161', 'ACH-001453', 'ACH-000642', 'ACH-000494', 'ACH-001854']


def load_model(ach):
    d = sio.loadmat(os.path.join(ROOT, 'models', 'sprintcore', ach + '.mat'),
                    squeeze_me=True, struct_as_record=False)
    return d['m']


def main():
    hdr = ('model', 'nRxn', 'nMet', 'nGene', 'fvaRows', 'orderMatch', 'setMatch')
    print('%-12s %7s %7s %7s %8s %11s %9s' % hdr)
    ok = True
    for i, a in enumerate(ALL):
        m = load_model(a)
        rx = [str(x) for x in np.atleast_1d(m.rxns)]
        mets = [str(x) for x in np.atleast_1d(m.mets)]
        genes = [str(x) for x in np.atleast_1d(m.genes)] if hasattr(m, 'genes') else []
        rows = xlsxio.read_sheets(os.path.join(FVA, a + '.xlsx'))['Sheet1']
        fr = [r[0] for r in rows[1:] if r and r[0] not in (None, '')]
        order_match = (rx == fr)
        set_match = (set(rx) == set(fr))
        ok &= order_match and set_match
        print('%-12s %7d %7d %7d %8d %11s %9s'
              % (a, len(rx), len(mets), len(genes), len(fr), order_match, set_match))
        if i == 0:
            fields = sorted(f for f in dir(m) if not f.startswith('_'))
            print('   model fields: %s' % ', '.join(fields))
            S = m.S
            print('   S: %s   lb range [%g, %g]   ub range [%g, %g]'
                  % (S.shape, np.min(m.lb), np.max(m.lb), np.min(m.ub), np.max(m.ub)))
            print('   distinct bound magnitudes: %s'
                  % np.unique(np.abs(np.concatenate([m.lb, m.ub])))[:10])

    print()
    if ok:
        print('PASS: every model matches its FVA sheet in both reaction set and order.')
    else:
        print('FAIL: at least one model does not line up with its FVA sheet. '
              'Do not proceed to sampling.')
        sys.exit(1)


if __name__ == '__main__':
    main()
