#!/usr/bin/env python3
"""
Verify the sc3gf1d component conversion against pyrocko engine synthetics.

For random moment tensors and azimuths, synthetics computed by the pyrocko
engine from the store (Z up, R away from the source, T at azimuth + 90 deg)
are compared with synthetics built from the converted sc3gf1d components
with the Minson & Dreger (2008) formulas. Both must agree to numerical
precision.

Usage: verify_conversion.py <store_dir> [--n 20] [--seed 1]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from pyrocko import gf, moment_tensor as pmt

from convert_to_sc3gf1d import standard_components, store_columns


def mt_synthetics(c, mt, az):
    """Z, R, T from sc3gf1d components c for NED tensor mt and azimuth az [deg]."""
    mxx, myy, mzz, mxy, mxz, myz = mt
    a = np.radians(az)
    c2, s2, ca, sa = np.cos(2 * a), np.sin(2 * a), np.cos(a), np.sin(a)

    def zr(p):
        return (mxx * (0.5 * c[p + 'SS'] * c2 - c[p + 'DD'] / 6 + c[p + 'EP'] / 3)
                + myy * (-0.5 * c[p + 'SS'] * c2 - c[p + 'DD'] / 6 + c[p + 'EP'] / 3)
                + mzz * (c[p + 'DD'] / 3 + c[p + 'EP'] / 3)
                + mxy * c[p + 'SS'] * s2 + mxz * c[p + 'DS'] * ca + myz * c[p + 'DS'] * sa)

    t = (c['TSS'] * (0.5 * (mxx - myy) * s2 - mxy * c2)
         + c['TDS'] * (mxz * sa - myz * ca))
    return zr('Z'), zr('R'), t


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('store_dir')
    ap.add_argument('--n', type=int, default=20, help='number of random cases')
    ap.add_argument('--seed', type=int, default=1)
    args = ap.parse_args()

    store_dir = Path(args.store_dir).resolve()
    engine = gf.LocalEngine(store_dirs=[str(store_dir)])
    store_id = gf.Store(str(store_dir)).config.id
    store = engine.get_store(store_id)
    cfg = store.config
    if cfg.component_scheme != 'elastic10':
        sys.exit(f"store uses '{cfg.component_scheme}', only 'elastic10' is supported")

    rng = np.random.default_rng(args.seed)
    depths = cfg.coords[0]
    dists = cfg.coords[1]
    worst = 1.0
    for _ in range(args.n):
        dep = float(rng.choice(depths))
        dist = float(rng.choice(dists[1:] if len(dists) > 1 else dists))
        az = float(rng.uniform(0, 360))
        m6 = pmt.MomentTensor(strike=rng.uniform(0, 360), dip=rng.uniform(5, 90),
                              rake=rng.uniform(-180, 180), scalar_moment=1.0).m6()
        m6 = m6 + 0.3 * np.array([1, 1, 1, 0, 0, 0])   # add an isotropic part
        src = gf.MTSource(lat=0, lon=0, depth=dep, m6=m6)
        a = np.radians(az)
        north, east = dist * np.cos(a), dist * np.sin(a)
        targets = [gf.Target(lat=0, lon=0, north_shift=north, east_shift=east,
                             azimuth=azi, dip=dip, interpolation='nearest_neighbor',
                             store_id=store_id, quantity='displacement')
                   for azi, dip in ((0.0, -90.0), (az, 0.0), (az + 90.0, 0.0))]
        ref = engine.process(src, targets).pyrocko_traces()

        cols, itmin = store_columns(store, dep, dist, cfg.ncomponents)
        comps = standard_components(cols)
        t = (itmin + np.arange(len(cols[0]))) / cfg.sample_rate
        # m6 = (mnn, mee, mdd, mne, mnd, med) = (xx, yy, zz, xy, xz, yz) in NED
        syn = mt_synthetics(comps, m6, az)

        ccs = []
        for name, r, y in zip('ZRT', ref, syn):
            yi = np.interp(r.get_xdata(), t, y)
            d = r.ydata
            norm = np.linalg.norm(d) * np.linalg.norm(yi)
            cc = np.dot(d, yi) / norm if norm > 0 else 1.0
            ccs.append(cc)
            worst = min(worst, cc)
        print(f"depth {dep / 1e3:6.1f} km  dist {dist / 1e3:7.1f} km  az {az:5.1f}  "
              f"CC Z {ccs[0]:+.5f}  R {ccs[1]:+.5f}  T {ccs[2]:+.5f}")

    print(f"worst CC {worst:+.5f}")
    sys.exit(0 if worst > 0.999 else 1)


if __name__ == '__main__':
    main()
