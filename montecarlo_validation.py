"""
montecarlo_validation.py
========================

Monte Carlo validation of the Keplerian collision-rate engine at reduced N.

Two independent checks against the analytic rate of collision_models.keplerian_model:

  A. Cube estimator (Liou et al., LEGEND): propagate the population with exact
     Keplerian dynamics, sample random epochs, partition space into cubes of
     side h, and for every pair found in the same cube accumulate
     sigma * v_rel / h^3. This converges to the kinetic-theory rate computed
     with the TRUE local densities and velocities, so it validates the
     analytic latitude integrals without assuming their form.

  B. Direct conjunction counting: exact two-body propagation of the same kind
     of population, detection of every close approach below an inflated
     capture radius R_c (sigma_c = pi R_c^2), using a radial prefilter
     (circular orbits: minimum possible separation >= |r1 - r2|) and a
     linearized minimum-distance refinement inside coarse time steps.
     This makes NO kinetic assumption at all: it is the ground truth the
     kinetic form must reproduce.

Both use a population drawn from the same element distributions as the
analytic model: inclination mix with uniform dispersion, RAAN ~ U(0, 2pi),
phase ~ U(0, 2pi), radius from the chosen radial profile, zero eccentricity.

Run:  python3 montecarlo_validation.py [cube|direct|all]
Results are appended to mc_results.json.
"""

from __future__ import annotations
import json
import sys
import time
import numpy as np

import collision_models as cm

RESULTS_FILE = "mc_results.json"


# ----------------------------------------------------------------------
# Population sampling and exact circular-orbit propagation
# ----------------------------------------------------------------------
def sample_population(inp: cm.Inputs, N: int, rng):
    """Draw orbital elements consistent with the analytic model assumptions."""
    d = cm.derive(inp)
    mix = np.array(inp.inclination_mix, dtype=float)
    w = mix[:, 1] / mix[:, 1].sum()
    pop_idx = rng.choice(len(w), size=N, p=w)
    di = np.deg2rad(inp.incl_dispersion_deg)
    inc = np.deg2rad(mix[pop_idx, 0]) + rng.uniform(-di, di, N)
    raan = rng.uniform(0.0, 2.0 * np.pi, N)
    u0 = rng.uniform(0.0, 2.0 * np.pi, N)
    if inp.radial_profile == "uniform_volume":
        # pdf ~ r^2 on [R_int, R_est]
        cdf = rng.uniform(0.0, 1.0, N)
        r = (d.R_int**3 + cdf * (d.R_est**3 - d.R_int**3)) ** (1.0 / 3.0)
    else:
        r0 = cm.R_EARTH + inp.shell_alt_km * 1e3
        wsh = inp.shell_width_km * 1e3
        r = rng.uniform(r0 - wsh / 2, r0 + wsh / 2, N)
    return {"inc": inc, "raan": raan, "u0": u0, "r": r,
            "n_mot": np.sqrt(cm.MU_EARTH / r**3), "v": np.sqrt(cm.MU_EARTH / r)}


def propagate(pop, t):
    """Exact two-body circular propagation. Returns positions and velocities,
    arrays (N, 3), at time t (scalar)."""
    u = pop["u0"] + pop["n_mot"] * t
    cu, su = np.cos(u), np.sin(u)
    cO, sO = np.cos(pop["raan"]), np.sin(pop["raan"])
    ci, si = np.cos(pop["inc"]), np.sin(pop["inc"])
    r = pop["r"]
    pos = np.stack([
        r * (cO * cu - sO * su * ci),
        r * (sO * cu + cO * su * ci),
        r * (su * si),
    ], axis=1)
    vmag = pop["v"]
    vel = np.stack([
        vmag * (-cO * su - sO * cu * ci),
        vmag * (-sO * su + cO * cu * ci),
        vmag * (cu * si),
    ], axis=1)
    return pos, vel


# ----------------------------------------------------------------------
# A. Cube estimator
# ----------------------------------------------------------------------
def cube_estimate(inp: cm.Inputs, N=2000, n_epochs=400, h_km=100.0, seed=1):
    rng = np.random.default_rng(seed)
    pop = sample_population(inp, N, rng)
    d = cm.derive(inp)
    h = h_km * 1e3
    rate_samples = np.zeros(n_epochs)
    n_pairs_tot = 0
    # spread epochs over many days so RAAN/phase geometry decorrelates via motion
    epochs = rng.uniform(0.0, 30 * 86400.0, n_epochs)
    for k, t in enumerate(epochs):
        pos, vel = propagate(pop, t)
        keys = np.floor(pos / h).astype(np.int64)
        # exact collision-free cube key (coordinates are within +-2^20 cells)
        kh = ((keys[:, 0] + 2**20) * 2**42
              + (keys[:, 1] + 2**20) * 2**21
              + (keys[:, 2] + 2**20))
        order = np.argsort(kh)
        khs = kh[order]
        rate = 0.0
        start = 0
        for j in range(1, len(khs) + 1):
            if j == len(khs) or khs[j] != khs[start]:
                if j - start >= 2:
                    idx = order[start:j]
                    for a in range(len(idx)):
                        for b in range(a + 1, len(idx)):
                            vrel = np.linalg.norm(vel[idx[a]] - vel[idx[b]])
                            rate += d.sigma * vrel / h**3
                            n_pairs_tot += 1
                start = j
        rate_samples[k] = rate
    mean = rate_samples.mean()
    se = rate_samples.std(ddof=1) / np.sqrt(n_epochs)
    return {"method": "cube", "N": N, "n_epochs": n_epochs, "h_km": h_km,
            "rate_per_s": float(mean), "rate_se": float(se),
            "E_per_year": float(mean * cm.YEAR), "E_se": float(se * cm.YEAR),
            "pairs_sampled": int(n_pairs_tot)}


# ----------------------------------------------------------------------
# B. Direct conjunction counting
# ----------------------------------------------------------------------
def direct_count(inp: cm.Inputs, N=1000, T_days=3.0, R_c_km=5.0,
                 dt_coarse=10.0, seed=2, t_offset_days=0.0):
    rng = np.random.default_rng(seed)
    pop = sample_population(inp, N, rng)
    R_c = R_c_km * 1e3

    # exact prefilter: circular orbits cannot approach closer than |r1 - r2|
    ii, jj = np.triu_indices(N, k=1)
    ok = np.abs(pop["r"][ii] - pop["r"][jj]) <= R_c
    pi, pj = ii[ok], jj[ok]
    n_elig = len(pi)

    # coarse threshold: max approach speed ~ 2 v_orb; within a coarse step a
    # pair can close by at most v_max * dt; flag generously.
    v_max = 2.05 * pop["v"].max()
    flag_dist = v_max * dt_coarse + 5 * R_c

    n_steps = int(T_days * 86400.0 / dt_coarse)
    t0 = t_offset_days * 86400.0
    events = []          # (time, i, j, miss_distance)
    last_event = {}
    for s in range(n_steps):
        t = t0 + s * dt_coarse
        pos, vel = propagate(pop, t)
        dr = pos[pi] - pos[pj]
        d2 = np.einsum("ij,ij->i", dr, dr)
        cand = np.where(d2 < flag_dist**2)[0]
        if len(cand) == 0:
            continue
        dv = vel[pi[cand]] - vel[pj[cand]]
        drc = dr[cand]
        dv2 = np.einsum("ij,ij->i", dv, dv)
        tstar = -np.einsum("ij,ij->i", drc, dv) / np.maximum(dv2, 1e-12)
        # minimum inside this coarse interval (centered convention: [0, dt))
        tstar_cl = np.clip(tstar, 0.0, dt_coarse)
        dmin2 = (np.einsum("ij,ij->i", drc, drc)
                 + 2.0 * tstar_cl * np.einsum("ij,ij->i", drc, dv)
                 + tstar_cl**2 * dv2)
        hits = np.where(dmin2 < R_c**2)[0]
        for hh in hits:
            a, b = int(pi[cand[hh]]), int(pj[cand[hh]])
            t_ev = t + float(tstar_cl[hh])
            key = (a, b)
            if key in last_event and t_ev - last_event[key] < 60.0:
                continue
            last_event[key] = t_ev
            events.append((t_ev, a, b, float(np.sqrt(dmin2[hh]))))
    return {"method": "direct", "N": N, "T_days": T_days, "R_c_km": R_c_km,
            "dt_coarse": dt_coarse, "n_eligible_pairs": int(n_elig),
            "n_events": len(events), "seed": seed,
            "events_per_year": len(events) * 365.25 / T_days,
            "event_times_days": [round(e[0] / 86400.0, 5) for e in events[:5000]],
            "miss_distances_km": [round(e[3] / 1e3, 3) for e in events[:5000]]}


def analytic_reference(inp: cm.Inputs, N, sigma_override=None):
    """Analytic Keplerian E[collisions/yr] for the reduced-N population."""
    kep = cm.keplerian_model(inp, n_quad=800)
    d = cm.derive(inp)
    E = kep["E_collisions"] * (N / inp.N) ** 2
    if sigma_override is not None:
        E *= sigma_override / d.sigma
    return E


def save(res):
    try:
        with open(RESULTS_FILE) as f:
            allres = json.load(f)
    except FileNotFoundError:
        allres = []
    allres.append(res)
    with open(RESULTS_FILE, "w") as f:
        json.dump(allres, f, indent=1)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    inp = cm.Inputs(n_mc_incl_pairs=300)
    if what in ("cube", "all"):
        t0 = time.time()
        res = cube_estimate(inp, N=2000, n_epochs=int(sys.argv[2]) if len(sys.argv) > 2 else 400)
        res["E_analytic"] = analytic_reference(inp, res["N"])
        res["runtime_s"] = round(time.time() - t0, 1)
        save(res)
        print(json.dumps({k: v for k, v in res.items() if k != "event_times_days"}, indent=1))
    if what in ("direct", "all"):
        t0 = time.time()
        days = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2
        off = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
        res = direct_count(inp, N=1000, T_days=days, R_c_km=5.0, seed=seed, t_offset_days=off)
        d = cm.derive(inp)
        res["E_analytic_per_year"] = analytic_reference(inp, res["N"], sigma_override=np.pi * (5e3) ** 2)
        res["runtime_s"] = round(time.time() - t0, 1)
        save(res)
        print(json.dumps({k: v for k, v in res.items()
                          if k not in ("event_times_days", "miss_distances_km")}, indent=1))
