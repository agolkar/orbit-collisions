"""
collision_models.py
===================

Collision-rate engine for LEO mega-constellations. Two models:

  1. Kinetic baseline ("Brownian"): Kessler (1978) particle-in-a-box.
     Satellites treated as molecules of a dilute gas in isotropic random
     motion inside a spherical shell. Rate = n * sigma * v_rel with a
     single scalar density and a single scalar relative speed.

  2. Keplerian kinetic model: same functional form rate = n * sigma * v_rel,
     but with the spatial density field n(r, beta) derived from the
     orbit-element distributions (Kessler 1981) and the relative-velocity
     structure derived from the admissible headings of circular orbits at
     each latitude. No Maxwellian assumption.

Conventions
-----------
* SI units internally (m, s, kg). Inputs in convenient units (km, deg, years).
* beta = geocentric latitude. i = inclination. A = heading angle measured
  from local East, A in [0, pi] (A > pi/2 for retrograde orbits).
* All rates are NATURAL (uncontrolled) rates unless multiplied by the
  avoidance failure fraction f_fail; residual = natural * f_fail.

Validity flags computed and reported:
* dilute regime: mean free path >> shell thickness, nu*T << 1 per object.
* thin-shell approximation quality.
* gravitational focusing negligible: v_esc(structure) << v_rel.

Author: developed with Claude for A. Golkar, TUM. June 2026.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
import numpy as np
from scipy import integrate

# ----------------------------------------------------------------------
# Physical constants
# ----------------------------------------------------------------------
MU_EARTH = 3.986004418e14        # m^3/s^2
R_EARTH = 6371.0e3               # m
YEAR = 3.156e7                   # s (value fixed in the problem statement)

# Starlink-like inclination mix: (inclination_deg, weight)
STARLINK_LIKE_MIX = ((43.0, 0.20), (53.0, 0.40), (70.0, 0.20), (97.6, 0.20))


# ----------------------------------------------------------------------
# Inputs (pure parameters) vs derived quantities: kept strictly separate.
# ----------------------------------------------------------------------
@dataclass
class Inputs:
    """All free parameters of both models. Nothing here is derived."""
    N: float = 80_000.0                  # number of satellites
    radiator_area_m2: float = 120.0      # A, one-sided radiator area
    shape_factor: float = 4.0            # sigma = shape_factor * A (4 = enveloping sphere)
    h_min_km: float = 500.0              # bottom of altitude band
    h_max_km: float = 800.0              # top of altitude band
    inclination_mix: tuple = STARLINK_LIKE_MIX   # ((deg, weight), ...)
    radial_profile: str = "uniform_volume"       # "uniform_volume" | "thin_shell"
    shell_alt_km: float = 650.0          # used if radial_profile == "thin_shell"
    shell_width_km: float = 10.0         # top-hat width, if thin_shell
    v_rel_iso_ms: float = 10.0e3         # isotropic-model relative speed (canonical 10 km/s)
    incl_dispersion_deg: float = 0.5     # half-width of uniform inclination spread per population.
                                         # Regularizes the log-divergent <n^2> at the latitude caps;
                                         # the RATE is insensitive to it, the F_spatial/F_velocity
                                         # decomposition depends on it logarithmically (see report).
    n_mc_incl_pairs: int = 800           # MC samples over (i, i') per population pair
    seed: int = 12345
    avoidance_failure: float = 1.0       # f_fail: fraction of encounters NOT avoided (1 = no CA)
    T_years: float = 1.0                 # exposure window
    # cascade-criterion parameters (parametric, see report)
    frag_per_collision: float = 1000.0   # lethal fragments per catastrophic collision
    frag_lifetime_years: float = 25.0    # mean fragment residence time in band
    frag_shape_factor: float = 1.0       # sigma_sat-frag = frag_shape_factor * A


@dataclass
class Derived:
    """Geometry and reference quantities computed once from Inputs."""
    R_int: float = 0.0       # m
    R_est: float = 0.0       # m
    R_bar: float = 0.0       # m
    V: float = 0.0           # m^3 (exact spherical shell volume)
    sigma: float = 0.0       # m^2
    n_mean: float = 0.0      # m^-3
    v_orb: float = 0.0       # m/s at R_bar
    T: float = 0.0           # s


def derive(inp: Inputs) -> Derived:
    d = Derived()
    d.R_int = R_EARTH + inp.h_min_km * 1e3
    d.R_est = R_EARTH + inp.h_max_km * 1e3
    d.R_bar = 0.5 * (d.R_int + d.R_est)
    d.V = (4.0 / 3.0) * np.pi * (d.R_est**3 - d.R_int**3)
    d.sigma = inp.shape_factor * inp.radiator_area_m2
    d.n_mean = inp.N / d.V
    d.v_orb = np.sqrt(MU_EARTH / d.R_bar)
    d.T = inp.T_years * YEAR
    return d


# ----------------------------------------------------------------------
# Model 1: kinetic baseline (particle in a box)
# ----------------------------------------------------------------------
def kinetic_model(inp: Inputs) -> dict:
    """Kessler 1978 particle-in-a-box. Returns all baseline outputs."""
    d = derive(inp)
    n, sig, vrel, T, N = d.n_mean, d.sigma, inp.v_rel_iso_ms, d.T, inp.N

    nu = n * sig * vrel                                  # s^-1, per satellite
    P1 = 1.0 - np.exp(-nu * T)                           # per-satellite probability
    E_total = 0.5 * N * nu * T                           # expected collisions, whole fleet
    mfp = 1.0 / (np.sqrt(2.0) * n * sig)                 # m, mean free path
    shell_thickness = d.R_est - d.R_int

    return {
        "model": "kinetic_baseline",
        "n_mean_m3": n,
        "sigma_m2": sig,
        "v_rel_ms": vrel,
        "nu_per_s": nu,
        "nu_per_year": nu * YEAR,
        "P1_per_window": P1,
        "E_collisions": E_total,
        "E_residual": E_total * inp.avoidance_failure,
        "mean_free_path_m": mfp,
        "mfp_over_shell_thickness": mfp / shell_thickness,
        "dilute_regime_ok": bool(mfp > 1e3 * shell_thickness and nu * T < 0.3),
        "V_m3": d.V,
        "v_orb_ms": d.v_orb,
    }


def required_shell_thickness_m(inp: Inputs, E_acc: float, exact: bool = True) -> float:
    """Invert the baseline: shell thickness Delta_r such that E_collisions <= E_acc.

    Thin-shell closed form: Delta_r = N^2 sigma v T / (8 pi R_bar^2 E_acc).
    With exact=True, solves the exact cubic V(R_out) = V_required keeping
    R_int fixed, since the thin-shell formula breaks for large Delta_r.
    """
    d = derive(inp)
    V_req = inp.N**2 * d.sigma * inp.v_rel_iso_ms * d.T / (2.0 * E_acc)
    if not exact:
        return V_req / (4.0 * np.pi * d.R_bar**2)
    R_out = (d.R_int**3 + 3.0 * V_req / (4.0 * np.pi)) ** (1.0 / 3.0)
    return R_out - d.R_int


# ----------------------------------------------------------------------
# Model 2: Keplerian kinetic model
# ----------------------------------------------------------------------
def i_effective(i_rad: float) -> float:
    """Maximum latitude reached by a circular orbit of inclination i."""
    return min(i_rad, np.pi - i_rad)


def latitude_pdf(beta, i_rad):
    """Time-fraction density p(beta) for one circular orbit, inclination i.
    p(beta) = (1/pi) cos(beta) / sqrt(sin^2 i_eff - sin^2 beta), |beta| < i_eff.
    Normalized: integral over (-i_eff, i_eff) = 1. Integrable singularity at
    the turning latitudes (the well-known density pile-up at beta = +-i)."""
    ie = i_effective(i_rad)
    s2 = np.sin(ie) ** 2 - np.sin(beta) ** 2
    out = np.where(s2 > 0, np.cos(beta) / (np.pi * np.sqrt(np.maximum(s2, 1e-300))), 0.0)
    return out


def heading_from_east(beta, i_rad):
    """Heading angle A in [0, pi] of the velocity, measured from local East,
    for a circular orbit of inclination i crossing latitude beta.
    cos A = cos i / cos beta ; the ascending pass has heading +A (northward
    component), the descending pass -A."""
    ie = i_effective(i_rad)
    s2 = np.sin(ie) ** 2 - np.sin(beta) ** 2
    sinA = np.sqrt(np.maximum(s2, 0.0)) / np.cos(beta)
    cosA = np.cos(i_rad) / np.cos(beta)
    return np.arctan2(sinA, cosA)


def pair_vrel_options(beta, i1_rad, i2_rad, v_orb):
    """The two admissible relative speeds at latitude beta between members of
    populations with inclinations i1, i2 (circular, same radius).
    Senses (ascending/descending) are equiprobable and independent, so each
    option carries probability 1/2:
       theta = |A1 - A2|  (same sense)   -> v_a
       theta =  A1 + A2   (opposite)     -> v_b
    v_rel = 2 v_orb sin(theta/2). Returns (v_a, v_b)."""
    A1 = heading_from_east(beta, i1_rad)
    A2 = heading_from_east(beta, i2_rad)
    v_a = 2.0 * v_orb * np.abs(np.sin(0.5 * (A1 - A2)))
    v_b = 2.0 * v_orb * np.sin(0.5 * (A1 + A2))
    return v_a, v_b


def _radial_factor(inp: Inputs, d: Derived) -> float:
    """I_r = integral g(r)^2 / (2 pi r^2) dr, where g is the radial pdf.
    uniform_volume: g = 3 r^2/(R2^3 - R1^3)  ->  I_r = 2/V (exact).
    thin_shell (top-hat width w at r0):       I_r ~ 1/(2 pi r0^2 w)."""
    if inp.radial_profile == "uniform_volume":
        return 2.0 / d.V
    if inp.radial_profile == "thin_shell":
        r0 = R_EARTH + inp.shell_alt_km * 1e3
        w = inp.shell_width_km * 1e3
        return 1.0 / (2.0 * np.pi * r0**2 * w)
    raise ValueError(inp.radial_profile)


_LEGGAUSS_CACHE: dict = {}


def _leggauss(n):
    if n not in _LEGGAUSS_CACHE:
        _LEGGAUSS_CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _LEGGAUSS_CACHE[n]


def _pair_latitude_integral(i1, i2, v_orb, with_velocity=True, n_quad=1200):
    """J_jk = integral over beta of p_1 p_2 <v_rel> / cos(beta),
    on (-b_max, b_max), b_max = min(i_eff1, i_eff2). If with_velocity is
    False, <v_rel> is replaced by 1 (used for the pure spatial factor).

    Substitution sin(beta) = sin(b_max) sin(phi) removes the endpoint
    singularity of the lower-inclination population; for equal effective
    inclinations the vanishing of v_rel at the caps keeps the integrand
    finite. Gauss-Legendre on phi in (0, pi/2), doubled for symmetry."""
    b_max = min(i_effective(i1), i_effective(i2))
    if b_max <= 0:
        return 0.0
    nodes, weights = _leggauss(n_quad)
    phi = 0.5 * (nodes + 1.0) * (np.pi / 2.0)
    wphi = weights * (np.pi / 4.0)
    sb = np.sin(b_max) * np.sin(phi)
    beta = np.arcsin(np.clip(sb, -1.0, 1.0))
    dbeta_dphi = np.sin(b_max) * np.cos(phi) / np.cos(beta)

    p1 = latitude_pdf(beta, i1)
    p2 = latitude_pdf(beta, i2)
    if with_velocity:
        v_a, v_b = pair_vrel_options(beta, i1, i2, v_orb)
        v_mean = 0.5 * (v_a + v_b)
    else:
        v_mean = 1.0
    integrand = p1 * p2 * v_mean / np.cos(beta) * dbeta_dphi
    integrand = np.where(np.isfinite(integrand), integrand, 0.0)
    return 2.0 * np.sum(wphi * integrand)   # factor 2: beta symmetry


def keplerian_model(inp: Inputs, n_quad=1200) -> dict:
    """Keplerian kinetic model. Rate
       R = 1/2 sigma I_r sum_jk N_j N_k J_jk
    with J_jk the latitude integral above (v evaluated at R_bar; the
    variation of v_orb across the band is < 2 percent and is declared as an
    approximation). Also returns the velocity / spatial decomposition:
       R_kep / R_kin = F_spatial * F_velocity
       F_spatial  = [integral (sum_k n_k)^2 dV] / (n_mean^2 V)
       F_velocity = rate-weighted mean v_rel / v_rel_iso
    """
    d = derive(inp)
    mix = [(np.deg2rad(ideg), w) for ideg, w in inp.inclination_mix]
    wsum = sum(w for _, w in mix)
    mix = [(i, w / wsum) for i, w in mix]
    I_r = _radial_factor(inp, d)
    di = np.deg2rad(inp.incl_dispersion_deg)
    rng = np.random.default_rng(inp.seed)

    n_pop = len(mix)
    J_rate = np.zeros((n_pop, n_pop))
    J_geom = np.zeros((n_pop, n_pop))
    for a, (i1, w1) in enumerate(mix):
        for b, (i2, w2) in enumerate(mix):
            if di > 0:
                # average over the inclination dispersion (uniform +-di);
                # regularizes the log-divergent same-inclination <n^2> term
                u1 = i1 + di * (2.0 * rng.random(inp.n_mc_incl_pairs) - 1.0)
                u2 = i2 + di * (2.0 * rng.random(inp.n_mc_incl_pairs) - 1.0)
                jr = jg = 0.0
                for x, y in zip(u1, u2):
                    jr += _pair_latitude_integral(x, y, d.v_orb, True, n_quad)
                    jg += _pair_latitude_integral(x, y, d.v_orb, False, n_quad)
                J_rate[a, b] = jr / inp.n_mc_incl_pairs
                J_geom[a, b] = jg / inp.n_mc_incl_pairs
            else:
                J_rate[a, b] = _pair_latitude_integral(i1, i2, d.v_orb, True, n_quad)
                J_geom[a, b] = _pair_latitude_integral(i1, i2, d.v_orb, False, n_quad)

    w = np.array([w for _, w in mix])
    Nk = inp.N * w
    S_rate = Nk @ J_rate @ Nk
    S_geom = Nk @ J_geom @ Nk
    # per-satellite rate for a member of population k
    nu_pop = d.sigma * I_r * (J_rate @ Nk)               # s^-1, by population

    R_kep = 0.5 * d.sigma * I_r * S_rate                 # collisions / s, whole fleet

    kin = kinetic_model(inp)
    R_kin = kin["E_collisions"] / d.T

    F_spatial = (I_r * S_geom) / (d.n_mean**2 * d.V)     # >= 1 by Cauchy-Schwarz for matched radial profile
    v_eff = S_rate / S_geom                              # rate-weighted mean relative speed
    F_velocity = v_eff / inp.v_rel_iso_ms

    nu = 2.0 * R_kep / inp.N                             # mean per-satellite rate
    E_total = R_kep * d.T
    return {
        "model": "keplerian",
        "E_collisions": E_total,
        "E_residual": E_total * inp.avoidance_failure,
        "nu_per_year_mean": nu * YEAR,
        "P1_per_window_mean": 1.0 - np.exp(-nu * d.T),
        "v_rel_effective_ms": v_eff,
        "F_spatial": F_spatial,
        "F_velocity": F_velocity,
        "ratio_kep_over_kin": E_total / kin["E_collisions"],
        "R_kin_check": R_kin * d.T,
        "v_orb_ms": d.v_orb,
        "nu_per_year_by_population": {f"{ideg:.1f}deg": float(nu * YEAR)
                                      for (ideg, _), nu in zip(inp.inclination_mix, nu_pop)},
    }


# ----------------------------------------------------------------------
# Distributions (for plots, the report and the future web app)
# ----------------------------------------------------------------------
def density_vs_latitude(inp: Inputs, beta_grid=None):
    """Number density (m^-3) vs latitude at r = R_bar, plus the uniform value."""
    d = derive(inp)
    if beta_grid is None:
        beta_grid = np.linspace(-np.pi / 2 * 0.999, np.pi / 2 * 0.999, 2001)
    mix = [(np.deg2rad(ideg), w) for ideg, w in inp.inclination_mix]
    wsum = sum(w for _, w in mix)
    if inp.radial_profile == "uniform_volume":
        g = 3.0 * d.R_bar**2 / (d.R_est**3 - d.R_int**3)
    else:
        g = 1.0 / (inp.shell_width_km * 1e3)
    n = np.zeros_like(beta_grid)
    for i, w in mix:
        n += inp.N * (w / wsum) * g * latitude_pdf(beta_grid, i) / (2.0 * np.pi * d.R_bar**2 * np.cos(beta_grid))
    return beta_grid, n, d.n_mean


def collision_rate_vs_latitude(inp: Inputs, n_beta=4000, m_disp=21):
    """dR/dbeta (collisions per second per radian of latitude), whole fleet.
    Averaged over the inclination dispersion (m_disp sub-bins per population)
    so the integrable cap singularities are smoothed and the curve can be
    integrated numerically."""
    d = derive(inp)
    di = np.deg2rad(inp.incl_dispersion_deg)
    sub = np.linspace(-di, di, m_disp) if (di > 0 and m_disp > 1) else np.array([0.0])
    mix = []
    wsum = sum(w for _, w in inp.inclination_mix)
    for ideg, w in inp.inclination_mix:
        for off in sub:
            mix.append((np.deg2rad(ideg) + off, (w / wsum) / len(sub)))
    I_r = _radial_factor(inp, d)
    beta = np.linspace(-np.pi / 2 * 0.999, np.pi / 2 * 0.999, n_beta)
    dRdb = np.zeros_like(beta)
    for i1, w1 in mix:
        for i2, w2 in mix:
            p1 = latitude_pdf(beta, i1)
            p2 = latitude_pdf(beta, i2)
            v_a, v_b = pair_vrel_options(beta, i1, i2, d.v_orb)
            term = (inp.N * w1) * (inp.N * w2) * p1 * p2 * 0.5 * (v_a + v_b) / np.cos(beta)
            dRdb += 0.5 * d.sigma * I_r * np.where(np.isfinite(term), term, 0.0)
    return beta, dRdb


def impact_velocity_distribution(inp: Inputs, n_beta=20000, bins=120):
    """Rate-weighted histogram of impact speeds (the true collision-velocity
    spectrum). Returns bin centers (m/s) and pdf normalized to 1."""
    d = derive(inp)
    mix = [(np.deg2rad(ideg), w) for ideg, w in inp.inclination_mix]
    wsum = sum(w for _, w in mix)
    mix = [(i, w / wsum) for i, w in mix]
    vs, ws = [], []
    for i1, w1 in mix:
        for i2, w2 in mix:
            b_max = min(i_effective(i1), i_effective(i2))
            phi = (np.arange(n_beta) + 0.5) / n_beta * (np.pi / 2)
            beta = np.arcsin(np.sin(b_max) * np.sin(phi))
            dbeta = np.sin(b_max) * np.cos(phi) / np.cos(beta) * (np.pi / 2 / n_beta)
            p1, p2 = latitude_pdf(beta, i1), latitude_pdf(beta, i2)
            v_a, v_b = pair_vrel_options(beta, i1, i2, d.v_orb)
            base = w1 * w2 * p1 * p2 / np.cos(beta) * dbeta
            for v_opt in (v_a, v_b):
                vs.append(v_opt)
                ws.append(0.5 * base * v_opt)   # rate weight includes v itself
    vs = np.concatenate(vs); ws = np.concatenate(ws)
    ok = np.isfinite(vs) & np.isfinite(ws)
    hist, edges = np.histogram(vs[ok], bins=bins, range=(0, 2.05 * d.v_orb), weights=ws[ok], density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, hist


# ----------------------------------------------------------------------
# Kessler-cascade boundary (simplified branching criterion, parametric)
# ----------------------------------------------------------------------
def cascade_criterion(inp: Inputs, E_natural_per_year: float) -> dict:
    """Branching number kappa = expected secondary catastrophic collisions
    caused by the fragment cloud of one collision, before the fragments decay:
        kappa = F * n_sat * sigma_frag * v_rel * tau_frag
    Avoidance does not apply to the (largely untracked) fragment flux.
    kappa < 1: self-limiting debris growth. kappa >= 1: runaway (cascade).
    The boundary is also expressed as a critical satellite count N_crit."""
    d = derive(inp)
    sigma_f = inp.frag_shape_factor * inp.radiator_area_m2
    tau = inp.frag_lifetime_years * YEAR
    kappa = inp.frag_per_collision * d.n_mean * sigma_f * inp.v_rel_iso_ms * tau
    N_crit = inp.N / kappa if kappa > 0 else np.inf
    return {
        "kappa": kappa,
        "N_crit_for_kappa_1": N_crit,
        "regime": "RUNAWAY (Kessler cascade)" if kappa >= 1 else "self-limiting",
        "residual_collisions_per_year": E_natural_per_year * inp.avoidance_failure,
    }


# ----------------------------------------------------------------------
# Self-test: anchor checks from the problem statement
# ----------------------------------------------------------------------
if __name__ == "__main__":
    inp = Inputs()
    d = derive(inp)
    kin = kinetic_model(inp)
    print("=== Geometry ===")
    print(f"V = {d.V:.3e} m^3   (anchor ~1.86e20)")
    print(f"v_orb(R_bar) = {d.v_orb:.1f} m/s (anchor ~7500)")
    print("=== Kinetic baseline ===")
    print(f"n = {kin['n_mean_m3']:.2e} m^-3      (anchor ~4.3e-16)")
    print(f"nu = {kin['nu_per_s']:.2e} s^-1      (anchor ~2.1e-9)")
    print(f"nu = {kin['nu_per_year']:.4f} /yr    (anchor ~0.065)")
    print(f"P1 = {kin['P1_per_window']*100:.2f} %/yr (anchor ~6.5)")
    print(f"E  = {kin['E_collisions']:.0f} /yr    (anchor ~2600)")
    print(f"mfp = {kin['mean_free_path_m']/1e3:.2e} km (anchor ~9e10)")
    print(f"dilute regime: {kin['dilute_regime_ok']}")
    for E_acc in (100.0, 10.0, 1.0):
        dr_thin = required_shell_thickness_m(inp, E_acc, exact=False) / 1e3
        dr_exact = required_shell_thickness_m(inp, E_acc, exact=True) / 1e3
        print(f"E_acc={E_acc:6.0f}/yr: Delta_r thin-shell {dr_thin:.3g} km, exact {dr_exact:.3g} km")
    print("=== Keplerian model (Starlink-like mix) ===")
    kep = keplerian_model(inp)
    for k, v in kep.items():
        print(f"{k}: {v}")
    print("=== Cascade criterion ===")
    print(cascade_criterion(inp, kep["E_collisions"]))
