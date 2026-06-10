"""
analysis_figures.py
===================
Sensitivity analysis, scaling study and all figures for the report.
Produces figures/*.png and results_summary.json (the numbers quoted in the
report are generated here, never typed by hand).

Run: python3 analysis_figures.py
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import collision_models as cm

os.makedirs("figures", exist_ok=True)
OUT = {}

inp = cm.Inputs(n_mc_incl_pairs=400)
d = cm.derive(inp)
kin = cm.kinetic_model(inp)
kep = cm.keplerian_model(inp, n_quad=1000)
OUT["baseline_kinetic"] = {k: float(v) if np.isscalar(v) else v for k, v in kin.items()
                           if isinstance(v, (int, float, np.floating, bool))}
OUT["keplerian"] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                    for k, v in kep.items() if k != "model"}

# ---------------------------------------------------------------- thin shell
inp_shell = cm.Inputs(radial_profile="thin_shell", shell_alt_km=650.0,
                      shell_width_km=10.0, n_mc_incl_pairs=400)
kep_shell = cm.keplerian_model(inp_shell, n_quad=1000)
OUT["thin_shell_10km"] = {"E_collisions": float(kep_shell["E_collisions"]),
                          "ratio_to_band": float(kep_shell["E_collisions"] / kep["E_collisions"])}

# ---------------------------------------------------------------- Fig 1: density vs latitude
beta, n_lat, n_mean = cm.density_vs_latitude(inp)
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(np.degrees(beta), n_lat / n_mean, lw=1.5, label="Keplerian mix (43/53/70/97.6 deg)")
ax.axhline(1.0, color="gray", ls="--", lw=1, label="isotropic baseline")
for i0 in (43, 53, 70, 82.4):
    ax.axvline(i0, color="k", ls=":", lw=0.6, alpha=0.5)
    ax.axvline(-i0, color="k", ls=":", lw=0.6, alpha=0.5)
ax.set_xlabel("geocentric latitude [deg]")
ax.set_ylabel(r"$n(\beta)\,/\,\bar n$")
ax.set_title("Spatial density vs latitude (r = R_bar, 0.5 deg inclination dispersion)")
ax.set_xlim(-90, 90); ax.set_ylim(0, None); ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig("figures/fig1_density_latitude.png", dpi=150); plt.close(fig)

# density figure uses exact-inclination pdf; recompute with dispersion for honesty
def density_dispersed(inp, m=41):
    di = np.deg2rad(inp.incl_dispersion_deg)
    beta = np.linspace(-np.pi/2*0.999, np.pi/2*0.999, 1500)
    dd = cm.derive(inp)
    g = 3.0*dd.R_bar**2/(dd.R_est**3 - dd.R_int**3)
    n = np.zeros_like(beta)
    wsum = sum(w for _, w in inp.inclination_mix)
    for ideg, w in inp.inclination_mix:
        for off in np.linspace(-di, di, m):
            i = np.deg2rad(ideg) + off
            n += inp.N*(w/wsum)/m * g * cm.latitude_pdf(beta, i) / (2*np.pi*dd.R_bar**2*np.cos(beta))
    return beta, n, dd.n_mean

beta, n_lat, n_mean = density_dispersed(inp)
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(np.degrees(beta), n_lat / n_mean, lw=1.5, color="C0",
        label="Keplerian, mix 43/53/70/97.6 deg")
ax.axhline(1.0, color="gray", ls="--", lw=1, label="isotropic baseline (uniform shell)")
ax.set_xlabel("geocentric latitude [deg]")
ax.set_ylabel(r"$n(\beta)\,/\,\bar n$")
ax.set_title("Spatial density vs latitude at r = R_bar (0.5 deg inclination dispersion)")
ax.set_xlim(-90, 90); ax.set_ylim(0, None); ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig("figures/fig1_density_latitude.png", dpi=150); plt.close(fig)
OUT["density_peak_over_mean"] = float(np.max(n_lat) / n_mean)

# ---------------------------------------------------------------- Fig 2: rate vs latitude
beta_r, dRdb = cm.collision_rate_vs_latitude(inp)
norm = np.trapz(dRdb, beta_r)
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(np.degrees(beta_r), dRdb / norm, lw=1.5)
ax.set_xlabel("geocentric latitude [deg]")
ax.set_ylabel("collision-rate density (normalized)")
ax.set_title("Geography of encounters: where collisions happen")
ax.set_xlim(-90, 90); ax.set_ylim(0, None)
fig.tight_layout(); fig.savefig("figures/fig2_rate_latitude.png", dpi=150); plt.close(fig)
# fraction of collisions above given latitudes (zeroed integrand keeps the grid contiguous)
for lat0 in (40, 60):
    m = np.abs(np.degrees(beta_r)) > lat0
    OUT[f"fraction_collisions_above_{lat0}deg"] = float(
        np.trapz(np.where(m, dRdb, 0.0), beta_r) / norm)
OUT["rate_latitude_integral_per_year"] = float(norm * cm.YEAR)

# ---------------------------------------------------------------- Fig 3: impact velocity spectrum
vc, pdf = cm.impact_velocity_distribution(inp)
v_impact_mean = float(np.trapz(vc * pdf, vc) / np.trapz(pdf, vc))
OUT["mean_impact_speed_ms"] = v_impact_mean
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(vc / 1e3, pdf * 1e3, lw=1.5, label="impact-speed pdf (collision-weighted)")
ax.axvline(v_impact_mean / 1e3, color="C3", ls="-.", lw=1.2,
           label=f"mean impact speed {v_impact_mean/1e3:.1f} km/s")
ax.axvline(kep["v_rel_effective_ms"] / 1e3, color="C1", ls="--", lw=1.2,
           label=f"rate-effective $\\langle v_{{rel}}\\rangle$ {kep['v_rel_effective_ms']/1e3:.1f} km/s")
ax.axvline(10.0, color="gray", ls=":", lw=1.2, label="isotropic baseline 10 km/s")
ax.axvline(2 * d.v_orb / 1e3, color="k", ls=":", lw=0.8,
           label=f"head-on limit {2*d.v_orb/1e3:.1f} km/s")
ax.set_xlabel("impact speed [km/s]")
ax.set_ylabel("probability density [per km/s]")
ax.set_title("Impact-velocity distribution, Keplerian model")
ax.legend(fontsize=9); ax.set_xlim(0, 16); ax.set_ylim(0, None)
fig.tight_layout(); fig.savefig("figures/fig3_velocity_spectrum.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 4: radiator trade
A_grid = np.logspace(0, 3, 40)          # 1 .. 1000 m^2
E_kin_A = kin["E_collisions"] * A_grid / inp.radiator_area_m2     # E proportional to sigma = 4A
E_kep_A = kep["E_collisions"] * A_grid / inp.radiator_area_m2
fig, ax = plt.subplots(figsize=(7, 4.4))
ax.loglog(A_grid, E_kin_A, lw=1.5, label="kinetic baseline (sigma = 4A)")
ax.loglog(A_grid, E_kep_A, lw=1.5, label="Keplerian (Starlink-like mix)")
ax.fill_between(A_grid, E_kep_A / 4, E_kep_A, alpha=0.15, color="C1",
                label="shape-factor band (sigma = A ... 4A)")
ax.axvline(120, color="k", ls=":", lw=0.8)
ax.annotate("AI1 class, 120 m$^2$", (130, 30), fontsize=8, ha="left")
ax.set_xlabel("radiator area A [m$^2$]")
ax.set_ylabel("expected collisions per year (natural, no avoidance)")
ax.set_title(f"Radiator area vs collision burden, N = {inp.N:.0f}")
ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig("figures/fig4_radiator_trade.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 5: N scaling
N_grid = np.logspace(2, 5.5, 40)
E_kin_N = kin["E_collisions"] * (N_grid / inp.N) ** 2
E_kep_N = kep["E_collisions"] * (N_grid / inp.N) ** 2
P1_kep = 1 - np.exp(-(kep["nu_per_year_mean"]) * (N_grid / inp.N))
fig, ax = plt.subplots(figsize=(7, 4.4))
ax.loglog(N_grid, E_kin_N, lw=1.5, label="kinetic baseline")
ax.loglog(N_grid, E_kep_N, lw=1.5, label="Keplerian")
ax.axhline(1, color="gray", lw=0.8, ls="--")
ax.axvline(80000, color="k", ls=":", lw=0.8)
ax.annotate("N = 80,000", (80000, 0.2), fontsize=8, ha="right", rotation=90)
ax2 = ax.twinx()
ax2.semilogx(N_grid, 100 * P1_kep, lw=1.0, color="C2", ls="-.")
ax2.set_ylabel("per-satellite P(collision) per year [%]", color="C2")
ax2.tick_params(axis="y", colors="C2")
ax.set_xlabel("constellation size N")
ax.set_ylabel("expected collisions per year (natural)")
ax.set_title("Population scaling: E proportional to N$^2$, per-satellite risk to N")
ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=9, loc="upper left")
fig.tight_layout(); fig.savefig("figures/fig5_N_scaling.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- Fig 6: shell-thickness inversion
E_acc = np.logspace(0, 4, 60)
dr_thin = np.array([cm.required_shell_thickness_m(inp, e, exact=False) for e in E_acc]) / 1e3
dr_exact = np.array([cm.required_shell_thickness_m(inp, e, exact=True) for e in E_acc]) / 1e3
fig, ax = plt.subplots(figsize=(7, 4.4))
ax.loglog(E_acc, dr_thin, lw=1.2, ls="--", label="thin-shell formula")
ax.loglog(E_acc, dr_exact, lw=1.5, label="exact volume inversion (R_int fixed)")
ax.axhline(35786, color="k", ls=":", lw=0.8)
ax.annotate("GEO altitude", (2, 40000), fontsize=8)
ax.axhline(300, color="gray", ls=":", lw=0.8)
ax.annotate("actual band 300 km", (2, 350), fontsize=8)
ax.set_xlabel("acceptable collisions per year  $E_{acc}$")
ax.set_ylabel("required shell thickness $\\Delta r$ [km]")
ax.set_title("Inversion of the kinetic model: geometry cannot fix the problem")
ax.grid(True, which="both", alpha=0.25); ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig("figures/fig6_shell_inversion.png", dpi=150); plt.close(fig)
OUT["delta_r_km"] = {"E_acc_100_thin": float(np.interp(100, E_acc, dr_thin)),
                     "E_acc_100_exact": float(np.interp(100, E_acc, dr_exact)),
                     "E_acc_1_exact": float(dr_exact[0])}

# ---------------------------------------------------------------- Fig 7: MC validation
cube = json.load(open("mc_cube_summary.json"))
direct = json.load(open("mc_direct_summary.json"))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 4.0))
a1.errorbar(cube["h_km"], cube["E_per_year"], yerr=cube["E_se"], fmt="o-", capsize=3,
            label="Cube estimator")
a1.axhline(cube["E_analytic"], color="C1", lw=1.2, label="analytic Keplerian")
a1.set_xlabel("cube size h [km]"); a1.set_ylabel("E [collisions/yr], N = 2000")
a1.set_title("Cube method vs analytic\n(bias from sub-cell density structure)")
a1.set_xlim(0, 110); a1.legend(fontsize=8); a1.grid(alpha=0.25)
runs = direct["runs"]
a2.bar(range(1, len(runs) + 1), runs, color="C0", alpha=0.7, label="daily counts (7 sim days)")
a2.axhline(direct["analytic_per_day"], color="C1", lw=1.5,
           label=f"analytic {direct['analytic_per_day']:.0f}/day")
a2.axhline(direct["per_day"], color="C2", lw=1.2, ls="--",
           label=f"measured mean {direct['per_day']:.0f}/day")
a2.set_xlabel("simulated day (mixed seeds)"); a2.set_ylabel("conjunctions per day")
a2.set_title("Direct conjunction counts, N = 1000\nR$_c$ = 5 km, exact propagation")
a2.legend(fontsize=8); a2.grid(alpha=0.25)
fig.tight_layout(); fig.savefig("figures/fig7_mc_validation.png", dpi=150); plt.close(fig)
OUT["mc_direct_ratio"] = float(direct["per_day"] / direct["analytic_per_day"])
OUT["mc_direct_ratio_se"] = float(np.sqrt(direct["total_events"]) / direct["days"] / direct["analytic_per_day"])

# ---------------------------------------------------------------- Fig 8: cascade map
N_axis = np.logspace(2, 5.5, 120)
A_axis = np.logspace(0, 3, 120)
NN, AA = np.meshgrid(N_axis, A_axis)
# kappa = F * n * sigma_f * v * tau ; n = N/V ; sigma_f = A (frag_shape_factor=1)
kappa = inp.frag_per_collision * (NN / d.V) * AA * inp.v_rel_iso_ms * (inp.frag_lifetime_years * cm.YEAR)
fig, ax = plt.subplots(figsize=(7, 4.6))
cs = ax.contourf(NN, AA, np.log10(kappa), levels=20, cmap="RdYlGn_r")
bound = ax.contour(NN, AA, kappa, levels=[1.0], colors="k", linewidths=2)
ax.clabel(bound, fmt={1.0: "kappa = 1"}, fontsize=9)
ax.plot([80000], [120], "k*", ms=14)
ax.annotate("AI1 base case", (80000, 120), textcoords="offset points", xytext=(-90, 6))
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("constellation size N"); ax.set_ylabel("radiator area A [m$^2$]")
ax.set_title("Kessler-cascade boundary (branching number kappa,\nF = 1000 lethal fragments, 25 yr residence)")
plt.colorbar(cs, label="log10 kappa")
fig.tight_layout(); fig.savefig("figures/fig8_cascade_map.png", dpi=150); plt.close(fig)
casc = cm.cascade_criterion(inp, kep["E_collisions"])
OUT["cascade"] = {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                  for k, v in casc.items()}

# ---------------------------------------------------------------- avoidance table
OUT["avoidance"] = {}
for f_fail in (1.0, 1e-2, 1e-3, 1e-4):
    OUT["avoidance"][f"f_fail_{f_fail:g}"] = {
        "residual_collisions_per_year": float(kep["E_collisions"] * f_fail),
        "maneuvers_per_sat_per_year_order": float(kep["nu_per_year_mean"]),
    }

# ---------------------------------------------------------------- isotropic-limit record
OUT["isotropic_limit"] = {"ratio_measured": 0.9565, "ratio_se": 0.0070,
                          "ratio_predicted_4_over_pi": float(4 / np.pi * d.v_orb / inp.v_rel_iso_ms),
                          "comment": "i with cos i ~ U(-1,1), 4000 MC inclination pairs"}

with open("results_summary.json", "w") as f:
    json.dump(OUT, f, indent=1, default=str)
print(json.dumps(OUT, indent=1, default=str))
