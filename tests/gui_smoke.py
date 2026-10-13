"""Click through the real wizard with hospital inputs and check the results.

Run from the repo root: python tests/gui_smoke.py
References: Theagenio INFINITY door sheet (H_PS) and the NCRP 151 / Table B.5a
direct-door values. Exits non-zero on any mismatch.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import shielding_wizard as w

VARS = ("_ack", "_single_energy", "_imrt_used", "_tbi_used", "_maze_exists",
        "_primary_reaches", "_primary_at_door")


def fill(step, vals):
    for k, v in vals.items():
        if k in VARS:
            getattr(step, k).set(v)
            if k in ("_imrt_used", "_tbi_used"):
                step._toggle()
        elif k in getattr(step, "_entries", {}):
            e = step._entries[k]
            e.delete(0, "end")
            e.insert(0, str(v))
            e._check()
        elif k in getattr(step, "_combos", {}):
            step._combos[k].set(str(v))
        else:
            raise KeyError(f"{type(step).__name__} has no input {k}")


def run(pages):
    app = w.ShieldingWizard()
    app.withdraw()
    app.start_radiotherapy_branch()
    for vals in pages:
        cur = app.steps[app.current_step_index]
        fill(cur, vals)
        if isinstance(cur, w.Step_RT_Checklist):
            cur._try_advance()
        else:
            app.next_step()
        if app.steps[app.current_step_index] is cur:
            raise RuntimeError(f"stuck on {type(cur).__name__}")
    data = dict(app.data)
    app.destroy()
    return data


def check(name, got, want, rel=0.01):
    ok = math.isclose(got, want, rel_tol=rel)
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {got:.4e} (expected {want:.4e})")
    return ok


ALB = dict(alpha_z_material="concrete", alpha_z_angle=0, alpha_z_scatter_angle=75,
           alpha0_material="concrete", alpha0_angle=0, alpha0_scatter_angle=75,
           alpha1_material="concrete", alpha1_angle=45, alpha1_scatter_angle=0)
SETUP = [
    dict(_ack=True),
    dict(_single_energy=True, E_max="18", therapies_per_day_max=50, days_per_week_max=5, Grays_max=4),
    dict(_imrt_used=True, p_imrt_min=0, C_I_min=1, p_imrt_max=1, C_I_max=2.5),
    dict(_tbi_used=False),
]

# Theagenio INFINITY maze door, 18 MV: W = 1000 (scatter), W_L = 2500 (leakage).
maze = run(SETUP + [
    dict(_maze_exists=True),
    dict(_primary_reaches=False),
    dict(P=2e-5, U_G=0.25, L_f=0.001, f=1.0, F=1600, h_e=3.6, w_1=0.15625, w_tongue=1.8, t=1.5,
         d_pp=3.0, d_r=7.2, d_sec=7.5, d_sca=1.0, d_L=6.5, d_z=5.7, d_zz=7.2, i_theta=0, s_theta=45,
         door_energy="1 MeV", **ALB),
    {},
])

# Theagenio secondary wall through the direct-door branch: 90 deg, 5 m, 18 MV.
direct = run(SETUP + [
    dict(_maze_exists=False),
    dict(_primary_at_door=False),
    dict(P=2e-4, T=0.25, L_f=0.001, F=1600, d_sca=1.0, d_sec=5.0, d_L=5.0, theta=90),
    {},
])

results = [
    check("W (scatter)", maze["W_max"], 1000.0, 1e-9),
    check("W_L (leakage)", maze["W_L_max"], 2500.0, 1e-9),
    check("H_PS vs hospital 4.224e-05", maze["H_PS"], 4.224e-5),
    check("H_ph vs hospital 1.273e-04 (+2% conservative)", maze["H_ph"], 1.273e-4, 0.05),
    check("direct door leakage barrier vs hospital 72.2 cm", direct["t_L"], 72.2, 0.03),
]
sys.exit(0 if all(results) else 1)
