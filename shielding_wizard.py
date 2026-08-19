"""
RT Bunker Door Shielding Calculation Wizard
===========================================
A step-by-step GUI for calculating radiation shielding in radiotherapy
and diagnostic radiology scenarios.

Reference standards used:
    1) NCRP Report No. 151 - Structural Shielding Design and Evaluation
       for Megavoltage X- and Gamma-Ray Radiotherapy Facilities
    2) Radiation Shielding for Diagnostic Radiology, 2nd Edition (Sutton, BIR)

If you are using a different standard, that is fine.
We may or may not call the Ordo Hereticus.
This work is designed for the specific calculations that suit
the author's location and clinical standards. Adapt at your own risk.
(Realistically, it is tailored to the point that starting from scratch may be simpler. Consider yourself warned.)

Workflow branches:
    SETUP    : Common steps - workload, IMRT correction, TBI, maze existence check.
               Every branch passes through here.
    M-branch : Linac photon, maze geometry, E_max < 10 MV.
    N-branch : Linac photon, maze geometry, E_max >= 10 MV.
               Extends M with neutron and capture gamma at door.
    D-branch : Linac photon, direct door, no maze.
    CT-branch: Diagnostic CT room shielding (Sutton/BIR). Pending.

Units contract.
    Distances        : m
    Areas            : m2
    Workloads        : Gy/week
    Design goal P    : Sv/week
    Dose rates Ddot0 : Gy/h
    IDR values       : Sv/h
    Field area F     : cm2  - NCRP normalisation. Do NOT convert to m2.
                        It has been done before. It did not end well.

Keep this in mind for proper calculations.
I will probably and hopefully not be the one irradiated if you fail this.
Unless you are building something that warrants its own exclusion zone.
If that is the case, I sincerely hope you are not using a wizard GUI for the shielding calculations.
If you are, have you tried therapy, a confession to your local priest, or Monte Carlo? I've heard they work better.

No external dependencies. Python stdlib only: math, tkinter.
If something breaks, the problem is in this file. At least it is easy to find.

Author: A severely underslept, overworked and severely caffeinated medical and radiation physicist. AKA notcosmonaut

"""


from __future__ import annotations
import math
import tkinter as tk
from tkinter import ttk, messagebox
from typing import Optional


# ==============================================================================
# TVL TABLE - Tenth Value Layers for primary barrier (concrete, cm)
# Source: NCRP Report No. 151, Table B.2
# ------------------------------------------------------------------------------
# TVL1  : first TVL - accounts for beam hardening in the initial layer
# TVLe  : equilibrium TVL - for all subsequent layers after the beam hardens
# Energy in MV. Co-60 is a special case (1.25 MeV photons), handled separately.
# Units: cm of ordinary concrete.
# Used by the D-branch barrier thickness calculation.
# The table carries more energies than the SETUP dropdown offers. That is
# deliberate: the dropdown is limited by the albedo tables, not by this one.
# ==============================================================================

_TVL_TABLE = [
    # (MV,  TVL1, TVLe)
    (4,    33,   28),
    (6,    34,   29),
    (10,   35,   31),
    (15,   36,   33),
    (18,   36,   34),
    (20,   36,   34),
    (25,   37,   35),
    (30,   37,   36),
]

# Co-60: monoenergetic 1.25 MeV photons, no spectrum to harden.
# TVL1 == TVLe for exactly that reason.
_CO60_TVL = (1.25, 21, 21)


def lookup_tvl(energy_mv: float) -> tuple[float, float]:
    """
    Return (TVL1, TVLe) in cm of concrete for a given beam energy in MV.

    Anything not in the table raises ValueError - linacs come in discrete
    energies, not a continuous spectrum. An energy not in this table means
    something went wrong upstream, not here.
    """
    if energy_mv <= 1.25:
        return _CO60_TVL[1], _CO60_TVL[2]

    for row in _TVL_TABLE:
        if energy_mv == row[0]:
            return row[1], row[2]

    raise ValueError(
        f"Energy {energy_mv} MV not found in TVL table. "
        f"Valid energies: {[row[0] for row in _TVL_TABLE]}"
    )


# Leakage and scattered-photon TVLs, Table B.5, NCRP 151, concrete only.
# The primary-beam TVL table is the WRONG table for the D-branch. Leaked and
# scattered photons are softer than the primary beam and attenuate faster,
# with their own TVLs, angle-dependent for scatter. Digitised at 6 and 18 MV
# only, concrete only. 10 MV is not covered -- if selected, the D-branch
# raises a clear error rather than silently reusing an adjacent energy's numbers.
_TVL_LEAKAGE_SCATTER_TABLE = {
    # energy_mv : (TVLl1, TVLle, TVLs20, TVLs90)  -- concrete, cm
    6:  (34, 29, 31.3, 17),
    18: (36, 34, 40.0, 19),
}


def lookup_tvl_leakage(energy_mv: float) -> tuple[float, float]:
    """Return (TVLl1, TVLle) in cm of concrete for head-leakage photons."""
    row = _TVL_LEAKAGE_SCATTER_TABLE.get(energy_mv)
    if row is None:
        raise ValueError(
            f"No leakage TVL data for {energy_mv} MV. Only 6 and 18 MV concrete "
            f"are digitised (Table B.5). 10 MV is a known gap, not a bug."
        )
    return row[0], row[1]


def lookup_tvl_scatter(energy_mv: float, scatter_angle_deg: float) -> float:
    """
    Return a single TVL in cm of concrete for scattered photons at a given
    scatter angle. Rounds to the nearest of the two measured angles (20, 90),
    midpoint at 55 deg.
    """
    row = _TVL_LEAKAGE_SCATTER_TABLE.get(energy_mv)
    if row is None:
        raise ValueError(
            f"No scattered-photon TVL data for {energy_mv} MV. Only 6 and 18 MV "
            f"concrete are digitised (Table B.5). 10 MV is a known gap, not a bug."
        )
    _, _, tvl_s20, tvl_s90 = row
    return tvl_s20 if scatter_angle_deg < 55.0 else tvl_s90


# ==============================================================================
# SCATTER FRACTION TABLE - Table B.4, NCRP 151
# ------------------------------------------------------------------------------
# Patient scatter fractions a(theta) at 1 m from a human-size phantom.
# Target-to-phantom distance: 1 m, field size: 400 cm2.
# Only four energies measured: 6, 10, 18, 24 MV. No Co-60.
# Structure: {scatter_angle: {energy_mv: fraction}}
# Values are direct scatter fractions, NOT x10^-3.
# Used for the patient scatter component in both the maze and direct-door paths.
# ==============================================================================

ALBEDO_SCATTER_TABLE_B4 = {
    #  angle : {6 MV,     10 MV,    18 MV,    24 MV}
    10:  {6: 1.04e-2, 10: 1.66e-2, 18: 1.42e-2, 24: 1.78e-2},
    20:  {6: 6.73e-3, 10: 5.79e-3, 18: 5.39e-3, 24: 6.32e-3},
    30:  {6: 2.77e-3, 10: 3.18e-3, 18: 2.53e-3, 24: 2.74e-3},
    45:  {6: 1.39e-3, 10: 1.35e-3, 18: 8.64e-4, 24: 8.30e-4},
    60:  {6: 8.24e-4, 10: 7.46e-4, 18: 4.24e-4, 24: 3.86e-4},
    90:  {6: 4.26e-4, 10: 3.81e-4, 18: 1.89e-4, 24: 1.74e-4},
    135: {6: 3.00e-4, 10: 3.02e-4, 18: 1.24e-4, 24: 1.20e-4},
    150: {6: 2.87e-4, 10: 2.74e-4, 18: 1.20e-4, 24: 1.13e-4},
}

# The only energies B.4 was measured at. Used to validate the dropdown selection.
_B4_ENERGIES = (6, 10, 18, 24)


def lookup_scatter_fraction_b4(scatter_angle: float, energy_mv: float) -> float:
    """
    Return patient scatter fraction a(theta) from NCRP 151 Table B.4.

    Valid energies: 6, 10, 18, 24 MV (dropdown upstream).
    Valid angle range: 10 to 150 degrees.

    Scatter angle rounded to nearest measured value, midpoints rounding up:
        < 15               -> 10
        15 to < 25         -> 20
        25 to < 37.5       -> 30
        37.5 to < 52.5     -> 45
        52.5 to < 75       -> 60
        75 to < 112.5      -> 90
        112.5 to < 142.5   -> 135
        142.5 to 150       -> 150

    Values are direct scatter fractions, NOT x10^-3. This table is for the
    patient scatter component only - do not use it for wall albedo.
    """
    # Past 150 degrees the data runs out. At that point, just use Monte Carlo.
    if scatter_angle < 10.0 or scatter_angle > 150.0:
        raise ValueError(
            f"Scatter angle {scatter_angle} deg is outside the measured range [10, 150]."
        )

    if energy_mv not in _B4_ENERGIES:
        raise ValueError(
            f"Energy {energy_mv} MV not in B.4 table. Valid energies: 6, 10, 18, 24 MV."
        )

    # Round scatter angle to nearest measured value (midpoints round up).
    if 10.0 <= scatter_angle < 15.0:
        angle_key = 10
    elif 15.0 <= scatter_angle < 25.0:
        angle_key = 20
    elif 25.0 <= scatter_angle < 37.5:
        angle_key = 30
    elif 37.5 <= scatter_angle < 52.5:
        angle_key = 45
    elif 52.5 <= scatter_angle < 75.0:
        angle_key = 60
    elif 75.0 <= scatter_angle < 112.5:
        angle_key = 90
    elif 112.5 <= scatter_angle < 142.5:
        angle_key = 135
    else:  # 142.5 to 150
        angle_key = 150

    return ALBEDO_SCATTER_TABLE_B4[angle_key][energy_mv]


# ==============================================================================
# POTATO'S PATIENT-SCATTER ALBEDO -- deliberate, isolated exception
# ------------------------------------------------------------------------------
# H_PS uses a SECOND, DIFFERENT lookup of Table B.4 than the one above. Potato
# ran the B.4 table (keyed angle-then-energy) through its wall-albedo lookup
# function, which does true bilinear interpolation over BOTH axes with no
# exact-match requirement. Since B.4's outer keys are angles and inner keys
# are energies, this reads the angle argument against the energy axis and
# vice versa -- semantically wrong, but it is what the verified reference
# case actually used, and it reproduces the confirmed thickness (0.71 cm)
# exactly. Our own nearest-neighbour lookup_scatter_fraction_b4 above cannot
# reproduce this number: it is a different algorithm (bucket rounding, not
# interpolation), not just a different table.
#
# This function exists ONLY to reproduce that one verified quantity. It is
# not a general-purpose replacement for lookup_albedo or lookup_scatter_
# fraction_b4, and should not be reused elsewhere. "No interpolation" stays
# the rule everywhere else in this file; this is the one deliberate,
# documented exception, kept because potato worked.
# ==============================================================================

def lookup_alpha_theta_potato(scatter_angle_deg: float, energy_mv: float) -> float:
    """
    Reproduces potato's H_PS scatter coefficient exactly: bilinear
    interpolation of ALBEDO_SCATTER_TABLE_B4 with angle and energy read
    against each other's axis. Verified: (10, 18) -> 0.007952.
    """
    table = ALBEDO_SCATTER_TABLE_B4
    energies = sorted(table.keys())
    angles = sorted(next(iter(table.values())).keys())

    lower_e = upper_e = None
    for e in energies:
        if energy_mv == e:
            lower_e = upper_e = e
            break
        if energy_mv > e:
            lower_e = e
        elif energy_mv < e and upper_e is None:
            upper_e = e
    if lower_e is None:
        lower_e = upper_e = energies[0]
    if upper_e is None:
        lower_e = upper_e = energies[-1]

    lower_a = upper_a = None
    for a in angles:
        if scatter_angle_deg == a:
            lower_a = upper_a = a
            break
        if scatter_angle_deg > a:
            lower_a = a
        elif scatter_angle_deg < a and upper_a is None:
            upper_a = a
    if lower_a is None:
        lower_a = upper_a = angles[0]
    if upper_a is None:
        lower_a = upper_a = angles[-1]

    def gv(e, a):
        return table[e][a]

    if lower_e == upper_e and lower_a == upper_a:
        return gv(lower_e, lower_a)
    elif lower_e == upper_e:
        v0, v1 = gv(lower_e, lower_a), gv(lower_e, upper_a)
        t = (scatter_angle_deg - lower_a) / (upper_a - lower_a)
        return v0 + t * (v1 - v0)
    elif lower_a == upper_a:
        v0, v1 = gv(lower_e, lower_a), gv(upper_e, lower_a)
        t = (energy_mv - lower_e) / (upper_e - lower_e)
        return v0 + t * (v1 - v0)
    else:
        v00, v01 = gv(lower_e, lower_a), gv(lower_e, upper_a)
        v10, v11 = gv(upper_e, lower_a), gv(upper_e, upper_a)
        te = (energy_mv - lower_e) / (upper_e - lower_e)
        ta = (scatter_angle_deg - lower_a) / (upper_a - lower_a)
        v0 = v00 + ta * (v01 - v00)
        v1 = v10 + ta * (v11 - v10)
        return v0 + te * (v1 - v0)


# ==============================================================================
# WALL ALBEDO TABLES - Differential dose albedo coefficients
# Source: NCRP 151, Tables B.8a-f
# ------------------------------------------------------------------------------
# Six tables: three materials (concrete, iron, lead) x two incident angles (0, 45).
# Structure: {energy_mv: {scatter_angle: coefficient}}
# The x10^-3 factor is baked in at storage - what comes out is ready to use.
# Incident angle is chosen by selecting which table you pass to lookup_albedo,
# so there is deliberately no incident-angle argument in that function.
# Used for maze wall scatter components (H_s, H_LS, H_PS, H_LT).
# ==============================================================================


# Ordinary concrete, normal incidence (0) - Table B.8a
ALBEDO_CONCRETE_0_TABLE = {
    # energy_mv : {scatter_angle: coefficient}
    1.25: {0: 7.0e-3, 30: 6.5e-3, 45: 6.0e-3, 60: 5.5e-3, 75: 3.8e-3},  # Co-60
    4:    {0: 6.7e-3, 30: 6.4e-3, 45: 5.8e-3, 60: 4.9e-3, 75: 3.1e-3},
    6:    {0: 5.3e-3, 30: 5.2e-3, 45: 4.7e-3, 60: 4.0e-3, 75: 2.7e-3},
    10:   {0: 4.3e-3, 30: 4.1e-3, 45: 3.8e-3, 60: 3.1e-3, 75: 2.1e-3},
    18:   {0: 3.4e-3, 30: 3.4e-3, 45: 3.0e-3, 60: 2.5e-3, 75: 1.6e-3},
    24:   {0: 3.2e-3, 30: 3.2e-3, 45: 2.8e-3, 60: 2.3e-3, 75: 1.5e-3},
    30:   {0: 3.0e-3, 30: 2.7e-3, 45: 2.6e-3, 60: 2.2e-3, 75: 1.5e-3},
}

# Ordinary concrete, 45 incidence - Table B.8b
# No Co-60 entry in NCRP 151 for this configuration.
ALBEDO_CONCRETE_45_TABLE = {
    # energy_mv : {scatter_angle: coefficient}
    4:    {0: 7.6e-3, 30: 8.5e-3, 45: 9.0e-3, 60: 9.2e-3, 75: 9.5e-3},
    6:    {0: 6.4e-3, 30: 7.1e-3, 45: 7.3e-3, 60: 7.7e-3, 75: 8.0e-3},
    10:   {0: 5.1e-3, 30: 5.7e-3, 45: 5.8e-3, 60: 6.0e-3, 75: 6.0e-3},
    18:   {0: 4.5e-3, 30: 4.6e-3, 45: 4.6e-3, 60: 4.3e-3, 75: 4.0e-3},
    24:   {0: 3.7e-3, 30: 3.9e-3, 45: 3.9e-3, 60: 3.7e-3, 75: 3.4e-3},
    30:   {0: 4.8e-3, 30: 5.0e-3, 45: 4.9e-3, 60: 4.0e-3, 75: 3.0e-3},
}

# Iron, normal incidence (0) - Table B.8c
ALBEDO_IRON_0_TABLE = {
    # energy_mv : {scatter_angle: coefficient}
    4:    {0: 6.0e-3, 30: 5.4e-3, 45: 5.1e-3, 60: 4.8e-3, 75: 3.1e-3},
    6:    {0: 5.5e-3, 30: 4.9e-3, 45: 4.7e-3, 60: 4.2e-3, 75: 2.8e-3},
    10:   {0: 5.0e-3, 30: 4.5e-3, 45: 4.3e-3, 60: 3.9e-3, 75: 2.5e-3},
    18:   {0: 5.1e-3, 30: 4.5e-3, 45: 4.3e-3, 60: 3.8e-3, 75: 2.4e-3},
    30:   {0: 5.5e-3, 30: 4.7e-3, 45: 4.4e-3, 60: 3.8e-3, 75: 2.3e-3},
}

# Iron, 45 incidence - Table B.8d
ALBEDO_IRON_45_TABLE = {
    # energy_mv : {scatter_angle: coefficient}
    4:    {0: 7.1e-3, 30: 8.1e-3, 45: 10.0e-3, 60: 10.6e-3, 75: 11.5e-3},
    6:    {0: 6.0e-3, 30: 7.0e-3, 45: 8.5e-3,  60: 9.0e-3,  75: 9.5e-3},
    10:   {0: 6.1e-3, 30: 6.8e-3, 45: 7.1e-3,  60: 7.2e-3,  75: 7.2e-3},
    18:   {0: 6.5e-3, 30: 6.4e-3, 45: 6.2e-3,  60: 6.0e-3,  75: 5.6e-3},
    30:   {0: 6.6e-3, 30: 6.5e-3, 45: 6.3e-3,  60: 5.5e-3,  75: 4.6e-3},
}

# Lead, normal incidence (0) - Table B.8e
ALBEDO_LEAD_0_TABLE = {
    # energy_mv : {scatter_angle: coefficient}
    4:    {0: 5.9e-3, 30: 5.2e-3, 45: 4.7e-3, 60: 4.2e-3, 75: 3.0e-3},
    6:    {0: 5.0e-3, 30: 4.5e-3, 45: 4.2e-3, 60: 3.8e-3, 75: 2.6e-3},
    10:   {0: 4.5e-3, 30: 3.9e-3, 45: 3.6e-3, 60: 3.2e-3, 75: 2.2e-3},
    18:   {0: 3.9e-3, 30: 3.4e-3, 45: 3.2e-3, 60: 2.8e-3, 75: 1.8e-3},
    30:   {0: 3.5e-3, 30: 3.0e-3, 45: 2.7e-3, 60: 2.4e-3, 75: 1.5e-3},
}

# Lead, 45 incidence - Table B.8f
ALBEDO_LEAD_45_TABLE = {
    # energy_mv : {scatter_angle: coefficient}
    4:    {0: 6.5e-3, 30: 7.6e-3, 45: 8.3e-3, 60: 8.6e-3, 75: 9.0e-3},
    6:    {0: 6.5e-3, 30: 6.8e-3, 45: 7.0e-3, 60: 7.3e-3, 75: 7.8e-3},
    10:   {0: 5.4e-3, 30: 5.8e-3, 45: 6.0e-3, 60: 5.9e-3, 75: 5.8e-3},
    18:   {0: 4.9e-3, 30: 5.0e-3, 45: 5.0e-3, 60: 4.8e-3, 75: 4.5e-3},
    30:   {0: 4.1e-3, 30: 4.2e-3, 45: 4.1e-3, 60: 3.7e-3, 75: 3.2e-3},
}


def lookup_albedo(table: dict, scatter_angle: float, energy_mv: float) -> float:
    """
    Return albedo coefficient for given material table, scatter angle, and beam energy.
    Incident angle is encoded by table choice - no argument needed here.
    Energy: exact match only, dropdown upstream.
    Scatter angle rounded to nearest table value, midpoints round up.
    Values returned are final coefficients - x10^-3 already baked in.
    """
    # Past 90 degrees the tables end. At that point, Monte Carlo is the honest answer.
    if scatter_angle < 0.0 or scatter_angle > 90.0:
        raise ValueError(
            f"Scatter angle {scatter_angle} deg is outside the supported range [0, 90]."
        )

    if energy_mv not in table:
        raise ValueError(
            f"Energy {energy_mv} MV not found in albedo table. "
            f"Valid energies: {sorted(table.keys())}"
        )

    # Round scatter angle to nearest measured value (midpoints round up).
    if 0.0 <= scatter_angle < 15.0:
        angle_key = 0
    elif 15.0 <= scatter_angle < 37.5:
        angle_key = 30
    elif 37.5 <= scatter_angle < 52.5:
        angle_key = 45
    elif 52.5 <= scatter_angle < 67.5:
        angle_key = 60
    elif 67.5 <= scatter_angle <= 90.0:
        angle_key = 75
    else:
        # Unreachable given the range check above. If you are reading this in a
        # stack trace: light a candle, some myrrh, and pray over the input data.
        # If it still fails, consult the names theory book - with less enthusiasm.
        raise ValueError(
            f"Scatter angle {scatter_angle} deg slipped through range validation."
        )

    return table[energy_mv][angle_key]


# ==============================================================================
# SOFT-SPECTRUM ALBEDO CONSTANTS
# ------------------------------------------------------------------------------
# Radiation on its second or later reflection has degraded to roughly 0.5 MeV
# whatever the beam's nominal energy, so those lookups fix the energy axis
# instead of tying it to E_max/E_min.
# ==============================================================================

# TODO: source NCRP 151 albedo value at true 0.5 MeV. Currently using the
# Co-60 (1.25 MeV) row as the nearest tabulated stand-in.
_SOFT_SPECTRUM_MV = 1.25

# Patient-scattered radiation arrives at the maze wall equally soft, so its
# wall albedo is one fixed value, not an energy-split pair.
# 0.5 MeV, 25 deg incidence, 0 deg reflection. Verified against a published
# door dose-rate reference.
# TODO: source full NCRP 0.5 MeV table.
ALPHA1_PATIENT_SCATTER = 0.0205


def lookup_albedo_soft(table: dict, scatter_angle: float) -> float:
    """
    Albedo for degraded-spectrum radiation. Energy fixed at the soft stand-in,
    not the beam energy. Tables without a Co-60 row fall back to their softest
    tabulated energy, which is still the least wrong number available.
    """
    energy = _SOFT_SPECTRUM_MV if _SOFT_SPECTRUM_MV in table else min(table.keys())
    return lookup_albedo(table, scatter_angle, energy)


# ==============================================================================
# IPEM IDR COMPLIANCE CHECK
# ------------------------------------------------------------------------------
# Two IPEM limits at the door, both in Sv/h:
#   IDR_total <= 7.5e-6   (instantaneous, through the shielding)
#   R8h       <= 0.5e-6   (time-averaged over an 8 h day, 40 h week)
# Rw = unshielded total dose rate x weekly beam-on hours x U, with U = 1 at a
# secondary barrier. R8h = Rw / 40 -- confirmed against three independent
# reference values, so resist the urge to rederive it.
# ==============================================================================

_IPEM_IDR_LIMIT = 7.5e-6   # Sv/h
_IPEM_R8H_LIMIT = 0.5e-6   # Sv/h


def ipem_idr_checks(idr_total: float, dr_total_unshielded: float,
                    weekly_beam_on_hours: float) -> dict:
    """Both IPEM door checks. Returns values, pass flags, and margins."""
    U = 1.0  # secondary barrier, not a choice
    Rw = dr_total_unshielded * weekly_beam_on_hours * U
    R8h = Rw / 40.0
    return {
        "IDR_check_value":  idr_total,
        "IDR_check_pass":   idr_total <= _IPEM_IDR_LIMIT,
        "IDR_check_margin": _IPEM_IDR_LIMIT - idr_total,
        "Rw":               Rw,
        "R8h":              R8h,
        "R8h_check_pass":   R8h <= _IPEM_R8H_LIMIT,
        "R8h_check_margin": _IPEM_R8H_LIMIT - R8h,
    }


# ==============================================================================
# GUI HELPER CLASSES
# Reused across steps. Everything else stays inline per step.
# ==============================================================================

class _Tooltip:
    """Hover tooltip. Shows explanatory label on mouse-enter, removes it on mouse-leave."""

    def __init__(self, widget: tk.Widget, text: str) -> None:
        self._widget = widget
        self._text = text
        self._tw: Optional[tk.Toplevel] = None
        # Bind to enter/leave so the tooltip appears and vanishes on hover.
        widget.bind("<Enter>", self._show)
        widget.bind("<Leave>", self._hide)

    def _show(self, _evt: tk.Event) -> None:
        # Position the tooltip just below and slightly right of the widget.
        # Absolute screen coordinates - Toplevel ignores frame hierarchy.
        x = self._widget.winfo_rootx() + 20
        y = self._widget.winfo_rooty() + self._widget.winfo_height() + 4
        self._tw = tw = tk.Toplevel(self._widget)
        tw.wm_overrideredirect(True)   # strip the window border and title bar
        tw.wm_geometry(f"+{x}+{y}")
        tk.Label(
            tw, text=self._text, justify="left",
            background="#ffffe0", relief="solid", borderwidth=1,
            font=("Segoe UI", 9), padx=6, pady=4,
        ).pack()

    def _hide(self, _evt: tk.Event) -> None:
        # No tooltip on screen? Nothing to do. Otherwise, send it to the void.
        if self._tw:
            self._tw.destroy()
            self._tw = None


class _ValidatedEntry(ttk.Entry):
    """
    Float entry with bounds validation on focus-out.
    Goes red with inline reason label on invalid input. No popups.
    Step validate() sweeps fields and blocks advancement if any are red.
    minv/maxv are inclusive. None means unbounded on that side.
    """

    def __init__(self, parent, minv: Optional[float] = None,
                 maxv: Optional[float] = None, label: str = "Value",
                 reason_widget: Optional[ttk.Label] = None, **kw) -> None:
        super().__init__(parent, **kw)
        self._minv = minv
        self._maxv = maxv
        self._label = label
        self._reason = reason_widget
        # Optimistic until proven otherwise. An untouched default is not "wrong",
        # it is merely "not yet right".
        self.is_valid = True
        self.bind("<FocusOut>", self._check)

    def _expected_text(self) -> str:
        """Describe the acceptable range in words a tired person can understand."""
        if self._minv is not None and self._maxv is not None:
            return f"{self._label} must be a number between {self._minv} and {self._maxv}."
        if self._minv is not None:
            return f"{self._label} must be a number >= {self._minv}."
        if self._maxv is not None:
            return f"{self._label} must be a number <= {self._maxv}."
        return f"{self._label} must be a number."

    def _check(self, _evt: tk.Event = None) -> bool:
        """
        Validate current contents. Returns True/False and sets self.is_valid.
        The step may call this directly during its sweep, or rely on focus-out.
        """
        try:
            val = float(self.get().strip())
        except ValueError:
            self._mark_invalid()
            return False

        if self._minv is not None and val < self._minv:
            self._mark_invalid()
            return False
        if self._maxv is not None and val > self._maxv:
            self._mark_invalid()
            return False

        self._mark_valid()
        return True

    def _mark_valid(self) -> None:
        self.is_valid = True
        self.configure(foreground="")          # default colour
        if self._reason is not None:
            self._reason.config(text="")        # clear the complaint

    def _mark_invalid(self) -> None:
        self.is_valid = False
        self.configure(foreground="red")        # the universal colour of "no"
        if self._reason is not None:
            self._reason.config(text=self._expected_text())  # explain yourself


INTRO_TEXT = (
    "Radiation Shielding Calculation Wizard\n\n"
    "Step-by-step shielding calculations for radiotherapy bunkers and diagnostic "
    "radiology rooms, following NCRP 151 and BIR/Sutton respectively.\n\n"
    "This is a verification and support tool. It carries none of the professional "
    "responsibility. That cross is yours to carry.\n\n"
    "Before you start:\n\n"
    "  -  Have your floor plans, workload data, and design goals ready.\n"
    "  -  Units are stated on every field. They are not suggestions.\n"
    "  -  Red fields block advancement and tell you why they are unhappy.\n\n"
    "Choose your branch below.\n"
)

# ==============================================================================
# WIZARD STEP BASE CLASS
# ==============================================================================

class WizardStep:
    """Base class for all wizard steps."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        self.parent = parent
        self.wizard = wizard
        self.frame: Optional[tk.Frame] = None
        self._data_keys_added: list[str] = []

    def show(self) -> None:
        if self.frame:
            self.frame.pack(fill="both", expand=True, padx=20, pady=20)

    def hide(self) -> None:
        if self.frame:
            self.frame.pack_forget()

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        return {}

    def apply_data(self, wizard_data: dict) -> None:
        data = self.get_data()
        self._data_keys_added = list(data.keys())
        wizard_data.update(data)

    def rollback_data(self, wizard_data: dict) -> None:
        # Undo data contributions when navigating back.
        for key in self._data_keys_added:
            wizard_data.pop(key, None)
        self._data_keys_added.clear()


# ==============================================================================
# STEP: BRANCH SELECTOR
# ==============================================================================

class Step_BranchSelector(WizardStep):
    """First screen. Intro text and branch choice. No inputs, no physics."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=20, pady=(20, 10))

        text = tk.Text(
            text_frame, wrap="word", height=22,
            relief="flat", background="#f5f5f5",
            font=("Segoe UI", 10), padx=12, pady=12,
        )
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        text.insert("1.0", INTRO_TEXT)
        text.configure(state="disabled")
        text.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=(10, 20))

        ttk.Label(
            btn_frame,
            text="What are you shielding today?",
            font=("Segoe UI", 10, "italic"),
        ).pack(pady=(0, 10))

        ttk.Button(
            btn_frame, text="Radiotherapy - Linac Bunker",
            width=35, command=wizard.start_radiotherapy_branch,
        ).pack(pady=4)

        ttk.Button(
            btn_frame, text="Diagnostic Radiology - CT Room",
            width=35, command=wizard.start_radiology_branch,
        ).pack(pady=4)


# ==============================================================================
# STEP: RT CHECKLIST
# ==============================================================================

RT_CHECKLIST_TEXT = (
    "Radiotherapy Branch - Before You Begin\n\n"
    "This branch covers linac bunker shielding following NCRP 151. It will walk "
    "you through workload estimation, IMRT and TBI corrections, maze geometry, "
    "and door shielding in that order.\n\n"
    "Have the following ready:\n\n"
    "  -  Beam energies, dose rate at isocentre Ddot0 [Gy/h], and workload "
    "estimates. If you do not have Ddot0, call the vendor.\n"
    "  -  Bunker floor plan with dimensions. You will need several distances.\n"
    "  -  Design goal P [Sv/week] and occupancy factors. Typical values are "
    "pre-filled but you are expected to verify them for your specific situation.\n"
    "  -  If IMRT is in use: MU ratios. If TBI is in use: treatment distances "
    "and weekly doses.\n\n"
    "The workflow branches based on your answers. Not every screen applies to "
    "every bunker.\n"
)


class Step_RT_Checklist(WizardStep):
    """RT branch orientation and legal acknowledgement gate."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=20, pady=(20, 10))

        text = tk.Text(
            text_frame, wrap="word", height=18,
            relief="flat", background="#f5f5f5",
            font=("Segoe UI", 10), padx=12, pady=12,
        )
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=text.yview)
        text.configure(yscrollcommand=scrollbar.set)
        text.insert("1.0", RT_CHECKLIST_TEXT)
        text.configure(state="disabled")
        text.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        self._ack = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            self.frame,
            text=(
                "I understand that this tool is a verification aid only. "
                "All calculations and their consequences remain the sole "
                "responsibility of the qualified person signing off the work."
            ),
            variable=self._ack,
        ).pack(anchor="w", padx=20, pady=(5, 15))

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=(0, 20))

        ttk.Button(btn_frame, text="Next", command=self._try_advance).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=wizard.prev_step).pack(side="left", padx=5)

    def _try_advance(self) -> None:
        # The checkbox is the only gate. No tick, no progress.
        if not self._ack.get():
            messagebox.showwarning(
                "Acknowledgement Required",
                "Please acknowledge the statement before proceeding."
            )
            return
        self.wizard.next_step()

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        return {}


# ==============================================================================
# SETUP STEP 1 : WORKLOAD PARAMETERS
# ------------------------------------------------------------------------------
# Collects machine workload data for both energies.
# W_base = therapies/day x Gy/therapy x days/week  [Gy/week]
# IMRT and TBI corrections applied in subsequent steps.
#
# Energy dropdown is limited to the energies present in EVERY lookup table
# used downstream (B.4 scatter fractions and all six B.8 albedo tables).
# That intersection is 6, 10 and 18 MV. Anything else would crash a lookup
# somewhere in the maze branch. 15 MV in particular is not tabulated in
# NCRP 151, so it is not offered. Better a missing option than a wrong number.
# ==============================================================================

_ENERGY_OPTIONS = ["6", "10", "18"]
_ENERGY_VALUES  = [6, 10, 18]


class Step_SETUP_Workload(WizardStep):
    """Workload parameters for both beam energies. Calculates W_base inline."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="Workload Parameters",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        ttk.Label(
            self.frame,
            text="Enter data for both beam energies in use.",
            font=("Segoe UI", 9),
        ).pack(pady=(0, 10))

        # Single-energy machines occupy the E_max slot so downstream branching
        # on E_max survives. The E_min slot is hidden and its workload zeroed.
        self._single_energy = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            self.frame, text="Single energy machine?",
            variable=self._single_energy, command=self._toggle_single,
        ).pack(anchor="w", padx=20, pady=(0, 5))

        cols = ttk.Frame(self.frame)
        cols.pack(fill="x", padx=20)

        self._col_min = ttk.LabelFrame(cols, text="Lower Energy (E_min)", padding=10)
        self._col_min.pack(side="left", fill="both", expand=True, padx=(0, 5))

        self._col_max = ttk.LabelFrame(cols, text="Higher Energy (E_max)", padding=10)
        self._col_max.pack(side="left", fill="both", expand=True, padx=(5, 0))

        self._entries = {}
        self._reasons = {}
        self._combos  = {}

        self._build_energy_column(self._col_min, suffix="min")
        self._build_energy_column(self._col_max, suffix="max")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def _build_energy_column(self, parent: ttk.LabelFrame, suffix: str) -> None:

        def row(label_text: str, key: str, minv=None, maxv=None,
                default=None, tooltip: str = ""):
            r = ttk.Frame(parent)
            r.pack(fill="x", pady=2)
            ttk.Label(r, text=label_text, width=22, anchor="w").pack(side="left")
            reason = ttk.Label(r, text="", foreground="red", font=("Segoe UI", 8))
            entry = _ValidatedEntry(
                r, minv=minv, maxv=maxv,
                label=label_text.rstrip(":"),
                reason_widget=reason, width=10,
            )
            if default is not None:
                entry.insert(0, str(default))
            else:
                # Mandatory blank - starts invalid, user cannot skip it.
                entry.is_valid = False
            entry.pack(side="left", padx=4)
            reason.pack(side="left")
            if tooltip:
                _Tooltip(entry, tooltip)
            self._entries[key] = entry
            self._reasons[key] = reason

        e_row = ttk.Frame(parent)
        e_row.pack(fill="x", pady=2)
        ttk.Label(e_row, text="Beam energy [MV]:", width=22, anchor="w").pack(side="left")
        combo = ttk.Combobox(
            e_row, values=_ENERGY_OPTIONS,
            state="readonly", width=12,
        )
        combo.current(0)
        combo.pack(side="left", padx=4)
        _Tooltip(combo, "Beam energy. Only energies tabulated in every NCRP 151 table used here are offered.")
        self._combos[f"E_{suffix}"] = combo

        row("Therapies/day:",  f"therapies_per_day_{suffix}", minv=0.0, default=30,
            tooltip="Number of patient treatments per day at this energy.")
        row("Days/week:",      f"days_per_week_{suffix}",     minv=0.0, maxv=7.0, default=5,
            tooltip="Treatment days per week [1-7].")
        row("Gy/therapy:",     f"Grays_{suffix}",             minv=0.0, default=2.0,
            tooltip="Absorbed dose per treatment session [Gy].")
        row("Ddot0 [Gy/h]:",   f"Ddot0_{suffix}",             minv=0.0,
            tooltip="Instantaneous dose rate at isocentre [Gy/h]. From vendor specs.")
        row("Beam-on [min/tx]:", f"beam_on_{suffix}",         minv=0.0, default=2.0,
            tooltip="Average beam-on time per treatment [minutes]. Feeds the IPEM time-averaged dose rate check.")

    def _toggle_single(self) -> None:
        # Hidden fields keep their defaults; get_data zeroes the slot anyway.
        if self._single_energy.get():
            self._col_min.pack_forget()
        else:
            self._col_min.pack(side="left", fill="both", expand=True,
                               padx=(0, 5), before=self._col_max)

    def validate(self) -> bool:
        # Force-check all entries including untouched ones.
        # A hidden E_min column is exempt - the user cannot fix what they cannot see.
        skip_min = self._single_energy.get()
        return all(
            entry._check() for key, entry in self._entries.items()
            if not (skip_min and key.endswith("_min"))
        )

    def get_data(self) -> dict:
        single = self._single_energy.get()
        data = {"single_energy": single}
        for key, entry in self._entries.items():
            if single and key.endswith("_min"):
                # Unused slot zeroed, not forced equal. E_min stays a valid
                # table energy so downstream lookups do not explode on it.
                data[key] = 0.0
            else:
                data[key] = float(entry.get().strip())
        for key, combo in self._combos.items():
            data[key] = _ENERGY_VALUES[combo.current()]

        # W_base: weekly photon workload before IMRT and TBI corrections.
        # Simple product - if this is wrong, check the inputs, not this line.
        data["W_base_min"] = (
            data["therapies_per_day_min"] *
            data["Grays_min"] *
            data["days_per_week_min"]
        )
        data["W_base_max"] = (
            data["therapies_per_day_max"] *
            data["Grays_max"] *
            data["days_per_week_max"]
        )

        # IPEM time-averaged check needs the weekly beam-on hours.
        data["weekly_beam_on_hours"] = (
            data["beam_on_min"] * data["therapies_per_day_min"] * data["days_per_week_min"]
            + data["beam_on_max"] * data["therapies_per_day_max"] * data["days_per_week_max"]
        ) / 60.0
        return data


# ==============================================================================
# SETUP STEP 2 : IMRT CORRECTION
# ------------------------------------------------------------------------------
# If IMRT is not used, all IMRT values are zeroed and W remains unchanged.
# If IMRT is used, MU ratios and fractions are collected here and C_I,
# F_imrt, W_conv, W_IMRT are calculated inline on advancement.
# ==============================================================================

class Step_SETUP_IMRT(WizardStep):
    """IMRT yes/no with collapsible input fields. Physics inline in get_data."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="IMRT Correction",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        self._imrt_used = tk.BooleanVar(value=False)

        radio_frame = ttk.Frame(self.frame)
        radio_frame.pack(anchor="w", padx=20, pady=5)
        ttk.Label(radio_frame, text="Is IMRT used?").pack(side="left", padx=(0, 15))
        ttk.Radiobutton(
            radio_frame, text="Yes", variable=self._imrt_used,
            value=True, command=self._toggle,
        ).pack(side="left", padx=5)
        ttk.Radiobutton(
            radio_frame, text="No", variable=self._imrt_used,
            value=False, command=self._toggle,
        ).pack(side="left", padx=5)

        # Collapsible input section - hidden until user selects Yes.
        self._imrt_frame = ttk.LabelFrame(
            self.frame, text="IMRT Parameters", padding=10
        )

        self._entries = {}
        self._reasons = {}

        cols = ttk.Frame(self._imrt_frame)
        cols.pack(fill="x")

        col_min = ttk.LabelFrame(cols, text="E_min", padding=8)
        col_min.pack(side="left", fill="both", expand=True, padx=(0, 5))

        col_max = ttk.LabelFrame(cols, text="E_max", padding=8)
        col_max.pack(side="left", fill="both", expand=True, padx=(5, 0))

        self._build_imrt_column(col_min, suffix="min")
        self._build_imrt_column(col_max, suffix="max")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def _build_imrt_column(self, parent: ttk.LabelFrame, suffix: str) -> None:

        def row(label_text: str, key: str, minv=None, maxv=None,
                default=0.0, tooltip: str = ""):
            r = ttk.Frame(parent)
            r.pack(fill="x", pady=2)
            ttk.Label(r, text=label_text, width=20, anchor="w").pack(side="left")
            reason = ttk.Label(r, text="", foreground="red", font=("Segoe UI", 8))
            entry = _ValidatedEntry(
                r, minv=minv, maxv=maxv,
                label=label_text.rstrip(":"),
                reason_widget=reason, width=10,
            )
            entry.insert(0, str(default))
            entry.pack(side="left", padx=4)
            reason.pack(side="left")
            if tooltip:
                _Tooltip(entry, tooltip)
            self._entries[key] = entry
            self._reasons[key] = reason

        row("IMRT fraction p:",  f"p_imrt_{suffix}",   minv=0.0, maxv=1.0,
            tooltip="Fraction of treatments delivered as IMRT [0-1].")
        row("MU_IMRT:",          f"MU_IMRT_{suffix}",  minv=0.0,
            tooltip="Monitor units per IMRT treatment.")
        row("MU_conv:",  f"MU_conv_{suffix}",  minv=0.001,
            tooltip="Monitor units per conventional treatment. Must be > 0.")

    def _toggle(self) -> None:
        # Show or hide the IMRT input frame based on radio selection.
        if self._imrt_used.get():
            self._imrt_frame.pack(fill="x", padx=20, pady=(0, 10))
        else:
            self._imrt_frame.pack_forget()

    def validate(self) -> bool:
        # If IMRT not used, nothing to validate.
        if not self._imrt_used.get():
            return True
        return all(entry._check() for entry in self._entries.values())

    def get_data(self) -> dict:
        if not self._imrt_used.get():
            # IMRT not used -- zero everything, W unchanged.
            return {
                "imrt_used":    False,
                "p_imrt_min":   0.0, "p_imrt_max":   0.0,
                "MU_IMRT_min":  0.0, "MU_IMRT_max":  0.0,
                "MU_conv_min":  0.0, "MU_conv_max":  0.0,
                # C_I is a MU ratio. No IMRT means no MU inflation, ratio of 1.
                # Numerically inert here (multiplied by zero workload) but a
                # ratio of nought would be a lie.
                "C_I_min":      1.0, "C_I_max":      1.0,
                "F_imrt_min":   1.0, "F_imrt_max":   1.0,
                "W_conv_min":   self.wizard.data["W_base_min"],
                "W_conv_max":   self.wizard.data["W_base_max"],
                "W_IMRT_min":   0.0, "W_IMRT_max":   0.0,
            }

        p_min  = float(self._entries["p_imrt_min"].get())
        p_max  = float(self._entries["p_imrt_max"].get())
        mu_i_min = float(self._entries["MU_IMRT_min"].get())
        mu_i_max = float(self._entries["MU_IMRT_max"].get())
        mu_c_min = float(self._entries["MU_conv_min"].get())
        mu_c_max = float(self._entries["MU_conv_max"].get())

        # C_I: IMRT MU ratio. Guards against zero conventional MU.
        C_I_min = mu_i_min / mu_c_min if mu_c_min > 0.0 else 0.0
        C_I_max = mu_i_max / mu_c_max if mu_c_max > 0.0 else 0.0

        # F_imrt: effective workload multiplier accounting for IMRT fraction.
        F_imrt_min = p_min * C_I_min + (1.0 - p_min)
        F_imrt_max = p_max * C_I_max + (1.0 - p_max)

        W_base_min = self.wizard.data["W_base_min"]
        W_base_max = self.wizard.data["W_base_max"]

        return {
            "imrt_used":   True,
            "p_imrt_min":  p_min,   "p_imrt_max":  p_max,
            "MU_IMRT_min": mu_i_min, "MU_IMRT_max": mu_i_max,
            "MU_conv_min": mu_c_min, "MU_conv_max": mu_c_max,
            "C_I_min":     C_I_min,  "C_I_max":     C_I_max,
            "F_imrt_min":  F_imrt_min, "F_imrt_max": F_imrt_max,
            "W_conv_min":  (1.0 - p_min) * W_base_min,
            "W_conv_max":  (1.0 - p_max) * W_base_max,
            "W_IMRT_min":  p_min * W_base_min,
            "W_IMRT_max":  p_max * W_base_max,
        }


# ==============================================================================
# SETUP STEP 3 : TBI CORRECTION
# ------------------------------------------------------------------------------
# If TBI is not used, W_tbi is zeroed and W is unchanged.
# If TBI is used, d_tbi and weekly doses are collected here.
# W_tbi = D_tbi * d_tbi^2  [Gy/week]
# Final workload assembly (W, W_L) happens in Step_SETUP_Maze after TBI.
# ==============================================================================

class Step_SETUP_TBI(WizardStep):
    """TBI yes/no with collapsible input fields. Physics inline in get_data."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="TBI Correction",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        self._tbi_used = tk.BooleanVar(value=False)

        radio_frame = ttk.Frame(self.frame)
        radio_frame.pack(anchor="w", padx=20, pady=5)
        ttk.Label(radio_frame, text="Is TBI used?").pack(side="left", padx=(0, 15))
        ttk.Radiobutton(
            radio_frame, text="Yes", variable=self._tbi_used,
            value=True, command=self._toggle,
        ).pack(side="left", padx=5)
        ttk.Radiobutton(
            radio_frame, text="No", variable=self._tbi_used,
            value=False, command=self._toggle,
        ).pack(side="left", padx=5)

        # Collapsible input section - hidden until user selects Yes.
        self._tbi_frame = ttk.LabelFrame(
            self.frame, text="TBI Parameters", padding=10
        )

        self._entries = {}
        self._reasons = {}

        def row(label_text: str, key: str, minv=None, maxv=None,
                default=None, tooltip: str = ""):
            r = ttk.Frame(self._tbi_frame)
            r.pack(fill="x", pady=2)
            ttk.Label(r, text=label_text, width=28, anchor="w").pack(side="left")
            reason = ttk.Label(r, text="", foreground="red", font=("Segoe UI", 8))
            entry = _ValidatedEntry(
                r, minv=minv, maxv=maxv,
                label=label_text.rstrip(":"),
                reason_widget=reason, width=10,
            )
            if default is not None:
                entry.insert(0, str(default))
            else:
                entry.is_valid = False
            entry.pack(side="left", padx=4)
            reason.pack(side="left")
            if tooltip:
                _Tooltip(entry, tooltip)
            self._entries[key] = entry
            self._reasons[key] = reason

        row("d_tbi [m]:",          "d_tbi",     minv=0.001,
            tooltip="Treatment distance for TBI, isocentre to patient midplane [m].")
        row("D_tbi_min [Gy/week]:", "D_tbi_min", minv=0.0,
            tooltip="Weekly TBI dose at E_min [Gy/week].")
        row("D_tbi_max [Gy/week]:", "D_tbi_max", minv=0.0,
            tooltip="Weekly TBI dose at E_max [Gy/week].")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def _toggle(self) -> None:
        if self._tbi_used.get():
            self._tbi_frame.pack(fill="x", padx=20, pady=(0, 10))
        else:
            self._tbi_frame.pack_forget()

    def validate(self) -> bool:
        if not self._tbi_used.get():
            return True
        return all(entry._check() for entry in self._entries.values())

    def get_data(self) -> dict:
        if not self._tbi_used.get():
            return {
                "tbi_used":  False,
                "d_tbi":     0.0,
                "D_tbi_min": 0.0, "D_tbi_max": 0.0,
                "W_tbi_min": 0.0, "W_tbi_max": 0.0,
            }

        d_tbi     = float(self._entries["d_tbi"].get())
        D_tbi_min = float(self._entries["D_tbi_min"].get())
        D_tbi_max = float(self._entries["D_tbi_max"].get())

        # W_tbi accounts for the extended SSD used in TBI.
        # Inverse square correction relative to isocentre distance.
        W_tbi_min = D_tbi_min * (d_tbi ** 2)
        W_tbi_max = D_tbi_max * (d_tbi ** 2)

        return {
            "tbi_used":  True,
            "d_tbi":     d_tbi,
            "D_tbi_min": D_tbi_min, "D_tbi_max": D_tbi_max,
            "W_tbi_min": W_tbi_min, "W_tbi_max": W_tbi_max,
        }


# ==============================================================================
# SETUP STEP 4 : MAZE AND FINAL WORKLOAD ASSEMBLY
# ------------------------------------------------------------------------------
# Last SETUP step before the workflow branches.
# Assembles final W and W_L from all previous corrections.
#
# W_L = W_conv + W_tbi + C_I * W_IMRT   (leakage workload, used in L/N branch)
# W   = W_conv + F_imrt * W_IMRT         (scatter workload, used in M/D branch)
# ==============================================================================

class Step_SETUP_Maze(WizardStep):
    """Maze yes/no. Assembles final W and W_L inline before branching."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="Bunker Configuration",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        self._maze_exists = tk.BooleanVar(value=True)

        radio_frame = ttk.Frame(self.frame)
        radio_frame.pack(anchor="w", padx=20, pady=10)
        ttk.Label(radio_frame, text="Does the bunker have a maze?").pack(side="left", padx=(0, 15))
        ttk.Radiobutton(
            radio_frame, text="Yes", variable=self._maze_exists,
            value=True,
        ).pack(side="left", padx=5)
        ttk.Radiobutton(
            radio_frame, text="No", variable=self._maze_exists,
            value=False,
        ).pack(side="left", padx=5)

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        d = self.wizard.data

        # Final workload assembly - all corrections now applied.
        # W_L used by leakage-dominated paths (head leakage, neutron branch).
        # W used by scatter-dominated paths (maze photon, direct door).
        W_L_min = d["W_conv_min"] + d["W_tbi_min"] + d["C_I_min"] * d["W_IMRT_min"]
        W_L_max = d["W_conv_max"] + d["W_tbi_max"] + d["C_I_max"] * d["W_IMRT_max"]
        W_min   = d["W_conv_min"] + d["F_imrt_min"] * d["W_IMRT_min"]
        W_max   = d["W_conv_max"] + d["F_imrt_max"] * d["W_IMRT_max"]

        return {
            "maze_exists": self._maze_exists.get(),
            "W_L_min": W_L_min,
            "W_L_max": W_L_max,
            "W_min":   W_min,
            "W_max":   W_max,
        }


# ==============================================================================
# M STEP 1 : PRIMARY BEAM CHECK
# ------------------------------------------------------------------------------
# Primary beam at the maze entrance is a design failure, not a calculation.
# Three strikes and the programme leaves. So does our patience.
# ==============================================================================

class Step_M_Primary(WizardStep):
    """Yes/no: does the primary beam reach the maze entrance. Escalates on repeat 'Yes'."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        # How many times the user has insisted 'Yes'. Patience is finite.
        self._yes_count = 0

        ttk.Label(
            self.frame, text="Primary Beam at Maze Entrance",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        self._primary_reaches = tk.BooleanVar(value=False)

        radio_frame = ttk.Frame(self.frame)
        radio_frame.pack(anchor="w", padx=20, pady=10)
        ttk.Label(
            radio_frame,
            text="Does the primary beam reach the maze entrance?",
        ).pack(side="left", padx=(0, 15))
        ttk.Radiobutton(
            radio_frame, text="Yes", variable=self._primary_reaches,
            value=True,
        ).pack(side="left", padx=5)
        ttk.Radiobutton(
            radio_frame, text="No", variable=self._primary_reaches,
            value=False,
        ).pack(side="left", padx=5)

        # Inline reason label, sits under the radios. Grows ruder on repetition.
        self._reason = ttk.Label(
            self.frame, text="", foreground="red",
            font=("Segoe UI", 10), wraplength=500, justify="left",
        )
        self._reason.pack(anchor="w", padx=20, pady=(0, 10))

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self._try_advance).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def _try_advance(self) -> None:
        if not self._primary_reaches.get():
            self._yes_count = 0
            self._reason.config(text="")
            self.wizard.next_step()
            return

        self._yes_count += 1

        if self._yes_count == 1:
            self._reason.config(text=(
                "The primary beam should never reach the maze entrance. "
                "Check your geometry. This is a design problem, not an input."
            ))
        elif self._yes_count == 2:
            self._reason.config(text=(
                "It still should not. If your beam reaches the maze entrance, "
                "the maze is in the wrong place. Move the wall, not the goalposts."
            ))
        else:
            # Third strike. The programme and I are both leaving.
            messagebox.showerror(
                "Enough",
                "You need Jesus and I want no part in this."
            )
            self.wizard.exit_wizard()

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        return {"primary_reaches_maze": self._primary_reaches.get()}


# ==============================================================================
# ALBEDO TABLE MAP
# Maps (material, incident_angle) to the relevant table. Replaces globals().
# ==============================================================================

_ALBEDO_TABLE_MAP = {
    ("concrete", 0):  ALBEDO_CONCRETE_0_TABLE,
    ("concrete", 45): ALBEDO_CONCRETE_45_TABLE,
    ("iron", 0):      ALBEDO_IRON_0_TABLE,
    ("iron", 45):     ALBEDO_IRON_45_TABLE,
    ("lead", 0):      ALBEDO_LEAD_0_TABLE,
    ("lead", 45):     ALBEDO_LEAD_45_TABLE,
}


# ==============================================================================
# M STEP 2 : MAZE GEOMETRY AND PHOTON PARAMETERS
# ------------------------------------------------------------------------------
# Collects architectural dimensions, distances, angles and albedo selections.
# Derives A0, A1, Az and d_h so the user never has to type an area.
# Every derivation below is an assumption. Verify against your own layout.
# ==============================================================================

class Step_M_Geometry(WizardStep):
    """Maze geometry inputs and albedo selections."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="Maze Geometry and Photon Parameters",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        # Scrollable container. The field count does not fit on one screen.
        canvas = tk.Canvas(self.frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.frame, orient="vertical", command=canvas.yview)
        scroll_content = ttk.Frame(canvas)

        scroll_content.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=scroll_content, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        self._entries = {}
        self._reasons = {}
        self._combos = {}

        # --- Design and beam parameters
        p_frame = ttk.LabelFrame(scroll_content, text="Design and Beam Parameters", padding=10)
        p_frame.pack(fill="x", padx=10, pady=5)

        self._add_row(p_frame, "Design goal P [Sv/week]:", "P", minv=1e-6, default="0.00006",
                      tooltip="Weekly dose limit at the external boundary.")
        self._add_row(p_frame, "Use factor U_G:", "U_G", minv=0.0, maxv=1.0, default="1.0",
                      tooltip="Fraction of workload directed toward the maze wall.")
        self._add_row(p_frame, "Head leakage ratio L_f:", "L_f", minv=0.0, maxv=0.01, default="0.001",
                      tooltip="Fraction of primary beam escaping the head. Typically 0.001.")
        self._add_row(p_frame, "Patient transmission f:", "f", minv=0.0, maxv=1.0, default="0.0",
                      tooltip="Fraction of primary beam transmitted through the patient.")
        self._add_row(p_frame, "Field size F [cm2]:", "F", minv=1.0, default="400.0",
                      tooltip="Field size at isocentre in cm2. NCRP normalises to 400 cm2.")

        # --- Architectural dimensions
        g_frame = ttk.LabelFrame(scroll_content, text="Architectural Dimensions", padding=10)
        g_frame.pack(fill="x", padx=10, pady=5)

        self._add_row(g_frame, "Maze ceiling height h_e [m]:", "h_e", minv=0.01, default="2.0",
                      tooltip="Height of the maze corridor ceiling.")
        self._add_row(g_frame, "Maze corridor width w_1 [m]:", "w_1", minv=0.01, default="1.7",
                      tooltip="Width of the main maze passageway.")
        self._add_row(g_frame, "Maze tongue width w_tongue [m]:", "w_tongue", minv=0.01, default="2.6",
                      tooltip="Width of the inner baffle wall facing the maze opening.")
        self._add_row(g_frame, "Inner wall thickness t [m]:", "t", minv=0.01, default="1.5",
                      tooltip="Physical thickness of the concrete maze baffle.")

        # --- Distances and angles
        d_frame = ttk.LabelFrame(scroll_content, text="Distances and Angles", padding=10)
        d_frame.pack(fill="x", padx=10, pady=5)

        self._add_row(d_frame, "Isocentre to opposing wall d_pp [m]:", "d_pp", minv=0.01,
                      tooltip="Perpendicular distance from isocentre to the primary scattering wall.")
        self._add_row(d_frame, "Primary wall to maze d_r [m]:", "d_r", minv=0.01,
                      tooltip="Scatter path from the primary wall to the maze inner wall.")
        self._add_row(d_frame, "Target to maze entrance d_sec [m]:", "d_sec", minv=0.01,
                      tooltip="Distance from isocentre to the inner maze wall.")
        self._add_row(d_frame, "Patient to maze entrance d_sca [m]:", "d_sca", minv=0.01, default="1.0",
                      tooltip="Patient scatter distance. NCRP standard is 1.0 m.")
        self._add_row(d_frame, "Target to door (direct) d_L [m]:", "d_L", minv=0.01,
                      tooltip="Direct line of sight from isocentre to door.")
        self._add_row(d_frame, "Maze centreline to door d_z [m]:", "d_z", minv=0.01,
                      tooltip="Distance along the maze centreline to the doorway.")
        self._add_row(d_frame, "Total secondary path d_zz [m]:", "d_zz", minv=0.01,
                      tooltip="Total distance along the maze path to the door.")
        self._add_row(d_frame, "Incident angle theta_i [deg]:", "i_theta", minv=0.0, maxv=89.0, default="0.0",
                      tooltip="Angle of primary beam incidence on the opposing wall.")
        self._add_row(d_frame, "Scatter angle theta_s [deg]:", "s_theta", minv=10.0, maxv=150.0, default="45.0",
                      tooltip="Patient scatter angle towards the maze entrance. Table B.4 range is 10 to 150 deg.")

        # --- Albedo selections
        a_frame = ttk.LabelFrame(scroll_content, text="Reflection and Albedo Parameters", padding=10)
        a_frame.pack(fill="x", padx=10, pady=5)

        _soft_note = (
            "Second-reflection radiation is soft (~0.5 MeV) regardless of beam "
            "energy, so this lookup uses a fixed energy. True 0.5 MeV data is "
            "not yet digitised; the Co-60 row (or softest tabulated energy) "
            "stands in for now."
        )
        self._add_combo_row(a_frame, "a_z material:", "alpha_z_material", ["concrete", "iron", "lead"],
                            tooltip=_soft_note)
        self._add_combo_row(a_frame, "a_z incident angle:", "alpha_z_angle", [0, 45],
                            tooltip=_soft_note)
        self._add_combo_row(a_frame, "a_z reflection angle:", "alpha_z_scatter_angle", [0, 30, 45, 60, 75],
                            tooltip=_soft_note)

        self._add_combo_row(a_frame, "a0 material:", "alpha0_material", ["concrete", "iron", "lead"])
        self._add_combo_row(a_frame, "a0 incident angle:", "alpha0_angle", [0, 45])
        self._add_combo_row(a_frame, "a0 reflection angle:", "alpha0_scatter_angle", [0, 30, 45, 60, 75])

        _a1_note = (
            "Used for head-leakage scatter only. Patient scatter now uses a "
            "fixed 0.5 MeV albedo (0.0205, interim value pending the full "
            "NCRP 0.5 MeV table) and ignores these dropdowns."
        )
        self._add_combo_row(a_frame, "a1 material:", "alpha1_material", ["concrete", "iron", "lead"],
                            tooltip=_a1_note)
        self._add_combo_row(a_frame, "a1 incident angle:", "alpha1_angle", [0, 45],
                            tooltip=_a1_note)
        self._add_combo_row(a_frame, "a1 reflection angle:", "alpha1_scatter_angle", [0, 30, 45, 60, 75],
                            tooltip=_a1_note)

        # --- Geometry guide
        s_frame = ttk.LabelFrame(scroll_content, text="Geometry Guide", padding=10)
        s_frame.pack(fill="x", padx=10, pady=5)
        self._draw_photon_schematic(s_frame)

        canvas.pack(side="left", fill="both", expand=True, padx=(10, 0))
        scrollbar.pack(side="right", fill="y", padx=(0, 10))

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=10)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def _add_row(self, parent: ttk.Frame, label_text: str, key: str,
                 minv=None, maxv=None, default="", tooltip="") -> None:
        r = ttk.Frame(parent)
        r.pack(fill="x", pady=2)
        ttk.Label(r, text=label_text, width=34, anchor="w").pack(side="left")
        reason = ttk.Label(r, text="", foreground="red", font=("Segoe UI", 8))
        entry = _ValidatedEntry(
            r, minv=minv, maxv=maxv,
            label=label_text.rstrip(":"),
            reason_widget=reason, width=12,
        )
        if default != "":
            entry.insert(0, str(default))
        else:
            entry.is_valid = False
        entry.pack(side="left", padx=4)
        reason.pack(side="left")
        if tooltip:
            _Tooltip(entry, tooltip)
        self._entries[key] = entry
        self._reasons[key] = reason

    def _add_combo_row(self, parent: ttk.Frame, label_text: str, key: str,
                       values: list, tooltip: str = "") -> None:
        r = ttk.Frame(parent)
        r.pack(fill="x", pady=2)
        ttk.Label(r, text=label_text, width=34, anchor="w").pack(side="left")
        combo = ttk.Combobox(r, values=values, state="readonly", width=10)
        combo.current(0)
        combo.pack(side="left", padx=4)
        if tooltip:
            _Tooltip(combo, tooltip)
        self._combos[key] = combo

    def validate(self) -> bool:
        return all(entry._check() for entry in self._entries.values())

    def get_data(self) -> dict:
        data = {}
        for key, entry in self._entries.items():
            data[key] = float(entry.get().strip())
        for key, combo in self._combos.items():
            val = combo.get()
            data[key] = int(val) if val.isdigit() else val

        # ----------------------------------------------------------------------
        # DERIVED GEOMETRY
        # These four are assumptions, not measurements. Check them against your
        # layout before trusting any number downstream.
        # ----------------------------------------------------------------------

        # Target sits one SAD behind isocentre, so target-to-wall is d_pp + 1.0 m.
        data["d_h"] = data["d_pp"] + 1.0

        # A0: primary field area where it lands on the opposing wall.
        # F is defined at the isocentre plane, 1 m from the target. The wall is
        # a further d_pp beyond the isocentre on the same divergent beam, so the
        # field side scales by d_pp, not d_h. Matches the Agios Andreas working
        # (40 cm x 5.95 / 1 = 238 cm). Scaling by d_h overestimated A0 by ~37%.
        field_side_at_iso = math.sqrt(data["F"]) / 100.0   # cm2 -> m side length
        data["A0"] = (field_side_at_iso * data["d_pp"]) ** 2

        # A1: cross-section of the maze tongue facing the room.
        data["A1"] = data["w_tongue"] * data["h_e"]

        # Az: cross-section of the maze corridor.
        data["Az"] = data["w_1"] * data["h_e"]
        return data

    def _draw_photon_schematic(self, parent: ttk.Frame) -> None:
        """Wall-scatter and leakage paths for the M/N-branch door dose."""
        ttk.Label(parent, text="Photon Geometry", font=("Segoe UI", 9, "bold")).pack(pady=(5, 0))
        cv = tk.Canvas(parent, width=620, height=380, bg="white", highlightthickness=0)
        cv.pack(pady=(0, 10))

        IN_L, IN_T, IN_R, IN_B = 34, 40, 562, 328
        TONGUE_L, TONGUE_R, TONGUE_T = 350, 412, 168

        # Maze passage narrowed to 2/3 of its former width
        OLD_CORR_W = IN_R - TONGUE_R
        CORR_W = round(OLD_CORR_W * 2 / 3)
        CORR_L = TONGUE_R
        CORR_R = TONGUE_R + CORR_W
        CL = (CORR_L + CORR_R) // 2

        # Outer bounding box adjustments
        OUT_L, OUT_T, OUT_B = 10, 10, 370
        # Make the right wall thickness match the top wall thickness (30px)
        WALL_THICKNESS = IN_T - OUT_T
        OUT_R = CORR_R + WALL_THICKNESS

        # Door shifted to the lowest border
        DOOR_X0, DOOR_X1 = CORR_L, CORR_R
        DOOR_Y_BOTTOM = OUT_B
        DOOR_Y_TOP = DOOR_Y_BOTTOM - 16
        DOOR = (CL, DOOR_Y_TOP)

        ISO = (185, 225)
        A0_PT = (245, IN_T)
        A1_PT = (CL, IN_T)

        # Point b at the tongue centreline crossed with Az mid-height
        B_PT = (CL, 101)

        RED, WALL, DARK, BLUE = "#cc0000", "#c4c4c4", "#3a3a3a", "#0044aa"

        def label(x, y, txt, fill="black", size=10):
            t = cv.create_text(x, y, text=txt, fill=fill, font=("Segoe UI", size, "bold"))
            x0, y0, x1, y1 = cv.bbox(t)
            cv.tag_lower(cv.create_rectangle(x0-3, y0-1, x1+3, y1+1, fill="white", outline=""), t)

        cv.create_rectangle(OUT_L, OUT_T, DOOR_X0, OUT_B, fill=WALL, outline="#333", width=2)
        cv.create_rectangle(DOOR_X1, OUT_T, OUT_R, OUT_B, fill=WALL, outline="#333", width=2)
        cv.create_rectangle(OUT_L, OUT_T, OUT_R, IN_T, fill=WALL, outline="#333", width=2)
        cv.create_line(OUT_L, OUT_T, OUT_L, OUT_B, fill="#333", width=2)
        cv.create_line(OUT_R, OUT_T, OUT_R, OUT_B, fill="#333", width=2)
        cv.create_line(OUT_L, OUT_B, DOOR_X0, OUT_B, fill="#333", width=2)
        cv.create_line(DOOR_X1, OUT_B, OUT_R, OUT_B, fill="#333", width=2)

        # Extend maze corridor white polygon down to DOOR_Y_BOTTOM
        cv.create_polygon(
            IN_L, IN_T, CORR_R, IN_T, CORR_R, DOOR_Y_BOTTOM,
            TONGUE_R, DOOR_Y_BOTTOM, TONGUE_R, TONGUE_T,
            TONGUE_L, TONGUE_T, TONGUE_L, IN_B, IN_L, IN_B,
            fill="white", outline="#333", width=2,
        )

        cv.create_rectangle(DOOR_X0, DOOR_Y_TOP, DOOR_X1, DOOR_Y_BOTTOM, fill="#9aa0b0", outline="#333", width=2)
        label(CL, DOOR_Y_BOTTOM - 8, "door")

        cv.create_rectangle(190, IN_T-17, 300, IN_T, fill=DARK, outline="")
        label(245, IN_T-28, "A0")
        cv.create_rectangle(CORR_L, IN_T-17, CORR_R, IN_T, fill=DARK, outline="")
        label(CL, IN_T-28, "A1")

        # Az and its label behind it (inside the right wall)
        cv.create_rectangle(CORR_R, 52, CORR_R+17, 150, fill=DARK, outline="")
        label(CORR_R+16, 101, "Az")

        ix, iy = ISO
        cv.create_line(TONGUE_L, 218, TONGUE_R, 218, fill=BLUE, width=2, arrow=tk.BOTH)
        label((TONGUE_L+TONGUE_R)//2, 206, "t", fill=BLUE)

        cv.create_line(A0_PT[0], IN_T, A0_PT[0], IN_T+70, fill="#666", width=1, dash=(4, 3))
        cv.create_arc(A0_PT[0]-42, IN_T-42, A0_PT[0]+42, IN_T+42,
                      start=252, extent=18, style="arc", outline=BLUE, width=2)
        label(A0_PT[0]-52, IN_T+30, "i_theta", fill=BLUE, size=9)

        cv.create_line(ix, iy, *A0_PT, fill=RED, width=3, arrow=tk.BOTH)
        label(196, 130, "d_h", fill=RED)

        # d_r terminating at point b
        cv.create_line(*A0_PT, *B_PT, fill=RED, width=3, arrow=tk.BOTH)
        label((A0_PT[0] + B_PT[0])//2 - 5, (A0_PT[1] + B_PT[1])//2 - 15, "d_r", fill=RED)

        cv.create_line(ix, iy, *A1_PT, fill=RED, width=3, arrow=tk.BOTH)
        label(310, 148, "d_sec", fill=RED)

        # d_zz extending straight down to the door
        cv.create_line(A1_PT[0], IN_T, CL, DOOR_Y_TOP, fill=RED, width=3, arrow=tk.LAST)
        label(CL+22, 230, "d_zz", fill=RED)

        # d_z going FROM door TO b
        cv.create_line(CL, DOOR_Y_TOP, B_PT[0], B_PT[1], fill=RED, width=2, dash=(6, 3), arrow=tk.LAST)
        label(CL-30, 230, "d_z", fill=RED)

        # d_L extending straight to the door
        cv.create_line(ix, iy, *DOOR, fill=RED, width=3, arrow=tk.BOTH)
        label(305, 305, "d_L", fill=RED)

        cv.create_arc(ix-56, iy-56, ix+56, iy+56, start=72, extent=-38,
                      style="arc", outline=BLUE, width=2)
        label(ix+70, iy-40, "s_theta", fill=BLUE, size=9)

        cv.create_line(ix-9, iy, ix+9, iy, fill="black", width=2)
        cv.create_line(ix, iy-9, ix, iy+9, fill="black", width=2)
        label(ix-34, iy+16, "iso", size=9)
        label(B_PT[0]-14, B_PT[1]-12, "b")


# ==============================================================================
# M STEP 3 : PHOTON DOSE AT THE DOOR
# ------------------------------------------------------------------------------
# NCRP 151 Eq 3.4 to 3.8, then H_ph = 2.64 * H0 and the required transmission.
# Albedo tables already carry x10^-3, so nothing is rescaled here.
#
# Design goal allocation:
#   E_max <  10 MV : photons are the only component at the door. B_ph = P / H_ph.
#   E_max >= 10 MV : the goal is shared with neutrons and capture gamma, each
#                    taking P/2 (Agios Andreas convention). B_ph = P / (2 H_ph).
# ==============================================================================

class Step_M_Calculations(WizardStep):
    """Photon scatter and leakage dose at the maze door."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        # Populated by show(), returned by get_data() so rollback can undo it.
        self._results: dict = {}

        ttk.Label(
            self.frame, text="Photon Dose at the Door",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 10))

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=10, pady=5)

        self.text_area = tk.Text(text_frame, wrap="word", height=20, width=80,
                                 font=("Consolas", 10))
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.text_area.yview)
        self.text_area.configure(yscrollcommand=scrollbar.set)
        self.text_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def show(self) -> None:
        super().show()
        self._execute_calculations()

    def _execute_calculations(self) -> None:
        d = self.wizard.data
        self._results = {}

        try:
            # Direct key access. A missing key means an upstream step failed to
            # store it, and that must surface as an error, not a plausible default.
            W_min, W_max = d["W_min"], d["W_max"]
            W_L_min, W_L_max = d["W_L_min"], d["W_L_max"]
            E_min, E_max = d["E_min"], d["E_max"]
            P, U_G, L_f, f_trans, F = d["P"], d["U_G"], d["L_f"], d["f"], d["F"]

            d_h, d_r, d_z, d_zz = d["d_h"], d["d_r"], d["d_z"], d["d_zz"]
            d_sec, d_sca, d_L = d["d_sec"], d["d_sca"], d["d_L"]
            A0, A1, Az = d["A0"], d["A1"], d["Az"]
            t_wall, i_theta, s_theta = d["t"], d["i_theta"], d["s_theta"]

            # Slant path through the inner maze wall for leakage transmission.
            t_s_cm = (t_wall / math.cos(math.radians(i_theta))) * 100.0
            tvl_ph = 0.6
            B_LT = 10.0 ** (-(t_s_cm / tvl_ph))

            # Patient scatter fraction. Potato's real formula, verified: uses
            # its bilinear interpolation on Table B.4 with the axes swapped,
            # not a clean nearest-neighbour lookup. See lookup_alpha_theta_potato.
            alpha_theta = lookup_alpha_theta_potato(s_theta, E_max)

            # Wall albedo. Tables already carry x10^-3, so no rescaling.
            a_z_tbl = _ALBEDO_TABLE_MAP[(d["alpha_z_material"], d["alpha_z_angle"])]
            a0_tbl = _ALBEDO_TABLE_MAP[(d["alpha0_material"], d["alpha0_angle"])]
            a1_tbl = _ALBEDO_TABLE_MAP[(d["alpha1_material"], d["alpha1_angle"])]

            # Second reflection sees a soft spectrum, not the beam energy.
            alpha_z = lookup_albedo(a_z_tbl, d["alpha_z_scatter_angle"], E_max)
            alpha0_max = lookup_albedo(a0_tbl, d["alpha0_scatter_angle"], E_max)
            alpha0_min = lookup_albedo(a0_tbl, d["alpha0_scatter_angle"], E_min)
            alpha1_max = lookup_albedo(a1_tbl, d["alpha1_scatter_angle"], E_max)
            alpha1_min = lookup_albedo(a1_tbl, d["alpha1_scatter_angle"], E_min)

            # Eq 3.5: scatter off the opposing wall, then off the maze wall.
            H_s = (
                (W_max * alpha0_max + W_min * alpha0_min)
                * U_G * A0 * alpha_z * Az
            ) / ((d_h * d_r * d_z) ** 2)

            # Eq 3.6: head leakage scattered off the maze wall.
            H_LS = (
                L_f * (W_L_max * alpha1_max + W_L_min * alpha1_min) * U_G * A1
            ) / ((d_sec * d_zz) ** 2)

            # Eq 3.7: patient scatter off the maze wall. Potato's real weekly
            # formula (functions.py, weekly_dose_patient_scatter) has NO F/400
            # term - that normalisation only appears in potato's separate
            # instantaneous DR_ps calc (Step_M_Bph), not here. Confirmed:
            # including it here overstated H_PS by exactly 4x (F/400=4)
            # against the verified reference.
            H_PS = (
                alpha_theta * (W_max * alpha1_max + W_min * alpha1_min)
                * U_G * A1
            ) / ((d_sca * d_sec * d_zz) ** 2)

            # Eq 3.8: leakage straight through the inner maze wall.
            H_LT = (L_f * (W_L_max + W_L_min) * U_G * B_LT) / (d_L ** 2)

            # Eq 3.4 then Eq 3.10. The 2.64 covers all gantry orientations.
            H0 = f_trans * H_s + H_LS + H_PS + H_LT
            H_ph = 2.64 * H0

            # Design goal allocation. See module banner.
            if E_max >= 10.0:
                B_ph = P / (2.0 * H_ph)
                goal_note = "P/2 (shared with neutron and capture gamma at E_max >= 10 MV)"
            else:
                B_ph = P / H_ph
                goal_note = "full P (photons are the only component below 10 MV)"

            self._results = {
                "B_LT": B_LT, "alpha_theta": alpha_theta, "alpha_z": alpha_z,
                "alpha0_max": alpha0_max, "alpha0_min": alpha0_min,
                "alpha1_max": alpha1_max, "alpha1_min": alpha1_min,
                "H_s": H_s, "H_LS": H_LS, "H_PS": H_PS, "H_LT": H_LT,
                "H0": H0, "H_ph": H_ph, "B_ph": B_ph,
            }

            report = [
                "=" * 65,
                "NCRP 151 PHOTON MAZE DOOR RESULTS",
                "=" * 65,
                f"Workload W_min / W_max          : {W_min:.1f} / {W_max:.1f} Gy/wk",
                f"Leakage workload W_L_max        : {W_L_max:.1f} Gy/wk",
                f"Beam energies E_min / E_max     : {E_min} / {E_max} MV",
                "-" * 65,
                "Dose equivalent components:",
                f"  Opposing wall scatter  H_s    : {H_s:.3e} Sv/wk",
                f"  Head leakage scatter   H_LS   : {H_LS:.3e} Sv/wk",
                f"  Patient scatter        H_PS   : {H_PS:.3e} Sv/wk",
                f"  Leakage transmission   H_LT   : {H_LT:.3e} Sv/wk",
                "-" * 65,
                f"Unshielded total       H0     : {H0:.3e} Sv/wk",
                f"All orientations       H_ph   : {H_ph:.3e} Sv/wk",
                f"Design goal            P      : {P:.3e} Sv/wk",
                f"Goal applied to photons       : {goal_note}",
                f"Required transmission  B_ph   : {B_ph:.3e}",
                "=" * 65,
            ]
            if E_max >= 10.0:
                report.append("E_max >= 10 MV. Neutron and capture gamma branch follows.")
            else:
                report.append("E_max < 10 MV. Photon-only completion follows.")

            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", "\n".join(report))

        except Exception as e:
            # Genuine calculation failure, not user input. Input is caught upstream.
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", f"Calculation failed:\n\n{type(e).__name__}: {e}")

    def validate(self) -> bool:
        # No results means the calculation failed. Do not let it advance.
        return bool(self._results)

    def get_data(self) -> dict:
        return dict(self._results)


# ==============================================================================
# M STEP 4 : REQUIRED LEAD THICKNESS AND INSTANTANEOUS DOSE RATE
# ------------------------------------------------------------------------------
# B_ph already computed in Step_M_Calculations. This step converts it to a
# lead thickness, then re-derives the dose RATE (not weekly dose) through
# that same thickness as a cross-check using Ddot0 instead of workload.
# TVL for photon scatter/leakage in lead is 0.6 cm throughout.
# ==============================================================================

class Step_M_Bph(WizardStep):
    """Lead thickness from B_ph, then instantaneous dose rate through it."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)
        self._results: dict = {}

        ttk.Label(
            self.frame, text="Required Lead Thickness",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 10))

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=10, pady=5)
        self.text_area = tk.Text(text_frame, wrap="word", height=18, width=80,
                                 font=("Consolas", 10))
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.text_area.yview)
        self.text_area.configure(yscrollcommand=scrollbar.set)
        self.text_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def show(self) -> None:
        super().show()
        self._execute()

    def _execute(self) -> None:
        d = self.wizard.data
        self._results = {}
        tvl_ph = 0.6

        try:
            B_ph = d["B_ph"]

            # Required lead thickness. B_ph > 1 means the unshielded dose
            # already meets the design goal, so no lead is needed.
            if B_ph <= 1.0:
                N_ph = -math.log10(B_ph)
                x_ph_lead = N_ph * tvl_ph
            else:
                x_ph_lead = 0.0

            # Instantaneous dose rate through x_ph_lead, using Ddot0 rather
            # than weekly workload.
            d_h, d_r, d_z = d["d_h"], d["d_r"], d["d_z"]
            d_sca, d_sec, d_zz, d_L = d["d_sca"], d["d_sec"], d["d_zz"], d["d_L"]

            if min(d_h, d_r, d_z, d_sca, d_sec, d_zz, d_L) <= 0.0:
                raise ValueError("All distance parameters must be > 0.")

            # Structurally mirrors H_s: pure double-reflection geometry with
            # Ddot0 in place of W. The old B_pr transmission term was a
            # site-specific artefact, not NCRP - removed.
            # No independently verified numeric target exists for this form yet.
            DRs = (
                (d["Ddot0_max"] * d["alpha0_max"] + d["Ddot0_min"] * d["alpha0_min"])
                * d["U_G"] * d["A0"] * d["alpha_z"] * d["Az"]
            ) / ((d_h * d_r * d_z) ** 2)

            # Matches the reverted H_PS: energy-split alpha1, not the soft
            # spectrum substitution. Potato's real formula, verified.
            DR_ps = (
                (d["Ddot0_max"] * d["alpha1_max"] + d["Ddot0_min"] * d["alpha1_min"])
                * d["alpha_theta"] * d["A1"] * (d["F"] / 400.0)
            ) / ((d_sca * d_sec * d_zz) ** 2)

            DR_ls = (
                (d["Ddot0_max"] * d["alpha1_max"] + d["Ddot0_min"] * d["alpha1_min"])
                * d["L_f"] * d["A1"]
            ) / ((d_sec * d_zz) ** 2)

            DR_lt = (
                (d["Ddot0_max"] + d["Ddot0_min"]) * d["B_LT"] * d["L_f"]
            ) / (d_L ** 2)

            DR_ph = DR_ps + d["f"] * DRs + DR_ls + DR_lt
            IDR_ph = DR_ph * (10.0 ** (-x_ph_lead / tvl_ph))

            self._results = {
                "x_ph_lead": x_ph_lead,
                "DRs": DRs, "DR_ps": DR_ps, "DR_ls": DR_ls, "DR_lt": DR_lt,
                "DR_ph": DR_ph, "IDR_ph": IDR_ph,
            }

            report = [
                "=" * 65,
                "REQUIRED LEAD THICKNESS AND DOSE RATE CHECK",
                "=" * 65,
                f"B_ph (required transmission)   : {B_ph:.3e}",
                f"Lead thickness x_ph_lead       : {x_ph_lead:.2f} cm",
                "-" * 65,
                "Instantaneous dose rate through x_ph_lead:",
                f"  Wall scatter          DRs   : {DRs:.3e} Sv/h",
                f"  Patient scatter        DR_ps : {DR_ps:.3e} Sv/h",
                f"  Head leakage scatter   DR_ls : {DR_ls:.3e} Sv/h",
                f"  Leakage transmission   DR_lt : {DR_lt:.3e} Sv/h",
                "-" * 65,
                f"Total dose rate         DR_ph : {DR_ph:.3e} Sv/h",
                f"Instantaneous dose rate IDR_ph: {IDR_ph:.3e} Sv/h",
                "=" * 65,
            ]
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", "\n".join(report))

        except Exception as e:
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", f"Calculation failed:\n\n{type(e).__name__}: {e}")

    def validate(self) -> bool:
        return bool(self._results)

    def get_data(self) -> dict:
        return dict(self._results)


# ==============================================================================
# M STEP 5 : RESULTS, PHOTON ONLY
# ------------------------------------------------------------------------------
# Terminal screen for the M-branch (E_max < 10 MV). Neutron and capture gamma
# do not apply below 10 MV, so lead at the door is the whole answer.
# Renders in show() so returning from Back reflects any changed inputs.
# ==============================================================================

class Step_M_Results(WizardStep):
    """Final photon-only result: required lead at the door."""

    # Lead TVL for the scatter/leakage photon spectrum at the door.
    _TVL_PB = 0.6
    # log10(2) - one HVL expressed as a fraction of a TVL.
    _HVL_FRACTION = 0.301

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="Required Door Shielding: Photons",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 10))

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=10, pady=5)
        self.text_area = tk.Text(text_frame, wrap="word", height=16, width=80,
                                 font=("Consolas", 10))
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.text_area.yview)
        self.text_area.configure(yscrollcommand=scrollbar.set)
        self.text_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Terminal step. Restart returns to the start of this branch, not step 0.
        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Restart", command=self.wizard.restart_workflow).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def show(self) -> None:
        super().show()
        self._render()

    def _render(self) -> None:
        d = self.wizard.data
        try:
            x_pb = d["x_ph_lead"]
            idr = d["IDR_ph"]
            margin = self._HVL_FRACTION * self._TVL_PB
            x_total = x_pb + margin

            report = [
                "=" * 65,
                "REQUIRED DOOR SHIELDING: PHOTONS",
                "=" * 65,
                f"E_max = {d['E_max']} MV, below the 10 MV neutron threshold.",
                "Neutron and capture gamma components do not apply.",
                "-" * 65,
                f"Calculated lead thickness      : {x_pb:.2f} cm Pb",
                f"Plus 1 HVL safety margin       : {margin:.2f} cm Pb",
                f"Recommended total              : {x_total:.2f} cm Pb",
                "-" * 65,
                f"Instantaneous dose rate IDR_ph : {idr:.3e} Sv/h",
                f"                               = {idr * 1e6:.2f} uSv/h",
            ]

            # IPEM checks. Photons are the whole story below 10 MV.
            chk = ipem_idr_checks(idr, d["DR_ph"], d["weekly_beam_on_hours"])
            report += [
                "-" * 65,
                "IPEM IDR COMPLIANCE",
                f"  IDR total (shielded)      : {idr*1e6:.3f} uSv/h "
                f"(limit 7.5) -- {'PASS' if chk['IDR_check_pass'] else 'FAIL'}",
                f"    Margin                  : {chk['IDR_check_margin']*1e6:+.3f} uSv/h",
                f"  Weekly beam-on time       : {d['weekly_beam_on_hours']:.2f} h",
                f"  Rw (weekly, unshielded)   : {chk['Rw']*1e6:.3f} uSv",
                f"  R8h (daily-averaged)      : {chk['R8h']*1e6:.3f} uSv/h "
                f"(limit 0.5) -- {'PASS' if chk['R8h_check_pass'] else 'FAIL'}",
                f"    Margin                  : {chk['R8h_check_margin']*1e6:+.3f} uSv/h",
            ]
            if not chk["IDR_check_pass"]:
                report.append(
                    f"  FAILED: instantaneous IDR exceeds 7.5 uSv/h by "
                    f"{-chk['IDR_check_margin']*1e6:.3f} uSv/h. Thicken the door."
                )
            if not chk["R8h_check_pass"]:
                report.append(
                    f"  FAILED: R8h exceeds 0.5 uSv/h by "
                    f"{-chk['R8h_check_margin']*1e6:.3f} uSv/h. Thicken the door."
                )

            report += [
                "=" * 65,
                "",
                "Verify against your own calculation before signing anything.",
            ]
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", "\n".join(report))

        except KeyError as e:
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", f"Missing upstream value: {e}")

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        return {}


# ==============================================================================
# N STEP 1 : NEUTRON AND CAPTURE GAMMA AT THE DOOR
# ------------------------------------------------------------------------------
# Only reached when E_max >= 10 MV. Extends the M-branch photon result with
# neutron and capture-gamma dose through the door.
#
# S_0 and S_1 are NOT re-asked here. S_0 = Az and S_1 = A1 from
# Step_M_Geometry describe the same maze opening. Re-typing them risked the
# two numbers drifting apart, so this step reads them from wizard.data.
#
# Modified Kersey method (Wu & McGinley), reproduced term by term against
# the NCRP 151 worked example (Varian 1800, Sec 7.1.11-7.1.12).
# Constant is 2.4e-15, not 2.4 * 10e-15 (potato's version evaluated to
# 1e-14, ten times too large -- Python does not read "10e-15" as "times ten
# to the minus fifteen"). The S0/S1 ratio carries NO square root - the
# worked example only reproduces without it. phi_A third term uses 1.3,
# confirmed by the same worked example.
#
# Design goal allocation: neutrons and capture gamma each take P/2, sharing
# the goal with photons (Agios Andreas convention). See Step_M_Calculations.
# ==============================================================================

class Step_N_Params(WizardStep):
    """Neutron and capture gamma inputs and physics, extending the M-branch."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)
        self._results: dict = {}

        ttk.Label(
            self.frame, text="Neutron and Capture Gamma Parameters",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        canvas = tk.Canvas(self.frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.frame, orient="vertical", command=canvas.yview)
        scroll = ttk.Frame(canvas)
        scroll.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scroll, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        self._entries = {}
        self._reasons = {}

        def row(label_text: str, key: str, minv=None, maxv=None, default="", tooltip=""):
            r = ttk.Frame(scroll)
            r.pack(fill="x", pady=2)
            ttk.Label(r, text=label_text, width=40, anchor="w").pack(side="left")
            reason = ttk.Label(r, text="", foreground="red", font=("Segoe UI", 8))
            entry = _ValidatedEntry(r, minv=minv, maxv=maxv,
                                    label=label_text.rstrip(":"),
                                    reason_widget=reason, width=14)
            if default != "":
                entry.insert(0, str(default))
            else:
                entry.is_valid = False
            entry.pack(side="left", padx=4)
            reason.pack(side="left")
            if tooltip:
                _Tooltip(entry, tooltip)
            self._entries[key] = entry
            self._reasons[key] = reason

        row("beta, head transmission factor:", "beta", minv=0.0, maxv=1.0,
            tooltip="1.0 for lead head shielding, 0.85 for tungsten.")
        row("Q_n, neutron source strength:", "Q_n", minv=0.0,
            tooltip="Neutrons per Gy of photons at isocentre. Vendor-specific, energy-dependent.")
        row("d_n_1, isocentre to maze entrance [m]:", "d_n_1", minv=0.001,
            tooltip="Distance from isocentre to the maze entrance plane.")
        row("d_n_2, maze entrance to door [m]:", "d_n_2", minv=0.001,
            tooltip="Distance from the maze entrance plane to the door.")
        row("K, capture gamma ratio:", "K", minv=0.0,
            tooltip="Ratio of capture gamma dose equivalent to total neutron fluence at the maze entrance.")
        row("w_1_n, maze width at room entrance [m]:", "w_1_n", minv=0.001,
            tooltip="Width of the maze corridor specifically where it meets the room. "
                    "Can differ from the M-branch corridor width if the maze narrows "
                    "or widens along its run - this is that measurement, taken again "
                    "on purpose, not derived from earlier inputs.")
        row("Rr_L, room length without maze [m]:", "Rr_L", minv=0.001,
            tooltip="Treatment room length, excluding the maze corridor.")
        row("Rr_W, room width [m]:", "Rr_W", minv=0.001)
        row("Rr_H, room height [m]:", "Rr_H", minv=0.001)

        note = ttk.Label(
            scroll, foreground="#555555", font=("Segoe UI", 8, "italic"),
            wraplength=520, justify="left",
            text=(
                "S_0 is taken from the maze geometry entered earlier (Az) and "
                "is not asked again here. S_1 uses the width entered above, "
                "not a value from the M-branch - see its tooltip for why."
            ),
        )
        note.pack(anchor="w", padx=4, pady=(8, 4))

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def validate(self) -> bool:
        return all(entry._check() for entry in self._entries.values())

    def get_data(self) -> dict:
        d = {key: float(entry.get().strip()) for key, entry in self._entries.items()}
        w = self.wizard.data

        # S_0 = Az (maze corridor cross-section at the room-side opening)
        # is genuinely the same opening as Step_M_Geometry's Az and safe
        # to reuse.
        # S_1 is NOT A1, and NOT the M-branch's w_1 either. A1 is the maze
        # TONGUE width facing the room (w_tongue x h_e) - a different wall.
        # The M-branch's own w_1 is the corridor width at whatever point
        # Az is measured, which need not be the same point as the maze's
        # width at the room entrance. A maze can genuinely narrow or widen
        # along its run, so this is asked for separately as w_1_n, not
        # derived. Confirmed against a verified reference case where the
        # two widths differ (1.25 m vs 1.7 m in the same bunker).
        S_0 = w["Az"]
        S_1 = d["w_1_n"] * w["h_e"]
        d.update({"S_0": S_0, "S_1": S_1})

        S_r = 2.0 * (d["Rr_L"] * d["Rr_W"] + d["Rr_L"] * d["Rr_H"] + d["Rr_W"] * d["Rr_H"])
        d["S_r"] = S_r

        TVD_n = 2.06 * math.sqrt(S_1)
        d["TVD_n"] = TVD_n

        # Modified Kersey, Wu & McGinley. Direct + wall-reflected + thermal.
        # Third coefficient 1.3, confirmed by the NCRP worked example.
        phi_A = (
            (d["beta"] * d["Q_n"]) / (4.0 * math.pi * (d["d_n_1"] ** 2))
            + (5.4 * d["beta"] * d["Q_n"]) / (2.0 * math.pi * S_r)
            + (1.3 * d["Q_n"]) / (2.0 * math.pi * S_r)
        )
        d["phi_A"] = phi_A

        # Plain S0/S1, no square root. The NCRP worked example only
        # reproduces without it - the thesis version was wrong.
        H_n_D = (
            2.4 * 10e-15 * math.sqrt(S_0 / S_1) * phi_A
            * (1.64 * (10.0 ** (-(d["d_n_2"] / 1.9))) + (10.0 ** (-(d["d_n_2"] / TVD_n))))
        )
        d["H_n_D"] = H_n_D

        E_min, E_max = w["E_min"], w["E_max"]

        # Both energies photoproduce when E_min clears the threshold too.
        if E_min >= 10.0:
            H_n = (w["W_L_max"] + w["W_L_min"]) * H_n_D
        else:
            H_n = w["W_L_max"] * H_n_D

        # P/2: goal shared with photons and capture gamma.
        B_n = w["P"] / (2.0 * H_n)
        d["H_n"], d["B_n"] = H_n, B_n

        if B_n > 1.0:
            x_n_PBE = x_n_concrete = 0.0
        else:
            N_n = -math.log10(B_n)
            x_n_PBE = N_n * 4.5
            x_n_concrete = N_n * 15.0
        d["x_n_PBE"], d["x_n_concrete"] = x_n_PBE, x_n_concrete

        if E_min >= 10.0:
            DR_n = (w["Ddot0_max"] + w["Ddot0_min"]) * H_n_D
        else:
            DR_n = w["Ddot0_max"] * H_n_D
        IDR_n = DR_n * 10.0 ** (-x_n_PBE / 4.5)
        d["DR_n"], d["IDR_n"] = DR_n, IDR_n

        # Capture gamma only above 10 MV. At exactly 10 MV, photoneutron
        # production is not yet significant enough to pose a capture-gamma
        # hazard (Deye and Young, 1977), so the whole chain is skipped.
        if E_max == 10.0:
            d["cg_skipped"] = True
        else:
            d["cg_skipped"] = False

            # NCRP 151 Sec 7.1.11, TVD fixed at 5.4 m.
            TVD_cg = 5.4
            h_phi = phi_A * d["K"] * 10.0 ** (-d["d_n_2"] / TVD_cg)
            d["TVD_cg"], d["h_phi"] = TVD_cg, h_phi

            if E_min >= 10.0:
                H_cg = (w["W_L_max"] + w["W_L_min"]) * h_phi
            else:
                H_cg = w["W_L_max"] * h_phi

            # P/2: goal shared with photons and neutrons.
            B_cg = w["P"] / (2.0 * H_cg)
            d["H_cg"], d["B_cg"] = H_cg, B_cg

            if B_cg > 1.0:
                x_cg_lead = x_cg_concrete = x_cg_steel = 0.0
            else:
                N_cg = -math.log10(B_cg)
                x_cg_lead = N_cg * 6.1
                x_cg_concrete = N_cg * 33.0
                x_cg_steel = N_cg * 10.0
            d["x_cg_lead"] = x_cg_lead
            d["x_cg_concrete"] = x_cg_concrete
            d["x_cg_steel"] = x_cg_steel

            if E_min >= 10.0:
                DR_cg = (w["Ddot0_max"] + w["Ddot0_min"]) * h_phi
            else:
                DR_cg = w["Ddot0_max"] * h_phi
            IDR_cg = DR_cg * 10.0 ** (-x_cg_lead / 6.1)
            d["DR_cg"], d["IDR_cg"] = DR_cg, IDR_cg

        self._results = d
        return d


# ==============================================================================
# N STEP 2 : RESULTS, PHOTON + NEUTRON + CAPTURE GAMMA
# ------------------------------------------------------------------------------
# Terminal screen for the N-branch (E_max >= 10 MV). Combines the M-branch
# photon result with neutron and capture gamma from Step_N_Params.
# Each component gets its own +1 HVL margin, in its own material.
# ==============================================================================

class Step_N_Results(WizardStep):
    """Final combined result: photon, neutron, and capture gamma at the door."""

    _HVL_FRACTION = 0.301   # log10(2), one HVL as a fraction of a TVL

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="Required Door Shielding: Photon + Neutron + Capture Gamma",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 10))

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=10, pady=5)
        self.text_area = tk.Text(text_frame, wrap="word", height=24, width=80,
                                 font=("Consolas", 10))
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.text_area.yview)
        self.text_area.configure(yscrollcommand=scrollbar.set)
        self.text_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Restart", command=self.wizard.restart_workflow).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def show(self) -> None:
        super().show()
        self._render()

    def _render(self) -> None:
        d = self.wizard.data
        try:
            x_ph = d["x_ph_lead"]
            idr_ph = d["IDR_ph"]

            x_n_pbe = d["x_n_PBE"]
            x_n_conc = d["x_n_concrete"]
            idr_n = d["IDR_n"]

            cg_skipped = d.get("cg_skipped", False)

            m = self._HVL_FRACTION
            ph_margin = m * 0.6
            n_pbe_margin = m * 4.5
            n_conc_margin = m * 15.0

            report = [
                "=" * 65,
                "REQUIRED DOOR SHIELDING: PHOTON + NEUTRON + CAPTURE GAMMA",
                "=" * 65,
                f"E_max = {d['E_max']} MV, at or above the 10 MV neutron threshold.",
                "Design goal P split as P/2 for each component.",
                "-" * 65,
                "PHOTONS",
                f"  Calculated                : {x_ph:.2f} cm Pb",
                f"  Plus 1 HVL margin         : {ph_margin:.2f} cm Pb",
                f"  Recommended total         : {x_ph + ph_margin:.2f} cm Pb",
                f"  IDR                       : {idr_ph:.3e} Sv/h ({idr_ph*1e6:.2f} uSv/h)",
                "-" * 65,
                "NEUTRONS",
                f"  Calculated (PBE)          : {x_n_pbe:.2f} cm PBE",
                f"  Plus 1 HVL margin         : {n_pbe_margin:.2f} cm PBE",
                f"  Recommended total         : {x_n_pbe + n_pbe_margin:.2f} cm PBE",
                f"  Equivalent (concrete)     : {x_n_conc:.2f} cm, "
                f"+{n_conc_margin:.2f} HVL = {x_n_conc + n_conc_margin:.2f} cm",
                f"  IDR                       : {idr_n:.3e} Sv/h ({idr_n*1e6:.2f} uSv/h)",
                "-" * 65,
                "CAPTURE GAMMA",
            ]

            if cg_skipped:
                report += [
                    "  Capture gamma not calculated at 10 MV per NCRP 151 -- in this",
                    "  energy range, photoneutron production is not yet significant",
                    "  enough to pose a capture-gamma radiological hazard",
                    "  (Deye and Young, 1977).",
                ]
                idr_cg = 0.0
                dr_cg = 0.0
            else:
                x_cg_pb = d["x_cg_lead"]
                x_cg_conc = d["x_cg_concrete"]
                x_cg_steel = d["x_cg_steel"]
                idr_cg = d["IDR_cg"]
                dr_cg = d["DR_cg"]
                cg_pb_margin = m * 6.1
                cg_conc_margin = m * 33.0
                cg_steel_margin = m * 10.0
                report += [
                    f"  Lead                      : {x_cg_pb:.2f} cm, "
                    f"+{cg_pb_margin:.2f} HVL = {x_cg_pb + cg_pb_margin:.2f} cm",
                    f"  Concrete                  : {x_cg_conc:.2f} cm, "
                    f"+{cg_conc_margin:.2f} HVL = {x_cg_conc + cg_conc_margin:.2f} cm",
                    f"  Steel                     : {x_cg_steel:.2f} cm, "
                    f"+{cg_steel_margin:.2f} HVL = {x_cg_steel + cg_steel_margin:.2f} cm",
                    f"  IDR                       : {idr_cg:.3e} Sv/h ({idr_cg*1e6:.2f} uSv/h)",
                ]
            
            # IPEM checks: shielded IDRs against the instantaneous limit,
            # unshielded rates times beam-on hours against the daily average.
            idr_total = idr_ph + d["IDR_n"] + idr_cg
            dr_total = d["DR_ph"] + d["DR_n"] + dr_cg
            chk = ipem_idr_checks(idr_total, dr_total, d["weekly_beam_on_hours"])

            report += [
                "-" * 65,
                "IPEM IDR COMPLIANCE",
                f"  IDR total (shielded)      : {idr_total*1e6:.3f} uSv/h "
                f"(limit 7.5) -- {'PASS' if chk['IDR_check_pass'] else 'FAIL'}",
                f"    Margin                  : {chk['IDR_check_margin']*1e6:+.3f} uSv/h",
                f"  Weekly beam-on time       : {d['weekly_beam_on_hours']:.2f} h",
                f"  Rw (weekly, unshielded)   : {chk['Rw']*1e6:.3f} uSv",
                f"  R8h (daily-averaged)      : {chk['R8h']*1e6:.3f} uSv/h "
                f"(limit 0.5) -- {'PASS' if chk['R8h_check_pass'] else 'FAIL'}",
                f"    Margin                  : {chk['R8h_check_margin']*1e6:+.3f} uSv/h",
            ]
            if not chk["IDR_check_pass"]:
                report.append(
                    f"  FAILED: instantaneous IDR exceeds 7.5 uSv/h by "
                    f"{-chk['IDR_check_margin']*1e6:.3f} uSv/h. Thicken the door."
                )
            if not chk["R8h_check_pass"]:
                report.append(
                    f"  FAILED: R8h exceeds 0.5 uSv/h by "
                    f"{-chk['R8h_check_margin']*1e6:.3f} uSv/h. Thicken the door."
                )

            report += [
                "=" * 65,
                "",
                "These three components are not simply additive across materials.",
                "The door construction (layered lead / PBE / lead, or equivalent)",
                "determines how they combine. Verify against your own calculation",
                "before signing anything.",
            ]
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", "\n".join(report))

        except KeyError as e:
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", f"Missing upstream value: {e}")

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        return {}


# ==============================================================================
# D STEP 0 : PRIMARY BEAM AT DOOR CHECK
# ------------------------------------------------------------------------------
# No maze means the door has nothing between it and the primary beam except
# whatever wall is in front of it. If the primary beam reaches the door
# directly, the barrier there needs primary-barrier thickness, not the
# secondary-barrier scatter/leakage physics this branch calculates.
# Same three-strike pattern as Step_M_Primary. Design failure, not an input.
# ==============================================================================

class Step_D_Primary(WizardStep):
    """Yes/no: does the primary beam reach the door directly. Escalates on repeat 'Yes'."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        self._yes_count = 0

        ttk.Label(
            self.frame, text="Primary Beam at Door",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        self._primary_at_door = tk.BooleanVar(value=False)

        radio_frame = ttk.Frame(self.frame)
        radio_frame.pack(anchor="w", padx=20, pady=10)
        ttk.Label(
            radio_frame,
            text="Does the primary beam reach the door directly?",
        ).pack(side="left", padx=(0, 15))
        ttk.Radiobutton(
            radio_frame, text="Yes", variable=self._primary_at_door,
            value=True,
        ).pack(side="left", padx=5)
        ttk.Radiobutton(
            radio_frame, text="No", variable=self._primary_at_door,
            value=False,
        ).pack(side="left", padx=5)

        self._reason = ttk.Label(
            self.frame, text="", foreground="red",
            font=("Segoe UI", 10), wraplength=500, justify="left",
        )
        self._reason.pack(anchor="w", padx=20, pady=(0, 10))

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self._try_advance).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def _try_advance(self) -> None:
        if not self._primary_at_door.get():
            self._yes_count = 0
            self._reason.config(text="")
            self.wizard.next_step()
            return

        self._yes_count += 1

        if self._yes_count == 1:
            self._reason.config(text=(
                "The primary beam should never reach the door. This branch "
                "only handles scatter and leakage. Check your geometry."
            ))
        elif self._yes_count == 2:
            self._reason.config(text=(
                "Still no. A door in the primary beam needs primary barrier "
                "thickness, which this branch does not calculate. Move the door "
                "or move the isocentre, not the goalposts."
            ))
        else:
            messagebox.showerror(
                "Enough",
                "You need Jesus and I want no part in this."
            )
            self.wizard.exit_wizard()

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        return {"primary_at_door": self._primary_at_door.get()}


# ==============================================================================
# D STEP 1 : DIRECT DOOR GEOMETRY
# ------------------------------------------------------------------------------
# No maze. Door sits in the path of patient-scattered and leakage radiation
# only (never the primary beam - see Step_D_Primary check before this step).
# NCRP 151 Sec 2.2: U = 1 for both scatter and leakage at a secondary
# barrier. Not a user input, not a choice.
# ==============================================================================

class Step_D_Geometry(WizardStep):
    """Direct-door geometry and design parameters. No maze, no albedo tables."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="Direct Door Geometry",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 5))

        self._entries = {}
        self._reasons = {}

        def row(parent_frame, label_text, key, minv=None, maxv=None, default="", tooltip=""):
            r = ttk.Frame(parent_frame)
            r.pack(fill="x", pady=2)
            ttk.Label(r, text=label_text, width=34, anchor="w").pack(side="left")
            reason = ttk.Label(r, text="", foreground="red", font=("Segoe UI", 8))
            entry = _ValidatedEntry(r, minv=minv, maxv=maxv,
                                    label=label_text.rstrip(":"),
                                    reason_widget=reason, width=12)
            if default != "":
                entry.insert(0, str(default))
            else:
                entry.is_valid = False
            entry.pack(side="left", padx=4)
            reason.pack(side="left")
            if tooltip:
                _Tooltip(entry, tooltip)
            self._entries[key] = entry
            self._reasons[key] = reason

        p_frame = ttk.LabelFrame(self.frame, text="Design Parameters", padding=10)
        p_frame.pack(fill="x", padx=20, pady=5)
        row(p_frame, "Design goal P [Sv/week]:", "P", minv=1e-6, default="0.00006",
            tooltip="Weekly dose limit at the door.")
        row(p_frame, "Occupancy factor T:", "T", minv=0.025, maxv=1.0, default="1.0",
            tooltip="Fraction of time the space beyond the door is occupied.")
        row(p_frame, "Head leakage ratio L_f:", "L_f", minv=0.0, maxv=0.01, default="0.001",
            tooltip="Fraction of primary beam escaping the head. NCRP default 0.001.")
        row(p_frame, "Field size F [cm2]:", "F", minv=1.0, default="400.0",
            tooltip="Field size at isocentre. NCRP normalises to 400 cm2.")

        d_frame = ttk.LabelFrame(self.frame, text="Distances and Angle", padding=10)
        d_frame.pack(fill="x", padx=20, pady=5)
        row(d_frame, "Target to patient d_sca [m]:", "d_sca", minv=0.01, default="1.0",
            tooltip="Distance from target to scattering surface. NCRP standard 1.0 m.")
        row(d_frame, "Patient to door d_sec [m]:", "d_sec", minv=0.01,
            tooltip="Distance from the patient to the door.")
        row(d_frame, "Isocentre to door d_L [m]:", "d_L", minv=0.01,
            tooltip="Direct distance from isocentre to the door, for leakage.")
        row(d_frame, "Scatter angle theta [deg]:", "theta", minv=10.0, maxv=150.0, default="45.0",
            tooltip="Angle between the beam axis and the patient-to-door line. Table B.4 range 10-150 deg.")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def validate(self) -> bool:
        return all(entry._check() for entry in self._entries.values())

    def get_data(self) -> dict:
        data = {key: float(entry.get().strip()) for key, entry in self._entries.items()}
        return data


# ==============================================================================
# D STEP 2 : DIRECT DOOR DOSE AND BARRIER THICKNESS
# ------------------------------------------------------------------------------
# NCRP 151 Sec 2.2, Eq 2.7 (scatter) and Eq 2.8 (leakage), U=1 throughout.
# Scatter and leakage barriers are sized SEPARATELY then combined by the
# two-source rule (Sec 2.3), not simply added as doses.
# Photons are the only component here, so the full design goal P applies.
#
# TVLs come from NCRP Table B.5 (leakage and scattered photons, concrete),
# not the primary-beam table: leaked and scattered photons are softer and
# attenuate faster. Scatter gets a single angle-dependent TVL; leakage keeps
# the TVL1/TVLe split. Digitised at 6 and 18 MV only - 10 MV raises a clear
# ValueError rather than borrowing a neighbour's numbers.
# ==============================================================================

class Step_D_Calculations(WizardStep):
    """Direct scatter and leakage dose at the door, combined by the two-source rule."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)
        self._results: dict = {}

        ttk.Label(
            self.frame, text="Direct Door Dose and Thickness",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 10))

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=10, pady=5)
        self.text_area = tk.Text(text_frame, wrap="word", height=22, width=80,
                                 font=("Consolas", 10))
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.text_area.yview)
        self.text_area.configure(yscrollcommand=scrollbar.set)
        self.text_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Next", command=self.wizard.next_step).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def show(self) -> None:
        super().show()
        self._execute()

    def _execute(self) -> None:
        d = self.wizard.data
        self._results = {}

        try:
            P, T, L_f, F = d["P"], d["T"], d["L_f"], d["F"]
            d_sca, d_sec, d_L, theta = d["d_sca"], d["d_sec"], d["d_L"], d["theta"]
            W_min, W_max = d["W_min"], d["W_max"]
            E_max = d["E_max"]

            # U = 1 by NCRP 151 Sec 2.2 for secondary barriers. Multiplying
            # by one is a no-op, so it is documented here rather than coded.

            # Scatter fraction from the same Table B.4 lookup the M-branch
            # uses. Same physical quantity, same table, different geometry.
            a = lookup_scatter_fraction_b4(theta, E_max)

            # Eq 2.7 rearranged for H. Both energies' workloads contribute.
            H_ps = (
                a * (W_max + W_min) * T * (F / 400.0)
            ) / ((d_sca ** 2) * (d_sec ** 2))

            # Eq 2.8 rearranged for H. L_f is the leakage ratio (NCRP 1e-3),
            # exposed as a field so a vendor-stated value can override it.
            H_L = (L_f * (W_max + W_min) * T) / (d_L ** 2)

            # Full P: photons are the only component at a direct door.
            B_ps = P / H_ps if H_ps > 0.0 else float("inf")
            B_L = P / H_L if H_L > 0.0 else float("inf")

            # Table B.5 TVLs, each component's own. 10 MV raises here - the
            # gap is stated loudly instead of papered over.
            TVL_s = lookup_tvl_scatter(E_max, theta)
            TVL_L1, TVL_Le = lookup_tvl_leakage(E_max)

            def thickness_scatter(B: float) -> tuple[float, float]:
                """Single TVL - that is all Table B.5 gives for scatter."""
                if B >= 1.0:
                    return 0.0, 0.0
                n = -math.log10(B)
                return n, n * TVL_s

            def thickness_leakage(B: float) -> tuple[float, float]:
                """n and thickness with the TVL1/TVLe split. Eq 2.2 / 2.3."""
                if B >= 1.0:
                    return 0.0, 0.0
                n = -math.log10(B)
                t = n * TVL_L1 if n <= 1.0 else TVL_L1 + (n - 1.0) * TVL_Le
                return n, t

            def transmission_scatter(t: float) -> float:
                if t <= 0.0:
                    return 1.0
                return 10.0 ** (-t / TVL_s)

            def transmission_leakage(t: float) -> float:
                """Inverse of Eq 2.3: first decade through TVL1, the rest through TVLe."""
                if t <= 0.0:
                    return 1.0
                if t <= TVL_L1:
                    return 10.0 ** (-t / TVL_L1)
                return 10.0 ** (-(1.0 + (t - TVL_L1) / TVL_Le))

            n_ps, t_ps = thickness_scatter(B_ps)
            n_L, t_L = thickness_leakage(B_L)

            # Two-source rule, NCRP 151 Sec 2.3. Barriers are not additive.
            # The comparison and the added HVL use the CONTROLLING component's
            # own TVL, since a shared one no longer exists.
            ctrl_tvl = TVL_Le if t_L >= t_ps else TVL_s
            if abs(t_ps - t_L) >= ctrl_tvl:
                t_final = max(t_ps, t_L)
                rule_note = "Difference >= 1 TVL. Larger barrier used as-is."
            else:
                hvl_ctrl = 0.301 * ctrl_tvl
                t_final = max(t_ps, t_L) + hvl_ctrl
                rule_note = f"Difference < 1 TVL. Added 1 HVL ({hvl_ctrl:.2f} cm) to the larger."

            # Instantaneous dose rate. D0 plays the role W plays above.
            # Each component decays through its own TVL, not a shared one.
            D0_max = d["Ddot0_max"]
            IDR_ps_unshielded = (D0_max * a * (F / 400.0)) / ((d_sca ** 2) * (d_sec ** 2))
            IDR_L_unshielded = (D0_max * L_f) / (d_L ** 2)

            IDR_ps = IDR_ps_unshielded * transmission_scatter(t_final)
            IDR_L = IDR_L_unshielded * transmission_leakage(t_final)
            IDR_total = IDR_ps + IDR_L

            self._results = {
                "a_scatter": a, "H_ps": H_ps, "H_L": H_L,
                "B_ps": B_ps, "B_L": B_L,
                "n_ps": n_ps, "n_L": n_L, "t_ps": t_ps, "t_L": t_L,
                "t_final": t_final,
                "TVL_s": TVL_s, "TVL_L1": TVL_L1, "TVL_Le": TVL_Le,
                "IDR_ps": IDR_ps, "IDR_L": IDR_L, "IDR_total": IDR_total,
            }

            report = [
                "=" * 65,
                "NCRP 151 DIRECT DOOR RESULTS (SEC 2.2, U = 1)",
                "=" * 65,
                f"Scatter fraction a(theta={theta} deg, E_max={E_max} MV): {a:.3e}",
                f"TVLs (Table B.5, concrete): scatter {TVL_s:.1f} cm, "
                f"leakage {TVL_L1:.1f}/{TVL_Le:.1f} cm",
                "-" * 65,
                f"H_ps (patient scatter, unshielded) : {H_ps:.3e} Sv/wk",
                f"H_L  (leakage, unshielded)         : {H_L:.3e} Sv/wk",
                "-" * 65,
                f"B_ps : {B_ps:.3e}   n_ps : {n_ps:.2f}   t_ps : {t_ps:.2f} cm",
                f"B_L  : {B_L:.3e}   n_L  : {n_L:.2f}   t_L  : {t_L:.2f} cm",
                "-" * 65,
                "Two-source rule (NCRP 151 Sec 2.3):",
                f"  {rule_note}",
                f"  Combined barrier thickness      : {t_final:.2f} cm",
                "-" * 65,
                "IDR through combined barrier (per-component TVLs):",
                f"  Scatter   : {IDR_ps:.3e} Sv/h",
                f"  Leakage   : {IDR_L:.3e} Sv/h",
                f"  Total     : {IDR_total:.3e} Sv/h ({IDR_total*1e6:.2f} uSv/h)",
                "=" * 65,
            ]
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", "\n".join(report))

        except Exception as e:
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", f"Calculation failed:\n\n{type(e).__name__}: {e}")

    def validate(self) -> bool:
        return bool(self._results)

    def get_data(self) -> dict:
        return dict(self._results)


# ==============================================================================
# D STEP 3 : RESULTS
# ------------------------------------------------------------------------------
# Terminal screen for the D-branch. The two-source rule already applied its
# own safety margin where needed, so no blanket +1 HVL is added again here.
# ==============================================================================

class Step_D_Results(WizardStep):
    """Final direct-door result: combined barrier thickness and IDR."""

    def __init__(self, parent: tk.Frame, wizard: 'ShieldingWizard') -> None:
        super().__init__(parent, wizard)
        self.frame = ttk.Frame(parent)

        ttk.Label(
            self.frame, text="Required Door Shielding: Direct Line of Sight",
            font=("Segoe UI", 12, "bold"),
        ).pack(pady=(15, 10))

        text_frame = ttk.Frame(self.frame)
        text_frame.pack(fill="both", expand=True, padx=10, pady=5)
        self.text_area = tk.Text(text_frame, wrap="word", height=16, width=80,
                                 font=("Consolas", 10))
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.text_area.yview)
        self.text_area.configure(yscrollcommand=scrollbar.set)
        self.text_area.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        btn_frame = ttk.Frame(self.frame)
        btn_frame.pack(pady=15)
        ttk.Button(btn_frame, text="Restart", command=self.wizard.restart_workflow).pack(side="left", padx=5)
        ttk.Button(btn_frame, text="Back", command=self.wizard.prev_step).pack(side="left", padx=5)

    def show(self) -> None:
        super().show()
        self._render()

    def _render(self) -> None:
        d = self.wizard.data
        try:
            report = [
                "=" * 65,
                "REQUIRED DOOR SHIELDING: DIRECT LINE OF SIGHT",
                "=" * 65,
                f"E_max = {d['E_max']} MV. No maze. Door sees patient scatter",
                "and head leakage directly, never the primary beam.",
                "-" * 65,
                f"Controlling component     : "
                f"{'leakage' if d['t_L'] >= d['t_ps'] else 'scatter'}",
                f"Combined barrier thickness : {d['t_final']:.2f} cm "
                f"(TVL_s={d['TVL_s']:.1f}, TVL_L1={d['TVL_L1']:.1f}, "
                f"TVL_Le={d['TVL_Le']:.1f})",
                "-" * 65,
                f"IDR at the door            : {d['IDR_total']:.3e} Sv/h "
                f"({d['IDR_total']*1e6:.2f} uSv/h)",
                "=" * 65,
                "",
                "TVLs are NCRP 151 Table B.5 (leakage and scattered photons,",
                "concrete), digitised at 6 and 18 MV only. 10 MV is a known gap.",
                "",
                "Verify against your own calculation before signing anything.",
            ]
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", "\n".join(report))
        except KeyError as e:
            self.text_area.delete("1.0", "end")
            self.text_area.insert("1.0", f"Missing upstream value: {e}")

    def validate(self) -> bool:
        return True

    def get_data(self) -> dict:
        return {}


# ==============================================================================
# NAVIGATION BAR LABELS
# ------------------------------------------------------------------------------
# One short label per step class for the bottom nav bar. Branch tag plus a
# number, matching the section banners above. Anything not listed falls back
# to the class name, which is ugly but at least honest.
# ==============================================================================

_NAV_LABELS = {
    "Step_BranchSelector":  "START",
    "Step_RT_Checklist":    "RT",
    "Step_SETUP_Workload":  "SETUP 1",
    "Step_SETUP_IMRT":      "SETUP 2",
    "Step_SETUP_TBI":       "SETUP 3",
    "Step_SETUP_Maze":      "SETUP 4",
    "Step_M_Primary":       "M 1",
    "Step_M_Geometry":      "M 2",
    "Step_M_Calculations":  "M 3",
    "Step_M_Bph":           "M 4",
    "Step_M_Results":       "M 5",
    "Step_N_Params":        "N 1",
    "Step_N_Results":       "N 2",
    "Step_D_Primary":       "D 0",
    "Step_D_Geometry":      "D 1",
    "Step_D_Calculations":  "D 2",
    "Step_D_Results":       "D 3",
}


# ==============================================================================
# WIZARD CONTROLLER
# ------------------------------------------------------------------------------
# Owns the data dictionary, the step list, the routing between them, and the
# bottom navigation bar. Nav bar shows one square per VISITED step. Clicking
# a square jumps back to that step and rolls back everything after it.
# Clicking START is a full restart, because changing branch invalidates all.
# ==============================================================================

class ShieldingWizard(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("RT Bunker Shielding Wizard")
        self.geometry("900x780")

        self.data: dict = {}
        self.current_step_index: int = 0
        self.step_history: list[int] = []
        self.steps: list[WizardStep] = []

        # Nav bar packed first at the bottom so the container fills the rest.
        self.nav_frame = ttk.Frame(self, padding=(8, 4))
        self.nav_frame.pack(side="bottom", fill="x")
        ttk.Separator(self, orient="horizontal").pack(side="bottom", fill="x")

        self.container = ttk.Frame(self)
        self.container.pack(side="top", fill="both", expand=True)

        self.initialize_wizard()

    # --- lifecycle -----------------------------------------------------------

    def initialize_wizard(self) -> None:
        """Full restart from the branch selector. Everything is forgotten."""
        for step in self.steps:
            step.hide()
            if step.frame:
                step.frame.destroy()

        self.data = {}
        self.steps = [Step_BranchSelector(self.container, self)]
        self.current_step_index = 0
        self.step_history = [0]
        self.show_current_step()

    def restart_workflow(self) -> None:
        """Restart at the start of the current branch. Branch choice is kept."""
        for step in self.steps[1:]:
            step.hide()
            if step.frame:
                step.frame.destroy()
        self.data = {}
        self.steps = self.steps[:1]
        self.current_step_index = 0
        self.step_history = [0]
        # Only the RT branch exists today. When CT arrives, remember which
        # branch was active and route accordingly.
        self.start_radiotherapy_branch()

    def show_current_step(self) -> None:
        for step in self.steps:
            step.hide()
        if 0 <= self.current_step_index < len(self.steps):
            self.steps[self.current_step_index].show()
        self._refresh_nav()

    # --- navigation bar ------------------------------------------------------

    def _refresh_nav(self) -> None:
        for child in self.nav_frame.winfo_children():
            child.destroy()
        for i, step in enumerate(self.steps):
            name = type(step).__name__
            text = _NAV_LABELS.get(name, name)
            is_current = (i == self.current_step_index)
            btn = tk.Button(
                self.nav_frame, text=text, width=8,
                relief="sunken" if is_current else "raised",
                bg="#dbe6f5" if is_current else "#f0f0f0",
                activebackground="#c9d9ef",
                font=("Segoe UI", 8, "bold" if is_current else "normal"),
                command=lambda idx=i: self._goto_step(idx),
            )
            btn.pack(side="left", padx=2)

    def _goto_step(self, index: int) -> None:
        """Jump back to a visited step, rolling back everything after it."""
        if index == 0:
            self.initialize_wizard()
            return
        if index >= self.current_step_index:
            return
        for step in reversed(self.steps[index + 1:]):
            step.rollback_data(self.data)
            step.hide()
            if step.frame:
                step.frame.destroy()
        self.steps = self.steps[:index + 1]
        while self.step_history and self.step_history[-1] != index:
            self.step_history.pop()
        if not self.step_history:
            self.step_history = [index]
        self.current_step_index = index
        self.show_current_step()

    # --- branch entry --------------------------------------------------------

    def start_radiotherapy_branch(self) -> None:
        self._ensure_step_exists(Step_RT_Checklist)
        self._advance_to_next_available()

    def start_radiology_branch(self) -> None:
        messagebox.showinfo("Pending", "CT branch is pending implementation (Sutton/BIR).")

    # --- routing -------------------------------------------------------------

    def next_step(self) -> None:
        current = self.steps[self.current_step_index]

        if not current.validate():
            return

        current.apply_data(self.data)

        if isinstance(current, Step_RT_Checklist):
            self._ensure_step_exists(Step_SETUP_Workload)
            self._advance_to_next_available()

        elif isinstance(current, Step_SETUP_Workload):
            self._ensure_step_exists(Step_SETUP_IMRT)
            self._advance_to_next_available()

        elif isinstance(current, Step_SETUP_IMRT):
            self._ensure_step_exists(Step_SETUP_TBI)
            self._advance_to_next_available()

        elif isinstance(current, Step_SETUP_TBI):
            self._ensure_step_exists(Step_SETUP_Maze)
            self._advance_to_next_available()

        elif isinstance(current, Step_SETUP_Maze):
            if self.data.get("maze_exists"):
                self._ensure_step_exists(Step_M_Primary)
            else:
                self._ensure_step_exists(Step_D_Primary)
            self._advance_to_next_available()

        # --- M / N branch
        elif isinstance(current, Step_M_Primary):
            self._ensure_step_exists(Step_M_Geometry)
            self._advance_to_next_available()

        elif isinstance(current, Step_M_Geometry):
            self._ensure_step_exists(Step_M_Calculations)
            self._advance_to_next_available()

        elif isinstance(current, Step_M_Calculations):
            self._ensure_step_exists(Step_M_Bph)
            self._advance_to_next_available()

        elif isinstance(current, Step_M_Bph):
            if self.data["E_max"] >= 10.0:
                self._ensure_step_exists(Step_N_Params)
            else:
                self._ensure_step_exists(Step_M_Results)
            self._advance_to_next_available()

        elif isinstance(current, Step_N_Params):
            self._ensure_step_exists(Step_N_Results)
            self._advance_to_next_available()

        # --- D branch
        elif isinstance(current, Step_D_Primary):
            self._ensure_step_exists(Step_D_Geometry)
            self._advance_to_next_available()

        elif isinstance(current, Step_D_Geometry):
            self._ensure_step_exists(Step_D_Calculations)
            self._advance_to_next_available()

        elif isinstance(current, Step_D_Calculations):
            self._ensure_step_exists(Step_D_Results)
            self._advance_to_next_available()

    def prev_step(self) -> None:
        """Go back and destroy later steps so they rebuild against fresh data."""
        if len(self.step_history) > 1:
            self.step_history.pop()
            prev_idx = self.step_history[-1]

            for step in reversed(self.steps[prev_idx + 1:]):
                step.rollback_data(self.data)
                step.hide()
                if step.frame:
                    step.frame.destroy()

            self.steps = self.steps[:prev_idx + 1]
            self.current_step_index = prev_idx
            self.show_current_step()

    # --- helpers -------------------------------------------------------------

    def _ensure_step_exists(self, step_class) -> None:
        for step in self.steps:
            if isinstance(step, step_class):
                return
        self.steps.append(step_class(self.container, self))

    def _advance_to_step(self, index: int) -> None:
        self.step_history.append(index)
        self.current_step_index = index
        self.show_current_step()

    def _advance_to_next_available(self) -> None:
        self._advance_to_step(self.current_step_index + 1)

    def exit_wizard(self) -> None:
        self.quit()
        self.destroy()


# ==============================================================================
# EXECUTION
# ==============================================================================

if __name__ == "__main__":
    app = ShieldingWizard()
    app.mainloop()