# RT Bunker Door Shielding Wizard

This tool partially automates radiotherapy bunker door shielding calculations following NCRP 151 methodology. It evaluates photon, neutron, and capture gamma attenuation for high-energy LINACs. It is designed strictly as a supplementary verification instrument to cross-check analytical estimations, not to replace professional judgement.

## Features
* Calculates shielding requirements for photons, neutrons, and capture gamma rays.
* Adjusts primary and leakage workloads for Intensity-Modulated Radiotherapy (IMRT) and Total Body Irradiation (TBI).
* Evaluates Instantaneous Dose Rates (IDR) and time-averaged equivalent doses against regulatory limits.
* Supports composite barrier designs (Concrete, Lead, Borated Polyethylene, Steel).

## Requirements & Usage
The programme is written in Python and utilises standard libraries (`tkinter`, `math`). 

Run the script from your terminal:
```bash
python shielding_wizard.py
