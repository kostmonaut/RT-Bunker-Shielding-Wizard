# RT Bunker Door Shielding Wizard

This tool partially automates radiotherapy bunker door shielding calculations following NCRP 151 methodology. It evaluates photon, neutron, and capture gamma attenuation for high-energy LINACs. It is designed strictly as a supplementary verification instrument to cross-check analytical estimations, not to replace professional judgement.

Originally developed as the practical implementation of the master's research thesis *COMPUTATIONAL MODEL FOR RADIOTHERAPY BUNKER DOOR NEUTRON SHIELDING CALCULATIONS* in medical physics.

---

## Features

* **Comprehensive Calculation:** Evaluates shielding requirements for photons, neutrons, and capture gamma rays.
* **Workload Modulation:** Adjusts primary and leakage workloads for Intensity-Modulated Radiotherapy (IMRT) and Total Body Irradiation (TBI).
* **Regulatory Compliance:** Evaluates Instantaneous Dose Rates (IDR) and time-averaged equivalent doses against established limits.
* **Material Support:** Accommodates composite barrier designs including concrete, lead, borated polyethylene, and steel.

---

## Requirements & Usage

**Requires Python 3.9 or newer.** The program is written in standard Python and utilizes built-in libraries (`tkinter`, `math`). No external dependencies are strictly required to run the source code.

### Running from Source
Execute the script directly from your terminal:
```bash
python shielding_wizard.py
