# Quantum-6502

**ZZ coupling is not a defect. It is a resource.**

A self-stabilizing quantum BIOS and OS interface that turns ZZ coupling drift  
into a resource — achieving up to ×2.66 effective T2 on IBM Heron hardware.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19029431.svg)](https://doi.org/10.5281/zenodo.19029431)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

---

## Key Results (IBM Heron r2, ibm_marrakesh, 2026-03-15)

| Metric | Value |
|---|---|
| GPS (Gate Phase Stabilization) | +0.726 ± 0.011 |
| z-score | 66.0σ (single condition) |
| z_combined (Holm correction) | 8.5σ (8 conditions) |
| ZZ drift range | 57× (0.063–3.60 kHz, 6 days) |
| QFEED+ net phase | 0.000 rad (confirmed) |
| Effective T2 (projected, Heron r3) | ×2.66 |

> **Note:** The z-score of 66.0σ represents peak statistical significance
> achieved under localized optimal drift conditions (τ=50µs, N=1, single condition).
> z_combined = 8.5σ after Holm correction across 8 conditions.
> Continuous verification across backends and qubit pairs is ongoing
> to characterize non-Markovian hardware fluctuations over extended periods.

---

## Architecture

```
┌─────────────────────────────────────────┐
│        Q-Script API  (qscript.py)       │
│  qs.status() / qs.run() / qs.explain()  │
├─────────────────────────────────────────┤
│   Quantum-6502 BIOS v1.2.2             │
│   (eedt_bios.py)                        │
│                                         │
│  Layer 0  Universal Pair Discovery      │
│  ULC      Ultra-Light Calibration       │
│  Layer 2  Operating Window Calculation  │
│  Layer 3  GO / DEGRADED / NO-GO         │
│  EEDT     QFEED+ Adaptive Feedforward   │
├─────────────────────────────────────────┤
│     Qiskit Runtime / IBM Backend        │
└─────────────────────────────────────────┘
```

### OS Analogy

| Classical OS | Quantum-6502 |
|---|---|
| Device detection | Layer 0: Universal Pair Discovery |
| Driver init | ULC: Auto-calibration |
| Scheduler | Layer 2: Window calculation |
| Kernel decision | Layer 3: GO/DEGRADED/NO-GO |
| Execution | EEDT + QFEED+ |

---

## Quick Start

### 1. Install

```bash
pip install -r requirements.txt
```

### 2. Set your IBM Quantum token

```bash
export IBM_QUANTUM_TOKEN="your_token_here"
```

Or edit `TOKEN` in `eedt_bios.py`:

```python
TOKEN = "your_token_here"
BACKEND_NAME = "ibm_marrakesh"   # or "ibm_kingston"
```

### 3. Run BIOS directly

```bash
python eedt_bios.py
```

The BIOS will automatically:
1. Discover the best qubit pair from calibration data (0 jobs)
2. Run Ultra-Light Calibration (9 circuits, 2 jobs)
3. Calculate operating window (τ*, ω_ZZ·T2, N*, φ_ff)
4. Execute EEDT with QFEED+ adaptive feedforward
5. Save results to CSV

### 4. Use Q-Script API

```python
import qscript

qs = qscript.QuantumOS(backend="ibm_marrakesh")
print(qs.status())

result = qs.run(goal="maximize_fidelity")
print(result)
# → {"gps": +0.41, "fidelity": 0.89, "pair": [95, 94], "verdict": "confirmed", ...}

# Why did the OS choose this pair?
print(qs.explain())

# Research mode — full physics disclosure
print(qs.debug())
```

---

## Core Concepts

### EEDT (Entanglement-Enhanced Dynamical Tracking)

ZZ coupling between superconducting qubits is normally treated as noise to be eliminated.  
EEDT inverts this: the deterministic phase accumulation from ZZ coupling is tracked  
and corrected in real time via Mid-Circuit Measurement (MCM) and feedforward.

**Operating window:** ω_ZZ · T2 ≈ 1–10  
**Optimal phase:** φ* ≈ 0.873 rad  
**Optimal storage time:** τ* = φ* / (2π × ν_ZZ_code)

### QFEED+ ($C3')

Extension of the feedforward instruction that compensates for ν_ZZ drift  
across any backend, any qubit pair.

```
φ_ff = φ_target + Δφ_comp
Δφ_comp = 2π × (ν_ZZ_actual - ν_ZZ_target) × (τ/N)
```

Confirmed: net phase = 0.000 rad on ibm_kingston Q0+Q1.

### Universal Pair Discovery (Layer 0)

Automatically scores all coupling-map pairs using calibration data (0 jobs):

```
score = (T2 / 500µs) × (1 - CZ_error × 20) × (1 - RO_error × 10)
```

Selects top-5 (exploitation) + median-5 (exploration) = up to 10 candidates.  
No-GO pairs are skipped instantly — no waiting.

### 369 Operating Modes

| Mode | N | τ* | Platform |
|---|---|---|---|
| Mode-3 | 3 | ~40µs | IBM standard |
| Mode-6 | 6 | ~28µs | High-frequency |
| Mode-9 | 9 | ~24µs | IQM QB20 range |

---

## Instruction Set ($C0–$C8 + $C3')

| Opcode | Instruction | Description |
|---|---|---|
| $C0 | QINIT | Qubit initialization |
| $C1 | QGATE | Gate application with ZZ correction |
| $C2 | QMCM | Mid-Circuit Measurement |
| $C3 | QFEED | Feedforward Rz (Lambert W embedded) |
| $C3' | QFEED+ | Adaptive feedforward with ν_ZZ compensation |
| $C4 | QSYNC | ZZ phase synchronization |
| $C5 | QERR | Separation Theorem N* check |
| $C6 | QWAIT | Poisson optimal wait |
| $C7 | QBRANCH | Conditional branch on measurement |
| $C8 | QMODE | 369 mode auto-configuration |

---

## File Structure

```
Quantum-6502/
├── eedt_bios.py          # BIOS v1.2.2 — Universal Pair Discovery
├── qscript.py            # Q-Script v0.1 — OS interface API
├── requirements.txt      # Python dependencies
└── README.md
```

---

## Platform Support

| Platform | Feedforward | EEDT | ω_ZZ·T2 | Status |
|---|---|---|---|---|
| IBM Heron r2 (ibm_marrakesh) | ✓ | ✓ | 5.91 (paper) | **GPS confirmed** |
| IBM Heron r2 (ibm_kingston) | ✓ | ✓ | TBD | QFEED+ confirmed |
| IBM Heron r3 (ibm_pittsburgh) | ✓ | ✓ (projected) | ~6.6 | Pending |
| IQM Garnet | ✗ | ✗ | >26 (window exceeded) | Dephasing-type ZZ |

---

## Theory

**Separation Theorem:** N ≤ N* = 5 for sign stability at τ = 50µs  
**Lambert W formula:** τ*(N) = τ_LW / W(α·N·exp(α)), α = 0.38  
**T2 extension:** T2_eff = T2 / (1 - ω_ZZ·T2/(1+ω_ZZ·T2) × GPS)  
**Projected on Heron r3:** T2_eff ≈ 930µs (×2.66, T2=350µs)

---

## Citation

```bibtex
@misc{okuda2026quantum6502,
  author    = {Okuda, Takeshi},
  title     = {A Self-Stabilizing Qubit Storage Unit with a Defined
               Operating Window on IBM Heron r2},
  year      = {2026},
  doi       = {10.5281/zenodo.19029431},
  publisher = {Zenodo}
}
```

---

## Author

**Takeshi Okuda**  
Independent Theoretical Contributor, Osaka, Japan  
GitHub: [github.com/okudat9](https://github.com/okudat9)  
DOI: [10.5281/zenodo.19029431](https://doi.org/10.5281/zenodo.19029431)

---

## License

CC BY 4.0 — Free to use, cite, and build upon with attribution.
