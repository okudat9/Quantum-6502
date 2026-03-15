# EEDT: Entanglement-Enhanced Dynamical Tracking

**ZZ Coupling as a Stochastic Synchronization Clock for Superconducting Qubits**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19029431.svg)](https://doi.org/10.5281/zenodo.19029431)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)

**Author:** Takeshi Okuda — Independent Quantum Researcher, Japan  
**Contact:** o93dice@gmail.com  
**Zenodo:** [10.5281/zenodo.19029431](https://doi.org/10.5281/zenodo.19029431) (v5, 2026-03-15)

---

## Core Idea

```
ZZ = noise      →  ZZ = stochastic synchronization clock
(standard view)     (EEDT view)
```

τ* = T₂ / (ω_ZZ · T₂ + 1)   φ* ≈ 1 (Poisson MLE optimum)

---

## Key Results

| | |
|---|---|
| Hardware | ibm_marrakesh Q94–Q95 (IBM Heron r2) |
| GPS (primary) | +0.037 ± 0.011, z = 3.35σ |
| GPS (ref₀, follow-up) | **+0.726 ± 0.011, z = 66.0σ** |
| Combined z | 8.5σ (8/10 Holm-corrected) |
| τ* agreement | <5% from analytic formula |
| ω_ZZ·T₂ operating window | 1–10 |

---

## Repository Structure

```
├── paper/                      # LaTeX source + compiled PDF (v5)
├── figures/                    # All 4 figures + generation scripts
├── experiments/
│   ├── ibm/                    # IBM Quantum experiment scripts
│   └── iqm/                    # IQM Garnet experiment scripts
└── quantum_6502/               # Quantum-6502 emulator
    ├── quantum_6502.py         # Official opcodes ($C0–$C7)
    ├── quantum_6502_illegal.py # Illegal opcodes ($D0–$D5)
    └── quantum_6502_v01.py     # Prototype v0.1 (full emulator)
```

---

## Quick Start

```bash
# Figures
pip install numpy matplotlib scipy
python figures/gen_fig4_v5.py

# Quantum-6502 emulator
pip install qiskit qiskit-aer scipy numpy
python quantum_6502/quantum_6502_v01.py

# IBM experiments (requires IBM Quantum credentials)
pip install qiskit qiskit-ibm-runtime scipy numpy
# Set token = "YOUR_IBM_QUANTUM_TOKEN" in script
python experiments/ibm/eedt_minimal_v2.py

# IQM experiments
pip install "iqm-client[qiskit]>=33.0,<34.0" numpy
python experiments/iqm/step1_zz_ramsey.py --token YOUR_IQM_TOKEN
```

---

## Quantum-6502

EEDT protocol as a MOS 6502-inspired instruction set:

| Opcode | Instruction | Function |
|--------|-------------|----------|
| $C0 | QINIT | Qubit initialization |
| $C1 | QGATE | H / X / CX / Rz |
| $C2 | QMCM | Mid-circuit measurement |
| $C3 | QFEED | EEDT feedforward Rz (Lambert W) |
| $C4 | QSYNC | ZZ phase synchronization |
| $C5 | QERR | Separation Theorem N* check |
| $C6 | QWAIT | Poisson optimal wait |
| $C7 | QBRANCH | Conditional branch |

**$QT = 40 µs** (time budget)  **QOVF** fires at N > N* = 5

---

## IBM Quantum Job IDs

**GPS v8 (2026-02):** `d6g3ue0ddp9c73cf3k60` `d6go67qthhns73916ocg` `d6gphf648nic73am518g` `d6gprum48nic73am5eg0`

**Follow-up (2026-03-15):** `d6r1j9i0q0ls73cstq10` `d6r1mki0q0ls73csttdg` `d6r1pujopkic73fil83g` `d6r1sm4u243c73a0tcl0` `d6r20vq0q0ls73csu890` `d6r234vr88ds73dd0i40` `d6r24u7r88ds73dd0k9g` `d6r27cropkic73filn50` `d6r2alfr88ds73dd0r8g` `d6r2dv3opkic73filv1g` `d6r2diq0q0ls73csumhg` `d6r3ij20q0ls73csvu3g`

---

## Citation

```bibtex
@misc{okuda2026eedt,
  author    = {Okuda, Takeshi},
  title     = {{EEDT}: Entanglement-Enhanced Dynamical Tracking
               for {ZZ}-Coupled Superconducting Qubits},
  year      = {2026},
  doi       = {10.5281/zenodo.19029431},
  publisher = {Zenodo},
  note      = {Version 5, 2026-03-15}
}
```

**License:** Code: AGPL-3.0 / Paper: CC BY 4.0  
*AI assistance (Claude, Anthropic) used for computation and LaTeX. All scientific decisions made by the author.*
