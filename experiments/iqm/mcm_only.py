"""
MCM バックアクション分離テスト — IQM Garnet QB12-QB17
====================================================
「MCM測定自体がdephasingを起こしているか」を確認

3回路 × 1000shots:
  REF:      MCMなし（ベースライン）
  MCM_only: MCMあり、Rzなし（バックアクションのみ）
  EEDT_neg: MCMあり、Rz(-phi)

MCM_only ≈ REF → MCM自体は問題なし → Rz角度の問題
MCM_only << REF → MCM自体がdephasing → 原理的限界

消費: 3回路 × 1000shots ≈ 4〜6 credits

実行:
  python mcm_only.py --token YOUR_TOKEN
"""

import argparse, json, sys
import numpy as np
from datetime import datetime
from pathlib import Path

QB_ANCILLA = 11   # QB12（T2=10.35µs）
QB_TARGET  = 16   # QB17（T2=13.04µs）
SHOTS      = 1000
NU_ZZ_HZ   = 35.86e3
TAU_US     = 3.0

OUTPUT_DIR = Path(f"mcm_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}")

def connect(token, device):
    from iqm.qiskit_iqm import IQMProvider
    url = f"https://cocos.resonance.meetiqm.com/{device}"
    provider = IQMProvider(url, token=token)
    backend  = provider.get_backend()
    print(f"[OK] {device}  qubits={backend.num_qubits}")
    return backend

def run(backend, qc, shots, label):
    from iqm.qiskit_iqm import transpile_to_IQM
    tc  = transpile_to_IQM(qc, backend)
    job = backend.run(tc, shots=shots)
    counts = job.result().get_counts()
    total  = sum(counts.values())
    print(f"  [{label}] {dict(sorted(counts.items()))}")
    return counts, total

def p0(counts, total, bit_index=0):
    if total == 0: return 0.5
    pos = -(bit_index + 1)
    return sum(v for k,v in counts.items()
               if len(k) > bit_index and k[pos]=='0') / total

def build_ref(qt, tau_us):
    """REF: MCMなし"""
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(qt + 1, 1)
    qc.h(qt)
    qc.delay(int(tau_us * 1000), qt, unit="ns")
    qc.h(qt)
    qc.measure(qt, 0)
    return qc

def build_mcm_only(qa, qt, tau_us):
    """MCMのみ: ancillaを測定するがRzは当てない"""
    from qiskit import QuantumCircuit
    nq = max(qa, qt) + 1
    qc = QuantumCircuit(nq, 2)
    qc.h(qa)
    qc.h(qt)
    qc.delay(int(tau_us * 1000), qa, unit="ns")
    qc.delay(int(tau_us * 1000), qt, unit="ns")
    qc.measure(qa, 0)  # MCM（Rzなし）
    qc.h(qt)
    qc.measure(qt, 1)
    return qc

def build_eedt(qa, qt, tau_us, nu_zz_hz, rz_sign=-1):
    """EEDT: MCM + Rz補正"""
    from qiskit import QuantumCircuit
    nq  = max(qa, qt) + 1
    phi = 2 * np.pi * nu_zz_hz * tau_us * 1e-6
    qc  = QuantumCircuit(nq, 2)
    qc.h(qa)
    qc.h(qt)
    qc.delay(int(tau_us * 1000), qa, unit="ns")
    qc.delay(int(tau_us * 1000), qt, unit="ns")
    qc.measure(qa, 0)
    qc.rz(rz_sign * phi, qt)
    qc.h(qt)
    qc.measure(qt, 1)
    return qc

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--token",  required=True)
    p.add_argument("--device", default="garnet")
    p.add_argument("--qa",  type=int, default=QB_ANCILLA)
    p.add_argument("--qt",  type=int, default=QB_TARGET)
    p.add_argument("--tau", type=float, default=TAU_US)
    args = p.parse_args()

    OUTPUT_DIR.mkdir(exist_ok=True)
    phi = 2 * np.pi * NU_ZZ_HZ * args.tau * 1e-6

    print("=" * 52)
    print("  MCM バックアクション分離テスト")
    print(f"  QB{args.qa+1}(ancilla) → QB{args.qt+1}(target)")
    print(f"  tau={args.tau}µs  phi={phi:.4f}")
    print(f"  3回路 × {SHOTS}shots  推定~4〜6 credits")
    print("=" * 52)

    backend = connect(args.token, args.device)

    # --- REF ---
    print("\n--- REF（MCMなし）---")
    cnt_r, tot_r = run(backend, build_ref(args.qt, args.tau),
                       SHOTS, "REF")
    F_ref = p0(cnt_r, tot_r, 0)
    print(f"  F_ref = {F_ref:.4f}")

    # --- MCM only ---
    print("\n--- MCM_only（Rzなし）---")
    cnt_m, tot_m = run(backend, build_mcm_only(args.qa, args.qt, args.tau),
                       SHOTS, "MCM_only")
    F_mcm = p0(cnt_m, tot_m, 1)
    gps_mcm = F_mcm - F_ref
    print(f"  F_mcm = {F_mcm:.4f}  GPS={gps_mcm:+.4f}")

    # --- EEDT Rz(-phi) ---
    print("\n--- EEDT Rz(-phi) ---")
    cnt_e, tot_e = run(backend, build_eedt(args.qa, args.qt,
                                            args.tau, NU_ZZ_HZ, -1),
                       SHOTS, "EEDT_neg")
    F_eedt = p0(cnt_e, tot_e, 1)
    gps_eedt = F_eedt - F_ref
    print(f"  F_eedt = {F_eedt:.4f}  GPS={gps_eedt:+.4f}")

    # --- 判定 ---
    print()
    print("=" * 52)
    print("  MCMバックアクション分析")
    print("=" * 52)
    print(f"  REF:      F={F_ref:.4f}")
    print(f"  MCM_only: F={F_mcm:.4f}  GPS={gps_mcm:+.4f}")
    print(f"  EEDT-:    F={F_eedt:.4f}  GPS={gps_eedt:+.4f}")
    print()

    mcm_loss = F_ref - F_mcm
    rz_loss  = F_mcm - F_eedt

    print(f"  MCM自体の損失:  {mcm_loss:+.4f}")
    print(f"  Rz補正の損失:   {rz_loss:+.4f}")
    print()

    if mcm_loss > 0.05:
        print("  ★ MCM自体がdephasingを起こしている")
        print("    → IQM GarnetのMCMバックアクションが大きい")
        print("    → feedforward原理的に困難")
        print("    → 論文: MCMバックアクションがEEDTの限界因子")
        verdict = "MCM_DEPHASING"
    elif mcm_loss < 0.01:
        print("  ★ MCM自体の影響は小さい")
        print("    → 問題はRz補正の角度/方向")
        print("    → 正しい角度を探せばGPS>0の可能性")
        verdict = "RZ_PROBLEM"
    else:
        print("  △ MCMとRzの両方が寄与")
        verdict = "MIXED"

    results = {
        "tau_us": args.tau, "phi": float(phi),
        "F_ref": float(F_ref), "F_mcm": float(F_mcm),
        "F_eedt": float(F_eedt),
        "gps_mcm": float(gps_mcm), "gps_eedt": float(gps_eedt),
        "mcm_loss": float(mcm_loss), "rz_loss": float(rz_loss),
        "verdict": verdict
    }
    with open(OUTPUT_DIR / "mcm_test.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[完了] verdict={verdict}")
    print(f"  結果: {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()
