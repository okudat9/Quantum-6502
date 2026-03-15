"""
EEDT 符号スキャン — IQM Garnet QB12-QB17
==========================================
3回路 × 1000shots で符号を確定する

REF:    補正なし（ベースライン）
EEDT-:  無条件 Rz(-phi) 全ショット
EEDT+:  無条件 Rz(+phi) 全ショット

EEDT- > REF → 正しい補正は Rz(-phi)
EEDT+ > REF → 正しい補正は Rz(+phi)
どちらも < REF → 無条件補正では不可、条件付き必要

消費: 3回路 × 1000shots ≈ 3〜5 credits

実行:
  python sign_scan.py --token YOUR_TOKEN
"""

import argparse, json, sys, os
import numpy as np
from datetime import datetime
from pathlib import Path

QB_ANCILLA = 11   # QB12（T2=10.35µs）
QB_TARGET  = 16   # QB17（T2=13.04µs）
SHOTS      = 1000
NU_ZZ_HZ   = 35.86e3  # 実測値
TAU_US     = 3.0       # τ*付近

OUTPUT_DIR = Path(f"sign_scan_{datetime.now().strftime('%Y%m%d_%H%M%S')}")

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
    res = job.result()
    counts = res.get_counts()
    total  = sum(counts.values())
    print(f"  [{label}] {dict(sorted(counts.items()))}")
    return counts, total

def p0(counts, total, bit_index=0):
    """bit_index番目のcbitが0の確率"""
    if total == 0: return 0.5
    pos = -(bit_index + 1)
    cnt = sum(v for k,v in counts.items()
              if len(k) > bit_index and k[pos] == '0')
    return cnt / total

def build_ref(qt, tau_us):
    """参照回路: 補正なし"""
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(qt + 1, 1)
    qc.h(qt)
    qc.delay(int(tau_us * 1000), qt, unit="ns")
    qc.h(qt)
    qc.measure(qt, 0)
    return qc

def build_eedt_uncond(qa, qt, tau_us, nu_zz_hz, rz_sign):
    """
    無条件フィードフォワード:
    rz_sign=+1: Rz(+phi) を全ショットに適用
    rz_sign=-1: Rz(-phi) を全ショットに適用
    """
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
    print("  EEDT 符号スキャン")
    print(f"  QB{args.qa+1}(ancilla) → QB{args.qt+1}(target)")
    print(f"  tau={args.tau}µs  phi={phi:.4f} rad")
    print(f"  3回路 × {SHOTS}shots  推定~3〜5 credits")
    print("=" * 52)

    backend = connect(args.token, args.device)

    results = {}

    # --- REF ---
    print("\n--- REF（補正なし）---")
    qc_ref = build_ref(args.qt, args.tau)
    cnt_r, tot_r = run(backend, qc_ref, SHOTS, "REF")
    F_ref = p0(cnt_r, tot_r, 0)
    se_ref = np.sqrt(F_ref*(1-F_ref)/SHOTS)
    print(f"  F_ref = {F_ref:.4f} ± {se_ref:.4f}")
    results['ref'] = {'F': float(F_ref), 'SE': float(se_ref)}

    # --- EEDT- : Rz(-phi) ---
    print("\n--- EEDT- : 無条件 Rz(-phi) ---")
    qc_neg = build_eedt_uncond(args.qa, args.qt, args.tau, NU_ZZ_HZ, -1)
    cnt_n, tot_n = run(backend, qc_neg, SHOTS, "EEDT-")
    F_neg = p0(cnt_n, tot_n, 1)
    se_neg = np.sqrt(F_neg*(1-F_neg)/SHOTS)
    gps_neg = F_neg - F_ref
    print(f"  F_eedt = {F_neg:.4f} ± {se_neg:.4f}")
    print(f"  GPS    = {gps_neg:+.4f}")
    results['eedt_neg'] = {'F': float(F_neg), 'SE': float(se_neg),
                            'GPS': float(gps_neg)}

    # --- EEDT+ : Rz(+phi) ---
    print("\n--- EEDT+ : 無条件 Rz(+phi) ---")
    qc_pos = build_eedt_uncond(args.qa, args.qt, args.tau, NU_ZZ_HZ, +1)
    cnt_p, tot_p = run(backend, qc_pos, SHOTS, "EEDT+")
    F_pos = p0(cnt_p, tot_p, 1)
    se_pos = np.sqrt(F_pos*(1-F_pos)/SHOTS)
    gps_pos = F_pos - F_ref
    print(f"  F_eedt = {F_pos:.4f} ± {se_pos:.4f}")
    print(f"  GPS    = {gps_pos:+.4f}")
    results['eedt_pos'] = {'F': float(F_pos), 'SE': float(se_pos),
                            'GPS': float(gps_pos)}

    # --- 判定 ---
    print()
    print("=" * 52)
    print("  符号判定")
    print("=" * 52)
    print(f"  REF:   F={F_ref:.4f}")
    print(f"  Rz(-): F={F_neg:.4f}  GPS={gps_neg:+.4f}")
    print(f"  Rz(+): F={F_pos:.4f}  GPS={gps_pos:+.4f}")
    print()

    if gps_neg > gps_pos and gps_neg > 0.01:
        verdict = "SIGN_NEGATIVE"
        print("  ★ Rz(-phi) が有効 → 正しい補正はRz(-phi)")
        print("    次: cc_prx で ancilla=0→Rz(-phi), ancilla=1→Rz(+phi)")
    elif gps_pos > gps_neg and gps_pos > 0.01:
        verdict = "SIGN_POSITIVE"
        print("  ★ Rz(+phi) が有効 → 正しい補正はRz(+phi)")
        print("    次: cc_prx で ancilla=0→Rz(+phi), ancilla=1→Rz(-phi)")
    elif max(gps_neg, gps_pos) > 0:
        verdict = "WEAK_POSITIVE"
        print("  △ どちらも正だが弱い → shotsを増やして確認")
    else:
        verdict = "BOTH_NEGATIVE"
        print("  ✗ 両方負 → 無条件補正では不可")
        print("    ZZがdephasingのため位相補正が無意味な可能性")

    print("=" * 52)

    results['verdict'] = verdict
    results['tau_us']  = args.tau
    results['phi']     = float(phi)
    results['nu_zz']   = NU_ZZ_HZ

    with open(OUTPUT_DIR / "sign_scan.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[完了] {OUTPUT_DIR}/sign_scan.json")

if __name__ == "__main__":
    main()
