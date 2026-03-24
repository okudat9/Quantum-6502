#!/usr/bin/env python3
"""
calib_monitor.py — IBM キャリブデータ収集スクリプト
GitHub Actions から定期実行される（量子時間ゼロ）

取得データ:
  T2, T1, CZ error, RO error → t2_history.csv に追記
"""

import os, csv, datetime
from pathlib import Path

TOKEN        = os.environ["IBM_QUANTUM_TOKEN"]
BACKEND_NAME = os.environ.get("BACKEND_NAME", "ibm_marrakesh")
TARGET_Q     = int(os.environ.get("TARGET_Q", "95"))
ANCILLA_Q    = int(os.environ.get("ANCILLA_Q", "94"))
CSV_PATH     = Path("data/t2_history.csv")

def main():
    from qiskit_ibm_runtime import QiskitRuntimeService

    print(f"[calib_monitor] {BACKEND_NAME} Q{TARGET_Q}+Q{ANCILLA_Q}")
    svc     = QiskitRuntimeService(
        channel="ibm_quantum_platform", token=TOKEN)
    backend = svc.backend(BACKEND_NAME)

    # ── キャリブ取得
    try:
        t2_tgt = backend.qubit_properties(TARGET_Q).t2  * 1e6
        t2_anc = backend.qubit_properties(ANCILLA_Q).t2 * 1e6
        t1_tgt = backend.qubit_properties(TARGET_Q).t1  * 1e6
    except Exception as e:
        print(f"[ERROR] qubit_properties: {e}"); return

    try:
        cz = backend.target['cz'][(TARGET_Q, ANCILLA_Q)].error
    except Exception:
        try:
            cz = backend.properties().gate_error(
                'cz', [TARGET_Q, ANCILLA_Q])
        except Exception:
            cz = None

    try:
        ro_tgt = backend.properties().readout_error(TARGET_Q)
        ro_anc = backend.properties().readout_error(ANCILLA_Q)
    except Exception:
        ro_tgt = ro_anc = None

    # ── ν_ZZ回帰予測（T2から）
    nu_pred = round(0.013053 * t2_tgt + 0.474, 3) if t2_tgt else None

    # ── 前回値との比較
    t2_chg_pct = None
    if CSV_PATH.exists():
        with open(CSV_PATH, 'r') as f:
            rows = list(csv.DictReader(f))
        if rows:
            prev_t2 = float(rows[-1].get("t2_tgt_us", 0) or 0)
            if prev_t2:
                t2_chg_pct = round(
                    (t2_tgt - prev_t2) / prev_t2 * 100, 2)

    # ── ドリフトアラート判定
    alert = ""
    if t2_chg_pct is not None:
        if abs(t2_chg_pct) >= 17.0:
            alert = "DIRECT_TRIGGER"
        elif abs(t2_chg_pct) >= 10.0:
            alert = "WATCH"

    # ── CSV追記
    row = {
        "timestamp":    datetime.datetime.utcnow().strftime(
                            "%Y-%m-%dT%H:%M:%SZ"),
        "backend":      BACKEND_NAME,
        "target_q":     TARGET_Q,
        "ancilla_q":    ANCILLA_Q,
        "t2_tgt_us":    round(t2_tgt, 1)  if t2_tgt  else "",
        "t2_anc_us":    round(t2_anc, 1)  if t2_anc  else "",
        "t1_tgt_us":    round(t1_tgt, 1)  if t1_tgt  else "",
        "cz_error":     round(cz, 7)      if cz      else "",
        "ro_tgt":       round(ro_tgt, 5)  if ro_tgt  else "",
        "ro_anc":       round(ro_anc, 5)  if ro_anc  else "",
        "nu_zz_pred":   nu_pred           if nu_pred else "",
        "t2_chg_pct":   t2_chg_pct        if t2_chg_pct is not None else "",
        "alert":        alert,
    }

    CSV_PATH.parent.mkdir(exist_ok=True)
    write_header = not CSV_PATH.exists()
    with open(CSV_PATH, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=row.keys())
        if write_header:
            w.writeheader()
        w.writerow(row)

    # ── 出力
    print(f"  T2(Q{TARGET_Q})  = {t2_tgt:.1f} µs"
          f"  ({t2_chg_pct:+.1f}%)" if t2_chg_pct else "")
    print(f"  T2(Q{ANCILLA_Q})  = {t2_anc:.1f} µs")
    print(f"  CZ        = {cz:.6f}" if cz else "  CZ = N/A")
    print(f"  nu_pred   = {nu_pred} kHz")
    if alert:
        print(f"  ⚠️  ALERT: {alert} (T2変化 {t2_chg_pct:+.1f}%)")
    print(f"  → 保存: {CSV_PATH}")

if __name__ == "__main__":
    main()
