"""
STEP1: ZZ Ramsey測定 — IQM Garnet専用
======================================
目的: IQM GarnetのQB20-QB19間の残留νZZを実測する
     「tunable couplerでもνZZ=0ではない」かを確認

消費: 10条件 × 2000shots = 20,000shots ≈ 2クレジット

実行:
  python step1_zz_ramsey.py --token YOUR_IQM_TOKEN

QB20(index=19) - QB19(index=18):
  T2(QB20) = 25.4us, T2(QB19) = 6.0us
  τ範囲: 0.5〜15µs（短めのT2に合わせる）
"""

import argparse
import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.optimize import curve_fit

# ─── 設定 ─────────────────────────────────────────
QB_ANCILLA  = 19   # QB20（T2最長）
QB_TARGET   = 18   # QB19（隣接）
SHOTS       = 2000

# τリスト: QB19のT2=6µsを考慮して短めに設定
# νZZ=1kHz想定なら1周期=1ms >> T2 → 振幅で判断
# νZZ=10kHz想定なら1周期=100µs >> T2 → 同様
# νZZ=100kHz想定なら1周期=10µs ≈ T2 → 縞が見える可能性
TAU_US_LIST = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 12.0, 15.0]

OUTPUT_DIR  = Path(f"step1_zz_{datetime.now().strftime('%Y%m%d_%H%M%S')}")

# ─── ZZ Ramsey 回路 ───────────────────────────────

def build_zz_circuit(qa: int, qt: int, tau_us: float):
    """
    ZZ Ramsey回路:
    両qubitを|+>に初期化 → τ自由進化 → X基底測定

    ZZ結合があると:
      P(00) + P(11) ≠ P(01) + P(10)
      ZZコリレーター C = P(00)+P(11) - P(01)-P(10) が振動

    ZZ=0なら C=0（独立）
    ZZ≠0なら Cがτに応じて変化
    """
    from qiskit import QuantumCircuit
    nq  = max(qa, qt) + 1
    qc  = QuantumCircuit(nq, 2)

    # |+> 初期化
    qc.h(qa)
    qc.h(qt)

    # τ自由進化
    tau_ns = int(tau_us * 1000)
    qc.delay(tau_ns, qa, unit="ns")
    qc.delay(tau_ns, qt, unit="ns")

    # X基底測定（H → Z測定 = X基底）
    qc.h(qa)
    qc.h(qt)

    qc.measure(qa, 0)  # cbit0 = QB_ANCILLA
    qc.measure(qt, 1)  # cbit1 = QB_TARGET

    return qc

# ─── 実行 ─────────────────────────────────────────

def connect(token: str, device: str):
    try:
        from iqm.qiskit_iqm import IQMProvider
    except ImportError:
        print("[ERROR] pip install 'iqm-client[qiskit]>=33.0,<34.0'")
        sys.exit(1)

    url = f"https://cocos.resonance.meetiqm.com/{device}"
    print(f"[接続] {url}")
    provider = IQMProvider(url, token=token)
    backend  = provider.get_backend()
    print(f"[OK] qubits={backend.num_qubits}  queue={0}")

    # 隣接確認
    cmap = list(backend.coupling_map)
    qa, qt = QB_ANCILLA, QB_TARGET
    adjacent = ([qa, qt] in cmap or [qt, qa] in cmap)
    print(f"[qubit] QB{qa+1}(index={qa}) - QB{qt+1}(index={qt}) "
          f"隣接={'✓' if adjacent else '✗ 要確認'}")
    return backend

def run(backend, circuit, shots: int, label: str):
    from iqm.qiskit_iqm import transpile_to_IQM
    tc     = transpile_to_IQM(circuit, backend)
    job    = backend.run(tc, shots=shots)
    result = job.result()
    counts = result.get_counts()
    total  = sum(counts.values())
    print(f"  [{label}]  shots={total}  "
          f"counts={dict(sorted(counts.items()))}")
    return counts, total

# ─── 解析 ─────────────────────────────────────────

def zz_correlator(counts: dict, total: int) -> float:
    """
    ZZコリレーター:
    C = P(00) + P(11) - P(01) - P(10)

    Qiskitのbitstring: 右端がcbit0(ancilla), 左端がcbit1(target)
    例: '01' → cbit0=1(ancilla=1), cbit1=0(target=0) → '10'状態

    C > 0: ZZ相関あり（同じ状態になりやすい）
    C = 0: 独立（ZZ=0）
    C振動: νZZ×τ で振動
    """
    if total == 0:
        return 0.0

    # Qiskit bitstring: '結果cbit1 結果cbit0' (左が高位)
    p = {k: v/total for k, v in counts.items()}

    # 各状態の確率
    p00 = sum(v for k, v in p.items()
              if len(k) >= 2 and k[-1] == '0' and k[-2] == '0')
    p01 = sum(v for k, v in p.items()
              if len(k) >= 2 and k[-1] == '1' and k[-2] == '0')
    p10 = sum(v for k, v in p.items()
              if len(k) >= 2 and k[-1] == '0' and k[-2] == '1')
    p11 = sum(v for k, v in p.items()
              if len(k) >= 2 and k[-1] == '1' and k[-2] == '1')

    C = p00 + p11 - p01 - p10
    return float(C)

def analyze(data: list):
    """νZZをフィットして判定"""
    taus = np.array([d["tau_us"] for d in data]) * 1e-6
    corr = np.array([d["C"] for d in data])

    print()
    print("  τ(µs)   C(ZZ相関)   解釈")
    print("  " + "-"*40)
    for d in data:
        bar   = "█" * int(abs(d["C"]) * 20)
        sign  = "+" if d["C"] > 0 else "-"
        print(f"  {d['tau_us']:5.1f}   {d['C']:+.4f}    {sign}{bar}")

    # フィット試行
    nu_est = 0.0
    fit_ok = False

    # 振動モデル（減衰込み）
    def model_osc(t, A, nu, phi, offset, T2eff):
        return A * np.exp(-t/T2eff) * np.cos(2*np.pi*nu*t + phi) + offset

    # まず振幅だけ見る（νZZが小さすぎると縞が見えない）
    amp = float(np.max(np.abs(corr)))
    print(f"\n  ZZ相関 最大振幅: {amp:.4f}")

    try:
        popt, pcov = curve_fit(
            model_osc, taus, corr,
            p0   = [0.1, 50e3, 0, 0, 10e-6],
            bounds = (
                [0,    100, -np.pi, -0.3, 1e-6],
                [1.0, 1e6,  np.pi,  0.3, 100e-6]
            ),
            maxfev = 20000
        )
        nu_est = abs(popt[1])
        perr   = np.sqrt(np.diag(pcov))
        fit_ok = True
        print(f"  フィット: νZZ = {nu_est/1e3:.2f} ± "
              f"{perr[1]/1e3:.2f} kHz")
        print(f"           振幅={popt[0]:.4f}  "
              f"T2eff={popt[4]*1e6:.1f}µs")
    except Exception as e:
        print(f"  フィット失敗({e})")
        print(f"  → 振幅 {amp:.4f} から粗推定")
        # 振幅からの粗推定（縞が見えない=低周波or高振幅減衰）
        if amp > 0.05:
            nu_est = 0  # 振幅大 → 振動周期がτ範囲より長い可能性
        else:
            nu_est = 0

    # 判定
    print()
    print("  " + "="*45)
    if amp < 0.02 and not fit_ok:
        print("  ★ 判定: C≈0 → νZZ ≈ 0")
        print("     IQMのtunable couplerがZZを完全抑制")
        print("     → EEDT FAIL（設計通り）")
        verdict = "FAIL_ZZ_ZERO"
    elif fit_ok and nu_est > 10e3:
        print(f"  ★ 判定: νZZ = {nu_est/1e3:.1f} kHz → ✓✓ EEDT動作可能")
        print(f"     T2=25.4µsでのτ* = "
              f"{25.4e-6/(2*np.pi*nu_est*25.4e-6+1)*1e6:.1f} µs")
        verdict = "PASS"
    elif fit_ok and nu_est > 1e3:
        print(f"  ★ 判定: νZZ = {nu_est/1e3:.1f} kHz → △ MARGINAL")
        print(f"     動作可能性あり、τ*を短くすれば届くかも")
        verdict = "MARGINAL"
    elif amp > 0.05:
        print(f"  ★ 判定: 相関あり(amp={amp:.3f})だが周波数不明")
        print(f"     τ範囲を広げて再測定推奨")
        verdict = "UNKNOWN_AMP"
    else:
        print(f"  ★ 判定: νZZ < 1kHz（検出限界以下）")
        print(f"     EEDT動作困難")
        verdict = "FAIL_LOW_ZZ"
    print("  " + "="*45)

    return nu_est, verdict, amp

def make_plot(data: list, nu_est: float, verdict: str):
    taus = np.array([d["tau_us"] for d in data])
    corr = np.array([d["C"] for d in data])
    se   = np.array([d["SE"] for d in data])

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(
        f"IQM Garnet ZZ Ramsey — QB20-QB19\n"
        f"νZZ={'%.1f kHz' % (nu_est/1e3) if nu_est > 0 else '≈0'}  "
        f"verdict={verdict}",
        fontsize=12)

    # Left: ZZ correlator vs τ
    ax = axes[0]
    ax.errorbar(taus, corr, yerr=se, fmt='o', color='teal',
                ms=8, capsize=4, label='ZZ correlator C')
    ax.axhline(0, color='gray', lw=0.8)
    if nu_est > 1e3:
        tau_fit = np.linspace(0, max(taus), 300) * 1e-6
        T2min   = 6e-6
        fit_curve = 0.1 * np.exp(-tau_fit/T2min) * \
                    np.cos(2*np.pi*nu_est*tau_fit)
        ax.plot(tau_fit*1e6, fit_curve, 'r--', lw=1.5,
                label=f'理論 νZZ={nu_est/1e3:.1f}kHz')
    ax.set_xlabel("τ (µs)", fontsize=11)
    ax.set_ylabel("C = P(00)+P(11)-P(01)-P(10)", fontsize=10)
    ax.set_title("ZZ Correlator", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # Right: |C| vs τ (振幅の単調減少確認)
    ax = axes[1]
    ax.plot(taus, np.abs(corr), 'o-', color='coral', ms=8)
    ax.axhline(0.02, color='gray', lw=1, ls='--',
               label='検出限界 0.02')
    ax.fill_between(taus, 0, 0.02, alpha=0.1, color='gray')
    ax.set_xlabel("τ (µs)", fontsize=11)
    ax.set_ylabel("|C|", fontsize=11)
    ax.set_title("ZZ振幅", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    out = OUTPUT_DIR / "zz_ramsey_result.png"
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"\n[プロット] {out}")

# ─── メイン ───────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="STEP1: ZZ Ramsey on IQM Garnet")
    p.add_argument("--token",  required=True, help="IQM API token")
    p.add_argument("--device", default="garnet",
                   help="garnet / garnet:mock")
    p.add_argument("--qa", type=int, default=QB_ANCILLA,
                   help=f"ancilla qubit index (default={QB_ANCILLA}=QB20)")
    p.add_argument("--qt", type=int, default=QB_TARGET,
                   help=f"target qubit index (default={QB_TARGET}=QB19)")
    args = p.parse_args()

    OUTPUT_DIR.mkdir(exist_ok=True)

    print("=" * 50)
    print("  STEP1: ZZ Ramsey測定")
    print(f"  QB{args.qa+1}(index={args.qa}) - "
          f"QB{args.qt+1}(index={args.qt})")
    print(f"  {len(TAU_US_LIST)}条件 × {SHOTS}shots = "
          f"{len(TAU_US_LIST)*SHOTS:,}shots")
    print(f"  推定消費: ~2クレジット")
    print("=" * 50)

    backend = connect(args.token, args.device)

    data = []
    print(f"\n測定開始 ({datetime.now().strftime('%H:%M:%S')})")

    for tau_us in TAU_US_LIST:
        qc = build_zz_circuit(args.qa, args.qt, tau_us)
        counts, total = run(backend, qc, SHOTS, f"τ={tau_us}µs")

        C  = zz_correlator(counts, total)
        SE = 2.0 / np.sqrt(total)  # 保守的な誤差推定

        data.append({
            "tau_us": tau_us,
            "C":      float(C),
            "SE":     float(SE),
            "counts": counts,
            "total":  total
        })

    # 保存
    with open(OUTPUT_DIR / "zz_raw.json", "w") as f:
        json.dump(data, f, indent=2)

    # 解析・判定
    print("\n" + "=" * 50)
    print("  解析結果")
    print("=" * 50)
    nu_est, verdict, amp = analyze(data)

    # プロット
    make_plot(data, nu_est, verdict)

    # Table 2用テキスト
    print()
    print("  [Table 2 記録用]")
    print(f"  IQM Garnet QB20-QB19 | "
          f"νZZ={'%.2f kHz' % (nu_est/1e3) if nu_est > 100 else '<0.1 kHz'} | "
          f"T2=25.4µs | 判定={verdict}")
    print()
    print(f"  結果ファイル: {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()
