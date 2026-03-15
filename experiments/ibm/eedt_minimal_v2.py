"""
EEDT 最小実験セット — IQM / IBM 統合版 (修正済み v2)
=====================================================
13条件 × 2000 shots = 26,000 shots
無料枠で動作可能

【インストール】Windows PowerShell / コマンドプロンプト:

  # 仮想環境の作成（推奨）
  python -m venv eedt_env
  eedt_env/Scripts/activate  # PowerShell: eedt_env/Scripts/Activate.ps1

  # IQM用（Resonanceダッシュボードの COMPATIBLE CLIENT LIBRARY VERSIONS を確認）
  pip install "iqm-client[qiskit]>=33.0,<34.0"
  pip install numpy matplotlib scipy

  # IBM用（追加で必要な場合）
  pip install qiskit-ibm-runtime

  注意: qiskit-iqm は廃止済み。インストールしないこと。
  もしインストール済みなら: pip uninstall qiskit-iqm

【使い方】
  # IQM Garnet（qubit自動選択）
  python eedt_minimal.py --backend iqm --token YOUR_IQM_TOKEN --device garnet

  # IQM Garnet（qubit手動指定）
  python eedt_minimal.py --backend iqm --token YOUR_IQM_TOKEN --device garnet --qa 2 --qt 5

  # IQM mock（無料テスト）
  python eedt_minimal.py --backend iqm --token YOUR_IQM_TOKEN --device garnet:mock

  # IBM（参考）
  python eedt_minimal.py --backend ibm --token YOUR_IBM_TOKEN --device ibm_marrakesh --qa 94 --qt 95

【IQM API Tokenの取得】
  https://resonance.meetiqm.com/ → Dashboard → Generate token
"""

import argparse
import json
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from scipy.optimize import curve_fit

# ═══════════════════════════════════════════════════════════
#  設定
# ═══════════════════════════════════════════════════════════
SHOTS       = 1000
NU_ZZ_REF   = 3.6e3    # Hz（ibm_marrakesh参照値、STEP1で実測更新）
T_MCM_NS    = 700      # ns（MCMデッドタイム）
T2_PASS_US  = 200      # µs（動作判定閾値）

TAU_ZZ_US   = [10, 20, 40, 60, 80, 100, 120, 150, 180, 200]  # STEP1 (skipped)
TAU_T2_US   = [10, 30, 60, 100, 150, 200, 250, 300]           # STEP2
TAU_SCAN_US = [60, 80]                                         # STEP3 continuation only
N_LIST      = [1, 2, 3, 4, 5]                                  # STEP4
TAU_NSCAN_US = 40  # STEP4固定τ（実測T2で再計算）

OUTPUT_DIR = Path(f"eedt_{datetime.now().strftime('%Y%m%d_%H%M%S')}")


# ═══════════════════════════════════════════════════════════
#  バックエンド接続
# ═══════════════════════════════════════════════════════════

def connect_iqm(token: str, device: str) -> tuple:
    """
    IQM Resonance 接続
    戻り値: (backend, qa, qt) — coupling_map から隣接qubitペアを自動選択
    """
    try:
        from iqm.qiskit_iqm import IQMProvider
    except ImportError:
        print("[ERROR] iqm-client[qiskit] がインストールされていません")
        print("  pip install 'iqm-client[qiskit]>=33.0,<34.0'")
        sys.exit(1)

    url = "https://resonance.meetiqm.com/"
    print(f"[IQM] 接続中: {url}  quantum_computer={device}")

    try:
        provider = IQMProvider(url, quantum_computer=device, token=token)
        backend = provider.get_backend()
    except Exception as e:
        print(f"[ERROR] 接続失敗: {e}")
        print("  → TOKEN と device名を確認してください")
        print("  → mock試験: --device garnet:mock")
        sys.exit(1)

    n_qubits = backend.num_qubits
    cmap = list(backend.coupling_map)
    print(f"[IQM] {device} 接続成功  qubits={n_qubits}")
    print(f"      coupling_map (一部): {cmap[:8]}...")

    # T2情報はIQM Resonanceカレンダー/ダッシュボードで確認
    # backend.properties() はIQMでは非対応のため使用しない
    print()
    print("  [!] T2/T1データはResonanceダッシュボードで確認してください")
    print(f"      https://resonance.meetiqm.com/")

    return backend, cmap


def auto_select_qubits(cmap: list, qa_hint: int, qt_hint: int) -> tuple:
    """
    coupling_map から隣接qubitペアを選択する
    ヒントが隣接していれば採用、そうでなければ最初のペアを使用
    """
    edges = [(a, b) for a, b in cmap]
    # ヒントが隣接しているか確認
    if [qa_hint, qt_hint] in cmap or [qt_hint, qa_hint] in cmap:
        print(f"[qubit] 指定ペア Q{qa_hint}-Q{qt_hint} は隣接しています ✓")
        return qa_hint, qt_hint

    # 最初の隣接ペアを採用
    if edges:
        qa, qt = edges[0]
        print(f"[qubit] 指定ペアQ{qa_hint}-Q{qt_hint}は非隣接 → "
              f"自動選択: Q{qa}-Q{qt}")
        return qa, qt

    # フォールバック
    print(f"[qubit] coupling_map 空 → デフォルト Q{qa_hint}-Q{qt_hint}")
    return qa_hint, qt_hint


def connect_ibm(token: str, device: str, instance: str) -> object:
    """IBM Quantum 接続"""
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
    except ImportError:
        print("[ERROR] qiskit-ibm-runtime がインストールされていません")
        print("  pip install qiskit-ibm-runtime")
        sys.exit(1)

    service = QiskitRuntimeService(
        channel="ibm_quantum", token=token, instance=instance)
    backend = service.backend(device)
    print(f"[IBM] {device} 接続成功")
    return backend


# ═══════════════════════════════════════════════════════════
#  回路実行（IQM / IBM 共通ラッパー）
# ═══════════════════════════════════════════════════════════

def run_iqm(backend, circuit, shots: int, label: str = "") -> tuple:
    """
    IQM専用実行ラッパー
    transpile_to_IQM を使用（Star アーキテクチャ対応）
    """
    from qiskit import transpile
    try:
        from iqm.qiskit_iqm import transpile_to_IQM
        tc = transpile_to_IQM(circuit, backend)
    except Exception:
        # フォールバック: 標準transpile
        tc = transpile(circuit, backend, optimization_level=1)

    job = backend.run(tc, shots=shots)
    result = job.result()
    counts = result.get_counts()
    total = sum(counts.values())
    print(f"    [{label}] shots={total}  "
          f"states={len(counts)}  "
          f"top={list(counts.items())[:3]}")
    return counts, total


def run_ibm(backend, circuit, shots: int, label: str = "") -> tuple:
    """IBM専用実行ラッパー（SamplerV2使用）"""
    from qiskit import transpile
    from qiskit_ibm_runtime import SamplerV2 as Sampler
    tc = transpile(circuit, backend, optimization_level=1)
    sampler = Sampler(mode=backend)
    job = sampler.run([tc], shots=shots)
    result = job.result()
    # SamplerV2 の結果形式
    counts = result[0].data.c.get_counts()
    total = sum(counts.values())
    print(f"    [{label}] shots={total}  states={len(counts)}")
    return counts, total


def run_job(backend, circuit, shots: int, label: str,
            use_iqm: bool = True) -> tuple:
    """バックエンド自動判定ラッパー"""
    if use_iqm:
        return run_iqm(backend, circuit, shots, label)
    else:
        return run_ibm(backend, circuit, shots, label)


# ═══════════════════════════════════════════════════════════
#  統計ヘルパー
# ═══════════════════════════════════════════════════════════

def p_zero(counts: dict, total: int, bit_index: int = 0) -> float:
    """
    指定した古典ビット (bit_index) が '0' の確率を返す

    Qiskit の bitstring は右端が古典レジスタの最低位ビット (cbit 0)
    例: '01' → cbit0=1, cbit1=0

    bit_index=0 → counts の各文字列の[-1]が '0' の確率
    bit_index=1 → counts の各文字列の[-2]が '0' の確率
    """
    if total == 0:
        return 0.5
    pos = -(bit_index + 1)  # 右端からのオフセット
    cnt = sum(v for k, v in counts.items()
              if len(k) >= (bit_index + 1) and k[pos] == "0")
    return cnt / total


def gps_and_se(counts_eedt: dict, total_e: int,
               counts_ref: dict, total_r: int,
               bit_index: int = 0) -> tuple:
    """GPS = F_EEDT - F_ref と標準誤差を返す"""
    p_e = p_zero(counts_eedt, total_e, bit_index)
    p_r = p_zero(counts_ref, total_r, bit_index)
    se = np.sqrt(p_e * (1 - p_e) / max(total_e, 1) +
                 p_r * (1 - p_r) / max(total_r, 1))
    return p_e - p_r, se


# ═══════════════════════════════════════════════════════════
#  回路ビルダー
# ═══════════════════════════════════════════════════════════

def _n_qubits(qa: int, qt: int) -> int:
    """回路に必要な最小 qubit 数"""
    return max(qa, qt) + 1


def build_zz_ramsey(qa: int, qt: int, tau_us: float):
    """
    ZZ Ramseyスペクトロスコピー回路
    両qubitを|+>に初期化してτ待機後に測定
    ZZ結合があると、測定結果に相関が現れる
    """
    from qiskit import QuantumCircuit
    nq = _n_qubits(qa, qt)
    qc = QuantumCircuit(nq, 2)
    qc.h(qa)
    qc.h(qt)
    tau_ns = int(tau_us * 1000)
    qc.delay(tau_ns, qa, unit="ns")
    qc.delay(tau_ns, qt, unit="ns")
    qc.h(qa)
    qc.h(qt)
    qc.measure(qa, 0)
    qc.measure(qt, 1)
    return qc


def build_t2_ramsey(qt: int, tau_us: float, detuning_hz: float = 300.0):
    """
    T₂ Ramsey 回路
    小さいデチューニングで縞を作り、減衰からT₂を抽出
    """
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(qt + 1, 1)
    qc.h(qt)
    tau_ns = int(tau_us * 1000)
    qc.delay(tau_ns, qt, unit="ns")
    # デチューニング相当の位相回転
    phase = 2 * np.pi * detuning_hz * tau_us * 1e-6
    qc.rz(phase, qt)
    qc.h(qt)
    qc.measure(qt, 0)
    return qc


def build_eedt(qa: int, qt: int, tau_us: float,
               nu_zz_hz: float, n_mcm: int = 1):
    """
    EEDT 回路: N回 MCM + フィードフォワード
    IQMでは条件分岐ゲートが非対応のため、無条件Rz補正（近似版）を使用
    完全版はIf-Else動的回路が使える環境が必要
    """
    from qiskit import QuantumCircuit
    nq = _n_qubits(qa, qt)
    n_cbits = n_mcm + 1  # MCM用 + target測定用
    qc = QuantumCircuit(nq, n_cbits)

    qc.h(qa)
    qc.h(qt)

    dt_us = tau_us / n_mcm
    omega_zz = 2 * np.pi * nu_zz_hz

    for i in range(n_mcm):
        dt_ns = int(dt_us * 1000)
        qc.delay(dt_ns, qa, unit="ns")
        qc.delay(dt_ns, qt, unit="ns")
        # Mid-circuit measurement（ancilla qubit）
        qc.measure(qa, i)
        # フィードフォワード（近似: 期待値補正）
        # 完全版: c_if(creg[i], 1) で回転量を切り替える
        theta = omega_zz * dt_us * 1e-6 * 0.5
        qc.rz(theta, qt)

    qc.h(qt)
    qc.measure(qt, n_mcm)
    return qc


def build_ref(qt: int, tau_us: float):
    """参照回路（フィードフォワードなし）"""
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(qt + 1, 1)
    qc.h(qt)
    tau_ns = int(tau_us * 1000)
    qc.delay(tau_ns, qt, unit="ns")
    qc.h(qt)
    qc.measure(qt, 0)
    return qc


# ═══════════════════════════════════════════════════════════
#  STEP 1: ZZ Ramsey
# ═══════════════════════════════════════════════════════════

def step1_zz_ramsey(backend, qa: int, qt: int,
                    use_iqm: bool) -> float:
    print("\n" + "─" * 55)
    print("STEP 1: ZZ Ramsey測定（νZZ確認）")
    print("─" * 55)

    data = []
    for tau_us in TAU_ZZ_US:
        qc = build_zz_ramsey(qa, qt, tau_us)
        counts, total = run_job(
            backend, qc, SHOTS, f"ZZ τ={tau_us}µs", use_iqm)
        if total == 0:
            continue

        # P(qt=0): cbit 1（qt が measure → cbit 1 に記録）
        p_qt0 = p_zero(counts, total, bit_index=1)
        # P(qa=0): cbit 0
        p_qa0 = p_zero(counts, total, bit_index=0)

        # ZZシグネチャ: |P(qt=0) - 0.5| × 2（理想ZZ=0なら0）
        zz_sig = abs(p_qt0 - 0.5) * 2

        data.append({
            "tau_us": tau_us,
            "p_qt0": float(p_qt0),
            "p_qa0": float(p_qa0),
            "zz_sig": float(zz_sig)
        })
        print(f"  τ={tau_us:3d}µs: P(qt=0)={p_qt0:.4f}  "
              f"P(qa=0)={p_qa0:.4f}  |ZZsig|={zz_sig:.4f}")

    _save("step1_zz_ramsey", data)
    nu_zz = _analyze_zz(data)
    return nu_zz


def _analyze_zz(data: list) -> float:
    if len(data) < 4:
        print("  [ZZ判定] データ不足")
        return 0.0

    taus = np.array([d["tau_us"] for d in data]) * 1e-6
    p0s  = np.array([d["p_qt0"] for d in data])

    nu_est = 0.0
    try:
        def model(t, A, nu, phi, offset):
            return A * np.cos(2 * np.pi * nu * t + phi) + offset
        popt, _ = curve_fit(
            model, taus, p0s,
            p0=[0.1, 2e3, 0, 0.5],
            bounds=([0, 10, -np.pi, 0.3], [0.5, 50e3, np.pi, 0.7]),
            maxfev=10000)
        nu_est = abs(popt[1])
        amp = abs(popt[0])
        print(f"\n  [ZZ推定] νZZ ≈ {nu_est/1e3:.3f} kHz  "
              f"振幅={amp:.4f}")
    except Exception as e:
        amp = float(np.max(np.abs(p0s - 0.5)))
        nu_est = amp * 8e3  # 粗い推定
        print(f"\n  [ZZ推定] フィット失敗({e})  "
              f"max|P-0.5|={amp:.4f}  推定νZZ≈{nu_est/1e3:.2f}kHz")

    print()
    if nu_est > 500:
        print(f"  ★ ZZ判定: νZZ={nu_est/1e3:.2f}kHz > 0.5kHz → ✓ EEDT候補")
    elif nu_est > 100:
        print(f"  ★ ZZ判定: νZZ={nu_est/1e3:.2f}kHz → △ MARGINAL")
    else:
        print(f"  ★ ZZ判定: νZZ≈0 → ✗ FAIL")
        print("     IQMのtunable couplerがZZ結合を設計上抑制している")
        print("     → Table 2にIQM FAIL（理論通り）として記録")
    return nu_est


# ═══════════════════════════════════════════════════════════
#  STEP 2: T₂ Ramsey
# ═══════════════════════════════════════════════════════════

def step2_t2(backend, qt: int, use_iqm: bool) -> float:
    print("\n" + "─" * 55)
    print("STEP 2: T₂ Ramsey測定")
    print("─" * 55)

    data = []
    for tau_us in TAU_T2_US:
        qc = build_t2_ramsey(qt, tau_us, detuning_hz=300)
        counts, total = run_job(
            backend, qc, SHOTS, f"T2 τ={tau_us}µs", use_iqm)
        p0 = p_zero(counts, total, bit_index=0)
        data.append({"tau_us": tau_us, "p0": float(p0)})
        print(f"  τ={tau_us:3d}µs: P(0)={p0:.4f}")

    _save("step2_t2", data)
    T2 = _analyze_t2(data)
    return T2


def _analyze_t2(data: list) -> float:
    if len(data) < 4:
        return 0.0

    taus = np.array([d["tau_us"] for d in data]) * 1e-6
    p0s  = np.array([d["p0"] for d in data])

    T2_est = 0.0
    try:
        def model(t, A, T2, f, phi):
            return A * np.exp(-t / T2) * np.cos(2 * np.pi * f * t + phi) + 0.5
        popt, _ = curve_fit(
            model, taus, p0s,
            p0=[0.4, 200e-6, 300, 0],
            bounds=([0, 1e-6, 0, -np.pi], [0.5, 2e-3, 5e3, np.pi]),
            maxfev=10000)
        T2_est = popt[1]
        print(f"\n  [T₂推定] T₂ ≈ {T2_est*1e6:.1f} µs")
    except Exception as e:
        print(f"\n  [T₂推定] フィット失敗({e})")

    print()
    if T2_est > T2_PASS_US * 1e-6:
        tau_star = T2_est / (2 * np.pi * NU_ZZ_REF * T2_est + 1)
        print(f"  ★ T₂判定: {T2_est*1e6:.0f}µs > {T2_PASS_US}µs → ✓ PASS")
        print(f"     予測τ*(νZZ=3.6kHz) ≈ {tau_star*1e6:.1f} µs")
    elif T2_est > 150e-6:
        print(f"  ★ T₂判定: {T2_est*1e6:.0f}µs → △ MARGINAL")
    elif T2_est > 0:
        print(f"  ★ T₂判定: {T2_est*1e6:.0f}µs < 150µs → ✗ FAIL")
    return T2_est


# ═══════════════════════════════════════════════════════════
#  STEP 3: τ-scan
# ═══════════════════════════════════════════════════════════

def step3_tau_scan(backend, qa: int, qt: int,
                   nu_zz_hz: float, use_iqm: bool) -> list:
    print("\n" + "─" * 55)
    print("STEP 3: τ-scan (N=1固定、φe^{-φ}山型確認)")
    print("─" * 55)

    # 参照回路を先に実行
    ref_p0 = {}
    for tau_us in TAU_SCAN_US:
        qc_ref = build_ref(qt, tau_us)
        counts_r, total_r = run_job(
            backend, qc_ref, SHOTS, f"REF τ={tau_us}µs", use_iqm)
        ref_p0[tau_us] = p_zero(counts_r, total_r, bit_index=0)

    data = []
    for tau_us in TAU_SCAN_US:
        qc = build_eedt(qa, qt, tau_us, nu_zz_hz, n_mcm=1)
        counts_e, total_e = run_job(
            backend, qc, SHOTS, f"EEDT τ={tau_us}µs", use_iqm)

        F_eedt = p_zero(counts_e, total_e, bit_index=1)  # target = cbit 1
        F_ref  = ref_p0[tau_us]
        gps    = F_eedt - F_ref
        se     = np.sqrt(F_eedt*(1-F_eedt)/SHOTS + F_ref*(1-F_ref)/SHOTS)
        phi    = 2 * np.pi * nu_zz_hz * tau_us * 1e-6

        data.append({
            "tau_us": tau_us, "phi": float(phi),
            "F_eedt": float(F_eedt), "F_ref": float(F_ref),
            "GPS": float(gps), "SE": float(se)
        })
        sign = "✓" if gps > 0 else "✗"
        print(f"  τ={tau_us:2d}µs  φ={phi:.3f}  "
              f"GPS={gps:+.4f}±{se:.4f}  {sign}")

    _save("step3_tau_scan", data)
    _analyze_tau_scan(data)
    return data


def _analyze_tau_scan(data: list):
    gps = {d["tau_us"]: d["GPS"] for d in data}
    peak = max((gps.get(t, -999) for t in [30, 40, 50]), default=-999)

    print()
    neg20   = gps.get(20, 0) < 0
    pos_mid = peak > 0.005
    decay   = gps.get(60, 999) < gps.get(40, -999)

    print(f"  τ=20µs GPS: {gps.get(20,0):+.4f}  "
          f"{'負 ✓' if neg20 else '正（想定外）'}")
    print(f"  τ=30-50µs peak: {peak:+.4f}  "
          f"{'正 ✓' if pos_mid else 'FAIL'}")
    print(f"  τ=60µs GPS: {gps.get(60,0):+.4f}  "
          f"{'減衰 ✓' if decay else '非単調'}")
    print()
    if neg20 and pos_mid and decay:
        print("  ★ τ-scan判定: 山型確認 → ✓✓✓ EEDT動作確定")
        print("     φe^{-φ} シグネチャが ibm_marrakesh と一致")
    elif pos_mid:
        print("  ★ τ-scan判定: ピーク正あり → △ 部分動作")
    else:
        print("  ★ τ-scan判定: 山型なし → ✗ FAIL")


# ═══════════════════════════════════════════════════════════
#  STEP 4: N-scan
# ═══════════════════════════════════════════════════════════

def step4_n_scan(backend, qa: int, qt: int,
                 nu_zz_hz: float, T2: float,
                 use_iqm: bool) -> list:
    # τ*を理論値から計算
    if nu_zz_hz > 0 and T2 > 0:
        tau_star = T2 / (2 * np.pi * nu_zz_hz * T2 + 1)
        tau_us = float(np.clip(round(tau_star * 1e6 / 10) * 10, 30, 50))
    else:
        tau_us = float(TAU_NSCAN_US)

    print("\n" + "─" * 55)
    print(f"STEP 4: N-scan (τ={tau_us:.0f}µs固定、符号反転確認)")
    print("─" * 55)

    # Ref（N=0）
    qc_ref = build_ref(qt, tau_us)
    counts_r, total_r = run_job(
        backend, qc_ref, SHOTS, "REF", use_iqm)
    F_ref = p_zero(counts_r, total_r, bit_index=0)

    data = []
    for N in N_LIST:
        qc = build_eedt(qa, qt, tau_us, nu_zz_hz, n_mcm=N)
        counts_e, total_e = run_job(
            backend, qc, SHOTS, f"N={N}", use_iqm)
        F_eedt = p_zero(counts_e, total_e, bit_index=N)  # target = cbit N
        gps    = F_eedt - F_ref
        se     = np.sqrt(F_eedt*(1-F_eedt)/SHOTS + F_ref*(1-F_ref)/SHOTS)
        f_th   = N * T_MCM_NS * 1e-9 / max(tau_us * 1e-6, 1e-9)

        data.append({
            "N": N, "tau_us": tau_us,
            "F_eedt": float(F_eedt), "F_ref": float(F_ref),
            "GPS": float(gps), "SE": float(se), "f_th": float(f_th)
        })
        sign = "✓" if gps > 0 else "✗"
        print(f"  N={N}: GPS={gps:+.4f}±{se:.4f}  "
              f"f_th={f_th:.1%}  {sign}")

    _save("step4_n_scan", data)
    _analyze_n_scan(data)
    return data


def _analyze_n_scan(data: list):
    gps = {d["N"]: d["GPS"] for d in data}
    n_star = max(gps, key=gps.get) if gps else None
    reversal = any(gps.get(n, 0) < 0 for n in [4, 5])

    print()
    print(f"  N*（GPS最大）= {n_star}  GPS_max={gps.get(n_star,0):+.4f}")
    print()
    if reversal and gps.get(1, 0) > 0:
        f_th_n5 = data[-1]["f_th"] if data else 0
        print(f"  ★ N-scan判定: 符号反転確認 → ✓✓ backactionモデル成立")
        print(f"     f_th(N=5)={f_th_n5:.1%}  (理論9.3%と比較)")
    elif gps.get(1, 0) > 0:
        print(f"  ★ N-scan判定: N=1で正GPS → △ 部分検証")
    else:
        print("  ★ N-scan判定: 正のGPSなし → ✗")


# ═══════════════════════════════════════════════════════════
#  プロット生成
# ═══════════════════════════════════════════════════════════

def make_plots(step3_data: list, step4_data: list, nu_zz: float):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"EEDT 最小実験セット — 検証結果  (νZZ={nu_zz/1e3:.2f} kHz)",
        fontsize=13)

    # ─── τ-scan ──────────────────────────────────────────
    ax = axes[0]
    if step3_data:
        taus = [d["tau_us"]  for d in step3_data]
        gps  = [d["GPS"]     for d in step3_data]
        se   = [d["SE"]      for d in step3_data]
        phi  = [d["phi"]     for d in step3_data]

        phi_th = np.linspace(0, max(phi) * 1.3, 300)
        g_th   = phi_th * np.exp(-phi_th) - np.exp(-1)
        scale  = max(abs(g) for g in gps) / max(abs(g_th)) if max(abs(g_th)) > 0 else 1
        g_th  *= scale

        ax.plot(phi_th, g_th, "k--", lw=1.5, alpha=0.6,
                label=r"$\propto \phi e^{-\phi}$ (理論)")
        ax.errorbar(phi, gps, yerr=se, fmt="o", color="steelblue",
                    ms=8, capsize=4, label="測定 GPS")
        ax.axhline(0, color="gray", lw=0.8)
        ax.axvline(1.0, color="red", lw=1, ls=":", alpha=0.5,
                   label="φ*=1 (最適)")
        ax.set_xlabel("φ = ωZZ × τ [rad]", fontsize=11)
        ax.set_ylabel("GPS", fontsize=11)
        ax.set_title("STEP 3: τ-scan (N=1)", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    # ─── N-scan ──────────────────────────────────────────
    ax = axes[1]
    if step4_data:
        ns   = [d["N"]   for d in step4_data]
        gps4 = [d["GPS"] for d in step4_data]
        se4  = [d["SE"]  for d in step4_data]
        tau_fix = step4_data[0]["tau_us"]
        n_star_th = tau_fix * 1e-6 / (0.093 * T_MCM_NS * 1e-9)

        ax.errorbar(ns, gps4, yerr=se4, fmt="s", color="coral",
                    ms=8, capsize=4, label="測定 GPS")
        ax.axhline(0, color="gray", lw=0.8)
        ax.axvline(n_star_th, color="red", lw=1, ls=":",
                   label=f"N*(f_th=9.3%) ≈ {n_star_th:.1f}")
        ax.set_xlabel("N (MCM回数)", fontsize=11)
        ax.set_ylabel("GPS", fontsize=11)
        ax.set_title(f"STEP 4: N-scan (τ={tau_fix:.0f}µs)", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)
        ax.set_xticks(ns)

    plt.tight_layout()
    out = OUTPUT_DIR / "eedt_result.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\n[プロット] {out}")


# ═══════════════════════════════════════════════════════════
#  最終レポート
# ═══════════════════════════════════════════════════════════

def final_report(nu_zz: float, T2: float,
                 step3_data: list, step4_data: list,
                 qa: int, qt: int, device_name: str):
    gps3 = {d["tau_us"]: d["GPS"] for d in step3_data}
    gps4 = {d["N"]: d["GPS"] for d in step4_data}
    peak3 = max((gps3.get(t, -999) for t in [30, 40, 50]), default=-999)

    c1 = "✓" if nu_zz  > 500 else ("△" if nu_zz > 100 else "✗")
    c2 = "✓" if T2     > 200e-6 else ("△" if T2 > 150e-6 else "✗")
    c3 = "✓" if peak3  > 0.005 else "✗"
    c4 = "✓" if any(gps4.get(n, 0) < 0 for n in [4, 5]) else (
         "△" if gps4.get(1, 0) > 0 else "✗")

    tau_star = T2 / (2*np.pi*max(nu_zz, 1)*T2+1) * 1e6 if T2 > 0 else 0

    print("\n" + "═" * 55)
    print("  EEDT 最小実験セット — 最終判定")
    print("═" * 55)
    print(f"  デバイス : {device_name}")
    print(f"  Qubit   : Q{qa}(ancilla) - Q{qt}(target)")
    print(f"  日時    : {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print()
    print(f"  条件① νZZ > 0.5kHz: {nu_zz/1e3:.3f} kHz  {c1}")
    print(f"  条件② T₂ > 200µs:  {T2*1e6:.1f} µs  {c2}")
    print(f"  条件③ GPS山型:      peak={peak3:.4f}  {c3}")
    print(f"  条件④ N符号反転:    "
          f"{'あり' if c4=='✓' else 'なし'}  {c4}")
    print()

    all_pass = all(x == "✓" for x in [c1, c2, c3])
    if all_pass:
        print("  ★★★ EEDT動作確定 → Quantum-6502完成度 95-100%")
    elif c1 == "✗":
        print("  ✗   νZZ≈0 → IQMのtunable couplerがZZを抑制")
        print("      → Table 2に「IQM FAIL（設計通り）」として記録")
        print("      → EEDTはZZ=0では動かないという理論が強化される")
    elif c1 == "✓" and c2 == "✓":
        print("  ★★  EEDT動作可能性あり → 完成度70-80%")
    else:
        print("  △   条件不足 → 別qubitペアで再試行を推奨")

    print()
    print("  [Table 2 追加行（論文用）]")
    verdict = "PASS" if all_pass else "FAIL"
    print(f"  {device_name} Q{qa}-Q{qt} | "
          f"νZZ={nu_zz/1e3:.2f}kHz | "
          f"T₂={T2*1e6:.0f}µs | "
          f"τ*≈{tau_star:.0f}µs | "
          f"判定={verdict}")
    print("═" * 55)

    _save("summary", {
        "device": device_name,
        "qa": qa, "qt": qt,
        "nu_zz_hz": float(nu_zz),
        "T2_us": float(T2 * 1e6),
        "tau_star_us": float(tau_star),
        "GPS_peak": float(peak3),
        "verdict": verdict,
        "total_shots": 13 * SHOTS,
        "timestamp": datetime.now().isoformat()
    })


# ═══════════════════════════════════════════════════════════
#  ユーティリティ
# ═══════════════════════════════════════════════════════════

def _save(name: str, data):
    OUTPUT_DIR.mkdir(exist_ok=True)
    path = OUTPUT_DIR / f"{name}.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


# ═══════════════════════════════════════════════════════════
#  メイン
# ═══════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="EEDT 最小実験セット（IQM/IBM統合版 v2）",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--backend", choices=["iqm", "ibm"],
                   required=True, help="バックエンド: iqm or ibm")
    p.add_argument("--token",   required=True,
                   help="API token")
    p.add_argument("--device",  default="garnet",
                   help="IQM: garnet / garnet:mock / emerald  "
                        "IBM: ibm_marrakesh 等")
    p.add_argument("--instance", default="ibm-q/open/main",
                   help="IBM instanceのみ使用")
    p.add_argument("--qa",  type=int, default=0,
                   help="ancilla qubit index（自動選択あり）")
    p.add_argument("--qt",  type=int, default=1,
                   help="target qubit index（自動選択あり）")
    p.add_argument("--skip",    default="",
                   help="スキップするSTEP番号: 例 '3,4'")
    return p.parse_args()


def main():
    args = parse_args()
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("═" * 55)
    print("  EEDT 最小実験セット  v2 (iqm-client[qiskit]対応)")
    print(f"  {13*SHOTS:,} shots  出力: {OUTPUT_DIR}")
    print("═" * 55)

    skip    = {int(s) for s in args.skip.split(",") if s.strip().isdigit()}
    use_iqm = (args.backend == "iqm")

    # ─── バックエンド接続 ─────────────────────────────
    if use_iqm:
        backend, cmap = connect_iqm(args.token, args.device)
        qa, qt = auto_select_qubits(cmap, args.qa, args.qt)
        device_name = f"IQM {args.device}"
    else:
        backend = connect_ibm(args.token, args.device, args.instance)
        qa, qt = args.qa, args.qt
        cmap = []
        device_name = args.device

    print(f"\n使用qubit: Q{qa}(ancilla) - Q{qt}(target)")
    print()

    # 結果保持
    nu_zz      = NU_ZZ_REF
    T2         = 0.0
    step3_data = []
    step4_data = []

    # ─── STEP 1: ZZ Ramsey ───────────────────────────
    if 1 not in skip:
        nu_zz = step1_zz_ramsey(backend, qa, qt, use_iqm)
    else:
        print("[STEP 1 スキップ]")

    # νZZ≈0 なら STEP3/4は意味なし（STEP2は継続）
    eedt_possible = (nu_zz > 100)

    # ─── STEP 2: T₂ ──────────────────────────────────
    if 2 not in skip:
        T2 = step2_t2(backend, qt, use_iqm)
    else:
        print("[STEP 2 スキップ]")

    if not eedt_possible:
        print("\n[INFO] νZZ≈0 のため STEP3/4 をスキップ")
        print("       EEDT は ZZ=0 では動作しない（理論通り）")
    else:
        # ─── STEP 3: τ-scan ──────────────────────────
        if 3 not in skip:
            step3_data = step3_tau_scan(
                backend, qa, qt, nu_zz, use_iqm)
        else:
            print("[STEP 3 スキップ]")

        # GPS>0 の場合のみ STEP4
        peak3 = max((d["GPS"] for d in step3_data), default=0)
        if 4 not in skip and peak3 > 0:
            step4_data = step4_n_scan(
                backend, qa, qt, nu_zz, T2, use_iqm)
        elif 4 not in skip:
            print("[STEP 4 スキップ] STEP3でGPS≤0のため")

    # ─── プロット ─────────────────────────────────────
    if step3_data or step4_data:
        make_plots(step3_data, step4_data, nu_zz)

    # ─── 最終レポート ──────────────────────────────────
    final_report(nu_zz, T2, step3_data, step4_data,
                 qa, qt, device_name)
    print(f"\n[完了] 結果: {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
