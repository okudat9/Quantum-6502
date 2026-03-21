"""
eedt_bios_v122.py
=================
Quantum-6502 BIOS v1.2.2 — Universal Pair Discovery Edition

人類初の量子OS。ZZカップリングをノイズではなくリソースとして扱う
自律的自己診断・自己適応型量子BIOSシステム。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
アーキテクチャ
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Layer 0  Universal Pair Discovery
           backend.properties()でキャリブデータ取得（ジョブ0）
           上位5ペア（搾取）+ 中央値±2ペア（探索）= 最大10候補
           過去CSV学習: 成功ペアのスコアを1.5倍に優遇

  ULC      Ultra-Light Calibration（9回路・2ジョブ）
           1qジョブ: RO×2 + T2×1
           2qジョブ: Ramsey×6（3τ × anc0/anc1）
           ν_ZZ推定: 2パラメータfit（A, ν_ZZ）T2s/delta/offset固定

  Layer 2  動作窓計算（τ* / ω_ZZ·T2 / N* / φ_ff）
           F_base推定にZZ振動項を含む（物理的に正確）

  Layer 3  GO / DEGRADED / NO-GO 判定
           NO-GO → 待たずに次候補ペアへ即移行

  EEDT     QFEED+回路（Baseline + Static Rz + 動的MCM）
           GPS_dynamic = GPS_raw - GPS_static が主指標

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
修正履歴 v1.2.1 → v1.2.2（合議録 第二次精査）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [③] ★★★ dict|None → Optional[dict]（Python 3.9互換性）
           from typing import Optional を追加

  [④] ★★★ EEDT回路: ancilla測定前にbarrier()を追加（退行バグ修正）
           delay完了後にmeasureが実行されることをコンパイラに保証
           IBM動的回路ガイドライン準拠（合議録#49 外部レビュー③）

  [⑤] ★   se()をモジュールレベルのgps_se()として移動
           テスト・再利用可能に。SHOTSを引数で受け取る形に変更

  [⑥] ★★  _get_gate_error: gate.qubitsの型安全化
           Qubit オブジェクト / int 両方に対応（Qiskitバージョン差を吸収）
           _to_int()ヘルパーで段階的フォールバック変換

修正履歴 v1.2.0 → v1.2.1（合議録 第一次精査）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  [A] ★★★ FFT削除・Ramsey fit を2パラメータ(A, ν_ZZ)に限定
  [B] ★★★ f_base_est にZZ振動項 cos(ω_ZZ·τ) を追加
  [C] ★★★ boot_decision内の重複ロジックを削除
  [D] ★★  get_p0: total=0のゼロ除算ガードを追加
  [E] ★★  score_pair: gate属性名をhasattr+fallbackで安全化
  [F] ★★★ ULC: 1q回路と2q回路を別jobに分割
  [G] ★★  Ramsey点を [10, 30, 90] に変更

著者: Takeshi Okuda（Independent, Osaka）
実装: Claude（Anthropic）
合議体: Dirac（理論）/ Woz（アーキテクチャ）/ Feynman（実機）
DOI: 10.5281/zenodo.19029431
GitHub: github.com/okudat9/Quantum-6502
更新: 2026-03-21 v1.2.2
"""

import math
import os
import csv
import datetime
import sys
import atexit
import glob
from typing import Optional
import numpy as np
from scipy.optimize import curve_fit
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, Session, SamplerV2 as Sampler

# ============================================================
# 設定（この4行を変えるだけで任意バックエンド・任意ペアに対応）
# ============================================================
TOKEN = os.environ.get("IBM_QUANTUM_TOKEN", "YOUR_IBM_TOKEN_HERE")
BACKEND_NAME = "ibm_marrakesh"

# ============================================================
# BIOS定数
# ============================================================

# ショット数
SHOTS_ULC_1Q  = 512    # Readout cal + T2 probe（1q、軽量）
SHOTS_ULC_2Q  = 1024   # Ramsey scan（2q）
SHOTS_EEDT    = 8192   # EEDT本実験

# ULC Ramseyスキャン点（修正G: 10/30/90µs）
# 10µs : ν_ZZ=50kHzで半周期 → 高周波検出
# 30µs : ν_ZZ=3kHz付近のτ*に近い → 中周波検出
# 90µs : T2減衰確認 + ν_ZZ=1kHz付近の1周期
RAMSEY_ULC_US = [10, 30, 90]

# T2簡易チェック用τ（T2の~1/2を狙う）
T2_PROBE_US = 150.0

# EEDTの理論基準値（3/15実測・論文確定値）
NU_ZZ_TARGET   = 2.778e3   # Hz feedforward式定義（= ν_ZZ_phys / 2）
NU_ZZ_PHYS_REF = 5.556e3   # Hz C(τ)振動周波数
PHI_STAR       = 0.873     # rad 最適動作点（実測・理論値と5%以内）
T2_REF_US      = 287.7     # µs フォールバック用参照値

# GO / DEGRADED / NO-GO しきい値
F_BASE_GO        = 0.75    # GOに必要な最低F_base推定値
F_BASE_DEGRADED  = 0.60    # DEGRADEDに必要な最低F_base推定値
NU_ZZ_MIN_KHZ    = 1.5     # ω_ZZ·T2動作窓: ν_ZZ下限
NU_ZZ_MAX_KHZ    = 6.0     # ω_ZZ·T2動作窓: ν_ZZ上限
OMEGA_T2_MIN     = 1.0     # 動作窓下限
OMEGA_T2_MAX     = 10.0    # 動作窓上限
FIT_QUALITY_WARN = 0.20    # ν_ZZ相対誤差 > 20% → POOR（Cramér-Rao下限から導出）
N_STAR_MAX       = 5       # Separation Theorem上限（τ=50µs基準）

# QFEED+補正可能判定ウィンドウ（修正済み: mod-2πベース）
# cos曲線の同じ山に乗っているか = φ_actual % 2π が PHI_STAR ± π/4 以内
QFEED_CORRECTION_WINDOW = math.pi / 4  # ±45度

# MCMオーバーヘッド（IBM Heron実測値 2025年）
MCM_FF_LATENCY_US = 0.600   # フィードフォワード古典処理遅延
MCM_READOUT_US    = 1.340   # MidCircuit測定時間
MCM_OVERHEAD_US   = MCM_FF_LATENCY_US + MCM_READOUT_US   # ≈ 1.94µs

# Layer 0 ペア選定パラメータ
TOP_N          = 5     # 上位ペア数（搾取）
MEDIAN_N       = 5     # 中央値周辺ペア数（探索）
HISTORY_BOOST  = 1.5   # 過去成功ペアのスコア倍率
CSV_HISTORY_GLOB = "eedt_bios_*.csv"   # 過去CSV検索パターン

BIOS_VERSION = "1.2.2"


# ============================================================
# ユーティリティ関数（モジュールレベル: from eedt_bios import * 可能）
# ============================================================

def us_to_dt(t_us: float, dt: float) -> int:
    """µs → dt単位（四捨五入）。round()でN回ループの系統誤差を防ぐ。"""
    return round((t_us * 1e-6) / dt)


def get_p0(pub, reg: str = 'c0'):
    """
    1bitレジスタからP(|0>)を安全に取得。

    修正D: total=0のゼロ除算ガードを追加。
    Returns: (p0: float, counts: dict)
    """
    if not hasattr(pub.data, reg):
        available = [a for a in dir(pub.data) if not a.startswith('_')]
        raise ValueError(
            f"レジスタ '{reg}' が存在しません。利用可能: {available}")
    counts = getattr(pub.data, reg).get_counts()
    total  = sum(counts.values())
    if total == 0:
        # ショット数0 → 測定不能、安全値0.5を返す
        return 0.5, counts
    c0 = counts.get('0', 0)
    return c0 / total, counts


def rem_correct(f_raw: float, M_inv) -> float:
    """
    REM（Readout Error Mitigation）補正。
    M_inv=Noneのとき（特異行列）は生値をそのまま返す。
    """
    if M_inv is None:
        return f_raw
    vec = np.array([f_raw, 1.0 - f_raw])
    corrected = np.clip(M_inv @ vec, 0.0, 1.0)
    s = corrected.sum()
    if s == 0:
        return f_raw
    corrected /= s
    return float(corrected[0])


def can_qfeed_correct(nu_zz_code: float, seg_us: float) -> bool:
    """
    QFEED+補正可能判定（修正C: 唯一の実装。boot_decision内から呼び出す）

    mod-2πベース判定:
      φ_actual = 2π × ν_ZZ_code × τ_seg
      補正可能 ⟺ |φ_actual % 2π - PHI_STAR| < π/4 （円周距離）

    背景: fez Q88の実験（ν_ZZ=25kHz, τ=48µs）でφ_actual=7.71rad。
          7.71 % 2π = 1.43rad → 目標0.873radとの差0.56rad > π/4
          → 正しく「補正不可」と判定できる。
    """
    phi_actual = 2.0 * math.pi * nu_zz_code * seg_us * 1e-6
    phi_mod    = phi_actual % (2.0 * math.pi)
    delta      = abs(phi_mod - PHI_STAR)
    delta      = min(delta, 2.0 * math.pi - delta)   # 円周上の最短距離
    return delta < QFEED_CORRECTION_WINDOW


def ramsey_model_2param(tau, A, nu_zz, T2s_fixed, delta_fixed, offset_fixed):
    """
    2パラメータRamseyモデル（修正A）。
    T2s / delta / offset を外側から固定して呼び出すことで
    3データ点に対して数学的に解ける（自由度 = 3 - 2 = 1）。

    使用例:
        def model(tau, A, nu_zz):
            return ramsey_model_2param(tau, A, nu_zz,
                                       T2s_est, 0.0, 0.0)
        popt, pcov = curve_fit(model, tau_arr, C_tau, ...)
    """
    return (A
            * np.cos(2.0 * math.pi * nu_zz * tau + delta_fixed)
            * np.exp(-tau / T2s_fixed)
            + offset_fixed)


def estimate_f_base(tau_us: float, T2_est_us: float,
                    nu_zz_code: float) -> float:
    """
    F_baseline推定（修正B: ZZ振動項を含む物理的に正確な式）。

    F_base = 0.5 × (1 + exp(-τ/T2) × cos(ω_ZZ × τ))
    ω_ZZ = 2π × ν_ZZ_code

    背景: 旧実装は cos 項を欠いており F_base を常に過大評価していた。
          ν_ZZ=3kHz, τ=50µs → cos(0.94rad)≈0.59 で補正が入る。
    """
    tau    = tau_us * 1e-6
    t2     = max(T2_est_us * 1e-6, 1e-9)   # ゼロ除算防止
    omega  = 2.0 * math.pi * nu_zz_code
    decay  = math.exp(-tau / t2)
    cosine = math.cos(omega * tau)
    return max(0.0, min(1.0, 0.5 * (1.0 + decay * cosine)))


def gps_se(fa: float, fb: float, shots: int) -> float:
    """
    GPS z-score用の標準誤差（誤差伝播）。
    SE = sqrt(fa*(1-fa)/n + fb*(1-fb)/n)
    モジュールレベルに置くことでテスト・再利用を可能にする。
    """
    return math.sqrt(fa * (1.0 - fa) / shots + fb * (1.0 - fb) / shots)


# ============================================================
# 回路ファクトリ関数（モジュールレベル）
# ============================================================

def make_ro_cal(state: int) -> QuantumCircuit:
    """Readoutキャリブ回路: state=0→|0>測定, state=1→|1>測定"""
    qr = QuantumRegister(1, 'q')
    cr = ClassicalRegister(1, 'c')
    qc = QuantumCircuit(qr, cr)
    if state == 1:
        qc.x(qr[0])
    qc.measure(qr[0], cr[0])
    qc.name = f"RO_cal_{state}"
    return qc


def make_t2_probe(tau_us: float, dt: float) -> QuantumCircuit:
    """T2簡易チェック回路（1qubit Ramsey、デチューニングなし）"""
    qr = QuantumRegister(1, 'q')
    cr = ClassicalRegister(1, 'c')
    qc = QuantumCircuit(qr, cr)
    qc.h(qr[0])
    qc.delay(us_to_dt(tau_us, dt), qr[0], unit='dt')
    qc.h(qr[0])
    qc.measure(qr[0], cr[0])
    qc.name = f"T2probe_{tau_us:.0f}us"
    return qc


def make_ramsey_anc0(tau_us: float, dt: float) -> QuantumCircuit:
    """
    Ramsey回路: ancilla |0⟩
    ancilla|0⟩ → target に -ν_ZZ_phys/2 のZZシフトが入る
    ancillaにも同じdelayを入れてZZ結合を対称に働かせる
    """
    qr = QuantumRegister(2, 'q')
    cr = ClassicalRegister(1, 'c0')
    qc = QuantumCircuit(qr, cr)
    qc.h(qr[0])
    qc.delay(us_to_dt(tau_us, dt), qr[0], unit='dt')
    qc.delay(us_to_dt(tau_us, dt), qr[1], unit='dt')
    qc.h(qr[0])
    qc.measure(qr[0], cr[0])
    qc.name = f"R0_{tau_us:.0f}us"
    return qc


def make_ramsey_anc1(tau_us: float, dt: float) -> QuantumCircuit:
    """
    Ramsey回路: ancilla |1⟩
    ancilla|1⟩ → target に +ν_ZZ_phys/2 のZZシフトが入る
    C(τ) = P0(anc=0) - P0(anc=1) の振動周波数 = ν_ZZ_phys
    """
    qr = QuantumRegister(2, 'q')
    cr = ClassicalRegister(1, 'c0')
    qc = QuantumCircuit(qr, cr)
    qc.x(qr[1])
    qc.h(qr[0])
    qc.delay(us_to_dt(tau_us, dt), qr[0], unit='dt')
    qc.delay(us_to_dt(tau_us, dt), qr[1], unit='dt')
    qc.h(qr[0])
    qc.measure(qr[0], cr[0])
    qc.name = f"R1_{tau_us:.0f}us"
    return qc


# ============================================================
# Layer 0: Universal Pair Discovery（ジョブ0: IBMアクセス不要）
# ============================================================

def load_csv_history(pattern: str = CSV_HISTORY_GLOB) -> set:
    """
    過去CSVから成功ペアを読み込む。
    成功条件: GPS_dynamic > 0 かつ z_score >= 1.0
    _nogo.csvにはGPS_dynamic列がないためサイレントスキップされる（正常動作）。
    Returns: set of (target_q, ancilla_q) as int tuples
    """
    successful_pairs = set()
    for path in glob.glob(pattern):
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if not row.get('target_q'):
                        continue
                    try:
                        gps_dyn = float(row.get('GPS_dynamic', 0) or 0)
                        z_score = float(row.get('z_score', 0) or 0)
                        tq = int(float(row['target_q']))
                        aq = int(float(row['ancilla_q']))
                        if gps_dyn > 0 and z_score >= 1.0:
                            successful_pairs.add((tq, aq))
                    except (ValueError, KeyError):
                        continue
        except Exception:
            continue
    return successful_pairs


def _get_gate_error(props, q1: int, q2: int) -> float:
    """
    2qubitゲートエラーを安全に取得（修正E + 修正⑥）。

    Qiskitバージョン差を吸収:
      - 旧 (0.x): gate.qubits が Qubit オブジェクトのリスト
      - 新 (1.x): gate.qubits が int のリスト
      どちらの形式でも int インデックスに変換して比較する。

    取得できない場合はペナルティ値0.01を返す（スコアリングは継続）。
    """
    target_set = frozenset([q1, q2])

    def _to_int(q) -> Optional[int]:
        """Qubit オブジェクトまたは int を安全に int に変換する。"""
        if isinstance(q, int):
            return q
        # Qubit オブジェクト: _index 属性を試みる
        for attr in ('_index', 'index'):
            val = getattr(q, attr, None)
            if val is not None:
                try:
                    return int(val)
                except (TypeError, ValueError):
                    pass
        # 最終手段: 文字列表現から数値を抽出 (例: "Qubit(QuantumRegister(2,'q'), 1)")
        try:
            s = str(q)
            # 末尾の ), N) パターンから数値を取り出す
            num = int(s.rstrip(')').split(',')[-1].strip())
            return num
        except (ValueError, IndexError):
            return None

    try:
        for gate in props.gates:
            # gate.gate / gate.name どちらの属性名にも対応（修正E）
            gate_name = (getattr(gate, 'gate', None)
                         or getattr(gate, 'name', ''))
            if gate_name not in ('cx', 'cz', 'ecr'):
                continue

            # qubitインデックスを int に変換して比較（修正⑥）
            indices = set()
            for q in gate.qubits:
                idx = _to_int(q)
                if idx is None:
                    indices = None   # 変換失敗 → このゲートをスキップ
                    break
                indices.add(idx)

            if indices is None or frozenset(indices) != target_set:
                continue

            for param in gate.parameters:
                if param.name == 'gate_error':
                    return float(param.value)

    except Exception:
        pass

    return 0.01   # 取得不可時のペナルティ値


def score_pair(props, q1: int, q2: int) -> float:
    """
    キャリブデータからペアスコアを計算（修正E）。
    スコア = T2 / 500µs × (1 - CZ×20) × (1 - RO×10)

    T2を主軸にCZエラー・ROエラーをペナルティとして乗じる。
    Returns: float ≥ 0（高いほど良いペア）
    """
    try:
        t2_q1 = props.qubit_property(q1, 'T2').value * 1e6   # µs
        t2_q2 = props.qubit_property(q2, 'T2').value * 1e6
        t2    = (t2_q1 + t2_q2) / 2.0

        ro_q1 = props.qubit_property(q1, 'readout_error').value
        ro_q2 = props.qubit_property(q2, 'readout_error').value
        ro    = (ro_q1 + ro_q2) / 2.0

        cz = _get_gate_error(props, q1, q2)

        score = (t2 / 500.0)
        score *= max(0.0, 1.0 - cz * 20.0)
        score *= max(0.0, 1.0 - ro * 10.0)
        return max(0.0, score)

    except Exception:
        return 0.0


def discover_pairs(backend) -> list:
    """
    Layer 0: キャリブデータからEEDT候補ペアを自動選定。

    アルゴリズム:
      1. backend.properties()で全qubitのT2/RO/CZを取得（ジョブ0）
      2. 全coupling_mapエッジをスコアリング（重複除去）
      3. 上位TOP_N（搾取）+ 中央値周辺MEDIAN_N（探索）を統合
      4. 過去CSV成功ペアのスコアをHISTORY_BOOST倍に優遇

    Returns: list of (target_q, ancilla_q, score) sorted by priority
    """
    print("\n┌─ Layer 0: Universal Pair Discovery ──────────────────")

    props    = backend.properties()
    cal_time = props.last_update_date
    try:
        age_h = ((datetime.datetime.now(cal_time.tzinfo) - cal_time)
                 .total_seconds() / 3600.0)
        print(f"│  キャリブデータ: {cal_time.strftime('%Y-%m-%d %H:%M')}"
              f" ({age_h:.1f}時間前)")
        if age_h > 12:
            print("│  ⚠ データが古い可能性（T2実態と乖離リスク）")
    except Exception:
        print("│  キャリブ時刻: 取得不可")

    # 過去成功ペアを読み込む
    history = load_csv_history()
    if history:
        print(f"│  過去成功ペア: {len(history)}件 → スコア×{HISTORY_BOOST}")

    cmap = backend.coupling_map

    # 全エッジをスコアリング（coupling_map非対称対応）
    seen   = set()
    scored = []
    for q1, q2 in cmap.get_edges():
        key = tuple(sorted([q1, q2]))
        if key in seen:
            continue
        seen.add(key)

        s = score_pair(props, key[0], key[1])
        if s <= 0.0:
            continue

        # 過去成功ペアのスコア優遇（順序不問）
        if (key in history or (key[1], key[0]) in history):
            s *= HISTORY_BOOST

        scored.append((key[0], key[1], s))

    if not scored:
        print("│  [警告] スコアリング可能なペアが見つかりません")
        print("└" + "─" * 54)
        return []

    scored.sort(key=lambda x: x[2], reverse=True)
    total = len(scored)
    print(f"│  結合ペア総数: {total}  スコアリング完了")

    # 上位TOP_N（搾取）
    top_pairs  = scored[:TOP_N]

    # 中央値周辺MEDIAN_N（探索）
    mid        = total // 2
    lo         = max(0, mid - 2)
    hi         = min(total, mid + 3)
    mid_pairs  = scored[lo:hi]

    # 重複除去して統合（上位優先）
    seen_keys  = set()
    candidates = []
    for entry in top_pairs + mid_pairs:
        key = (entry[0], entry[1])
        if key not in seen_keys:
            seen_keys.add(key)
            candidates.append(entry)

    # 表示
    print(f"│  上位{TOP_N}ペア   : "
          f"{[(c[0], c[1]) for c in top_pairs[:TOP_N]]}")
    print(f"│  中央値周辺: "
          f"{[(c[0], c[1]) for c in mid_pairs]}")
    print(f"│  候補合計  : {len(candidates)}ペア")
    print("└" + "─" * 54)

    return candidates


# ============================================================
# ULC: Ultra-Light Calibration（修正A, D, F, G）
# ============================================================

def run_ulc(backend, session, target_q: int, ancilla_q: int,
            dt: float) -> Optional[dict]:
    """
    Ultra-Light Calibration（9回路・2ジョブ）。

    修正F: 1q回路と2q回路を別jobに分割。
           異なるlayoutを1ジョブに混在させるとIBMが拒否するため。
    修正G: Ramseyを[10, 30, 90]µsに変更。
           旧46µsはν_ZZ=50kHzで7周目 → フィット不能。
    修正A: Ramseyフィットを2パラメータ(A, ν_ZZ)に限定。
           T2s=T2_est固定・delta=0固定・offset=0固定
           3点データで自由度1 → 数学的に解ける。

    Returns: dict（診断結果）または None（結合なし/エラー）
    """
    print(f"\n┌─ ULC: Q{target_q}+Q{ancilla_q}"
          + "─" * max(0, 46 - len(str(target_q)) - len(str(ancilla_q))))

    # coupling_mapバリデーション
    cmap = backend.coupling_map
    connected = (cmap.graph.has_edge(target_q, ancilla_q)
                 or cmap.graph.has_edge(ancilla_q, target_q))
    if not connected:
        print(f"│  ✗ Q{target_q}-Q{ancilla_q} は直接結合なし → スキップ")
        print("└" + "─" * 54)
        return None

    # パスマネージャ（1q / 2qは分離）
    pm_1q = generate_preset_pass_manager(
        backend=backend, optimization_level=1,
        initial_layout=[target_q])
    pm_2q = generate_preset_pass_manager(
        backend=backend, optimization_level=1,
        initial_layout=[target_q, ancilla_q])

    # ─── Job 1: 1q回路（RO×2 + T2×1）─── 修正F
    circs_1q = [
        pm_1q.run(make_ro_cal(0)),
        pm_1q.run(make_ro_cal(1)),
        pm_1q.run(make_t2_probe(T2_PROBE_US, dt)),
    ]

    # ─── Job 2: 2q回路（Ramsey 3τ × anc0/anc1 = 6回路）─── 修正F
    circs_2q = []
    for tau in RAMSEY_ULC_US:
        circs_2q.append(pm_2q.run(make_ramsey_anc0(tau, dt)))
        circs_2q.append(pm_2q.run(make_ramsey_anc1(tau, dt)))

    # 2ジョブを並列submit
    samp_1q = Sampler(mode=session)
    samp_1q.options.default_shots = SHOTS_ULC_1Q
    samp_2q = Sampler(mode=session)
    samp_2q.options.default_shots = SHOTS_ULC_2Q

    job_1q = samp_1q.run(circs_1q)
    job_2q = samp_2q.run(circs_2q)
    print(f"│  1q Job: {job_1q.job_id()}")
    print(f"│  2q Job: {job_2q.job_id()}")
    print("│  結果待機中（2ジョブ並列）...")

    try:
        result_1q = job_1q.result(timeout=600)
        result_2q = job_2q.result(timeout=600)
    except Exception as e:
        print(f"│  ✗ ジョブ失敗: {e}")
        print("└" + "─" * 54)
        return None

    # ── Layer 1-A: Readout解析 ──────────────────────────────
    p00, _ = get_p0(result_1q[0], 'c')
    p01, _ = get_p0(result_1q[1], 'c')
    p10, p11 = 1.0 - p00, 1.0 - p01
    M = np.array([[p00, p01], [p10, p11]])
    try:
        M_inv  = np.linalg.inv(M)
        ro_ok  = True
    except np.linalg.LinAlgError:
        M_inv  = None
        ro_ok  = False
        print("│  [警告] 混同行列が特異 → REM補正を無効化")

    ro_quality = p00 - p01
    ro_status  = ("GOOD"     if ro_quality > 0.8
                  else ("DEGRADED" if ro_quality > 0.5
                  else  "POOR"))
    print(f"│  RO品質: {ro_quality:.3f} ({ro_status})")

    # ── Layer 1-B: T2簡易推定 ──────────────────────────────
    f_probe, _ = get_p0(result_1q[2], 'c')
    val = min(2.0 * f_probe - 1.0, 0.9999)
    if 0.1 < val < 0.999:
        T2_est_us = -T2_PROBE_US / math.log(val)
        t2_valid  = True
    else:
        # コヒーレンス喪失またはval範囲外 → フォールバック
        T2_est_us = T2_REF_US * 0.5
        t2_valid  = False

    print(f"│  T2推定: {T2_est_us:.1f} µs"
          + (" (フォールバック)" if not t2_valid else ""))

    # ── Layer 1-C: Ramsey解析（修正A: 2パラメータfit）──────
    tau_arr = np.array(RAMSEY_ULC_US) * 1e-6
    C_list  = []
    for i in range(len(RAMSEY_ULC_US)):
        f0, _ = get_p0(result_2q[2 * i],     'c0')
        f1, _ = get_p0(result_2q[2 * i + 1], 'c0')
        C_list.append(f0 - f1)
    C_tau = np.array(C_list)

    fit_ok           = False
    fit_quality      = "UNKNOWN"
    fit_quality_warn = False
    nu_zz_phys       = NU_ZZ_PHYS_REF
    nu_zz_code       = NU_ZZ_TARGET
    nu_zz_err        = float('nan')
    T2s_ramsey       = T2_est_us * 1e-6

    try:
        T2s_fixed = max(T2_est_us * 1e-6, 5e-6)

        # 修正A: 2パラメータfit（A, ν_ZZ）
        # T2s=固定・delta=0・offset=0 → 3点で自由度1 → 解ける
        def model_2p(tau, A, nu_zz):
            return ramsey_model_2param(tau, A, nu_zz,
                                       T2s_fixed, 0.0, 0.0)

        # 初期値: 振幅は実測C_tauのmax、ν_ZZは複数候補で試す
        A_init    = max(abs(C_tau.max()), abs(C_tau.min()), 0.05)
        best_popt = None
        best_pcov = None
        best_res  = float('inf')

        # ν_ZZ初期候補: 1〜50kHzを対数スケールで7点試す（局所最小回避）
        nu_candidates = np.logspace(
            math.log10(1e3), math.log10(50e3), 7)

        for nu_init in nu_candidates:
            try:
                popt, pcov = curve_fit(
                    model_2p, tau_arr, C_tau,
                    p0=[A_init, nu_init],
                    bounds=([0.0, 0.5e3], [1.0, 80e3]),
                    maxfev=5000)
                residual = np.sum((model_2p(tau_arr, *popt) - C_tau) ** 2)
                if residual < best_res:
                    best_res  = residual
                    best_popt = popt
                    best_pcov = pcov
            except RuntimeError:
                continue

        if best_popt is None:
            raise RuntimeError("全初期値でフィット失敗")

        nu_zz_phys = best_popt[1]
        nu_zz_code = nu_zz_phys / 2.0
        nu_zz_err  = (math.sqrt(best_pcov[1, 1])
                      if (best_pcov is not None
                          and np.isfinite(best_pcov[1, 1]))
                      else float('inf'))
        fit_ok = True

        # フィット品質評価
        rel_err = nu_zz_err / nu_zz_code if nu_zz_code > 0 else float('inf')
        if rel_err > FIT_QUALITY_WARN:
            fit_quality      = "POOR"
            fit_quality_warn = True
        else:
            fit_quality = "GOOD"

        print(f"│  ν_ZZ_phys: {nu_zz_phys/1e3:.3f} ± "
              f"{nu_zz_err/1e3:.3f} kHz  [{fit_quality}]")
        print(f"│  ν_ZZ_code: {nu_zz_code/1e3:.3f} kHz")

        # QFEED+補正可能性チェック（情報表示のみ）
        tau_star_est = PHI_STAR / (2.0 * math.pi * nu_zz_code) * 1e6
        seg_est      = tau_star_est   # N=1想定
        qfeed_check  = can_qfeed_correct(nu_zz_code, seg_est)
        print(f"│  QFEED+補正: {'✓ 可能' if qfeed_check else '✗ 不可（窓外）'}"
              f"  τ*≈{tau_star_est:.1f}µs")

    except RuntimeError as e:
        print(f"│  Ramseyフィット失敗: {e}")
        print(f"│  ⚠ フィット失敗 → QFEED+補正の信頼性が担保できません")
        print(f"│  → fit_quality=FAILED として boot_decision で NO-GO へ移行します")
        print(f"│  （未知の物理的相互作用・異常ノイズの兆候として記録します）")
        fit_quality      = "FAILED"
        fit_quality_warn = True
        # fit_quality_warn=True により boot_decision は GO に昇格しない。
        # フォールバック値（NU_ZZ_TARGET）で位相補正を試みると
        # 実態と乖離した位相を強要しコヒーレンスを破壊するリスクがある。
        # → このペアは DEGRADED / NO-GO として次候補へ移行する。
        # 合議体 Dirac: 「ノイズが大きすぎると切り捨てるのではなく、
        # 既存理論と衝突した兆候として診断CSVに記録する」

    print("└" + "─" * 54)

    return {
        "target_q":        target_q,
        "ancilla_q":       ancilla_q,
        "nu_zz_phys":      nu_zz_phys,
        "nu_zz_code":      nu_zz_code,
        "nu_zz_err":       nu_zz_err,
        "T2_est_us":       T2_est_us,
        "t2_valid":        t2_valid,
        "T2s_ramsey":      T2s_ramsey,
        "ro_quality":      ro_quality,
        "ro_status":       ro_status,
        "ro_ok":           ro_ok,
        "M_inv":           M_inv,
        "fit_ok":          fit_ok,
        "fit_quality":     fit_quality,
        "fit_quality_warn": fit_quality_warn,
        "C_tau":           C_list,
        "pm_2q":           pm_2q,
        "job_1q_id":       job_1q.job_id(),
        "job_2q_id":       job_2q.job_id(),
    }


# ============================================================
# Layer 2: 動作窓計算
# ============================================================

def compute_window(ulc: dict) -> dict:
    """
    ULC結果から動作窓パラメータを計算。

    τ*  = PHI_STAR / (2π × ν_ZZ_code)   ← EEDT最適貯蔵時間
    ω_ZZ·T2 ≈ 1〜10                       ← EEDT動作窓定義
    N*  ≤ N_STAR_MAX = 5                  ← Separation Theorem
    """
    nu_zz_code  = ulc["nu_zz_code"]
    T2_est_us   = ulc["T2_est_us"]

    tau_star_us  = PHI_STAR / (2.0 * math.pi * nu_zz_code) * 1e6
    omega_zz_T2  = 2.0 * math.pi * nu_zz_code * T2_est_us * 1e-6
    in_window    = OMEGA_T2_MIN <= omega_zz_T2 <= OMEGA_T2_MAX
    nu_ok        = NU_ZZ_MIN_KHZ <= nu_zz_code / 1e3 <= NU_ZZ_MAX_KHZ

    # N*: T2 / (2 × τ*) で推定（上限はSeparation Theorem）
    N_star = min(max(1, int(T2_est_us / (2.0 * tau_star_us))), N_STAR_MAX)

    return {
        "tau_star_us": tau_star_us,
        "omega_zz_T2": omega_zz_T2,
        "in_window":   in_window,
        "nu_ok":       nu_ok,
        "N_star":      N_star,
    }


# ============================================================
# Layer 3: GO / DEGRADED / NO-GO 判定（修正B, C）
# ============================================================

def boot_decision(ulc: dict, win: dict) -> dict:
    """
    GO / DEGRADED / NO-GO を判定してEEDT実行パラメータを確定。

    修正B: estimate_f_base()でZZ振動項を考慮（物理的に正確）
    修正C: can_qfeed_correct()を直接呼び出し（重複ロジック削除）

    GO条件:
      F_base_est ≥ 0.75  かつ  nu_ok  かつ  fit_ok  かつ  fit_quality GOOD
    DEGRADED条件（N=1縮退）:
      F_base_est ≥ 0.60  かつ  nu_ok
    NO-GO:
      それ以外 → 次候補ペアへ即移行（待機なし）
    """
    nu_zz_code       = ulc["nu_zz_code"]
    T2_est_us        = ulc["T2_est_us"]
    nu_ok            = win["nu_ok"]
    fit_ok           = ulc["fit_ok"]
    fit_quality_warn = ulc["fit_quality_warn"]
    tau_star_us      = win["tau_star_us"]

    # F_base推定（修正B: ZZ振動項を含む）
    f_base_est = estimate_f_base(tau_star_us, T2_est_us, nu_zz_code)

    # 判定
    if (f_base_est >= F_BASE_GO
            and nu_ok
            and fit_ok
            and not fit_quality_warn):
        boot_mode = "GO"
        N_run     = min(win["N_star"], 3)
        tau_run   = tau_star_us
    elif f_base_est >= F_BASE_DEGRADED and nu_ok:
        boot_mode = "DEGRADED"
        N_run     = 1
        tau_run   = tau_star_us
    else:
        boot_mode = "NO-GO"
        N_run     = 1
        tau_run   = tau_star_us

    # φ_ff計算（セグメント単位）
    seg_us        = tau_run / N_run
    phi_ff_target = 2.0 * math.pi * NU_ZZ_TARGET * seg_us * 1e-6
    delta_phi     = (2.0 * math.pi
                     * (nu_zz_code - NU_ZZ_TARGET)
                     * seg_us * 1e-6)

    # 修正C: can_qfeed_correct()を唯一の実装として呼び出す
    qfeed_ok     = can_qfeed_correct(nu_zz_code, seg_us)
    phi_ff_total = phi_ff_target + delta_phi if qfeed_ok else phi_ff_target

    return {
        "boot_mode":     boot_mode,
        "N_run":         N_run,
        "tau_run":       tau_run,
        "seg_us":        seg_us,
        "phi_ff_target": phi_ff_target,
        "delta_phi":     delta_phi,
        "phi_ff_total":  phi_ff_total,
        "qfeed_ok":      qfeed_ok,
        "f_base_est":    f_base_est,
    }


# ============================================================
# EEDT実験（Baseline + Static Rz + QFEED+）
# ============================================================

def run_eedt(backend, session, ulc: dict, win: dict,
             dec: dict, dt: float) -> dict:
    """
    EEDT本実験。3条件を同一ジョブで測定:
      1. Baseline    : MCMなし・補正なし
      2. Static Rz   : MCMなし・静的位相補正のみ（対照実験）
      3. EEDT QFEED+ : MCMあり・動的feedforward補正

    GPS_dynamic = GPS_raw - GPS_static
    = 動的MCMの純粋な寄与（readout error / 静的位相補正の影響を排除）
    """
    target_q     = ulc["target_q"]
    ancilla_q    = ulc["ancilla_q"]
    M_inv        = ulc["M_inv"]
    ro_ok        = ulc["ro_ok"]
    pm_2q        = ulc["pm_2q"]
    tau_run      = dec["tau_run"]
    N_run        = dec["N_run"]
    phi_ff_total = dec["phi_ff_total"]
    qfeed_ok     = dec["qfeed_ok"]

    # MCMオーバーヘッド補正（delay = τ_run/N - MCM_overhead）
    seg_us_raw    = tau_run / N_run
    seg_corrected = max(seg_us_raw - MCM_OVERHEAD_US, 0.5)
    tau_eff_us    = seg_corrected * N_run + MCM_OVERHEAD_US * N_run
    SEG_DT        = us_to_dt(seg_corrected, dt)

    print(f"\n[EEDT] {dec['boot_mode']} | Q{target_q}+Q{ancilla_q}")
    print(f"  τ={tau_run:.1f}µs  τ_eff≈{tau_eff_us:.2f}µs  N={N_run}")
    print(f"  φ_ff={phi_ff_total:.4f}rad  "
          f"QFEED+={'✓' if qfeed_ok else '✗（通常QFEED使用）'}")

    # ── Baseline回路 ──────────────────────────────────────────
    qr_b = QuantumRegister(2, 'q')
    cr_b = ClassicalRegister(1, 'c0')
    qc_b = QuantumCircuit(qr_b, cr_b)
    qc_b.h(qr_b[0])
    qc_b.delay(us_to_dt(tau_run, dt), qr_b[0], unit='dt')
    qc_b.h(qr_b[0])
    qc_b.measure(qr_b[0], cr_b[0])
    qc_b.name = "Baseline"

    # ── Static Rz回路（対照実験）────────────────────────────
    # 同じ位相補正をMCMなしで静的に適用
    # GPS_dynamic = GPS_raw - GPS_static で「動的MCMの純粋寄与」を分離
    qr_s = QuantumRegister(2, 'q')
    cr_s = ClassicalRegister(1, 'c0')
    qc_s = QuantumCircuit(qr_s, cr_s)
    qc_s.h(qr_s[0])
    qc_s.delay(us_to_dt(tau_run, dt), qr_s[0], unit='dt')
    qc_s.rz(-phi_ff_total * N_run, qr_s[0])   # N_run分の合計位相を一括補正
    qc_s.h(qr_s[0])
    qc_s.measure(qr_s[0], cr_s[0])
    qc_s.name = "Static_Rz"

    # ── EEDT QFEED+回路 ──────────────────────────────────────
    # $C0 QINIT: |+⟩準備
    # $C2 QMCM × N_run: ancilla測定（MidCircuit）
    # $C3' QFEED+: ancilla=1のときのみ Rz(-φ_ff_total) 適用
    # $C4 QSYNC: 最終測定
    qr_e   = QuantumRegister(2, 'q')
    cr_anc = ClassicalRegister(N_run, 'anc')
    cr_tgt = ClassicalRegister(1, 'tgt')
    qc_e   = QuantumCircuit(qr_e, cr_anc, cr_tgt)
    qc_e.h(qr_e[0])
    for i in range(N_run):
        qc_e.delay(SEG_DT, qr_e[0], unit='dt')
        # barrier: delayが完了してからmeasureが実行されることをコンパイラに保証
        # IBM動的回路ガイドライン準拠（合議録#49 外部レビュー③）
        qc_e.barrier(qr_e[0], qr_e[1])
        qc_e.measure(qr_e[1], cr_anc[i])
        clbit = cr_anc[i]   # if_testにはClbitを明示
        with qc_e.if_test((clbit, 1)):
            qc_e.rz(-phi_ff_total, qr_e[0])
    qc_e.h(qr_e[0])
    qc_e.measure(qr_e[0], cr_tgt[0])
    qc_e.name = "EEDT_QFEED+"

    # トランスパイル & 実行
    isa_b = pm_2q.run(qc_b)
    isa_s = pm_2q.run(qc_s)
    isa_e = pm_2q.run(qc_e)

    sampler = Sampler(mode=session)
    sampler.options.default_shots = SHOTS_EEDT
    job_e = sampler.run([isa_b, isa_s, isa_e])
    print(f"  EEDT Job: {job_e.job_id()}")

    result_e = job_e.result(timeout=900)

    # 測定結果取得
    F_base_raw,   cnt_b = get_p0(result_e[0], 'c0')
    F_static_raw, cnt_s = get_p0(result_e[1], 'c0')
    F_eedt_raw,   cnt_e = get_p0(result_e[2], 'tgt')

    # REM補正（参考値）
    F_base_rem   = rem_correct(F_base_raw,   M_inv) if ro_ok else F_base_raw
    F_static_rem = rem_correct(F_static_raw, M_inv) if ro_ok else F_static_raw
    F_eedt_rem   = rem_correct(F_eedt_raw,   M_inv) if ro_ok else F_eedt_raw

    # GPS計算
    GPS_raw     = F_eedt_raw   - F_base_raw
    GPS_static  = F_static_raw - F_base_raw
    GPS_dynamic = GPS_raw - GPS_static   # 動的MCMの純粋寄与（主指標）
    GPS_rem     = F_eedt_rem   - F_base_rem

    # z-score（誤差伝播: モジュールレベルのgps_se()を使用）
    SE       = gps_se(F_base_raw, F_eedt_raw,   SHOTS_EEDT)
    SE_s     = gps_se(F_base_raw, F_static_raw, SHOTS_EEDT)
    z        = GPS_raw    / SE   if SE   > 0 else 0.0
    z_static = GPS_static / SE_s if SE_s > 0 else 0.0

    return {
        "job_id":        job_e.job_id(),
        "F_base_raw":    F_base_raw,
        "F_static_raw":  F_static_raw,
        "F_eedt_raw":    F_eedt_raw,
        "F_base_rem":    F_base_rem,
        "F_static_rem":  F_static_rem,
        "F_eedt_rem":    F_eedt_rem,
        "GPS_raw":       GPS_raw,
        "GPS_static":    GPS_static,
        "GPS_dynamic":   GPS_dynamic,
        "GPS_rem":       GPS_rem,
        "SE":            SE,
        "z":             z,
        "z_static":      z_static,
        "tau_eff_us":    tau_eff_us,
        "depth_base":    isa_b.depth(),
        "depth_static":  isa_s.depth(),
        "depth_eedt":    isa_e.depth(),
    }


# ============================================================
# 結果表示
# ============================================================

def print_result(ulc: dict, win: dict, dec: dict, eedt: dict,
                 attempt: int, n_total: int) -> tuple:
    """結果を整形して表示。(verdict, dyn_verdict)を返す。"""
    z           = eedt["z"]
    GPS_dynamic = eedt["GPS_dynamic"]

    # EEDT全体の判定
    if z >= 3.0:
        verdict = "✅ 3σ以上 → 統計的に有意"
    elif z >= 2.0:
        verdict = "✓✓ 2σ以上 → 有望"
    elif z >= 1.0:
        verdict = "✓  1σ以上 → 正の傾向"
    else:
        verdict = "△  1σ未満 → 条件継続確認"

    # 動的MCM寄与の判定
    if GPS_dynamic > 0 and z >= 2.0:
        dyn_verdict = "★ 動的MCMが静的補正を上回る（本物の動的効果）"
    elif GPS_dynamic > 0:
        dyn_verdict = "△ 動的有利だが有意差未確定"
    elif GPS_dynamic < 0:
        dyn_verdict = "✗ 静的補正の方が高い → MCMオーバーヘッドが利得を消している"
    else:
        dyn_verdict = "― 差なし"

    print()
    print("=" * 62)
    print(f"  Quantum-6502 BIOS v{BIOS_VERSION} — 実験結果")
    print("=" * 62)
    print(f"  日時     : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  試行     : {attempt}/{n_total}ペア")
    print(f"  対象     : Q{ulc['target_q']}+Q{ulc['ancilla_q']}")
    print(f"  Boot Mode: {dec['boot_mode']}")
    print()
    print("  ── 診断 ──────────────────────────────────────────")
    print(f"  ν_ZZ_code  : {ulc['nu_zz_code']/1e3:.4f} kHz"
          f"  [{ulc['fit_quality']}]")
    print(f"  T2推定     : {ulc['T2_est_us']:.1f} µs")
    print(f"  RO品質     : {ulc['ro_quality']:.3f}"
          f"  ({ulc['ro_status']})")
    print(f"  ω_ZZ×T2   : {win['omega_zz_T2']:.2f}"
          f"  {'[IN]' if win['in_window'] else '[OUT]'}")
    print(f"  F_base推定 : {dec['f_base_est']:.4f}")
    print(f"  φ_ff_total : {dec['phi_ff_total']:.4f} rad"
          f"  QFEED+={'✓' if dec['qfeed_ok'] else '✗'}")
    print()
    print("  ── 3条件比較（REM補正なし・主指標）─────────────")
    print(f"  F_baseline : {eedt['F_base_raw']:.4f}")
    print(f"  F_static   : {eedt['F_static_raw']:.4f}"
          f"  GPS_static ={eedt['GPS_static']:+.4f}"
          f"  z={eedt['z_static']:.2f}σ")
    print(f"  F_EEDT     : {eedt['F_eedt_raw']:.4f}"
          f"  GPS_eedt   ={eedt['GPS_raw']:+.4f}"
          f"  z={eedt['z']:.2f}σ")
    print(f"  GPS_dynamic: {eedt['GPS_dynamic']:+.4f}"
          f"  ← 動的MCMの純粋寄与（主指標）")
    print()
    print("  ── REM補正あり（参考）────────────────────────────")
    print(f"  F_base(REM): {eedt['F_base_rem']:.4f}")
    print(f"  F_EEDT(REM): {eedt['F_eedt_rem']:.4f}")
    print(f"  GPS(REM)   : {eedt['GPS_rem']:+.4f}")
    print()
    print(f"  判定（EEDT）    : {verdict}")
    print(f"  判定（動的寄与）: {dyn_verdict}")
    print()
    print(f"  参照: 3/15 marrakesh  GPS=+0.726  z=66.0σ")
    print("=" * 62)

    return verdict, dyn_verdict


# ============================================================
# CSV保存
# ============================================================

def save_csv_full(ulc: dict, win: dict, dec: dict, eedt: dict,
                  verdict: str, dyn_verdict: str,
                  attempt: int) -> str:
    """EEDT成功時の完全CSVを保存。"""
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fn = f"eedt_bios_{ts}.csv"
    with open(fn, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow([
            "bios_version", "attempt", "timestamp", "backend",
            "target_q", "ancilla_q", "boot_mode",
            "ro_quality", "ro_status", "fit_quality",
            "nu_zz_code_kHz", "T2_est_us", "tau_star_us",
            "tau_run_us", "tau_eff_us", "N_run",
            "omega_zz_T2", "N_star", "qfeed_ok",
            "phi_ff_total", "F_base_est",
            "F_baseline_raw", "F_static_raw", "F_eedt_raw",
            "GPS_static", "GPS_raw", "GPS_dynamic",
            "F_baseline_rem", "F_static_rem", "F_eedt_rem", "GPS_rem",
            "SE", "z_static", "z_score",
            "depth_base", "depth_static", "depth_eedt",
            "verdict", "dynamic_verdict",
            "ulc_job_1q", "ulc_job_2q", "eedt_job",
        ])
        w.writerow([
            BIOS_VERSION, attempt,
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            BACKEND_NAME, ulc["target_q"], ulc["ancilla_q"],
            dec["boot_mode"],
            round(ulc["ro_quality"], 4), ulc["ro_status"], ulc["fit_quality"],
            round(ulc["nu_zz_code"] / 1e3, 6), round(ulc["T2_est_us"], 1),
            round(win["tau_star_us"], 2),
            round(dec["tau_run"], 2), round(eedt["tau_eff_us"], 2),
            dec["N_run"],
            round(win["omega_zz_T2"], 3), win["N_star"], dec["qfeed_ok"],
            round(dec["phi_ff_total"], 6), round(dec["f_base_est"], 4),
            round(eedt["F_base_raw"], 6), round(eedt["F_static_raw"], 6),
            round(eedt["F_eedt_raw"], 6),
            round(eedt["GPS_static"], 6), round(eedt["GPS_raw"], 6),
            round(eedt["GPS_dynamic"], 6),
            round(eedt["F_base_rem"], 6), round(eedt["F_static_rem"], 6),
            round(eedt["F_eedt_rem"], 6), round(eedt["GPS_rem"], 6),
            round(eedt["SE"], 6), round(eedt["z_static"], 4),
            round(eedt["z"], 4),
            eedt["depth_base"], eedt["depth_static"], eedt["depth_eedt"],
            verdict, dyn_verdict,
            ulc["job_1q_id"], ulc["job_2q_id"], eedt["job_id"],
        ])
        # ULC Ramseyデータ
        w.writerow([])
        w.writerow(["--- ULC Ramsey C(tau) ---"])
        w.writerow(["tau_us", "C_tau"])
        for tau, c in zip(RAMSEY_ULC_US, ulc["C_tau"]):
            w.writerow([tau, round(c, 6)])
    print(f"\n[BIOS] CSV保存: {fn}")
    return fn


def save_csv_nogo(ulc: dict, win: dict, attempt: int) -> str:
    """NO-GO時の簡易CSVを保存。load_csv_historyはこれをスキップする（正常）。"""
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fn = f"eedt_bios_{ts}_nogo.csv"
    with open(fn, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow([
            "bios_version", "attempt", "timestamp", "backend",
            "target_q", "ancilla_q", "boot_mode",
            "nu_zz_code_kHz", "T2_est_us", "omega_zz_T2",
            "ro_quality", "fit_quality",
            "ulc_job_1q", "ulc_job_2q",
        ])
        w.writerow([
            BIOS_VERSION, attempt,
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            BACKEND_NAME, ulc["target_q"], ulc["ancilla_q"], "NO-GO",
            round(ulc["nu_zz_code"] / 1e3, 6),
            round(ulc["T2_est_us"], 1),
            round(win["omega_zz_T2"], 3),
            round(ulc["ro_quality"], 4), ulc["fit_quality"],
            ulc["job_1q_id"], ulc["job_2q_id"],
        ])
    return fn


# ============================================================
# メインループ
# ============================================================

if __name__ == "__main__":

    ts_boot = datetime.datetime.now()
    print("=" * 62)
    print("  ██████╗ ██╗ ██████╗ ███████╗")
    print("  ██╔══██╗██║██╔═══██╗██╔════╝")
    print("  ██████╔╝██║██║   ██║███████╗")
    print("  ██╔══██╗██║██║   ██║╚════██║")
    print("  ██████╔╝██║╚██████╔╝███████║")
    print("  ╚═════╝ ╚═╝ ╚═════╝ ╚══════╝")
    print(f"  Quantum-6502 BIOS v{BIOS_VERSION}")
    print(f"  Boot: {ts_boot.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Backend: {BACKEND_NAME}")
    print("=" * 62)

    # ── IBM接続 ──────────────────────────────────────────────
    service = QiskitRuntimeService(
        channel="ibm_quantum_platform", token=TOKEN)
    backend = service.backend(BACKEND_NAME)
    DT = backend.dt
    print(f"  接続OK: {backend.name}  dt={DT*1e9:.4f}ns\n")

    # ── Layer 0: ペア候補選定（ジョブ0）──────────────────────
    candidates = discover_pairs(backend)
    if not candidates:
        print("[BIOS] 候補ペアなし → 終了")
        sys.exit(0)

    # ── Session開始 ──────────────────────────────────────────
    session = Session(backend=backend)
    atexit.register(session.close)
    atexit.register(print, "\n[BIOS] Session終了")
    print(f"  Session: {session.session_id}\n")

    # ── 候補ペアを順番に試す（待たずに次へ）─────────────────
    found = False
    for attempt, (tq, aq, score) in enumerate(candidates, 1):
        print(f"\n{'━' * 62}")
        print(f"  候補 {attempt}/{len(candidates)}: "
              f"Q{tq}+Q{aq}  スコア={score:.3f}")
        print(f"{'━' * 62}")

        # ULC（9回路・2ジョブ）
        ulc = run_ulc(backend, session, tq, aq, DT)
        if ulc is None:
            continue   # 結合なし → 即次候補

        # Layer 2: 動作窓計算
        win = compute_window(ulc)
        print(f"  τ*={win['tau_star_us']:.1f}µs  "
              f"ω·T2={win['omega_zz_T2']:.2f}  "
              f"N*={win['N_star']}  "
              f"窓内={'✓' if win['in_window'] else '✗'}")

        # Layer 3: 判定
        dec = boot_decision(ulc, win)
        print(f"  → {dec['boot_mode']}"
              f"  F_base推定={dec['f_base_est']:.3f}")

        if dec["boot_mode"] == "NO-GO":
            fn = save_csv_nogo(ulc, win, attempt)
            print(f"  NO-GO CSV: {fn}")
            print("  次候補へ即移行（待機なし）")
            continue

        # GO / DEGRADED → EEDT実行
        try:
            eedt = run_eedt(backend, session, ulc, win, dec, DT)
        except Exception as e:
            print(f"  [ERROR] EEDT実行失敗: {e}")
            print("  次候補へ即移行")
            continue

        verdict, dyn_verdict = print_result(
            ulc, win, dec, eedt, attempt, len(candidates))
        save_csv_full(ulc, win, dec, eedt, verdict, dyn_verdict, attempt)
        found = True
        break   # 1ペア成功したら終了（全ペア実行したい場合はbreakを削除）

    if not found:
        print(f"\n[BIOS] 全{len(candidates)}候補がNO-GO → 終了")
        print("  次回: 新アカウント取得またはキャリブ回復後に再実行")
        sys.exit(0)
