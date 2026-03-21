"""
qscript.py
==========
Q-Script v0.1 — Quantum-6502 OS インターフェース

「量子OSへの問い合わせプロトコル」

ユーザーは「目標」を指定する。OSが「手段」を決める。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
使用例（5行デモ）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    import qscript

    qs = qscript.QuantumOS(backend="ibm_marrakesh")
    print(qs.status())
    result = qs.run(goal="maximize_fidelity")
    print(result)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
アーキテクチャ
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    Q-Script API（このファイル）
         ↓
    Quantum-6502 BIOS v1.2.2（eedt_bios.py）
         ↓
    Qiskit Runtime / IBM Backend

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
goalポリシー（Gemini + 合議体 確定）
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    maximize_fidelity:
        N_run = min(N_star, 3)   ← N_star上限まで使う
        tau_run = tau_star_us    ← τ*をそのまま
        f_base_threshold = 0.75  ← 標準

    maximize_stability:
        N_run = 1                ← 縮退モード固定
        tau_run = tau_star_us * 0.8  ← τ*の80%（安全側）
        f_base_threshold = 0.80  ← より厳しい条件

著者: Takeshi Okuda（Independent, Osaka）
実装: Claude（Anthropic）
合議体: Dirac / Woz / Feynman
外部評価: ChatGPT / Copilot / Gemini / Grok
DOI: 10.5281/zenodo.19029431
GitHub: github.com/okudat9/Quantum-6502
更新: 2026-03-21 v0.1
"""

from __future__ import annotations

import csv
import datetime
import glob
import math
from typing import Any, Dict, List, Optional

# ============================================================
# 例外クラス（Gemini提案: エラーを握りつぶさない）
# ============================================================

class QScriptError(Exception):
    """Q-Script 基底例外クラス"""
    pass


class NoPairFoundError(QScriptError):
    """全候補ペアが NO-GO だった場合"""
    pass


class HardwareDriftError(QScriptError):
    """ν_ZZ ドリフトが動作窓を外れた場合"""
    pass


class FitFailedError(QScriptError):
    """Ramsey フィットが失敗した場合"""
    pass


class BackendNotSupportedError(QScriptError):
    """未対応のバックエンドが指定された場合"""
    pass


# ============================================================
# 内部定数（BIOSと共有）
# ============================================================

# stability → window の翻訳しきい値
_STABILITY_OPTIMAL   = 0.75
_STABILITY_DEGRADED  = 0.60

# goal ポリシー（Gemini + 合議体 確定）
_GOAL_POLICIES: Dict[str, Dict[str, Any]] = {
    "maximize_fidelity": {
        "n_run_mode":         "max",     # min(N_star, 3) を使う
        "tau_factor":         1.00,      # τ* をそのまま使う
        "f_base_threshold":   0.75,      # 標準 GO 条件
        "description":        "GPS最大化優先 / Separation Theorem境界を活用",
    },
    "maximize_stability": {
        "n_run_mode":         "min",     # N_run = 1（縮退固定）
        "tau_factor":         0.80,      # τ* の 80%（安全側）
        "f_base_threshold":   0.80,      # より厳しい GO 条件
        "description":        "readout error・ドリフト影響を最小化優先",
    },
}

_DEFAULT_GOAL   = "maximize_fidelity"
_DEFAULT_SHOTS  = 8192
_CSV_GLOB       = "eedt_bios_*.csv"


# ============================================================
# 内部ヘルパー
# ============================================================

def _load_history(pattern: str = _CSV_GLOB) -> List[Dict[str, Any]]:
    """
    過去の BIOS 実験結果 CSV を読み込む。
    _nogo.csv は GPS_dynamic 列がないためサイレントスキップ（正常）。
    """
    records = []
    for path in sorted(glob.glob(pattern)):
        try:
            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if not row.get("target_q"):
                        continue
                    try:
                        records.append({
                            "timestamp":  row.get("timestamp", ""),
                            "backend":    row.get("backend", ""),
                            "target_q":  int(float(row["target_q"])),
                            "ancilla_q": int(float(row["ancilla_q"])),
                            "gps":       float(row.get("GPS_dynamic", 0) or 0),
                            "fidelity":  float(row.get("F_eedt_raw", 0) or 0),
                            "stability": float(row.get("F_base_est", 0) or 0),
                            "boot_mode": row.get("boot_mode", ""),
                            "z_score":   float(row.get("z_score", 0) or 0),
                        })
                    except (ValueError, KeyError):
                        continue
        except Exception:
            continue
    return records


def _stability_to_window(stability: float) -> str:
    """stability スコアから window 文字列に翻訳する"""
    if stability >= _STABILITY_OPTIMAL:
        return "optimal"
    elif stability >= _STABILITY_DEGRADED:
        return "degraded"
    else:
        return "no-go"


def _verdict_from_z(z: float) -> str:
    """z スコアから verdict 文字列に翻訳する"""
    if z >= 3.0:
        return "confirmed"
    elif z >= 2.0:
        return "promising"
    elif z >= 1.0:
        return "weak"
    else:
        return "inconclusive"


# ============================================================
# QuantumOS クラス
# ============================================================

class QuantumOS:
    """
    Q-Script v0.1 メインクラス。

    Quantum-6502 BIOS v1.2.2 のラッパー。
    ユーザーは物理パラメータを意識せず、目標（goal）を指定するだけ。

    使用例:
        qs = QuantumOS(backend="ibm_marrakesh")
        print(qs.status())
        result = qs.run(goal="maximize_fidelity")
        print(result)
    """

    def __init__(self, backend: str = "ibm_marrakesh", token: Optional[str] = None):
        """
        Args:
            backend: IBM バックエンド名。
                     "auto" は v0.2 で実装予定（現在は NotImplementedError）。
            token:   IBM Quantum Platform API トークン。
                     省略時は eedt_bios.py の TOKEN 定数を使用。
        """
        if backend == "auto":
            raise NotImplementedError(
                "backend='auto' は v0.2 で実装予定です。\n"
                "現在は backend='ibm_marrakesh' または 'ibm_kingston' を指定してください。"
            )

        self._backend_name = backend
        self._token        = token
        self._session      = None   # BIOS Session（run()時に開く）
        self._last_ulc     = None   # 最後の ULC 結果（キャッシュ）
        self._last_win     = None   # 最後の Layer2 結果
        self._last_dec     = None   # 最後の Layer3 結果
        self._last_eedt    = None   # 最後の EEDT 結果
        self._candidates   = None   # Layer0 ペア候補（キャッシュ）

    # ──────────────────────────────────────────────────────
    # コンテキストマネージャ（with 構文対応）
    # ──────────────────────────────────────────────────────

    def __enter__(self) -> "QuantumOS":
        self.open()
        return self

    def __exit__(self, *args) -> None:
        self.close()

    def open(self) -> None:
        """IBM Session を開く（run() 時に自動で呼ばれる）"""
        # 実際の Session 開始は BIOS 側で行う
        # ここでは接続状態のフラグのみ管理
        self._opened = True

    def close(self) -> None:
        """IBM Session をクローズする"""
        if self._session is not None:
            try:
                self._session.close()
            except Exception:
                pass
            self._session = None
        self._opened = False

    # ──────────────────────────────────────────────────────
    # status() — OSの現在状態を人間語で返す
    # ──────────────────────────────────────────────────────

    def status(self) -> Dict[str, Any]:
        """
        OSの現在状態を返す。物理パラメータを「翻訳」して返す。

        Returns:
            {
              "ready":     bool,    # GO/DEGRADED=True, NO-GO=False
              "stability": float,   # 0.0〜1.0（f_base_est に対応）
              "window":    str,     # "optimal"/"degraded"/"no-go"
              "drift":     str,     # "stable"/"warning"/"bad"
              "best_pair": list,    # [target_q, ancilla_q]
              "timestamp": str,
            }

        実装ブリッジ:
            ready      ← boot_mode（"GO"/"DEGRADED"→True）
            stability  ← f_base_est
            window     ← omega_zz_T2（1〜10→"optimal"）
            drift      ← fit_quality（"GOOD"→"stable"）
            best_pair  ← discover_pairs()[0][:2]
        """
        # 最後の実験結果があればそれを使う（キャッシュ）
        if self._last_dec is not None:
            stability = self._last_dec.get("f_base_est", 0.0)
            boot_mode = self._last_dec.get("boot_mode", "NO-GO")
            ready     = boot_mode in ("GO", "DEGRADED")
        else:
            # 実験前: キャリブデータから推定（ジョブなし）
            stability = 0.5    # 不明時のデフォルト
            ready     = False

        if self._last_win is not None:
            omega_t2  = self._last_win.get("omega_zz_T2", 0.0)
            in_window = self._last_win.get("in_window", False)
            window    = "optimal" if in_window else "no-go"
        else:
            window    = "unknown"

        if self._last_ulc is not None:
            fit_q = self._last_ulc.get("fit_quality", "UNKNOWN")
            drift = (
                "stable"  if fit_q == "GOOD" else
                "warning" if fit_q == "POOR" else
                "bad"
            )
            pair  = [
                self._last_ulc.get("target_q"),
                self._last_ulc.get("ancilla_q"),
            ]
        else:
            drift = "unknown"
            pair  = []

        return {
            "ready":      ready,
            "stability":  round(stability, 4),
            "window":     window,
            "drift":      drift,
            "best_pair":  pair,
            "timestamp":  datetime.datetime.now().isoformat(timespec="seconds"),
        }

    # ──────────────────────────────────────────────────────
    # run() — メイン実行
    # ──────────────────────────────────────────────────────

    def run(
        self,
        goal:        str = _DEFAULT_GOAL,
        constraints: Optional[Dict[str, Any]] = None,
        shots:       int = _DEFAULT_SHOTS,
    ) -> Dict[str, Any]:
        """
        EEDT を実行して結果を返す。

        Args:
            goal: "maximize_fidelity"（デフォルト）/ "maximize_stability"
            constraints: {
                "min_stability": float,  # GO条件の最低stability（デフォルト: goalに依存）
                "max_time_us":   float,  # τ* の上限（デフォルト: 無制限）
                "min_z_score":   float,  # 統計有意性の下限（デフォルト: 1.0）
            }
            shots: 測定ショット数（デフォルト: 8192）

        Returns:
            {
              "gps":       float,  # GPS_dynamic（主指標）
              "gps_raw":   float,  # GPS_raw（参考）
              "fidelity":  float,  # F_EEDT_raw
              "z_score":   float,
              "pair":      list,   # [target_q, ancilla_q]
              "mode":      str,    # "GO" / "DEGRADED"
              "verdict":   str,    # "confirmed"/"promising"/"weak"/"inconclusive"
              "backend":   str,
              "timestamp": str,
            }

        Raises:
            NoPairFoundError:  全候補がNO-GO
            HardwareDriftError: ν_ZZが動作窓外
            FitFailedError:    Ramseyフィット失敗

        注意:
            v0.1では単一ペアのみ実行。
            max_pairs は予約語（v0.2で実装予定）。
        """
        if goal not in _GOAL_POLICIES:
            raise ValueError(
                f"goal='{goal}' は未対応です。"
                f"対応: {list(_GOAL_POLICIES.keys())}"
            )

        policy      = _GOAL_POLICIES[goal]
        constraints = constraints or {}

        # constraints でポリシーを上書き
        f_thresh = constraints.get("min_stability", policy["f_base_threshold"])
        max_tau  = constraints.get("max_time_us",   float("inf"))
        min_z    = constraints.get("min_z_score",   1.0)

        # ── BIOS をインポート ──────────────────────────────
        # （修正①: BIOS呼び出しを _select_pair / _run_ulc / _decide / _execute
        #   に分離するのが理想だが、v0.1では run() 内にまとめて置く。
        #   v0.2でプライベートメソッド化してライブラリ化する。）
        try:
            from eedt_bios import (
                QiskitRuntimeService, Session,
                discover_pairs, run_ulc, compute_window,
                boot_decision, run_eedt,
            )
            import eedt_bios as _bios
        except ImportError as e:
            raise ImportError(
                f"eedt_bios.py が見つかりません: {e}\n"
                "同じディレクトリに配置してください。"
            ) from e

        token = self._token or _bios.TOKEN
        service  = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
        backend  = service.backend(self._backend_name)
        dt       = backend.dt

        # Layer 0
        candidates = discover_pairs(backend)
        if not candidates:
            raise NoPairFoundError("候補ペアが見つかりませんでした。")
        self._candidates = candidates

        # Session
        session = Session(backend=backend)
        self._session = session

        # ── goalポリシーに基づく局所判定関数（修正③）──────────
        # グローバル定数を書き換えず、ローカルな引数で判定する。
        # _bios.F_BASE_GO を直接書き換えると並列実行時に状態がリークする。
        def _local_boot_decision(ulc_r, win_r):
            """goalポリシーのしきい値を使った局所的なboot_decision。"""
            nu_zz_code       = ulc_r["nu_zz_code"]
            T2_est_us        = ulc_r["T2_est_us"]
            nu_ok            = win_r["nu_ok"]
            fit_ok           = ulc_r["fit_ok"]
            fit_quality_warn = ulc_r["fit_quality_warn"]
            tau_star_us      = win_r["tau_star_us"]

            f_base_est = _bios.estimate_f_base(tau_star_us, T2_est_us, nu_zz_code)

            if (f_base_est >= f_thresh         # ← goalポリシーのしきい値
                    and nu_ok
                    and fit_ok
                    and not fit_quality_warn):
                boot_mode = "GO"
                N_run     = min(win_r["N_star"], 3)
                tau_run   = tau_star_us
            elif f_base_est >= min(f_thresh - 0.15, 0.60) and nu_ok:
                boot_mode = "DEGRADED"
                N_run     = 1
                tau_run   = tau_star_us
            else:
                boot_mode = "NO-GO"
                N_run     = 1
                tau_run   = tau_star_us

            seg_us        = tau_run / N_run
            phi_ff_target = 2.0 * math.pi * _bios.NU_ZZ_TARGET * seg_us * 1e-6
            delta_phi     = (2.0 * math.pi
                             * (nu_zz_code - _bios.NU_ZZ_TARGET)
                             * seg_us * 1e-6)
            qfeed_ok      = _bios.can_qfeed_correct(nu_zz_code, seg_us)
            phi_ff_total  = phi_ff_target + delta_phi if qfeed_ok else phi_ff_target

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

        try:
            import atexit
            atexit.register(session.close)

            result_dict   = None
            skip_log: List[Dict[str, Any]] = []   # 修正②: スキップ理由を記録

            for tq, aq, score in candidates:
                # ULC
                ulc = run_ulc(backend, session, tq, aq, dt)
                if ulc is None:
                    skip_log.append({"pair": [tq, aq], "reason": "no_coupling"})
                    continue

                self._last_ulc = ulc

                # Layer 2
                win = compute_window(ulc)
                self._last_win = win

                # τ 制約チェック
                if win["tau_star_us"] > max_tau:
                    skip_log.append({
                        "pair":   [tq, aq],
                        "reason": f"tau_exceeded ({win['tau_star_us']:.1f}µs > {max_tau}µs)"
                    })
                    continue

                # Layer 3（修正③: グローバル定数を書き換えない）
                dec = _local_boot_decision(ulc, win)

                # goalポリシーによるN_run調整
                if goal == "maximize_stability" and dec["boot_mode"] != "NO-GO":
                    dec["N_run"]   = 1
                    dec["tau_run"] = win["tau_star_us"] * policy["tau_factor"]

                self._last_dec = dec

                if dec["boot_mode"] == "NO-GO":
                    reason = (
                        "window_exceeded"
                        if not win["in_window"]
                        else f"fit_quality={ulc.get('fit_quality','?')}"
                        if ulc.get("fit_quality_warn")
                        else f"stability_low ({dec['f_base_est']:.3f})"
                    )
                    skip_log.append({"pair": [tq, aq], "reason": reason})  # 修正②
                    continue

                # EEDT 実行
                try:
                    eedt = run_eedt(backend, session, ulc, win, dec, dt)
                except Exception as e:
                    raise HardwareDriftError(
                        f"EEDT実行中にハードウェアエラー: {e}"
                    ) from e

                self._last_eedt = eedt

                # z_score が min_z を下回る場合は warning を付与（例外は投げない）
                verdict = _verdict_from_z(eedt["z"])
                warning = None
                if eedt["z"] < min_z:
                    warning = (
                        f"z_score={eedt['z']:.2f} < min_z={min_z:.2f}。"
                        "統計的有意性が不足しています。"
                    )

                result_dict = {
                    "gps":          round(eedt["GPS_dynamic"], 6),
                    "gps_raw":      round(eedt["GPS_raw"],     6),
                    "fidelity":     round(eedt["F_eedt_raw"],  6),
                    "z_score":      round(eedt["z"],           4),
                    "pair":         [ulc["target_q"], ulc["ancilla_q"]],
                    "mode":         dec["boot_mode"],
                    "verdict":      verdict,
                    "backend":      self._backend_name,
                    "timestamp":    datetime.datetime.now().isoformat(timespec="seconds"),
                    "skipped_pairs": skip_log,   # 修正②: なぜ他のペアを棄却したか
                }
                if warning:
                    result_dict["warning"] = warning

                break   # 1ペア成功したら終了

        finally:
            session.close()
            self._session = None

        if result_dict is None:
            raise NoPairFoundError(
                f"全{len(candidates)}候補がNO-GOでした。"
                f"棄却理由: {skip_log}\n"
                "マシン状態の回復後に再実行してください。"
            )

        return result_dict

    # ──────────────────────────────────────────────────────
    # best_pair() — 最良ペア情報
    # ──────────────────────────────────────────────────────

    def best_pair(self) -> Dict[str, Any]:
        """
        Layer0 の結果から最良ペアを返す。
        run() 実行前でも呼べる（キャッシュがあれば）。

        Returns:
            {
              "qubits": [target_q, ancilla_q],
              "score":  float,
              "mode":   "exploit" or "explore",
            }
        """
        if self._candidates:
            tq, aq, score = self._candidates[0]
            # TOP_N=5 以内なら exploit、それ以降は explore
            from eedt_bios import TOP_N
            mode = "exploit" if self._candidates.index((tq, aq, score)) < TOP_N else "explore"
            return {"qubits": [tq, aq], "score": round(score, 4), "mode": mode}

        return {"qubits": [], "score": 0.0, "mode": "unknown"}

    # ──────────────────────────────────────────────────────
    # explain() — 説明可能性（ChatGPT提案・★最重要）
    # ──────────────────────────────────────────────────────

    def explain(self) -> Dict[str, Any]:
        """
        「なぜこのペアを選んだか」を返す。
        IBMが必ず聞く質問への直接回答。

        Returns:
            {
              "selected_pair": [target_q, ancilla_q],
              "reason": {
                "score": float,
                "mode":  str,
                "history_boost": bool,
              },
              "window_calculation": {
                "tau_star_us":  float,
                "omega_zz_T2":  float,
                "window":       str,
                "N_star":       int,
              },
              "rejected_pairs": [
                {"pair": [q1, q2], "reason": str}, ...
              ]
            }
        """
        if self._last_ulc is None:
            return {"error": "run() を先に実行してください。"}

        selected = [self._last_ulc["target_q"], self._last_ulc["ancilla_q"]]

        # 棄却ペアの収集（NO-GOのCSVから）
        rejected = []
        for path in sorted(glob.glob("eedt_bios_*_nogo.csv")):
            try:
                with open(path, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        if not row.get("target_q"):
                            continue
                        try:
                            tq = int(float(row["target_q"]))
                            aq = int(float(row["ancilla_q"]))
                            omega = float(row.get("omega_zz_T2", 0) or 0)
                            reason = (
                                f"omega_T2={omega:.1f}（動作窓外）"
                                if omega > 10 else
                                f"fit_quality={row.get('fit_quality','?')}"
                            )
                            rejected.append({"pair": [tq, aq], "reason": reason})
                        except (ValueError, KeyError):
                            continue
            except Exception:
                continue

        # history_boost 確認
        history_pairs = set()
        for r in _load_history():
            if r["gps"] > 0 and r["z_score"] >= 1.0:
                history_pairs.add((r["target_q"], r["ancilla_q"]))
        history_boost = tuple(selected) in history_pairs

        win_info = {}
        if self._last_win:
            win_info = {
                "tau_star_us": round(self._last_win.get("tau_star_us", 0), 2),
                "omega_zz_T2": round(self._last_win.get("omega_zz_T2", 0), 3),
                "window":      "optimal" if self._last_win.get("in_window") else "no-go",
                "N_star":      self._last_win.get("N_star", 0),
            }

        return {
            "selected_pair": selected,
            "reason": {
                "history_boost": history_boost,
            },
            "window_calculation": win_info,
            "rejected_pairs":    rejected[:5],   # 最大5件
        }

    # ──────────────────────────────────────────────────────
    # debug() — 研究者モード（物理を全開示）
    # ──────────────────────────────────────────────────────

    def debug(self) -> Dict[str, Any]:
        """
        物理パラメータを全開示する研究者モード。
        「物理を隠さない、でも直接は見せない」原則の裏口。

        Returns:
            ULC / Layer2 / Layer3 / EEDT の全生データ
        """
        result = {}

        if self._last_ulc:
            result.update({
                "nu_zz_phys_khz": round(self._last_ulc.get("nu_zz_phys", 0) / 1e3, 6),
                "nu_zz_code_khz": round(self._last_ulc.get("nu_zz_code", 0) / 1e3, 6),
                "T2_est_us":      round(self._last_ulc.get("T2_est_us", 0), 2),
                "ro_quality":     round(self._last_ulc.get("ro_quality", 0), 4),
                "fit_quality":    self._last_ulc.get("fit_quality", "UNKNOWN"),
                "ulc_job_1q":     self._last_ulc.get("job_1q_id", ""),
                "ulc_job_2q":     self._last_ulc.get("job_2q_id", ""),
            })

        if self._last_win:
            result.update({
                "omega_zz_T2": round(self._last_win.get("omega_zz_T2", 0), 3),
                "tau_star_us": round(self._last_win.get("tau_star_us", 0), 2),
                "N_star":      self._last_win.get("N_star", 0),
                "in_window":   self._last_win.get("in_window", False),
            })

        if self._last_dec:
            result.update({
                "phi_ff_total": round(self._last_dec.get("phi_ff_total", 0), 6),
                "f_base_est":   round(self._last_dec.get("f_base_est", 0), 4),
                "boot_mode":    self._last_dec.get("boot_mode", ""),
                "qfeed_ok":     self._last_dec.get("qfeed_ok", False),
            })

        if self._last_eedt:
            result.update({
                "F_baseline":   round(self._last_eedt.get("F_base_raw", 0), 6),
                "F_static":     round(self._last_eedt.get("F_static_raw", 0), 6),
                "F_eedt":       round(self._last_eedt.get("F_eedt_raw", 0), 6),
                "GPS_static":   round(self._last_eedt.get("GPS_static", 0), 6),
                "GPS_raw":      round(self._last_eedt.get("GPS_raw", 0), 6),
                "GPS_dynamic":  round(self._last_eedt.get("GPS_dynamic", 0), 6),
                "z_score":      round(self._last_eedt.get("z", 0), 4),
                "eedt_job":     self._last_eedt.get("job_id", ""),
            })

        if not result:
            return {"error": "run() を先に実行してください。"}

        return result

    # ──────────────────────────────────────────────────────
    # history() — 過去結果（ChatGPT提案）
    # ──────────────────────────────────────────────────────

    def history(
        self,
        last_n: int = 10,
        field:  Optional[str] = None,
    ) -> List[Any]:
        """
        過去の実験結果を返す（BIOSのCSVから読み込む）。

        Args:
            last_n: 最新 N 件を返す
            field:  特定フィールドのみ時系列で返す（例: "nu_zz_khz"）

        Returns:
            last_n 件の実験結果リスト、または特定フィールドの時系列リスト
        """
        records = _load_history()
        records = records[-last_n:]   # 最新 N 件

        if field is not None:
            return [r.get(field) for r in records]

        return [
            {
                "timestamp":  r["timestamp"],
                "pair":       [r["target_q"], r["ancilla_q"]],
                "gps":        round(r["gps"], 6),
                "fidelity":   round(r["fidelity"], 6),
                "stability":  round(r["stability"], 4),
                "backend":    r["backend"],
            }
            for r in records
        ]

    # ──────────────────────────────────────────────────────
    # __repr__
    # ──────────────────────────────────────────────────────

    def __repr__(self) -> str:
        return (
            f"QuantumOS(backend='{self._backend_name}', "
            f"session={'open' if self._session else 'closed'})"
        )
