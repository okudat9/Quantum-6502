"""
======================================================================
  Quantum-6502 Prototype v0.1
  EEDT-Native Hybrid Quantum Processor Emulator
  
  Council: Dirac (理論) + Woz (アーキテクチャ) + Feynman (実機)
  Author : Takeshi Okuda (okudat9) — Independent Quantum Researcher
  Base   : EEDT v19, Zenodo 10.5281/zenodo.18899504
  Env    : Windows / Python 3.10+ / Qiskit 2.x / Qiskit-Aer
======================================================================

命令セット ($C0–$C7):
  $C0  QINIT   qb, state       — qubit初期化 |0> or |+>
  $C1  QGATE   qb, gate        — H / X / CX / Rz ゲート適用
  $C2  QMCM    qb, creg        — Mid-Circuit Measurement → clbit
  $C3  QFEED   qb, creg        — Feedforward Rz (EEDT GPS core)
  $C4  QSYNC   tau_mode        — ZZ位相同期 (Lambert W τ* 計算)
  $C5  QERR    qb, N           — Separation Theorem N* チェック
  $C6  QWAIT   mode            — Poisson最適待機
  $C7  QBRANCH creg, label     — 測定結果条件分岐

タイムバジェット: $QT レジスタ (40µs = τ* window)
QOVF フラグ: バジェット超過 → Separation Theorem違反と等価
"""

import numpy as np
from scipy.special import lambertw
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
from qiskit.circuit.library import RZGate
import warnings, sys
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────
#  EEDT Hardware Constants (ibm_marrakesh Q94-Q95)
# ─────────────────────────────────────────────
TAU_CR        = 31.3e-6   # Critical time τ_CR [s]
T2            = 261e-6    # Coherence time T2 [s]
ALPHA         = 0.38      # EEDT scaling coefficient α
N_STAR_MAX    = 5         # Separation Theorem N* upper bound
TAU_STAR_MIN  = 30e-6     # Quantum-6502 window min [s]
TAU_STAR_MAX  = 50e-6     # Quantum-6502 window max [s]
GPS_EXPECTED  = +0.037    # Expected GPS gain (±0.011)
ZZ_DRIFT_KHZ  = (0.063, 3.60)  # Observed ZZ drift range [kHz]
CHI2_DOF      = 0.11      # χ²/dof from hardware fit

# ─────────────────────────────────────────────
#  Opcode Table (Woz: $C0–$C7 illegal opcode空間)
# ─────────────────────────────────────────────
OPCODE = {
    'QINIT':   0xC0, 'QGATE':   0xC1,
    'QMCM':    0xC2, 'QFEED':   0xC3,
    'QSYNC':   0xC4, 'QERR':    0xC5,
    'QWAIT':   0xC6, 'QBRANCH': 0xC7,
}

# Instruction cost table [µs] — Feynman's budget
COST_US = {
    'QINIT':   1.0,  'QGATE':   0.1,
    'QMCM':    2.0,  'QFEED':   0.1,
    'QSYNC':   3.0,  'QERR':    2.0,
    'QWAIT':   0.5,  'QBRANCH': 2.5,
}

BUDGET_US = 40.0   # τ* window budget ≈ 38µs (+ margin)

RESET  = "\033[0m"
BOLD   = "\033[1m"
CYAN   = "\033[96m"
GREEN  = "\033[92m"
YELLOW = "\033[93m"
RED    = "\033[91m"
MAGENTA= "\033[95m"
BLUE   = "\033[94m"


# ══════════════════════════════════════════════
#  EEDT 数学コア (Dirac's Lambert W engine)
# ══════════════════════════════════════════════
class EEDTCore:
    """
    EEDT理論の数学コア。
    Lambert W公式によるτ*計算、GPS利得計算、
    Separation Theorem N*チェックを実装。
    """

    @staticmethod
    def tau_star(N: int, nu_zz_khz: float = 1.0) -> float:
        """
        Lambert W公式によるτ*計算 (Dirac's QSYNC_A core)
        
        τ*(ν_ZZ, N) = τ_CR / W(α · N · exp(α))
        
        ZZ driftを引数に取る適応型 (QSYNC_A)。
        """
        arg = ALPHA * N * np.exp(ALPHA)
        W_val = float(lambertw(arg, k=0).real)
        tau_s = TAU_CR / W_val
        # ZZ drift補正: ν_ZZ が大きいほどτ*を短縮
        drift_factor = 1.0 / (1.0 + 0.01 * nu_zz_khz)
        return tau_s * drift_factor

    @staticmethod
    def gps_gain(p_with: float, p_without: float) -> float:
        """
        GPS (Gate Phase Stabilization) 利得
        GPS = P(correct|EEDT) - P(correct|no EEDT)
        """
        return p_with - p_without

    @staticmethod
    def phi_opt(tau_us: float, nu_zz_khz: float) -> float:
        """
        EEDT feedforward位相 φ* = 2π · ν_ZZ · τ
        φ*=1 (rad) が最適動作点 (Separation Property)
        """
        return 2.0 * np.pi * nu_zz_khz * 1e3 * tau_us * 1e-6

    @staticmethod
    def separation_check(N: int) -> bool:
        """
        Separation Theorem: N ≤ N* = 5
        違反 → 符号反転リスク (QOVF相当)
        """
        return N <= N_STAR_MAX

    @staticmethod
    def poisson_mle(counts: list) -> float:
        """
        Poisson MLE: λ̂ = mean(counts)
        waiting-time isomorphism: inter-MCM intervals → Poisson process
        """
        return float(np.mean(counts)) if counts else 0.0


# ══════════════════════════════════════════════
#  Quantum-6502 CPU
# ══════════════════════════════════════════════
class Quantum6502:
    """
    Quantum-6502 仮想プロセッサ
    
    クラシックレジスタ: A, X, Y, SP, PC
    量子レジスタ      : Q-reg (Qiskit QuantumRegister)
    特殊レジスタ      : $QT (タイムバジェット), $QN (MCMカウンタ)
    ステータスフラグ  : QOVF (Quantum Overflow), C (Carry), Z (Zero)
    """

    def __init__(self, n_qubits: int = 4, shots: int = 1024):
        # ── クラシックレジスタ ──
        self.A  = 0x00   # Accumulator
        self.X  = 0x00   # Index X
        self.Y  = 0x00   # Index Y
        self.SP = 0xFF   # Stack Pointer
        self.PC = 0x0000 # Program Counter

        # ── 量子リソース ──
        self.n_qubits = n_qubits
        self.shots    = shots
        self.simulator = AerSimulator()

        # ── 特殊レジスタ ──
        self.QT = BUDGET_US      # Time Budget [µs]
        self.QN = 0              # MCM実行カウンタ ($QN)

        # ── フラグ ──
        self.QOVF = False   # Quantum Overflow (N*超過 / バジェット枯渇)
        self.C    = False   # Carry
        self.Z    = False   # Zero

        # ── 現在のQiskit回路 ──
        self.qr   = QuantumRegister(n_qubits, 'q')
        self.cr   = ClassicalRegister(n_qubits, 'c')
        self.circ = QuantumCircuit(self.qr, self.cr)

        # ── EEDTコア ──
        self.eedt = EEDTCore()

        # ── 実行ログ ──
        self.exec_log   = []
        self.phase_log  = []   # feedforward位相履歴
        self.mcm_results= []   # MCM結果列

        # ── ZZ drift (Adaptive QSYNC用) ──
        self.nu_zz_khz = 1.0   # 現在のZZカップリング [kHz]

        print(f"{BOLD}{CYAN}━━━ Quantum-6502 CPU RESET ━━━{RESET}")
        print(f"  Qubits  : {n_qubits}")
        print(f"  Budget  : {BUDGET_US:.1f} µs (τ* window)")
        print(f"  N*_max  : {N_STAR_MAX}  (Separation Theorem)")
        print(f"  Shots   : {shots}")
        print(f"  τ_CR    : {TAU_CR*1e6:.1f} µs")
        print(f"  T2      : {T2*1e6:.0f} µs")
        print()

    # ─────────────────────────────────────────
    #  内部ユーティリティ
    # ─────────────────────────────────────────
    def _consume_budget(self, instr: str) -> bool:
        """タイムバジェット消費。QOVF判定付き。"""
        cost = COST_US.get(instr, 0.0)
        self.QT -= cost
        if self.QT < 0:
            self.QOVF = True
            print(f"  {RED}[QOVF] Budget exhausted at {instr}! "
                  f"Remaining: {self.QT:.2f} µs{RESET}")
            return False
        return True

    def _log(self, instr: str, detail: str = ""):
        opcode = OPCODE.get(instr, 0xFF)
        cost   = COST_US.get(instr, 0.0)
        entry  = f"  ${opcode:02X}  {instr:<8}  QT={self.QT:5.1f}µs  {detail}"
        self.exec_log.append(entry)
        print(entry)

    def _reset_circuit(self):
        """新しいQiskit回路を生成。"""
        self.qr   = QuantumRegister(self.n_qubits, 'q')
        self.cr   = ClassicalRegister(self.n_qubits, 'c')
        self.circ = QuantumCircuit(self.qr, self.cr)

    # ─────────────────────────────────────────
    #  命令実装
    # ─────────────────────────────────────────
    def QINIT(self, qb: int, state: str = '0'):
        """
        $C0 QINIT — qubit初期化
        state: '0'→|0>, '+'→|+>, '1'→|1>
        """
        if not self._consume_budget('QINIT'): return
        if state == '+':
            self.circ.h(self.qr[qb])
            detail = f"q[{qb}] ← |+>"
        elif state == '1':
            self.circ.x(self.qr[qb])
            detail = f"q[{qb}] ← |1>"
        else:
            detail = f"q[{qb}] ← |0>"
        self._log('QINIT', detail)

    def QGATE(self, qb: int, gate: str, qb2: int = None, angle: float = None):
        """
        $C1 QGATE — 量子ゲート適用
        gate: 'H','X','Y','Z','CX','Rz'
        """
        if not self._consume_budget('QGATE'): return
        g = gate.upper()
        if   g == 'H':  self.circ.h(self.qr[qb]);             detail = f"H  q[{qb}]"
        elif g == 'X':  self.circ.x(self.qr[qb]);             detail = f"X  q[{qb}]"
        elif g == 'Y':  self.circ.y(self.qr[qb]);             detail = f"Y  q[{qb}]"
        elif g == 'Z':  self.circ.z(self.qr[qb]);             detail = f"Z  q[{qb}]"
        elif g == 'CX' and qb2 is not None:
            self.circ.cx(self.qr[qb], self.qr[qb2]);          detail = f"CX q[{qb}],q[{qb2}]"
        elif g == 'RZ' and angle is not None:
            self.circ.rz(angle, self.qr[qb]);                 detail = f"Rz({angle:.4f}) q[{qb}]"
        else:
            detail = f"UNKNOWN gate {gate}"
        self._log('QGATE', detail)

    def QMCM(self, qb: int, creg: int = None):
        """
        $C2 QMCM — Mid-Circuit Measurement
        測定結果→クラシックレジスタ $creg
        MCMカウンタ QN++ / Separation Theorem チェック
        """
        if not self._consume_budget('QMCM'): return
        if creg is None: creg = qb
        self.circ.measure(self.qr[qb], self.cr[creg])
        self.QN += 1
        # Separation Theorem チェック
        if not self.eedt.separation_check(self.QN):
            self.QOVF = True
            print(f"  {RED}[QOVF] N={self.QN} > N*={N_STAR_MAX}: "
                  f"Separation Theorem violated!{RESET}")
        self._log('QMCM', f"measure q[{qb}]→c[{creg}]  QN={self.QN}")

    def QFEED(self, qb: int, creg: int, nu_zz_khz: float = None,
              tau_us: float = 38.0):
        """
        $C3 QFEED — EEDT Feedforward Rz
        
        φ* = 2π · ν_ZZ · τ を計算し、測定結果に応じてRz適用。
        これがEEDTのGPSゲインの実体。
        Lambert W公式内包 (Dirac's requirement)。
        """
        if not self._consume_budget('QFEED'): return
        if nu_zz_khz is None:
            nu_zz_khz = self.nu_zz_khz

        # Lambert W → τ* 計算
        tau_star_us = self.eedt.tau_star(
            max(1, self.QN), nu_zz_khz) * 1e6

        # φ* = 2π · ν_ZZ · τ
        phi = self.eedt.phi_opt(tau_us, nu_zz_khz)
        self.phase_log.append(phi)

        # 条件付きRz: c_if (creg==1) → Rz(φ*)
        with self.circ.if_test((self.cr[creg], 1)):
            self.circ.rz(phi, self.qr[qb])

        detail = (f"Rz(φ*={phi:.4f}rad) q[{qb}] | c[{creg}]  "
                  f"τ*={tau_star_us:.1f}µs  ν_ZZ={nu_zz_khz:.3f}kHz")
        self._log('QFEED', detail)

    def QSYNC(self, tau_mode: str = 'adaptive'):
        """
        $C4 QSYNC — ZZ位相同期
        
        mode 'static'   : 固定τ* (QSYNC_S, 高速)
        mode 'adaptive' : Lambert W + ZZ drift補正 (QSYNC_A)
        """
        if not self._consume_budget('QSYNC'): return
        if tau_mode == 'adaptive':
            tau_s = self.eedt.tau_star(max(1, self.QN), self.nu_zz_khz)
            mode_str = f"ADAPTIVE τ*={tau_s*1e6:.1f}µs"
        else:
            tau_s = TAU_CR
            mode_str = f"STATIC   τ*={tau_s*1e6:.1f}µs"

        in_window = TAU_STAR_MIN <= tau_s <= TAU_STAR_MAX
        window_str = (f"{GREEN}IN Quantum-6502 window ✓{RESET}"
                      if in_window else
                      f"{YELLOW}OUT of window ⚠{RESET}")
        self._log('QSYNC', f"{mode_str}  {window_str}")

    def QERR(self, qb: int, N: int = None):
        """
        $C5 QERR — Separation Theorem Error Mitigation
        N回MCMをタイムバジェット内でループ実行する前に
        N* チェックをかける安全ガード。
        """
        if not self._consume_budget('QERR'): return
        check_N = N if N is not None else self.QN
        ok = self.eedt.separation_check(check_N)
        status = (f"{GREEN}SAFE N={check_N}≤{N_STAR_MAX}{RESET}"
                  if ok else
                  f"{RED}RISK N={check_N}>{N_STAR_MAX} — sign flip risk!{RESET}")
        self._log('QERR', status)

    def QWAIT(self, mode: str = 'poisson'):
        """
        $C6 QWAIT — Poisson最適待機
        inter-MCM intervals → Poisson MLE → λ̂ = mean
        """
        if not self._consume_budget('QWAIT'): return
        lam = self.eedt.poisson_mle(self.mcm_results)
        self._log('QWAIT', f"Poisson MLE λ̂={lam:.4f}  mode={mode}")

    def QBRANCH(self, creg: int, label: str, expected: int = 1):
        """
        $C7 QBRANCH — 測定結果条件分岐 (ログのみ、制御フロー記録)
        """
        if not self._consume_budget('QBRANCH'): return
        self._log('QBRANCH', f"if c[{creg}]=={expected} → {label}")

    # ─────────────────────────────────────────
    #  回路実行 & GPS利得計算
    # ─────────────────────────────────────────
    def execute(self, label: str = "") -> dict:
        """現在の回路をAerで実行し、カウント辞書を返す。"""
        job    = self.simulator.run(self.circ, shots=self.shots)
        result = job.result()
        counts = result.get_counts(self.circ)
        if label:
            print(f"\n  {BLUE}[RUN] {label}: {counts}{RESET}")
        return counts

    def measure_gps(self, shots: int = None) -> dict:
        """
        EEDT GPS利得を測定。
        with EEDT (QFEED有り) と without (QFEED無し) を比較。
        返り値: {'gps': float, 'p_with': float, 'p_without': float}
        """
        if shots is None: shots = self.shots

        # ── WITHOUT EEDT ──
        qr0 = QuantumRegister(2, 'q')
        cr0 = ClassicalRegister(2, 'c')
        c0  = QuantumCircuit(qr0, cr0)
        c0.h(qr0[0])
        c0.cx(qr0[0], qr0[1])
        c0.measure(qr0, cr0)
        r0     = self.simulator.run(c0, shots=shots).result().get_counts(c0)
        p0_00  = r0.get('00', 0) / shots   # |00> 期待値
        p0_11  = r0.get('11', 0) / shots   # |11> 期待値
        p_without = p0_00 + p0_11          # Bell state正解確率

        # ── WITH EEDT (QMCM → QFEED パイプライン) ──
        qr1 = QuantumRegister(2, 'q')
        cr1 = ClassicalRegister(2, 'c')
        c1  = QuantumCircuit(qr1, cr1)
        c1.h(qr1[0])
        c1.cx(qr1[0], qr1[1])
        # MCM on qubit 0
        c1.measure(qr1[0], cr1[0])
        # Feedforward Rz on qubit 1 based on MCM result
        phi = self.eedt.phi_opt(38.0, self.nu_zz_khz)
        with c1.if_test((cr1[0], 1)):
            c1.rz(phi, qr1[1])
        c1.measure(qr1, cr1)
        r1     = self.simulator.run(c1, shots=shots).result().get_counts(c1)
        p1_00  = r1.get('00', 0) / shots
        p1_11  = r1.get('11', 0) / shots
        p_with = p1_00 + p1_11

        gps = self.eedt.gps_gain(p_with, p_without)
        return {
            'gps':       gps,
            'p_with':    p_with,
            'p_without': p_without,
            'counts_without': r0,
            'counts_with':    r1,
        }

    def status(self):
        """CPUステータス表示。"""
        flag_qovf = f"{RED}QOVF{RESET}" if self.QOVF else f"{GREEN}----{RESET}"
        flag_c    = f"C={'1' if self.C else '0'}"
        flag_z    = f"Z={'1' if self.Z else '0'}"
        budget_bar = int(self.QT / BUDGET_US * 20)
        bar = f"[{'█'*budget_bar}{'░'*(20-budget_bar)}]"
        pct = self.QT / BUDGET_US * 100

        print(f"\n  {BOLD}── CPU Status ──────────────────────────{RESET}")
        print(f"  A={self.A:02X}  X={self.X:02X}  Y={self.Y:02X}  "
              f"SP={self.SP:02X}  PC={self.PC:04X}")
        print(f"  Flags: {flag_qovf}  {flag_c}  {flag_z}")
        print(f"  $QT (Budget): {bar} {self.QT:5.1f}µs / {BUDGET_US:.0f}µs  ({pct:.0f}%)")
        print(f"  $QN (MCM count): {self.QN}  (N*={N_STAR_MAX})")
        print(f"  ν_ZZ: {self.nu_zz_khz:.3f} kHz  "
              f"(drift range: {ZZ_DRIFT_KHZ[0]}–{ZZ_DRIFT_KHZ[1]} kHz)")
        print()


# ══════════════════════════════════════════════
#  Demo Programs
# ══════════════════════════════════════════════

def demo_qmcm_qfeed_pipeline():
    """
    Demo 1: QMCM → QFEED コアパイプライン
    
    EEDTのGPSゲインの「2命令マイクロプログラム」表現。
    合議結論: 「最も実装価値が高い最初のターゲット」
    """
    print(f"\n{BOLD}{MAGENTA}{'='*60}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 1: QMCM → QFEED Pipeline{RESET}")
    print(f"{BOLD}{MAGENTA}  (EEDT GPS Core as 2-instruction Microprogram){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*60}{RESET}\n")

    cpu = Quantum6502(n_qubits=4, shots=2048)
    cpu.nu_zz_khz = 1.5   # ZZ drift設定

    print(f"{BOLD}── Program: Bell State + EEDT Stabilization ──{RESET}")
    print(f"  (Woz: opcodes $C0→$C1→$C1→$C4→$C2→$C3→$C5)")
    print()

    # ── プログラム実行 ──
    cpu.QINIT(0, '0')                      # $C0: q[0] ← |0>
    cpu.QINIT(1, '0')                      # $C0: q[1] ← |0>
    cpu.QGATE(0, 'H')                      # $C1: H q[0]
    cpu.QGATE(0, 'CX', qb2=1)             # $C1: CX q[0],q[1] → Bell |Φ+>
    cpu.QSYNC('adaptive')                  # $C4: Lambert W τ* sync
    cpu.QMCM(0, creg=0)                   # $C2: MCM q[0] → c[0]
    cpu.QFEED(1, creg=0,                   # $C3: feedforward Rz q[1]
              nu_zz_khz=cpu.nu_zz_khz,
              tau_us=38.0)
    cpu.QERR(1, N=cpu.QN)                 # $C5: Separation check
    cpu.QWAIT('poisson')                   # $C6: Poisson wait

    cpu.status()

    # ── GPS利得測定 ──
    print(f"{BOLD}── GPS Gain Measurement ──{RESET}")
    gps_result = cpu.measure_gps()

    gps   = gps_result['gps']
    p_w   = gps_result['p_with']
    p_wo  = gps_result['p_without']
    color = GREEN if gps > 0 else RED

    print(f"\n  P(correct | WITHOUT EEDT) = {p_wo:.4f}")
    print(f"  P(correct |    WITH EEDT) = {p_w:.4f}")
    print(f"  {color}{BOLD}GPS = {gps:+.4f}{RESET}  "
          f"(expected ≈ +{GPS_EXPECTED:.3f})")

    sigma = abs(gps) / 0.011
    print(f"  Significance ≈ {sigma:.1f}σ  (EEDT paper: z=8.5σ)")
    print(f"  Counts (without): {gps_result['counts_without']}")
    print(f"  Counts (with):    {gps_result['counts_with']}")

    return gps_result


def demo_zz_drift_scan():
    """
    Demo 2: ZZ Drift スキャン
    
    ν_ZZ を 0.063〜3.60 kHz (57× range) でスキャンし、
    τ* と φ* の変化を観測。
    「新物理の兆候」: φ*=1 からの乖離がEEDT崩壊の予兆。
    """
    print(f"\n{BOLD}{MAGENTA}{'='*60}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 2: ZZ Drift Scan (τ* adaptive tracking){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*60}{RESET}\n")

    eedt = EEDTCore()
    nu_range = np.linspace(ZZ_DRIFT_KHZ[0], ZZ_DRIFT_KHZ[1], 10)
    N_test   = 3  # N* = 3 (within Separation Theorem)

    print(f"  {'ν_ZZ [kHz]':>12}  {'τ* [µs]':>10}  "
          f"{'φ* [rad]':>10}  {'In window':>10}  {'φ*≈1?':>8}")
    print(f"  {'─'*12}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*8}")

    for nu in nu_range:
        tau_s  = eedt.tau_star(N_test, nu) * 1e6
        phi    = eedt.phi_opt(tau_s, nu)
        in_win = TAU_STAR_MIN*1e6 <= tau_s <= TAU_STAR_MAX*1e6
        phi_ok = 0.8 <= phi <= 1.2   # φ*=1 ± 20%

        win_str = f"{GREEN}✓{RESET}" if in_win else f"{YELLOW}⚠{RESET}"
        phi_str = f"{GREEN}✓{RESET}" if phi_ok  else f"{RED}✗{RESET}"
        print(f"  {nu:>12.3f}  {tau_s:>10.2f}  "
              f"{phi:>10.4f}  {win_str:>18}  {phi_str:>16}")

    print(f"\n  {YELLOW}Note: φ* >> 1 at high ν_ZZ → potential 'new physics' signal{RESET}")
    print(f"  {YELLOW}(EEDT theory predicts GPS flip if φ* ≫ 2π){RESET}")


def demo_separation_theorem_stress():
    """
    Demo 3: Separation Theorem 境界テスト
    
    N = 1..8 でQMCMを繰り返し、N* = 5 でQOVFが
    発火することをシミュレート。
    Feynman: 「N*超過 = QOVF = 符号反転リスク」
    """
    print(f"\n{BOLD}{MAGENTA}{'='*60}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 3: Separation Theorem Stress Test{RESET}")
    print(f"{BOLD}{MAGENTA}  (N* boundary: QOVF fires at N > {N_STAR_MAX}){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*60}{RESET}\n")

    cpu = Quantum6502(n_qubits=4, shots=512)

    print(f"  Running QMCM loop N=1..8  (N*={N_STAR_MAX}):\n")

    for i in range(8):
        if cpu.QOVF:
            print(f"  {RED}  └─ QOVF active: halting loop (simulates CPU trap){RESET}")
            break
        cpu.QMCM(0, creg=0)
        cpu.QERR(0, N=cpu.QN)

    cpu.status()


def demo_quantum_6502_full():
    """
    Demo 4: Quantum-6502 フルプログラム
    
    Bell state生成 → ZZ drift適応 → EEDT安定化
    → GPS測定 → QBRANCH（結果分岐）
    全命令セット ($C0–$C7) を使ったデモ。
    """
    print(f"\n{BOLD}{MAGENTA}{'='*60}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 4: Full Quantum-6502 Program{RESET}")
    print(f"{BOLD}{MAGENTA}  (All opcodes $C0–$C7){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*60}{RESET}\n")

    cpu = Quantum6502(n_qubits=4, shots=1024)
    cpu.nu_zz_khz = 0.8   # 低ZZ drift (安定動作域)

    print(f"{BOLD}── Quantum-6502 Assembly ──{RESET}")
    print(f"  ; Bell state + EEDT feedback + conditional branch")
    print()

    cpu.QINIT(0, '+')                         # $C0  QINIT  q0, |+>
    cpu.QINIT(1, '0')                         # $C0  QINIT  q1, |0>
    cpu.QGATE(0, 'CX', qb2=1)                # $C1  QGATE  CX q0,q1
    cpu.QSYNC('adaptive')                     # $C4  QSYNC  adaptive
    cpu.QMCM(0, creg=0)                      # $C2  QMCM   q0 → c0
    cpu.QFEED(1, creg=0,
              nu_zz_khz=cpu.nu_zz_khz,
              tau_us=38.0)                    # $C3  QFEED  q1,c0
    cpu.QMCM(1, creg=1)                      # $C2  QMCM   q1 → c1
    cpu.QERR(1, N=cpu.QN)                    # $C5  QERR   N check
    cpu.QWAIT('poisson')                      # $C6  QWAIT  poisson
    cpu.QBRANCH(creg=1, label='EEDT_PASS',
                expected=0)                   # $C7  QBRANCH c1==0 → PASS

    cpu.status()

    # GPS測定
    print(f"{BOLD}── Final GPS Measurement ──{RESET}")
    g = cpu.measure_gps()
    gps = g['gps']
    c = GREEN if gps > 0 else RED
    print(f"\n  {c}{BOLD}GPS = {gps:+.4f}{RESET}  "
          f"({'PASS ✓' if gps > 0 else 'FAIL ✗'})")

    return cpu


# ══════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════
if __name__ == '__main__':
    print(f"\n{BOLD}{CYAN}{'━'*60}{RESET}")
    print(f"{BOLD}{CYAN}  Quantum-6502 Prototype v0.1{RESET}")
    print(f"{BOLD}{CYAN}  EEDT-Native Hybrid Quantum Processor{RESET}")
    print(f"{BOLD}{CYAN}  Council: Dirac + Woz + Feynman{RESET}")
    print(f"{BOLD}{CYAN}{'━'*60}{RESET}")

    print(f"\n{BOLD}EEDT Constants:{RESET}")
    print(f"  τ_CR  = {TAU_CR*1e6:.1f} µs  |  T2 = {T2*1e6:.0f} µs")
    print(f"  α     = {ALPHA}          |  GPS_expected = +{GPS_EXPECTED:.3f}")
    print(f"  N*    = {N_STAR_MAX}              |  τ* window = [{TAU_STAR_MIN*1e6:.0f}, {TAU_STAR_MAX*1e6:.0f}] µs")
    print(f"  ZZ drift = {ZZ_DRIFT_KHZ[0]}–{ZZ_DRIFT_KHZ[1]} kHz (57× range, 6-day measurement)")

    # ── Lambert W τ* preview ──
    print(f"\n{BOLD}Lambert W τ*(N) Preview:{RESET}")
    eedt = EEDTCore()
    for N in [1, 2, 3, 4, 5]:
        tau_s = eedt.tau_star(N) * 1e6
        in_w  = TAU_STAR_MIN*1e6 <= tau_s <= TAU_STAR_MAX*1e6
        mark  = f"{GREEN}← Quantum-6502 window{RESET}" if in_w else ""
        print(f"  N={N}: τ* = {tau_s:.2f} µs  {mark}")

    # ── Run demos ──
    r1 = demo_qmcm_qfeed_pipeline()
    demo_zz_drift_scan()
    demo_separation_theorem_stress()
    demo_quantum_6502_full()

    print(f"\n{BOLD}{CYAN}{'━'*60}{RESET}")
    print(f"{BOLD}{GREEN}  Quantum-6502 Prototype: All demos complete ✓{RESET}")
    print(f"{BOLD}{CYAN}{'━'*60}{RESET}\n")
