"""
Quantum-6502 Prototype v0.1
EEDT-Native Hybrid Quantum Processor Emulator

Council: Dirac (理論) + Woz (アーキテクチャ) + Feynman (実機)
Author : Takeshi Okuda (okudat9) — Independent Quantum Researcher
Base   : EEDT v19, Zenodo 10.5281/zenodo.18899504
Env    : Windows / Python 3.10+ / Qiskit 2.x / Qiskit-Aer
"""

import numpy as np
from scipy.special import lambertw
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
import warnings
warnings.filterwarnings('ignore')

TAU_CR        = 31.3e-6
T2            = 261e-6
ALPHA         = 0.38
N_STAR_MAX    = 5
TAU_STAR_MIN  = 30e-6
TAU_STAR_MAX  = 50e-6
GPS_EXPECTED  = +0.037
ZZ_DRIFT_KHZ  = (0.063, 3.60)
CHI2_DOF      = 0.11

OPCODE = {
    'QINIT':0xC0,'QGATE':0xC1,'QMCM':0xC2,'QFEED':0xC3,
    'QSYNC':0xC4,'QERR':0xC5,'QWAIT':0xC6,'QBRANCH':0xC7,
}
COST_US = {
    'QINIT':1.0,'QGATE':0.1,'QMCM':2.0,'QFEED':0.1,
    'QSYNC':3.0,'QERR':2.0,'QWAIT':0.5,'QBRANCH':2.5,
}
BUDGET_US = 40.0

RESET="\033[0m"; BOLD="\033[1m"; CYAN="\033[96m"
GREEN="\033[92m"; YELLOW="\033[93m"; RED="\033[91m"
MAGENTA="\033[95m"; BLUE="\033[94m"

class EEDTCore:
    @staticmethod
    def tau_star(N,nu_zz_khz=1.0):
        arg=ALPHA*N*np.exp(ALPHA)
        W_val=float(lambertw(arg,k=0).real)
        tau_s=TAU_CR/W_val
        drift_factor=1.0/(1.0+0.01*nu_zz_khz)
        return tau_s*drift_factor

    @staticmethod
    def gps_gain(p_with,p_without): return p_with-p_without

    @staticmethod
    def phi_opt(tau_us,nu_zz_khz):
        return 2.0*np.pi*nu_zz_khz*1e3*tau_us*1e-6

    @staticmethod
    def separation_check(N): return N<=N_STAR_MAX

    @staticmethod
    def poisson_mle(counts):
        return float(np.mean(counts)) if counts else 0.0

class Quantum6502:
    def __init__(self,n_qubits=4,shots=1024):
        self.A=0x00; self.X=0x00; self.Y=0x00
        self.SP=0xFF; self.PC=0x0000
        self.n_qubits=n_qubits; self.shots=shots
        self.simulator=AerSimulator()
        self.QT=BUDGET_US; self.QN=0
        self.QOVF=False; self.C=False; self.Z=False
        self.qr=QuantumRegister(n_qubits,'q')
        self.cr=ClassicalRegister(n_qubits,'c')
        self.circ=QuantumCircuit(self.qr,self.cr)
        self.eedt=EEDTCore()
        self.exec_log=[]; self.phase_log=[]; self.mcm_results=[]
        self.nu_zz_khz=1.0
        print(f"{BOLD}{CYAN}━━━ Quantum-6502 CPU RESET ━━━{RESET}")
        print(f"  Qubits:{n_qubits}  Budget:{BUDGET_US:.1f}µs  N*:{N_STAR_MAX}  Shots:{shots}")
        print(f"  τ_CR:{TAU_CR*1e6:.1f}µs  T2:{T2*1e6:.0f}µs\n")

    def _consume_budget(self,instr):
        cost=COST_US.get(instr,0.0); self.QT-=cost
        if self.QT<0:
            self.QOVF=True
            print(f"  {RED}[QOVF] Budget exhausted at {instr}! {self.QT:.2f}µs{RESET}")
            return False
        return True

    def _log(self,instr,detail=""):
        opcode=OPCODE.get(instr,0xFF)
        entry=f"  ${opcode:02X}  {instr:<8}  QT={self.QT:5.1f}µs  {detail}"
        self.exec_log.append(entry); print(entry)

    def _reset_circuit(self):
        self.qr=QuantumRegister(self.n_qubits,'q')
        self.cr=ClassicalRegister(self.n_qubits,'c')
        self.circ=QuantumCircuit(self.qr,self.cr)

    def QINIT(self,qb,state='0'):
        if not self._consume_budget('QINIT'): return
        if state=='+': self.circ.h(self.qr[qb]); detail=f"q[{qb}] ← |+>"
        elif state=='1': self.circ.x(self.qr[qb]); detail=f"q[{qb}] ← |1>"
        else: detail=f"q[{qb}] ← |0>"
        self._log('QINIT',detail)

    def QGATE(self,qb,gate,qb2=None,angle=None):
        if not self._consume_budget('QGATE'): return
        g=gate.upper()
        if g=='H': self.circ.h(self.qr[qb]); detail=f"H q[{qb}]"
        elif g=='X': self.circ.x(self.qr[qb]); detail=f"X q[{qb}]"
        elif g=='Y': self.circ.y(self.qr[qb]); detail=f"Y q[{qb}]"
        elif g=='Z': self.circ.z(self.qr[qb]); detail=f"Z q[{qb}]"
        elif g=='CX' and qb2 is not None:
            self.circ.cx(self.qr[qb],self.qr[qb2]); detail=f"CX q[{qb}],q[{qb2}]"
        elif g=='RZ' and angle is not None:
            self.circ.rz(angle,self.qr[qb]); detail=f"Rz({angle:.4f}) q[{qb}]"
        else: detail=f"UNKNOWN {gate}"
        self._log('QGATE',detail)

    def QMCM(self,qb,creg=None):
        if not self._consume_budget('QMCM'): return
        if creg is None: creg=qb
        self.circ.measure(self.qr[qb],self.cr[creg])
        self.QN+=1
        if not self.eedt.separation_check(self.QN):
            self.QOVF=True
            print(f"  {RED}[QOVF] N={self.QN}>{N_STAR_MAX}: Separation Theorem violated!{RESET}")
        self._log('QMCM',f"measure q[{qb}]→c[{creg}]  QN={self.QN}")

    def QFEED(self,qb,creg,nu_zz_khz=None,tau_us=38.0):
        if not self._consume_budget('QFEED'): return
        if nu_zz_khz is None: nu_zz_khz=self.nu_zz_khz
        tau_star_us=self.eedt.tau_star(max(1,self.QN),nu_zz_khz)*1e6
        phi=self.eedt.phi_opt(tau_us,nu_zz_khz)
        self.phase_log.append(phi)
        with self.circ.if_test((self.cr[creg],1)):
            self.circ.rz(phi,self.qr[qb])
        self._log('QFEED',f"Rz(φ*={phi:.4f}rad) q[{qb}]|c[{creg}]  τ*={tau_star_us:.1f}µs  ν_ZZ={nu_zz_khz:.3f}kHz")

    def QSYNC(self,tau_mode='adaptive'):
        if not self._consume_budget('QSYNC'): return
        if tau_mode=='adaptive':
            tau_s=self.eedt.tau_star(max(1,self.QN),self.nu_zz_khz)
            mode_str=f"ADAPTIVE τ*={tau_s*1e6:.1f}µs"
        else:
            tau_s=TAU_CR; mode_str=f"STATIC τ*={tau_s*1e6:.1f}µs"
        in_window=TAU_STAR_MIN<=tau_s<=TAU_STAR_MAX
        win=f"{GREEN}IN window ✓{RESET}" if in_window else f"{YELLOW}OUT of window ⚠{RESET}"
        self._log('QSYNC',f"{mode_str}  {win}")

    def QERR(self,qb,N=None):
        if not self._consume_budget('QERR'): return
        check_N=N if N is not None else self.QN
        ok=self.eedt.separation_check(check_N)
        s=(f"{GREEN}SAFE N={check_N}≤{N_STAR_MAX}{RESET}" if ok
           else f"{RED}RISK N={check_N}>{N_STAR_MAX} sign flip!{RESET}")
        self._log('QERR',s)

    def QWAIT(self,mode='poisson'):
        if not self._consume_budget('QWAIT'): return
        lam=self.eedt.poisson_mle(self.mcm_results)
        self._log('QWAIT',f"Poisson MLE λ̂={lam:.4f}  mode={mode}")

    def QBRANCH(self,creg,label,expected=1):
        if not self._consume_budget('QBRANCH'): return
        self._log('QBRANCH',f"if c[{creg}]=={expected} → {label}")

    def execute(self,label=""):
        job=self.simulator.run(self.circ,shots=self.shots)
        result=job.result(); counts=result.get_counts(self.circ)
        if label: print(f"\n  {BLUE}[RUN] {label}: {counts}{RESET}")
        return counts

    def measure_gps(self,shots=None):
        if shots is None: shots=self.shots
        qr0=QuantumRegister(2,'q'); cr0=ClassicalRegister(2,'c')
        c0=QuantumCircuit(qr0,cr0)
        c0.h(qr0[0]); c0.cx(qr0[0],qr0[1]); c0.measure(qr0,cr0)
        r0=self.simulator.run(c0,shots=shots).result().get_counts(c0)
        p_without=(r0.get('00',0)+r0.get('11',0))/shots

        qr1=QuantumRegister(2,'q'); cr1=ClassicalRegister(2,'c')
        c1=QuantumCircuit(qr1,cr1)
        c1.h(qr1[0]); c1.cx(qr1[0],qr1[1])
        c1.measure(qr1[0],cr1[0])
        phi=self.eedt.phi_opt(38.0,self.nu_zz_khz)
        with c1.if_test((cr1[0],1)): c1.rz(phi,qr1[1])
        c1.measure(qr1,cr1)
        r1=self.simulator.run(c1,shots=shots).result().get_counts(c1)
        p_with=(r1.get('00',0)+r1.get('11',0))/shots
        gps=self.eedt.gps_gain(p_with,p_without)
        return {'gps':gps,'p_with':p_with,'p_without':p_without,
                'counts_without':r0,'counts_with':r1}

    def status(self):
        flag_qovf=f"{RED}QOVF{RESET}" if self.QOVF else f"{GREEN}----{RESET}"
        budget_bar=max(0,int(self.QT/BUDGET_US*20))
        bar=f"[{'█'*budget_bar}{'░'*(20-budget_bar)}]"
        pct=self.QT/BUDGET_US*100
        print(f"\n  {BOLD}── CPU Status ──{RESET}")
        print(f"  A={self.A:02X} X={self.X:02X} Y={self.Y:02X} SP={self.SP:02X} PC={self.PC:04X}")
        print(f"  Flags:{flag_qovf} C={'1' if self.C else '0'} Z={'1' if self.Z else '0'}")
        print(f"  $QT:{bar} {self.QT:5.1f}µs/{BUDGET_US:.0f}µs ({pct:.0f}%)")
        print(f"  $QN:{self.QN} (N*={N_STAR_MAX})  ν_ZZ:{self.nu_zz_khz:.3f}kHz\n")

def demo_qmcm_qfeed_pipeline():
    print(f"\n{BOLD}{MAGENTA}{'='*55}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 1: QMCM→QFEED Pipeline (EEDT GPS Core){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*55}{RESET}\n")
    cpu=Quantum6502(n_qubits=4,shots=2048); cpu.nu_zz_khz=1.5
    cpu.QINIT(0,'0'); cpu.QINIT(1,'0')
    cpu.QGATE(0,'H'); cpu.QGATE(0,'CX',qb2=1)
    cpu.QSYNC('adaptive'); cpu.QMCM(0,creg=0)
    cpu.QFEED(1,creg=0,nu_zz_khz=cpu.nu_zz_khz,tau_us=38.0)
    cpu.QERR(1,N=cpu.QN); cpu.QWAIT('poisson')
    cpu.status()
    print(f"{BOLD}── GPS Gain ──{RESET}")
    g=cpu.measure_gps()
    gps=g['gps']; c=GREEN if gps>0 else RED
    print(f"  P(without EEDT)={g['p_without']:.4f}")
    print(f"  P(with    EEDT)={g['p_with']:.4f}")
    print(f"  {c}{BOLD}GPS={gps:+.4f}{RESET}  (expected≈+{GPS_EXPECTED:.3f})")
    print(f"  Significance≈{abs(gps)/0.011:.1f}σ")
    return g

def demo_zz_drift_scan():
    print(f"\n{BOLD}{MAGENTA}{'='*55}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 2: ZZ Drift Scan (τ* adaptive){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*55}{RESET}\n")
    eedt=EEDTCore()
    nu_range=np.linspace(ZZ_DRIFT_KHZ[0],ZZ_DRIFT_KHZ[1],10)
    print(f"  {'ν_ZZ[kHz]':>10}  {'τ*[µs]':>8}  {'φ*[rad]':>8}  {'window':>8}  {'φ*≈1':>6}")
    print(f"  {'─'*10}  {'─'*8}  {'─'*8}  {'─'*8}  {'─'*6}")
    for nu in nu_range:
        tau_s=eedt.tau_star(3,nu)*1e6
        phi=eedt.phi_opt(tau_s,nu)
        in_win=TAU_STAR_MIN*1e6<=tau_s<=TAU_STAR_MAX*1e6
        phi_ok=0.8<=phi<=1.2
        ws=f"{GREEN}✓{RESET}" if in_win else f"{YELLOW}⚠{RESET}"
        ps=f"{GREEN}✓{RESET}" if phi_ok else f"{RED}✗{RESET}"
        print(f"  {nu:>10.3f}  {tau_s:>8.2f}  {phi:>8.4f}  {ws:>16}  {ps:>14}")

def demo_separation_theorem_stress():
    print(f"\n{BOLD}{MAGENTA}{'='*55}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 3: Separation Theorem Stress (N*={N_STAR_MAX}){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*55}{RESET}\n")
    cpu=Quantum6502(n_qubits=4,shots=512)
    for i in range(8):
        if cpu.QOVF:
            print(f"  {RED}└─ QOVF: halting (CPU trap){RESET}"); break
        cpu.QMCM(0,creg=0); cpu.QERR(0,N=cpu.QN)
    cpu.status()

def demo_quantum_6502_full():
    print(f"\n{BOLD}{MAGENTA}{'='*55}{RESET}")
    print(f"{BOLD}{MAGENTA}  Demo 4: Full Quantum-6502 ($C0–$C7){RESET}")
    print(f"{BOLD}{MAGENTA}{'='*55}{RESET}\n")
    cpu=Quantum6502(n_qubits=4,shots=1024); cpu.nu_zz_khz=0.8
    cpu.QINIT(0,'+'); cpu.QINIT(1,'0')
    cpu.QGATE(0,'CX',qb2=1); cpu.QSYNC('adaptive')
    cpu.QMCM(0,creg=0)
    cpu.QFEED(1,creg=0,nu_zz_khz=cpu.nu_zz_khz,tau_us=38.0)
    cpu.QMCM(1,creg=1); cpu.QERR(1,N=cpu.QN)
    cpu.QWAIT('poisson'); cpu.QBRANCH(creg=1,label='EEDT_PASS',expected=0)
    cpu.status()
    g=cpu.measure_gps(); gps=g['gps']; c=GREEN if gps>0 else RED
    print(f"  {c}{BOLD}GPS={gps:+.4f}  ({'PASS ✓' if gps>0 else 'FAIL ✗'}){RESET}")
    return cpu

if __name__=='__main__':
    print(f"\n{BOLD}{CYAN}{'━'*55}{RESET}")
    print(f"{BOLD}{CYAN}  Quantum-6502 Prototype v0.1{RESET}")
    print(f"{BOLD}{CYAN}  EEDT-Native Hybrid Quantum Processor{RESET}")
    print(f"{BOLD}{CYAN}{'━'*55}{RESET}")
    print(f"\n{BOLD}EEDT Constants:{RESET}")
    print(f"  τ_CR={TAU_CR*1e6:.1f}µs  T2={T2*1e6:.0f}µs  α={ALPHA}  N*={N_STAR_MAX}")
    print(f"  GPS_expected=+{GPS_EXPECTED:.3f}  ZZ drift={ZZ_DRIFT_KHZ[0]}–{ZZ_DRIFT_KHZ[1]}kHz")
    print(f"\n{BOLD}Lambert W τ*(N) Preview:{RESET}")
    eedt=EEDTCore()
    for N in [1,2,3,4,5]:
        tau_s=eedt.tau_star(N)*1e6
        in_w=TAU_STAR_MIN*1e6<=tau_s<=TAU_STAR_MAX*1e6
        mark=f"{GREEN}← Quantum-6502 window{RESET}" if in_w else ""
        print(f"  N={N}: τ*={tau_s:.2f}µs  {mark}")
    demo_qmcm_qfeed_pipeline()
    demo_zz_drift_scan()
    demo_separation_theorem_stress()
    demo_quantum_6502_full()
    print(f"\n{BOLD}{CYAN}{'━'*55}{RESET}")
    print(f"{BOLD}{GREEN}  Quantum-6502: All demos complete ✓{RESET}")
    print(f"{BOLD}{CYAN}{'━'*55}{RESET}\n")
