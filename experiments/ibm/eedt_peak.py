"""
EEDT τ*精密化スキャン
τ=80, 85, 90µs / N=1
τ*=70µsのピーク確定 → φ* ∝ 1/T2の精密検証
Q94-Q95 / 2026-03-15
"""
import numpy as np, json, os
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

TS  = datetime.now().strftime("%Y%m%d_%H%M%S")
DIR = f"eedt_peak_{TS}"
os.makedirs(DIR, exist_ok=True)

service = QiskitRuntimeService(
    channel  = "ibm_quantum_platform",
    token    = "YOUR_IBM_QUANTUM_TOKEN",
    instance = "YOUR_IBM_QUANTUM_CRN"
)
backend = service.backend("ibm_marrakesh")
print(f"[{TS}] Connected")

ANC, TGT = 94, 95
NU_ZZ    = 3.6
SHOTS    = 8000
omega_zz = 2*np.pi*NU_ZZ*1e3

pm = generate_preset_pass_manager(
    backend=backend, optimization_level=1,
    initial_layout=[ANC, TGT]
)

TAU_NEW = [80, 85, 90]
circs   = []

for tau_us in TAU_NEW:
    # WITH EEDT N=1
    qr = QuantumRegister(2,'q'); cr = ClassicalRegister(2,'c')
    c  = QuantumCircuit(qr, cr)
    c.h(qr[0]); c.h(qr[1])
    c.delay(int(tau_us*1000), qr[0], unit='ns')
    c.measure(qr[0], cr[0])
    with c.if_test((cr[0],1)):
        c.rz(-omega_zz*tau_us*1e-6, qr[1])
    c.h(qr[1]); c.measure(qr[1], cr[1])
    circs.append(('with', tau_us, c))

    # REF (ancilla |+>)
    qr2 = QuantumRegister(2,'q'); cr2 = ClassicalRegister(1,'c')
    c2  = QuantumCircuit(qr2, cr2)
    c2.h(qr2[0]); c2.h(qr2[1])
    c2.delay(int(tau_us*1000), qr2[1], unit='ns')
    c2.h(qr2[1]); c2.measure(qr2[1], cr2[0])
    circs.append(('ref', tau_us, c2))

print(f"  {len(circs)}回路 (τ=80,85,90µs × with/ref)")

isa = [pm.run(c) for _,_,c in circs]
job = SamplerV2(mode=backend).run(isa, shots=SHOTS)
print(f"  Job: {job.job_id()}")
with open(f"{DIR}/jobid.txt","w") as f: f.write(job.job_id())

print("  Waiting...")
res = job.result()

pairs = {}
for i,(kind,tau_us,_) in enumerate(circs):
    cts = res[i].data.c.get_counts()
    if kind == 'with':
        p0 = sum(v for k,v in cts.items() if k[0]=='0') / SHOTS
    else:
        p0 = sum(v for k,v in cts.items() if k[-1]=='0') / SHOTS
    pairs.setdefault(tau_us,{})[kind] = p0

sigma   = 0.011
new_res = {}
for tau_us in TAU_NEW:
    d = pairs.get(tau_us,{})
    if 'with' in d and 'ref' in d:
        gps = d['with']-d['ref']
        new_res[tau_us] = {'gps':gps,'z':gps/sigma,
                           'p_with':d['with'],'p_ref':d['ref']}

# 全データ統合
all_data = {
     5: 0.022,  10: 0.058,  15: 0.113,
    20: 0.183,  30: 0.326,  33: 0.371,
    50: 0.588,  70: 0.684, 100: 0.250,
}
for t,v in new_res.items():
    all_data[t] = v['gps']

print(f"\n  τ(µs)   φ=ωτ    GPS      判定")
print(f"  ─────  ──────  ───────  ────")
peak_gps, peak_tau = 0, 0
for t in sorted(all_data.keys()):
    g   = all_data[t]
    phi = omega_zz*t*1e-6
    src = '★NEW' if t in TAU_NEW else ''
    mk  = '◀PEAK' if g > peak_gps else ''
    if g > peak_gps:
        peak_gps, peak_tau = g, t
    print(f"  {t:>5}  {phi:.3f}  {g:+.4f}  {src} {mk}")

# τ* 精密値
print(f"\n── τ*精密化 ──")
print(f"  観測ピーク: τ*≈{peak_tau}µs  GPS={peak_gps:.4f}")
phi_star_obs = omega_zz*peak_tau*1e-6
print(f"  φ* = {phi_star_obs:.4f}")

# 理論値との比較
T2_today   = 133.0
tau_theory = T2_today*1e-6/(omega_zz*T2_today*1e-6+1)*1e6
phi_theory = omega_zz*tau_theory*1e-6
print(f"  τ*(理論, T2=133µs) = {tau_theory:.1f}µs  φ*={phi_theory:.4f}")
print(f"  τ*(論文, T2=261µs) = 38.0µs  φ*=0.860")
print(f"  τ*比(実測): {peak_tau/38.0:.2f}×  (T2比: {261/133:.2f}×)")

# Lambert W式との一致確認
print(f"\n  τ* = T2/(ω_ZZ·T2+1)")
print(f"  T2=261µs → τ*={261e-6/(omega_zz*261e-6+1)*1e6:.1f}µs (論文: 38µs)")
print(f"  T2=133µs → τ*={133e-6/(omega_zz*133e-6+1)*1e6:.1f}µs (今日観測: {peak_tau}µs)")

json.dump({
    'anc':ANC,'tgt':TGT,'nu_zz_khz':NU_ZZ,'t2_us':T2_today,
    'peak_tau_us':peak_tau,'peak_gps':peak_gps,
    'phi_star_obs':phi_star_obs,
    'tau_star_theory':tau_theory,
    'new_results':new_res,
    'all_data':all_data
}, open(f"{DIR}/peak_results.json","w"), indent=2)

print(f"\n  保存: {DIR}/  完了")
