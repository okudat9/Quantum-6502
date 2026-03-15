"""
EEDT τ延長スキャン Q94-Q95
τ=70, 100µs / N=1のみ / 4回路
山型ピーク確認用
"""
import numpy as np, json, os
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

TS  = datetime.now().strftime("%Y%m%d_%H%M%S")
DIR = f"eedt_tauscan_{TS}"
os.makedirs(DIR, exist_ok=True)

service = QiskitRuntimeService(
    channel  = "ibm_quantum_platform",
    token    = "YOUR_IBM_QUANTUM_TOKEN",
    instance = "YOUR_IBM_QUANTUM_CRN"
)
backend = service.backend("ibm_marrakesh")
print(f"[{TS}] Connected")

ANC, TGT  = 94, 95
NU_ZZ     = 3.6    # kHz
SHOTS     = 8000
omega_zz  = 2*np.pi*NU_ZZ*1e3

pm = generate_preset_pass_manager(
    backend=backend, optimization_level=1,
    initial_layout=[ANC, TGT]
)

# 今日の既存データ（前のジョブ）
prev = {
    20: {'gps':0.182875, 'z':16.6},
    30: {'gps':0.326000, 'z':29.6},
    33: {'gps':0.371000, 'z':33.7},
    50: {'gps':0.588375, 'z':53.5},
}

# 新規スキャン: τ=70, 100µs / N=1のみ
TAU_NEW = [70, 100]
N = 1
circs = []

for tau_us in TAU_NEW:
    # WITH
    qr = QuantumRegister(2,'q'); cr = ClassicalRegister(2,'c')
    c  = QuantumCircuit(qr, cr)
    c.h(qr[0]); c.h(qr[1])
    c.delay(int(tau_us*1000), qr[0], unit='ns')
    c.measure(qr[0], cr[0])
    phi_ff = omega_zz * tau_us*1e-6
    with c.if_test((cr[0], 1)):
        c.rz(-phi_ff, qr[1])
    c.h(qr[1]); c.measure(qr[1], cr[1])
    circs.append(('with', tau_us, c))

    # REF (ancilla NOT in |+>, just delay)
    qr2 = QuantumRegister(2,'q'); cr2 = ClassicalRegister(1,'c')
    c2  = QuantumCircuit(qr2, cr2)
    c2.h(qr2[0]); c2.h(qr2[1])
    c2.delay(int(tau_us*1000), qr2[1], unit='ns')
    c2.h(qr2[1]); c2.measure(qr2[1], cr2[0])
    circs.append(('ref', tau_us, c2))

isa   = [pm.run(c) for _,_,c in circs]
job   = SamplerV2(mode=backend).run(isa, shots=SHOTS)
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

sigma = 0.011
new_results = {}
for tau_us in TAU_NEW:
    d = pairs.get(tau_us,{})
    if 'with' in d and 'ref' in d:
        gps = d['with'] - d['ref']
        z   = gps / sigma
        new_results[tau_us] = {'gps':gps,'z':z,
                                'p_with':d['with'],'p_ref':d['ref']}

# ── 全データ表示 ──
print(f"\n  τ(µs)   GPS      z      判定")
print(f"  ─────  ───────  ─────  ──────")
all_tau = sorted(list(prev.keys()) + TAU_NEW)
peak_gps, peak_tau = 0, 0
for t in all_tau:
    if t in prev:
        g,z = prev[t]['gps'], prev[t]['z']
        src = '(前回)'
    elif t in new_results:
        g,z = new_results[t]['gps'], new_results[t]['z']
        src = '(今回)'
    else:
        continue
    mark = '★PEAK?' if g > peak_gps else ''
    if g > peak_gps:
        peak_gps, peak_tau = g, t
    ok = '★' if g>0 and z>2 else '✗'
    print(f"  {t:>5}   {g:+.4f}  {z:+.1f}  {ok} {src} {mark}")

# ── 山型判定 ──
print(f"\n── 山型解析 ──")
if TAU_NEW[0] in new_results and TAU_NEW[1] in new_results:
    g70  = new_results[70]['gps']
    g100 = new_results[100]['gps']
    g50  = prev[50]['gps']

    if g70 < g50 and g100 < g70:
        print(f"  → 山型確認 ✓  ピーク: τ∈[50,70]µs")
        print(f"  → τ*が33µsより大きい（T2劣化でτ*シフト）")
        shape = 'peaked'
    elif g70 > g50:
        print(f"  → τ=70でも上昇中 → ピーク未到達")
        print(f"  → τ*>70µs の可能性（T2劣化と矛盾？）")
        shape = 'rising'
    else:
        print(f"  → 平坦またはノイズ範囲")
        shape = 'flat'

    # τ*(T2劣化版)理論値
    T2_new = 133.0
    tau_star_theory = T2_new*1e-6/(omega_zz*T2_new*1e-6+1)*1e6
    print(f"\n  T2={T2_new}µs時のτ*理論値: {tau_star_theory:.1f}µs")
    print(f"  実測ピーク付近: τ≈{peak_tau}µs (GPS={peak_gps:.3f})")

    json.dump({
        'anc':ANC,'tgt':TGT,'nu_zz_khz':NU_ZZ,'t2_us':T2_new,
        'tau_star_theory':tau_star_theory,
        'shape':shape,'peak_tau':peak_tau,'peak_gps':peak_gps,
        'new':new_results,'prev':prev
    }, open(f"{DIR}/tauscan_results.json","w"), indent=2)

print(f"\n  保存: {DIR}/  完了")
