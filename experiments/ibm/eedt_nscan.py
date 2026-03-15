"""
EEDT Nスキャン τ=70µs固定
Separation Theorem T2劣化下での成立確認
N=1〜6 / Q94-Q95 / 2026-03-15
"""
import numpy as np, json, os
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

TS  = datetime.now().strftime("%Y%m%d_%H%M%S")
DIR = f"eedt_nscan_{TS}"
os.makedirs(DIR, exist_ok=True)

service = QiskitRuntimeService(
    channel  = "ibm_quantum_platform",
    token    = "YOUR_IBM_QUANTUM_TOKEN",
    instance = "YOUR_IBM_QUANTUM_CRN"
)
backend = service.backend("ibm_marrakesh")
print(f"[{TS}] Connected")

ANC, TGT  = 94, 95
TAU_US    = 70       # τ*ピーク付近
NU_ZZ     = 3.6      # kHz
SHOTS     = 8000
omega_zz  = 2*np.pi*NU_ZZ*1e3
N_LIST    = [1, 2, 3, 4, 5, 6]

print(f"  Q{ANC}-Q{TGT}  τ={TAU_US}µs  N=1〜6")
print(f"  論文N*=5 → T2劣化下で変化するか？")

pm = generate_preset_pass_manager(
    backend=backend, optimization_level=1,
    initial_layout=[ANC, TGT]
)

circs = []
for N in N_LIST:
    # WITH EEDT
    qr = QuantumRegister(2,'q')
    cr = ClassicalRegister(N+1,'c')
    c  = QuantumCircuit(qr, cr)
    c.h(qr[0]); c.h(qr[1])
    dt = int(TAU_US*1000/N)
    for k in range(N):
        c.delay(dt, qr[0], unit='ns')
        c.measure(qr[0], cr[k])
        phi_ff = omega_zz * TAU_US*1e-6 * (k+1)/N
        with c.if_test((cr[k], 1)):
            c.rz(-phi_ff, qr[1])
    c.h(qr[1]); c.measure(qr[1], cr[N])
    circs.append(('with', N, c))

    # REF（N=1と共通で1本でよいが揃えて取る）
    qr2 = QuantumRegister(2,'q')
    cr2 = ClassicalRegister(1,'c')
    c2  = QuantumCircuit(qr2, cr2)
    c2.h(qr2[0]); c2.h(qr2[1])
    c2.delay(int(TAU_US*1000), qr2[1], unit='ns')
    c2.h(qr2[1]); c2.measure(qr2[1], cr2[0])
    circs.append(('ref', N, c2))

isa = [pm.run(c) for _,_,c in circs]
job = SamplerV2(mode=backend).run(isa, shots=SHOTS)
print(f"  Job: {job.job_id()}")
with open(f"{DIR}/jobid.txt","w") as f: f.write(job.job_id())

print("  Waiting...")
res = job.result()

pairs = {}
for i,(kind,N,_) in enumerate(circs):
    cts = res[i].data.c.get_counts()
    if kind == 'with':
        p0 = sum(v for k,v in cts.items() if k[0]=='0') / SHOTS
    else:
        p0 = sum(v for k,v in cts.items() if k[-1]=='0') / SHOTS
    pairs.setdefault(N,{})[kind] = p0

sigma   = 0.011
results = []
print(f"\n  N    GPS      z      判定  (論文比較)")
print(f"  ─  ───────  ─────  ──────────────────")

prev_gps = None
n_star   = None
for N in N_LIST:
    d = pairs.get(N,{})
    if 'with' in d and 'ref' in d:
        gps = d['with'] - d['ref']
        z   = gps / sigma
        # 符号反転検出
        if prev_gps is not None and prev_gps > 0 and gps < 0 and n_star is None:
            n_star = N
        prev_gps = gps
        paper = {1:'+3.35σ',2:'+3.00σ',3:'+2.68σ',4:'+1.45σ',5:'+2.62σ',6:'-2.16σ'}
        ok = '★' if gps>0 and z>2 else ('△' if gps>0 else f'✗ ← 符号反転!')
        print(f"  {N}  {gps:+.4f}  {z:+.1f}  {ok}  (論文:{paper.get(N,'?')})")
        results.append({'N':N,'gps':gps,'z':z,'p_with':d['with'],'p_ref':d['ref']})

print(f"\n── Separation Theorem 解析 ──")
if n_star:
    print(f"  N* = {n_star-1} (N={n_star}で符号反転)")
    if n_star-1 < 5:
        print(f"  → N*が5→{n_star-1}に低下: T2劣化でSeparation Thm境界がシフト")
        print(f"  → 新発見: N* ∝ T2 の可能性")
    elif n_star-1 == 5:
        print(f"  → N*=5維持: Separation ThmはT2劣化に対してロバスト ✓")
    else:
        print(f"  → N*>{n_star-1}: Separation Thm境界が拡大（予想外）")
else:
    if all(r['gps']>0 for r in results):
        print(f"  符号反転なし: N*>6 (T2劣化でN*拡大？)")
    else:
        print(f"  全条件GPS<0: T2={133}µs では全域で動作不可")

# 理論的なN*予測
T_MCM  = 0.7   # µs
alpha  = 0.38
n_star_theory = int(TAU_US / (alpha * T_MCM * 5))
print(f"  理論N*(τ=70µs): N ≤ τ/(5αT_MCM) = {n_star_theory}")

json.dump({
    'anc':ANC,'tgt':TGT,'tau_us':TAU_US,
    'nu_zz_khz':NU_ZZ,'t2_us':133.0,
    'n_star_observed':n_star,
    'n_star_theory':n_star_theory,
    'results':results
}, open(f"{DIR}/nscan_results.json","w"), indent=2)

print(f"\n  保存: {DIR}/  完了")
