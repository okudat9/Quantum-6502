"""
EEDT GPS 正しいnu_ZZ=2.778kHzでの再測定
τ=70µs / N=1,2,3 / Q94-Q95
今日の全GPSは3.6kHz仮定（30%過補正）だった
正しいフィードフォワード角での真値取得
"""
import numpy as np, json, os
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

TS  = datetime.now().strftime("%Y%m%d_%H%M%S")
DIR = f"eedt_corrected_{TS}"
os.makedirs(DIR, exist_ok=True)

service = QiskitRuntimeService(
    channel  = "ibm_quantum_platform",
    token    = "YOUR_IBM_QUANTUM_TOKEN",
    instance = "YOUR_IBM_QUANTUM_CRN"
)
backend = service.backend("ibm_marrakesh")
print(f"[{TS}] Connected")

ANC, TGT  = 94, 95
TAU_US    = 70
SHOTS     = 8000

NU_ZZ_OLD = 3.600   # 今日使っていた（論文値）
NU_ZZ_NEW = 2.778   # 実測値（Job8 Ramsey結果）

omega_old = 2*np.pi*NU_ZZ_OLD*1e3
omega_new = 2*np.pi*NU_ZZ_NEW*1e3

print(f"\n  τ={TAU_US}µs  Q{ANC}-Q{TGT}")
print(f"  ν_ZZ(旧): {NU_ZZ_OLD} kHz → φ_ff={omega_old*TAU_US*1e-6:.4f}rad")
print(f"  ν_ZZ(新): {NU_ZZ_NEW} kHz → φ_ff={omega_new*TAU_US*1e-6:.4f}rad")
print(f"  過補正率: {NU_ZZ_OLD/NU_ZZ_NEW:.3f}×")

pm = generate_preset_pass_manager(
    backend=backend, optimization_level=1,
    initial_layout=[ANC, TGT]
)

N_LIST = [1, 2, 3]
circs  = []

for N in N_LIST:
    # WITH EEDT (正しいν_ZZ)
    qr = QuantumRegister(2,'q'); cr = ClassicalRegister(N+1,'c')
    c  = QuantumCircuit(qr, cr)
    c.h(qr[0]); c.h(qr[1])
    dt = int(TAU_US*1000/N)
    for k in range(N):
        c.delay(dt, qr[0], unit='ns')
        c.measure(qr[0], cr[k])
        phi_ff = omega_new * TAU_US*1e-6 * (k+1)/N  # ← 正しい角度
        with c.if_test((cr[k],1)):
            c.rz(-phi_ff, qr[1])
    c.h(qr[1]); c.measure(qr[1], cr[N])
    circs.append(('with_new', N, c))

    # REF
    qr2 = QuantumRegister(2,'q'); cr2 = ClassicalRegister(1,'c')
    c2  = QuantumCircuit(qr2, cr2)
    c2.h(qr2[0]); c2.h(qr2[1])
    c2.delay(int(TAU_US*1000), qr2[1], unit='ns')
    c2.h(qr2[1]); c2.measure(qr2[1], cr2[0])
    circs.append(('ref', N, c2))

print(f"\n  {len(circs)}回路送信中...")
isa = [pm.run(c) for _,_,c in circs]
job = SamplerV2(mode=backend).run(isa, shots=SHOTS)
print(f"  Job: {job.job_id()}")
with open(f"{DIR}/jobid.txt","w") as f: f.write(job.job_id())

print("  Waiting...")
res = job.result()

pairs = {}
for i,(kind,N,_) in enumerate(circs):
    cts = res[i].data.c.get_counts()
    if kind == 'with_new':
        p0 = sum(v for k,v in cts.items() if k[0]=='0') / SHOTS
    else:
        p0 = sum(v for k,v in cts.items() if k[-1]=='0') / SHOTS
    pairs.setdefault(N,{})[kind] = p0

# 旧値（過補正）との比較
old_gps = {1:0.684125, 2:0.553125, 3:0.336000}

sigma   = 0.011
results = []
print(f"\n  N   GPS(旧φ)  GPS(正φ)  差      改善?")
print(f"  ─  ─────────  ─────────  ──────  ─────")
for N in N_LIST:
    d = pairs.get(N,{})
    if 'with_new' in d and 'ref' in d:
        gps_new = d['with_new'] - d['ref']
        gps_old = old_gps.get(N, float('nan'))
        diff    = gps_new - gps_old
        better  = '↑改善' if gps_new > gps_old else ('↓低下' if gps_new < gps_old else '変化なし')
        print(f"  {N}   {gps_old:+.4f}    {gps_new:+.4f}    {diff:+.4f}  {better}")
        results.append({
            'N':N,'gps_old_phi':gps_old,'gps_new_phi':gps_new,
            'diff':diff,'p_with':d['with_new'],'p_ref':d['ref']
        })

# 最適φの解析
print(f"\n── フィードフォワード角の最適化 ──")
print(f"  φ(旧) = ω×3.6kHz × τ = {omega_old*TAU_US*1e-6:.4f} rad")
print(f"  φ(新) = ω×2.778kHz × τ = {omega_new*TAU_US*1e-6:.4f} rad")
print(f"  φ*(論文) = 0.855 rad (ϕ*=1/A)")
print(f"  φ(新)に対するφ* = {omega_new*TAU_US*1e-6:.4f} → {'過補正' if omega_new*TAU_US*1e-6 > 1.0 else '適正'}")

json.dump({
    'anc':ANC,'tgt':TGT,'tau_us':TAU_US,
    'nu_zz_old':NU_ZZ_OLD,'nu_zz_new':NU_ZZ_NEW,
    'overcorrection_factor':NU_ZZ_OLD/NU_ZZ_NEW,
    'results':results
}, open(f"{DIR}/corrected_results.json","w"), indent=2)

print(f"\n  保存: {DIR}/  完了")
