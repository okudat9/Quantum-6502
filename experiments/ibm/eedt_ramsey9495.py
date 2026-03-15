"""
ZZ Ramsey Q94-Q95 ν_ZZ現在値確認
論文値3.6kHz(9日前)から変化しているか？
τ*ズレの原因解明 / 2026-03-15
"""
import numpy as np, json, os
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

TS  = datetime.now().strftime("%Y%m%d_%H%M%S")
DIR = f"eedt_ramsey9495_{TS}"
os.makedirs(DIR, exist_ok=True)

service = QiskitRuntimeService(
    channel  = "ibm_quantum_platform",
    token    = "YOUR_IBM_QUANTUM_TOKEN",
    instance = "YOUR_IBM_QUANTUM_CRN"
)
backend = service.backend("ibm_marrakesh")
print(f"[{TS}] Connected")

ANC, TGT = 94, 95
SHOTS    = 8000

pm = generate_preset_pass_manager(
    backend=backend, optimization_level=1,
    initial_layout=[ANC, TGT]
)

# 論文τ_CR=31.3µsを基準に符号反転が見えるτ点を選ぶ
# ν_ZZ=3.6kHz → 周期=277µs → 1/4=69µs付近で反転
# ν_ZZが変化していれば反転点がずれる
TAU_R = [10, 30, 50, 70, 100, 140, 180]
circs = []

for tau_us in TAU_R:
    for anc_st in [0, 1]:
        qr = QuantumRegister(2,'q'); cr = ClassicalRegister(1,'c')
        c  = QuantumCircuit(qr, cr)
        if anc_st == 1: c.x(qr[0])
        c.h(qr[1])
        c.delay(int(tau_us*1000), qr[1], unit='ns')
        c.h(qr[1])
        c.measure(qr[1], cr[0])
        circs.append((tau_us, anc_st, c))

print(f"  {len(circs)}回路 (τ=7点×anc|0>/|1>)")

isa = [pm.run(c) for _,_,c in circs]
job = SamplerV2(mode=backend).run(isa, shots=SHOTS)
print(f"  Job: {job.job_id()}")
with open(f"{DIR}/jobid.txt","w") as f: f.write(job.job_id())

print("  Waiting...")
res = job.result()

freqs = {}
for i,(tau_us,anc_st,_) in enumerate(circs):
    cts = res[i].data.c.get_counts()
    p0  = cts.get('0',0)/SHOTS
    freqs.setdefault(tau_us,{})[anc_st] = p0

corr = {t: d[0]-d[1] for t,d in freqs.items() if 0 in d and 1 in d}

print(f"\n  τ(µs)   C(τ)=P0(0)-P0(1)   バー")
print(f"  ─────  ─────────────────  ─────────────────────")
for t in sorted(corr):
    c_val = corr[t]
    bar   = '█'*int(abs(c_val)*40)
    sgn   = '+' if c_val >= 0 else '-'
    print(f"  {t:>5}  {c_val:+.4f}              {bar}")

# ν_ZZ推定（符号反転から）
tau_list  = sorted(corr.keys())
nu_zz_est = None
for i in range(len(tau_list)-1):
    t1,t2 = tau_list[i], tau_list[i+1]
    if corr[t1]*corr[t2] < 0:
        tc    = t1+(t2-t1)*abs(corr[t1])/(abs(corr[t1])+abs(corr[t2]))
        nu_zz_est = 1000/(4*tc)
        print(f"\n  ★ 符号反転 τ≈{tc:.1f}µs → ν_ZZ≈{nu_zz_est:.3f} kHz")
        break

if nu_zz_est is None:
    print(f"\n  符号反転なし → ν_ZZ<1.8kHz or >5kHz")
    nu_zz_est = None

# 論文値との比較
print(f"\n── ν_ZZ比較 ──")
print(f"  論文値(3/6):  3.600 kHz  T2=261µs")
if nu_zz_est:
    ratio = nu_zz_est/3.6
    print(f"  今日の測定:  {nu_zz_est:.3f} kHz  T2=133µs")
    print(f"  変化率: {ratio:.2f}×")

    # τ*再計算
    omega_new  = 2*np.pi*nu_zz_est*1e3
    T2_today   = 133.0
    tau_new    = T2_today*1e-6/(omega_new*T2_today*1e-6+1)*1e6
    phi_new    = omega_new*tau_new*1e-6
    print(f"\n  ν_ZZ={nu_zz_est:.3f}kHz, T2=133µs での τ*(理論)={tau_new:.1f}µs")
    print(f"  φ* = {phi_new:.4f}")
    print(f"  実測ピーク: τ*≈70µs")
    print(f"  → 理論と実測の一致度: {abs(tau_new-70)/70*100:.0f}%ズレ")
else:
    print(f"  今日の測定: 符号反転なし（測定範囲外）")

json.dump({
    'anc':ANC,'tgt':TGT,
    'corr':{str(k):v for k,v in corr.items()},
    'nu_zz_khz_paper':3.6,
    'nu_zz_khz_today':nu_zz_est,
    't2_today_us':133.0
}, open(f"{DIR}/ramsey_results.json","w"), indent=2)

print(f"\n  保存: {DIR}/  完了")
