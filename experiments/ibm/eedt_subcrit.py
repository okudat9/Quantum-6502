"""
EEDT Subcritical境界スキャン
τ=5, 10, 15µs / N=1
論文: τ=20µsでGPS<0 → 今日は+0.183（逆転）
→ どこで符号が変わるか確定する
Q94-Q95 / 2026-03-15
"""
import numpy as np, json, os
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

TS  = datetime.now().strftime("%Y%m%d_%H%M%S")
DIR = f"eedt_subcrit_{TS}"
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

# 既存データ（今日取得済み）
prev = {
    20: +0.182875,
    30: +0.326000,
    33: +0.371000,
    50: +0.588375,
    70: +0.684125,
   100: +0.249875,
}

# 新規: τ=5, 10, 15µs（subcritical候補）
TAU_NEW = [5, 10, 15]
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

    # REF (ancilla |+>、論文と同定義)
    qr2 = QuantumRegister(2,'q'); cr2 = ClassicalRegister(1,'c')
    c2  = QuantumCircuit(qr2, cr2)
    c2.h(qr2[0]); c2.h(qr2[1])
    c2.delay(int(tau_us*1000), qr2[1], unit='ns')
    c2.h(qr2[1]); c2.measure(qr2[1], cr2[0])
    circs.append(('ref', tau_us, c2))

print(f"  回路数: {len(circs)}  (τ=5,10,15µs × with/ref)")

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
        gps = d['with'] - d['ref']
        new_res[tau_us] = {
            'gps':gps,'z':gps/sigma,
            'p_with':d['with'],'p_ref':d['ref']
        }

# ── 全τ表示（subcritical境界を探す）──
all_gps = {**{t:{'gps':g} for t,g in prev.items()},
           **{t:v for t,v in new_res.items()}}

print(f"\n  τ(µs)   φ=ωτ    GPS      z      判定  論文比較")
print(f"  ─────  ──────  ───────  ─────  ────  ────────────")

sign_cross = None
prev_gps_v = None
prev_tau   = None

for tau_us in sorted(all_gps.keys()):
    d    = all_gps[tau_us]
    gps  = d['gps']
    phi  = omega_zz * tau_us * 1e-6
    z    = d.get('z', gps/sigma)
    src  = '(新)' if tau_us in new_res else '(前)'
    ok   = '★' if gps>0 and abs(z)>2 else ('✗' if gps<0 else '△')

    # 論文の値
    paper = {5:'?',10:'?',15:'?',20:'-3.7σ',30:'+2.7σ',50:'+3.1σ'}
    pm_str = paper.get(tau_us, '')

    print(f"  {tau_us:>5}  {phi:.3f}  {gps:+.4f}  {z:+.1f}  {ok} {src}  {pm_str}")

    # 符号反転検出
    if prev_gps_v is not None:
        if prev_gps_v < 0 and gps > 0:
            sign_cross = (prev_tau, tau_us)
        elif prev_gps_v > 0 and gps < 0:
            sign_cross = (prev_tau, tau_us)
    prev_gps_v = gps
    prev_tau   = tau_us

# ── 結論 ──
print(f"\n── Subcritical境界解析 ──")
print(f"  論文: τ=20µs(φ=0.45) → GPS=-0.041  subcritical確認")
print(f"  今日: τ=20µs(φ=0.45) → GPS=+0.183  ← 逆転！")

if sign_cross:
    t1,t2 = sign_cross
    print(f"\n  符号反転: τ∈[{t1},{t2}]µs")
    phi1 = omega_zz*t1*1e-6
    phi2 = omega_zz*t2*1e-6
    print(f"  φ境界: [{phi1:.3f},{phi2:.3f}]")
    print(f"  → T2劣化でsubcritical境界がφ<{phi1:.2f}側にシフト")
    print(f"  → D0（ベースライン）がT2劣化で正にシフトした証拠")
else:
    min_gps = min(all_gps[t]['gps'] for t in TAU_NEW)
    if min_gps > 0:
        print(f"\n  τ=5µsでもGPS>0 → subcritical境界はτ<5µs")
        print(f"  → T2劣化でD0が完全に正にシフト")
        print(f"  → 論文のsubcritical regime（GPS<0）が消失")
    print(f"  (符号反転なし)")

# φ* 論文との比較
tau_star_paper = 38.0
tau_star_today = 70.0
phi_star_paper = omega_zz * tau_star_paper * 1e-6
phi_star_today = omega_zz * tau_star_today * 1e-6
print(f"\n  φ*(論文): {phi_star_paper:.3f}  τ*=38µs  T2=261µs")
print(f"  φ*(今日): {phi_star_today:.3f}  τ*=70µs  T2=133µs")
print(f"  → φ*がほぼ2倍: T2半減でτ*は約2倍にシフト")

json.dump({
    'anc':ANC,'tgt':TGT,'nu_zz_khz':NU_ZZ,'t2_us':133.0,
    'sign_cross_tau':sign_cross,
    'new_results':new_res,
    'prev_results':prev,
    'conclusion':{
        'phi_star_paper':phi_star_paper,
        'phi_star_today':phi_star_today,
        'tau_star_shift_us':[tau_star_paper, tau_star_today]
    }
}, open(f"{DIR}/subcrit_results.json","w"), indent=2)

print(f"\n  保存: {DIR}/  完了")
