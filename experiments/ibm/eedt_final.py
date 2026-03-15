"""
EEDT GPS ref0確定測定
τ=50µs(τ*) / ν_ZZ=2.778kHz(実測) / N=1,2,3
ancilla|0>基準での真のEEDT利得確定
昼食後・ハードウェア安定確認込み
"""
import numpy as np, json, os
from datetime import datetime
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

TS  = datetime.now().strftime("%Y%m%d_%H%M%S")
DIR = f"eedt_final_{TS}"
os.makedirs(DIR, exist_ok=True)

service = QiskitRuntimeService(
    channel  = "ibm_quantum_platform",
    token    = "YOUR_IBM_QUANTUM_TOKEN",
    instance = "YOUR_IBM_QUANTUM_CRN"
)
backend = service.backend("ibm_marrakesh")
print(f"[{TS}] Connected: ibm_marrakesh")

ANC, TGT = 94, 95
NU_ZZ    = 2.778   # kHz 実測値（Job8 Ramsey）
TAU_US   = 50      # τ* 確定値
SHOTS    = 8000
omega_zz = 2*np.pi*NU_ZZ*1e3

phi_star = omega_zz * TAU_US*1e-6
print(f"\n  Q{ANC}-Q{TGT}  ν_ZZ={NU_ZZ}kHz  τ={TAU_US}µs")
print(f"  φ = {phi_star:.4f} rad  (論文φ*=0.855)")

pm = generate_preset_pass_manager(
    backend=backend, optimization_level=1,
    initial_layout=[ANC, TGT]
)

N_LIST = [1, 2, 3]
circs  = []

for N in N_LIST:
    # ── WITH EEDT ──
    qr = QuantumRegister(2,'q'); cr = ClassicalRegister(N+1,'c')
    c  = QuantumCircuit(qr, cr)
    c.h(qr[0]); c.h(qr[1])
    dt = int(TAU_US*1000/N)
    for k in range(N):
        c.delay(dt, qr[0], unit='ns')
        c.measure(qr[0], cr[k])
        phi_ff = omega_zz * TAU_US*1e-6 * (k+1)/N
        with c.if_test((cr[k],1)):
            c.rz(-phi_ff, qr[1])
    c.h(qr[1]); c.measure(qr[1], cr[N])
    circs.append(('with', N, c))

    # ── REF_plus (ancilla |+>, 論文定義) ──
    qr2 = QuantumRegister(2,'q'); cr2 = ClassicalRegister(1,'c')
    c2  = QuantumCircuit(qr2, cr2)
    c2.h(qr2[0]); c2.h(qr2[1])
    c2.delay(int(TAU_US*1000), qr2[1], unit='ns')
    c2.h(qr2[1]); c2.measure(qr2[1], cr2[0])
    circs.append(('ref_plus', N, c2))

    # ── REF_zero (ancilla |0>, 真のベースライン) ──
    qr3 = QuantumRegister(2,'q'); cr3 = ClassicalRegister(1,'c')
    c3  = QuantumCircuit(qr3, cr3)
    # ancilla は |0>のまま（hなし）
    c3.h(qr3[1])
    c3.delay(int(TAU_US*1000), qr3[1], unit='ns')
    c3.h(qr3[1]); c3.measure(qr3[1], cr3[0])
    circs.append(('ref_zero', N, c3))

print(f"\n  {len(circs)}回路 (N=1,2,3 × with/ref+/ref0)")

isa = [pm.run(c) for _,_,c in circs]
job = SamplerV2(mode=backend).run(isa, shots=SHOTS)
print(f"  Job: {job.job_id()}")
with open(f"{DIR}/jobid.txt","w") as f: f.write(job.job_id())

print("  Waiting...")
res = job.result()

# ── 結果取得 ──
data = {}
for i,(kind,N,_) in enumerate(circs):
    cts = res[i].data.c.get_counts()
    if kind == 'with':
        p0 = sum(v for k,v in cts.items() if k[0]=='0') / SHOTS
    else:
        p0 = sum(v for k,v in cts.items() if k[-1]=='0') / SHOTS
    data.setdefault(N,{})[kind] = p0

sigma   = 0.011
results = []

print(f"\n  N   p_with  p_ref+  p_ref0  GPS(+)   GPS(0)   z(0)")
print(f"  ─  ──────  ──────  ──────  ───────  ───────  ─────")

for N in N_LIST:
    d   = data[N]
    pw  = d.get('with',     float('nan'))
    rp  = d.get('ref_plus', float('nan'))
    rz  = d.get('ref_zero', float('nan'))
    gp  = pw - rp
    gz  = pw - rz
    z_gz = gz / sigma
    print(f"  {N}  {pw:.4f}  {rp:.4f}  {rz:.4f}  "
          f"{gp:+.4f}   {gz:+.4f}   {z_gz:+.1f}")
    results.append({
        'N':N,'p_with':pw,'p_ref_plus':rp,'p_ref_zero':rz,
        'gps_plus':gp,'gps_zero':gz,'z_zero':z_gz
    })

# ── ハードウェア安定確認 ──
print(f"\n── ハードウェア安定確認 ──")
p_ref_prev = 0.158   # Job10のp_ref(τ=50µs)
p_ref_now  = data.get(1,{}).get('ref_plus', float('nan'))
drift      = abs(p_ref_now - p_ref_prev)
stable     = drift < 0.05
print(f"  p_ref+(前回Job10): {p_ref_prev:.4f}")
print(f"  p_ref+(今回):      {p_ref_now:.4f}")
print(f"  ドリフト: {drift:.4f}  {'✓ 安定' if stable else '⚠ ドリフト継続'}")

# ── 最終サマリ ──
gz_n1 = next((r['gps_zero'] for r in results if r['N']==1), None)
gp_n1 = next((r['gps_plus'] for r in results if r['N']==1), None)

print(f"\n── 最終確定値（τ*=50µs, ν_ZZ=2.778kHz） ──")
if gz_n1 is not None:
    print(f"  GPS(ref+, N=1) = {gp_n1:+.4f}  z={gp_n1/sigma:+.1f}σ")
    print(f"  GPS(ref0, N=1) = {gz_n1:+.4f}  z={gz_n1/sigma:+.1f}σ  ← 真のEEDT利得")
    print(f"  ZZ拷問寄与    = {gp_n1-gz_n1:+.4f} ({(gp_n1-gz_n1)/gp_n1*100:.0f}%)")
    if gz_n1 > 0:
        print(f"\n  ★ 真のEEDT利得確定: GPS(ref0)>0 ✓")
        print(f"  φ={phi_star:.3f}≈φ*=0.855 で論文理論と一致")
    else:
        print(f"\n  ✗ GPS(ref0)<0: T2劣化でEEDT効果消失")

json.dump({
    'backend':'ibm_marrakesh','anc':ANC,'tgt':TGT,
    'nu_zz_khz':NU_ZZ,'tau_us':TAU_US,'phi':phi_star,
    'hardware_stable':stable,'p_ref_drift':float(drift),
    'results':results,
    'note':'Final measurement: true EEDT gain with correct nu_ZZ and ref0'
}, open(f"{DIR}/final_results.json","w"), indent=2)

print(f"\n  保存: {DIR}/  完了")
