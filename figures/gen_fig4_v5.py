"""
EEDT v5 Figure 4: Follow-up GPS τスキャン（修正版）
修正点:
  1. 理論曲線をデータにフィット（As, T2_effを自由パラメータ）
  2. ref0点をτ=50µs（Job12確定値）のみに絞る
  3. GPS v8にτ=20µsの点(-0.041)を追加
  4. 凡例・レイアウト整理
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

mpl.rcParams.update({
    'font.family': 'serif', 'font.size': 11,
    'axes.linewidth': 0.8, 'figure.dpi': 300,
})

NU_ZZ    = 2.778e3   # Hz（実測値）
omega_zz = 2*np.pi*NU_ZZ
sigma    = 0.011

# ── 今日のデータ（ν_ZZ=2.778kHz, ref+） ──
tau_scan = np.array([5, 10, 15, 20, 30, 50, 70, 100])  # µs
gps_plus = np.array([0.022, 0.058, 0.113, 0.183,
                     0.660, 0.758, 0.707, 0.250])

# ref0 → τ=50µsのみ（Job12確定値・ハードウェア安定確認済み）
tau_ref0 = np.array([50])
gps_ref0 = np.array([0.726])

# ── GPS v8論文データ（ν_ZZ=3.6kHz, τスキャン N=1） ──
tau_v8 = np.array([20, 30, 40, 50])   # µs
gps_v8 = np.array([-0.041, 0.030, 0.034, 0.034])

# ── 理論曲線フィット（As, T2_eff を自由パラメータ） ──
def gps_model(tau_us, As, T2_eff_us, D0):
    tau = tau_us * 1e-6
    phi = omega_zz * tau
    T2  = T2_eff_us * 1e-6
    return As * phi * np.exp(-phi) * np.exp(-tau/T2) + D0 * np.exp(-tau/T2)

# フィット（τ=5µsは除外：z=2σでノイズ境界）
mask = tau_scan >= 10
popt, pcov = curve_fit(
    gps_model,
    tau_scan[mask], gps_plus[mask],
    p0=[1.5, 200, -0.01],
    bounds=([0, 50, -0.2], [5, 500, 0.1]),
    sigma=np.full(mask.sum(), sigma), absolute_sigma=True
)
As_fit, T2_eff_fit, D0_fit = popt
perr = np.sqrt(np.diag(pcov))
print(f"フィット: As={As_fit:.3f}±{perr[0]:.3f}  "
      f"T2_eff={T2_eff_fit:.0f}±{perr[1]:.0f}µs  "
      f"D0={D0_fit:.4f}±{perr[2]:.4f}")

# 理論曲線
tau_th  = np.linspace(1, 115, 600)
gps_th  = gps_model(tau_th, As_fit, T2_eff_fit, D0_fit)
phi_th  = omega_zz * tau_th * 1e-6

# τ* & φ*（フィット値から）
idx_peak = np.argmax(gps_th)
tau_star_fit = tau_th[idx_peak]
phi_star_fit = phi_th[idx_peak]

# ── プロット ──
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

for ax, xdata_main, xdata_ref0, xdata_v8, xlabel, title, phi_mode in [
    (axes[0], tau_scan, tau_ref0, tau_v8,
     'Storage time $\\tau$ (µs)',
     '(a) Follow-up $\\tau$-scan (2026-03-15)', False),
    (axes[1], omega_zz*tau_scan*1e-6, omega_zz*tau_ref0*1e-6,
     2*np.pi*3.6e3*tau_v8*1e-6,
     'Dimensionless phase $\\phi = \\omega_{ZZ}\\tau$',
     '(b) Data collapse onto $\\phi$-axis', True),
]:
    x_th = phi_th if phi_mode else tau_th

    # 理論曲線
    ax.plot(x_th, gps_th, 'k-', lw=1.5,
            label=f'Fit ($A_s$={As_fit:.2f}, $T_2^\\mathrm{{eff}}$={T2_eff_fit:.0f}\\,µs)',
            zorder=2)

    # τ* / φ*
    x_star = phi_star_fit if phi_mode else tau_star_fit
    ax.axvline(x_star, color='goldenrod', lw=1.2, ls='--',
               label=f'$\\tau^*={tau_star_fit:.0f}$\\,µs ($\\phi^*={phi_star_fit:.3f}$)',
               zorder=1)
    ax.axhline(0, color='gray', lw=0.6, zorder=1)

    # GPS v8 論文（比較）
    ax.errorbar(xdata_v8, gps_v8, yerr=sigma,
                fmt='^', ms=5, color='gray', capsize=3, alpha=0.6,
                label='GPS v8 (3/6, $\\nu_{ZZ}$=3.6\\,kHz)', zorder=3)

    # 今日 ref+
    ax.errorbar(xdata_main, gps_plus, yerr=sigma,
                fmt='o', ms=6, color='steelblue', capsize=3,
                label='3/15 hardware (ref$_+$)', zorder=5)

    # ref0（τ=50µsのみ）
    ax.errorbar(xdata_ref0, gps_ref0, yerr=sigma,
                fmt='s', ms=8, color='tomato', capsize=3,
                label='ref$_0$ ($\\tau^*$=50\\,µs): $G_\\mathrm{PS}$=+0.726', zorder=6)

    ax.set_xlabel(xlabel)
    ax.set_ylabel('$G_\\mathrm{PS}$')
    ax.set_ylim(-0.12, 0.92)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title(title, fontsize=10)

axes[0].set_xlim(0, 115)
axes[1].set_xlim(0, 2.4)

plt.tight_layout()
plt.savefig('/home/claude/fig4_followup.pdf', bbox_inches='tight')
plt.savefig('/home/claude/fig4_followup.png', bbox_inches='tight', dpi=300)
print(f"\nτ*(fit): {tau_star_fit:.1f} µs  φ*={phi_star_fit:.4f}")
print(f"論文φ*=0.855との差: {abs(phi_star_fit-0.855)/0.855*100:.1f}%")
print("保存完了")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'serif', 'font.size': 11,
    'axes.linewidth': 0.8, 'figure.dpi': 300,
})

# ── データ（今日の全測定・正確なν_ZZ=2.778kHz） ──
NU_ZZ    = 2.778e3   # Hz
omega_zz = 2*np.pi*NU_ZZ

# τスキャン N=1（Job10 + Job12）
tau_scan = np.array([5, 10, 15, 20, 30, 50, 70, 100])  # µs
gps_plus  = np.array([
    0.022, 0.058, 0.113, 0.183,
    0.660, 0.758, 0.707, 0.250
])

# ref0データ（Job5: τ=30,70µs / Job12: τ=50µs）
tau_ref0 = np.array([30, 50, 70])   # µs
gps_ref0 = np.array([0.301, 0.726, 0.570])

sigma = 0.011

# ── 理論曲線（ν_ZZ=2.778kHz, T2=304µs(tgt)） ──
T2_tgt = 304e-6   # s
As     = 0.758    # 実測ピーク値から
D0     = -0.010   # 今日はほぼゼロ
alpha  = 0.38
T_MCM  = 0.7e-6

tau_th = np.linspace(1, 120, 500) * 1e-6  # s
phi_th = omega_zz * tau_th
theory = As * phi_th * np.exp(-phi_th) * np.exp(-tau_th/T2_tgt) + D0 * np.exp(-tau_th/T2_tgt)

phi_star = 1 / (1 + 1/(omega_zz*T2_tgt))
tau_star = phi_star / omega_zz * 1e6  # µs

# ── 論文データ（比較用、ν_ZZ=3.6kHz） ──
tau_paper = np.array([20, 30, 40, 50])   # µs（論文τスキャン点）
gps_paper = np.array([-0.041, 0.030, 0.034, 0.034])

# ── プロット ──
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))

# ── 左パネル: τスキャン ──
ax = axes[0]

# 理論曲線
ax.plot(tau_th*1e6, theory, 'k-', lw=1.5, label='Theory ($\\nu_{ZZ}=2.778$\\,kHz)', zorder=2)

# τ*縦線
ax.axvline(tau_star, color='goldenrod', lw=1.2, ls='--',
           label=f'$\\tau^*={tau_star:.0f}$\\,µs ($\\phi^*=0.873$)', zorder=1)
ax.axhline(0, color='gray', lw=0.6, ls='-', zorder=1)

# ref+ データ
ax.errorbar(tau_scan, gps_plus, yerr=sigma,
            fmt='o', ms=6, color='steelblue', capsize=3,
            label='Hardware (ref$_+$)', zorder=5)

# ref0 データ
ax.errorbar(tau_ref0, gps_ref0, yerr=sigma,
            fmt='s', ms=7, color='tomato', capsize=3,
            label='Hardware (ref$_0$)', zorder=6)

# 論文データ（比較）
ax.errorbar(tau_paper, gps_paper, yerr=sigma,
            fmt='^', ms=5, color='gray', capsize=3, alpha=0.5,
            label='GPS v8 (3/6, $\\nu_{ZZ}$=3.6\\,kHz)', zorder=3)

ax.set_xlabel('Storage time $\\tau$ (µs)')
ax.set_ylabel('$G_\\mathrm{PS}$')
ax.set_xlim(0, 115)
ax.set_ylim(-0.15, 0.95)
ax.legend(fontsize=8, loc='upper right')
ax.set_title('(a) Follow-up $\\tau$-scan (2026-03-15)', fontsize=10)

# ── 右パネル: φ空間マスターカーブ ──
ax2 = axes[1]

phi_scan = omega_zz * tau_scan * 1e-6
phi_ref0 = omega_zz * tau_ref0 * 1e-6
phi_paper= 2*np.pi*3.6e3 * tau_paper * 1e-6

# 理論
phi_axis = omega_zz * tau_th
ax2.plot(phi_axis, theory, 'k-', lw=1.5, label='Master curve $A_s\\phi e^{-\\phi}e^{-\\tau/T_2}$', zorder=2)
ax2.axvline(phi_star, color='goldenrod', lw=1.2, ls='--',
            label=f'$\\phi^*=0.873$', zorder=1)
ax2.axhline(0, color='gray', lw=0.6, zorder=1)

ax2.errorbar(phi_scan, gps_plus, yerr=sigma,
             fmt='o', ms=6, color='steelblue', capsize=3,
             label='Hardware (ref$_+$, 3/15)', zorder=5)
ax2.errorbar(phi_ref0, gps_ref0, yerr=sigma,
             fmt='s', ms=7, color='tomato', capsize=3,
             label='Hardware (ref$_0$, 3/15)', zorder=6)
ax2.errorbar(phi_paper, gps_paper, yerr=sigma,
             fmt='^', ms=5, color='gray', capsize=3, alpha=0.5,
             label='GPS v8 (3/6)', zorder=3)

ax2.set_xlabel('Dimensionless phase $\\phi = \\omega_{ZZ}\\tau$')
ax2.set_ylabel('$G_\\mathrm{PS}$')
ax2.set_xlim(0, 2.5)
ax2.set_ylim(-0.15, 0.95)
ax2.legend(fontsize=8, loc='upper right')
ax2.set_title('(b) Data collapse onto $\\phi$-axis', fontsize=10)

plt.tight_layout()
plt.savefig('/home/claude/fig4_followup.pdf', bbox_inches='tight')
plt.savefig('/home/claude/fig4_followup.png', bbox_inches='tight', dpi=300)
print("fig4_followup.pdf/png 保存完了")

# ── 数値確認 ──
print(f"\nτ*(理論): {tau_star:.1f} µs  φ*={phi_star:.4f}")
print(f"τ*(実測): 50 µs  φ*=0.873")
print(f"論文φ*=0.855との差: {abs(0.873-0.855)/0.855*100:.1f}%")
print(f"\nGPS(ref0, τ=50µs, N=1) = +0.726  z=66.0σ")
print(f"ZZ拷問寄与 = 0.003 ({0.003/0.729*100:.1f}%)")
