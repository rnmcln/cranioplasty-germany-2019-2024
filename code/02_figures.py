"""
Figures from output/analysis_results.json and data/input/drg_cost_b20.csv.

Main figures
    figure1_volumes          annual discharge episodes by material group; autologous share
    figure2_mortality        in-hospital mortality by year and pooled, Wilson 95% CI
    figure3_infection        coded infection (T81.4, T85.72), pooled, Wilson 95% CI
Supplementary figures
    figureS1_drg_cost        whole-DRG per-case cost, B20 family, 2022-2024
    figureS2_case_mix        length of stay, PCCL 0, male and long-stay proportions by year

Each figure is written as vector PDF and 600 dpi PNG to figures/.
Run from the repository root:  python code/02_figures.py
"""
import csv
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

MM = 1 / 25.4
avail = {f.name for f in font_manager.fontManager.ttflist}
FONT = next((f for f in ('Arial', 'Helvetica', 'Liberation Sans', 'Nimbus Sans', 'DejaVu Sans') if f in avail), 'sans-serif')
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': [FONT], 'font.size': 7,
    'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.direction': 'out', 'ytick.direction': 'out',
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#3B4252',
    'xtick.color': '#3B4252', 'ytick.color': '#3B4252', 'axes.labelcolor': '#2E3440', 'text.color': '#2E3440',
    'lines.linewidth': 0.7, 'lines.markersize': 3.2, 'legend.frameon': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.dpi': 600, 'figure.dpi': 150, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
})

YEARS = [2019, 2020, 2021, 2022, 2023, 2024]
G = ['Autologous', 'Simple alloplastic', 'CAD']
LABEL = {'Autologous': 'Autologous', 'Simple alloplastic': 'Simple alloplastic', 'CAD': 'CAD'}
COL = {'Autologous': '#2E3440', 'Simple alloplastic': '#8C939E', 'CAD': '#5E81AC'}   # greys + one muted accent
MK = {'Autologous': 'o', 'Simple alloplastic': 's', 'CAD': '^'}
FILL = {'Autologous': '#2E3440', 'Simple alloplastic': 'white', 'CAD': 'white'}
OFF = {'Autologous': -0.17, 'Simple alloplastic': 0.0, 'CAD': 0.17}

res = json.load(open('output/analysis_results.json'))
os.makedirs('figures', exist_ok=True)


def save(fig, name):
    fig.savefig(f'figures/{name}.pdf', metadata={'CreationDate': None, 'Producer': None})
    fig.savefig(f'figures/{name}.png')
    plt.close(fig)


def panel(ax, letter):
    ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontsize=8, va='bottom', ha='left')


def year_axis(ax):
    ax.set_xticks(YEARS)
    ax.set_xlim(2018.6, 2024.4)


def spread(values, gap):
    """Vertical positions for direct labels so that neighbouring labels are at least `gap` apart."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    pos = [values[i] for i in order]
    for _ in range(50):
        moved = False
        for i in range(1, len(pos)):
            if pos[i] - pos[i - 1] < gap:
                shift = (gap - (pos[i] - pos[i - 1])) / 2
                pos[i - 1] -= shift
                pos[i] += shift
                moved = True
        if not moved:
            break
    out = [0.0] * len(values)
    for k, i in enumerate(order):
        out[i] = pos[k]
    return out


def end_labels(ax, ends, gap, x=2024.25, size=6.3):
    ys = spread([ends[g] for g in G], gap)
    for g, yl in zip(G, ys):
        ax.text(x, yl, LABEL[g], color=COL[g], va='center', ha='left', fontsize=size)


def point(ax, x, y, g, **kw):
    ax.plot(x, y, marker=MK[g], color=COL[g], mfc=FILL[g], mec=COL[g], mew=0.6, ls='none', zorder=3, **kw)


# ------------------------------------------------------------------ Figure 1
fig, axs = plt.subplots(1, 2, figsize=(180 * MM, 62 * MM), gridspec_kw={'width_ratios': [1.35, 1], 'wspace': 0.42})
ax = axs[0]
for g in G:
    v = [res['volumes'][str(y)][g] for y in YEARS]
    ax.plot(YEARS, v, color=COL[g], lw=0.7, zorder=2)
    point(ax, YEARS, v, g)
end_labels(ax, {g: res['volumes']['2024'][g] for g in G}, 250)
year_axis(ax)
ax.set_ylim(0, 3600)
ax.set_yticks(range(0, 3601, 1000))
ax.set_yticklabels([f'{t:,}' for t in range(0, 3601, 1000)])
ax.set_ylabel('Discharge episodes (n)')
panel(ax, 'A')

ax = axs[1]
s = [res['autologous_share_pct'][str(y)] for y in YEARS]
ax.plot(YEARS, s, color=COL['Autologous'], lw=0.7)
point(ax, YEARS, s, 'Autologous')
for y, v in zip(YEARS, s):
    ax.text(y, v + 1.4, f'{v:.1f}', ha='center', va='bottom', fontsize=6, color='#4C566A')
year_axis(ax)
ax.set_ylim(0, 50)
ax.set_ylabel('Autologous share of cranial\nreconstruction (%)')
panel(ax, 'B')
save(fig, 'figure1_volumes')

# ------------------------------------------------------------------ Figure 2
fig, axs = plt.subplots(1, 2, figsize=(180 * MM, 62 * MM), gridspec_kw={'width_ratios': [1.45, 1], 'wspace': 0.6})
ax = axs[0]
for g in G:
    x = [y + OFF[g] for y in YEARS]
    est = [res['mortality_by_year'][str(y)][g]['pct_ci'] for y in YEARS]
    for xi, (p, lo, hi) in zip(x, est):
        ax.plot([xi, xi], [lo, hi], color=COL[g], lw=0.6, zorder=2)
    point(ax, x, [e[0] for e in est], g)
    ax.text(2024.42, {'Autologous': 3.45, 'Simple alloplastic': 2.55, 'CAD': 1.2}[g], LABEL[g], color=COL[g], va='center', ha='left', fontsize=6.3)
ax.set_xticks(YEARS)
ax.set_xlim(2018.6, 2024.4)
ax.set_ylim(0, 4.6)
ax.set_ylabel('In-hospital mortality (%)')
panel(ax, 'A')

ax = axs[1]
ys = {g: i for i, g in enumerate(reversed(G))}
for g in G:
    m = res['mortality_pooled'][g]
    p, lo, hi = m['pct_ci']
    ax.plot([lo, hi], [ys[g]] * 2, color=COL[g], lw=0.7)
    point(ax, [p], [ys[g]], g)
    ax.text(4.5, ys[g], f'{m["deaths"]:,}/{m["N"]:,}', va='center', ha='right', fontsize=6, color='#4C566A')
ax.set_yticks(list(ys.values()))
ax.set_yticklabels([LABEL[g] for g in ys])
ax.set_ylim(-0.6, 2.6)
ax.set_xlim(0, 4.5)
ax.set_xlabel('In-hospital mortality, 2019–2024 (%)')
ax.spines['left'].set_visible(False)
ax.tick_params(axis='y', length=0)
panel(ax, 'B')
save(fig, 'figure2_mortality')

# ------------------------------------------------------------------ Figure 3
fig, axs = plt.subplots(1, 2, figsize=(150 * MM, 46 * MM), sharey=True, gridspec_kw={'wspace': 0.12})
spec = [('T81.4', 'Post-procedural infection, T81.4 (%)', 5.0), ('T85.72', 'Nervous-system implant infection, T85.72 (%)', 1.5)]
ys = {g: i for i, g in enumerate(reversed(G))}
for ax, (code, xl, xmax), letter in zip(axs, spec, 'AB'):
    for g in G:
        c = res['infection_pooled'][code][g]
        p, lo, hi = c['pct_ci']
        ax.plot([lo, hi], [ys[g]] * 2, color=COL[g], lw=0.7)
        point(ax, [p], [ys[g]], g)
        ax.text(xmax, ys[g] + 0.22, f'{c["k"]}/{c["N"]:,}', va='bottom', ha='right', fontsize=6, color='#4C566A')
    ax.set_xlim(0, xmax)
    if xmax < 2:
        ax.set_xticks([0, 0.5, 1.0, 1.5])
        ax.set_xticklabels(['0', '0.5', '1.0', '1.5'])
    ax.set_xlabel(xl)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    panel(ax, letter)
axs[0].set_yticks(list(ys.values()))
axs[0].set_yticklabels([LABEL[g] for g in ys])
axs[0].set_ylim(-0.6, 2.75)
save(fig, 'figure3_infection')

# ------------------------------------------------------------------ Figure S1
cost = list(csv.DictReader(open('data/input/drg_cost_b20.csv', newline='')))
cats = ['B20A', 'B20B', 'B20C', 'B20D', 'B20E', 'B20 family']
ticklab = ['B20A', 'B20B', 'B20C', 'B20D', 'B20E', 'B20 family\n(case-weighted)']
cyears = [2022, 2023, 2024]
shade = {2022: '#B0B6BF', 2023: '#6B7380', 2024: '#2E3440'}
cmk = {2022: 'o', 2023: 's', 2024: '^'}
fig, ax = plt.subplots(figsize=(120 * MM, 55 * MM))
for j, y in enumerate(cyears):
    rows = {r['drg']: r for r in cost if int(r['year']) == y}
    xs = [i + (j - 1) * 0.2 for i in range(len(cats))]
    vals = [float(rows[c]['mean_cost_eur']) for c in cats]
    ax.plot(xs, vals, ls='none', marker=cmk[y], color=shade[y], mfc='white' if y != 2024 else shade[y], mew=0.6)
    ax.text(xs[0], vals[0] + 900, str(y), rotation=90, fontsize=5.5, ha='center', va='bottom', color=shade[y])
ax.axvline(len(cats) - 1.5, color='#D8DEE9', lw=0.5)
ax.set_xticks(range(len(cats)))
ax.set_xticklabels(ticklab)
ax.set_xlim(-0.5, len(cats) - 0.5)
ax.set_ylim(0, 18000)
ax.set_yticks(range(0, 18001, 4000))
ax.set_yticklabels([f'{t:,}' for t in range(0, 18001, 4000)])
ax.set_ylabel('Mean per-case cost (EUR)')
save(fig, 'figureS1_drg_cost')

# ------------------------------------------------------------------ Figure S2
fig, axs = plt.subplots(2, 2, figsize=(170 * MM, 105 * MM), gridspec_kw={'wspace': 0.45, 'hspace': 0.45})
spec = [('mean_los', 'Mean length of stay (days)', (0, 25)), ('pccl0_pct', 'PCCL 0 (%)', (0, 70)),
        ('male_pct', 'Male (%)', (0, 70)), ('lang_pct', 'Long-stay episodes (%)', (0, 20))]
for ax, (f, yl, lim), letter in zip(axs.flat, spec, 'ABCD'):
    for g in G:
        v = [res['case_mix'][str(y)][g][f] for y in YEARS]
        ax.plot(YEARS, v, color=COL[g], lw=0.7)
        point(ax, YEARS, v, g)
    end_labels(ax, {g: res['case_mix']['2024'][g][f] for g in G}, 0.075 * (lim[1] - lim[0]), size=6)
    year_axis(ax)
    ax.set_xlim(2018.6, 2025.9)
    ax.set_ylim(*lim)
    ax.set_ylabel(yl)
    panel(ax, letter)
save(fig, 'figureS2_case_mix')
print('Wrote figures to figures/ (font:', FONT + ')')
