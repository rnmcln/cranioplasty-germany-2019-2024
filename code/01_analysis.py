"""
Cranioplasty in Germany, 2019-2024: material-group analysis of aggregate
inpatient data (InEK DatenBrowser, section 21 KHEntgG).

Input  : data/input/mortality.csv, complications.csv (T81.4, T85.72, T85.81), demographics.csv,
         overlap.csv, secondary_diagnoses_autologous_2023.csv,
         drg_mix_autologous_2023.csv
Output : output/analysis_results.json, output/tables/*.csv

Material groups
    Autologous          OPS 5-020.2   (reinsertion of the patient's own bone flap)
    Simple alloplastic  OPS 5-020.70  without a concurrent CAD code (derived)
    CAD                 OPS 5-020.71/.72/.74/.75 (computer-aided-design implants)
    Alloplastic total   any of 5-020.70/.71/.72/.74/.75 (descriptive parent group)

The DatenBrowser returns the number of discharge episodes carrying at least one
of the selected codes. The simple alloplastic group is therefore obtained as
alloplastic total minus CAD, which equals the episodes coded 5-020.70 and no CAD
code. Counts are subtracted directly; attribute means and percentages are
derived by case-weighted subtraction.

Episodes coded with both 5-020.2 and a 5-020.7 code (overlap.csv) belong to the
autologous group and to an alloplastic group. A sensitivity analysis removes them from
both and recomputes rates and risk ratios over every combination of values that the
suppressed overlap cells (1-4 episodes each) can take.

Statistics (Python standard library only, so that results do not depend on
library versions): Wilson 95% confidence intervals; risk ratios with
log-method 95% confidence intervals; Pearson chi-square test without continuity
correction (Fisher's exact test if any cell is below 5); Holm adjustment across
the three mortality contrasts; Cochran-Armitage test for linear trend in annual
mortality.

Run from the repository root:  python code/01_analysis.py
"""
import csv
import json
import math
import os
from itertools import combinations

IN = 'data/input/'
OUT = 'output/'
YEARS = [2019, 2020, 2021, 2022, 2023, 2024]
GROUPS = ['Autologous', 'Simple alloplastic', 'CAD']
ALLO = 'Alloplastic total'
SRC_ALLO = 'Broad alloplastic'   # label of the alloplastic total in the input files
Z = 1.959963984540054            # standard normal 0.975 quantile


# ---------------------------------------------------------------- input
def load(fn):
    with open(IN + fn, newline='') as f:
        return list(csv.DictReader(f))


def num(v):
    return None if v in ('', 'None', None) else float(v)


mort = load('mortality.csv')
comp = load('complications.csv')
demo = load('demographics.csv')
for r in mort:
    r['year'], r['N'], r['deaths'] = int(r['year']), int(r['N']), int(r['deaths'])
for r in comp:
    for k in ('year', 'N', 't85_72_n', 't81_4_n'):
        r[k] = int(r[k])
for r in demo:
    r['year'] = int(r['year'])
    for k in ('N', 'male_pct', 'female_pct', 'mean_los', 'lang_pct', 'pccl0_pct', 'adult_pct'):
        r[k] = num(r[k])


def row(rows, year, tier):
    hits = [r for r in rows if r['year'] == year and r['tier'] == tier]
    if len(hits) != 1:
        raise ValueError(f'expected one row for {tier} {year}, found {len(hits)}')
    return hits[0]


# ---------------------------------------------------------------- statistics
def wilson(k, n):
    p = k / n
    d = 1 + Z ** 2 / n
    c = (p + Z ** 2 / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z ** 2 / (4 * n ** 2)) / d
    return [round(100 * p, 2), round(100 * (c - h), 2), round(100 * (c + h), 2)]


def risk_ratio(k1, n1, k2, n2):
    rr = (k1 / n1) / (k2 / n2)
    se = math.sqrt(1 / k1 - 1 / n1 + 1 / k2 - 1 / n2)
    return rr, rr * math.exp(-Z * se), rr * math.exp(Z * se)


def chi2_p(k1, n1, k2, n2):
    a, b, c, d = k1, n1 - k1, k2, n2 - k2
    n = a + b + c + d
    x2 = n * (a * d - b * c) ** 2 / ((a + b) * (c + d) * (a + c) * (b + d))
    return math.erfc(math.sqrt(x2 / 2))          # upper tail of chi-square, 1 df


def fisher_p(k1, n1, k2, n2):
    m, n_tot = k1 + k2, n1 + n2
    lo, hi = max(0, m - n2), min(m, n1)
    pmf = {x: math.comb(n1, x) * math.comb(n2, m - x) / math.comb(n_tot, m) for x in range(lo, hi + 1)}
    obs = pmf[k1]
    return min(1.0, sum(v for v in pmf.values() if v <= obs * (1 + 1e-7)))


def test(k1, n1, k2, n2):
    if min(k1, n1 - k1, k2, n2 - k2) < 5:
        return fisher_p(k1, n1, k2, n2), 'Fisher'
    return chi2_p(k1, n1, k2, n2), 'chi-square'


def holm(p):
    order = sorted(range(len(p)), key=lambda i: p[i])
    adj, running = [0.0] * len(p), 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[i]))
        adj[i] = running
    return adj


def cochran_armitage(ks, ns, xs):
    N, K = sum(ns), sum(ks)
    p = K / N
    xbar = sum(n * x for n, x in zip(ns, xs)) / N
    t = sum(k * (x - xbar) for k, x in zip(ks, xs))
    v = p * (1 - p) * sum(n * (x - xbar) ** 2 for n, x in zip(ns, xs))
    z = t / math.sqrt(v)
    return z, math.erfc(abs(z) / math.sqrt(2))


def sig(x, digits=10):
    """Round to significant digits so that stored values do not depend on platform libm."""
    return float(f'{x:.{digits}g}')


def fmt_p(p):
    return '<0.001' if p < 0.001 else f'{p:.3f}'


# ---------------------------------------------------------------- group values
def deaths(y, g):
    if g == 'Simple alloplastic':
        b, c = row(mort, y, SRC_ALLO), row(mort, y, 'CAD')
        return b['deaths'] - c['deaths'], b['N'] - c['N']
    if g == ALLO:
        r = row(mort, y, SRC_ALLO)
        return r['deaths'], r['N']
    r = row(mort, y, g)
    return r['deaths'], r['N']


def infection(y, g, key):
    if g == 'Simple alloplastic':
        b, c = row(comp, y, SRC_ALLO), row(comp, y, 'CAD')
        return b[key] - c[key], b['N'] - c['N']
    r = row(comp, y, SRC_ALLO if g == ALLO else g)
    return r[key], r['N']


def attribute(y, g, field):
    if g == 'Simple alloplastic':
        b, c = row(demo, y, SRC_ALLO), row(demo, y, 'CAD')
        if b[field] is None or c[field] is None:
            return None
        n = b['N'] - c['N']
        return round((b[field] * b['N'] - c[field] * c['N']) / n, 1)
    r = row(demo, y, SRC_ALLO if g == ALLO else g)
    return None if r[field] is None else round(r[field], 1)


# ---------------------------------------------------------------- consistency of inputs
for y in YEARS:
    for rows in (mort, comp, demo):
        assert row(rows, y, SRC_ALLO)['N'] >= row(rows, y, 'CAD')['N']
    assert row(mort, y, 'Autologous')['N'] == row(comp, y, 'Autologous')['N'] == row(demo, y, 'Autologous')['N']
    assert row(mort, y, 'CAD')['N'] == row(comp, y, 'CAD')['N'] == row(demo, y, 'CAD')['N']
    assert row(mort, y, SRC_ALLO)['N'] == row(comp, y, SRC_ALLO)['N'] == row(demo, y, SRC_ALLO)['N']
    for key in ('t85_72_n', 't81_4_n'):
        assert row(comp, y, SRC_ALLO)[key] >= row(comp, y, 'CAD')[key]
    assert row(mort, y, SRC_ALLO)['deaths'] >= row(mort, y, 'CAD')['deaths']

ALLG = GROUPS + [ALLO]
res = {'groups': GROUPS, 'years': YEARS}

# ---------------------------------------------------------------- volumes
vol = {y: {g: deaths(y, g)[1] for g in ALLG} for y in YEARS}
res['volumes'] = vol
res['volumes_total_2019_2024'] = {g: sum(vol[y][g] for y in YEARS) for g in ALLG}
res['volume_change_2019_2024_pct'] = {g: round(100 * (vol[2024][g] - vol[2019][g]) / vol[2019][g], 1) for g in ALLG}
res['autologous_share_pct'] = {y: round(100 * vol[y]['Autologous'] / (vol[y]['Autologous'] + vol[y][ALLO]), 1) for y in YEARS}

# ---------------------------------------------------------------- mortality
res['mortality_by_year'] = {y: {g: dict(zip(('deaths', 'N'), deaths(y, g)), pct_ci=wilson(*deaths(y, g))) for g in ALLG} for y in YEARS}
pool = {g: (sum(deaths(y, g)[0] for y in YEARS), sum(deaths(y, g)[1] for y in YEARS)) for g in ALLG}
res['mortality_pooled'] = {g: {'deaths': pool[g][0], 'N': pool[g][1], 'pct_ci': wilson(*pool[g])} for g in ALLG}
res['cad_lowest_every_year'] = all(
    res['mortality_by_year'][y]['CAD']['pct_ci'][0] < min(res['mortality_by_year'][y][g]['pct_ci'][0] for g in ('Autologous', 'Simple alloplastic'))
    for y in YEARS)
res['autologous_above_simple_years'] = [y for y in YEARS if deaths(y, 'Autologous')[0] / deaths(y, 'Autologous')[1] > deaths(y, 'Simple alloplastic')[0] / deaths(y, 'Simple alloplastic')[1]]

pairs = list(combinations(GROUPS, 2))
raw, info = [], []
for a, b in pairs:
    p, t = test(*pool[a], *pool[b])
    raw.append(p)
    info.append((a, b, risk_ratio(*pool[a], *pool[b]), t))
res['mortality_pairwise'] = {}
for (a, b, rr, t), p, pa in zip(info, raw, holm(raw)):
    res['mortality_pairwise'][f'{a} vs {b}'] = {'RR': round(rr[0], 2), 'RR_CI': [round(rr[1], 2), round(rr[2], 2)],
                                                'p_raw': sig(p), 'p_holm': sig(pa), 'p_holm_fmt': fmt_p(pa), 'test': t}
res['mortality_trend'] = {}
for g in GROUPS:
    z, p = cochran_armitage([deaths(y, g)[0] for y in YEARS], [deaths(y, g)[1] for y in YEARS], YEARS)
    res['mortality_trend'][g] = {'z': round(z, 3), 'p': round(p, 3), 'annual_pct': [res['mortality_by_year'][y][g]['pct_ci'][0] for y in YEARS]}

# ---------------------------------------------------------------- coded infection
res['infection_by_year'] = {}
res['infection_pooled'] = {}
for code, key in (('T81.4', 't81_4_n'), ('T85.72', 't85_72_n')):
    res['infection_by_year'][code] = {y: {g: dict(zip(('k', 'N'), infection(y, g, key)), pct_ci=wilson(*infection(y, g, key))) for g in ALLG} for y in YEARS}
    pv = {g: (sum(infection(y, g, key)[0] for y in YEARS), sum(infection(y, g, key)[1] for y in YEARS)) for g in ALLG}
    block = {g: {'k': pv[g][0], 'N': pv[g][1], 'pct_ci': wilson(*pv[g])} for g in ALLG}
    for a, b in (('CAD', 'Autologous'), ('Simple alloplastic', 'Autologous'), ('CAD', 'Simple alloplastic')):
        rr = risk_ratio(*pv[a], *pv[b])
        p, t = test(*pv[a], *pv[b])
        block[f'{a} vs {b}'] = {'RR': round(rr[0], 2), 'RR_CI': [round(rr[1], 2), round(rr[2], 2)], 'p': sig(p), 'p_fmt': fmt_p(p), 'test': t}
    res['infection_pooled'][code] = block

# ---------------------------------------------------------------- case mix
FIELDS = ['male_pct', 'mean_los', 'lang_pct', 'pccl0_pct', 'adult_pct']
res['case_mix'] = {y: {g: dict(N=vol[y][g], **{f: attribute(y, g, f) for f in FIELDS}) for g in ALLG} for y in YEARS}

# ---------------------------------------------------------------- secondary diagnoses, autologous 2023
sdx = load('secondary_diagnoses_autologous_2023.csv')
res['secondary_dx_autologous_2023'] = [{'icd10gm': r['icd10gm'], 'label': r['label'], 'n': int(r['n']), 'N': int(r['N']),
                                        'pct': round(100 * int(r['n']) / int(r['N']), 1)} for r in sdx]
assert all(int(r['N']) == row(mort, 2023, 'Autologous')['N'] for r in sdx)

# ---------------------------------------------------------------- DRG mix, autologous 2023 (context for cost data)
res['drg_mix_autologous_2023'] = {r['drg_family']: {'n': int(r['n']), 'N': int(r['N']), 'pct': round(100 * int(r['n']) / int(r['N']), 1)}
                                  for r in load('drg_mix_autologous_2023.csv')}

# ---------------------------------------------------------------- overlap and exclusive-group sensitivity analysis
# Episodes coded with both 5-020.2 and a 5-020.7 code were counted by combining the procedure
# field (5-020.2) with the further-procedure field (5-020.7 codes). A cell below the display
# threshold is recorded as 's'. The DatenBrowser reports an empty selection separately
# ("keine Faelle"), so a suppressed cell is at least 1; the smallest count displayed in the
# extraction is 5 (complications.csv), so a suppressed cell is at most 4.
SUPP = (1, 4)
ovl = load('overlap.csv')


def cell(v):
    return set(range(SUPP[0], SUPP[1] + 1)) if v == 's' else {int(v)}


def totals(col_all, col_cad):
    """Feasible pooled totals (x_all, x_cad) with x_cad <= x_all in every year."""
    feas = {(0, 0)}
    for r in ovl:
        yr = {(a, c) for a in cell(r[col_all]) for c in cell(r[col_cad]) if c <= a}
        feas = {(s[0] + a, s[1] + c) for s in feas for a, c in yr}
    return feas


for r in ovl:
    r['year'] = int(r['year'])
    for k in ('auto_N', 'overlap_alloplastic_N', 'overlap_cad_N'):
        r[k] = int(r[k])
    assert r['auto_N'] == row(mort, r['year'], 'Autologous')['N']
    assert r['overlap_cad_N'] <= r['overlap_alloplastic_N']
assert sorted(r['year'] for r in ovl) == YEARS

O_all = sum(r['overlap_alloplastic_N'] for r in ovl)
O_cad = sum(r['overlap_cad_N'] for r in ovl)
O_simple = O_all - O_cad
tot = res['volumes_total_2019_2024']
res['overlap'] = {
    'by_year': {r['year']: {'autologous_and_alloplastic': r['overlap_alloplastic_N'], 'autologous_and_CAD': r['overlap_cad_N'],
                            'autologous_and_simple': r['overlap_alloplastic_N'] - r['overlap_cad_N']} for r in ovl},
    'total': {'autologous_and_alloplastic': O_all, 'autologous_and_CAD': O_cad, 'autologous_and_simple': O_simple},
    'share_pct': {'Autologous': round(100 * O_all / tot['Autologous'], 1), 'Simple alloplastic': round(100 * O_simple / tot['Simple alloplastic'], 1),
                  'CAD': round(100 * O_cad / tot['CAD'], 1)},
    'suppressed_cell_bounds': list(SUPP),
}
excl_N = {'Autologous': tot['Autologous'] - O_all, 'Simple alloplastic': tot['Simple alloplastic'] - O_simple, 'CAD': tot['CAD'] - O_cad}


def excl_events(pooled_k, xa, xc):
    return {'Autologous': pooled_k['Autologous'] - xa, 'Simple alloplastic': pooled_k['Simple alloplastic'] - (xa - xc), 'CAD': pooled_k['CAD'] - xc}


def sens_block(pooled_k, feas, contrasts):
    out = {'overlap_events_range': {'autologous_and_alloplastic': [min(f[0] for f in feas), max(f[0] for f in feas)],
                                    'autologous_and_CAD': [min(f[1] for f in feas), max(f[1] for f in feas)]},
           'N': excl_N, 'pct_range': {}, 'contrasts': {}}
    scen = [excl_events(pooled_k, xa, xc) for xa, xc in sorted(feas)]
    for g in GROUPS:
        pcts = [wilson(s[g], excl_N[g])[0] for s in scen]
        out['pct_range'][g] = [min(pcts), max(pcts)]
    for a, b in contrasts:
        rrs = [risk_ratio(s[a], excl_N[a], s[b], excl_N[b]) for s in scen]
        ps = [test(s[a], excl_N[a], s[b], excl_N[b])[0] for s in scen]
        out['contrasts'][f'{a} vs {b}'] = {'RR_range': [round(min(r[0] for r in rrs), 2), round(max(r[0] for r in rrs), 2)],
                                           'CI_outer': [round(min(r[1] for r in rrs), 2), round(max(r[2] for r in rrs), 2)],
                                           'p_max': sig(max(ps)), 'p_max_fmt': fmt_p(max(ps))}
    return out


mk = {g: res['mortality_pooled'][g]['deaths'] for g in GROUPS}
res['sensitivity_exclusive'] = {'mortality': sens_block(mk, totals('overlap_alloplastic_deaths', 'overlap_cad_deaths'), pairs)}
for code, ca, cc in (('T81.4', 'overlap_alloplastic_t81_4', 'overlap_cad_t81_4'), ('T85.72', 'overlap_alloplastic_t85_72', 'overlap_cad_t85_72')):
    kk = {g: res['infection_pooled'][code][g]['k'] for g in GROUPS}
    res['sensitivity_exclusive'][code] = sens_block(kk, totals(ca, cc), (('CAD', 'Autologous'), ('Simple alloplastic', 'Autologous'), ('CAD', 'Simple alloplastic')))

# ---------------------------------------------------------------- T85.81, descriptive only
# T85.81 (other complications of nervous-system implants) is below the display threshold in most
# group-years ('s', 1-4 episodes). It is reported descriptively: queried groups as displayed, and
# the derived simple alloplastic count as the range compatible with the alloplastic total and CAD.
res['t85_81_by_year'] = {}
for y in YEARS:
    a, b, c = (row(comp, y, t)['t85_81_n'] for t in ('Autologous', SRC_ALLO, 'CAD'))
    simple = sorted({x - z for x in cell(b) for z in cell(c) if z <= x})
    assert simple, f'T85.81 {y}: CAD cannot exceed the alloplastic total'
    res['t85_81_by_year'][y] = {'Autologous': a, ALLO: b, 'CAD': c, 'Simple alloplastic_range': [simple[0], simple[-1]]}
res['t85_81_max_per_group_year'] = max(max(SUPP[1] if d[g] == 's' else int(d[g]) for g in ('Autologous', 'CAD')) for d in res['t85_81_by_year'].values())
res['t85_81_max_per_group_year'] = max([res['t85_81_max_per_group_year']] + [d['Simple alloplastic_range'][1] for d in res['t85_81_by_year'].values()])

# ---------------------------------------------------------------- write
os.makedirs(OUT + 'tables', exist_ok=True)
with open(OUT + 'analysis_results.json', 'w') as f:
    json.dump(res, f, indent=1, sort_keys=True)

with open(OUT + 'tables/mortality_by_year.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['year', 'group', 'deaths', 'N', 'pct', 'ci_low', 'ci_high'])
    for y in YEARS:
        for g in ALLG:
            m = res['mortality_by_year'][y][g]
            w.writerow([y, g, m['deaths'], m['N'], *m['pct_ci']])
with open(OUT + 'tables/infection_by_year.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['code', 'year', 'group', 'k', 'N', 'pct', 'ci_low', 'ci_high'])
    for code in ('T81.4', 'T85.72'):
        for y in YEARS:
            for g in ALLG:
                m = res['infection_by_year'][code][y][g]
                w.writerow([code, y, g, m['k'], m['N'], *m['pct_ci']])
with open(OUT + 'tables/case_mix.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['year', 'group', 'N'] + FIELDS)
    for y in YEARS:
        for g in ALLG:
            c = res['case_mix'][y][g]
            w.writerow([y, g, c['N']] + [c[k] for k in FIELDS])

# ---------------------------------------------------------------- console summary
print('Pooled in-hospital mortality, 2019-2024')
for g in ALLG:
    m = res['mortality_pooled'][g]
    print(f'  {g:20s} {m["deaths"]:4d}/{m["N"]:6d}  {m["pct_ci"][0]:.2f}% ({m["pct_ci"][1]:.2f}-{m["pct_ci"][2]:.2f})')
for k, v in res['mortality_pairwise'].items():
    print(f'  {k:34s} RR {v["RR"]:.2f} ({v["RR_CI"][0]:.2f}-{v["RR_CI"][1]:.2f})  Holm p {v["p_holm_fmt"]}')
for code in ('T81.4', 'T85.72'):
    print(f'Coded infection {code}, pooled 2019-2024')
    for g in ALLG:
        c = res['infection_pooled'][code][g]
        print(f'  {g:20s} {c["k"]:4d}/{c["N"]:6d}  {c["pct_ci"][0]:.2f}% ({c["pct_ci"][1]:.2f}-{c["pct_ci"][2]:.2f})')
print('Volume change 2019-2024 (%):', res['volume_change_2019_2024_pct'])
print('Overlap (autologous and alloplastic code in one episode):', res['overlap']['total'], res['overlap']['share_pct'])
for k, v in res['sensitivity_exclusive']['mortality']['contrasts'].items():
    print(f'  exclusive groups, mortality {k:34s} RR {v["RR_range"][0]:.2f}-{v["RR_range"][1]:.2f}  CI {v["CI_outer"][0]:.2f}-{v["CI_outer"][1]:.2f}')
print('Wrote', OUT + 'analysis_results.json')
