"""
Internal-consistency and external-benchmark checks.

1. Every group-year count is non-negative and every derived simple-alloplastic
   count (alloplastic total minus CAD) is non-negative.
2. Rates, confidence intervals and risk ratios in output/analysis_results.json are
   recomputed here with an independent implementation and must agree.
3. Overlap counts are consistent with the group counts and the exclusive-group
   denominators equal the pooled counts minus the overlap.
4. Annual autologous (5-020.2) and CAD (5-020.71/.72/.74/.75) counts are compared
   with the national procedure counts published by Hirsch et al. (2026), Table 1,
   https://doi.org/10.3390/brainsci16050444, for 2019-2023.

Exit status 0 = all checks passed.
Run from the repository root:  python code/03_validate.py
"""
import csv
import json
import math
import sys

res = json.load(open('output/analysis_results.json'))
YEARS = [2019, 2020, 2021, 2022, 2023, 2024]
G = ['Autologous', 'Simple alloplastic', 'CAD']
fails = []


def check(cond, msg):
    print(('  ok    ' if cond else '  FAIL  ') + msg)
    if not cond:
        fails.append(msg)


def load(fn):
    with open('data/input/' + fn, newline='') as f:
        return list(csv.DictReader(f))


mort, comp = load('mortality.csv'), load('complications.csv')

print('1. Non-negativity and additivity')
for y in YEARS:
    m = {r['tier']: r for r in mort if int(r['year']) == y}
    c = {r['tier']: r for r in comp if int(r['year']) == y}
    sd = int(m['Broad alloplastic']['deaths']) - int(m['CAD']['deaths'])
    sn = int(m['Broad alloplastic']['N']) - int(m['CAD']['N'])
    s1 = int(c['Broad alloplastic']['t81_4_n']) - int(c['CAD']['t81_4_n'])
    s2 = int(c['Broad alloplastic']['t85_72_n']) - int(c['CAD']['t85_72_n'])
    check(min(sd, sn, s1, s2) >= 0, f'{y}: simple alloplastic deaths={sd}, N={sn}, T81.4={s1}, T85.72={s2}')
    v = res['volumes'][str(y)]
    check(v['Simple alloplastic'] + v['CAD'] == v['Alloplastic total'], f'{y}: simple + CAD = alloplastic total')

print('2. Independent recomputation')
q = 1.959963984540054
for g in G:
    m = res['mortality_pooled'][g]
    k, n = m['deaths'], m['N']
    p = k / n
    # Wilson interval written as the roots of the score-test quadratic
    a = 1 + q * q / n
    b = -(2 * p + q * q / n)
    cc = p * p
    lo = (-b - math.sqrt(b * b - 4 * a * cc)) / (2 * a)
    hi = (-b + math.sqrt(b * b - 4 * a * cc)) / (2 * a)
    check(abs(100 * lo - m['pct_ci'][1]) < 0.006 and abs(100 * hi - m['pct_ci'][2]) < 0.006, f'mortality CI {g}')
    check(sum(res['mortality_by_year'][str(y)][g]['deaths'] for y in YEARS) == k, f'pooled deaths = sum of annual deaths, {g}')
for key, v in res['mortality_pairwise'].items():
    a_, b_ = key.split(' vs ')
    A, B = res['mortality_pooled'][a_], res['mortality_pooled'][b_]
    rr = (A['deaths'] / A['N']) / (B['deaths'] / B['N'])
    check(abs(round(rr, 2) - v['RR']) < 1e-9 and v['RR_CI'][0] < v['RR'] < v['RR_CI'][1], f'mortality RR {key}')
for code in ('T81.4', 'T85.72'):
    blk = res['infection_pooled'][code]
    for key in [k for k in blk if ' vs ' in k]:
        a_, b_ = key.split(' vs ')
        rr = (blk[a_]['k'] / blk[a_]['N']) / (blk[b_]['k'] / blk[b_]['N'])
        check(abs(round(rr, 2) - blk[key]['RR']) < 1e-9, f'{code} RR {key}')

print('2b. Overlap of autologous and alloplastic codes')
ov = res['overlap']
for y in YEARS:
    o = ov['by_year'][str(y)]
    check(o['autologous_and_CAD'] + o['autologous_and_simple'] == o['autologous_and_alloplastic'], f'{y}: overlap CAD + simple = overlap total')
    check(o['autologous_and_alloplastic'] <= res['volumes'][str(y)]['Autologous'] and o['autologous_and_CAD'] <= res['volumes'][str(y)]['CAD'], f'{y}: overlap within group counts')
for g in G:
    n_ex = res['sensitivity_exclusive']['mortality']['N'][g]
    o = {'Autologous': ov['total']['autologous_and_alloplastic'], 'Simple alloplastic': ov['total']['autologous_and_simple'], 'CAD': ov['total']['autologous_and_CAD']}[g]
    check(n_ex == res['volumes_total_2019_2024'][g] - o, f'exclusive N {g} = pooled N minus overlap')
for key, v in res['sensitivity_exclusive']['mortality']['contrasts'].items():
    check(v['RR_range'][0] <= v['RR_range'][1] and v['CI_outer'][0] < v['RR_range'][0], f'sensitivity RR range ordered, {key}')

print('2c. T85.81 (descriptive)')
for y in YEARS:
    d = res['t85_81_by_year'][str(y)]
    lo, hi = d['Simple alloplastic_range']
    check(0 <= lo <= hi, f'{y}: T85.81 simple alloplastic range {lo}-{hi} is feasible')

print('3. External benchmark (Hirsch et al. 2026, Table 1)')
published = {'Autologous': {2019: 3039, 2020: 2497, 2021: 2088, 2022: 2082, 2023: 2129},
             'CAD': {2019: 1112, 2020: 1050, 2021: 1022, 2022: 1134, 2023: 1284}}
for y, n in published['Autologous'].items():
    check(res['volumes'][str(y)]['Autologous'] == n, f'{y}: autologous {res["volumes"][str(y)]["Autologous"]} vs published {n}')
diffs = {y: n - res['volumes'][str(y)]['CAD'] for y, n in published['CAD'].items()}
check(all(0 <= d <= 10 for d in diffs.values()), f'CAD differences (published minus present): {diffs}')

print('\nVALIDATION', 'PASSED' if not fails else f'FAILED ({len(fails)} checks)')
sys.exit(1 if fails else 0)
