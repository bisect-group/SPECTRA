#!/usr/bin/env python3
"""finalize_stage3_mcmc_and_plots.py
Runs automatically after Stage 2B (80-chain RHMC MCMC sampling on the 20 newly
rebuilt SPECTRA_ME models) and Stage 2A (20-model FVA + Stage 1 1,479 models) finish:
  1. Merges all 80 RHMC MCMC chains (4 chains x 20 models, 2,500 ESS/chain -> ~10,000 ESS),
     computing Gelman-Rubin R-hat, ESS, two-sided Mann-Whitney U on |v| and signed v,
     Cliff's delta, and BH FDR -> saves sampling_per_reaction.csv & sampling_summary.csv.
  2. Computes reaction-level and pathway-level corroboration between the new FVA FSR/FEA
     calls and the new RHMC MCMC flux distributions -> saves step4_rxn_fsr_mcmc_concordance.csv
     and step4_new_fsr_pathways_with_mcmc_corroboration.csv.
  3. Exports updated multi-sheet Data_S4, Data_S5, and Data_S6 Excel workbooks.
  4. Renders updated Fig. 5d (bubble_fea_orange_sleek_new.pdf / .png) with 1.2 pt
     Burnt Orange (#E66100) MCMC corroboration rings.
"""
import json
import os
import zipfile
from xml.sax.saxutils import escape
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import mannwhitneyu, hypergeom

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

ROOT_Q44 = '/scratch/pavan/q_4_4'
RERUN_ROOT = '/scratch/pavan/q_4_4/rerun_69d2281'
CHAINS_DIR = os.path.join(RERUN_ROOT, 'sampling', 'chains')
FVA_DIR = os.path.join(RERUN_ROOT, 'fva')
OUT_DIR = os.path.join(RERUN_ROOT, 'results')
FIG_DIR = os.path.join(RERUN_ROOT, 'figures')
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'supp_data'), exist_ok=True)

CONTEXTS = ['Brain', 'Kidney', 'Liver', 'Lung', 'Eye',
            'Pancreas', 'Prostate', 'Breast', 'Ovary', 'Biliary Tract']
CANCER = ['ACH-000376', 'ACH-000649', 'ACH-000739', 'ACH-000681', 'ACH-001295',
          'ACH-000108', 'ACH-000177', 'ACH-000019', 'ACH-000001', 'ACH-001961']
NORMAL = ['ACH-003150', 'ACH-003143', 'ACH-003155', 'ACH-000064', 'ACH-002466',
          'ACH-003161', 'ACH-001453', 'ACH-000642', 'ACH-000494', 'ACH-001854']

NCHAIN = 4
DELTA_CUT = 0.33
DELTA_GRID = [0.147, 0.33, 0.474, 0.6]
FDR_CUT = 0.05
CHUNK = 400


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum.accumulate((ranked * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out


def load_chains(ach):
    mats, rxns, counts = [], None, []
    for c in range(1, NCHAIN + 1):
        f = os.path.join(CHAINS_DIR, f'{ach}_chain{c}.mat')
        if not os.path.exists(f):
            raise FileNotFoundError(f)
        d = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
        s = np.asarray(d['samples'], dtype=np.float32)
        r = [str(x) for x in np.atleast_1d(d['rxns'])]
        if rxns is None:
            rxns = r
        elif r != rxns:
            raise ValueError(f'reaction order differs between chains of {ach}')
        mats.append(s)
        counts.append(s.shape[1])
    return rxns, np.hstack(mats), counts


def split_rhat(mats):
    n = min(m.shape[1] for m in mats)
    x = np.stack([m[:, :n] for m in mats], axis=0).astype(np.float64)
    m_chain = x.mean(axis=2)
    v_chain = x.var(axis=2, ddof=1)
    W = v_chain.mean(axis=0)
    B = n * m_chain.var(axis=0, ddof=1)
    var_hat = (n - 1) / n * W + B / n
    with np.errstate(divide='ignore', invalid='ignore'):
        rh = np.sqrt(var_hat / W)
    rh[~np.isfinite(rh)] = 1.0
    return rh


def ess_1d(x):
    n = x.size
    x = x - x.mean()
    v = np.dot(x, x) / n
    if v <= 0:
        return float(n)
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n] / (n * v)
    s, t = 0.0, 1
    while t + 1 < n:
        pair = ac[t] + ac[t + 1]
        if pair <= 0:
            break
        s += pair
        t += 2
    return float(min(n, n / (1 + 2 * s)))


def chain_ess(mats, idx):
    out = np.zeros(len(idx))
    for m in mats:
        for j, i in enumerate(idx):
            out[j] += ess_1d(m[i, :].astype(np.float64))
    return out


def thin_to(x, k):
    n = x.shape[1]
    if k >= n:
        return x
    return x[:, np.linspace(0, n - 1, k).astype(int)]


def merge_mcmc_all_pairs():
    print("=== 1. Merging 80 RHMC MCMC Chains & Running Distribution Tests ===", flush=True)
    frames, summary = [], []
    for ctx, c, n in zip(CONTEXTS, CANCER, NORMAL):
        print(f"  [{ctx}] {c} vs {n} ...", flush=True)
        rxn_c, Xc, cnt_c = load_chains(c)
        rxn_n, Xn, cnt_n = load_chains(n)
        mats_c = np.split(Xc, np.cumsum(cnt_c)[:-1], axis=1)
        mats_n = np.split(Xn, np.cumsum(cnt_n)[:-1], axis=1)
        rhat_c = split_rhat(mats_c)
        rhat_n = split_rhat(mats_n)

        common = sorted(set(rxn_c) & set(rxn_n))
        common = [r for r in common if not (r.startswith('EX_') or r.startswith('sink_') or r.startswith('DM_'))]
        ic = np.array([rxn_c.index(r) for r in common])
        inn = np.array([rxn_n.index(r) for r in common])

        ess_c = chain_ess(mats_c, ic)
        ess_n = chain_ess(mats_n, inn)
        n_eff = np.clip(np.minimum(ess_c, ess_n).astype(int), 10, min(Xc.shape[1], Xn.shape[1]))

        recs = []
        for start in range(0, len(common), CHUNK):
            sl = slice(start, min(start + CHUNK, len(common)))
            A = np.abs(Xc[ic[sl], :]).astype(np.float64)
            B = np.abs(Xn[inn[sl], :]).astype(np.float64)
            Asig = Xc[ic[sl], :].astype(np.float64)
            Bsig = Xn[inn[sl], :].astype(np.float64)
            for j in range(A.shape[0]):
                gi = start + j
                k = int(n_eff[gi])
                a = thin_to(A[j][None, :], k)[0]
                b = thin_to(B[j][None, :], k)[0]
                asg = thin_to(Asig[j][None, :], k)[0]
                bsg = thin_to(Bsig[j][None, :], k)[0]
                if np.ptp(a) == 0 and np.ptp(b) == 0 and a[0] == b[0]:
                    p, delta = 1.0, 0.0
                else:
                    u, p = mannwhitneyu(a, b, alternative='two-sided')
                    delta = 2.0 * u / (a.size * b.size) - 1.0
                if np.ptp(asg) == 0 and np.ptp(bsg) == 0 and asg[0] == bsg[0]:
                    p_s, delta_s = 1.0, 0.0
                else:
                    us, p_s = mannwhitneyu(asg, bsg, alternative='two-sided')
                    delta_s = 2.0 * us / (asg.size * bsg.size) - 1.0
                recs.append({
                    'context': ctx, 'rxn': common[gi], 'n_eff': k,
                    'ess_cancer': ess_c[gi], 'ess_normal': ess_n[gi],
                    'rhat_cancer': rhat_c[ic[gi]], 'rhat_normal': rhat_n[inn[gi]],
                    'median_abs_cancer': float(np.median(a)), 'median_abs_normal': float(np.median(b)),
                    'mean_abs_cancer': float(a.mean()), 'mean_abs_normal': float(b.mean()),
                    'median_signed_cancer': float(np.median(asg)), 'median_signed_normal': float(np.median(bsg)),
                    'p_abs': float(p), 'cliffs_delta_abs': float(delta),
                    'p_signed': float(p_s), 'cliffs_delta_signed': float(delta_s),
                })
        df = pd.DataFrame(recs)
        df['fdr_abs'] = bh_fdr(df.p_abs.values)
        df['fdr_signed'] = bh_fdr(df.p_signed.values)
        df['call_up'] = (df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs >= DELTA_CUT)
        df['call_down'] = (df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs <= -DELTA_CUT)
        frames.append(df)

        row = {'context': ctx, 'cancer': c, 'normal': n, 'n_common': len(df),
               'n_up': int(df.call_up.sum()), 'n_down': int(df.call_down.sum()),
               'n_fdr_sig': int((df.fdr_abs < FDR_CUT).sum()),
               'max_rhat': float(max(df.rhat_cancer.max(), df.rhat_normal.max())),
               'frac_rhat_gt_1p05': float(((df.rhat_cancer > 1.05) | (df.rhat_normal > 1.05)).mean()),
               'min_ess': float(min(df.ess_cancer.min(), df.ess_normal.min())),
               'median_ess': float(np.median(np.r_[df.ess_cancer, df.ess_normal]))}
        for g in DELTA_GRID:
            row[f'up_d{g:.3f}'] = int(((df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs >= g)).sum())
            row[f'down_d{g:.3f}'] = int(((df.fdr_abs < FDR_CUT) & (df.cliffs_delta_abs <= -g)).sum())
        summary.append(row)
        print(f"    common={len(df)} | FDR<0.05={row['n_fdr_sig']} | up={row['n_up']} | dn={row['n_down']} | max Rhat={row['max_rhat']:.4f}", flush=True)

    allf = pd.concat(frames, ignore_index=True)
    allf.to_csv(os.path.join(OUT_DIR, 'sampling_per_reaction.csv'), index=False)
    sm = pd.DataFrame(summary)
    sm.to_csv(os.path.join(OUT_DIR, 'sampling_summary.csv'), index=False)
    return allf, sm


def save_xlsx(path, sheets_dict):
    def col_letter(idx):
        s, n = '', idx + 1
        while n > 0:
            n, r = divmod(n - 1, 26)
            s = chr(65 + r) + s
        return s

    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as zf:
        sheet_names = list(sheets_dict.keys())
        overrides = ''.join(
            f'<Override PartName="/xl/worksheets/sheet{i+1}.xml" '
            f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            for i in range(len(sheet_names))
        )
        zf.writestr('[Content_Types].xml',
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                    '<Default Extension="xml" ContentType="application/xml"/>'
                    '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                    + overrides + '</Types>')
        zf.writestr('_rels/.rels',
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                    '</Relationships>')
        wb_sheets = ''.join(f'<sheet name="{escape(nm)}" sheetId="{i+1}" r:id="rId{i+1}"/>' for i, nm in enumerate(sheet_names))
        zf.writestr('xl/workbook.xml',
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
                    f'<sheets>{wb_sheets}</sheets></workbook>')
        wb_rels = ''.join(f'<Relationship Id="rId{i+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i+1}.xml"/>' for i in range(len(sheet_names)))
        zf.writestr('xl/_rels/workbook.xml.rels',
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    f'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">{wb_rels}</Relationships>')
        for i, nm in enumerate(sheet_names):
            row_xml = []
            for r_idx, row in enumerate(sheets_dict[nm]):
                cells = []
                for c_idx, val in enumerate(row):
                    ref = f"{col_letter(c_idx)}{r_idx+1}"
                    if val is None or (isinstance(val, float) and (np.isnan(val) or np.isinf(val))):
                        continue
                    elif isinstance(val, (int, float, np.integer, np.floating)) and not isinstance(val, bool):
                        cells.append(f'<c r="{ref}"><v>{float(val):.12g}</v></c>')
                    else:
                        cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{escape(str(val))}</t></is></c>')
                row_xml.append(f'<row r="{r_idx+1}">{"".join(cells)}</row>')
            zf.writestr(f'xl/worksheets/sheet{i+1}.xml',
                        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                        f'<sheetData>{"".join(row_xml)}</sheetData></worksheet>')


def run_corroboration_and_export(sampl_df):
    print("\n=== 2. Running FSR + Upper-Tail FEA + RHMC MCMC Corroboration ===", flush=True)
    S_rec = sio.loadmat(os.path.join(ROOT_Q44, 'Results_BuildDepMapModels.mat'),
                        squeeze_me=True, struct_as_record=False)['model']
    u_rxns = [str(x) for x in np.atleast_1d(S_rec.rxns)]
    u_names = [str(x) for x in np.atleast_1d(S_rec.rxnNames)]
    u_gr = [str(x) if not (isinstance(x, np.ndarray) and x.size == 0) else '' for x in np.atleast_1d(S_rec.grRules)]
    u_subs = [str(s.item()) if isinstance(s, np.ndarray) and s.size == 1 else str(s) for s in np.atleast_1d(S_rec.subSystems)]
    u_mets = [str(x) for x in np.atleast_1d(S_rec.mets)]
    S_csc = S_rec.S.tocsc()
    lb_u = np.asarray(S_rec.lb, dtype=float).ravel()
    ub_u = np.asarray(S_rec.ub, dtype=float).ravel()

    rxn_meta = {}
    for j, r in enumerate(u_rxns):
        p0, p1 = S_csc.indptr[j], S_csc.indptr[j + 1]
        subs_p, prod_p = [], []
        for mi, val in zip(S_csc.indices[p0:p1], S_csc.data[p0:p1]):
            tok = f"{abs(val):g} {u_mets[mi]}" if abs(abs(val) - 1.0) > 1e-9 else u_mets[mi]
            if val < 0:
                subs_p.append(tok)
            elif val > 0:
                prod_p.append(tok)
        arrow = ' <=> ' if (lb_u[j] < 0 and ub_u[j] > 0) else ' -> '
        rxn_meta[r] = (u_names[j], ' + '.join(subs_p) + arrow + ' + '.join(prod_p), u_gr[j], u_subs[j])

    nb = json.load(open('/scratch/pavan/spectra/Visualisations/cancer_spectra_pathway.ipynb'))
    code = ''.join(nb['cells'][0]['source']).replace('\u2011', '-')
    ns = {'pd': pd, 'np': np}
    exec(code, ns)
    df_cur = ns['df'].copy()
    cmap = {
        'Glioblastoma': 'Brain', 'Renal Clear Cell \n Carcinoma': 'Kidney',
        'Liver cancer': 'Liver', 'Lung Adenocarcinoma': 'Lung',
        'Retinoblastoma': 'Eye', 'Adenosquamous Carcinoma \n (Pancreas)': 'Pancreas',
        'Prostate Small \n Cell Carcinoma': 'Prostate', 'Invasive Breast \n Carcinoma': 'Breast',
        'High-Grade Serous \n Ovarian Cancer': 'Ovary', 'Gallbladder \n Cancer': 'Biliary Tract',
    }
    df_cur['context'] = df_cur['cancer'].map(cmap)
    ev_map = {(r['context'], r['pathway']): r['Evidence level'] for _, r in df_cur.iterrows()}
    ev_map.update({
        ('Biliary Tract', 'Heme synthesis'): 'Uncertain',
        # Newly enriched pathways under rebuilt SPECTRA_ME (commit 69d2281), literature-curated:
        ('Brain', 'Keratan sulfate degradation'): 'Component-level',
        ('Brain', 'Chondroitin sulfate degradation'): 'Component-level',
        ('Brain', 'Fatty acid synthesis'): 'Pathway-level',
        ('Brain', 'Histidine metabolism'): 'Component-level',
        ('Brain', 'Hyaluronan metabolism'): 'Pathway-level',
        ('Brain', 'Tryptophan metabolism'): 'Pathway-level',
        ('Kidney', 'Tryptophan metabolism'): 'Pathway-level',
        ('Kidney', 'Heme synthesis'): 'Pathway-level',
        ('Kidney', 'Nucleotide interconversion'): 'Pathway-level',
        ('Liver', 'Aminoacyl-tRNA biosynthesis'): 'Pathway-level',
        ('Liver', 'Purine synthesis'): 'Pathway-level',
        ('Liver', 'Protein degradation'): 'Component-level',
        ('Liver', 'Protein assembly'): 'Component-level',
        ('Liver', 'Hyaluronan metabolism'): 'Pathway-level',
        ('Liver', 'Citric acid cycle'): 'Pathway-level',
        ('Lung', 'Keratan sulfate degradation'): 'Component-level',
        ('Lung', 'Chondroitin sulfate degradation'): 'Component-level',
        ('Lung', 'Fatty acid synthesis'): 'Pathway-level',
        ('Lung', 'Eicosanoid metabolism'): 'Pathway-level',
        ('Lung', 'Hyaluronan metabolism'): 'Pathway-level',
        ('Eye', 'Lysine metabolism'): 'Pathway-level',
        ('Eye', 'Fatty acid synthesis'): 'Pathway-level',
        ('Eye', 'Histidine metabolism'): 'Pathway-level',
        ('Pancreas', 'N-glycan degradation'): 'Pathway-level',
        ('Pancreas', 'Glyoxylate and dicarboxylate metabolism'): 'Component-level',
        ('Prostate', 'Keratan sulfate degradation'): 'Uncertain',
        ('Prostate', 'Chondroitin sulfate degradation'): 'Component-level',
        ('Prostate', 'Purine synthesis'): 'Pathway-level',
        ('Prostate', 'Histidine metabolism'): 'Component-level',
        ('Prostate', 'Glycerophospholipid metabolism'): 'Pathway-level',
        ('Prostate', 'Hyaluronan metabolism'): 'Pathway-level',
        ('Prostate', 'Pyrimidine synthesis'): 'Pathway-level',
        ('Ovary', 'Keratan sulfate degradation'): 'Component-level',
        ('Ovary', 'Chondroitin sulfate degradation'): 'Component-level',
        ('Ovary', 'Heparan sulfate degradation'): 'Pathway-level',
        ('Ovary', 'Chondroitin synthesis'): 'Component-level',
        ('Ovary', 'Hyaluronan metabolism'): 'Pathway-level',
        ('Ovary', 'Sphingolipid metabolism'): 'Pathway-level',
        ('Ovary', 'Glyoxylate and dicarboxylate metabolism'): 'Component-level',
        ('Biliary Tract', 'Eicosanoid metabolism'): 'Pathway-level',
        ('Biliary Tract', 'Histidine metabolism'): 'Component-level',
        ('Biliary Tract', 'Hyaluronan metabolism'): 'Component-level',
        ('Biliary Tract', 'Valine, leucine, and isoleucine metabolism'): 'Component-level',
    })

    rxn_conc_rows = []
    pw_conc_rows = []
    up_sheets, dn_sheets, fea_sheets = {}, {}, {}

    for ctx, c_id, n_id in zip(CONTEXTS, CANCER, NORMAL):
        fc = pd.read_csv(os.path.join(FVA_DIR, c_id + '.csv')).set_index('rxn')
        fn = pd.read_csv(os.path.join(FVA_DIR, n_id + '.csv')).set_index('rxn')
        c_subs = pd.Series([rxn_meta[r][3] for r in fc.index], index=fc.index)
        common = [r for r in fc.index.intersection(fn.index)
                  if not (r.startswith('EX_') or r.startswith('sink_') or r.startswith('DM_'))]
        sub_s = sampl_df[sampl_df['context'] == ctx].set_index('rxn')

        up_rows = [['Upregulated_disease', 'fsr_Upregulated_disease', 'Upregulated_disease_names',
                    'Upregulated_disease_formulae', 'Upregulated_disease_genes', 'Upregulated_disease_subsystems',
                    'minFlux_cancer', 'maxFlux_cancer', 'minFlux_normal', 'maxFlux_normal',
                    'mcmc_median_abs_cancer', 'mcmc_median_abs_normal', 'mcmc_cliffs_delta_abs', 'mcmc_fdr_abs', 'mcmc_concordant']]
        dn_rows = [['Downregulated_disease', 'fsr_Downregulated_disease', 'Down_disease_names',
                    'Down_disease_formulae', 'Down_disease_genes', 'Down_disease_subsystems',
                    'minFlux_cancer', 'maxFlux_cancer', 'minFlux_normal', 'maxFlux_normal',
                    'mcmc_median_abs_cancer', 'mcmc_median_abs_normal', 'mcmc_cliffs_delta_abs', 'mcmc_fdr_abs', 'mcmc_concordant']]

        up_set, dn_set, clean2b = set(), set(), set()
        n_mcmc_conc = 0

        for r in common:
            min_c, max_c = float(fc.loc[r, 'minimum']), float(fc.loc[r, 'maximum'])
            min_n, max_n = float(fn.loc[r, 'minimum']), float(fn.loc[r, 'maximum'])
            span_c, span_n = max_c - min_c, max_n - min_n
            if span_n <= 0 or not np.isfinite(span_c / span_n):
                continue
            fsr = span_c / span_n
            is_up = (max_c > max_n) and (2.0 <= fsr < 20000.0)
            is_dn = (max_c < max_n) and (0.01 <= fsr <= 0.5)
            if not (is_up or is_dn):
                continue

            if is_up:
                up_set.add(r)
                if (max_c - max_n) >= (min_n - min_c) - 1e-9:
                    clean2b.add(r)
            else:
                dn_set.add(r)
                clean2b.add(r)

            rname, rform, rgr, rsub = rxn_meta.get(r, ('', '', '', ''))
            if r in sub_s.index:
                sr = sub_s.loc[r]
                m_c, m_n = float(sr['median_abs_cancer']), float(sr['median_abs_normal'])
                cd, fdr_m = float(sr['cliffs_delta_abs']), float(sr['fdr_abs'])
                conc = bool((fdr_m < 0.05) and ((cd > 0) if is_up else (cd < 0)))
            else:
                m_c, m_n, cd, fdr_m, conc = np.nan, np.nan, np.nan, np.nan, False
            if conc:
                n_mcmc_conc += 1

            rv = [r, fsr, rname, rform, rgr, rsub, min_c, max_c, min_n, max_n, m_c, m_n, cd, fdr_m, str(conc)]
            if is_up:
                up_rows.append(rv)
            else:
                dn_rows.append(rv)

        up_sheets[ctx] = up_rows
        dn_sheets[ctx] = dn_rows
        diff_rxns = list(up_set | dn_set)
        tot_fsr = len(diff_rxns)
        rxn_conc_rows.append({
            'context': ctx, 'cancer': c_id, 'normal': n_id,
            'fsr_up': len(up_set), 'fsr_dn': len(dn_set), 'fsr_tot': tot_fsr,
            'mcmc_concordant_n': n_mcmc_conc,
            'mcmc_concordant_pct': round(100.0 * n_mcmc_conc / max(tot_fsr, 1), 2),
        })

        # FEA on diff_rxns and on clean2b
        M = len(c_subs)
        all_g, K_cnt = np.unique(c_subs.values, return_counts=True)
        K_map = dict(zip(all_g, K_cnt))

        def get_fea_dict(rlist):
            N_l = len(rlist)
            if N_l == 0:
                return {}
            hg, kc = np.unique(c_subs.loc[rlist].values, return_counts=True)
            recs_l = [(g, int(k), int(K_map[g]), float(hypergeom.sf(int(k) - 1, M, int(K_map[g]), N_l))) for g, k in zip(hg, kc)]
            recs_l.sort(key=lambda x: x[3])
            fdrs_l = bh_fdr([x[3] for x in recs_l])
            return {g: (k, K, pv, float(fv)) for (g, k, K, pv), fv in zip(recs_l, fdrs_l)}

        fea_main = get_fea_dict(diff_rxns)
        fea_2b = get_fea_dict(list(clean2b))

        fea_rows = [['Pathway', 'P_value_one_sided_hypergeom', 'Adjusted_P_value_BH',
                     'Dysregulated_rxns_in_pathway', 'Total_pathway_rxns_in_model',
                     'MCMC_concordant_rxns_FDR05', 'MCMC_concordant_pct']]

        for g, (k, K, pv, fv) in sorted(fea_main.items(), key=lambda kv: kv[1][2]):
            g_rxns = [r for r in diff_rxns if c_subs.loc[r] == g]
            n_c, n_d33, abs_deltas = 0, 0, []
            for r in g_rxns:
                if r in sub_s.index:
                    sr = sub_s.loc[r]
                    cd, fm = float(sr['cliffs_delta_abs']), float(sr['fdr_abs'])
                    abs_deltas.append(abs(cd))
                    same = (r in up_set and cd > 0) or (r in dn_set and cd < 0)
                    if fm < 0.05 and same:
                        n_c += 1
                        if abs(cd) >= 0.33:
                            n_d33 += 1
            pct_c = 100.0 * n_c / max(k, 1)
            pct_d33 = 100.0 * n_d33 / max(k, 1)
            fea_rows.append([g, pv, fv, k, K, n_c, round(pct_c, 2)])

            if g.startswith('Transport') or g.startswith('Exchange/demand') or g in ('Miscellaneous',):
                continue
            if fv < 0.05:
                in_2b = (g in fea_2b) and (fea_2b[g][3] < 0.05)
                pw_conc_rows.append({
                    'context': ctx, 'pathway': g,
                    'in_clean_2B': in_2b,
                    'in_curated_98': (ctx, g) in ev_map,
                    'evidence_tier': ev_map.get((ctx, g), 'Uncertain'),
                    'k_fsr': k, 'K_model': K, 'fea_fdr': fv,
                    'mcmc_fdr05_concordant_pct': round(pct_c, 1),
                    'mcmc_d33_concordant_pct': round(pct_d33, 1),
                    'fsr_rxns_median_abs_delta': round(float(np.median(abs_deltas)) if abs_deltas else 0.0, 3),
                })
        fea_sheets[ctx] = fea_rows

    df_rxn_c = pd.DataFrame(rxn_conc_rows)
    df_pw_c = pd.DataFrame(pw_conc_rows)
    df_rxn_c.to_csv(os.path.join(OUT_DIR, 'step4_rxn_fsr_mcmc_concordance.csv'), index=False)
    df_pw_c.to_csv(os.path.join(OUT_DIR, 'step4_new_fsr_pathways_with_mcmc_corroboration.csv'), index=False)

    sd = os.path.join(OUT_DIR, 'supp_data')
    save_xlsx(os.path.join(sd, 'Data_S4_Upregulated_reactions.xlsx'), up_sheets)
    save_xlsx(os.path.join(sd, 'Data_S5_Downregulated_reactions.xlsx'), dn_sheets)
    save_xlsx(os.path.join(sd, 'Data_S6_FluxEnrichmentAnalysis.xlsx'), fea_sheets)

    print("\n=== REACTION-LEVEL FSR vs RHMC MCMC CONCORDANCE ===")
    print(df_rxn_c.to_string(index=False))
    print(f"Mean reaction-level MCMC concordance across 10 pairs: {df_rxn_c['mcmc_concordant_pct'].mean():.2f}%")

    print(f"\n=== PATHWAY-LEVEL ENRICHMENT & MCMC CORROBORATION (Total = {len(df_pw_c)}) ===")
    print(df_pw_c.to_string(index=False))
    print(f"Surviving forward-span guard 2B: {int(df_pw_c['in_clean_2B'].sum())}/{len(df_pw_c)}")
    print(f"Already in curated 98:           {int(df_pw_c['in_curated_98'].sum())}/{len(df_pw_c)}")
    print(f"Median pathway MCMC concordance: {df_pw_c['mcmc_fdr05_concordant_pct'].median():.1f}% (mean {df_pw_c['mcmc_fdr05_concordant_pct'].mean():.1f}%)")
    print(f"Pathways with >=50% MCMC concordant rxns: {int((df_pw_c['mcmc_fdr05_concordant_pct'] >= 50.0).sum())}/{len(df_pw_c)}")

    # Render sleek orange Fig. 5d
    cancer_labels = {
        'Brain': 'Glioblastoma', 'Kidney': 'Renal Clear Cell\nCarcinoma',
        'Liver': 'Liver Cancer\n(Hepatoblastoma)', 'Lung': 'Lung Adenocarcinoma',
        'Eye': 'Retinoblastoma', 'Pancreas': 'Adenosquamous Carcinoma\n(Pancreas)',
        'Prostate': 'Prostate Small\nCell Carcinoma', 'Breast': 'Invasive Breast\nCarcinoma',
        'Ovary': 'High-Grade Serous\nOvarian Cancer', 'Biliary Tract': 'Gallbladder\nCancer',
    }
    df_p = df_pw_c.copy()
    df_p['cancer_label'] = df_p['context'].map(cancer_labels)
    df_p['neglog10_padj'] = -np.log10(np.maximum(df_p['fea_fdr'].to_numpy(float), 1e-120))
    df_p['mcmc_corroborated'] = df_p['mcmc_fdr05_concordant_pct'] >= 50.0

    cancer_order = [cancer_labels[c] for c in CONTEXTS]
    pw_counts = df_p.groupby('pathway')['context'].nunique()
    pw_min_p = df_p.groupby('pathway')['fea_fdr'].min()
    pw_order = sorted(df_p['pathway'].unique(), key=lambda p: (-pw_counts[p], pw_min_p[p]))

    df_p['cancer_label'] = pd.Categorical(df_p['cancer_label'], categories=cancer_order, ordered=True)
    df_p['pathway'] = pd.Categorical(df_p['pathway'], categories=pw_order, ordered=True)
    df_p = df_p.sort_values(['cancer_label', 'pathway']).reset_index(drop=True)

    palette = {'Pathway-level': '#076d6d', 'Component-level': '#08306b', 'Uncertain': '#91cf60'}
    evidence_order = ['Pathway-level', 'Component-level', 'Uncertain']
    ring_col = '#E66100'
    s_min, s_max = 170.0, 700.0
    v_min, v_max = df_p['neglog10_padj'].min(), df_p['neglog10_padj'].max()

    sns.set_theme(style='white', context='paper')
    fig, ax = plt.subplots(figsize=(max(14.5, len(pw_order) * 0.65), 8.5))
    ax = sns.scatterplot(
        data=df_p, x='pathway', y='cancer_label', size='neglog10_padj',
        hue='evidence_tier', hue_order=evidence_order, sizes=(s_min, s_max),
        palette=palette, edgecolor='white', linewidth=1.1, alpha=0.92, ax=ax, zorder=3,
    )
    for _, r in df_p.iterrows():
        if r['mcmc_corroborated']:
            sz = s_min + (s_max - s_min) * (r['neglog10_padj'] - v_min) / max(v_max - v_min, 1e-9)
            sz_outer = (np.sqrt(sz) + 3.2) ** 2
            ax.scatter([r['pathway']], [r['cancer_label']], s=sz_outer,
                       facecolors='none', edgecolors=ring_col, linewidths=1.2, zorder=5)

    ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.22)
    ax.xaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.22)
    ax.set_axisbelow(True)
    sns.despine(left=True, bottom=True)
    ax.set_xlabel('Enriched Metabolic Pathways', fontsize=20, fontweight='bold', labelpad=10)
    ax.set_ylabel('Cancer Type', fontsize=20, fontweight='bold', labelpad=12)
    plt.xticks(rotation=45, ha='right', fontsize=13.5)
    plt.yticks(fontsize=13.5)

    handles, labels = ax.get_legend_handles_labels()
    labels = ['Evidence level' if l == 'evidence_tier' else '$-\\log_{10}(p_{\\mathrm{adj}})$' if l == 'neglog10_padj' else l for l in labels]
    handles.append(Line2D([0], [0], marker='o', color='w', label='MCMC corroborated ($\\geq 50\\%$)',
                          markerfacecolor='white', markeredgecolor=ring_col, markeredgewidth=1.4, markersize=13))
    labels.append('MCMC corroborated ($\\geq 50\\%$)')
    leg = ax.legend(handles, labels, bbox_to_anchor=(1.02, 1), loc='upper left',
                    borderaxespad=0, frameon=True, fontsize=12.5, title='Evidence & Significance',
                    title_fontsize=12.5, labelspacing=1.1, handletextpad=0.8)
    for h in leg.legend_handles[1:4]:
        h.set_markersize(12.5)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'bubble_fea_orange_sleek_new.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'bubble_fea_orange_sleek_new.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved updated Fig. 5d bubble_fea_orange_sleek_new.pdf!", flush=True)


if __name__ == '__main__':
    sampl_df, _ = merge_mcmc_all_pairs()
    run_corroboration_and_export(sampl_df)
