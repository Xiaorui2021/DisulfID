# -*- coding: utf-8 -*-
"""
DisulfID v1.0 - Disulfide Linkage Identification System (5-step integrated)

Workflow: input MGF file + peptide sequence -> one-click pipeline:
  1) Combination/species database (exactly 4 Cys, 2 disulfide bonds)
  2) Theoretical b/y fragment library (whole-molecule fragmentation, linking-ion
     model, loop-internal suppression incl. Cys itself, fragment charge-state
     m/z window 200-1200 by default, GUI-adjustable)
  3) MGF parsing (MS1 RT/mz/charge/Area/M + MS2 fragments)
  4) Library search & scoring (precursor n=0/1/2 rule + dual-channel fragments
     + key-ion weighting + dedup)
  5) Judgment of disulfide linkage (beads / globular / ribbon /
     linear) + report with key evidence interpretation
Output: GUI table of identification per chromatographic peak + judgment;
      exports disulfide_judgment.csv / disulfide_report.md /
      disulfide_report.html and ms2_annotated/ annotated spectra.

Dependencies: PyQt6 (GUI), matplotlib (annotated spectra, optional)
Usage:
  GUI:        python disulfide_analyzer.py
  CLI:        python disulfide_analyzer.py --cli <mgf file> [sequence] [out dir]
"""

import bisect
import csv
import os
import re
import sys

try:
    sys.stdout.reconfigure(errors='replace')
except Exception:
    pass

# ===================== Core constants =====================
AA_MONO = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}
H2O = 18.010565
PROTON = 1.007276466812
DISULFIDE_LOSS_PER_BOND = -2 * 1.007825
DEFAULT_MODIFICATION = -0.98
MZ_MIN, MZ_MAX = 200.0, 1200.0   # fragment charge-state m/z window (GUI-adjustable)
BRIDGES = {
    'beads': [(0, 1), (2, 3)],
    'globular': [(0, 2), (1, 3)],
    'ribbon': [(0, 3), (1, 2)],
}
BRIDGE_INFO = {
    'beads': ('beads', 'C1-C2, C3-C4'),
    'globular': ('globular', 'C1-C3, C2-C4'),
    'ribbon': ('ribbon', 'C1-C4, C2-C3'),
    'linear': ('Linear peptide', 'no disulfide'),
}

def bridge_disp(b):
    """Internal bridge value -> display string."""
    return 'linear' if str(b) == 'linear' else str(b)
DEFAULT_SEQ = "ACQEFCTYERDCASFDFDPSSNKC"

# ===================== Step 1-2: Theoretical fragment library =====================
def calc_mass(seq):
    return sum(AA_MONO[aa] for aa in seq) + H2O

def charge_states(mass, lo=MZ_MIN, hi=MZ_MAX):
    states = []
    z = 1
    while True:
        mz = (mass + z * PROTON) / z
        if mz < lo:
            break
        if mz <= hi:
            states.append((z, round(mz, 4)))
        z += 1
    return states

def parse_sequence(raw):
    seq, mod = raw, 0.0
    for op, cl in (('（', '）'), ('(', ')')):
        if op in seq:
            head, rest = seq.split(op, 1)
            tag = rest.split(cl, 1)[0].strip()
            try:
                mod = float(tag)
            except ValueError:
                mod = 0.0
            seq = head
            break
    return seq, mod

def get_c_fragments(seq):
    cuts = [0]
    for i, aa in enumerate(seq):
        if aa in ('R', 'K'):
            cuts.append(i + 1)
    cuts.append(len(seq))
    frags = []
    for a in range(len(cuts)):
        for b in range(a + 1, len(cuts)):
            s, e = cuts[a], cuts[b]
            if 'C' in seq[s:e]:
                frags.append((seq[s:e], s, e))
    return frags

def is_non_overlapping(combo):
    for i in range(len(combo)):
        for j in range(i + 1, len(combo)):
            s1, e1 = combo[i][1], combo[i][2]
            s2, e2 = combo[j][1], combo[j][2]
            if not (e1 <= s2 or e2 <= s1):
                return False
    return True

def get_free_fragments(seq):
    cuts = [0]
    for i, aa in enumerate(seq):
        if aa in ('R', 'K'):
            cuts.append(i + 1)
    cuts.append(len(seq))
    seen, out = set(), []
    for a in range(len(cuts)):
        for b in range(a + 1, len(cuts)):
            fseq = seq[cuts[a]:cuts[b]]
            if 'C' in fseq or fseq in seen:
                continue
            seen.add(fseq)
            out.append(fseq)
    return out

def make_species(combo, bridge, seq):
    c_pos = [i for i, aa in enumerate(seq) if aa == 'C']
    bonds_abs = [(c_pos[a], c_pos[b]) for a, b in bridge]
    frag_of_pos = {}
    for fid, (fseq, s, e) in enumerate(combo):
        for p in c_pos:
            if s <= p < e:
                frag_of_pos[p] = fid
    parent = list(range(len(combo)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx
    for pa, pb in bonds_abs:
        union(frag_of_pos[pa], frag_of_pos[pb])
    groups = {}
    for fid in range(len(combo)):
        groups.setdefault(find(fid), []).append(fid)
    species_list = []
    for root, fids in groups.items():
        chains = []
        for fid in sorted(fids):
            fseq, s, e = combo[fid]
            local_cys = sorted(p - s for p in c_pos if s <= p < e)
            chains.append({'fid': fid, 'seq': fseq, 's': s, 'e': e,
                           'local_cys': local_cys})
        fset = set(fids)
        bonds = []
        for pa, pb in bonds_abs:
            fa, fb = frag_of_pos[pa], frag_of_pos[pb]
            if fa in fset and fb in fset:
                bonds.append((fa, fb, pa - combo[fa][1], pb - combo[fb][1]))
        species_list.append({'chains': chains, 'bonds': bonds})
    return species_list

def bond_map_display(seq, bset):
    out, c_no = [], 0
    for r, aa in enumerate(seq):
        if aa == 'C':
            c_no += 1
            label = 'C%d' % c_no
        else:
            label = aa
        out.append('%s%s' % (label, '✗' if r in bset else '✓'))
    return ''.join(out)

def species_ions(sp, combo, ctx):
    """Fragment a single species as whole molecule -> [(fid,type,n,seq,attach,nbond,mass)] + skipped + maps."""
    rows, skipped, maps = [], [], []
    chain_mass = {ch['fid']: calc_mass(ch['seq']) for ch in sp['chains']}
    nb_total = len(sp['bonds'])
    species_m = round(sum(chain_mass.values())
                      + nb_total * DISULFIDE_LOSS_PER_BOND, 4)
    intra, inter = {}, {}
    for fa, fb, la, lb in sp['bonds']:
        if fa == fb:
            intra.setdefault(fa, []).append((la, lb))
        else:
            inter.setdefault(fa, []).append((fb, la))
            inter.setdefault(fb, []).append((fa, lb))
    blocked = {}
    for fid, bonds in intra.items():
        bset = set()
        for p, q in bonds:
            bset.update(range(p, q))       # loop-internal suppression: includes the Cys itself
        blocked[fid] = bset
    for fid, bset in blocked.items():
        ch = next(c for c in sp['chains'] if c['fid'] == fid)
        maps.append((fid, bond_map_display(ch['seq'], bset)))
    def carried(fid, cys_in):
        seen, frontier = {fid}, [fid]
        while frontier:
            cur = frontier.pop()
            for other, pos in inter.get(cur, ()):
                if cur == fid and pos not in cys_in:
                    continue
                if other not in seen:
                    seen.add(other)
                    frontier.append(other)
        return seen
    for ch in sp['chains']:
        fid, seq, L = ch['fid'], ch['seq'], len(ch['seq'])
        cset = set(ch['local_cys'])
        bset = blocked.get(fid, set())
        for i in range(L - 1):
            if i in bset:
                skipped.append((fid, i, 'loop-internal suppressed'))
                continue
            pre = {r for r in cset if r <= i}
            suf = {r for r in cset if r > i}
            cb = carried(fid, pre)
            cy = carried(fid, suf)
            if len(cb & cy) > 1:
                skipped.append((fid, i, 'macro-loop suppressed'))
                continue
            def bonds_on(side_chains, fid_in):
                n = 0
                for fa, fb, la, lb in sp['bonds']:
                    a_in = (fa == fid and fid_in(la)) or (fa != fid and fa in side_chains)
                    b_in = (fb == fid and fid_in(lb)) or (fb != fid and fb in side_chains)
                    if a_in and b_in:
                        n += 1
                return n
            nb_b = bonds_on(cb, lambda x: x <= i)
            nb_y = bonds_on(cy, lambda x: x > i)
            b_mass = round(sum(AA_MONO[a] for a in seq[:i + 1])
                           + sum(chain_mass[c] for c in cb if c != fid)
                           + nb_b * DISULFIDE_LOSS_PER_BOND, 4)
            y_mass = round(sum(AA_MONO[a] for a in seq[i + 1:]) + H2O
                           + sum(chain_mass[c] for c in cy if c != fid)
                           + nb_y * DISULFIDE_LOSS_PER_BOND, 4)
            attach_b = '+'.join(str(combo[c][0]) for c in sorted(cb - {fid}))
            attach_y = '+'.join(str(combo[c][0]) for c in sorted(cy - {fid}))
            rows.append((fid, 'b', i + 1, seq[:i + 1], attach_b, nb_b, b_mass))
            rows.append((fid, 'y', L - (i + 1), seq[i + 1:], attach_y,
                         nb_y, y_mass))
            assert abs(b_mass + y_mass - species_m) < 0.01, \
                ('Mass conservation failed: %.4f + %.4f != %.4f (%s)'
                 % (b_mass, y_mass, species_m, ctx))
    return rows, skipped, maps

def combinations2(c_frags):
    """Equivalent to itertools.combinations for k=1..4."""
    from itertools import combinations
    for k in range(1, 5):
        for combo in combinations(c_frags, k):
            yield combo

def build_library(raw_seq, msgs=None, mz_min=MZ_MIN, mz_max=MZ_MAX):
    """Step 1-2: sequence -> species table + ion index. Returns (species_rows, ions_index, combos_info)."""
    if msgs is None:
        msgs = []
    seq, mod = parse_sequence(raw_seq)
    c_pos = [i for i, aa in enumerate(seq) if aa == 'C']
    if len(c_pos) != 4:
        raise ValueError('Sequence contains %d Cys; this software handles '
                         'exactly 4 Cys (2 disulfide bonds)' % len(c_pos))
    for aa in seq:
        if aa not in AA_MONO:
            raise ValueError('Non-standard amino acid: %s' % aa)
    c_frags = get_c_fragments(seq)
    combos = []
    for combo in combinations2(c_frags):
        if not is_non_overlapping(combo):
            continue
        if sum(f[0].count('C') for f in combo) == 4:
            combos.append(combo)
    if not combos:
        raise ValueError('No valid digestion combination (non-overlapping with '
                         'total C = 4)')
    msgs.append('Sequence %s | length %d aa | Cys positions(0-based) %s | '
                'modification %g' % (seq, len(seq), c_pos, mod))
    msgs.append('Valid digestion combos %d x 3 bridge hypotheses; fragment '
                'charge-state m/z range %.0f-%.0f'
                % (len(combos), mz_min, mz_max))
    species_rows, ions_index = [], {}
    for ci, combo in enumerate(combos, 1):
        combo_mass = round(sum(calc_mass(f[0]) for f in combo)
                           + 2 * DISULFIDE_LOSS_PER_BOND + mod, 4)
        for bridge_name, bridge in BRIDGES.items():
            sl = make_species(combo, bridge, seq)
            for si, sp in enumerate(sl, 1):
                nb = len(sp['bonds'])
                sp_mass = round(sum(calc_mass(ch['seq']) for ch in sp['chains'])
                                + nb * DISULFIDE_LOSS_PER_BOND, 4)
                chains_str = ' + '.join(ch['seq'] for ch in sp['chains'])
                def _clabel(pos):
                    return 'C%d' % (c_pos.index(pos) + 1)
                bonds_str = '; '.join(
                    '%s-%s' % (_clabel(combo[_fa][1] + _la),
                               _clabel(combo[_fb][1] + _lb))
                    for _fa, _fb, _la, _lb in sp['bonds'])
                ion_rows, skipped, maps = species_ions(sp, combo,
                                                       'combo%d %s' % (ci, bridge_name))
                entry = {'ions': {}}
                for fid, itype, n, iseq, attach, nb_i, mass in ion_rows:
                    ion_id = (ci, bridge_name, si, fid, itype, iseq,
                              attach or '-')
                    if ion_id not in entry['ions']:
                        entry['ions'][ion_id] = {
                            'id': ion_id, 'fid': fid, 'type': itype,
                            'seq': iseq, 'carry': attach or '-',
                            'nbond': nb_i, 'M': mass, 'mz': []}
                    for z, mz in charge_states(mass, mz_min, mz_max):
                        entry['ions'][ion_id]['mz'].append((z, mz))
                entry['total'] = len(entry['ions'])
                entry['total_w'] = sum(ion_weight(entry['ions'][i])
                                       for i in entry['ions'])
                entry['by_mass'] = {}
                for ion in entry['ions'].values():
                    entry['by_mass'].setdefault(ion['M'], []).append(ion)
                entry['sorted_masses'] = sorted(entry['by_mass'].keys())
                ions_index[(ci, bridge_name, si)] = entry
                species_rows.append({'combo': ci, 'bridge': bridge_name,
                                     'sp': si, 'compo': chains_str,
                                     'bonds': bonds_str, 'nb': nb,
                                     'M': sp_mass,
                                     'keymap': '; '.join('chain%d: %s' % (f, m)
                                                         for f, m in maps) or '-'})
    free_frags = get_free_fragments(seq)
    if free_frags:
        msgs.append('Free single-chain (non-disulfide) fragments: %s'
                    % ', '.join(free_frags))
        for fi, fseq in enumerate(free_frags, 1):
            fm = round(calc_mass(fseq), 4)
            species_rows.append({'combo': 0, 'bridge': 'linear', 'sp': fi,
                                 'compo': fseq, 'bonds': '-', 'nb': 0,
                                 'M': fm, 'keymap': '-'})
            entry = {'ions': {}}
            for i in range(len(fseq) - 1):
                b_m = round(sum(AA_MONO[a] for a in fseq[:i + 1]), 4)
                y_m = round(sum(AA_MONO[a] for a in fseq[i + 1:]) + H2O, 4)
                for itype, iseq, mass in (('b', fseq[:i + 1], b_m),
                                          ('y', fseq[i + 1:], y_m)):
                    ion_id = (0, 'linear', fi, 0, itype, iseq, '-')
                    if ion_id not in entry['ions']:
                        entry['ions'][ion_id] = {
                            'id': ion_id, 'fid': 0, 'type': itype,
                            'seq': iseq, 'carry': '-', 'nbond': 0,
                            'M': mass, 'mz': []}
                    for z, mz in charge_states(mass, mz_min, mz_max):
                        entry['ions'][ion_id]['mz'].append((z, mz))
            entry['total'] = len(entry['ions'])
            entry['total_w'] = sum(1.0 for _ in entry['ions'])
            entry['by_mass'] = {}
            for ion in entry['ions'].values():
                entry['by_mass'].setdefault(ion['M'], []).append(ion)
            entry['sorted_masses'] = sorted(entry['by_mass'].keys())
            ions_index[(0, 'linear', fi)] = entry
    msgs.append('Theoretical library: species %d, ion groups '
                '(combo x bridge x sp) %d, total ions %d'
                % (len(species_rows), len(ions_index),
                   sum(e['total'] for e in ions_index.values())))
    return species_rows, ions_index, combos

# ===================== Step 3: MGF parsing =====================
def infer_frag_charge(peaks, ppm_lim):
    n = len(peaks)
    zs = [1] * n
    for i in range(n):
        Mi = peaks[i][0]
        cutoff = Mi * ppm_lim * 1e-6
        best_d = None
        for j in range(i + 1, n):
            d = peaks[j][0] - Mi
            if d >= 1.02:
                break
            if best_d is None or d < best_d:
                best_d = d
        if best_d is None:
            continue
        if abs(best_d - 1.0073) <= cutoff:
            zs[i] = 1
        elif abs(best_d - 0.5036) <= cutoff:
            zs[i] = 2
        elif abs(best_d - 0.3357) <= cutoff:
            zs[i] = 3
    return zs

def parse_mgf_file(path, params, msgs=None):
    """Step 3: parse MGF -> (specs, frag_map)
       specs: {title, rt_min, mz, charge, scan, area, M}
       frag_map: key(scan or title) -> [{mz,intensity,z,neutral}]"""
    if msgs is None:
        msgs = []
    a = 1.0073
    ppm_limit = params['ppm_charge'] if 'ppm_charge' in params else 5
    cutoff = params['cutoff']
    ms1_filter = params['ms1_filter']
    sp = None
    spectra = []
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line == 'BEGIN IONS':
                sp = {'title': '', 'rt_sec': None, 'pepmass': None,
                      'area': None, 'scan': None, 'charge': None, 'sign': '+',
                      'peaks': []}
                continue
            if sp is None:
                continue
            if line == 'END IONS':
                spectra.append(sp)
                sp = None
                continue
            if line.startswith('TITLE='):
                sp['title'] = line.split('=', 1)[1].strip()
            elif line.startswith('RTINSECONDS='):
                try:
                    sp['rt_sec'] = float(line.split('=', 1)[1].strip())
                except ValueError:
                    sp['rt_sec'] = None
            elif line.startswith('PEPMASS='):
                parts = line.split('=', 1)[1].split()
                try:
                    sp['pepmass'] = float(parts[0])
                except (ValueError, IndexError):
                    sp['pepmass'] = None
                if len(parts) > 1:
                    try:
                        sp['area'] = float(parts[1])
                    except ValueError:
                        sp['area'] = None
            elif line.startswith('SCANS='):
                m_ = re.search(r'\d+', line.split('=', 1)[1])
                if m_:
                    sp['scan'] = int(m_.group(0))
            elif line.startswith('CHARGE='):
                val = line.split('=', 1)[1].strip()
                m_ = re.search(r'\d+', val)
                if m_:
                    sp['charge'] = int(m_.group(0))
                if val.endswith('+'):
                    sp['sign'] = '+'
                elif val.endswith('-'):
                    sp['sign'] = '-'
            else:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        inten = float(parts[1])
                        if inten > cutoff:
                            sp['peaks'].append((mz, inten))
                    except ValueError:
                        pass
    specs, frag_map = [], {}
    n_skip_area = 0
    for sp in spectra:
        if sp['pepmass'] is None or sp['area'] is None:
            n_skip_area += 1
            continue
        if sp['area'] < ms1_filter:
            n_skip_area += 1
            continue
        n = float(sp['charge'] or 0.0)
        if sp['scan'] is None:
            m_ = re.search(r'scan[=: ](\d+)', sp['title'], re.I)
            sp['scan'] = int(m_.group(1)) if m_ else ''
        rt_min = round((sp['rt_sec'] or 0.0) / 60, 2)
        m = round(sp['pepmass'], 4)
        M = round(m * n - a * n, 4)
        peaks = sorted(sp['peaks'], key=lambda x: x[0])
        zs = infer_frag_charge(peaks, ppm_limit)
        frs = []
        for (mz, inten), z in zip(peaks, zs):
            mi = round(mz * z - a * (z - 1), 4)
            frs.append({'mz': round(mz, 4), 'intensity': inten,
                        'z': float(z), 'neutral': mi})
        key = str(sp['scan']) if sp['scan'] else (sp['title'] + '|')
        frag_map[key] = frs
        specs.append({'title': sp['title'], 'rt_min': rt_min, 'mz': m,
                      'charge': '%g%s' % (n, sp['sign']), 'scan': sp['scan'],
                      'area': sp['area'], 'M': M})
    msgs.append('MGF parsed: spectra %d (Area>=%.0f, skipped %d no/weak '
                'intensity) | MS2 peaks %d'
                % (len(specs), ms1_filter, n_skip_area,
                   sum(len(v) for v in frag_map.values())))
    return specs, frag_map

# ===================== Step 4: Library search and scoring =====================
KEY_ION_WEIGHT, LOOP_ION_WEIGHT, FREE_ION_WEIGHT = 3.0, 2.0, 1.0
N_OFFSETS = (0, 1, 2)
DEDUP_RT_TOL, DEDUP_MASS_TOL = 1.0, 0.02

def ion_weight(ion):
    if ion['carry'] and ion['carry'] != '-':
        # linking ions: both b and y ions are topology-discriminating evidence
        return KEY_ION_WEIGHT
    if ion['nbond'] and ion['nbond'] > 0:
        return LOOP_ION_WEIGHT
    return FREE_ION_WEIGHT

def precursor_match(theo_m, exp_m, ppm, offsets):
    if theo_m is None or exp_m is None:
        return True, '-'
    tol = exp_m * ppm * 1e-6
    diff = abs(theo_m - exp_m)
    for n in offsets:
        if abs(diff - n) <= tol:
            return True, n
    return False, None

def score_candidate(spec_frags, ions_idx, ppm2):
    total_theo = ions_idx['total']
    total_w = ions_idx['total_w']
    by_mass, sorted_m = ions_idx['by_mass'], ions_idx['sorted_masses']
    matched_ion, matched_frag = set(), set()
    matched_int = 0.0
    details = []
    for fi, fr in enumerate(spec_frags):
        hit, channel = None, None
        if fr['neutral'] is not None:
            tol = fr['neutral'] * ppm2 * 1e-6
            i0 = bisect.bisect_left(sorted_m, fr['neutral'] - tol)
            i1 = bisect.bisect_right(sorted_m, fr['neutral'] + tol)
            for k in range(i0, i1):
                for ion in by_mass[sorted_m[k]]:
                    if ion['id'] not in matched_ion:
                        hit, channel = ion, 'neutral'
                        break
                if hit:
                    break
        if hit is None and fr['mz'] is not None:
            tol = fr['mz'] * ppm2 * 1e-6
            for ion in ions_idx['ions'].values():
                if ion['id'] in matched_ion:
                    continue
                for z, mz in ion['mz']:
                    if abs(mz - fr['mz']) <= tol:
                        hit, channel = ion, 'mz'
                        break
                if hit:
                    break
        if hit is not None:
            matched_ion.add(hit['id'])
            matched_frag.add(fi)
            matched_int += fr['intensity'] or 0.0
            if channel == 'neutral':
                dd = abs(hit['M'] - fr['neutral'])
                ppm = dd / fr['neutral'] * 1e6 if fr['neutral'] else 0.0
            else:
                mz_hit = min(hit['mz'], key=lambda x: abs(x[1] - fr['mz']))
                dd = abs(mz_hit[1] - fr['mz'])
                ppm = dd / fr['mz'] * 1e6 if fr['mz'] else 0.0
            details.append({'fi': fi, 'fr': fr, 'ion': hit,
                            'channel': channel, 'dd': round(dd, 4),
                            'ppm': round(ppm, 2), 'weight': ion_weight(hit)})
    matched_w = 0.0
    matched_key = 0
    for ion_id in matched_ion:
        w = ion_weight(ions_idx['ions'][ion_id])
        matched_w += w
        if w > FREE_ION_WEIGHT:
            matched_key += 1
    total_int = sum((fr['intensity'] or 0.0) for fr in spec_frags)
    cov_ion_w = matched_w / total_w if total_w > 0 else 0.0
    cov_ion_raw = len(matched_ion) / total_theo if total_theo else 0.0
    cov_int = matched_int / total_int if total_int > 0 else 0.0
    stat = {'n_theo': total_theo, 'n_match_ion': len(matched_ion),
            'n_match_key': matched_key, 'cov_ion': round(cov_ion_w, 4),
            'cov_ion_raw': round(cov_ion_raw, 4), 'cov_int': round(cov_int, 4),
            'score': round(100 * (0.7 * cov_ion_w + 0.3 * cov_int), 2)}
    return stat, details

def dedup_results(rows, rt_tol, mass_tol):
    cands = [r for r in rows if r['n_offset'] != '-']
    others = [r for r in rows if r['n_offset'] == '-']
    keep, used = [], [False] * len(cands)
    # Dedup cluster representative: rows that CAN discriminate connectivity
    # (n_disc_ev > 0, i.e. loop-ion / single-Cys y linking evidence exists)
    # outrank degenerate rows (n_disc_ev = 0). Within the same group (both can
    # discriminate, or both cannot), prefer the higher-scoring row, then
    # smaller scan id. This keeps e.g. scan 1102 (ribbon, explicit loop-ion
    # evidence) over a same-APEX scan with higher score but degenerate globular/ribbon,
    # while among spectra that can both discriminate it picks the best score.
    order = sorted(range(len(cands)),
                   key=lambda i: (0 if (fnum(cands[i].get('n_disc_ev'), 0.0)
                                        or 0.0) > 0 else 1,
                                  -(fnum(cands[i]['score'], 0.0) or 0.0),
                                  str(cands[i]['scan'])))
    for i in order:
        if used[i]:
            continue
        ri = cands[i]
        rt_i = fnum(ri['RT_min'], 0.0) or 0.0
        m_i = fnum(ri['exp_neutral_M'], 0.0) or 0.0
        cluster = [i]
        for j in order:
            if used[j] or j == i:
                continue
            rj = cands[j]
            rt_j = fnum(rj['RT_min'], 0.0) or 0.0
            m_j = fnum(rj['exp_neutral_M'], 0.0) or 0.0
            if abs(rt_i - rt_j) <= rt_tol and abs(m_i - m_j) <= mass_tol:
                cluster.append(j)
        for j in cluster:
            used[j] = True
        keep.append(ri)
    keep.extend(others)
    return keep

def fnum(x, dflt=None):
    try:
        return float(x)
    except (ValueError, TypeError):
        return dflt

def inum(x, dflt=None):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return dflt

def run_search(species_rows, ions_index, specs, frag_map, params, msgs=None):
    """Step 4: per-spectrum scoring + dedup -> (rows_results, rows_matches, spec_store, stats)."""
    if msgs is None:
        msgs = []
    ppm1 = params['ppm1']
    ppm2 = params['ppm2']
    rows_results, rows_matches = [], []
    spec_store = {}
    tot_spec = tot_cand = tot_hit = 0
    for sp in specs:
        key = str(sp['scan']) if sp['scan'] else (sp['title'] + '|')
        spec_frags = frag_map.get(key, [])
        cands = []
        for spc in species_rows:
            ok, n_off = precursor_match(spc['M'], sp['M'], ppm1, N_OFFSETS)
            if ok:
                cands.append((spc, n_off))
        if not cands:
            rows_results.append({
                'scan': sp['scan'], 'TITLE': sp['title'],
                'RT_min': sp['rt_min'] if sp['rt_min'] is not None else '',
                'precursor_mz': sp['mz'] if sp['mz'] is not None else '',
                'precursor_charge': sp['charge'],
                'precursor_Area': sp['area'] if sp['area'] is not None else '',
                'exp_neutral_M': sp['M'] if sp['M'] is not None else '',
                'n_offset': '-', 'combo_id': '-', 'bridge': '-', 'species': '-',
                'composition': '-', 'ss_bonds': '-', 'theo_species_mass': '-',
                'Δmass_Da': '-', 'Δppm': '-', 'n_theo_ions': 0,
                'n_matched_ion': 0, 'n_matched_key': 0,
                'cov_ion': 0, 'cov_ion_raw': 0, 'intensity_coverage': 0,
                'score': 0, 'n_disc_ev': 0,
                'conclusion': 'No candidate (precursor mass failed n=0/1/2 rule)'})
            tot_spec += 1
            continue
        scored = []
        for spc, n_off in cands:
            key_t = (spc['combo'], spc['bridge'], spc['sp'])
            if key_t not in ions_index:
                continue
            stat, details = score_candidate(spec_frags, ions_index[key_t], ppm2)
            diff = abs(spc['M'] - sp['M']) if (spc['M'] is not None
                                               and sp['M'] is not None) else None
            scored.append((stat, details, spc, n_off, diff))
        scored.sort(key=lambda x: x[0]['score'], reverse=True)
        for rank, (stat, details, spc, n_off, diff) in enumerate(scored, 1):
            dppm = (diff / sp['M'] * 1e6) if (diff is not None and sp['M']) else None
            # Count b-/y-ion pairing discrimination evidence for this candidate
            # (loop ions closing intra-chain disulfides + single-Cys b-/y-linking
            # ions) so that dedup can prefer spectra with explicit connectivity
            # evidence over higher-scoring but degenerate ones.
            compo_parts_i = [p.strip() for p in spc['compo'].split('+')]
            disc_ev = []
            for dt in details:
                ion, fr = dt['ion'], dt['fr']
                fid = ion['fid']
                cseq = (compo_parts_i[fid] if fid < len(compo_parts_i) else '')
                if (str(ion.get('carry') or '-') != '-'
                        or inum(ion.get('nbond'), 0) > 0):
                    mdict = {'chain_idx': ion['fid'], 'chain_seq': cseq,
                             'ion_type': ion['type'], 'ion_seq': ion['seq'],
                             'carries': ion['carry'],
                             'nbond': ion['nbond']}
                    ok, why = check_ion_consistency(mdict, compo_parts_i,
                                                    parse_bonds(spc['bonds']))
                    disc_ev.append((mdict, ok, why))
            _, _, disc_notes = infer_connectivity_from_evidence(
                disc_ev, compo_parts_i)
            n_disc_ev = len(disc_notes)
            store_key = (str(sp['scan']), spc['combo'], spc['bridge'], spc['sp'])
            spec_store.setdefault(store_key, {'frags': spec_frags,
                                              'details': None})
            spec_store[store_key]['details'] = details
            rows_results.append({
                'scan': sp['scan'], 'TITLE': sp['title'],
                'RT_min': sp['rt_min'] if sp['rt_min'] is not None else '',
                'precursor_mz': sp['mz'] if sp['mz'] is not None else '',
                'precursor_charge': sp['charge'],
                'precursor_Area': sp['area'] if sp['area'] is not None else '',
                'exp_neutral_M': sp['M'] if sp['M'] is not None else '',
                'n_offset': n_off, 'combo_id': spc['combo'],
                'bridge': spc['bridge'], 'species': spc['sp'],
                'composition': spc['compo'], 'ss_bonds': spc['bonds'],
                'theo_species_mass': spc['M'] if spc['M'] is not None else '',
                'Δmass_Da': ('%.4f' % diff) if diff is not None else '-',
                'Δppm': ('%.2f' % dppm) if dppm is not None else '-',
                'n_theo_ions': stat['n_theo'],
                'n_matched_ion': stat['n_match_ion'],
                'n_matched_key': stat['n_match_key'],
                'cov_ion': stat['cov_ion'],
                'cov_ion_raw': stat['cov_ion_raw'],
                'intensity_coverage': stat['cov_int'],
                'score': stat['score'], 'n_disc_ev': n_disc_ev,
                'conclusion': 'Candidate (rank %d)' % rank})
            compo_parts = [p.strip() for p in spc['compo'].split('+')]
            for dt in details:
                ion, fr = dt['ion'], dt['fr']
                fid = ion['fid']
                chain_seq = compo_parts[fid] if fid < len(compo_parts) else ''
                rows_matches.append({
                    'scan': sp['scan'], 'TITLE': sp['title'],
                    'combo_id': spc['combo'], 'bridge': spc['bridge'],
                    'species': spc['sp'], 'chain_idx': ion['fid'],
                    'chain_seq': chain_seq,
                    'ion_type': ion['type'], 'ion_seq': ion['seq'],
                    'carries': ion['carry'],
                    'nbond': ion['nbond'],
                    'theo_neutral_mass': ion['M'],
                    'match_channel': dt['channel'],
                    'weight': dt['weight'],
                    'exp_frag_mz': fr['mz'] if fr['mz'] is not None else '',
                    'exp_charge_z': fr['z'] if fr['z'] is not None else '-',
                    'exp_intensity': fr['intensity'],
                    'Δ_Da': dt['dd'], 'Δ_ppm': dt['ppm']})
        tot_spec += 1
        tot_cand += len(scored)
        tot_hit += sum(1 for s, _, _, _, _ in scored if s['n_match_ion'] > 0)
        if scored:
            best = scored[0][0]
            msgs.append('Spectrum scan=%s M=%.4f z=%s | candidates %d | best '
                        'combo%d-%s-sp%d score %.2f key ions %d/%d n_offset %s'
                        % (sp['scan'], sp['M'] or 0.0, sp['charge'],
                           len(scored), scored[0][2]['combo'],
                           scored[0][2]['bridge'], scored[0][2]['sp'],
                           best['score'], best['n_match_key'], best['n_theo'],
                           scored[0][3]))
    n_before = len(rows_results)
    rt_tol = params.get('dedup_rt', DEDUP_RT_TOL)
    mass_tol = params.get('dedup_mass', DEDUP_MASS_TOL)
    results_all = list(rows_results)          # every candidate row, no dedup
    matches_all = list(rows_matches)          # every matched ion, no dedup
    rows_results = dedup_results(rows_results, rt_tol, mass_tol)
    kept_keys = set()
    for r in rows_results:
        if r['n_offset'] != '-':
            kept_keys.add((str(r['scan']), inum(r['combo_id']),
                           r['bridge'], inum(r['species'])))
            r['conclusion'] = 'Identified (kept after dedup)'
    rows_matches = [m for m in rows_matches
                    if (str(m['scan']), inum(m['combo_id']), m['bridge'],
                        inum(m['species'])) in kept_keys]
    def sort_key(r):
        s = r['scan']
        try:
            s = int(s)
        except (ValueError, TypeError):
            s = 10 ** 9
        return (s, -(fnum(r['score'], 0.0) or 0.0))
    rows_results.sort(key=sort_key)
    msgs.append('Dedup (RT+/-%.1f min and M+/-%.3f Da): %d -> %d rows | '
                'spectra %d | candidates %d | hits %d'
                % (rt_tol, mass_tol, n_before,
                   len(rows_results), tot_spec, tot_cand, tot_hit))
    return (results_all, matches_all, rows_results, rows_matches,
            spec_store, {'tot_spec': tot_spec, 'tot_cand': tot_cand,
                         'tot_hit': tot_hit})

# ===================== Step 5: Judgment and report =====================
STRONG_SCORE, MID_SCORE = 60.0, 40.0
STRONG_KEY_MIN, MID_KEY_MIN = 3, 1

def parse_bonds(bond_str):
    out = []
    for seg in re.split(r'[;；,，、]', str(bond_str)):
        m_ = re.match(r'C(\d+)\s*[-–]\s*C(\d+)', seg.strip())
        if m_:
            out.append((int(m_.group(1)), int(m_.group(2))))
    return out

def chain_cys_map(compo_parts_list):
    per_chain, global_of = [], []
    for seq in compo_parts_list:
        pos = [p for p, aa in enumerate(seq) if aa == 'C']
        per_chain.append(pos)
        for p in pos:
            global_of.append((len(per_chain) - 1, p))
    return per_chain, global_of

def compo_parts(compo):
    return [p.strip() for p in str(compo).replace('+', '+').split('+') if p.strip()]

def check_ion_consistency(ion, parts, bonds):
    """Verify matched ion consistency with hypothesized disulfides"""
    carry = str(ion.get('carries') or '-')
    nbond = inum(ion.get('nbond'), 0)
    itype = str(ion.get('ion_type') or 'b')
    n = len(str(ion.get('ion_seq') or ''))
    label = '%s%d' % (itype, n)
    chain_seq = str(ion.get('chain_seq') or '')
    if chain_seq in parts:
        fid = parts.index(chain_seq)
    else:
        fid = inum(ion.get('chain_idx'), -1)
        if fid < 0 or fid >= len(parts):
            return True, 'ion %s: chain not located, skipped' % label
    _, global_of = chain_cys_map(parts)
    if carry != '-':
        carried_idxs = []
        for cseq in str(carry).split('+'):
            cseq = cseq.strip()
            if cseq in parts:
                carried_idxs.append(parts.index(cseq))
        fid_c = [g + 1 for g, (ci, _p) in enumerate(global_of) if ci == fid]
        ok_pairs, miss = [], []
        for ci2 in carried_idxs:
            c2 = [g + 1 for g, (ci, _p) in enumerate(global_of) if ci == ci2]
            hit = any((a in fid_c and b in c2) or (b in fid_c and a in c2)
                      for a, b in bonds)
            if hit:
                ok_pairs.append('chain%d' % ci2)
            else:
                miss.append('chain%d' % ci2)
        if ok_pairs and not miss:
            return True, ('linking ion %s: chain%d carries %s, consistent with '
                          'hypothesis' % (label, fid, '+'.join(ok_pairs)))
        return False, ('linking ion %s: chain%d carries %s, but no matching '
                       'inter-chain disulfide in hypothesis'
                       % (label, fid, '+'.join(miss) if miss else '?'))
    if nbond and nbond > 0:
        fid_c = [g + 1 for g, (ci, _p) in enumerate(global_of) if ci == fid]
        intra = [b for b in bonds if b[0] in fid_c and b[1] in fid_c]
        if intra:
            return True, ('loop ion %s(-2H): chain%d intra-chain disulfide %s '
                          'consistent'
                          % (label, fid,
                             ', '.join('C%d-C%d' % b for b in intra)))
        return False, ('loop ion %s(-2H): chain%d has no intra-chain disulfide, '
                       'inconsistent' % (label, fid))
    return True, ('free ion %s: no disulfide info, not verified' % label)

def infer_connectivity_from_evidence(evidence, parts):
    """Discriminate beads/globular/ribbon using b- and y-ion pairing evidence.

    Two independent evidence types (either one is sufficient to assign pairing):
      (a) Loop ions: a b or y ion that closes an intra-chain disulfide
          (nbond>0) and whose sequence contains exactly 2 Cys ->
          these 2 Cys form the intra-chain disulfide pair;
      (b) Linking ions: a b or y ion whose sequence contains exactly 1 Cys and which
          carries exactly one single-Cys fragment -> the 2 Cys form the
          inter-chain disulfide pair.
    Pairing map: {C1-C2, C3-C4}->beads, {C1-C3, C2-C4}->globular,
                 {C1-C4, C2-C3}->ribbon; every possible Cys pair maps to
                 exactly one hypothesis, so a single observed pair assigns the
                 connectivity uniquely.
    Returns (concl, pairs, notes); concl in {'beads','globular','ribbon','conflict',None};
    None = no pairing evidence (degenerate).
    """
    per_chain, global_of = chain_cys_map(parts)
    pairs, notes = [], []

    def gnum(ci, p):
        return next((g + 1 for g, (ci2, p2) in enumerate(global_of)
                     if ci2 == ci and p2 == p), None)

    for m, ok, why in evidence:
        if not ok:
            continue
        # Both b and y ions participate in pairing discrimination
        chain_seq = str(m.get('chain_seq') or '')
        ion_seq = str(m.get('ion_seq') or '')
        if chain_seq not in parts or not ion_seq:
            continue
        ci = parts.index(chain_seq)
        idx = chain_seq.find(ion_seq)
        if idx < 0:
            continue
        local_c = [p for p in per_chain[ci] if idx <= p < idx + len(ion_seq)]
        nb = inum(m.get('nbond'), 0)
        carry = str(m.get('carries') or '-')
        if nb > 0 and len(local_c) == 2:
            # loop ion: closed intra-chain disulfide inside the ion
            gs = sorted(g for g in (gnum(ci, p) for p in local_c)
                        if g is not None)
            if len(gs) == 2:
                pair = tuple(gs)
                pairs.append(pair)
                notes.append('loop ion %s%d (%s, chain%d) closes intra-chain '
                             'disulfide C%d-C%d'
                             % (m.get('ion_type','?'), len(ion_seq), ion_seq, ci,
                                pair[0], pair[1]))
        elif carry not in ('-', '') and len(local_c) == 1:
            # linking ion: single Cys in ion + single-Cys carried fragment
            carried = [c.strip() for c in str(carry).split('+')
                       if c.strip() in parts]
            if len(carried) != 1:
                continue
            cj = parts.index(carried[0])
            cy = [g for g in (gnum(cj, p) for p in per_chain[cj])
                  if g is not None]
            if len(cy) != 1:
                continue
            gx = gnum(ci, local_c[0])
            if gx is None:
                continue
            pair = tuple(sorted((gx, cy[0])))
            pairs.append(pair)
            notes.append('linking %s ion %s%d (%s, chain%d) contains C%d and '
                         'carries single-Cys fragment of chain%d (C%d) -> '
                         'inter-chain C%d-C%d'
                         % (m.get('ion_type','?'), m.get('ion_type','?'),
                            len(ion_seq), ion_seq, ci, gx, cj, cy[0],
                            pair[0], pair[1]))
    if not pairs:
        return None, [], notes
    ps = set(pairs)
    hyp_map = {'beads': {(1, 2), (3, 4)},
               'globular': {(1, 3), (2, 4)},
               'ribbon': {(1, 4), (2, 3)}}
    supported = [h for h, hp in hyp_map.items() if ps & hp]
    if len(supported) == 1:
        return supported[0], sorted(ps), notes
    return 'conflict', sorted(ps), notes

def evidence_strength(score, n_key, n_conflict):
    if n_conflict:
        return 'Weak (conflict)'
    if score >= STRONG_SCORE and n_key >= STRONG_KEY_MIN:
        return 'Strong'
    if score >= MID_SCORE and n_key >= MID_KEY_MIN:
        return 'Medium'
    return 'Weak'

def sample_of(row):
    t = str(row.get('TITLE') or '-')
    return t.split('.')[0] if t != '-' else '-'

def judge_and_report(rows_results, rows_matches, msgs=None,
                     rows_all=None, matches_all=None):
    """Step 5: consolidated judgment + report generation.
    Returns (judge_rows, md_text, html_text, sample_notes, judge_all);
    judge_rows = deduped (one row per compound peak); judge_all = all-scan view
    (one row per scan, no RT/M dedup; scans with no candidate still listed)."""
    if msgs is None:
        msgs = []
    match_map = {}
    for m in rows_matches:
        k = (str(m.get('scan')), inum(m.get('combo_id')),
             m.get('bridge'), inum(m.get('species')))
        match_map.setdefault(k, []).append(m)

    def _judge_loop(idents, mmap):
        out, report = [], []
        for i, r in enumerate(idents, 1):
            scan = r.get('scan')
            combo = inum(r.get('combo_id'))
            bridge = r.get('bridge')
            sp = inum(r.get('species'))
            compo = r.get('composition')
            bond_str = r.get('ss_bonds')
            bonds = parse_bonds(bond_str)
            parts = compo_parts(compo)
            score = fnum(r.get('score'), 0.0) or 0.0
            n_key = inum(r.get('n_matched_key'), 0)
            m_key = mmap.get((str(scan), combo, bridge, sp), [])
            key_ions = [m for m in m_key
                        if (str(m.get('carries') or '-') != '-'
                            or inum(m.get('nbond'), 0) > 0)]
            evidence, n_conflict = [], 0
            for m in key_ions:
                ok, why = check_ion_consistency(m, parts, bonds)
                if not ok:
                    n_conflict += 1
                evidence.append((m, ok, why))
            for m in m_key:
                if m not in key_ions:
                    evidence.append((m, True,
                                     'free ion %s%d: no disulfide info'
                                     % (m.get('ion_type'),
                                        len(str(m.get('ion_seq') or '')))))
            # Connectivity discrimination: loop ions (closed intra-chain
            # disulfide, 2 Cys inside) + single-Cys b-/y-linking ions
            disc, disc_pairs, disc_notes = infer_connectivity_from_evidence(
                evidence, parts)
            strength = evidence_strength(score, len(key_ions), n_conflict)
            bname, bpat = BRIDGE_INFO.get(bridge, (bridge, bond_str))
            bond_display = bond_str
            if bridge == 'linear':
                concl = 'Linear peptide (no disulfide)'
            elif bridge in ('beads', 'globular', 'ribbon'):
                if disc is None:
                    if bridge == 'beads':
                        concl = '%s (%s): %s' % (bridge, bname, bpat or bond_str)
                    else:
                        concl = ('globular/ribbon (2 disulfides between chains confirmed; '
                                 'no loop-ion / single-Cys y-ion pairing evidence, '
                                 'crossover/nested indistinguishable)')
                        strength = 'Weak (degenerate)'
                        bond_display = 'globular/ribbon'
                elif disc == 'conflict':
                    concl = 'Connectivity conflict: ' + ', '.join(
                        'C%d-C%d' % p for p in disc_pairs)
                    strength = 'Weak (conflict)'
                    n_conflict += 1
                    bond_display = '; '.join('C%d-C%d' % p for p in disc_pairs)
                else:
                    dname, dpat = BRIDGE_INFO[disc]
                    concl = '%s (%s): %s' % (disc, dname, dpat)
                    bond_display = dpat
                    if disc != bridge:
                        concl += (' [corrected by ion-pairing discrimination; '
                                  'initial %s]' % bridge)
            else:
                concl = '%s (%s): %s' % (bridge, bname, bpat or bond_str)
            if n_conflict:
                concl += ' [evidence conflict, review needed]'
            out.append({'idx': i, 'scan': scan, 'sample': sample_of(r),
                        'RT_min': r.get('RT_min'),
                        'precursor_Area': r.get('precursor_Area'),
                        'exp_neutral_M': r.get('exp_neutral_M'),
                        'combo_id': combo, 'bridge': bridge, 'species': sp,
                        'composition': compo, 'ss_bonds': bond_display,
                        'score': score,
                        'cov_ion': fnum(r.get('cov_ion'), 0.0) or 0.0,
                        'n_matched_ion': inum(r.get('n_matched_ion'), 0),
                        'n_matched_key': n_key, 'n_conflict': n_conflict,
                        'n_disc': len(disc_pairs),
                        'disc_pairs': '; '.join(
                            'C%d-C%d' % p for p in disc_pairs)
                        if disc_pairs else '-',
                        'disc_notes': '; '.join(disc_notes)
                        if disc_notes else '-',
                        'strength': strength, 'conclusion': concl})
            report.append((r, out[-1], evidence))
        return out, report

    def _allscan_judge(rows_all, matches_all):
        """Every scan gets exactly one row (best candidate; no-candidate scans
        included with a 'No candidate' conclusion), no RT/M dedup."""
        mmap_all = {}
        for m in matches_all:
            k = (str(m.get('scan')), inum(m.get('combo_id')),
                 m.get('bridge'), inum(m.get('species')))
            mmap_all.setdefault(k, []).append(m)
        scan_order, seen = [], set()
        best, nocand = {}, {}
        for r in rows_all:
            s = str(r.get('scan'))
            if s not in seen:
                seen.add(s)
                scan_order.append(s)
            if r.get('n_offset') == '-':
                nocand.setdefault(s, r)
            else:
                sc = fnum(r.get('score'), 0.0) or 0.0
                if s not in best or sc > (fnum(best[s].get('score'), 0.0) or 0.0):
                    best[s] = r
        all_judge, _ = _judge_loop([best[s] for s in scan_order if s in best],
                                   mmap_all)
        jmap = {str(j['scan']): j for j in all_judge}
        merged = []
        for s in scan_order:
            if s in jmap:
                merged.append(jmap[s])
            elif s in nocand:
                r = nocand[s]
                merged.append({
                    'idx': 0, 'scan': r.get('scan'),
                    'sample': sample_of(r), 'RT_min': r.get('RT_min'),
                    'precursor_Area': r.get('precursor_Area'),
                    'exp_neutral_M': r.get('exp_neutral_M'),
                    'combo_id': '-', 'bridge': '-', 'species': '-',
                    'composition': '-', 'ss_bonds': '-', 'score': 0.0,
                    'cov_ion': 0.0, 'n_matched_ion': 0,
                    'n_matched_key': 0, 'n_conflict': 0,
                    'n_disc': 0, 'disc_pairs': '-', 'disc_notes': '-',
                    'strength': '-',
                    'conclusion': 'No candidate (precursor mass failed n=0/1/2 rule)'})
        for k, j in enumerate(merged, 1):
            j['idx'] = k
        return merged

    idents = [r for r in rows_results if r.get('n_offset') != '-']
    judge_all = []
    if not idents:
        md_text = 'No identifications.'
        html_text = ''
        sample_notes = {}
        judge_all = (_allscan_judge(rows_all, matches_all)
                     if rows_all is not None and matches_all is not None else [])
        return [], md_text, html_text, sample_notes, judge_all
    judge_rows, report_rows = _judge_loop(idents, match_map)
    by_sample = {}
    for j in judge_rows:
        by_sample.setdefault(j['sample'], []).append(j)
    sample_notes = {}
    for sname, rows in by_sample.items():
        h1s = [j for j in rows if j['bridge'] == 'beads']
        h1_bonds = set()
        for j in h1s:
            h1_bonds.update(b for b in parse_bonds(j['ss_bonds']))
        notes = []
        if len(h1s) >= 2 and len(h1_bonds) >= 2:
            notes.append('MS1 dual-species evidence: intra-chain loop and '
                         'inter-chain dimer species both identified (%s), '
                         'beads (C1-C2, C3-C4) supported'
                         % (', '.join('sp%d' % j['species'] for j in h1s)))
        if any(j['bridge'] in ('globular', 'ribbon') for j in rows):
            notes.append('Single-complex species (globular/ribbon) present: globular/ribbon '
                         'discrimination needs y-ion pairing evidence (loop ions '
                         'closing an intra-chain disulfide, or single-Cys linking '
                         'ions; b-ions at low m/z end unreliable, not used)')
        if not notes:
            notes.append('No sample-level cross evidence; judge each species '
                         'independently')
        sample_notes[sname] = notes
    judge_all = (_allscan_judge(rows_all, matches_all)
                 if rows_all is not None and matches_all is not None else [])
    # ---- report text (Markdown) ----
    md = ['# DisulfID Disulfide Linkage Identification Report', '',
          '- Identified species: %d; with evidence conflict: %d'
          % (len(judge_rows), sum(1 for j in judge_rows if j['n_conflict'])),
          '- Strength thresholds: Strong(score>=%.0f and key ions>=%d), '
          'Medium(score>=%.0f and key ions>=%d)'
          % (STRONG_SCORE, STRONG_KEY_MIN, MID_SCORE, MID_KEY_MIN), '',
          '## 1. Summary', '',
          '| Idx | Sample | RT(min) | Area | M(Da) | Bridge | S-S bonds | Score | '
          'Key ions | Conflict | Strength | Conclusion |',
          '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for j in judge_rows:
        md.append('| %s | %s | %s | %s | %s | %s | %s | %.2f | %d | %d | %s | %s |'
                  % (j['idx'], j['sample'], j['RT_min'], j['precursor_Area'],
                     j['exp_neutral_M'], bridge_disp(j['bridge']), j['ss_bonds'],
                     j['score'], j['n_matched_key'], j['n_conflict'],
                     j['strength'], j['conclusion'].replace('|', '\\|')))
    if judge_all:
        md += ['', '## 2. All Scans Identification (no dedup)', '',
               '| Idx | Scan | Sample | RT(min) | Area | M(Da) | Bridge | '
               'S-S bonds | Score | Key ions | Conflict | Strength | Conclusion |',
               '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
        for j in judge_all:
            md.append('| %d | %s | %s | %s | %s | %s | %s | %s | %.2f | %d | %d | %s | %s |'
                      % (j['idx'], j['scan'], j['sample'], j['RT_min'],
                         j['precursor_Area'], j['exp_neutral_M'],
                         bridge_disp(j['bridge']), j['ss_bonds'], j['score'],
                         j['n_matched_key'], j['n_conflict'], j['strength'],
                         str(j['conclusion']).replace('|', '\\|')))
    md += ['', '## 3. Sample-Level Judgment', '']
    for sname, notes in sample_notes.items():
        md.append('**%s**: %s' % (sname, '; '.join(notes)))
    md += ['', '## 4. Key Evidence Details', '']
    for i, (row_src, j, evidence) in enumerate(report_rows, 1):
        md.append('### Species %d: scan %s | RT %s min | M %s Da | Area %s\n'
                  % (i, j['scan'], j['RT_min'], j['exp_neutral_M'],
                     j['precursor_Area']))
        md.append('- Verdict: **%s** (%s)' % (j['conclusion'], j['strength']))
        md.append('- Score %.2f = 100 x (0.7 x weighted ion coverage %.0f%% '
                  '+ 0.3 x intensity coverage %.0f%%)'
                  % (j['score'], j['cov_ion'] * 100,
                     fnum(row_src.get('intensity_coverage'), 0.0) or 0.0))
        md.append('- MS1 evidence: precursor M=%s, n_offset=%s, d=%.4f Da (%.1f ppm)'
                  % (j['exp_neutral_M'], row_src.get('n_offset'),
                     fnum(row_src.get('Δmass_Da'), 0.0) or 0.0,
                     fnum(row_src.get('Δppm'), 0.0) or 0.0))
        md.append('- Key ion evidence:')
        md.append('  | Ion | Seq | Carried | #S-S | Theo M | exp m/z | z | '
                  'dDa | dppm | W | Check |')
        md.append('  |---|---|---|---|---|---|---|---|---|---|---|')
        for m, ok, why in evidence:
            md.append('  | %s%d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |'
                      % (m.get('ion_type'), len(str(m.get('ion_seq') or '')),
                         m.get('ion_seq'), m.get('carries'),
                         m.get('nbond'), m.get('theo_neutral_mass'),
                         m.get('exp_frag_mz'), m.get('exp_charge_z'),
                         m.get('Δ_Da'), m.get('Δ_ppm'), m.get('weight'),
                         ('OK' if ok else '**Conflict**')))
        md.append('- Evidence interpretation:')
        for m, ok, why in evidence:
            md.append('  - %s%s' % (why, '' if ok else ' (**not consistent**'))
        if str(j['bridge']) in ('beads', 'globular', 'ribbon'):
            if j.get('disc_notes') and j.get('disc_notes') != '-':
                md.append('- Connectivity discrimination (loop ions closing '
                          'intra-chain disulfides / single-Cys b-/y-linking ions):')
                for n in str(j['disc_notes']).split('; '):
                    md.append('  - %s' % n.strip())
                md.append('  - Discrimination result: **%s**'
                          % j.get('disc_pairs'))
            else:
                md.append('- Connectivity discrimination: no loop-ion / '
                          'single-Cys b-/y-linking pairing evidence, connectivity '
                          'indistinguishable (use different enzyme / partial '
                          'digestion for single-Cys fragments)')
        md.append('- Overall conclusion: **%s**\n' % j['conclusion'])
    # ---- report text (HTML) ----
    hbody = ['<div style="font-family:Roboto,PingFang SC,Segoe UI,Arial,sans-serif;'
             'max-width:960px;margin:0 auto;color:#1A1B1C;">',
             '<h2>DisulfID Disulfide Linkage Identification Report</h2>',
             '<p style="color:#6B7280;font-size:13px;">Identified species %d · '
             'evidence conflicts %d</p>'
             % (len(judge_rows), sum(1 for j in judge_rows if j['n_conflict'])),
             '<h3>1. Summary</h3>',
             '<table border="1" cellspacing="0" cellpadding="4" '
             'style="border-collapse:collapse;font-size:12px;border-color:#E4E3DD;">',
             '<tr style="background:#F4F3EE;"><th>Idx</th><th>Sample</th><th>RT</th>'
             '<th>Area</th><th>M</th><th>Bridge</th><th>S-S bonds</th><th>Score</th>'
             '<th>Key ions</th><th>Conflict</th><th>Strength</th><th>Conclusion</th></tr>']
    for j in judge_rows:
        hbody.append('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>'
                     '<td>%s</td><td>%s</td><td>%.2f</td><td>%d</td><td>%d</td>'
                     '<td>%s</td><td>%s</td></tr>'
                     % (j['idx'], _esc(j['sample']), _esc(j['RT_min']),
                        _esc(j['precursor_Area']), _esc(j['exp_neutral_M']),
                        bridge_disp(j['bridge']), _esc(j['ss_bonds']), j['score'],
                        j['n_matched_key'], j['n_conflict'], j['strength'],
                        _esc(j['conclusion'])))
    if judge_all:
        hbody.append('<h3>2. All Scans Identification (no dedup)</h3>')
        hbody.append('<p style="color:#6B7280;font-size:13px;">Every scan, '
                     'one row (best candidate), no RT/M dedup.</p>')
        hbody.append('<table border="1" cellspacing="0" cellpadding="4" '
                     'style="border-collapse:collapse;font-size:12px;'
                     'border-color:#E4E3DD;">'
                     '<tr style="background:#F4F3EE;"><th>Idx</th><th>Scan</th>'
                     '<th>Sample</th><th>RT</th><th>Area</th><th>M</th>'
                     '<th>Bridge</th><th>S-S bonds</th><th>Score</th>'
                     '<th>Key ions</th><th>Conflict</th><th>Strength</th>'
                     '<th>Conclusion</th></tr>')
        for j in judge_all:
            hbody.append('<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td>'
                         '<td>%s</td><td>%s</td><td>%s</td><td>%s</td>'
                         '<td>%.2f</td><td>%d</td><td>%d</td><td>%s</td>'
                         '<td>%s</td></tr>'
                         % (j['idx'], _esc(j['scan']), _esc(j['sample']),
                            _esc(j['RT_min']), _esc(j['precursor_Area']),
                            _esc(j['exp_neutral_M']), bridge_disp(j['bridge']),
                            _esc(j['ss_bonds']), j['score'],
                            j['n_matched_key'], j['n_conflict'], j['strength'],
                            _esc(j['conclusion'])))
        hbody.append('</table>')
    hbody.append('<h3>3. Sample-Level Judgment</h3>')
    for sname, notes in sample_notes.items():
        hbody.append('<p style="font-size:13px;"><b>%s</b>: %s</p>'
                     % (_esc(sname), '; '.join(_esc(n) for n in notes)))
    hbody.append('<h3>4. Key Evidence Details</h3>')
    for i, (row_src, j, evidence) in enumerate(report_rows, 1):
        hbody.append('<h4>Species %d: scan %s | RT %s min | M %s Da</h4>'
                     % (i, j['scan'], j['RT_min'], j['exp_neutral_M']))
        hbody.append('<p style="font-size:13px;">Verdict: <b>%s</b> (%s)<br/>'
                     'Score <b>%.2f</b><br/>MS1: M=%s, n_offset=%s, d=%.4f Da</p>'
                     % (_esc(j['conclusion']), j['strength'], j['score'],
                        _esc(j['exp_neutral_M']), row_src.get('n_offset'),
                        fnum(row_src.get('Δmass_Da'), 0.0) or 0.0))
        hbody.append('<table border="1" cellspacing="0" cellpadding="3" '
                     'style="border-collapse:collapse;font-size:11.5px;'
                     'border-color:#E4E3DD;"><tr style="background:#F4F3EE;">'
                     '<th>Ion</th><th>Seq</th><th>Carried</th><th>#S-S</th>'
                     '<th>Theo M</th><th>exp m/z</th><th>z</th><th>dDa</th>'
                     '<th>dppm</th><th>W</th><th>Check</th></tr>')
        for m, ok, why in evidence:
            hbody.append('<tr><td>%s%d</td><td>%s</td><td>%s</td><td>%s</td>'
                         '<td>%s</td><td>%s</td><td>%s</td><td>%s</td>'
                         '<td>%s</td><td>%s</td><td>%s</td></tr>'
                         % (m.get('ion_type'), len(str(m.get('ion_seq') or '')),
                            _esc(m.get('ion_seq')), _esc(m.get('carries')),
                            m.get('nbond'), m.get('theo_neutral_mass'),
                            m.get('exp_frag_mz'), m.get('exp_charge_z'),
                            m.get('Δ_Da'), m.get('Δ_ppm'), m.get('weight'),
                            ('OK' if ok else
                             '<span style="color:#EA6668;">Conflict</span>')))
        hbody.append('</table>')
        for m, ok, why in evidence:
            hbody.append('<p style="font-size:12.5px;margin:2px 0;">&#8226; %s%s</p>'
                         % (_esc(why), '' if ok else
                            ' (<b style="color:#EA6668;">not consistent</b>)'))
        if str(j['bridge']) in ('beads', 'globular', 'ribbon'):
            if j.get('disc_notes') and j.get('disc_notes') != '-':
                hbody.append('<p style="font-size:12.5px;margin:2px 0;">'
                             '<b style="color:#1f4e79;">Connectivity discrimination'
                             ' (loop ions closing intra-chain disulfides / '
                             'single-Cys b-/y-linking ions):</b><br/>%s</p>'
                             % '<br/>'.join(
                                 _esc(n.strip()) for n in
                                 str(j['disc_notes']).split('; ')))
            else:
                hbody.append('<p style="font-size:12.5px;margin:2px 0;">'
                             '<b style="color:#b91c1c;">Connectivity '
                             'discrimination:</b> no loop-ion / single-Cys '
                             'b-/y-linking pairing evidence, connectivity '
                             'indistinguishable (use different enzyme / partial '
                             'digestion for single-Cys fragments)</p>')
        hbody.append('<p style="font-size:13px;font-weight:600;">Overall '
                     'conclusion: %s</p>' % _esc(j['conclusion']))
    hbody.append('</div>')
    html_text = ('<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                 '<title>DisulfID Identification Report</title>\n</head>\n'
                 '<body style="margin:16px;background:#FFFFFF;">\n%s\n</body>\n</html>'
                 % '\n'.join(hbody))
    msgs.append('Judgment: %d species, %d with evidence conflict'
                % (len(judge_rows),
                   sum(1 for j in judge_rows if j['n_conflict'])))
    if judge_all:
        msgs.append('All-scan view: %d scans (%d identified, %d no candidate)'
                    % (len(judge_all),
                       sum(1 for j in judge_all if j['n_matched_key'] or
                           j.get('combo_id') != '-'),
                       sum(1 for j in judge_all
                           if str(j['conclusion']).startswith('No candidate'))))
    return judge_rows, '\n'.join(md), html_text, sample_notes, judge_all

def _esc(s):
    s = str(s)
    return (s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
             .replace('"', '&quot;'))

# ===================== Annotated MS2 spectrum =====================
CHAIN_COLORS = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
                '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
UNMATCH_COLOR = '#c0c0c0'

def plot_annotated_ms2(frags, details, compo_parts, row, out_png):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    if not frags:
        return False
    # normalization: base = highest peak within plotting range (200-1200);
    # if no peak in range, use global highest peak
    in_range = [(fr['intensity'] or 0.0) for fr in frags
                if fr['mz'] is not None and 200.0 <= fr['mz'] <= 1200.0]
    base = max(in_range) if in_range else max(
        (fr['intensity'] or 0.0) for fr in frags)
    base = base or 1.0
    mz_all = [fr['mz'] for fr in frags if fr['mz'] is not None]
    if not mz_all:
        return False
    matched_fi = {dt['fi'] for dt in details}
    fig, ax = plt.subplots(figsize=(13, 5.5), dpi=150)
    for fi, fr in enumerate(frags):
        mz = fr['mz']
        if mz is None:
            continue
        rel = (fr['intensity'] or 0.0) / base * 100.0
        color = UNMATCH_COLOR
        if fi in matched_fi:
            dt = next(d for d in details if d['fi'] == fi)
            ion = dt['ion']
            fid = ion['fid']
            color = CHAIN_COLORS[fid % len(CHAIN_COLORS)]
            label = '%s%d' % (ion['type'], len(ion['seq']))
            zlab = ('z=%g' % dt['fr']['z']) if (dt['fr']['z'] is not None) else 'z=?'
            ax.annotate(label, xy=(mz, rel), xytext=(mz, rel + 3),
                        ha='center', va='bottom', fontsize=7.5, color=color,
                        fontweight='bold',
                        arrowprops=dict(arrowstyle='-', color=color, lw=0.6))
            ax.text(mz, rel - 2, '%.4f  %s' % (mz, zlab), ha='center',
                    va='top', fontsize=6.5, color='#444444')
        ax.plot([mz, mz], [0, rel], color=color, lw=1.1,
                solid_capstyle='round', zorder=2)
    ax.set_xlim(200, 1200)
    ax.set_ylim(0, 112)
    ax.set_xlabel('m/z', fontsize=10)
    ax.set_ylabel('Relative Abundance (%)', fontsize=10)
    ax.tick_params(labelsize=8)
    ax.set_title('DisulfID | scan %s | RT %.2f min | precursor %.4f (%s) | M %.4f | '
                 'score %.2f | combo%s %s sp%s | %s'
                 % (row['scan'], fnum(row.get('RT_min'), 0.0) or 0.0,
                    fnum(row.get('precursor_mz'), 0.0) or 0.0,
                    row.get('precursor_charge'),
                    fnum(row.get('exp_neutral_M'), 0.0) or 0.0,
                    row.get('score'), row.get('combo_id'),
                    bridge_disp(row.get('bridge')),
                    row.get('species'),
                    str(row.get('composition')).replace(' + ', '/')),
                 fontsize=9, color='#1A1B1C')
    ss = str(row.get('ss_bonds') or '').strip()
    ss_text = ('linear peptide (no S-S)' if ss in ('-', '')
               else 'S-S bonds: %s' % ss)
    ax.text(0.985, 0.97, ss_text, transform=ax.transAxes, fontsize=8.5,
            ha='right', va='top', color='#d62728', fontweight='bold')
    handles = []
    for i, part in enumerate(compo_parts):
        handles.append(plt.Line2D([0], [0],
                                  color=CHAIN_COLORS[i % len(CHAIN_COLORS)],
                                  lw=3, label=part))
    handles.append(plt.Line2D([0], [0], color=UNMATCH_COLOR, lw=3,
                              label='unmatched'))
    ax.legend(handles=handles, fontsize=8, loc='upper left', framealpha=0.9,
              ncol=2)
    fig.tight_layout()
    fig.savefig(out_png)
    # also export vector PDF: every peak/label/axis is an independent path/text,
    # editable in Adobe Illustrator / Inkscape, suitable for publication
    try:
        out_pdf = os.path.splitext(out_png)[0] + '.pdf'
        fig.savefig(out_pdf)
    except Exception:
        pass
    plt.close(fig)
    return True

# ===================== Pipeline orchestration (GUI and CLI) =====================
def run_pipeline(mgf_path, raw_seq, out_dir, params, msgs):
    """Five-step pipeline, returns result dict (for GUI/CLI display and export)."""
    msgs.append('[1/5] Building theoretical library ...')
    species_rows, ions_index, combos = build_library(
        raw_seq, msgs,
        mz_min=params.get('mz_min', MZ_MIN),
        mz_max=params.get('mz_max', MZ_MAX))
    msgs.append('[2/5] Theoretical b/y fragment library built')
    msgs.append('[3/5] Parsing MGF: %s' % os.path.basename(mgf_path))
    specs, frag_map = parse_mgf_file(mgf_path, params, msgs)
    if not specs:
        raise ValueError('No usable spectra after MGF parsing (no PEPMASS '
                         'intensity or Area < %.0f)' % params['ms1_filter'])
    msgs.append('[4/5] Library search & scoring ...')
    (results_all, matches_all, rows_results, rows_matches,
     spec_store, stats) = run_search(species_rows, ions_index, specs,
                                     frag_map, params, msgs)
    msgs.append('[5/5] Judgment & report generation ...')
    judge_rows, md_text, html_text, sample_notes, judge_all = judge_and_report(
        rows_results, rows_matches, msgs,
        rows_all=results_all, matches_all=matches_all)
    # Backfill judgment conclusion / S-S display into result rows: for the
    # degenerate case show "globular/ribbon" instead of the initial hypothesis so it is
    # not misread as confirmed; MS2 annotation label stays in sync
    judge_map = {}
    for j in judge_rows:
        judge_map[(str(j.get('scan')), inum(j.get('combo_id')),
                   j.get('bridge'), inum(j.get('species')))] = j
    for r in rows_results:
        if r.get('n_offset') == '-':
            continue
        j = judge_map.get((str(r['scan']), inum(r.get('combo_id')),
                           r.get('bridge'), inum(r.get('species'))))
        if j:
            r['ss_bonds'] = j.get('ss_bonds', r.get('ss_bonds'))
            r['conclusion'] = j.get('conclusion', r.get('conclusion'))
    # annotated spectra
    n_plot = 0
    if params.get('plot', True):
        plot_dir = os.path.join(out_dir, 'ms2_annotated')
        try:
            os.makedirs(plot_dir, exist_ok=True)
            for r in rows_results:
                if r.get('n_offset') == '-':
                    continue
                if (fnum(r.get('score'), 0.0) or 0.0) <= params.get(
                        'plot_score_min', 10.0):
                    continue
                skey = (str(r['scan']), inum(r.get('combo_id')),
                        r.get('bridge'), inum(r.get('species')))
                store = spec_store.get(skey)
                if not store or not store['details']:
                    continue
                compo_parts = [p.strip() for p in str(r.get('composition')).split('+')]
                scan_str = str(r['scan']).replace('/', '_').replace('\\', '_')
                fname = ('annotated_scan%s_combo%s_%s_sp%s.png'
                         % (scan_str, r.get('combo_id'),
                            bridge_disp(r.get('bridge')),
                            r.get('species')))
                out_png = os.path.join(plot_dir, fname)
                if plot_annotated_ms2(store['frags'], store['details'],
                                      compo_parts, r, out_png):
                    n_plot += 1
        except Exception as e:
            msgs.append('Plot error (ignorable): %s' % e)
    msgs.append('Done: spectra %d | identified species %d | annotated spectra %d'
                % (stats['tot_spec'], len(judge_rows), n_plot))
    return {'species': species_rows, 'ions_index': ions_index,
            'results': rows_results, 'results_all': results_all,
            'matches': rows_matches, 'matches_all': matches_all,
            'judge': judge_rows, 'judge_all': judge_all,
            'md': md_text, 'html': html_text,
            'notes': sample_notes, 'stats': stats, 'n_plot': n_plot}

def export_outputs(payload, out_dir):
    """Export judgment.csv / report.md / report.html (spectra in ms2_annotated/)"""
    CSV_HEADER_MAP = {
        'idx': 'Idx', 'scan': 'Scan', 'sample': 'Sample', 'RT_min': 'RT(min)',
        'precursor_Area': 'Area', 'exp_neutral_M': 'M(Da)', 'combo_id': 'Combo',
        'bridge': 'Bridge', 'species': 'Sp', 'composition': 'Composition',
        'ss_bonds': 'S-S bonds', 'score': 'Score', 'cov_ion': 'IonCov',
        'n_matched_ion': 'MatchedIons', 'n_matched_key': 'KeyIons',
        'n_conflict': 'Conflict', 'n_disc': 'nDisc',
        'disc_pairs': 'DiscPairs',
        'disc_notes': 'DiscNotes', 'strength': 'Strength', 'conclusion': 'Conclusion',
    }
    if payload.get('judge'):
        p = os.path.join(out_dir, 'disulfide_judgment.csv')
        with open(p, 'w', newline='', encoding='utf-8-sig') as fh:
            keys = list(payload['judge'][0].keys())
            headers = [CSV_HEADER_MAP.get(k, k) for k in keys]
            w = csv.DictWriter(fh, fieldnames=headers)
            w.writeheader()
            for j in payload['judge']:
                w.writerow({CSV_HEADER_MAP.get(k, k): j.get(k, '')
                            for k in keys})
    if payload.get('judge_all'):
        p = os.path.join(out_dir, 'disulfide_judgment_allscans.csv')
        with open(p, 'w', newline='', encoding='utf-8-sig') as fh:
            keys = list(payload['judge_all'][0].keys())
            headers = [CSV_HEADER_MAP.get(k, k) for k in keys]
            w = csv.DictWriter(fh, fieldnames=headers)
            w.writeheader()
            for j in payload['judge_all']:
                w.writerow({CSV_HEADER_MAP.get(k, k): j.get(k, '')
                            for k in keys})
    p = os.path.join(out_dir, 'disulfide_report.md')
    with open(p, 'w', encoding='utf-8') as fh:
        fh.write(payload['md'] or 'No identifications.')
    p = os.path.join(out_dir, 'disulfide_report.html')
    with open(p, 'w', encoding='utf-8') as fh:
        fh.write(payload['html'] or '<p>No identifications.</p>')
    return os.path.join(out_dir, 'disulfide_report.html')

# ===================== Command-line mode =====================
def run_cli(mgf_path, raw_seq, out_dir):
    msgs = []
    print('=' * 72)
    print('DisulfID - Disulfide Linkage Identification (CLI mode)')
    print('=' * 72)
    payload = run_pipeline(mgf_path, raw_seq, out_dir,
                           {'ppm1': 20, 'ppm2': 20, 'cutoff': 100,
                            'ms1_filter': 10000, 'plot': True,
                            'plot_score_min': 10.0,
                            'dedup_rt': 1.0, 'dedup_mass': 0.02,
                            'mz_min': MZ_MIN, 'mz_max': MZ_MAX}, msgs)
    for m in msgs:
        print(m)
    report = export_outputs(payload, out_dir)
    print('-' * 72)
    if payload['judge']:
        print('%-5s %-8s %-8s %-12s %-16s %7s %-16s %s'
              % ('Idx', 'Bridge', 'RT', 'M', 'S-S bonds', 'Score',
                 'Strength', 'Conclusion'))
        for j in payload['judge']:
            print('%-5d %-8s %-8s %-12s %-16s %7.2f %-16s %s'
                  % (j['idx'], bridge_disp(j['bridge']), j['RT_min'],
                     j['exp_neutral_M'], j['ss_bonds'], j['score'],
                     j['strength'], j['conclusion']))
    print('Exported: %s / disulfide_judgment.csv / disulfide_report.md'
          % report)

# ===================== PyQt6 graphical interface =====================
try:
    from PyQt6.QtCore import Qt, QThread, pyqtSignal
    from PyQt6.QtGui import (QPixmap, QIcon, QPainter, QPen, QColor, QBrush)
    from PyQt6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QGroupBox, QLabel, QLineEdit, QPushButton, QFileDialog, QTableWidget,
        QTableWidgetItem, QHeaderView, QTextEdit, QTabWidget, QSplitter,
        QProgressBar, QFormLayout, QMessageBox, QComboBox)
    PYQT_OK = True
except ImportError:
    PYQT_OK = False

if PYQT_OK:
    class PipelineWorker(QThread):
        log_sig = pyqtSignal(str)
        done_sig = pyqtSignal(bool, str, object)
        def __init__(self, mgf_path, raw_seq, out_dir, params, parent=None):
            super().__init__(parent)
            self.mgf_path, self.raw_seq, self.out_dir, self.params = (
                mgf_path, raw_seq, out_dir, params)
        def run(self):
            msgs = []
            try:
                payload = run_pipeline(self.mgf_path, self.raw_seq,
                                       self.out_dir, self.params, msgs)
                for m in msgs:
                    self.log_sig.emit(m)
                self.done_sig.emit(True, '', payload)
            except Exception as e:
                import traceback
                tb = traceback.format_exc()
                for m in msgs:
                    self.log_sig.emit(m)
                self.log_sig.emit('Error: %s' % e)
                self.log_sig.emit(tb)
                self.done_sig.emit(False, '%s' % e, None)

    class MainWindow(QMainWindow):
        RESULT_COLS = ['scan', 'RT_min', 'precursor_Area', 'exp_neutral_M',
                       'n_offset', 'bridge', 'ss_bonds', 'score', 'n_matched_key',
                       'strength', 'conclusion']
        RESULT_HEADERS = ['Scan', 'RT(min)', 'Area', 'M(Da)', 'n_offset',
                          'Bridge', 'S-S bonds', 'Score', 'Key ions',
                          'Strength', 'Conclusion']
        EVID_COLS = ['ion', 'seq', 'carries', 'nbond', 'theo_m', 'exp_mz', 'z',
                     'delta_da', 'delta_ppm', 'weight', 'check']
        EVID_HEADERS = ['Ion', 'Seq', 'Carried', '#S-S', 'Theo M', 'exp m/z',
                        'z', 'dDa', 'dppm', 'W', 'Check']
        def __init__(self):
            super().__init__()
            self.setWindowTitle('DisulfID - Disulfide Linkage Identification')
            self.resize(1180, 860)
            self.setWindowIcon(self._build_icon())
            self._build_menu()
            self.payload = None
            self._build_ui()
            self.worker = None

        # ---- Application icon (programmatic drawing, no external icon file) ----
        @staticmethod
        def _build_icon():
            pix = QPixmap(64, 64)
            pix.fill(Qt.GlobalColor.transparent)
            qp = QPainter(pix)
            qp.setRenderHint(QPainter.RenderHint.Antialiasing)
            pen_b = QPen(QColor('#2b6cb0'), 6)
            qp.setPen(pen_b)
            qp.drawLine(18, 18, 46, 18)   # upper peptide chain
            qp.drawLine(18, 46, 46, 46)   # lower peptide chain
            qp.drawLine(18, 18, 18, 46)   # left U-turn
            qp.setPen(QPen(QColor('#d9a413'), 3))
            qp.drawLine(30, 18, 30, 46)   # S-S bridge 1
            qp.drawLine(40, 18, 40, 46)   # S-S bridge 2
            qp.setPen(QPen(QColor('#d9a413'), 1.2))
            qp.setBrush(QBrush(QColor('#f6d743')))
            for (x, y) in ((30, 27), (30, 37), (40, 27), (40, 37)):
                qp.drawEllipse(x - 3, y - 3, 6, 6)
            qp.end()
            return QIcon(pix)

        def _build_menu(self):
            mb = self.menuBar()
            help_menu = mb.addMenu('&Help')
            act = help_menu.addAction('About DisulfID')
            act.triggered.connect(self._about)

        def _about(self):
            QMessageBox.about(
                self, 'About DisulfID',
                '<h3 style="color:#1f4e79;">DisulfID v1.0</h3>'
                '<p><b>Disulfide Linkage Identification System</b></p>'
                '<p>Enzymatic digestion + Orbitrap MS: automatic disulfide '
                'linkage identification for single-chain peptides (2 '
                'disulfide bonds).</p>'
                '<p>Pipeline: theoretical library -> fragment library -> '
                'MGF parsing -> search & scoring -> judgment & report</p>'
                '<p style="color:#6b7280;">Bridge hypotheses: beads '
                '(C1-C2, C3-C4) | globular (C1-C3, C2-C4) | ribbon '
                '(C1-C4, C2-C3)</p>')

        def _build_ui(self):
            central = QWidget(self)
            self.setCentralWidget(central)
            root = QVBoxLayout(central)

            # ---- title bar ----
            title = QLabel('DisulfID - Disulfide Linkage Identification')
            title.setStyleSheet('font-size:15px; font-weight:600; color:#1f4e79;')
            root.addWidget(title)

            # ---- input area ----
            g_in = QGroupBox('Input')
            f1 = QFormLayout(g_in)
            row1 = QHBoxLayout()
            self.edit_mgf = QLineEdit()
            self.edit_mgf.setPlaceholderText('Select MGF file ...')
            b1 = QPushButton('Browse...')
            b1.clicked.connect(self._pick_mgf)
            row1.addWidget(self.edit_mgf, 1)
            row1.addWidget(b1)
            self.edit_seq = QLineEdit(DEFAULT_SEQ)
            row2 = QHBoxLayout()
            self.edit_out = QLineEdit()
            self.edit_out.setPlaceholderText('Output directory (default: MGF folder)')
            b2 = QPushButton('Browse...')
            b2.clicked.connect(self._pick_out)
            row2.addWidget(self.edit_out, 1)
            row2.addWidget(b2)
            f1.addRow('MGF file', row1)
            f1.addRow('Peptide sequence', self.edit_seq)
            f1.addRow('Output dir', row2)
            # advanced parameters
            g_adv = QGroupBox('Advanced parameters')
            fa = QFormLayout(g_adv)
            self.e_ppm1 = QLineEdit('20')
            self.e_ppm2 = QLineEdit('20')
            self.e_cutoff = QLineEdit('100')
            self.e_ms1 = QLineEdit('10000')
            self.e_plot = QLineEdit('10')
            self.e_rt = QLineEdit('1.0')
            self.e_mass = QLineEdit('0.02')
            self.e_mzmin = QLineEdit('200')
            self.e_mzmax = QLineEdit('1200')
            for w in (self.e_ppm1, self.e_ppm2, self.e_cutoff, self.e_ms1,
                      self.e_plot, self.e_rt, self.e_mass,
                      self.e_mzmin, self.e_mzmax):
                w.setFixedWidth(70)
            fa.addRow('MS1 ppm', self.e_ppm1)
            fa.addRow('MS2 ppm', self.e_ppm2)
            fa.addRow('MS2 peak intensity cutoff', self.e_cutoff)
            fa.addRow('MS1 Area cutoff', self.e_ms1)
            fa.addRow('Plot score threshold', self.e_plot)
            fa.addRow('Dedup RT tol (min)', self.e_rt)
            fa.addRow('Dedup mass tol (Da)', self.e_mass)
            fa.addRow('Fragment m/z min', self.e_mzmin)
            fa.addRow('Fragment m/z max', self.e_mzmax)
            g_adv.setMaximumWidth(300)
            row_top = QHBoxLayout()
            row_top.addWidget(g_in, 1)
            row_top.addWidget(g_adv)
            root.addLayout(row_top)

            # ---- run buttons ----
            row_btn = QHBoxLayout()
            self.btn_run = QPushButton('Run Identification')
            self.btn_run.clicked.connect(self._start)
            self.btn_export = QPushButton('Export Report')
            self.btn_export.setEnabled(False)
            self.btn_export.clicked.connect(self._export)
            self.btn_open = QPushButton('Open Output Folder')
            self.btn_open.setEnabled(False)
            self.btn_open.clicked.connect(self._open_out)
            self.progress = QProgressBar()
            self.progress.setRange(0, 0)
            self.progress.setVisible(False)
            row_btn.addWidget(self.btn_run)
            row_btn.addWidget(self.btn_export)
            row_btn.addWidget(self.btn_open)
            row_btn.addWidget(self.progress, 1)
            root.addLayout(row_btn)

            # ---- results area ----
            split = QSplitter(Qt.Orientation.Vertical)
            self.table = QTableWidget(0, len(self.RESULT_COLS))
            self.table.setHorizontalHeaderLabels(self.RESULT_HEADERS)
            self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
            self.table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
            self.table.horizontalHeader().setSectionResizeMode(
                QHeaderView.ResizeMode.Interactive)
            self.table.itemSelectionChanged.connect(self._on_select)

            # ---- results view switch: peak-level (dedup) vs all scans ----
            view_row = QHBoxLayout()
            view_row.addWidget(QLabel('Results view:'))
            self.view_sel = QComboBox()
            self.view_sel.addItems(['Peak-level (dedup)',
                                    'All scans (no dedup)'])
            self.view_sel.setMinimumWidth(220)
            self.view_sel.currentIndexChanged.connect(self._reload_view)
            view_row.addWidget(self.view_sel)
            view_row.addStretch(1)
            view_host = QWidget()
            view_host.setLayout(view_row)
            split.addWidget(view_host)
            split.addWidget(self.table)
            tabs = QTabWidget()
            self.ev_table = QTableWidget(0, len(self.EVID_COLS))
            self.ev_table.setHorizontalHeaderLabels(self.EVID_HEADERS)
            self.ev_table.horizontalHeader().setSectionResizeMode(
                QHeaderView.ResizeMode.Interactive)
            self.txt_why = QTextEdit()
            self.txt_why.setReadOnly(True)
            self.lbl_plot = QLabel('(select a result row to show the annotated '
                                   'MS2 spectrum; score must be above the plot '
                                   'threshold)')
            self.lbl_plot.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.lbl_plot.setMinimumHeight(300)
            tabs.addTab(self.ev_table, 'Key Evidence')
            tabs.addTab(self.txt_why, 'Evidence Interpretation')
            tabs.addTab(self.lbl_plot, 'MS2 Annotated Spectrum')
            split.addWidget(tabs)
            split.setSizes([360, 400])
            root.addWidget(split, 1)

            # ---- judgment + log ----
            tabs2 = QTabWidget()
            self.txt_judge = QTextEdit()
            self.txt_judge.setReadOnly(True)
            self.txt_log = QTextEdit()
            self.txt_log.setReadOnly(True)
            tabs2.addTab(self.txt_judge, 'Judgment')
            tabs2.addTab(self.txt_log, 'Log')
            tabs2.setMaximumHeight(170)
            root.addWidget(tabs2)

        def _pick_mgf(self):
            p, _ = QFileDialog.getOpenFileName(
                self, 'Select MGF file', '', 'Mass spectra (*.mgf);;All files (*.*)')
            if p:
                self.edit_mgf.setText(p)
                if not self.edit_out.text().strip():
                    self.edit_out.setText(os.path.dirname(p))
        def _pick_out(self):
            d = QFileDialog.getExistingDirectory(self, 'Select output directory',
                                                 self.edit_out.text().strip() or '')
            if d:
                self.edit_out.setText(d)
        def _open_out(self):
            d = self.edit_out.text().strip()
            if d and os.path.isdir(d):
                os.startfile(d)   # Windows only

        def _params(self):
            def f(w, dflt):
                try:
                    return float(w.text().strip())
                except ValueError:
                    return dflt
            return {'ppm1': f(self.e_ppm1, 20), 'ppm2': f(self.e_ppm2, 20),
                    'cutoff': f(self.e_cutoff, 100),
                    'ms1_filter': f(self.e_ms1, 10000),
                    'plot': True, 'plot_score_min': f(self.e_plot, 10),
                    'dedup_rt': f(self.e_rt, 1.0),
                    'dedup_mass': f(self.e_mass, 0.02),
                    'mz_min': f(self.e_mzmin, MZ_MIN),
                    'mz_max': f(self.e_mzmax, MZ_MAX)}

        def _start(self):
            mgf = self.edit_mgf.text().strip()
            seq = self.edit_seq.text().strip()
            if not mgf:
                QMessageBox.warning(self, 'Notice', 'Please select an MGF file first')
                return
            if not os.path.exists(mgf):
                QMessageBox.warning(self, 'Notice', 'MGF file not found:\n%s' % mgf)
                return
            if not seq:
                QMessageBox.warning(self, 'Notice', 'Please enter the peptide sequence')
                return
            out_dir = self.edit_out.text().strip() or os.path.dirname(mgf)
            try:
                os.makedirs(out_dir, exist_ok=True)
            except Exception as e:
                QMessageBox.warning(self, 'Notice', 'Output directory unavailable:\n%s' % e)
                return
            self.btn_run.setEnabled(False)
            self.btn_export.setEnabled(False)
            self.progress.setVisible(True)
            self.txt_log.clear()
            self.table.setRowCount(0)
            self.ev_table.setRowCount(0)
            self.txt_why.clear()
            self.txt_judge.clear()
            self.lbl_plot.setText('(running...)')
            self.payload = None
            self.worker = PipelineWorker(mgf, seq, out_dir, self._params(), self)
            self.worker.log_sig.connect(lambda s: self.txt_log.append(s))
            self.worker.done_sig.connect(self._on_done)
            self.worker.start()

        def _on_done(self, ok, err, payload):
            self.progress.setVisible(False)
            self.btn_run.setEnabled(True)
            if not ok:
                QMessageBox.critical(self, 'Run failed', str(err))
                return
            self.payload = payload
            self.view_sel.setCurrentIndex(0)
            self._reload_view()
            self._fill_judge(payload)
            self.btn_export.setEnabled(True)
            self.btn_open.setEnabled(True)
            n_id = len(payload['judge'])
            QMessageBox.information(
                self, 'Done',
                'Identification finished: spectra %d | identified species %d | '
                'annotated spectra %d\nSelect a result row to view key evidence '
                'and the annotated MS2 spectrum; click "Export Report" to export.'
                % (payload['stats']['tot_spec'], n_id, payload['n_plot']))

        def _reload_view(self):
            if not self.payload:
                return
            if self.view_sel.currentIndex() == 0:
                self._fill_table(self.payload['results'])
            else:
                self._fill_table(self.payload.get('judge_all', []))

        def _fill_table(self, results):
            self._rows = results
            self.table.setRowCount(0)
            for r in results:
                row = self.table.rowCount()
                self.table.insertRow(row)
                vals = [r.get(c, '') for c in self.RESULT_COLS]
                for ci, v in enumerate(vals):
                    item = QTableWidgetItem('' if v is None else str(v))
                    if ci == self.RESULT_COLS.index('score'):
                        try:
                            item.setText('%.2f' % float(v))
                        except (ValueError, TypeError):
                            pass
                    if ci == self.RESULT_COLS.index('bridge'):
                        item.setForeground(Qt.GlobalColor.darkBlue)
                    self.table.setItem(row, ci, item)

        def _on_select(self):
            rows = self.table.selectionModel().selectedRows()
            if not rows or not self.payload or not hasattr(self, '_rows'):
                return
            row_i = rows[0].row()
            if row_i >= len(self._rows):
                return
            r = self._rows[row_i]
            if (r.get('n_offset') == '-'
                    or str(r.get('combo_id')) in ('-', '')):
                self.ev_table.setRowCount(0)
                self.txt_why.setText('No candidate for this spectrum (precursor '
                                     'mass failed n=0/1/2 rule)')
                self.lbl_plot.setText('(no spectrum)')
                return
            key = (str(r['scan']), inum(r.get('combo_id')), r.get('bridge'),
                   inum(r.get('species')))
            matches = [m for m in (self.payload.get('matches_all')
                                   or self.payload['matches']) if (
                str(m.get('scan')) == str(r['scan'])
                and inum(m.get('combo_id')) == inum(r.get('combo_id'))
                and m.get('bridge') == r.get('bridge')
                and inum(m.get('species')) == inum(r.get('species')))]
            self.ev_table.setRowCount(0)
            for m in matches:
                row = self.ev_table.rowCount()
                self.ev_table.insertRow(row)
                ilabel = '%s%d' % (m.get('ion_type'),
                                   len(str(m.get('ion_seq') or '')))
                vals = [ilabel, m.get('ion_seq'), m.get('carries'),
                        m.get('nbond'), m.get('theo_neutral_mass'),
                        m.get('exp_frag_mz'), m.get('exp_charge_z'),
                        m.get('Δ_Da'), m.get('Δ_ppm'), m.get('weight'),
                        '—']
                for ci, v in enumerate(vals):
                    self.ev_table.setItem(row, ci, QTableWidgetItem(
                        '' if v is None else str(v)))
            why = []
            for j in self.payload['judge']:
                if (str(j.get('scan')) == str(r['scan'])
                        and str(j.get('combo_id')) == str(r.get('combo_id'))
                        and j.get('bridge') == r.get('bridge')):
                    why.append('Verdict: %s (%s)' % (j['conclusion'], j['strength']))
            if not why:
                why.append('Verdict: %s (%s)'
                           % (r.get('conclusion'), r.get('n_matched_key')))
            why.append('')
            for m in matches:
                carry = m.get('carries')
                is_key = (carry not in ('-', None)) or inum(
                    m.get('nbond'), 0) > 0
                if is_key:
                    why.append('Key ion %s%d (%s): theo %.4f -> exp m/z %s '
                               '(z=%s), d=%.4f Da (%.2f ppm), carried=%s'
                               % (m.get('ion_type'),
                                  len(str(m.get('ion_seq') or '')),
                                  m.get('ion_seq'),
                                  fnum(m.get('theo_neutral_mass'), 0.0) or 0.0,
                                  m.get('exp_frag_mz'), m.get('exp_charge_z'),
                                  fnum(m.get('Δ_Da'), 0.0) or 0.0,
                                  fnum(m.get('Δ_ppm'), 0.0) or 0.0, carry))
            self.txt_why.setText('\n'.join(why))
            # annotated spectrum
            scan_str = str(r['scan']).replace('/', '_').replace('\\', '_')
            png = os.path.join(self.edit_out.text().strip() or
                               os.path.dirname(self.edit_mgf.text().strip()),
                               'ms2_annotated',
                               'annotated_scan%s_combo%s_%s_sp%s.png'
                               % (scan_str, r.get('combo_id'),
                                  bridge_disp(r.get('bridge')),
                                  r.get('species')))
            if os.path.exists(png):
                pix = QPixmap(png)
                if not pix.isNull():
                    self.lbl_plot.setPixmap(pix.scaled(
                        self.lbl_plot.width(),
                        max(300, self.lbl_plot.height()),
                        Qt.AspectRatioMode.KeepAspectRatio,
                        Qt.TransformationMode.SmoothTransformation))
                else:
                    self.lbl_plot.setText('(spectrum load failed)')
            else:
                self.lbl_plot.setText('(score below plot threshold, no spectrum)')

        def _fill_judge(self, payload):
            jd, notes = payload['judge'], payload['notes']
            lines = ['=== DisulfID Judgment ===', '']
            lines.append('Spectra %d | candidates %d | hits %d | identified %d | '
                         'spectra %d'
                         % (payload['stats']['tot_spec'],
                            payload['stats']['tot_cand'],
                            payload['stats']['tot_hit'],
                            len(jd), payload['n_plot']))
            lines.append('')
            for sname, ns in notes.items():
                lines.append('[Sample %s]' % sname)
                lines.extend('  - ' + n for n in ns)
            lines.append('')
            if jd:
                lines.append('%-4s %-8s %-10s %-11s %-16s %7s %-16s %s'
                             % ('#', 'Bridge', 'RT', 'M', 'S-S bonds', 'Score',
                                'Strength', 'Conclusion'))
                for j in jd:
                    lines.append('%-4d %-8s %-10s %-11s %-16s %7.2f %-16s %s'
                                 % (j['idx'], bridge_disp(j['bridge']), j['RT_min'],
                                    j['exp_neutral_M'], j['ss_bonds'],
                                    j['score'], j['strength'], j['conclusion']))
            self.txt_judge.setText('\n'.join(lines))

        def _export(self):
            if not self.payload:
                return
            d = QFileDialog.getExistingDirectory(
                self, 'Select export directory', self.edit_out.text().strip() or '')
            if not d:
                return
            try:
                report = export_outputs(self.payload, d)
                QMessageBox.information(self, 'Export success',
                                        'Exported:\n%s\ndisulfide_judgment.csv\n'
                                        'disulfide_report.md\nms2_annotated/'
                                        % report)
            except Exception as e:
                QMessageBox.critical(self, 'Export failed', str(e))

def main():
    if len(sys.argv) > 1 and sys.argv[1] == '--cli':
        if len(sys.argv) < 3:
            print('Usage: python disulfide_analyzer.py --cli <mgf file> '
                  '[peptide sequence] [output dir]')
            sys.exit(1)
        mgf = sys.argv[2]
        seq = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_SEQ
        out_dir = sys.argv[4] if len(sys.argv) > 4 else os.path.dirname(mgf) or '.'
        run_cli(mgf, seq, out_dir)
        return
    if not PYQT_OK:
        print('[Notice] PyQt6 not installed; falling back to CLI mode.')
        print('  Usage: python disulfide_analyzer.py --cli <mgf file> '
              '[peptide sequence] [output dir]')
        print('  Install GUI deps: pip install PyQt6 matplotlib')
        sys.exit(1)
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
