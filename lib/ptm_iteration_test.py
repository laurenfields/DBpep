# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 15:31:32 2022

@author: lawashburn
"""

from itertools import permutations
import pandas as pd
import time
start = time.time()
#sequence = ['ELMAPEN']
sequence = ['ELEEEEEMAPEN']
amidation = True
max_modifications = 3

mods = ['E(Pyro-glu','M(Oxidation)']

mod_dict = {'E':'E(Pyro-glu)',
            'M':'M(Oxidation)'}

modded_aas = []

for a in mods:
    modded_aas.append(a[0])

aas_present = []
aas_present_index = []

aas_absent = []
aas_absent_index = []

for b in sequence:
    seq_len = len(b)
    for c in range(0,seq_len):
        aa_ID = b[c]
        if aa_ID in modded_aas:
            aas_present.append(aa_ID)
            aas_present_index.append(str(c))
            aas_absent.append(aa_ID)
            aas_absent_index.append(str(c))
        else:
            aas_absent.append(aa_ID)
            aas_absent_index.append(str(c))

number_mods = len(aas_present)

seq_log = pd.DataFrame()
seq_log['Original Residue'] = aas_present
seq_log['Index'] = aas_present_index
seq_log['Res+Index'] = seq_log[['Original Residue', 'Index']].apply(lambda x: ''.join(x), axis=1)
seq_w_ind = seq_log['Res+Index'].values.tolist()

seq_absent_log = pd.DataFrame()
seq_absent_log['Original Residue'] = aas_absent
seq_absent_log['Index'] = aas_absent_index



for d in range(0,number_mods):
    seq_w_ind.append('X'+(str(d)))

aa_combination_no_dups = []

if number_mods > max_modifications:
    aa_combinations = permutations(seq_w_ind,max_modifications)
    aa_combo_list = list(aa_combinations)
    for a in aa_combo_list:
            if a not in aa_combination_no_dups:
                aa_combination_no_dups.append(a)
else:
    aa_combinations = permutations(seq_w_ind,number_mods)
    aa_combo_list = list(aa_combinations)
    for a in aa_combo_list:
            if a not in aa_combination_no_dups:
                aa_combination_no_dups.append(a)

restyle_aa_list = []
restyle_aa_ind_list = []
for z in aa_combination_no_dups:
    restyle_aa = []
    restyle_aa_ind = []
    for y in z:
        restyle_aa.append(y[0])
        restyle_aa_ind.append(y[1:])
    restyle_aa_list.append(restyle_aa)
    restyle_aa_ind_list.append(restyle_aa_ind)

unmodded_seq = pd.DataFrame()
unmodded_seq['Original Residue'] = aas_absent
unmodded_seq['Index'] = aas_absent_index

all_modded_seqs = []
for k in range(0,len(restyle_aa_list)):
    new_seq_log = pd.DataFrame()
    new_seq_log['New Residue'] = restyle_aa_list[k]
    new_seq_log['Index'] = restyle_aa_ind_list[k]
    new_seq_log_filtered = new_seq_log[new_seq_log['New Residue'] != 'X']
    ptm_applied_seq_log = new_seq_log_filtered.replace({'New Residue':mod_dict})
    
    merge_seq_log = unmodded_seq.merge(ptm_applied_seq_log, on='Index', how='left')
    ptm_applied_seq_log = merge_seq_log.sort_values(by='Index')
    ptm_applied_seq_log['New Residue'] = ptm_applied_seq_log['New Residue'].replace('', pd.NA).fillna(ptm_applied_seq_log['Original Residue'])

    modified_seq = ptm_applied_seq_log['New Residue'].values.tolist()
    modified_seq_format = "".join([str(item) for item in modified_seq])
    all_modded_seqs.append(modified_seq_format)

all_modded_seqs_nodups = []
for l in all_modded_seqs:
    if l not in all_modded_seqs_nodups:
        all_modded_seqs_nodups.append(l)

terminal_seqs = []
if amidation == True:
    for m in all_modded_seqs_nodups:
        terminal_seqs.append(m)
        amidated_pep = m + '(Amidated)'
        terminal_seqs.append(amidated_pep)
else:
    for m in all_modded_seqs_nodups:
        terminal_seqs.append(m)

print(terminal_seqs)
print(len(terminal_seqs))
end = time.time()
print('Analysis complete. Time elapsed:',(end - start),'s')