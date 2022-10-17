# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:38:03 2022

@author: lawashburn
"""

import csv
import pandas as pd
import re
import os

output_folder = r"C:\Users\lawashburn\Documents\DBpep_v2\DDA_PTM_results" #folder in which all output directories will be generated
raw_converter_path =  r"C:\Users\lawashburn\Documents\DBpep_v2\DDA_results\Brain1_2_formatted.txt" #path to formatted RawConverter output
db_path = r"D:\DIApep\20220816_test\database\cNP_db_20220725_database.csv"
sample_name = 'DIA_20220505_Exp15_20220816'
promex_cutoff = -10
precursor_error_cutoff = 20 #ppm
fragment_error_cutoff = 0.02
precursor_charges = [1,2,3,4,5,6,7,8]
fragment_charges = [1,2,3,4]
h_mass = 1.00784


peptide_report_output = output_folder+'\\peptide_reports'
if not os.path.exists(peptide_report_output):
    os.makedirs(peptide_report_output)

### Theoretical fragment calculator ###
proton_mass = 1.00727646688
charge = 1 #fragment charge


H = 1.0078250352
O = 15.99491463
C = 12.0000000
N = 14.003074
P = 30.973762
S = 31.9720707

aa_masses = {
    'G' : C*2  + H*3  + N   + O,
    'A' : C*3  + H*5  + N   + O,
    'S' : C*3  + H*5  + N   + O*2,
    'P' : C*5  + H*7  + N   + O,
    'V' : C*5  + H*9  + N   + O,
    'T' : C*4  + H*7  + N   + O*2,
    'C' : C*3  + H*5  + N   + O   + S,
    'L' : C*6  + H*11 + N   + O,
    'I' : C*6  + H*11 + N   + O,
    'N' : C*4  + H*6  + N*2 + O*2,
    'D' : C*4  + H*5  + N   + O*3,
    'Q' : C*5  + H*8  + N*2 + O*2,
    'K' : C*6  + H*12 + N*2 + O ,
    'E' : C*5  + H*7  + N   + O*3 ,
    'M' : C*5  + H*9  + N   + O   + S ,
    'H' : C*6  + H*7  + N*3 + O ,
    'F' : C*9  + H*9  + N   + O ,
    'R' : C*6  + H*12 + N*4 + O ,
    'Y' : C*9  + H*9  + N   + O*2 ,
    'W' : C*11 + H*10 + N*2 + O ,
    'O' : C*5  + H*12 + N*2 + O*2,
    'C(Pyro-glu)' : C*3  + H * 2 + O + S,
    'Q(Pyro-glu)' : C*5  + H*5  + N + O*2,
    'E(Pyro-glu)' : C*5  + H*4 + O*3,
    'M(Oxidation)' : C*5  + H*9  + N   + O*2   + S,
    'Y(Sulfo)' :  C*9  + H*9  + N   + O*5 + S
    }

#database_input = pd.read_csv(database_path)
#sequence = database_input['Sequence'].values.tolist()

termini = {'Standard' : H * 2 + O,
'(Amidated)' : N + H * 3}

PTMs = {'C(Pyro-glu)' : H * -3 - N,
        'Q(Pyro-glu)' : H * -3 - N,
        'E(Pyro-glu)' : H * -3 - N,
        'M(Oxidation)' : O,
        'Y(Sulfo)' : S + O * 3
        }

adducts = {
    'H2O' : H * 2 + O,
    'NH3' : N + H * 3}

def check_termini_Key(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        return termini['Standard']

def check_PTM_Key(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        return 0

def check_termini(pot_mod_peptide):
    if '(' in pot_mod_peptide:
        term_start = (pot_mod_peptide.rindex('('))
        termini_ID = pot_mod_peptide[term_start:]
        termini_mass_change = check_termini_Key(termini, termini_ID)
        return termini_mass_change
    else:
        return termini['Standard']

def check_PTM(pot_mod_peptide):
    number_of_mods = pot_mod_peptide.count('(')
    if number_of_mods > 0:
        current_peptide = []
        mass_change_collection = []
        current_peptide.append(pot_mod_peptide)
        for a in range(1,number_of_mods):
            peptide_less_mod = current_peptide[-1]
            ptm_start = (peptide_less_mod.index('('))-1
            ptm_end = (peptide_less_mod.index(')'))+1
            ptm_ID = peptide_less_mod[ptm_start:ptm_end]
            ptm_mass_change = check_PTM_Key(PTMs, ptm_ID)
            mass_change_collection.append(ptm_mass_change)
            peptide_less_mod2 = peptide_less_mod[:ptm_start] + peptide_less_mod[ptm_end:]
            current_peptide.append(peptide_less_mod2)
            
        ptm_mass_change = sum(mass_change_collection)
        return ptm_mass_change
    else:
        ptm_mass_change = 0
        return ptm_mass_change


def list_of_residues(pot_mod_peptide):

    list_of_res = []
    
    pep_update = []
    pep_update.append(pot_mod_peptide)
    
    no_mods = pot_mod_peptide.count('(')
    if no_mods > 1:
        for c in range(0,no_mods+1):
            pep_of_interest = pep_update[-1]
            if '(' in pep_of_interest:
                first_ptm_start = pep_of_interest.index('(')
                first_ptm_end = pep_of_interest.index(')')

                first_residues = pep_of_interest[:(first_ptm_start-1)]
                for a in first_residues:
                    list_of_res.append(a)
                ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                list_of_res.append(ptm_residue)
                remaining_pep = pep_of_interest[(first_ptm_end+1):]
                pep_update.append(remaining_pep)
                
            else:
                for d in pep_of_interest:
                    list_of_res.append(d)
    elif no_mods == 1:
        for c in range(0,1):
            pep_of_interest = pep_update[-1]
            if '(' in pep_of_interest:
                first_ptm_start = pep_of_interest.index('(')
                first_ptm_end = pep_of_interest.index(')')
                
                if first_ptm_start == 1:
                    ptm_residue =  pep_of_interest[0] + (pep_of_interest[(first_ptm_start):(first_ptm_end+1)])
                    list_of_res.append(ptm_residue)
                    remaining_pep = pep_of_interest[(first_ptm_end+1):]
                    for d in remaining_pep:
                        list_of_res.append(d)

                if first_ptm_start != 1:
                    first_residues = pep_of_interest[:(first_ptm_start-1)]
                    for a in first_residues:
                        list_of_res.append(a)
                    ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                    list_of_res.append(ptm_residue)
                    remaining_pep = pep_of_interest[(first_ptm_end+1):]
                    for d in remaining_pep:
                        list_of_res.append(d)
                
            else:
                for d in pep_of_interest:
                    list_of_res.append(d)
    
    elif no_mods == 0:
        for c in pot_mod_peptide:
            list_of_res.append(c)
    return list_of_res

### end of theoretical fragment calculator ###

raw_converter = pd.read_csv(raw_converter_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge","null"])

raw_converter = raw_converter[raw_converter['charge'] != 0]

db = pd.read_csv(db_path)

precursor_mz = []
precursor_z = []
precursor_scan = []

rawconv_mz = raw_converter['MS2'].values.tolist()
rawconv_z = raw_converter['precursor_charge'].values.tolist()
rawconv_scan = raw_converter['scan_number'].values.tolist()

for b in rawconv_mz:
    precursor_mz.append(b)

for d in rawconv_z:
    precursor_z.append(d)

for e in rawconv_scan:
    precursor_scan.append(e)

exp_precursor = pd.DataFrame()
exp_precursor['Precursor actual m/z'] = precursor_mz
exp_precursor['Precursor actual charge'] = precursor_z
exp_precursor['Precursor scan'] = precursor_scan


exp_precursor = exp_precursor.drop_duplicates() 
exp_precursor['Monoisotopic Mass'] =  ((exp_precursor['Precursor actual m/z']) * (exp_precursor['Precursor actual charge']))-(h_mass*(exp_precursor['Precursor actual charge']))

if len(db)<1: #throws an error if database file is empty
    raise ValueError('Database file is empty')

precursor_temp_cutoff = precursor_error_cutoff/1000

precursor_amm_actual_mz = []
precursor_amm_actual_z = []
precursor_amm_actual_scan = []
precursor_amm_actual_monoiso = []
precursor_amm_theoretical_sequence = []
precursor_amm_theoretical_monoiso = []
precursor_amm_err = []

for a in precursor_charges:
    db_sorted = db.sort_values(by = 'Monoisotopic Mass')
    db_sorted = db_sorted.rename(columns={'Monoisotopic Mass':'Theoretical Monoisotopic Mass'})

    exp_precursor_sorted = exp_precursor.sort_values(by = 'Monoisotopic Mass')
    exp_precursor_sorted = exp_precursor_sorted.rename(columns={'Monoisotopic Mass':'Actual Monoisotopic Mass'})

    exp_precursor_z_filter = exp_precursor_sorted[exp_precursor_sorted['Precursor actual charge'] == a]

    merge_match = pd.merge_asof(exp_precursor_z_filter,db_sorted, left_on='Actual Monoisotopic Mass', right_on='Theoretical Monoisotopic Mass',
                                tolerance = precursor_temp_cutoff, allow_exact_matches=True)

    merge_match_filtered = merge_match.dropna(subset=['Sequence','Theoretical Monoisotopic Mass'])
    merge_match_filtered['Precursor error (ppm)'] = ((abs((merge_match_filtered['Theoretical Monoisotopic Mass'])-(merge_match_filtered['Actual Monoisotopic Mass'])))/
                                                     (merge_match_filtered['Theoretical Monoisotopic Mass'])) * 1E6
    merge_match_filtered2 = merge_match_filtered[merge_match_filtered['Precursor error (ppm)'] <= precursor_error_cutoff]

    actual_mz = merge_match_filtered2['Precursor actual m/z'].values.tolist()
    actual_z = merge_match_filtered2['Precursor actual charge'].values.tolist()
    actual_scan = merge_match_filtered2['Precursor scan'].values.tolist()
    actual_monoiso = merge_match_filtered2['Actual Monoisotopic Mass'].values.tolist()
    theoretical_sequence = merge_match_filtered2['Sequence'].values.tolist()
    theoretical_monoiso = merge_match_filtered2['Theoretical Monoisotopic Mass'].values.tolist()
    err = merge_match_filtered2['Precursor error (ppm)'].values.tolist()
    
    for f in actual_mz:
        precursor_amm_actual_mz.append(f)
    for g in actual_z:
        precursor_amm_actual_z.append(g)
    for h in actual_scan:
        precursor_amm_actual_scan.append(h)
    for i in actual_monoiso:
        precursor_amm_actual_monoiso.append(i)
    for j in theoretical_sequence:
        precursor_amm_theoretical_sequence.append(j)
    for k in theoretical_monoiso:
        precursor_amm_theoretical_monoiso.append(k)
    for l in err:
        precursor_amm_err.append(l)

precursor_amm_results = pd.DataFrame()
precursor_amm_results['Precursor Actual m/z'] = precursor_amm_actual_mz
precursor_amm_results['Precursor Actual z'] = precursor_amm_actual_z
precursor_amm_results['Scan'] = precursor_amm_actual_scan
precursor_amm_results['Precursor Actual Monoisotopic'] = precursor_amm_actual_monoiso
precursor_amm_results['Sequence'] = precursor_amm_theoretical_sequence
precursor_amm_results['Precursor Theoretical Monoisotopic'] = precursor_amm_theoretical_monoiso
precursor_amm_results['Precursor error (ppm)'] = precursor_amm_err

secondary_amm_prep = pd.merge(raw_converter,precursor_amm_results,left_on=['MS2','scan_number','precursor_charge'], right_on=['Precursor Actual m/z','Scan','Precursor Actual z'])
secondary_amm_prep_clean = secondary_amm_prep.drop(['MS2','scan_number','precursor_charge','null','resolution','intensity'],axis=1)
secondary_amm_prep_clean = secondary_amm_prep_clean.rename(columns={'charge':'Fragment actual charge','m/z':'Fragment actual m/z'})

# =============================================================================
# file_path = output_folder + '\\secondary_amm_prep.csv'
# with open(file_path,'w',newline='') as filec:
#         writerc = csv.writer(filec)
#         secondary_amm_prep_clean.to_csv(filec,index=False)   
# =============================================================================

candidate_sequences_raw = secondary_amm_prep_clean['Sequence'].values.tolist()

candidate_sequences = []

final_report_storage_scan = []
final_report_storage_seq_coverage = []
final_report_storage_peptide = []

for m in candidate_sequences_raw:
    if m not in candidate_sequences:
        candidate_sequences.append(m)

for peptide in candidate_sequences:
        ###pull theoretical masses###
        plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide)
        res_list_for_fragment = list_of_residues(peptide)
        mass_of_residues = []
        for residue in plain_peptide:
            residue_mass = aa_masses[residue]
            mass_of_residues.append(residue_mass)
        peptide_mass = (sum(mass_of_residues)) + check_termini(peptide) + check_PTM(peptide)
        mass_to_charge = (peptide_mass + (proton_mass * charge))/charge

        num_ions = len(plain_peptide)-1

        b_ions = []
        b_ion_name = []
        
        y_ions = []
        y_ion_name = []
        
        for a in range(0,num_ions):
            
            residue_identity = res_list_for_fragment[a]
            if len(b_ions) == 0:
                ion_mass = aa_masses[residue_identity]
                ion_mz = ion_mass + proton_mass
                b_ions.append(ion_mz)
                ion_name = 'b' + str(a+1)
                b_ion_name.append(ion_name)
            
            elif len(b_ions) > 0:
                ion_mass = (aa_masses[residue_identity]) + b_ions[-1]
                b_ions.append(ion_mass)
                ion_name = 'b' + str(a+1)
                b_ion_name.append(ion_name)
        
        for b in (range(0,num_ions)):
            residue_identity = res_list_for_fragment[b]
            if len(y_ions) == 0:
                ion_mass = mass_to_charge - aa_masses[residue_identity]
                y_ions.append(ion_mass)
                ion_name = 'y' + str((num_ions-b))
                y_ion_name.append(ion_name)
            elif len(y_ions) > 0:
                ion_mass = y_ions[-1] - aa_masses[residue_identity]
                y_ions.append(ion_mass)
                ion_name = 'y' + str((num_ions-b))
                y_ion_name.append(ion_name)
            
        y_ion_name.append('MH')
        y_ions.append(mass_to_charge)
            
        b_ions_report = pd.DataFrame()
        b_ions_report['ion'] = b_ion_name
        b_ions_report['mass'] = b_ions
        
        b_ions_water_adduct = pd.DataFrame()
        b_ions_water_adduct['ion'] = b_ions_report['ion'] + '-H2O'
        b_ions_water_adduct['mass'] = b_ions_report['mass'] - adducts['H2O']
        
        b_ions_ammonia_adduct = pd.DataFrame()
        b_ions_ammonia_adduct['ion'] = b_ions_report['ion'] + '-NH3'
        b_ions_ammonia_adduct['mass'] = b_ions_report['mass'] - adducts['NH3']
        
        y_ions_report = pd.DataFrame()
        y_ions_report['ion'] = y_ion_name
        y_ions_report['mass'] = y_ions
        
        y_ions_ammonia_adduct = pd.DataFrame()
        y_ions_ammonia_adduct['ion'] = y_ions_report['ion'] + '-NH3'
        y_ions_ammonia_adduct['mass'] = y_ions_report['mass'] - adducts['NH3']

        y_ions_water_adduct = pd.DataFrame()
        y_ions_water_adduct['ion'] = y_ions_report['ion'] + '-H2O'
        y_ions_water_adduct['mass'] = y_ions_report['mass'] - adducts['H2O']
        
        ion_report = pd.DataFrame()
        ion_report = ion_report.append(b_ions_report)
        ion_report = ion_report.append(y_ions_report)
        ion_report = ion_report.append(b_ions_water_adduct)
        ion_report = ion_report.append(b_ions_ammonia_adduct)
        ion_report = ion_report.append(y_ions_ammonia_adduct)
        ion_report = ion_report.append(y_ions_water_adduct)
        ion_report = ion_report.drop_duplicates()
        ion_report = ion_report.rename(columns={'mass':'Fragment theoretical Monoisotopic Mass'})
        ### end of pulling theoretical fragments ###

        filtered_secondary_amm = secondary_amm_prep_clean[secondary_amm_prep_clean['Sequence'] == peptide] #filter amm from before for the sequence we are looking at

        filtered_secondary_amm['Fragment Actual Monoisotopic Mass'] = (filtered_secondary_amm['Fragment actual m/z'] * filtered_secondary_amm['Fragment actual charge']) - (h_mass*filtered_secondary_amm['Fragment actual charge'])
        
        scans_present_raw = filtered_secondary_amm['Scan'].values.tolist()
        
        scans_present = []
        for z in scans_present_raw:
            if z not in scans_present:
                scans_present.append(z)
        
        for y in scans_present:
            final_report_storage_scan.append(y)
            scans_filtered_secondary_amm = filtered_secondary_amm[filtered_secondary_amm['Scan'] == y] #filter to look at just one scan
            scans_filtered_secondary_amm = scans_filtered_secondary_amm.sort_values(by='Fragment Actual Monoisotopic Mass')
            ion_report = ion_report.sort_values(by='Fragment theoretical Monoisotopic Mass')
            prelim_fragment_matches = pd.merge_asof(scans_filtered_secondary_amm,ion_report,left_on='Fragment Actual Monoisotopic Mass',
                                                    right_on='Fragment theoretical Monoisotopic Mass', tolerance=fragment_error_cutoff,allow_exact_matches=True)
            merge_fragment_match_filtered = prelim_fragment_matches.dropna(subset=['ion','Fragment theoretical Monoisotopic Mass'])
            
            merge_fragment_match_filtered['Fragment error (Da)'] = merge_fragment_match_filtered['Fragment Actual Monoisotopic Mass'] - merge_fragment_match_filtered['Fragment theoretical Monoisotopic Mass']
            merge_fragment_match_filtered = merge_fragment_match_filtered[merge_fragment_match_filtered['Fragment error (Da)'] <= fragment_error_cutoff]
            
            if len(merge_fragment_match_filtered)>0:
                file_path = peptide_report_output+'\\'+peptide+'_'+str(y)+'_prelim_fragment_matches.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        merge_fragment_match_filtered.to_csv(filec,index=False) 

            merge_fragment_match_filtered['ion_format'] = merge_fragment_match_filtered['ion'].str.extract('(\d+)', expand=False)
            ion_report['ion_format'] = ion_report['ion'].str.extract('(\d+)', expand=False)
            
            merge_fragment_match_filtered_ion = merge_fragment_match_filtered.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
            ion_report_filtered_ion = ion_report.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
            
            seq_coverage = (len(merge_fragment_match_filtered_ion)/len(ion_report_filtered_ion))*100
            final_report_storage_seq_coverage.append(seq_coverage)
            final_report_storage_peptide.append(peptide)
#%%
final_report = pd.DataFrame()
final_report['Neuropeptide'] = final_report_storage_peptide
final_report['Scan'] = final_report_storage_scan
final_report['Sequence Coverage'] = final_report_storage_seq_coverage
final_report = final_report[final_report['Sequence Coverage'] > 0]

final_report = final_report.sort_values(by='Sequence Coverage',ascending=False)
final_report_seq = final_report.drop_duplicates(subset=['Neuropeptide'])
final_report_scan = final_report_seq.drop_duplicates(subset=['Scan'])

file_path = output_folder + '\\final_report.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        final_report_scan.to_csv(filec,index=False) 
#%%
# =============================================================================
# 
# final_report_max = pd.merge(final_report,final_report_scan,on=['Neuropeptide','Scan'],how='left')
# print(final_report_max)
# file_path = output_folder + '\\final_report_max.csv'
# with open(file_path,'w',newline='') as filec:
#         writerc = csv.writer(filec)
#         final_report_max.to_csv(filec,index=False) 
# missing_scans = final_report_max[final_report_max['Sequence Coverage_y'].isna()]
# 
# banned_scans = final_report_scan['Scan'].values.tolist()
# banned_seqs = final_report_scan['Neuropeptide'].values.tolist()
# 
# print(missing_scans)
# missing_scans_filter1 = missing_scans[~missing_scans['Scan'].isin(banned_scans)]
# missing_scans_filter2 = missing_scans_filter1[~missing_scans_filter1['Neuropeptide'].isin(banned_seqs)]
# print(missing_scans_filter2)
# =============================================================================
# =============================================================================
# 
# 
# sequences = db['Sequence'].values.tolist()
# sequence_store = []
# charge_store = []
# target_mz_store = []
# query_mz_store = []
# ppm_err_store = []
# 
# query = topFD_combined_out
# 
# for b in sequences:
#     db2 = db[db['Sequence'] == b]
#     for c in charges:
#         target_mz = db2[str(c)]
#         query2 = query[query['Precursor actual charge'] == c]
#         monomz = query2['Precursor actual m/z'].values.tolist()
#         for d in target_mz:
#             for e in monomz:
#                 ppm_err = ((abs(d-e))/d) * 1E6
#                 
#                 sequence_store.append(b)
#                 charge_store.append(c)
#                 target_mz_store.append(d)
#                 query_mz_store.append(e)
#                 ppm_err_store.append(ppm_err)
#                 
# results_df = pd.DataFrame()
# results_df['Sequence'] = sequence_store
# results_df['Precursor actual charge'] = charge_store  
# results_df['Theoretical precursor m/z'] = target_mz_store
# results_df['Precursor actual m/z'] = query_mz_store
# results_df['Precursor error (ppm)'] = ppm_err_store
# results_df = results_df[results_df['Precursor error (ppm)'] <= error_cutoff]
# results_df = results_df.sort_values(by='Precursor error (ppm)')
# results_df = results_df.drop_duplicates(subset=['Sequence','Precursor actual charge'],keep='first').reset_index(drop=True)
# 
# complete_results_df = pd.merge(results_df, query,  how='left', left_on=['Precursor actual m/z','Precursor actual charge'], right_on = ['Precursor actual m/z','Precursor actual charge'])
# 
# complete_results_df = complete_results_df.sort_values(by='Precursor error (ppm)')
# complete_results_df = complete_results_df.drop_duplicates(subset=['Sequence','Precursor actual charge'],keep='first').reset_index(drop=True)
# 
# file_path = prelim_AMM_out_file_name_path
# with open(file_path,'w',newline='') as filec:
#         writerc = csv.writer(filec)
#         complete_results_df.to_csv(filec,index=False)        
# prelim_AMM_out = complete_results_df
# 
# =============================================================================
