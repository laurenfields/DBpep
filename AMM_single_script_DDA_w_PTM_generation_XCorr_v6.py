# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:25:25 2023

@author: lawashburn
"""

import csv
import pandas as pd
import re
import os
from itertools import permutations
from Bio.SeqIO.FastaIO import SimpleFastaParser
import random
import collections
import time
import numpy as np
from scipy import signal
from datetime import datetime
start = time.time()

##User input##
output_parent_directory = r"C:\Users\lawashburn\Documents\DBpep_v2\XCorr_Opt\XCorr_validation\20230118\v5_no_drop_dups" #folder in which all output directories will be generated
db_path = r"C:\Users\lawashburn\Desktop\ALC50_Mass_Search_Files\duplicate_removed_crustacean_database_validated_formatted20220725.fasta" #database fasta path
base_file_path = r"C:\Users\lawashburn\Documents\DBpep_v2\results_log\formatted_MS2"

precursor_error_cutoff = 20 #ppm
fragment_error_cutoff = 0.02
precursor_charges = [1,2]
fragment_charges = [1,2]
h_mass = 1.00784
bin_size = 0.02
number_of_steps = 7
min_seq_coverage = 25
normalization=50

amidation = True
oxidation_M_status = True
pyroglu_E_status = True
pyroglu_Q_status = True
sulfo_Y_status = True
max_modifications = 1

raw_converter_path_input = [r"C:\Users\lawashburn\Documents\DBpep_v2\results_log\formatted_MS2\SG2_formatted.txt"]

##Definition storage

### generating database of mods selected ###
mods = []
mod_dict = {}
if oxidation_M_status == True:
    mods.append('M(Oxidation)')
    mod_dict['M'] = 'M(Oxidation)'
if pyroglu_E_status == True:
    mods.append('E(Pyro-glu)')
    mod_dict['E'] = 'E(Pyro-glu)'
if pyroglu_Q_status == True:
    mods.append('Q(Pyro-glu)')
    mod_dict['Q'] = 'Q(Pyro-glu)'
if sulfo_Y_status == True:
    mods.append('Y(Sulfo)')
    mod_dict['Y'] = 'Y(Sulfo)'
else:
    pass

modded_aas = []
for a in mods:
    modded_aas.append(a[0])

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
        for a in range(0,number_of_mods):
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
### start of monoisotopic mass calculator ###
def monoisotopic_mass_calculator(peptide_from_fasta):
        plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide_from_fasta)
        res_list_for_fragment = list_of_residues(peptide_from_fasta)
        mass_of_residues = []
        for residue in plain_peptide:
            residue_mass = aa_masses[residue]
            mass_of_residues.append(residue_mass)
        peptide_mass = (sum(mass_of_residues)) + check_termini(peptide_from_fasta) + check_PTM(peptide_from_fasta)
        mass_to_charge = (peptide_mass + (proton_mass * charge))/charge
        return mass_to_charge
### end of monoisotopic mass calculator ###

## Start of theoretical spectra generator ##
def theoretical_spectra_generator(peptide_to_check):
    ###pull theoretical masses###
    plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide_to_check) #removes any modification for mass calculations
    res_list_for_fragment = list_of_residues(peptide_to_check)
    mass_of_residues = []
    for residue in plain_peptide:
        residue_mass = aa_masses[residue]
        mass_of_residues.append(residue_mass)
    peptide_mass = (sum(mass_of_residues)) + check_termini(peptide_to_check) + check_PTM(peptide_to_check) #calculated MH mass
    mass_to_charge = (peptide_mass + (proton_mass * charge))/charge #translates the MH mass to m/z for each charge of interest

    num_ions = len(plain_peptide)-1 #number of expected fragment ions is one less than the number of AAs in the sequence

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
    return ion_report

### start of convert fasta to dataframe with monoisotopic masses ###
fasta_to_df = []
with open(db_path) as fasta_file:  # Will close handle cleanly
    for title, sequence in SimpleFastaParser(fasta_file):
        fasta_to_df.append(sequence)
um_database_list=fasta_to_df

final_seq_list = []
for b in fasta_to_df:
    aas_present = []
    aas_present_index = []
    aas_absent = []
    aas_absent_index = []
    
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
    
    if amidation == True:
        for m in all_modded_seqs_nodups:
            final_seq_list.append(m)
            amidated_pep = m + '(Amidated)'
            final_seq_list.append(amidated_pep)
    else:
        for m in all_modded_seqs_nodups:
            final_seq_list.append(m)
            
    ###start of mass calculations ###
    fasta_monoiso_mass = []
    pg_filtered_list = []
    for sequence in final_seq_list:
        if '(Pyro-glu)' in sequence:
            count_pg = sequence.count('(Pyro-glu')
            if count_pg == 1:
                if sequence[1] == '(':
                    pg_filtered_list.append(sequence)
                else:
                    pass
            else:
                    pass
        else:
            pg_filtered_list.append(sequence)
    finalized_mod_list = []  
    for seq in pg_filtered_list:
        if ')(' not in seq:
            finalized_mod_list.append(seq)
        else:
            pass
    for sequence in finalized_mod_list:
        mono_mass = monoisotopic_mass_calculator(sequence)
        fasta_monoiso_mass.append(mono_mass)
    
    db = pd.DataFrame()
    db['Sequence'] = finalized_mod_list
    db['Monoisotopic Mass'] = fasta_monoiso_mass
    
    ### end of convert fasta to dataframe with monoisotopic masses ###




###end of PTM generation###

for raw_converter_path in raw_converter_path_input:
    ##storage for results##
    final_peptide = []
    final_scan = []
    final_seq_coverage = []
    final_xcorr = []
    ##
    raw_file_sample_name1 = raw_converter_path.replace(base_file_path,'')
    raw_file_sample_name2 = raw_file_sample_name1.replace('_formatted','')
    raw_file_sample_name3 = raw_file_sample_name2.replace('\\','')
    sample_name = raw_file_sample_name3.replace('.txt','')
    print(sample_name)
    
    output_folder = output_parent_directory+'\\'+sample_name
    print(output_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    ### generate output folder ###
    peptide_report_output = output_folder+'\\peptide_reports'
    if not os.path.exists(peptide_report_output):
        os.makedirs(peptide_report_output)

    def param_log_export():
    
        out_path_report = 'Ouput directory path: ' + output_folder
        rawconverter_path_report = 'RawConverter output file path: ' + raw_converter_path    
        db_path_report = 'Database path: ' + db_path
        sample_ID_report = 'Sample name: ' + sample_name
        p_err_report = 'Precursor error threshold (ppm): ' + str(precursor_error_cutoff)
        f_err_report = 'Fragment error threshold (Da): ' + str(fragment_error_cutoff)
        max_p_charge_report = 'Maximum precursor charge: +' + str(precursor_charges[-1])
        max_f_charge_report = 'Maximum fragment charge: +' + str(fragment_charges[-1])
        bin_size_report = 'Bin size: ' + str(bin_size)
        no_steps_report = 'Number of steps: ' + str(number_of_steps)
        normalization_report = 'Normalized intensity : ' + str(normalization)    
        max_modifications_report = 'Maximum number of modifications : ' + str(max_modifications)    
        amidation_modifications_report = 'Amidation? : ' + str(amidation) 
        oxidation_M_status_modifications_report = 'Oxidation on M? : ' + str(oxidation_M_status) 
        pyroglu_E_status_modifications_report = 'Pyro-glu on E? : ' + str(pyroglu_E_status) 
        pyroglu_Q_status_modifications_report = 'Pyro-glu on Q? : ' + str(pyroglu_Q_status) 
        sulfo_Y_status_modifications_report = 'Sulfo on Y? : ' + str(sulfo_Y_status) 
        #elapsed_time_report = 'Time elapsed : ' + str(elapsed_time) 
        
        param_file_entries = [out_path_report,rawconverter_path_report,db_path_report,sample_ID_report,p_err_report,f_err_report,max_p_charge_report,max_f_charge_report,
                              bin_size_report,no_steps_report,normalization_report,max_modifications_report,amidation_modifications_report,oxidation_M_status_modifications_report,
                              pyroglu_E_status_modifications_report,pyroglu_Q_status_modifications_report,sulfo_Y_status_modifications_report]
    
        param_file_path = output_folder + '\\' + 'parameter_file.txt'
        with open(param_file_path,'a') as f:
            f.writelines('\n'.join(param_file_entries))
    
    ### start of importing spectra output from RawConverter ###
    raw_converter = pd.read_csv(raw_converter_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge","null"])
    raw_converter['charge'] = raw_converter['charge'].replace(to_replace=0,value=1) #assume z=0 is z=1
    
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
    
    ### end of importing spectra output from RawConverter ###
    
    if len(db)<1: #throws an error if database file is empty
        raise ValueError('Database file is empty')
    
    precursor_temp_cutoff = precursor_error_cutoff/100 #rough estimate of ppm to Da to minimize extra search space
    
    ### start of precursor AMM ###
    
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
    secondary_amm_prep_clean = secondary_amm_prep.drop(['MS2','scan_number','precursor_charge','null','resolution'],axis=1)
    secondary_amm_prep_clean = secondary_amm_prep_clean.rename(columns={'charge':'Fragment actual charge','m/z':'Fragment actual m/z'})
    
    candidate_sequences_raw = secondary_amm_prep_clean['Sequence'].values.tolist()
    
    candidate_sequences = []
    
    final_report_storage_scan = []
    final_report_storage_seq_coverage = []
    final_report_storage_peptide = []
    
    for m in candidate_sequences_raw:
        if m not in candidate_sequences:
            candidate_sequences.append(m)
    ### end of precursor AMM ###
    
    for peptide in candidate_sequences:
        ion_report = theoretical_spectra_generator(peptide)
    
        ### start of fragment AMM ###
        
        filtered_secondary_amm = secondary_amm_prep_clean[secondary_amm_prep_clean['Sequence'] == peptide] #filter amm from before for the sequence we are looking at
    
        filtered_secondary_amm['Fragment Actual Monoisotopic Mass'] = (filtered_secondary_amm['Fragment actual m/z'] * filtered_secondary_amm['Fragment actual charge']) - (h_mass*filtered_secondary_amm['Fragment actual charge'])     
    
        scans_present_raw = filtered_secondary_amm['Scan'].values.tolist()
        
        scans_present = []
        for z in scans_present_raw:
            if z not in scans_present:
                scans_present.append(z)
        
        for y in scans_present:
            scan_to_report = y
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
                final_seq_coverage.append(seq_coverage)
                final_report_storage_peptide.append(peptide)
                
                scan_filtered_score = merge_fragment_match_filtered[merge_fragment_match_filtered['Scan'] == y]
                scan_filtered_pep = scan_filtered_score['Sequence'].values.tolist()
            else:
                final_seq_coverage.append(0)
            actual_peak_intensities = raw_converter[raw_converter['scan_number'] == y]
            
            ###theoretical spectra for XCorr
            
            theoretical_mass_list_z1 = ion_report['Fragment theoretical Monoisotopic Mass'].values.tolist()
            
            theoretical_mass_list_z_all = []
            
            for zz in theoretical_mass_list_z1:
                for yy in fragment_charges:
                    new_mz = (zz + (h_mass*yy))/yy
                    theoretical_mass_list_z_all.append(new_mz)
            
            theo_mass_list = []
            theo_mass_list_intensity = []
            
            for dd in theoretical_mass_list_z_all:
                theo_mass_list.append(dd)
                theo_mass_list_intensity.append(50)
                
                theo_minus_flank1 = dd - bin_size
                theo_plus_flank1 = dd + bin_size
                
                theo_mass_list.append(theo_minus_flank1)
                theo_mass_list_intensity.append(25)
                
                theo_mass_list.append(theo_plus_flank1)
                theo_mass_list_intensity.append(25)
                
                theo_minus_flank2 = theo_minus_flank1 - bin_size
                theo_mass_list.append(theo_minus_flank2)
                theo_mass_list_intensity.append(10)
                
                theo_plus_flank2 = theo_plus_flank1 + bin_size
                theo_mass_list.append(theo_plus_flank2)
                theo_mass_list_intensity.append(10)
            
            theoretical_spectra_w_flanks = pd.DataFrame()
            theoretical_spectra_w_flanks['theoretical m/z'] = theo_mass_list
            theoretical_spectra_w_flanks['theoretical m/z'] = theoretical_spectra_w_flanks['theoretical m/z'].round(4)
            theoretical_spectra_w_flanks['theoretical intensity'] = theo_mass_list_intensity
            
            experimental_spectra = pd.DataFrame()
            experimental_spectra['experimental m/z'] = actual_peak_intensities['m/z']
            experimental_spectra['experimental non-normalized intensity'] = actual_peak_intensities['intensity']
            
            exp_max_int = experimental_spectra['experimental non-normalized intensity'].max()
            
            experimental_spectra['experimental intensity'] = (experimental_spectra['experimental non-normalized intensity']/(float(exp_max_int)))*50
            experimental_spectra['experimental intensity'] = experimental_spectra['experimental intensity'].round(2)
            
            range_eval = []
            
            range_eval.append(float(theoretical_spectra_w_flanks['theoretical m/z'].min()))
            range_eval.append(float(theoretical_spectra_w_flanks['theoretical m/z'].max()))
            range_eval.append(float(experimental_spectra['experimental m/z'].max()))
            range_eval.append(float(experimental_spectra['experimental m/z'].min()))
            xcorr_start_whole = (min(range_eval)) - 10
            xcorr_end_whole = (max(range_eval)) + 10
            xcorr_start = int(xcorr_start_whole)
            xcorr_end = int(xcorr_end_whole)
            
            interval_bins = []
            
            for aa in np.arange(xcorr_start,xcorr_end,bin_size):
                interval_bins.append(aa)
            
            bin_tolerance = bin_size
            
            xcorr_array = pd.DataFrame()
            xcorr_array['Bin #'] = interval_bins
            xcorr_array['Bin #'] = xcorr_array['Bin #'].round(2)
            
            xcorr_array = xcorr_array.sort_values(by='Bin #')
            theoretical_spectra_w_flanks = theoretical_spectra_w_flanks.sort_values(by='theoretical m/z')
            experimental_spectra = experimental_spectra.sort_values(by='experimental m/z')
            
            xcorr_array = pd.merge_asof(xcorr_array,experimental_spectra,left_on='Bin #',
                                                    right_on='experimental m/z', tolerance=bin_tolerance,allow_exact_matches=True)
            
            xcorr_array = pd.merge_asof(xcorr_array,theoretical_spectra_w_flanks,left_on='Bin #',
                                                    right_on='theoretical m/z', tolerance=bin_tolerance,allow_exact_matches=True)
            
            
            
            #xcorr_array = xcorr_array.drop_duplicates(subset=['theoretical m/z','experimental m/z','experimental intensity','theoretical intensity'])
            xcorr_array = xcorr_array.replace(np.nan, 0)

            folder_name = str(sample_name) + '_spectral_output'
            output_folder_spectra_compare = output_folder+'\\'+folder_name
            if not os.path.exists(output_folder_spectra_compare):
                os.makedirs(output_folder_spectra_compare)    
            theoretical_sample_title = str(peptide) + '_' + str(y) + '_theoretical_experimental_spectra.csv'
            file_path = output_folder_spectra_compare + '\\' + theoretical_sample_title
            with open(file_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    xcorr_array.to_csv(filec,index=False) 
            
            theoretical_sample_title2 = str(peptide) + '_' + str(y) + '_actual_peak_intensities.csv'
            file_path = output_folder_spectra_compare + '\\' + theoretical_sample_title2
            with open(file_path,'w',newline='') as filec:
                    writerc = csv.writer(filec)
                    actual_peak_intensities.to_csv(filec,index=False) 
            
            theoretical_array = xcorr_array['experimental intensity'].values.tolist()
            experimental_array = xcorr_array['theoretical intensity'].values.tolist()
            
            correlation_ctrl_array = signal.correlate(theoretical_array,experimental_array)
            correlation_score = correlation_ctrl_array.max()
            
            final_peptide.append(peptide)
            final_scan.append(y)
            final_xcorr.append(correlation_score)

    final_results_table = pd.DataFrame()
    final_results_table['Peptide'] = final_peptide
    final_results_table['Scan'] = final_scan
    final_results_table['Sequence coverage'] = final_seq_coverage
    final_results_table['XCorr'] = final_xcorr
    
    file_path = output_folder + '\\final_report.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            final_results_table.to_csv(filec,index=False) 