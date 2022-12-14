# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 09:50:47 2022

@author: lawashburn
"""

from tkinter import *
from PIL import ImageTk, Image
from tkinter.filedialog import askopenfilename
import sys
from itertools import islice
import subprocess
from subprocess import Popen, PIPE
from textwrap import dedent
from tkinter import messagebox
from tkinter import filedialog
import os
import pickle
from tkinter import ttk
import threading
import traceback
import time
import psutil

def hide_all_frames():
    logo_frame.pack_forget()
    about_click_frame.pack_forget()
    #modules_click_frame.pack_forget()
    files_click_frame.pack_forget()
    run_click_frame.pack_forget()
    help_click_frame.pack_forget()
### Definitions for file selection/user input decisions ###
def set_path_rawconverter_file_field():
    path_rawconverter_file = askopenfilename(filetypes=[("MS2 Files","*.MS2")],multiple=False) 
    input_rawconverter_file_text.set(path_rawconverter_file)
def get_rawconverter_file_path(): 
    """ Function provides the database full file path."""
    return path_rawconverter_file
def set_path_database_field():
    path_db = askopenfilename(filetypes=[("FASTA Files","*.fasta")]) 
    input_db_text.set(path_db)
def get_database_path(): 
    """ Function provides the database full file path."""
    return path_db
def set_path_out_dir_field():
    path_out_dir = filedialog.askdirectory() 
    out_dir_var.set(path_out_dir)
def Run_Search():
    try:
        subprocess.run(['python.exe','AMM_single_script_DDA_w_PTM_generation_XCorr_v2.py'],check=True)
        messagebox.showinfo("Message","Analysis results have been exported")
    except FileNotFoundError:
        messagebox.showinfo('Error','One or more input entries is invalid. See documentation.')
    except Exception:
        messagebox.showerror('Error','One or more input entries is invalid. See documentation.')
    
    # time.sleep(5)
    # os.system('python AMM_single_script_DDA_w_PTM_generation_XCorr_v2.py')
def param_log_export():

    sample_name_var_report = sample_name_var.get()
    sample_name_report = 'Sample name: ' + sample_name_var_report
    
    input_rawconverter_file_text_var_report = input_rawconverter_file_text.get()
    input_rawconverter_file_text_report = 'Spectra path: ' + input_rawconverter_file_text_var_report
    
    input_db_text_var_report = input_db_text.get()
    input_db_text_report = 'Database path: ' + input_db_text_var_report
    
    prec_err_var_report = prec_err_var.get()
    prec_err_report = 'Precursor error threshold (ppm): ' + str(prec_err_var_report)
    
    frag_err_var_report = frag_err_var.get()
    frag_err_report = 'Fragment error threshold (Da): ' + str(frag_err_var_report)
    
    prec_z_var_report = prec_z_var.get()
    prec_z_report = 'Maximum precursor charge: ' + str(prec_z_var_report)
    
    frag_z_var_report = frag_z_var.get()
    frag_z_report = 'Maximum fragment charge: ' + str(frag_z_var_report)
    
    amid_check_var_report = amid_check_var.get()
    if amid_check_var_report == 0:
        amid_check_report = 'C-terminal amidation: False'
    if amid_check_var_report == 1:
        amid_check_report = 'C-terminal amidation: True'
    
    oxM_check_var_report = oxM_check_var.get()
    if oxM_check_var_report == 0:
        oxM_check_report = 'Oxidation on M: False'
    if oxM_check_var_report == 1:
        oxM_check_report = 'Oxidation on M: True'

    PG_E_check_var_report = PG_E_check_var.get()
    if PG_E_check_var_report == 0:
        PG_E_check_report = 'Pyro-glu on E: False'
    if PG_E_check_var_report == 1:
        PG_E_check_report = 'Pyro-glu on E: True'
    
    PG_Q_check_var_report = PG_Q_check_var.get()
    if PG_Q_check_var_report == 0:
        PG_Q_check_report = 'Pyro-glu on Q: False'
    if PG_Q_check_var_report == 1:
        PG_Q_check_report = 'Pyro-glu on Q: True'
    
    sulf_Y_check_var_report = sulf_Y_check_var.get()
    if sulf_Y_check_var_report == 0:
        sulf_Y_check_report = 'Sulfation on Y: False'
    if sulf_Y_check_var_report == 1:
        sulf_Y_check_report = 'Sulfation on Y: True'

    max_mods_var_report = max_mods_var.get()
    max_mods_report = 'Maximum PTMs per neuropeptide: ' + str(max_mods_var_report)
    
    FDR_thresh_var_report = FDR_thresh_var.get()
    FDR_thresh_report = 'FDR threshold (%): ' + str(FDR_thresh_var_report)
    
    bin_size_var_report = bin_size_var.get()
    bin_size_report = 'Bin size: ' + str(bin_size_var_report)
    
    bin_steps_var_report = bin_steps_var.get()
    bin_steps_report = 'Number of bin steps: ' + str(bin_steps_var_report)
    
    seq_cov_thresh_var_report = seq_cov_thresh_var.get()
    seq_cov_thresh_report = 'Minimum sequence coverage threshold (%): ' + str(seq_cov_thresh_var_report)
    
    out_dir_var_report = out_dir_var.get()
    out_dir_report = 'Output directory path: ' + out_dir_var_report

    param_file_entries = [sample_name_report,input_rawconverter_file_text_report,input_db_text_report,prec_err_report,frag_err_report,
                          prec_z_report,frag_z_report,amid_check_report,oxM_check_report,PG_E_check_report,PG_Q_check_report,
                          sulf_Y_check_report,max_mods_report,FDR_thresh_report,bin_size_report,bin_steps_report,seq_cov_thresh_report,out_dir_report]

    param_file_path = out_dir_var_report + '\\' + sample_name_var_report + '_parameter_file.txt'
    with open(param_file_path,'a') as f:
        f.writelines('\n'.join(param_file_entries))

    sample_name_pick_path = 'sample_name.pkl'
    rawconverter_path_pick_path = 'rawconverter_path.pkl'
    db_path_pick_path = 'db_path.pkl'
    prec_err_pick_path = 'prec_err.pkl'
    frag_err_pick_path = 'frag_err.pkl'
    prec_z_pick_path = 'prec_z.pkl'
    frag_z_pick_path = 'frag_z.pkl'
    amid_pick_path = 'amidation.pkl'
    pyroglu_E_pick_path = 'pyroglu_E.pkl'
    pyroglu_Q_pick_path = 'pyroglu_Q.pkl'    
    oxo_M_pick_path = 'oxidation_M.pkl'    
    sulfo_Y_pick_path = 'sulfation_Y.pkl' 
    max_PTMs_pick_path = 'max_PTMs.pkl'  
    fdr_pick_path = 'FDR.pkl' 
    bin_size_pick_path = 'bin_size.pkl'
    bin_step_pick_path = 'bin_steps.pkl'
    seq_cov_pick_path = 'sequence_coverage.pkl'  
    output_dir_pick_path = 'output_directory.pkl'    
    
     #A new file will be created
    with open(sample_name_pick_path, 'wb') as file_e:
        pickle.dump(sample_name_var_report, file_e)
    with open(rawconverter_path_pick_path, 'wb') as file_f:
        pickle.dump(input_rawconverter_file_text_var_report, file_f)
    with open(db_path_pick_path, 'wb') as file_g:
        pickle.dump(input_db_text_var_report, file_g)      
    with open(prec_err_pick_path, 'wb') as file_h:
        pickle.dump(prec_err_var_report, file_h)        
    with open(frag_err_pick_path, 'wb') as file_i:
        pickle.dump(frag_err_var_report, file_i)        
    with open(prec_z_pick_path, 'wb') as file_j:
        pickle.dump(prec_z_var_report, file_j)        
    with open(frag_z_pick_path, 'wb') as file_k:
        pickle.dump(frag_z_var_report, file_k)
    with open(amid_pick_path, 'wb') as file_l:
        pickle.dump(amid_check_var_report, file_l)        
    with open(pyroglu_E_pick_path, 'wb') as file_m:
        pickle.dump(PG_E_check_var_report, file_m)     
    with open(pyroglu_Q_pick_path, 'wb') as file_n:
        pickle.dump(PG_Q_check_var_report, file_n)  
    with open(oxo_M_pick_path, 'wb') as file_o:
        pickle.dump(oxM_check_var_report, file_o)        
    with open(sulfo_Y_pick_path, 'wb') as file_p:
        pickle.dump(sulf_Y_check_var_report, file_p)       
    with open(max_PTMs_pick_path, 'wb') as file_q:
        pickle.dump(max_mods_var_report, file_q)       
    with open(fdr_pick_path, 'wb') as file_q:
        pickle.dump(FDR_thresh_var_report, file_q) 
    with open(bin_size_pick_path, 'wb') as file_q:
        pickle.dump(bin_size_var_report, file_q) 
    with open(bin_step_pick_path, 'wb') as file_q:
        pickle.dump(bin_steps_var_report, file_q) 
    with open(seq_cov_pick_path, 'wb') as file_s:
        pickle.dump(seq_cov_thresh_var_report, file_s)        
    with open(output_dir_pick_path, 'wb') as file_t:
        pickle.dump(out_dir_var_report, file_t)    
def combine_funcs(*funcs):
    def combined_func(*args, **kwargs):
        for f in funcs:
            f(*args, **kwargs)
    return combined_func

### GUI appearance settings ###
Font_tuple = ("Corbel Light", 20)
root = Tk()
root.title('Neuropeptide Database Search')
root.iconbitmap(r"hypep_icon.ico")
root.geometry('525x670')

### Variable storage ###
sample_name_var = StringVar()
sample_name_var.set('test')

input_rawconverter_file_text =  StringVar()
input_rawconverter_file_text.set('C:/Users/lawashburn/Documents/DBpep_v2/results_log/2021_0817_Brain_3.ms2')

input_db_text = StringVar()
input_db_text.set('C:/Users/lawashburn/Desktop/ALC50_Mass_Search_Files/duplicate_removed_crustacean_database_validated_formatted20220725.fasta')

prec_err_var = StringVar()
prec_err_var.set('20')

frag_err_var = StringVar()
frag_err_var.set('0.02')

prec_z_var = IntVar()
prec_z_var.set(8)

frag_z_var = IntVar()
frag_z_var.set(4)

amid_check_var = IntVar()
oxM_check_var = IntVar()
PG_E_check_var = IntVar()
PG_Q_check_var = IntVar()
sulf_Y_check_var = IntVar()

max_mods_var = IntVar()
max_mods_var.set(2)

FDR_thresh_var = StringVar()
FDR_thresh_var.set('1.0')

bin_size_var = StringVar()
bin_size_var.set('0.02')

bin_steps_var = IntVar()
bin_steps_var.set(7)

seq_cov_thresh_var = StringVar()
seq_cov_thresh_var.set('0.0')

out_dir_var = StringVar()
out_dir_var.set('C:/Users/lawashburn/Documents/DBpep_v2/results_log/20221201/v30')
###

###

#Add modules window
# modules_click_frame= Canvas(root, width= 460, height= 200)
# modules_click_frame.pack(side=TOP)

modules_title_frame= Canvas(root, width= 450, height= 55)
modules_title_frame.pack(side=TOP)
modules_title_frame.create_text(260, 25, text="Database Search Settings", fill="#2F4FAA", font=(Font_tuple,14),justify=CENTER)
modules_title_frame.pack()

sample_name_settings_frame = Canvas(root, width= 450, height= 55)
sample_name_settings_frame.pack(side=TOP)
sample_name_entry_text_frame = Canvas(sample_name_settings_frame, width= 103, height= 50)
sample_name_entry_text_frame.pack(side=LEFT)
sample_name_entry_text_frame.create_text(53, 25, text="Sample Name", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
sample_name_entry_text_frame.pack()
sample_name_entry_choice_frame = Canvas(sample_name_settings_frame, width= 300, height= 50)
sample_name_entry_choice_frame.pack(side=RIGHT)
sample_name_browse_frame = Canvas(sample_name_entry_choice_frame, width= 300, height= 50)
sample_name_browse_frame.pack(side=TOP)
sample_name_browse_entry = Entry(sample_name_browse_frame, textvariable = sample_name_var, width = 50)
sample_name_browse_entry.pack(side=LEFT)

rawconverter_file_entry_frame = Canvas(root, width= 450, height= 55)
rawconverter_file_entry_frame.pack(side=TOP)
rawconverter_file_entry_text_frame = Canvas(rawconverter_file_entry_frame, width= 103, height= 50)
rawconverter_file_entry_text_frame.pack(side=LEFT)
rawconverter_file_entry_text_frame.create_text(53, 25, text="RawConverter\nPre-processed\nSpectra (.MS2)", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
rawconverter_file_entry_text_frame.pack()
rawconverter_file_entry_choice_frame = Canvas(rawconverter_file_entry_frame, width= 300, height= 50)
rawconverter_file_entry_choice_frame.pack(side=RIGHT)
rawconverter_file_browse_frame = Canvas(rawconverter_file_entry_choice_frame, width= 300, height= 50)
rawconverter_file_browse_frame.pack(side=TOP)
rawconverter_file_choice_browse_entry = Entry(rawconverter_file_browse_frame, textvariable = input_rawconverter_file_text, width = 41)
rawconverter_file_choice_browse_entry.pack(side=LEFT)
rawconverter_file_choice_browse = Button(rawconverter_file_browse_frame, text = "Browse", command = set_path_rawconverter_file_field)
rawconverter_file_choice_browse.pack(side=RIGHT)

db_entry_frame = Canvas(root, width= 450, height= 55)
db_entry_frame.pack(side=TOP)
db_entry_text_frame = Canvas(db_entry_frame, width= 103, height= 50)
db_entry_text_frame.pack(side=LEFT)
db_entry_text_frame.create_text(53, 25, text="Database (.fasta)", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
db_entry_text_frame.pack()
db_entry_choice_frame = Canvas(db_entry_frame, width= 300, height= 50)
db_entry_choice_frame.pack(side=RIGHT)
db_choice_browse_entry = Entry(db_entry_choice_frame, textvariable = input_db_text, width = 41)
db_choice_browse_entry.pack(side=LEFT)
db_choice_browse = Button(db_entry_choice_frame, text = "Browse", command = set_path_database_field)
db_choice_browse.pack(side=RIGHT)

precursor_err_frag_err_settings_frame = Canvas(root, width= 450, height= 55)
precursor_err_frag_err_settings_frame.pack(side=TOP)
prec_err_entry_text_frame = Canvas(precursor_err_frag_err_settings_frame, width= 103, height= 50)
prec_err_entry_text_frame.pack(side=LEFT)
prec_err_entry_text_frame.create_text(53, 25, text="Precursor error \ntolerance (ppm)", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
prec_err_entry_text_frame.pack()
prec_err_entry_choice_frame = Canvas(precursor_err_frag_err_settings_frame, width= 300, height= 50)
prec_err_entry_choice_frame.pack(side=LEFT)
prec_err_browse_frame = Canvas(prec_err_entry_choice_frame, width= 300, height= 50)
prec_err_browse_frame.pack(side=TOP)
prec_err_browse_entry = Entry(prec_err_browse_frame, textvariable = prec_err_var, width = 15)
prec_err_browse_entry.pack(side=LEFT)
frag_err_entry_text_frame = Canvas(precursor_err_frag_err_settings_frame, width= 103, height= 50)
frag_err_entry_text_frame.pack(side=LEFT)
frag_err_entry_text_frame.create_text(53, 25, text="Fragment error \ntolerance (Da)", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
frag_err_entry_text_frame.pack()
frag_err_entry_choice_frame = Canvas(precursor_err_frag_err_settings_frame, width= 300, height= 50)
frag_err_entry_choice_frame.pack(side=LEFT)
frag_err_browse_frame = Canvas(frag_err_entry_choice_frame, width= 300, height= 50)
frag_err_browse_frame.pack(side=TOP)
frag_err_browse_entry = Entry(frag_err_browse_frame, textvariable = frag_err_var, width = 15)
frag_err_browse_entry.pack(side=LEFT)

prec_z_frag_z_settings_frame = Canvas(root, width= 450, height= 55)
prec_z_frag_z_settings_frame.pack(side=TOP)
prec_z_entry_text_frame = Canvas(prec_z_frag_z_settings_frame, width= 103, height= 50)
prec_z_entry_text_frame.pack(side=LEFT)
prec_z_entry_text_frame.create_text(53, 25, text="Max precursor\ncharge", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
prec_z_entry_text_frame.pack()
prec_z_entry_choice_frame = Canvas(prec_z_frag_z_settings_frame, width= 300, height= 50)
prec_z_entry_choice_frame.pack(side=LEFT)
prec_z_browse_frame = Canvas(prec_z_entry_choice_frame, width= 300, height= 50)
prec_z_browse_frame.pack(side=TOP)
prec_z_browse_entry = Entry(prec_z_browse_frame, textvariable = prec_z_var, width = 15)
prec_z_browse_entry.pack(side=LEFT)
frag_z_entry_text_frame = Canvas(prec_z_frag_z_settings_frame, width= 103, height= 50)
frag_z_entry_text_frame.pack(side=LEFT)
frag_z_entry_text_frame.create_text(53, 25, text="Max fragment\ncharge", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
frag_z_entry_text_frame.pack()
frag_z_entry_choice_frame = Canvas(prec_z_frag_z_settings_frame, width= 300, height= 50)
frag_z_entry_choice_frame.pack(side=LEFT)
frag_z_browse_frame = Canvas(frag_z_entry_choice_frame, width= 300, height= 50)
frag_z_browse_frame.pack(side=TOP)
frag_z_browse_entry = Entry(frag_z_browse_frame, textvariable = frag_z_var, width = 15)
frag_z_browse_entry.pack(side=LEFT)

ptm_choice_frame = Canvas(root,width=200,height=25)
ptm_choice_frame.pack(side=TOP)
ptm_choice_entry_text_frame = Canvas(ptm_choice_frame, width= 200, height= 30)
ptm_choice_entry_text_frame.pack(side=LEFT)
ptm_choice_entry_text_frame.create_text(100, 25, text="Variable PTM Selection", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
ptm_choice_entry_text_frame.pack(side=TOP)
ptm_available_frame = Canvas(ptm_choice_frame, width= 103, height= 50)
ptm_available_frame.pack(side=RIGHT)
ptm_available_frame1 = Canvas(ptm_available_frame, width= 103, height= 50)
ptm_available_frame1.pack(side=TOP)
ptm_available_frame2 = Canvas(ptm_available_frame, width= 103, height= 50)
ptm_available_frame2.pack(side=TOP)
amid_check = Checkbutton(ptm_available_frame1, text="C-terminal Amidation", variable=amid_check_var,font=(Font_tuple,8),onvalue=1,offvalue=0)
amid_check.pack(side=LEFT)
OxM_check = Checkbutton(ptm_available_frame1, text="Oxidation on M", variable=oxM_check_var,font=(Font_tuple,8),onvalue=1,offvalue=0)
OxM_check.pack(side=LEFT)
PG_E_check = Checkbutton(ptm_available_frame1, text="Pyro-glu on E", variable=PG_E_check_var,font=(Font_tuple,8),onvalue=1,offvalue=0)
PG_E_check.pack(side=LEFT)
PG_Q_check = Checkbutton(ptm_available_frame2, text="Pyro-glu on Q", variable=PG_Q_check_var,font=(Font_tuple,8),onvalue=1,offvalue=0)
PG_Q_check.pack(side=LEFT)
sulf_Y_check = Checkbutton(ptm_available_frame2, text="Sulfation on Y", variable=sulf_Y_check_var,font=(Font_tuple,8),onvalue=1,offvalue=0)
sulf_Y_check.pack(side=LEFT)

max_mods_FDR_thresh_settings_frame = Canvas(root, width= 450, height= 55)
max_mods_FDR_thresh_settings_frame.pack(side=TOP)
max_mods_entry_text_frame = Canvas(max_mods_FDR_thresh_settings_frame, width= 103, height= 50)
max_mods_entry_text_frame.pack(side=LEFT)
max_mods_entry_text_frame.create_text(53, 25, text="Max # PTMs per\nneuropeptide", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
max_mods_entry_text_frame.pack()
max_mods_entry_choice_frame = Canvas(max_mods_FDR_thresh_settings_frame, width= 300, height= 50)
max_mods_entry_choice_frame.pack(side=LEFT)
max_mods_browse_frame = Canvas(max_mods_entry_choice_frame, width= 300, height= 50)
max_mods_browse_frame.pack(side=TOP)
max_mods_browse_entry = Entry(max_mods_browse_frame, textvariable = max_mods_var, width = 15)
max_mods_browse_entry.pack(side=LEFT)
FDR_thresh_entry_text_frame = Canvas(max_mods_FDR_thresh_settings_frame, width= 103, height= 50)
FDR_thresh_entry_text_frame.pack(side=LEFT)
FDR_thresh_entry_text_frame.create_text(53, 25, text="FDR threshold (%)", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
FDR_thresh_entry_text_frame.pack()
FDR_thresh_entry_choice_frame = Canvas(max_mods_FDR_thresh_settings_frame, width= 300, height= 50)
FDR_thresh_entry_choice_frame.pack(side=LEFT)
FDR_thresh_browse_frame = Canvas(FDR_thresh_entry_choice_frame, width= 300, height= 50)
FDR_thresh_browse_frame.pack(side=TOP)
FDR_thresh_browse_entry = Entry(FDR_thresh_browse_frame, textvariable = FDR_thresh_var, width = 15)
FDR_thresh_browse_entry.pack(side=LEFT)

bin_size_bin_steps_settings_frame = Canvas(root, width= 450, height= 55)
bin_size_bin_steps_settings_frame.pack(side=TOP)
bin_size_entry_text_frame = Canvas(bin_size_bin_steps_settings_frame, width= 103, height= 50)
bin_size_entry_text_frame.pack(side=LEFT)
bin_size_entry_text_frame.create_text(53, 25, text="Bin size", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
bin_size_entry_text_frame.pack()
bin_size_entry_choice_frame = Canvas(bin_size_bin_steps_settings_frame, width= 300, height= 50)
bin_size_entry_choice_frame.pack(side=LEFT)
bin_size_browse_frame = Canvas(bin_size_entry_choice_frame, width= 300, height= 50)
bin_size_browse_frame.pack(side=TOP)
bin_size_browse_entry = Entry(bin_size_browse_frame, textvariable = bin_size_var, width = 15)
bin_size_browse_entry.pack(side=LEFT)
bin_steps_entry_text_frame = Canvas(bin_size_bin_steps_settings_frame, width= 103, height= 50)
bin_steps_entry_text_frame.pack(side=LEFT)
bin_steps_entry_text_frame.create_text(53, 25, text="Number of steps", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
bin_steps_entry_text_frame.pack()
bin_steps_entry_choice_frame = Canvas(bin_size_bin_steps_settings_frame, width= 300, height= 50)
bin_steps_entry_choice_frame.pack(side=LEFT)
bin_steps_browse_frame = Canvas(bin_steps_entry_choice_frame, width= 300, height= 50)
bin_steps_browse_frame.pack(side=TOP)
bin_steps_browse_entry = Entry(bin_steps_browse_frame, textvariable = bin_steps_var, width = 15)
bin_steps_browse_entry.pack(side=LEFT)

seq_cov_thresh_settings_frame = Canvas(root, width= 450, height= 55)
seq_cov_thresh_settings_frame.pack(side=TOP)
seq_cov_thresh_entry_text_frame = Canvas(seq_cov_thresh_settings_frame, width= 103, height= 50)
seq_cov_thresh_entry_text_frame.pack(side=LEFT)
seq_cov_thresh_entry_text_frame.create_text(53, 25, text="Mimimum sequence\ncoverage threshold\n(optional)", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
seq_cov_thresh_entry_text_frame.pack()
seq_cov_thresh_entry_choice_frame = Canvas(seq_cov_thresh_settings_frame, width= 300, height= 50)
seq_cov_thresh_entry_choice_frame.pack(side=RIGHT)
seq_cov_thresh_browse_frame = Canvas(seq_cov_thresh_entry_choice_frame, width= 300, height= 50)
seq_cov_thresh_browse_frame.pack(side=TOP)
seq_cov_thresh_browse_entry = Entry(seq_cov_thresh_browse_frame, textvariable = seq_cov_thresh_var, width = 41)
seq_cov_thresh_browse_entry.pack(side=LEFT)

out_dir_entry_frame = Canvas(root, width= 450, height= 55)
out_dir_entry_frame.pack(side=TOP)
out_dir_entry_text_frame = Canvas(out_dir_entry_frame, width= 103, height= 50)
out_dir_entry_text_frame.pack(side=LEFT)
out_dir_entry_text_frame.create_text(53, 25, text="Output directory", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
out_dir_entry_text_frame.pack()
out_dir_entry_choice_frame = Canvas(out_dir_entry_frame, width= 300, height= 50)
out_dir_entry_choice_frame.pack(side=RIGHT)
out_dir_choice_browse_entry = Entry(out_dir_entry_choice_frame, textvariable = out_dir_var, width = 45)
out_dir_choice_browse_entry.pack(side=LEFT)
out_dir_choice_browse = Button(out_dir_entry_choice_frame, text = "Browse", command = set_path_out_dir_field)
out_dir_choice_browse.pack(side=RIGHT)

run_click_frame= Canvas(root, width= 450, height= 525)
run_click_frame.pack(side=TOP)

# save_button = Button(run_click_frame, width=450,height=1,text='Run Analysis',fg='white',relief='flat',borderwidth=5, 
#                     bg='#2F4FAA',font=(Font_tuple,15),command = threading.Thread(target=combine_funcs(param_log_export,Run_Search)).start)

save_button = Button(run_click_frame, width=450,height=1,text='Run Analysis',fg='white',relief='flat',borderwidth=5, 
                    bg='#2F4FAA',font=(Font_tuple,15),command = lambda: threading.Thread(target=combine_funcs(param_log_export,Run_Search)).start())

save_button.pack(side=BOTTOM)


root.mainloop()
