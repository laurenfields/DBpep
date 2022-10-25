# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 10:16:13 2022

@author: lawashburn
"""

from tkinter import *
from PIL import ImageTk, Image
from tkinter.filedialog import askopenfilename
import sys
from itertools import islice
from subprocess import Popen, PIPE
from textwrap import dedent
from tkinter import messagebox
from tkinter import filedialog
import os
import pickle
from tkinter import ttk
import threading
import time

def hide_all_frames():
    logo_frame.pack_forget()
    about_click_frame.pack_forget()
    #modules_click_frame.pack_forget()
    files_click_frame.pack_forget()
    run_click_frame.pack_forget()
    help_click_frame.pack_forget()
def about_frame():
    hide_all_frames()
    about_click_frame.pack(side = RIGHT)
def add_modules():
    hide_all_frames()
    modules_click_frame.pack(side = RIGHT)  
def select_files():
    hide_all_frames()
    files_click_frame.pack(side = RIGHT)
def run_hypep():
    hide_all_frames()
    run_click_frame.pack(side = RIGHT)
def help_form():
    hide_all_frames()
    help_click_frame.pack(side = RIGHT)

#GUI appearance settings
Font_tuple = ("Corbel Light", 20)
root = Tk()
root.title('HyPep 1.0')
root.iconbitmap(r"hypep_icon.ico")
root.geometry('700x525')
frame= Frame(root, width= 20, height= 700, bg = '#2F4FAA')
frame.pack(side=LEFT,  anchor=N)
logo_frame= Frame(root, width= 700, height= 700, bg = 'purple')
logo_frame.pack(side=RIGHT, anchor=N)
###
im1 = Image.open(r"About.png")
resized_im1 = im1.resize((218,50))
button1 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=about_frame)
ph_im1 = ImageTk.PhotoImage(resized_im1) # <----------
button1.config(image=ph_im1)
button1.pack(fill=X, expand=1, anchor=N)

im2 = Image.open(r"Add Modules2.png")
resized_im2 = im2.resize((218,50))
button2 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=add_modules)
ph_im2 = ImageTk.PhotoImage(resized_im2) # <----------
button2.config(image=ph_im2)
button2.pack(fill=X, expand=1, anchor=N)

im3 = Image.open(r"Select Files.png")
resized_im3 = im3.resize((218,50))
button3 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA' ,command=select_files)
ph_im3 = ImageTk.PhotoImage(resized_im3) # <----------
button3.config(image=ph_im3)
button3.pack(fill=X, expand=1, anchor=N)


button4 = Canvas(frame, width=218,height=50, bg='#2F4FAA',highlightbackground = '#2F4FAA' )
button4.pack(fill=X, expand=1, anchor=N)

im5 = Image.open(r"Run Analysis.png")
resized_im5 = im5.resize((218,50))
button5 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=run_hypep)
ph_im5 = ImageTk.PhotoImage(resized_im5) # <----------
button5.config(image=ph_im5)
button5.pack(fill=X, expand=1, anchor=N)


button6 = Canvas(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',highlightbackground = '#2F4FAA' )
button6.pack(fill=X, expand=1, anchor=N, pady=50)

im7 = Image.open(r"Help.png")
resized_im7 = im7.resize((218,50))
button7 = Button(frame, width=218,height=50,relief='flat',borderwidth=5, bg='#2F4FAA',command=help_form)
ph_im7 = ImageTk.PhotoImage(resized_im7) # <----------
button7.config(image=ph_im7)
button7.pack(fill=X, expand=1, anchor=N)


logo = Image.open(r"hypep_logo.png")
resized_logo = logo.resize((500,500))
logo_pos = Label(logo_frame, width=500,height=500,relief='flat',borderwidth=5)
logo2 = ImageTk.PhotoImage(resized_logo) # <----------
logo_pos.config(image=logo2)
logo_pos.pack(expand= YES,fill= BOTH, anchor=CENTER)

#about window
about_click_frame= Canvas(root, width= 700, height= 500)
about_click_frame.pack()
title_text_box= Canvas(about_click_frame, width= 500, height= 100)
title_text_box.pack(anchor=N)
title_text_box.create_text(250, 50, text="About HyPep", fill="#2F4FAA", font=(Font_tuple,25),justify=CENTER)
title_text_box.pack(anchor=N)
body_text_box= Canvas(about_click_frame, width= 500, height= 500)
body_text_box.pack(anchor=CENTER)
body_text_box.create_text(225, 100, text="HyPep is a hybrid neuropeptide \nidentification software composed of a sequence homology search \nalgorithm and an accurate \nmass matching algorithm", 
                              fill="#2F4FAA", font=Font_tuple, justify=LEFT, width = 400)
body_text_box.pack(expand= YES,fill= BOTH, anchor=CENTER)
#Two main pages: DB search setup and select files, also run page

def param_file_input_save():
    #these "get" values are not for variable storage, just for checking if the input values are valid
    #This is the command for when the user selects "save preferences" on the file selection input page
    space_sample_name = input_sample_file_text.get()

def combine_funcs(*funcs):
    def combined_func(*args, **kwargs):
        for f in funcs:
            f(*args, **kwargs)
    return combined_func


def param_log_export():

    db_path = input_db_text.get()
    db_path_report = 'Database path: ' + db_path

    rawconverter_path = input_rawconverter_file_text.get()
    rawconverter_path_report = 'RawConverter output file path: ' + rawconverter_path

    out_path = input_out_dir_text.get()
    out_path_report = 'Ouput directory path: ' + out_path

    param_file_entries = [db_path_report,
                         rawconverter_path_report,
                          out_path_report]

#select files page, three input areas- .fasta, output directory, and spectra import
path_db = StringVar()
input_db_text = StringVar()
input_rawconverter_file_text =  StringVar()
input_out_dir_text = StringVar()
files_click_frame= Canvas(root, width= 450, height= 525)
files_click_frame.pack(side=TOP)

def set_path_database_field():
    path_db = askopenfilename(filetypes=[("FASTA Files","*.fasta")]) 
    input_db_text.set(path_db)
def get_database_path(): 
    """ Function provides the database full file path."""
    return path_db

def set_path_out_dir_field():
    path_out_dir = filedialog.askdirectory() 
    input_out_dir_text.set(path_out_dir)
def get_out_dir_path(): 
    """ Function provides the database full file path."""
    return path_out_dir

def set_path_rawconverter_file_field():
    path_rawconverter_file = askopenfilename(filetypes=[("MS2 Files","*.MS2")],multiple=False) 
    input_rawconverter_file_text.set(path_rawconverter_file)
def get_rawconverter_file_path(): 
    """ Function provides the database full file path."""
    return path_rawconverter_file
def RunHyPep():
    time.sleep(5)
    os.system('python command_center.py')
def retrievedata():
    ''' get data stored '''
    global list_data
    list_data = []
    try:
      with open("save.txt", "r", encoding="utf-8") as file:
       for f in file:
        listbox.insert(tk.END, f.strip())
        list_data.append(f.strip())
        print(list_data)
    except:
        pass

def reload_data():
    listbox.delete(0, tk.END)
    for d in list_data:
        listbox.insert(0, d)
def add_item(event=1):
    global list_data
    if input_rawconverter_file_text.get() != "":
        listbox.insert(END, input_rawconverter_file_text.get())
        test = input_rawconverter_file_text.get()
        list_data.append(input_rawconverter_file_text.get())
        input_rawconverter_file_text.set("")
def delete():
    global list_data
    listbox.delete(0, END)
    list_data = []
def delete_selected():

    try:
        selected = listbox.get(listbox.curselection())
        listbox.delete(listbox.curselection())
        list_data.pop(list_data.index(selected))
        # reload_data()
        # # listbox.selection_clear(0, END)
        listbox.selection_set(0)
        listbox.activate(0)
        listbox.event_generate("&lt;&lt;ListboxSelect>>")
        print(listbox.curselection())
    except:
        pass
files_title_text_box= Canvas(files_click_frame, width= 450, height= 75)
files_title_text_box.pack(side=TOP)
files_title_text_box.create_text(250, 25, text="Select Files", fill="#2F4FAA", font=(Font_tuple,20),justify=CENTER)
files_title_text_box.pack(side=TOP)

check_postion = Canvas(files_click_frame, width= 450, height= 375)
check_postion.pack(side=TOP)

db_entry_frame = Canvas(check_postion, width= 450, height= 55)
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

rawconverter_file_entry_frame = Canvas(check_postion, width= 450, height= 55)
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

rawconverter_file_options_frame = Canvas(rawconverter_file_entry_choice_frame, width= 300, height= 50)
rawconverter_file_options_frame.pack(side=BOTTOM)
button = Button(rawconverter_file_options_frame, text="Add Item", command=(add_item))
button.pack(side=LEFT)
button_delete_selected = Button(rawconverter_file_options_frame,text="Delete Selected", command=delete_selected)
button_delete_selected.pack(side=RIGHT)
button_delete = Button(rawconverter_file_options_frame,text="Clear all", command=delete)
button_delete.pack(side=RIGHT)


rawconverter_file_list_frame = Canvas(rawconverter_file_entry_choice_frame, width= 300, height= 50)
rawconverter_file_list_frame.pack(side=BOTTOM)

listbox = Listbox(rawconverter_file_list_frame, width = 50)
listbox.pack()
rawconverter_file_entry_frame.bind("&lt;Return>", add_item)
retrievedata()
out_dir_entry_frame = Canvas(check_postion, width= 450, height= 55)
out_dir_entry_frame.pack(side=TOP)
out_dir_entry_text_frame = Canvas(out_dir_entry_frame, width= 103, height= 50)
out_dir_entry_text_frame.pack(side=LEFT)
out_dir_entry_text_frame.create_text(53, 25, text="Output directory", fill="#2F4FAA", font=(Font_tuple,8),justify=CENTER)
out_dir_entry_text_frame.pack()
out_dir_entry_choice_frame = Canvas(out_dir_entry_frame, width= 300, height= 50)
out_dir_entry_choice_frame.pack(side=RIGHT)
out_dir_choice_browse_entry = Entry(out_dir_entry_choice_frame, textvariable = input_out_dir_text, width = 41)
out_dir_choice_browse_entry.pack(side=LEFT)
out_dir_choice_browse = Button(out_dir_entry_choice_frame, text = "Browse", command = set_path_out_dir_field)
out_dir_choice_browse.pack(side=RIGHT)

save_button_frame= Canvas(files_click_frame, width= 450, height= 100)
save_button_frame.pack(side=BOTTOM)
save_button = Button(save_button_frame, width=450,height=5,text='Save preferences',fg='white',relief='flat',borderwidth=5, 
                    bg='#2F4FAA',font=(Font_tuple,15),command=param_file_input_save)
save_button.pack(side=BOTTOM)



#Run page
run_click_frame= Canvas(root, width= 450, height= 525)

pb = ttk.Progressbar(
    run_click_frame,
    orient='horizontal',
    mode='indeterminate',
    length=280
)
pb.pack(side=BOTTOM)
save_button = Button(run_click_frame, width=450,height=1,text='Run Analysis!',fg='white',relief='flat',borderwidth=5, 
                    bg='#2F4FAA',font=(Font_tuple,15),command = threading.Thread(target=combine_funcs(param_log_export,pb.start,RunHyPep)).start)
save_button.pack(side=BOTTOM)

#Help page
help_click_frame= Canvas(root, width= 450, height= 525)
help_click_frame.pack(side=TOP)
help_text_box= Canvas(help_click_frame, width= 500, height= 100)
help_text_box.pack(anchor=N)
help_text_box.create_text(250, 50, text="Help", fill="#2F4FAA", font=(Font_tuple,25),justify=CENTER)
help_text_box.pack(anchor=N)
help_body_text_box= Canvas(help_click_frame, width= 500, height= 500)
help_body_text_box.pack(anchor=CENTER)
help_body_text_box.create_text(225, 100, text="Please direct all HyPep questions to Lauren Fields, lawashburn@wisc.edu", 
                              fill="#2F4FAA", font=Font_tuple, justify=LEFT, width = 400)
help_body_text_box.pack(expand= YES,fill= BOTH, anchor=CENTER)
root.mainloop()
