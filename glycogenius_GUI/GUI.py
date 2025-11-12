# GlycoGenius: Glycomics Data Analysis Tool
# Copyright (C) 2023 by Hector Franco Loponte
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or 
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. It is accessible within the program files
# or by typing 'license' after running it stand-alone in the terminal
# by typing 'glycogenius'. If not, see <https://www.gnu.org/licenses/>.

global gg_version, GUI_version
gg_version = '1.3.1'
GUI_version = '1.1.2'

from PIL import Image, ImageTk
from tkinter import messagebox
import tkinter as tk
import pathlib
import psutil
import os
import tempfile
import platform

# Get working directory
global current_dir, ico_image
current_dir = pathlib.Path(__file__).parent.resolve()
        
# Load icon to be used in the program
ico_path = os.path.join(current_dir, "Assets/gg_icon.ico")
ico_image = Image.open(ico_path)
    
def start_splash():
    global splash_screen

    exec_check_folder = os.path.join(tempfile.gettempdir())
    
    this_process_id = os.getpid()
    this_process = psutil.Process(this_process_id)
    this_process_ppid = this_process.ppid()
    if f"gg_{this_process_ppid}.txt" not in os.listdir(exec_check_folder):
        with open(os.path.join(exec_check_folder, f"gg_{this_process_id}.txt"), 'w') as f:
            f.write("Glycogenius has run")
            f.close()
    
        # Check if GG folder is writeable
        global gg_draw_glycans_path
        try:
            gg_draw_glycans_path = os.path.join(current_dir, "Assets/glycans")
            os.makedirs(gg_draw_glycans_path, exist_ok=True)
            with open(os.path.join(gg_draw_glycans_path, f"test.txt"), 'w') as f:
                f.write("Glycogenius can access this folder")
                f.close()
            os.remove(os.path.join(gg_draw_glycans_path, f"test.txt"))
        except Exception:
            messagebox.showwarning("Warning", "No permission to access GlycoGenius folder to save glycans figures in. Using temporary folder instead. This may lead to the need of building the figures from scratch again in a future activation of GG Draw. GlycoGenius will still work as intended.")
            gg_draw_glycans_path = os.path.join(tempfile.gettempdir(), "glycans_gg")
            os.makedirs(gg_draw_glycans_path, exist_ok=True)
            
        splash_screen = tk.Tk()
        splash_screen.withdraw()
        icon = ImageTk.PhotoImage(ico_image)
        splash_screen.iconphoto(False, icon)
        splash_screen.overrideredirect(True)
        splash_screen.resizable(False, False)
        splash_screen.attributes("-topmost", True)
    
        # Get the platform at which GG is running on to set whether splash logo can have transparent background
        curr_os = platform.system()
        if curr_os == "Windows":
            splash_screen.attributes("-transparentcolor", "white")
            
        splash_screen.grab_set()
        
        splash_image = Image.open(os.path.join(current_dir, "Assets/splash.png"))
        splash_size = splash_image.size
        splash_image = splash_image.resize((int(splash_size[0]/1.5), int(splash_size[1]/1.5)))
        tk_splash = ImageTk.PhotoImage(splash_image)
        
        splash_screen_label = tk.Label(splash_screen, bg="white", image=tk_splash)
        splash_screen_label.pack(pady=0, padx=0)
        
        splash_screen.update_idletasks()
        splash_screen.deiconify()
        window_width = splash_screen.winfo_width()
        window_height = splash_screen.winfo_height()
        screen_width = splash_screen.winfo_screenwidth()
        screen_height = splash_screen.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        splash_screen.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")

        # Update the GUI manually
        splash_screen.update()

# This has to come before importing everything, as these imports can take some time!
start_splash()

from glycogenius.Modules.core import main as run_glycogenius
from glycogenius.Modules import Execution_Functions, General_Functions, Config_Handler
from tkinter import Menu, ttk, filedialog, colorchooser
from tkinter.scrolledtext import ScrolledText
from ttkwidgets import CheckboxTreeview
from pyteomics import mass, mzxml, mzml
from itertools import product
from scipy.stats import gaussian_kde
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.colors import LogNorm, PowerNorm
from matplotlib.ticker import ScalarFormatter
from mpl_scatter_density import ScatterDensityArtist
from collections import defaultdict
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import threading
import matplotlib
import shutil
import importlib
import sys
import ctypes
import traceback
import multiprocessing
import concurrent.futures
import copy
import datetime
import zipfile
import dill
import random
import math

# Windows only
if platform.system() == 'Windows':
    import winreg as reg

# For File Editor
from reportlab.lib.pagesizes import letter, A4
from reportlab.pdfgen import canvas as pdf_canvas
from reportlab.lib.utils import ImageReader
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle

# For GG Draw
from PIL import ImageDraw

# Fix multiprocessing on linux and presumably mac
if platform.system() != 'Windows':
    multiprocessing.set_start_method('spawn', force=True)

# All the settings necessary to make a Glycogenius run
global min_max_monos, min_max_hex, min_max_hn, min_max_hexnac, min_max_fuc, min_max_sia, min_max_ac, min_max_ac, min_max_gc, min_max_ua, forced, min_max_proton, max_charges, reducing_end_tag, internal_standard, permethylated, lactonized_ethyl_esterified, reduced, fast_iso, high_res, number_cores, multithreaded_analysis, exp_lib_name, min_samples, lyase_digested, custom_monosaccharides, ion_mode_selected

custom_glycans_list = [False, '']
min_max_monos = [5, 22]
min_max_hex = [3, 10]
min_max_hn = [0, 0]
min_max_hexnac = [2, 8]
min_max_xyl = [0, 0]
min_max_sia = [0, 4]
min_max_fuc = [0, 2]
min_max_ac = [0, 4]
min_max_gc = [0, 0]
min_max_ua = [0, 0]
custom_monosaccharides = []
forced = 'n_glycans'
min_max_proton = [1, 3]
custom_adducts = []
max_charges = 3
ion_mode_selected = False
reducing_end_tag = 0.0
permethylated = False
reduced = False
lactonized_ethyl_esterified = False
min_max_sulfation = [0, 0]
min_max_phosphorylation = [0, 0]
lyase_digested = False
fast_iso = True
high_res = False
internal_standard = '0.0'
imp_exp_library = [False, False]
exp_lib_name = ''
library_path = ''
only_gen_lib = False

multithreaded_analysis = True
number_cores = 'all'
analyze_ms2 = [False, True, False]
reporter_ions = []
tolerance = ['ppm', 20]
ret_time_interval = [0.0, 999.0, 0.2]
rt_tolerance_frag = 0.2
min_isotopologue_peaks = 2
min_ppp = [False, 10]
close_peaks = [False, 5]
align_chromatograms = False
percentage_auc = 0.01
min_samples = 0
max_ppm = [-10, 10]
iso_fit_score = 0.9
curve_fit_score = 0.9
s_to_n = 3
fill_gaps = (False, 50, 0.2, False)
custom_noise = [False, []]
samples_path = ''
save_path = ''
plot_metaboanalyst = [False, '']
compositions = True
iso_fittings = False
reanalysis = False
reanalysis_path = ''
output_plot_data = False

verbose = False

global samples_list, samples_names
samples_list = []
samples_names = []

samples_dropdown_options = []
reducing_end_boolean = False

global gg_file_name
library_name = ""
gg_file_name = ""

global former_reanalysis_list
former_reanalysis_list = []

reducing_end_tags = {
    '2-AB': 'C7H8N2', 
    '2-AA': 'C7H7N1O1',
    'ProA': 'C13H21N3', 
    'GirP': 'C7H7N3',
    'APTS': 'C16H11N1O8S3'
    }

forced_classes = {
    'None': 'none', 
    'N-Glycans': 'n_glycans', 
    'O-Glycans': 'o_glycans', 
    'GAGs': 'gags'
    }

global which_progress_bar
which_progress_bar = 'ms1'

global qc_dist_opened
qc_dist_opened = False

global compare_samples_opened, former_alignments
compare_samples_opened = False
former_alignments = []

global max_spectrum_opened, maximum_spectra, glycans_list_quickcheck
max_spectrum_opened = False
maximum_spectra = {}
glycans_list_quickcheck = {}
glycans_list_quickcheck_save = {}

global quick_trace_opened, quick_traces_all, quick_traces_list_save
quick_trace_opened = False
quick_traces_all = {}
quick_traces_list_save = []

suppressed_prints = ["remaining parameters.", "set 'only_gen_lib' to False and input", "If you wish to analyze files,", "Press Enter to exit.", "Close the window or press CTRL+C to exit."]

bpc = {}

global sample_groups
sample_groups = {}

global last_xlims_chrom, last_ylims_chrom
last_xlims_chrom = None
last_ylims_chrom = None

global ms1_bound, ms1_binds, ms2_bound, ms2_binds, chromatogram_bound, chromatogram_binds
ms1_bound = False
ms1_binds = []
ms2_bound = False
ms2_binds = []
chromatogram_bound = False
chromatogram_binds = []

global clicked_spectra, current_line_ruler, distance_ruler_text, ruler_lines_list, selected_ruler
clicked_spectra = None
current_line_ruler = None
distance_ruler_text = None
selected_ruler = None
ruler_lines_list = []

global gg_draw_on, gg_draw_list, gg_draw_selected
gg_draw_on = False
gg_draw_list = []
gg_draw_selected = None

global former_selected_sample
former_selected_sample = None

global zoom_lock_chromatogram, zoom_lock_spectrum
zoom_lock_chromatogram = {'on': False, 'x_axis': (0, 0), 'y_axis': (0, 0)}
zoom_lock_spectrum = {'on': False, 'x_axis': (0, 0), 'y_axis': (0, 0)}

global key_states, selected_chromatograms, last_plotted
last_plotted = None
key_states = {}
selected_chromatograms = []

global original_lines
original_lines = {}

global gg_draw_zoom_scale_value_float
gg_draw_zoom_scale_value_float = 0.2

colors = [
    "#FF0000",  # Red
    "#0000FF",  # Blue
    "#00FF00",  # Green
    "#FFFF00",  # Yellow
    "#800080",  # Purple
    "#00FFFF",  # Cyan
    "#FFA500",  # Orange
    "#32CD32",  # Lime Green
    "#FFD700",  # Gold
    "#008080",  # Teal
    "#000080",  # Navy
    "#808000",  # Olive
    "#FF8C00",  # Dark Orange
    "#006400",  # Dark Green
    "#D2691E"   # Chocolate
]
color_number = 0

list_font_size = 10
list_font_size_smaller = 8
button_font_size = 10
big_button_size = (int(button_font_size), int(button_font_size*2))

global standard_library_parameters
standard_library_parameters = {
                               'N-glycans': {'monos': [3, 22], 'h': [1, 10], 'hn': [0, 0], 'n': [2, 8], 'x': [0, 0], 'f': [0, 2], 'sia': [0, 4], 'ua': [0, 0], 's': [0, 4], 'g': [0, 0], 'class': 'n_glycans', 'lyase': False},
                               'O-glycans': {'monos': [1, 10], 'h': [0, 8], 'hn': [0, 0], 'n': [1, 8], 'x': [0, 0], 'f': [0, 2], 'sia': [0, 4], 'ua': [0, 0], 's': [0, 4], 'g': [0, 0], 'class': 'o_glycans', 'lyase': False},
                               'GAGs': {'monos': [1, 6], 'h': [0, 0], 'hn': [0, 3], 'n': [0, 3], 'x': [0, 0], 'f': [0, 0], 'sia': [0, 0], 'ua': [0, 3], 's': [0, 0], 'g': [0, 0], 'class': 'gags', 'lyase': True}
                               }

# Classes
class gg_archive:
    '''A class to manage access to .gg files.
    
    Attributes
    ----------
    path : string
        The path to the .gg file.
        
    temp_path : string
        The path to the temporary folder to store data for access.
        
    Methods
    -------
    check_gg_archive
        Checks the contents of the gg_archive to see if it is actually a .gg file.
        
    extract
        Extracts the full content of the .gg file to the temporary folder.
        
    extract_file(file, origin)
        Extracts a single file from the .gg file to the temporary folder.
        
    list_samples
        Lists the samples in the .gg file.
        
    list_chromatograms(file_number)
        Returns a list of all the names of the chromatograms available for a given file number.
        
    get_rt_array(file_number)
        Returns the RT array for this sample.
        
    get_chromatogram(file_number, chromatogram_name, chromatogram_type)
        Extracts the target chromatogram from the target file and the target type to the temporary folder and returns its data.
        
    get_results_table
        Returns the results table of this .gg file analysis.
        
    get_version
        Retrieves the GlycoGenius version at which the sample was analyzed in.
        
    get_metadata
        Returns the metadata from the .gg file.
        
    get_isotopic_fittings
        Returns the isotopic fittings from the .gg file.
        
    get_curve_fittings
        Returns the curve fittings from the .gg file.
        
    close
        Unloads the data from the .gg file and cleans up the temporary folder.
    '''
    def __init__(self, path, temp_path = ''):
        '''Initializes the object, which contains the path to the .gg file and the path to the temporary folder. If not path to the temporary folder is provided, creates it by itself'''
        self.path = path
        
        if temp_path != '':
            self.temp_path = temp_path
        else:
            self.begin_time = datetime.datetime.now()
            self.begin_time_string = str(self.begin_time)[2:4]+str(self.begin_time)[5:7]+str(self.begin_time)[8:10]+"_"+str(self.begin_time)[11:13]+str(self.begin_time)[14:16]+str(self.begin_time)[17:19]
            while "gg_"+self.begin_time_string in os.listdir(tempfile.gettempdir()):
                self.begin_time_string = f"{self.begin_time_string[:-2]}{int(self.begin_time_string[-2:])+1}"
                if int(self.begin_time_string[-2:]) > 59:
                    self.begin_time_string = f"{self.begin_time_string[:-4]}{int(self.begin_time_string[-4:-2])+1}00"
            self.temp_path = os.path.join(tempfile.gettempdir(), "gg_"+self.begin_time_string)
            os.makedirs(self.temp_path, exist_ok=True)
            
        self.check_gg_archive()
            
    def check_gg_archive(self):
        '''Checks the contents of the gg_archive to see if it is actually a .gg file'''
        try:
            with zipfile.ZipFile(self.path, 'r') as zip_ref:
                self._file_list = zip_ref.namelist()
                if 'results' not in self._file_list or 'eics_list' not in self._file_list or 'isotopic_fittings' not in self._file_list or 'curve_fittings' not in self._file_list or 'metadata' not in self._file_list:
                    self.close_gg()
                zip_ref.close()
        except:
            self.close_gg()
        
    def extract(self):
        '''Extracts all the content from the .gg using zipfile.'''
        with zipfile.ZipFile(self.path, 'r') as zip_ref:
            zip_ref.extractall(self.temp_path)
            zip_ref.close()
    
    def extract_file(self, file, origin = ''):
        '''Extracts a single file from the .gg file or another zipfile.'''
        if origin == '':
            origin = self.path
        with zipfile.ZipFile(origin, 'r') as zip_ref:
            zip_ref.extract(file, self.temp_path)
            zip_ref.close()
        
    def list_chromatograms(self, file_number):
        '''Extracts the eics list from the .gg file and returns the chromatograms names.'''
        try:
            with open(os.path.join(self.temp_path, f'eics_list'), 'rb') as f:
                self._chromatograms_list = dill.load(f)
                f.close()
        except:
            self.extract_file(f'eics_list')
            with open(os.path.join(self.temp_path, f'eics_list'), 'rb') as f:
                self._chromatograms_list = dill.load(f)
                f.close()
            
        return self._chromatograms_list[int(file_number)][1:]
        
    def list_samples(self):
        '''Extracts the results table, acquires the samples list from it and returns it.'''
        self._results = self.get_results_table()
        
        self._samples_metadata = self._results[1]
        
        self._samples_list = {}
        
        for index, i in enumerate(self._samples_metadata['File_Name']):
            self._samples_list[index] = i
            
        return self._samples_list
        
    def get_rt_array(self, file_number):
        '''Extracts the rt array file and returns it.'''
        try:
            with open(os.path.join(self.temp_path, f'{file_number}_RTs'), 'rb') as f:
                _rt_array = dill.load(f)
                f.close()
        except:
            try:
                self.extract_file(f'{file_number}_RTs', os.path.join(self.temp_path, f'{file_number}_eics'))
            except:
                self.extract_file(f'{file_number}_eics')
                self.extract_file(f'{file_number}_RTs', os.path.join(self.temp_path, f'{file_number}_eics'))
                
            with open(os.path.join(self.temp_path, f'{file_number}_RTs'), 'rb') as f:
                _rt_array = dill.load(f)
                f.close()
            
        return _rt_array
            
    def get_chromatogram(self, file_number, chromatogram_name, chromatogram_type):
        '''Extracts the corresponding chromatogram file from the .gg file and returns it.'''
        try:
            with open(os.path.join(self.temp_path, f'{file_number}_{chromatogram_type}_{chromatogram_name}'), 'rb') as f:
                _chromatogram = dill.load(f)
                f.close()
        except:
            try:
                self.extract_file(f'{file_number}_{chromatogram_type}_{chromatogram_name}', os.path.join(self.temp_path, f'{file_number}_eics'))
            except:
                self.extract_file(f'{file_number}_eics')
                self.extract_file(f'{file_number}_{chromatogram_type}_{chromatogram_name}', os.path.join(self.temp_path, f'{file_number}_eics'))
            with open(os.path.join(self.temp_path, f'{file_number}_{chromatogram_type}_{chromatogram_name}'), 'rb') as f:
                _chromatogram = dill.load(f)
                f.close()

        return _chromatogram
            
    def get_fragments_library_ms2_scores(self):
        ''''''
        try:
            with open(os.path.join(self.temp_path, 'fragments_library'), 'rb') as f:
                _fragments_library = dill.load(f)
                f.close()
                
            with open(os.path.join(self.temp_path, 'spectra_score'), 'rb') as f:
                _spectra_score = dill.load(f)
                f.close()
        except:
            try:
                self.extract_file('fragments_library')
                self.extract_file('spectra_score')
                
                with open(os.path.join(self.temp_path, 'fragments_library'), 'rb') as f:
                    _fragments_library = dill.load(f)
                    f.close()
                    
                with open(os.path.join(self.temp_path, 'spectra_score'), 'rb') as f:
                    _spectra_score = dill.load(f)
                    f.close()
            except:
                _fragments_library = {}
                _spectra_score = {}
                
        return _fragments_library, _spectra_score
        
    def get_results_table(self):
        '''Extracts the results table file from the .gg file and returns it.'''
        try:
            with open(os.path.join(self.temp_path, 'results'), 'rb') as f:
                self._results = dill.load(f)
                f.close()
        except:
            self.extract_file('results')
            with open(os.path.join(self.temp_path, 'results'), 'rb') as f:
                self._results = dill.load(f)
                f.close()
            
        return self._results
        
    def get_version(self):
        '''Retrieves the GlycoGenius version at which the sample was analyzed in.'''
        self._results = self.get_results_table()
        if len(self._results) == 4:
            self._version = self._results[3]
        else:
            self._version = self._results[2]
            
        return self._version
        
    def get_metadata(self):
        '''Extracts the metadata file from the .gg file and returns it.'''
        try:
            with open(os.path.join(self.temp_path, 'metadata'), 'rb') as f:
                self._metadata = dill.load(f)
                f.close()
        except:
            self.extract_file('metadata')
            with open(os.path.join(self.temp_path, 'metadata'), 'rb') as f:
                self._metadata = dill.load(f)
                f.close()
                
        return self._metadata
        
    def get_isotopic_fittings(self):
        '''Extracts the isotopic fittings file from the .gg file and returns it.'''
        try:
            with open(os.path.join(self.temp_path, 'isotopic_fittings'), 'rb') as f:
                self._isotopic_fittings = dill.load(f)
                f.close()
        except:
            self.extract_file('isotopic_fittings')
            with open(os.path.join(self.temp_path, 'isotopic_fittings'), 'rb') as f:
                self._isotopic_fittings = dill.load(f)
                f.close()
                
        return self._isotopic_fittings
        
    def get_curve_fittings(self):
        '''Extracts the curve fittings file from the .gg file and returns it.'''
        try:
            with open(os.path.join(self.temp_path, 'curve_fittings'), 'rb') as f:
                self._curve_fittings = dill.load(f)
                f.close()
        except:
            self.extract_file('curve_fittings')
            with open(os.path.join(self.temp_path, 'curve_fittings'), 'rb') as f:
                self._curve_fittings = dill.load(f)
                f.close()
            
        return self._curve_fittings
        
    def close_gg(self):
        '''Unloads the variables within this object and cleans up the temporary folder'''
        shutil.rmtree(self.temp_path)
        
        self.__dict__.clear()
        
        self._is_closed = True

    def __repr__(self):
        '''Provide a string representation of the object, displaying a message if it's been closed.'''
        if getattr(self, '_is_closed', False):
            return f"File is unloaded."
        else:
            return f"Temporary folder: {self.temp_path}.\n.gg file path: {self.path}"
            
class ImageGallery(tk.Frame):
    def __init__(self, parent, image_paths=None, columns=3, cell_size=(120, 120)):
        super().__init__(parent)

        self.columns = columns  # Number of columns in the grid layout
        self.cell_width, self.cell_height = cell_size  # Size of each grid cell

        # Create a canvas to hold the images
        self.canvas = tk.Canvas(self)
        self.scroll_y = ttk.Scrollbar(self, orient=tk.VERTICAL, command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )
        
        self.current_os = platform.system()
        if self.current_os == "Windows":
            self.scrollable_frame.bind("<MouseWheel>", self.scroll_off_bar)
        else:
            self.scrollable_frame.bind("<Button-4>", self.scroll_off_bar)    # Linux Scroll Up
            self.scrollable_frame.bind("<Button-5>", self.scroll_off_bar)
        
        self.selected_image_path = ''

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scroll_y.set)

        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.scroll_y.pack(side=tk.RIGHT, fill=tk.Y)

        # Initialize with image paths if provided
        self.image_paths = image_paths or []
        self.image_frames = []  # Keep track of image frames
        self.selected_frame = None  # Track the selected frame
        self.load_images(self.image_paths)
        
    def scroll_off_bar(self, event):
        """Handles mouse wheel scrolling based on OS."""
        if event.delta:  # Windows & macOS
            self.canvas.yview_scroll(-1 * (event.delta // 120), "units")
        else:  # Linux (event.num is used instead of event.delta)
            if event.num == 4:  # Scroll up
                self.canvas.yview_scroll(-1, "units")
            elif event.num == 5:  # Scroll down
                self.canvas.yview_scroll(1, "units")

    def load_images(self, image_paths):
        # Clear existing images
        for frame in self.image_frames:
            frame.destroy()
        self.image_frames.clear()

        # Create a grid with a uniform cell size
        for idx, path in enumerate(image_paths):
            # Create a cell frame with a border and fixed size
            cell_frame = tk.Frame(self.scrollable_frame, width=self.cell_width, height=self.cell_height,
                                  borderwidth=2, relief="solid", bg="white")
            cell_frame.grid_propagate(False)  # Prevent resizing of the frame
            cell_frame.grid(row=idx // self.columns, column=idx % self.columns, padx=0, pady=0)  # No padding

            # Load the image and resize it to fit inside the cell
            img = Image.open(path).convert("RGBA")  # Convert to RGBA for transparency handling
            img.thumbnail((self.cell_width - 20, self.cell_height - 20))  # Resize for the cell

            # Create a transparent background and paste the image with alpha transparency
            transparent_bg = Image.new("RGBA", (self.cell_width - 20, self.cell_height - 20), (255, 255, 255, 0))
            transparent_bg.paste(img, (0, 0), img)

            photo = ImageTk.PhotoImage(transparent_bg)

            # Create a label for the image and center it in the cell
            label = tk.Label(cell_frame, image=photo, bg="white")
            label.image = photo  # Keep a reference to avoid garbage collection
            label.place(relx=0.5, rely=0.5, anchor="center")  # Center the image in the cell

            # Bind a click event to both the frame and the label (image)
            cell_frame.bind("<Button-1>", lambda e, idx=idx: self.select_image(idx))
            label.bind("<Button-1>", lambda e, idx=idx: self.select_image(idx))
            
            if self.current_os == "Windows":
                cell_frame.bind("<MouseWheel>", self.scroll_off_bar)
                label.bind("<MouseWheel>", self.scroll_off_bar)
            else:
                cell_frame.bind("<Button-4>", self.scroll_off_bar)    # Linux Scroll Up
                cell_frame.bind("<Button-5>", self.scroll_off_bar)
                label.bind("<Button-4>", self.scroll_off_bar)    # Linux Scroll Up
                label.bind("<Button-5>", self.scroll_off_bar)

            self.image_frames.append(cell_frame)  # Store the frame for future reference

    def select_image(self, idx):
        # Reset the background color of the previously selected frame
        if self.selected_frame:
            try:
                self.selected_frame.config(bg="white")
            except:
                pass

        # Set the background of the new selected frame to light blue
        self.selected_frame = self.image_frames[idx]
        self.selected_frame.config(bg="lightblue")

        # Print the selected image path (optional)
        self.selected_image_path = self.image_paths[idx]
        
    def get_selected(self):
        return self.selected_image_path

    def update_images(self, new_image_paths):
        self.image_paths = new_image_paths
        self.load_images(self.image_paths)

class CustomToolbar(NavigationToolbar2Tk):
    def __init__(self, canvas, window):
        super().__init__(canvas, window)
        self.canvas = canvas
        self.ax = self.canvas.figure.gca()
        self.ruler_on = False
        self.press_event = None
        self.line = None
        self.distance_text = None
        
        # Add custom 'Ruler' button to the toolbar
        self.add_ruler_button()

    def add_ruler_button(self):
        """Add a ruler button to the toolbar."""
        self.ruler_button = tk.Button(self, image=photo_ruler, relief=tk.FLAT, command=self.toggle_ruler)
        self.ruler_button.pack(side=tk.LEFT, padx=2, pady=2)
        ToolTip(self.ruler_button, "Measure the distance between peaks, vertically for m/z or horizontally for retention/migration time.")

    def toggle_ruler(self):
        """Toggle the ruler tool on or off and ensure exclusivity with pan/zoom."""
        self.ruler_on = not self.ruler_on
        if self.ruler_on:
            self.ruler_button.config(relief=tk.SUNKEN)
            # Deactivate pan and zoom tools
            self._button_click('pan')  # Disable pan
            self._button_click('zoom')  # Disable zoom
            # Connect mouse events for ruler
            self._id_press = self.canvas.mpl_connect('button_press_event', self.on_mouse_press)
            self._id_motion = self.canvas.mpl_connect('motion_notify_event', self.on_mouse_drag)
            self._id_release = self.canvas.mpl_connect('button_release_event', self.on_mouse_release)
        else:
            self.ruler_button.config(relief=tk.FLAT)
            # Disconnect mouse events when ruler is off
            self.canvas.mpl_disconnect(self._id_press)
            self.canvas.mpl_disconnect(self._id_motion)
            self.canvas.mpl_disconnect(self._id_release)
            self.clear_ruler()

    def _button_click(self, tool):
        """Deactivate the specified toolbar tool (e.g., pan, zoom)."""
        # The toolbar has pan and zoom buttons, so we need to turn them off
        if tool == 'pan' and self.mode == 'pan/zoom':
            super().pan()  # This turns off pan if it’s already active
        if tool == 'zoom' and self.mode == 'zoom rect':
            super().zoom()  # This turns off zoom if it’s already active

    def pan(self, *args):
        """Override the pan method to deactivate ruler when pan is activated."""
        if self.ruler_on:
            self.toggle_ruler()  # Deactivate ruler when pan is activated
        super().pan(*args)

    def zoom(self, *args):
        """Override the zoom method to deactivate ruler when zoom is activated."""
        if self.ruler_on:
            self.toggle_ruler()  # Deactivate ruler when zoom is activated
        super().zoom(*args)

    def on_mouse_press(self, event):
        """Store the initial click event."""
        if event.inaxes != self.ax:
            return
        self.press_event = event
        self.clear_ruler()

    def on_mouse_drag(self, event):
        """Draw a line constrained to horizontal or vertical movement."""
        if event.inaxes != self.ax or self.press_event is None:
            return

        self.clear_ruler()  # Clear any previous lines or text

        # Determine drag direction and draw the line accordingly
        x_range = self.ax.get_xlim()[1] - self.ax.get_xlim()[0]
        y_range = self.ax.get_ylim()[1] - self.ax.get_ylim()[0]
        if abs(event.xdata - self.press_event.xdata)/(x_range) > abs(event.ydata - self.press_event.ydata)/y_range:
            # Horizontal line from click to current mouse position
            self.line, = self.ax.plot([self.press_event.xdata, event.xdata],
                                      [self.press_event.ydata, self.press_event.ydata],
                                      color='#94f731', linestyle='-')
            distance = abs(event.xdata - self.press_event.xdata)
        else:
            # Vertical line from click to current mouse position
            self.line, = self.ax.plot([self.press_event.xdata, self.press_event.xdata],
                                      [self.press_event.ydata, event.ydata],
                                      color='#94f731', linestyle='-')
            distance = abs(event.ydata - self.press_event.ydata)

        # Display the distance as a text on the plot
        self.distance_text = self.ax.text(0, 1.05, f'Distance: {distance:.2f}',
                                          transform=self.ax.transAxes,
                                          fontsize=12, color='black',
                                          verticalalignment='top')

        self.canvas.draw()

    def on_mouse_release(self, event):
        """Reset the press event and finalize the line."""
        self.press_event = None

    def clear_ruler(self):
        """Clear any existing lines or distance text."""
        if self.line is not None:
            self.line.remove()
            self.line = None
        if self.distance_text is not None:
            self.distance_text.remove()
            self.distance_text = None

        self.canvas.draw()

class TreeviewTooltip:
    def __init__(self, tree, delay=500):
        self.tree = tree
        self.delay = delay  # delay in milliseconds
        self.tipwindow = None
        self.after_id = None
        self.last_item = ""

        self.tree.bind("<Motion>", self.on_motion)
        self.tree.bind("<Leave>", self.on_leave)

    def on_motion(self, event):
        item = self.tree.identify_row(event.y)
        if item != self.last_item:
            self.last_item = item
            self.cancel_scheduled_tooltip()
            self.schedule_tooltip(item, event.x_root + 10, event.y_root + 10)

    def on_leave(self, event=None):
        self.cancel_scheduled_tooltip()
        self.hide_tooltip()

    def schedule_tooltip(self, item, x, y):
        if item:
            self.after_id = self.tree.after(self.delay, lambda: self.show_tooltip(item, x, y))

    def cancel_scheduled_tooltip(self):
        if self.after_id:
            self.tree.after_cancel(self.after_id)
            self.after_id = None
        self.hide_tooltip()

    def show_tooltip(self, item, x, y):
        self.hide_tooltip()
        text = self.tree.item(item, "text")
        self.tipwindow = tw = tk.Toplevel(self.tree)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(tw, text=text, background="white", justify=tk.LEFT, relief="solid", borderwidth=1, wraplength=200)
        label.pack()

    def hide_tooltip(self):
        if self.tipwindow:
            self.tipwindow.destroy()
            self.tipwindow = None

class ToolTip:
    '''This allows to create the mouse hover tooltips'''
    def __init__(self, widget, text, delay=750):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.tip_window = None
        self.entered = False
        self.widget.bind("<Enter>", self.on_enter)
        self.widget.bind("<Leave>", self.on_leave)
        
        # Store tooltip in widget dictionary
        self.widget.tooltip = self  # Assign tooltip to widget

    def on_enter(self, event):
        self.entered = True
        self.widget.after(self.delay, self.show_tooltip)

    def on_leave(self, event):
        self.entered = False
        self.hide_tooltip()

    def show_tooltip(self):
        if self.entered:
            if not self.tip_window:
                # Get screen width
                screen_width = self.widget.winfo_screenwidth()

                # Get the mouse position
                x = self.widget.winfo_pointerx()+15
                y = self.widget.winfo_pointery()+15

                # Calculate width of the tooltip window to check if it goes off the screen
                tooltip_width = 200  # Width of the tooltip (wraplength)
                if x + tooltip_width > screen_width:  # Tooltip goes off the screen
                    x = self.widget.winfo_pointerx()-15 - tooltip_width  # Move to the left side

                self.tip_window = tk.Toplevel(self.widget)
                self.tip_window.attributes("-topmost", True)
                self.tip_window.wm_overrideredirect(True)
                self.tip_window.wm_geometry(f"+{x}+{y}")
                label = tk.Label(self.tip_window, text=self.text, justify=tk.LEFT,
                                 background="white", relief=tk.SOLID, borderwidth=1, wraplength=200)
                label.pack(ipadx=1)

    def hide_tooltip(self):
        if self.tip_window:
            self.tip_window.destroy()
            self.tip_window = None           
            
class TextRedirector_Gen_Lib(object): 
    '''This redirects the text from the sys.out to the text widget for the library generation progress window'''               
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, s):
        found = False
        for i in suppressed_prints:
            if i in s:
                found = True
                break
        if "File name is" in s:
            global library_name
            ok_lib_gen_button.config(state=tk.NORMAL)
            library_name = s.split("'")[1].split(".")[0]
        if not found:
            self.widget.config(state=tk.NORMAL)
            self.widget.insert(tk.END, s, (self.tag,))
            self.widget.config(state=tk.DISABLED)
            self.widget.see(tk.END)  # Autoscroll to the end
        
    def flush(self):
        pass
                    
class TextRedirector_Run_Analysis(object):
    '''This redirects the text from the sys.out to the text widget for the analysis progress window'''         
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, s):
        global which_progress_bar, progress_bar_run_analysis, progress_bar_ms1_label, progress_bar_run_analysis2, progress_bar_ms2_label, multithreaded_analysis, number_cores
        found = False
        for i in suppressed_prints:
            if i in s:
                found = True
                break
        if "MS1" in s:
            which_progress_bar = 'ms1'
        if "MS2" in s:
            which_progress_bar = 'ms2'
        if "MS1 tracing done" in s:
            progress_bar_run_analysis["value"] = 100
            progress_bar_ms1_label.config(text=f"MS1 Analysis Progress: {100}%")
        if "MS2 analysis done" in s:
            progress_bar_run_analysis2["value"] = 100
            progress_bar_ms2_label.config(text=f"MS2 Analysis Progress: {100}%")
        if "Finished" in s:
            ok_run_analysis_button.config(state=tk.NORMAL)
        if "Traced " in s or "Analyzed glycan " in s:
            progress_value = "%.1f" % round(((int(s.split(" ")[-1].split("/")[0]))/int(s.split(" ")[-1].split("/")[1]))*100, 1)
            if which_progress_bar == 'ms1':
                progress_bar_run_analysis["value"] = progress_value
                progress_bar_ms1_label.config(text=f"MS1 Analysis Progress: {progress_value}%")
            elif which_progress_bar == 'ms2':
                progress_bar_run_analysis2["value"] = progress_value
                progress_bar_ms2_label.config(text=f"MS2 Analysis Progress: {progress_value}%")
            run_analysis.update_idletasks()
        if "File name is" in s:
            global gg_file_name
            gg_file_name = s.split("'")[1]
        if not found:
            self.widget.config(state=tk.NORMAL)
            self.widget.insert(tk.END, s, (self.tag,))
            self.widget.config(state=tk.DISABLED)
            self.widget.see(tk.END)  # Autoscroll to the end
        
    def flush(self):
        pass
        
class TextRedirector_save_result(object):  
    '''This redirects the text from the sys.out to the text widget for the save results progress window'''                       
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, s):
        found = False
        for i in suppressed_prints:
            if i in s:
                found = True
                break
        if "Finished" in s:
            ok_save_result_button.config(state=tk.NORMAL)
        if not found:
            self.widget.config(state=tk.NORMAL)
            self.widget.insert(tk.END, s, (self.tag,))
            self.widget.config(state=tk.DISABLED)
            self.widget.see(tk.END)  # Autoscroll to the end
        
    def flush(self):
        pass
        
# Functions to manage file association with .gg files
def check_file_association(extension, file_type):
    """
    Check if a file association already exists for the given extension and file type.
    """
    try:
        # Check if the extension key exists
        extension_key = f"SOFTWARE\\Classes\\{extension}"
        with reg.OpenKey(reg.HKEY_CURRENT_USER, extension_key) as key:
            associated_file_type, _ = reg.QueryValueEx(key, "")
            if associated_file_type != file_type:
                return False, 1
        
        # Check if the file type key exists and the associated script is valid
        command_key = f"SOFTWARE\\Classes\\{file_type}\\shell\\open\\command"
        with reg.OpenKey(reg.HKEY_CURRENT_USER, command_key) as key:
            command, _ = reg.QueryValueEx(key, "")
            # Extract the .exe path from the command string (strip quotes and arguments)
            if len(command.split('" "')) > 2:
                exe_path = command.split('" "')[1].strip('"')
            else:
                exe_path = command.split('" "')[0].strip('"')
            
            if not os.path.exists(exe_path):
                return False, 2
            if exe_path != os.path.abspath(sys.argv[0]):
                return False, 3

        return True, 0
    except FileNotFoundError:
        return False, 4
        
def remove_file_association(extension, file_type):
    """
    Remove the file association for the given extension and file type.
    """
    try:
        # Define registry keys
        extension_key = f"SOFTWARE\\Classes\\{extension}"
        filetype_key = f"SOFTWARE\\Classes\\{file_type}"

        # Delete the extension key
        try:
            reg.DeleteKey(reg.HKEY_CURRENT_USER, extension_key)
        except FileNotFoundError:
            error_window(f"Registry key not found: {extension_key}")

        # Delete the file type key (recursively)
        def delete_recursively(key_path):
            try:
                with reg.OpenKey(reg.HKEY_CURRENT_USER, key_path, 0, reg.KEY_READ) as key:
                    subkeys_count, _, _ = reg.QueryInfoKey(key)
                    for i in range(subkeys_count):
                        subkey_name = reg.EnumKey(key, 0)  # Always delete the first subkey
                        delete_recursively(f"{key_path}\\{subkey_name}")
                reg.DeleteKey(reg.HKEY_CURRENT_USER, key_path)
            except FileNotFoundError:
                error_window(f"Registry key not found: {key_path}")

        delete_recursively(filetype_key)
    except Exception as e:
        error_window(f"Failed to remove file association: {e}")

def create_file_association(file_extension, file_type, description):
    try:
        # Paths
        python_executable = sys.executable
        script_path = os.path.abspath(sys.argv[0])  # Path to your Python script

        # Command to execute the file (using pythonw for no terminal window)
        if script_path.split(".")[-1] == "exe":
            command = f'"{script_path}" "%1"'
        else:
            command = f'"{python_executable}" "{script_path}" "%1"'

        # Registry keys
        extension_key = f"SOFTWARE\\Classes\\{file_extension}"
        filetype_key = f"SOFTWARE\\Classes\\{file_type}"
        command_key = f"{filetype_key}\\shell\\open\\command"
        default_icon_key = f"{filetype_key}\\DefaultIcon"

        # Create the extension key and associate it with the file type
        with reg.CreateKey(reg.HKEY_CURRENT_USER, extension_key) as key:
            reg.SetValue(key, "", reg.REG_SZ, file_type)

        # Create the file type key and set its description
        with reg.CreateKey(reg.HKEY_CURRENT_USER, filetype_key) as key:
            reg.SetValue(key, "", reg.REG_SZ, description)

        # Set the command to execute when the file is opened
        with reg.CreateKey(reg.HKEY_CURRENT_USER, command_key) as key:
            reg.SetValue(key, "", reg.REG_SZ, command)

        # Set the default icon for the file type
        with reg.CreateKey(reg.HKEY_CURRENT_USER, default_icon_key) as key:
            reg.SetValue(key, "", reg.REG_SZ, ico_path)
        
        confirmation_window(f"File association with {file_extension} files created succesfully!")
        
    except Exception as e:
        error_window(f"Failed to create file association: {e}")
        
# Custom square root function that handles negative values
def custom_sqrt(x):
    return np.where(x >= 0, np.sqrt(x), -np.sqrt(-x))

# Custom square function (inverse of sqrt)
def custom_square(x):
    return np.where(x >= 0, np.square(x), -np.square(-x))
            
def der_sialics_combo(mono, s_value):
    combinations = []
    
    if mono == 'S':
        # Loop over all possible values of Am
        for am in range(s_value + 1):
            e = s_value - am  # E is the remaining part that sums with Am to S
            combinations.append({'Am': am, 'E': e})
    
    if mono == 'G':
        # Loop over all possible values of Am
        for am in range(s_value + 1):
            e = s_value - am  # E is the remaining part that sums with Am to S
            combinations.append({'AmG': am, 'EG': e})
    
    return combinations
    
def replace_dict_keys(original_dict, original_key, replacement_dict):
    # Check if 'S' is in the original dict
    if original_key in original_dict:
        s_value = original_dict[original_key]
        del original_dict[original_key]
        
        for new_key, new_value in replacement_dict.items():
            original_dict[new_key] = new_value
    
    return original_dict
        
def make_total_glycans_df(df1,
                          df2,
                          percentage_auc = percentage_auc,
                          rt_tolerance = ret_time_interval[2]):
    '''This is used here for the real-time alignment when comparing samples. DEPRECATED'''
    global s_to_n, iso_fit_score, curve_fit_score
    sn = s_to_n
    df1 = copy.deepcopy(df1)
    df1_refactor = []
    for i_i, i in enumerate(df2["Sample_Number"]): #QCs cutoff
        df1_refactor.append({"Glycan" : [], "Adduct" : [], "mz" : [], "RT" : [], "AUC" : [], "PPM" : [], "S/N" : [], "Iso_Fitting_Score" : [], "Curve_Fitting_Score" : []})
        for j_j, j in enumerate(df1[i_i]["Adduct"]): 
            temp_rt = df1[i_i]["RT"][j_j]
            temp_auc = df1[i_i]["AUC"][j_j]
            temp_ppm = df1[i_i]["PPM"][j_j]
            temp_sn = df1[i_i]["S/N"][j_j]
            temp_fit = df1[i_i]["Iso_Fitting_Score"][j_j]
            temp_curve = df1[i_i]["Curve_Fitting_Score"][j_j]
            to_remove = []
            to_remove_glycan = []
            to_remove_adduct = []
            for k_k, k in enumerate(temp_sn):
                if df1[i_i]["Glycan"][j_j] != "Internal Standard":
                    if float(k) < sn:
                        to_remove.append(k_k)
                        to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                        to_remove_adduct.append(j)
                        continue
                    if float(temp_fit[k_k]) < iso_fit_score:
                        to_remove.append(k_k)
                        to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                        to_remove_adduct.append(j)
                        continue
                    if float(temp_curve[k_k]) < curve_fit_score:
                        to_remove.append(k_k)
                        to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                        to_remove_adduct.append(j)
                        continue
                    if float(temp_ppm[k_k]) < max_ppm[0] or float(temp_ppm[k_k]) > max_ppm[1]:
                        to_remove.append(k_k)
                        to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                        to_remove_adduct.append(j)
                        continue
            if len(to_remove) != 0:
                to_remove.reverse()
                to_remove_glycan.reverse()
                to_remove_adduct.reverse()
                for k_k in to_remove:
                    del temp_rt[k_k]
                    del temp_auc[k_k]
                    del temp_ppm[k_k]
                    del temp_sn[k_k]
                    del temp_fit[k_k]
                    del temp_curve[k_k]
            to_remove = [] #second pass to remove based on % of remained peaks
            to_remove_glycan = []
            to_remove_adduct = []  
            for k_k, k in enumerate(temp_sn): 
                if max(temp_auc) == 0.0:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
                if float(temp_auc[k_k]/max(temp_auc)) <= percentage_auc:
                    to_remove.append(k_k)
                    to_remove_glycan.append(df1[i_i]["Glycan"][j_j])
                    to_remove_adduct.append(j)
                    continue
            if len(to_remove) != 0:
                to_remove.reverse()
                to_remove_glycan.reverse()
                to_remove_adduct.reverse()
                for k_k in to_remove:
                    del temp_rt[k_k]
                    del temp_auc[k_k]
                    del temp_ppm[k_k]
                    del temp_sn[k_k]
                    del temp_fit[k_k]
                    del temp_curve[k_k]
            df1[i_i]["RT"][j_j] = str(temp_rt)[1:-1]
            df1[i_i]["AUC"][j_j] = str(temp_auc)[1:-1]
            df1[i_i]["PPM"][j_j] = str(temp_ppm)[1:-1]
            df1[i_i]["S/N"][j_j] = str(temp_sn)[1:-1]
            df1[i_i]["Iso_Fitting_Score"][j_j] = str(temp_fit)[1:-1]
            df1[i_i]["Curve_Fitting_Score"][j_j] = str(temp_curve)[1:-1]
    
    to_remove = []
    to_remove_glycan = []
    to_remove_adduct = []
    for i_i, i in enumerate(df1[0]["Adduct"]):
        glycan_good = False
        for j_j, j in enumerate(df2["Sample_Number"]):
            check = df1[j_j]["RT"][i_i]
            if len(check) != 0:
                glycan_good = True
                break
        if not glycan_good:
            to_remove.append(i_i)
            to_remove_glycan.append(df1[j_j]["Glycan"][i_i])
            to_remove_adduct.append(i)
    if len(to_remove) != 0:
        to_remove.reverse()
        to_remove_glycan.reverse()
        to_remove_adduct.reverse()
        for j_j, j in enumerate(df2["Sample_Number"]):
            for i_index, i_i in enumerate(to_remove):
                del df1[j_j]["Glycan"][i_i]
                del df1[j_j]["Adduct"][i_i]
                del df1[j_j]["mz"][i_i]
                del df1[j_j]["RT"][i_i]
                del df1[j_j]["AUC"][i_i]
                del df1[j_j]["PPM"][i_i]
                del df1[j_j]["S/N"][i_i]
                del df1[j_j]["Iso_Fitting_Score"][i_i]
                del df1[j_j]["Curve_Fitting_Score"][i_i]
                            
    for i_i, i in enumerate(df1): #final arrangement for standard results print
        for j_j, j in enumerate(df1[i_i]["Adduct"]):
            for k_k, k in enumerate(df1[i_i]["RT"][j_j].split(", ")):
                if k != "":
                    df1_refactor[i_i]["Glycan"].append(df1[i_i]["Glycan"][j_j])
                    df1_refactor[i_i]["Adduct"].append(df1[i_i]["Adduct"][j_j])
                    df1_refactor[i_i]["mz"].append(df1[i_i]["mz"][j_j])
                    df1_refactor[i_i]["RT"].append(float(k))
                    df1_refactor[i_i]["AUC"].append(float(df1[i_i]["AUC"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["PPM"].append(float(df1[i_i]["PPM"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["S/N"].append(float(df1[i_i]["S/N"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["Iso_Fitting_Score"].append(float(df1[i_i]["Iso_Fitting_Score"][j_j].split(", ")[k_k]))
                    df1_refactor[i_i]["Curve_Fitting_Score"].append(float(df1[i_i]["Curve_Fitting_Score"][j_j].split(", ")[k_k]))        
    total_dataframes = [] #total glycans AUC dataframe
    for i_i, i in enumerate(df1_refactor):
        total_dataframes.append({"Glycan": [], "RT": [], "AUC": []})
        current_glycan = ""
        RTs = []
        AUCs = []
        for j_j in range(len(i["Glycan"])):
            j = i["Glycan"][j_j]
            found = False
            if j == current_glycan:
                for k_k, k in enumerate(RTs):
                    if abs(i["RT"][j_j] - k) <= rt_tolerance:
                        RTs[k_k] = (k+i["RT"][j_j])/2
                        AUCs[k_k] = AUCs[k_k]+i["AUC"][j_j]
                        found = True
                        break
                if not found:
                    RTs.append(i["RT"][j_j])
                    AUCs.append(i["AUC"][j_j])
            if j != current_glycan:
                if j_j != 0:
                    for k_k, k in enumerate(RTs):
                            total_dataframes[i_i]["Glycan"].append(current_glycan)
                            total_dataframes[i_i]["RT"].append(k)
                            total_dataframes[i_i]["AUC"].append(AUCs[k_k])
                    RTs = []
                    AUCs = []
                current_glycan = j
                RTs.append(i["RT"][j_j])
                AUCs.append(i["AUC"][j_j])
            if j_j == len(i["Glycan"])-1:
                for k_k, k in enumerate(RTs):
                    total_dataframes[i_i]["Glycan"].append(current_glycan)
                    total_dataframes[i_i]["RT"].append(k)
                    total_dataframes[i_i]["AUC"].append(AUCs[k_k])
                RTs = []
                AUCs = [] #total glycans AUC dataframe
    
    arranged_total_dataframes = []
    for i_i, i in enumerate(total_dataframes):
        arranged_total_dataframes.append({'Glycan' : [], 'RT' : [], 'AUC' : []})
        list_of_glycans = []
        glycans_dict_rt = {}
        glycans_dict_AUC = {}
        current_glycan = ''
        for j_j, j in enumerate(i['Glycan']):
            if j != current_glycan:
                current_glycan = j
                list_of_glycans.append(j)
                glycans_dict_rt[j] = []
                glycans_dict_AUC[j] = []
                glycans_dict_rt[j].append(i['RT'][j_j])
                glycans_dict_AUC[j].append(i['AUC'][j_j])
            else:
                glycans_dict_rt[j].append(i['RT'][j_j])
                glycans_dict_AUC[j].append(i['AUC'][j_j])
        if "Internal Standard" in list_of_glycans:
            highest = glycans_dict_AUC["Internal Standard"].index(max(glycans_dict_AUC["Internal Standard"]))
            to_remove = []
            for j_j, j in enumerate(glycans_dict_AUC["Internal Standard"]):
                if j_j != highest:
                    to_remove.append(j_j)
            for j in sorted(to_remove, reverse = True):
                del glycans_dict_rt["Internal Standard"][j]
                del glycans_dict_AUC["Internal Standard"][j]
        list_of_glycans = sorted(list_of_glycans)
        for j_j, j in enumerate(list_of_glycans):
            for k in range(len(glycans_dict_rt[j])):
                arranged_total_dataframes[i_i]['Glycan'].append(j)
            zipped = zip(glycans_dict_rt[j], glycans_dict_AUC[j])
            zipped = sorted(zipped)
            current_RTs, current_AUCs = zip(*zipped)
            arranged_total_dataframes[i_i]['RT'] += list(current_RTs)
            arranged_total_dataframes[i_i]['AUC'] += list(current_AUCs)
    total_dataframes = arranged_total_dataframes
    return total_dataframes    
        
def pre_process(samples_list):
    '''This function starts the loading of mzml and mzxml files.'''
    global samples_names
    sample_info = {}
    
    results = []
    bad = False
    cpu_number = (os.cpu_count())-2 if os.cpu_count() < 60 else 60
    if cpu_number <= 0:
        cpu_number = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_number) as executor:
        for i_i, i in enumerate(samples_list):
            result = executor.submit(pre_process_one_sample, i, samples_names[i_i])
            results.append(result)
            
        for index, i in enumerate(results):
            result = i.result()
            if result == 'bad':
                bad = True
                break
            sample_info[result[1]] = result[0]
            results[index] = ''
            
    global processed_data
    if bad:
        processed_data = 'bad'
    else:
        processed_data = sample_info
    
def pre_process_one_sample(sample, sample_name):
    '''This is the effective part of pre_process function.'''
    
    def get_nested(dictionary, keys, fallback_keys = []):
        '''Helper function to retrieve nested keys safely.'''
        original_dict = copy.deepcopy(dictionary)
        for key in keys:
            if (type(key) == str and key in dictionary) or type(key) == int:
                dictionary = dictionary[key]
            else:
                return get_nested(original_dict, fallback_keys)
        return dictionary
    
    data = {}
    data['ms2'] = {}
    bpc = []
    rt_array = []
    ms1_array = []
    time_unit = ''
    file_type = sample.split('.')[-1].lower()
    if file_type == 'mzml':
        access = mzml.MzML(sample)
        rt_key = 'scanList'
        rt_subkey = 'scan start time'
        ms_key = 'ms level'
        intensity_key = 'base peak intensity'
        precursor_key = ['precursorList', 'precursor', 0, 'isolationWindow', 'isolation window target m/z']
        # Some MzML have isolation window, some don't, so fallback takes that into account
        fallback_precursor_key = ['precursorList', 'precursor', 0, 'selectedIonList', 'selectedIon', 0, 'selected ion m/z']
    elif file_type == 'mzxml':
        access = mzxml.MzXML(sample)
        rt_key = 'retentionTime'
        rt_subkey = None
        ms_key = 'msLevel'
        intensity_key = 'basePeakIntensity'
        precursor_key = ['precursorMz']
        
    try:
        # Determine time unit
        if file_type == 'mzml':
            if float(access[-1][rt_key]['scan'][0][rt_subkey]) > 300:
                time_unit = 'seconds'
            else:
                time_unit = 'minutes'
        elif file_type == 'mzxml':
            time_unit = 'minutes'
            
        last_ms1 = ''
        max_mz = 0
        for i_i, i in enumerate(access):
            ms_level = i[ms_key]
            current_rt = float(i[rt_key]['scan'][0][rt_subkey]) if file_type == 'mzml' else float(i[rt_key])

            if ms_level == 2:
                data['ms2'][last_ms1][current_rt] = [get_nested(i, precursor_key, fallback_precursor_key)] if file_type == 'mzml' else [k['precursorMz'] for k in get_nested(i, precursor_key)]

            elif ms_level == 1:
                ms1_array.append(i_i)
                rt_array.append(current_rt)

                if current_rt != last_ms1:
                    last_ms1 = current_rt
                    data['ms2'][current_rt] = {}

                # Get base peak intensity
                intensity_value = i.get(intensity_key, np.max(i['intensity array']) if len(i['intensity array']) > 0 else 0.0)
                bpc.append(float(intensity_value))

            # Track max m/z
            if len(i['m/z array']) > 0:
                temp_max = np.max(i['m/z array'])
                if temp_max > max_mz:
                    max_mz = temp_max

        # Populate data dictionary
        data.update({
            'time_unit': time_unit,
            'rt_array': rt_array,
            'bpc': bpc,
            'file_type': file_type,
            'access': access,
            'file_path': sample,
            'max_mz': max_mz,
            'ms1_array': ms1_array
        })
        
        return data, sample_name
    except:
        error_window(f"Something went wrong when loading the file {sample}. Check if it is an MzML/MzXML file. If it is, it might be corrupted.")
        return 'bad'

def analyze_fraction(access, ms_level, min_max_index, decimal_places = 2):
    '''This is the effective part of the Maximum Intensity Spectrum calculations.'''
    maximum = defaultdict(list)
    for i_i, i in enumerate(access[min_max_index[0]:min_max_index[1]]):
        if i[ms_level] == 1:
            mz_array = i['m/z array']
            int_array = i['intensity array']
            for j_j, j in enumerate(mz_array):
                value = round(j, decimal_places)
                maximum[value].append(int_array[j_j])
            if i_i%500 == 0:
                for j in maximum:
                    maximum[j] = [max(maximum[j])]
    for i in maximum:
        maximum[i] = [max(maximum[i])]
    return maximum

def calculate_ambiguities(df1):
    '''This function is used to calculate the number of ambiguities in a per sample basis.'''
    for i_i, i in enumerate(df1):
        i['Ambiguity'] = []
        for j in i['Glycan']:
            i['Ambiguity'].append([])
        for j_j, j in enumerate(i['Glycan']):
            glycan_j = j+'_'+i['Adduct'][j_j]
            for k_k, k in enumerate(i['Glycan'][j_j+1:]):
                k_k = j_j+k_k+1
                glycan_k = k+'_'+i['Adduct'][k_k]
                if j != k and i['mz'][j_j] == i['mz'][k_k]:
                    i['Ambiguity'][j_j].append(i['Glycan'][k_k]+'_'+i['Adduct'][k_k])
                    i['Ambiguity'][k_k].append(i['Glycan'][j_j]+'_'+i['Adduct'][j_j])
        for j_j, j in enumerate(i['Ambiguity']):
            if len(j) > 0:
                i['Ambiguity'][j_j] = ', '.join(j)
            else:
                i['Ambiguity'][j_j] = 'No'
    
def load_reanalysis(reanalysis_path):
    '''This function loads the .gg file.'''
    global glycans_per_sample, ms2_spectra_glycans, samples_dropdown_options, df1, df2, gg_file, parameters_gg, isotopic_fittings, curve_fittings, fragments_library_ms2_scores, gg_analysis_tolerance, adduct_combos_dict
    
    # Loads back data that had already been loaded before in this run or load from anew
    gg_file = gg_archive(reanalysis_path)
    
    # Tries to read the data
    try:
        # Load results table
        results = gg_file.get_results_table()
        if len(results) == 4:
            df1, df2, fragments_df = results[:3]
        else:
            df1, df2 = results[:2]
        
        # Load metadata
        parameters_gg = gg_file.get_metadata()
        
        if type(parameters_gg) != dict:
            error_window(f"This analysis file was created in a version of GlycoGenius before 1.3.0, and thus can only be opened in that or older versions.")
            return "Error"
        
    except:
        error_window(f"Something went wrong when loading the reanalysis file. Check if it is a .gg file. If it is, it might be corrupted.")
        return "Error"
    
    # Add the filename to the dropdown menu
    samples_dropdown_options = df2['File_Name']
    samples_dropdown['values'] = samples_dropdown_options
    
    # Compile the adduct combos into a cohesive dictionary for use later on
    adduct_combos_dict = {General_Functions.comp_to_formula(adduct[0]): {'charges': adduct[1], 'comp': adduct[0]} for adduct in parameters_gg['adduct combos']}
    
    # Create the dictionary to store data and calculate ambiguities
    glycans_per_sample = {}
    calculate_ambiguities(df1)
    
    # Create a dictionary to keep track of what glycans annotated a given MS2 spectrum
    ms2_spectra_glycans = {}
    
    # Goes through the data to store on the dictionary
    for i_i, i in enumerate(df1):
        # Add the sample to the glycans_per_sample dictionary
        glycans_per_sample[samples_dropdown_options[i_i]] = {} #glycans_per_sample = {'sample_0' : {}}
        
        # Add the sample to the ms2 spectra glycans dictionary
        ms2_spectra_glycans[samples_dropdown_options[i_i]] = {}
        
        # Start a variable to keep track of the last glycan analyzed
        last_glycan = ''
        
        # Go through the glycans in that particular sample
        for k_k, k in enumerate(i['Glycan']): #glycan
            # If the glycans is different from the last glycan, create it's spot in the dictionary
            if k != last_glycan:
                last_glycan = k
                glycans_per_sample[samples_dropdown_options[i_i]][k] = {} #glycans_per_sample = {'sample_0' : {'glycan_0': {}}}
                
            # Add the per-adduct data to the glycan, if it hasn't been done before for it
            if i['Adduct'][k_k] not in glycans_per_sample[samples_dropdown_options[i_i]][k]:
                # If there's no data, ignore it
                if len(i["S/N"][k_k]) == 1 and i["S/N"][k_k][0] == 0.0:
                    continue
                    
                # Add mz, retention times, aucs, PPM errors, iso scores, curve scores, sns and ambiguities of the adduct
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]] = {'mz' : i['mz'][k_k]} #glycans_per_sample = {'sample_0' : {'glycan_0': {'adduct_0' : {'mz' : 9999.99, 'peaks' : [peak_1, peak_2]}}}}
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['peaks'] = i['RT'][k_k]
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['auc'] = i['AUC'][k_k]
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ppm'] = i['PPM'][k_k]
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['iso'] = i['Iso_Fitting_Score'][k_k]
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['curve'] = i['Curve_Fitting_Score'][k_k]
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['sn'] = i["S/N"][k_k]
                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ambiguity'] = i["Ambiguity"][k_k]

                # If there's MS2 data, go through it
                if len(results) == 4:
                    # Create the dictionary entry
                    glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'] = {}
                    
                    # Initialize the variable to keep track of the last retention time checked
                    last_rt = ''
                    
                    # Go through the fragment data, glycans list
                    for l_l, l in enumerate(fragments_df[i_i]['Glycan']):
                        # If the current glycan+adduct matches the glycan+adduct that was being added to the glycans per sample dictionary, check it out
                        if l == k and fragments_df[i_i]['Adduct'][l_l] == i['Adduct'][k_k]:
                            
                            # If it's the first time this rt is being checked, create its entry and add the data
                            if fragments_df[i_i]['RT'][l_l] != last_rt:
                                # Create the base structure ([[fragment_mz], [fragment_intensity], [fragment_name], % TIC explained)
                                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'][fragments_df[i_i]['RT'][l_l]] = [[],[],[], fragments_df[i_i]['% TIC explained'][l_l]]
                                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'][fragments_df[i_i]['RT'][l_l]][0].append(fragments_df[i_i]['Fragment_mz'][l_l])
                                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'][fragments_df[i_i]['RT'][l_l]][1].append(fragments_df[i_i]['Fragment_Intensity'][l_l])
                                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'][fragments_df[i_i]['RT'][l_l]][2].append(fragments_df[i_i]['Fragment'][l_l])
                                last_rt = fragments_df[i_i]['RT'][l_l]
                                
                            # Otherwise just add the data
                            else:
                                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'][fragments_df[i_i]['RT'][l_l]][0].append(fragments_df[i_i]['Fragment_mz'][l_l])
                                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'][fragments_df[i_i]['RT'][l_l]][1].append(fragments_df[i_i]['Fragment_Intensity'][l_l])
                                glycans_per_sample[samples_dropdown_options[i_i]][k][i['Adduct'][k_k]]['ms2'][fragments_df[i_i]['RT'][l_l]][2].append(fragments_df[i_i]['Fragment'][l_l])
        
        # If MS2 is present, catalog the retention times and what glycans were annotated to it
        if len(results) == 4:
            for index, glycan in enumerate(fragments_df[i_i]['Glycan']):
                if fragments_df[i_i]['RT'][index] not in ms2_spectra_glycans[samples_dropdown_options[i_i]].keys():
                    ms2_spectra_glycans[samples_dropdown_options[i_i]][fragments_df[i_i]['RT'][index]] = {}
                if glycan not in ms2_spectra_glycans[samples_dropdown_options[i_i]][fragments_df[i_i]['RT'][index]].keys():
                    ms2_spectra_glycans[samples_dropdown_options[i_i]][fragments_df[i_i]['RT'][index]][glycan] = []
                if fragments_df[i_i]['Adduct'][index] not in ms2_spectra_glycans[samples_dropdown_options[i_i]][fragments_df[i_i]['RT'][index]][glycan]:
                    ms2_spectra_glycans[samples_dropdown_options[i_i]][fragments_df[i_i]['RT'][index]][glycan].append(fragments_df[i_i]['Adduct'][index])
    
    # Get isotopic fittings and curve fittings from file
    isotopic_fittings = gg_file.get_isotopic_fittings()
    curve_fittings = gg_file.get_curve_fittings()
    fragments_library_ms2_scores = gg_file.get_fragments_library_ms2_scores()
    gg_analysis_tolerance = parameters_gg['analysis settings']['mass tolerance']

def error_window(text):
    '''This function makes an error window with the text.'''
    messagebox.showerror("Error", text)
    
def warning_window(text):
    response = messagebox.askokcancel("Warning", text)
    return response
    
def confirmation_window(text):
    messagebox.showinfo("Confirmation", text)
        
def kill_concurrent_futures():
    '''This function finishes off all instances of GG multithreading executions.'''
    processes = psutil.pids()
    gui_id = os.getpid()
    this_process = psutil.Process(gui_id)
    this_process_name = this_process.name()
    for p in processes:
        try:
            process = psutil.Process(p)
            name = process.name()
            ppid = process.ppid()
            if (name == 'python.exe' or name == this_process_name) and ppid == gui_id:
                process.terminate()
        except:
            continue
            
def clean_temp_folder():
    '''Cleans the temporary folder off gg stuff'''
    processes = psutil.pids()
    gui_id = os.getpid()
    this_process = psutil.Process(gui_id)
    this_process_name = this_process.name()
    counter = 0
    for p in processes:
        try:
            process = psutil.Process(p)
            name = process.name()
            command_line = process.cmdline()
            for i in command_line:
                if "glycogenius" in i.lower():
                    counter+= 1
        except:
            pass
    if counter <= 1:
        general_temp_folder = os.path.join(tempfile.gettempdir())
        temp_dir_list = os.listdir(general_temp_folder)
        for i in temp_dir_list:
            if i.split("_")[0] == 'gg':
                if os.path.isdir(os.path.join(general_temp_folder, i)):
                    try:
                        shutil.rmtree(os.path.join(general_temp_folder, i))
                    except:
                        pass
    else:
        try:
            shutil.rmtree(temp_folder)
        except:
            pass
            
def annotate_ms2_spectrum(file, spectrum, fragments_dict, tolerance, target_glycan_comp):
    ''''''
    # Superscripts for pretty fragment name
    superscripts = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹', '+': '⁺', '-': '⁻', '=': '⁼', '(': '⁽', ')': '⁾', 'n': 'ⁿ', 'i': 'ⁱ'}
    
    # Annotated fragments list
    annotated_fragments = []
    
    # If empty array, return empty annotations
    if len(file[spectrum]['intensity array']) == 0:
        return annotated_fragments
    
    # Create indexed fragments dictionary
    indexed_fragments = {}
    for fragment, frag_data in fragments_dict.items():
        for adduct, mz in frag_data['Adducts_mz'].items():
            indexed_fragments[mz] = indexed_fragments.get(mz, []) + [f"{fragment}_{adduct}"]
    indexed_fragments = dict(sorted(indexed_fragments.items()))
    indexed_fragments_list = list(indexed_fragments.keys())
    
    # The array summed intensities information to calculate % TIC assigned
    total_array_intensity = sum(file[spectrum]['intensity array'])
    
    # The last peak information
    former_peak_mz = 0
    former_peak_intensity = 0
    former_peak_identified_mz = 0
    
    # The maximum intensity of the array and intensity cutoff base for array trimming
    max_int = max(file[spectrum]['intensity array'])
    intensity_cutoff = 0
    
    # Trim extremely large arrays
    max_size = 10000
    if len(file[spectrum]['intensity array']) > max_size:
        
        while np.sum(file[spectrum]['intensity array'] > intensity_cutoff) > max_size:
            intensity_cutoff += max_int*0.01
            
        mask = file[spectrum]['intensity array'] > intensity_cutoff
        mz_array = file[spectrum]['m/z array'][mask]
        int_array = file[spectrum]['intensity array'][mask]
        
    else:
        mz_array = file[spectrum]['m/z array']
        int_array = file[spectrum]['intensity array']
        
    # Check the spectrum peak by peak
    for mz_peak_index, mz_peak in enumerate(mz_array):
        
        # Moving intensity threshold to ignore minuscule peaks that are between isotopologues
        if int_array[mz_peak_index] < former_peak_intensity*0.05:
            continue
            
        # Check whether this mz peak is part of the isotopic envelope of the last checked peak. Checks for possibly singly, doubly or triply charged
        if (
            abs(mz_peak-(former_peak_mz+General_Functions.h_mass)) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
            or abs(mz_peak-(former_peak_mz+(General_Functions.h_mass/2))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
            or abs(mz_peak-(former_peak_mz+(General_Functions.h_mass/3))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak)
            ): 
            # And here it checks if it is part of the isotopic envelope of the last IDENTIFIED peak. If it is, it reduces its intensity from the total array intensity so to better represent the % TIC annotated
            if (
                abs(mz_peak-(former_peak_identified_mz+General_Functions.h_mass)) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
                or abs(mz_peak-(former_peak_identified_mz+(General_Functions.h_mass/2))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak) 
                or abs(mz_peak-(former_peak_identified_mz+(General_Functions.h_mass/3))) < General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak)
                ):
                former_peak_identified_mz = mz_peak
                total_array_intensity -= int_array[mz_peak_index]
                
            former_peak_mz = mz_peak
            continue
            
        former_peak_mz = mz_peak
        former_peak_intensity = int_array[mz_peak_index]
        
        # Checks whether this peak matches the mz of a fragment
        fragment_id = General_Functions.binary_search_with_tolerance(indexed_fragments_list, mz_peak, 0, len(indexed_fragments_list)-1, General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz_peak))
        
        # If no matching fragment is found, continue to the next mz peak
        if fragment_id == -1:
            continue
            
        # Identify the fragment mz
        identified_fragment_mz = indexed_fragments_list[fragment_id]
        
        # Otherwise, fetch the list of fragments that match the mz and, if the user wants to limit the output to the composition of precursor, filter it
        possible_fragments = indexed_fragments[identified_fragment_mz]
        
        if target_glycan_comp:
            to_remove = []
            for index, fragment in enumerate(possible_fragments):
                for mono, amount in fragments_dict[fragment.split("_")[0]]['Monos_composition'].items():
                    if mono != 'T' and amount > target_glycan_comp.get(mono, 0):
                        to_remove.append(index)
                        break
            
            for index in sorted(to_remove, reverse=True):
                del possible_fragments[index]
        
        # If after filtering, no possible fragments are found, move on to the next peak
        if len(possible_fragments) == 0:
            continue
            
        # Mark this as a successful identification and save the mz as the last identified
        former_peak_identified_mz = mz_peak 
        
        # Format the fragments name to look prettier
        fragment_name_list = []
        for fragment in possible_fragments:
            fragment_name, fragment_adduct = fragment.split("_")
            
            adduct_comp_frag = adduct_combos_dict[fragment_adduct]['comp']
            adduct_charge_frag = adduct_combos_dict[fragment_adduct]['charges']
            
            adduct_str = ""
            for o in adduct_comp_frag:
                polarity = '+' if adduct_comp_frag[o] > 0 else ''
                adduct_str += f"{polarity}{adduct_comp_frag[o]}{o}"
                
            superscript_polarity = superscripts['+'] if adduct_charge_frag > 0 else superscripts['-']
            
            fragment_name_list.append(f"{fragment_name}[M{adduct_str}]{superscript_polarity}{superscripts[str(abs(adduct_charge_frag))]}")
            
        fragment_name = "/".join(fragment_name_list)
            
        annotated_fragments.append(['Unknown Glycan', 'Unknown Adduct', fragment_name, mz_peak, int_array[mz_peak_index], 0.0, 0.0, total_array_intensity])
        
    # Update the total intensities after it's done with the array
    for fragment in annotated_fragments:
        fragment[7] = total_array_intensity
            
    return annotated_fragments
    
def on_closing():
    '''This function is used to remove the function from the Close Window button (x button).'''
    return
        
def handle_selection(event):
    '''This function handles the selection of the samples dropdown menu.'''
    global current_data, ax, canvas, ax_spec, canvas_spec, samples_dropdown, chromatograms_list, selected_item, compare_samples_button, loading_files, filter_list, samples_list, spectra_plot_frame, no_spectra_loaded_label, ms1_bound, ms1_binds, ms2_bound, ms2_binds, chromatogram_bound, chromatogram_binds, rt_label, precursor_label, gg_draw_on, former_selected_sample, last_plotted, selected_chromatograms, previous_spectrum_button, next_spectrum_button, annotate_spectrum_button
    
    selected_item = samples_dropdown.get()
    if selected_item == former_selected_sample:
        return
        
    former_selected_sample = selected_item
    
    # Deactivate the buttons to cycle through MS2 spectra
    previous_spectrum_button.config(state=tk.DISABLED)
    next_spectrum_button.config(state=tk.DISABLED)
    annotate_spectrum_button.config(state=tk.DISABLED)
    
    # Clear glycans selected
    last_plotted = None
    selected_chromatograms = []
    
    # Unbind MS1 plot
    if ms1_bound:
        for i in ms1_binds:
            canvas_spec.mpl_disconnect(i)
        ms1_bound = False
        ms1_binds = []
        
    # Unbind MS2 plot
    if ms2_bound:
        for i in ms2_binds:
            canvas_spec_ms2.mpl_disconnect(i)
        ms2_bound = False
        ms2_binds = []
        
    # Unbind chromatogram plot
    if chromatogram_bound:
        for i in chromatogram_binds:
            canvas.mpl_disconnect(i)
        chromatogram_bound = False
        chromatogram_binds = []
        
    # Clear labels
    rt_label.config(text="")
    precursor_label.config(text="")
    coordinate_label_spec.config(text='')
    
    compare_samples_button.config(state=tk.DISABLED)
    ToolTip(samples_dropdown, f"{selected_item}")
    if len(reanalysis_path) > 0 and len(glycans_per_sample[selected_item]) > 0:
        check_qc_dist_button.config(state=tk.NORMAL)
    else:
        check_qc_dist_button.config(state=tk.DISABLED)
        
    if len(samples_list) > 0:
        if processed_data == 'bad':
            loading_files.destroy()
            samples_dropdown['values'] = []
            samples_dropdown.set('')
            chromatograms_list.delete(*chromatograms_list.get_children())
            samples_list = []
            return
        else:
            if selected_item in processed_data:
                current_data = processed_data[selected_item]
    populate_treeview()
    clear_plot(ax, canvas)
    clear_plot(ax_spec, canvas_spec)
    clear_plot(ax_spec_ms2, canvas_spec_ms2)
    
    
    check_folder = os.path.join(tempfile.gettempdir())
    
    this_process_id = os.getpid()
    
    if f"mzml_window_{this_process_id}.txt" in os.listdir(check_folder):
        mzml_window_start(change_sample = selected_item)
    
    no_spectra_loaded_label_ms2.config(text=f"No MS2 spectra selected. If your raw data file\ncontains MS2 data, click on a purple diamond on\nan MS1 spectra or if you've analyzed your data\nwith MS2 analysis included, click on 'MS2'\nbesides a glycan peak retention time in the\nglycan's list.")
    no_spectra_loaded_label_ms2.place(relx=0.50, rely=0.45)
    
    if len(reanalysis_path) > 0 and selected_item not in samples_names:
        no_spectra_loaded_label.config(text=f"No MzML or MzXML file loaded for this sample.\n To visualize spectra for this sample, load\n{selected_item}.MzML or\n{selected_item}.MzXML\nin the Select Files menu.")
        no_spectra_loaded_label.place(relx=0.50, rely=0.45)
    else:
        no_spectra_loaded_label.config(text="")
        no_spectra_loaded_label.place(relx=0.01, rely=0.9)
        
def handle_selection_vigilance():
    check_folder = os.path.join(tempfile.gettempdir())
    this_process_id = os.getpid()
    if f"mzml_window_{this_process_id}.txt" in os.listdir(check_folder):
        handle_selection(None)
        main_window.after(100, handle_selection_vigilance)
                        
def check_peak_quality(**kwargs):
    '''max_ppm, iso_fit_score, curve_fit_score, s_to_n'''
    count = 0
    if kwargs['ppm'] > max_ppm[1] or kwargs['ppm'] < max_ppm[0]:
        count+= 1
    if kwargs['iso'] < iso_fit_score:
        count+= 1
    if kwargs['curve'] < curve_fit_score:
        count+= 1
    if kwargs['sn'] < s_to_n:
        count+= 1
        
    return count
    
def populate_treeview():
    '''This function populates the glycans list with the glycans from the .gg file.'''
    global selected_item, chromatograms_list, filter_list
    
    # Save the current view and opened rows
    view_top, view_bottom = chromatograms_list.yview()
    len_list = len(chromatograms_list.get_children())
    open_close_state = []
    counter = 0
    for item_id in chromatograms_list.get_children():
        open_close_state.append(chromatograms_list.item(item_id, "open"))
        counter+= 1
        for sub_item_id in chromatograms_list.get_children(item_id):
            open_close_state.append(chromatograms_list.item(sub_item_id, "open"))
            counter+= 1
    
    # Clean the list
    chromatograms_list.delete(*chromatograms_list.get_children())
    
    # Grab the query
    filter_query = [query.lower() for query in filter_list.get().split("+") if query != "Filter the list of glycans..."]
    
    # Initialize the counters
    good_count = 0
    average_count = 0
    bad_count = 0
    ambiguity_count = 0
    noted_ambiguities = set()
    
    # Add the BPC
    if len(samples_list) > 0:
        if selected_item in processed_data:
            if len(filter_query) == 0:
                chromatograms_list.insert("", "end", text="Base Peak Chromatogram/Electropherogram")
            else:
                for query in filter_query:
                    if query in "base peak chromatogram/electropherogram;bpc;bpe":
                        chromatograms_list.insert("", "end", text="Base Peak Chromatogram/Electropherogram")
                        break
    
    # Add and filter the glycans
    if len(reanalysis_path) > 0:
        sample_index = list(glycans_per_sample.keys()).index(selected_item)
        for glycan, adducts in glycans_per_sample[selected_item].items():
            if len(adducts) == 0:
                continue
                
            temp_adducts_quality = {}
            ambiguity = False
            glycan_value = ''
            
            queriables = [glycan]
            
            for adduct, adduct_data in adducts.items():
                queriables.append(adduct)
                queriables.append(str(adduct_data['mz']))
                
                if adduct_data.get('ms2', ''):
                    queriables.append('ms2')
                
                # Find ambiguity
                if adduct_data['ambiguity'] != 'No':
                    ambiguity = True
                    
                temp_peaks_quality = []
                
                for index, rt in enumerate(adduct_data['peaks']):
                    queriables.append(str(rt))
                    temp_peaks_quality.append(check_peak_quality(
                                                                 ppm = adduct_data['ppm'][index],
                                                                 iso = adduct_data['iso'][index],
                                                                 curve = adduct_data['curve'][index],
                                                                 sn = adduct_data['sn'][index]
                                                                 )
                                              )
                temp_adducts_quality[adduct] = temp_peaks_quality
            
            glycan_qualities = [item for sublist in temp_adducts_quality.values() for item in sublist]
            if 0 in glycan_qualities:
                glycan_tag = 'good'
            elif 1 in glycan_qualities:
                glycan_tag = 'warning'
            else:
                glycan_tag = 'bad'
                
            queriables.append(glycan_tag)
            if ambiguity:
                queriables.append('ambiguity')
                glycan_value = "   ♦"
                
            queriables = ";".join(queriables).lower()
            
            query_not_found = False
            for query in filter_query:
                if query not in queriables:
                    query_not_found = True
                    break
            if query_not_found:
                continue
            
            if glycan_tag == 'good':
                good_count += 1
            elif glycan_tag == 'warning':
                average_count += 1
            else:
                bad_count += 1
                
            if glycan_value == "   ♦":
                if glycan not in noted_ambiguities:
                    ambiguity_count += 1
                noted_ambiguities.add(glycan)
                noted_ambiguities.update(
                                         k.split("_")[0]
                                         for k in glycans_per_sample[selected_item][glycan][list(glycans_per_sample[selected_item][glycan].keys())[0]]['ambiguity'].split(", ")
                                        )
            
            # Start adding the glycans, adducts and peaks to the list
            glycan_row = chromatograms_list.insert("", "end", text=glycan, value=(glycan_value) if ambiguity else '', tags=(glycan_tag, ))
            
            for adduct, adduct_data in temp_adducts_quality.items():
                # Check the adduct quality
                if 0 in adduct_data:
                    adduct_tag = 'good'
                elif 1 in adduct_data:
                    adduct_tag = 'warning'
                else:
                    adduct_tag = 'bad'
                
                # Insert the adduct row
                adduct_row = chromatograms_list.insert(glycan_row, "end", text=f"{adduct} - {glycans_per_sample[selected_item][glycan][adduct]['mz']}", tags=(adduct_tag, ))
                
                for index, peak_data in enumerate(adduct_data):
                    # Check peak quality
                    if peak_data == 0:
                        peak_tag = 'good'
                    elif peak_data == 1:
                        peak_tag = 'warning'
                    else:
                        peak_tag = 'bad'
                        
                    # Default peak value
                    peak_value = ''
                    
                    if 'ms2' in glycans_per_sample[selected_item][glycan][adduct].keys():
                        # Lower and upper peak boundaries
                        lower_boundary = curve_fittings[sample_index][f"{glycan}+{adduct}_{glycans_per_sample[selected_item][glycan][adduct]['peaks'][index]}_RTs"][0]
                        upper_boundary = curve_fittings[sample_index][f"{glycan}+{adduct}_{glycans_per_sample[selected_item][glycan][adduct]['peaks'][index]}_RTs"][-1]
                        
                        # Targets list
                        targets = np.array(list(glycans_per_sample[selected_item][glycan][adduct]['ms2'].keys()))
                        
                        # Filtered array
                        if len(targets[(targets >= lower_boundary) & (targets <= upper_boundary)]) > 0:
                            peak_value = 'MS2'
                        
                    # Insert peak row
                    peak_row = chromatograms_list.insert(adduct_row, "end", text=glycans_per_sample[selected_item][glycan][adduct]['peaks'][index], values=(peak_value), tags=(peak_tag, ))
                    
        chromatograms_qc_numbers.config(text=f"Compositions Quality:\n        Good: {good_count}    Average: {average_count}    Bad: {bad_count}\n        Ambiguities: {ambiguity_count}")
    else:        
        chromatograms_qc_numbers.config(text=f"Compositions Quality:\n        Good: {0}    Average: {0}    Bad: {0}\n        Ambiguities: {0}")
        
    # If it has the same length as before, go to the same position and expand the same rows
    if len(chromatograms_list.get_children()) == len_list:
        counter = 0
        for item_id in chromatograms_list.get_children():
            if counter > len(open_close_state)-1:
                break
            chromatograms_list.item(item_id, open=open_close_state[counter])
            counter+= 1
            for sub_item_id in chromatograms_list.get_children(item_id):
                if counter > len(open_close_state)-1:
                    break
                chromatograms_list.item(sub_item_id, open=open_close_state[counter])
                counter+= 1
        chromatograms_list.yview_moveto(view_top)
    
def clear_plot(ax_here, canvas_here):
    '''This function clears the specified plot and redraws it.'''
    ax_here.clear()
    canvas_here.draw()
    
def on_right_click_plot(event, ax_here, canvas_here, clean_plot):
    '''This function shows the menu containing the Save Image button, which cleans the plot of markers and allows to save a high resolution image of it.'''
    over_x = event.y > ax_here.bbox.y1 or event.y < ax_here.bbox.y0
    over_y = event.x > ax_here.bbox.x1 or event.x < ax_here.bbox.x0
    
    # Get the pointer position relative to the screen
    screen_x, screen_y = main_window.winfo_pointerxy()
    
    if over_x or over_y:
        return
        
    if not over_x and not over_y:  # Right mouse button
        # Open a context menu
        context_menu = tk.Menu(main_window, tearoff=0)
        context_menu.add_command(label="Save Image", command=lambda: save_image(event, canvas_here, ax_here, clean_plot))
        if len([line.get_label() for line in ax_here.get_lines() if "MS2" in line.get_label()]) > 0:
            if len([[line.get_label(), line.get_ydata()[0]] for line in ax_spec_ms2.lines if isinstance(line, matplotlib.lines.Line2D) and line.get_marker() == '*']) > 0:
                context_menu.add_command(label="Copy Annotations to Clipboard", command=save_annotations)
        context_menu.post(screen_x, screen_y)
        
def save_annotations():
    ''''''
    labels_intensity = [[line.get_label(), line.get_ydata()[0]] for line in ax_spec_ms2.lines if isinstance(line, matplotlib.lines.Line2D) and line.get_marker() == '*']
    
    main_window.clipboard_clear()
    main_window.clipboard_append(f"Fragment\tmz\tIntensity\n")
    for data in labels_intensity:
        label, intensity = data
        label, mz = label.split("\n")
        main_window.clipboard_append(f"{label}\t{mz}\t{intensity}\n")
        
        
def save_image(event, canvas_here, ax_here, clean_plot):
    '''This is the working function of the on_right_click_plot function. It saves an image of the right clicked plot.'''
    
    def remove_marker():
        removed_markers = []
        artists = ax_here.get_children()
        for artist in artists:
            if isinstance(artist, matplotlib.lines.Line2D) and (artist.get_markersize() == 4 or artist.get_linestyle() == '--'):
                removed_markers.append(artist)
                artist.remove()
        return removed_markers
    
    removed_markers = []
    if clean_plot:
        removed_markers = remove_marker()
    canvas_here.draw()
    
    # Open a save dialog window
    file_dialog = tk.Toplevel()
    file_dialog.attributes("-topmost", True)
    file_dialog.withdraw()
    file_dialog.grab_set()
    
    file_path = filedialog.asksaveasfilename(defaultextension=".svg", filetypes=[("SVG files", "*.svg"), ("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("All files", "*.*")])
    if file_path:
        if file_path.lower().endswith('.svg'):
            canvas_here.figure.savefig(file_path, format='svg', dpi=600)
        else:
            canvas_here.figure.savefig(file_path, dpi=600)
        
    for artist in removed_markers:
        ax_here.add_line(artist)
    canvas_here.draw()
        
    file_dialog.destroy()
    
def on_mouse_press_spectra_plot(event, ax_here):
    global clicked_spectra
    
    if clicked_spectra == None:
        clicked_spectra = event
    
def on_mouse_release_spectra_plot(event, ax_here):
    global clicked_spectra, current_line_ruler, distance_ruler_text
    
    clicked_spectra = None
    
    if current_line_ruler != None and distance_ruler_text != None:
        ruler_lines_list.append(current_line_ruler+[[distance_ruler_text]])
        current_line_ruler = None
        distance_ruler_text = None
        
def select_ruler(event):
    global ruler_lines_list, selected_ruler, spectrum_ruler_on
    
    if spectrum_ruler_on:
        for ruler_number, annotation in enumerate(ruler_lines_list):
            cont, _ = annotation[3][0].contains(event)
            if selected_ruler != None:
                selected_ruler[0][3][0].set_bbox(None)
                selected_ruler = None
            if cont:
                selected_ruler = [annotation, ruler_number]
                selected_ruler[0][3][0].set_bbox(dict(boxstyle='square', fc='none', ec='black', lw=1, linestyle=':'))
                break
    else:
        if selected_ruler != None:
            selected_ruler[0][3][0].set_bbox(None)
            selected_ruler = None
            
    canvas_spec.draw_idle()
    canvas_spec_ms2.draw_idle()
        
def on_pan(event, ax_here, canvas_here, type_coordinate):
    '''This function handles the panning of chromatogram and spectrum plots.'''
    global panning_enabled, current_line_ruler, distance_ruler_text, ruler_lines_list, gg_draw_on, gg_draw_list, gg_draw_selected
    
    over_x = event.y > ax_here.bbox.y1 or event.y < ax_here.bbox.y0
    over_y = event.x > ax_here.bbox.x1 or event.x < ax_here.bbox.x0
    
    if event.name == 'button_press_event':
        if (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            # Ignore panning when the shift key is pressed
            return
        panning_enabled = True
        on_pan.last_x, on_pan.last_y = event.x, event.y
    elif event.name == 'button_release_event':
        panning_enabled = False
        on_pan.last_x, on_pan.last_y = None, None
    elif event.name == 'motion_notify_event':
        if panning_enabled and event.button == 1:
            if on_pan.last_x is None or on_pan.last_y is None:
                return
            if over_x:
                dx = (event.x - on_pan.last_x) * 1  # Adjust the panning speed here
                x_min, x_max = ax_here.get_xlim()
                x_scale = (x_max - x_min) / ax_here.bbox.width
                ax_here.set_xlim(x_min - dx * x_scale, x_max - dx * x_scale)
                
                # Automatically adjust y-axis based on the most intense y-value in the visible region
                x_data = ax_here.get_lines()[0].get_xdata()
                y_data = ax_here.get_lines()[0].get_ydata()
                visible_indices = np.where((x_data >= x_min) & (x_data <= x_max))
                if len(visible_indices[0]) > 0:
                    if type_coordinate == 'chromatogram':
                        lines = ax_here.get_lines()
                        max_value = float('-inf')
                        for line in lines:
                            y_check_data = line.get_ydata()
                            try:
                                line_max_value = np.max(y_check_data[visible_indices])
                                max_value = max(max_value, line_max_value)
                            except:
                                pass
                        if max_value != float('-inf'):
                            max_y_value = max_value
                    else:
                        max_y_value = np.max(y_data[visible_indices])
                    ax_here.set_ylim(0, max_y_value + 0.1 * max_y_value)  # Adjust y-axis limit
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
            elif over_y:
                dy = (event.y - on_pan.last_y) * 1  # Adjust the panning speed here
                y_min, y_max = ax_here.get_ylim()
                y_scale = (y_max - y_min) / ax_here.bbox.height
                ax_here.set_ylim(y_min - dy * y_scale, y_max - dy * y_scale)
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
            else:
            
                if gg_draw_on:
                    if gg_draw_selected != None:
                        gg_draw_list[gg_draw_selected][1].xybox = (event.xdata, event.ydata)
                        gg_draw_list[gg_draw_selected][3] = [event.xdata, event.ydata]
                        return
                            
                if not spectrum_ruler_on or type_coordinate != 'spectra':
                    dx = (event.x - on_pan.last_x) * 1  # Adjust the panning speed here
                    x_min, x_max = ax_here.get_xlim()
                    x_scale = (x_max - x_min) / ax_here.bbox.width
                    ax_here.set_xlim(x_min - dx * x_scale, x_max - dx * x_scale)
                    dy = (event.y - on_pan.last_y) * 1  # Adjust the panning speed here
                    y_min, y_max = ax_here.get_ylim()
                    y_scale = (y_max - y_min) / ax_here.bbox.height
                    ax_here.set_ylim(y_min - dy * y_scale, y_max - dy * y_scale)
                    if type_coordinate == 'spectra':
                        annotate_top_y_values(ax_here, canvas_here)
                else:
                    if clicked_spectra != None:
                        if current_line_ruler == None:
                            # Horizontal line from click to current mouse position
                            horizontal_line = ax_here.plot([clicked_spectra.xdata, event.xdata], [event.ydata, event.ydata], color='blue', linestyle='-', zorder = 0)
                            
                            # Vertical lines of the ruler
                            vertical_line1 = ax_here.plot([clicked_spectra.xdata, clicked_spectra.xdata], [event.ydata, 0], color='blue', linestyle='-', zorder = 0)
                            vertical_line2 = ax_here.plot([event.xdata, event.xdata], [event.ydata, 0], color='blue', linestyle='-', zorder = 0)
                            
                            current_line_ruler = [horizontal_line, vertical_line1, vertical_line2]
                        else:
                            current_line_ruler[0][0].set_data([clicked_spectra.xdata, event.xdata], [event.ydata, event.ydata])
                            current_line_ruler[1][0].set_data([clicked_spectra.xdata, clicked_spectra.xdata], [event.ydata, 0])
                            current_line_ruler[2][0].set_data([event.xdata, event.xdata], [event.ydata, 0])
                        distance = abs(event.xdata - clicked_spectra.xdata)
                            
                        # Distance annotation of the ruler
                        x_line = current_line_ruler[0][0].get_xdata()
                        y_line = current_line_ruler[0][0].get_ydata()
                        
                        y_limits = ax_here.get_ylim()
                        if distance_ruler_text == None:
                            distance_ruler_text = ax_here.annotate(f'{distance:.4f}', xy=(x_line[0]+(x_line[1]-x_line[0])/2, y_line[1]+((y_limits[1]-y_limits[0])*0.02)), ha='center', fontsize=9, clip_on=True)
                        else:
                            distance_ruler_text.set_x(x_line[0]+(x_line[1]-x_line[0])/2)
                            distance_ruler_text.set_y(y_line[1]+((y_limits[1]-y_limits[0])*0.02))
                            distance_ruler_text.set_text(f'{distance:.4f}')
                    
            on_pan.last_x, on_pan.last_y = event.x, event.y
global panning_enabled
panning_enabled = False
on_pan.last_x, on_pan.last_y = None, None

def on_pan_right_click(event, ax_here, canvas_here, type_coordinate):
    '''This function handles the right-click panning of chromatograms and spectrum plots.'''
    global panning_enabled_right_click
    over_x = event.y > ax_here.bbox.y1 or event.y < ax_here.bbox.y0
    over_y = event.x > ax_here.bbox.x1 or event.x < ax_here.bbox.x0
    
    if event.name == 'button_press_event':
        panning_enabled_right_click = True
        on_pan_right_click.last_x, on_pan_right_click.last_y = event.x, event.y
    elif event.name == 'button_release_event':
        panning_enabled_right_click = False
        on_pan_right_click.last_x, on_pan_right_click.last_y = None, None
    elif event.name == 'motion_notify_event':
        if panning_enabled_right_click and event.button == 3:
            if on_pan_right_click.last_x is None or on_pan_right_click.last_y is None:
                return
            if over_x:
                dx = (event.x - on_pan_right_click.last_x) * 1.2  # Adjust the panning speed here
                x_min, x_max = ax_here.get_xlim()
                x_scale = (x_max - x_min) / ax_here.bbox.width
                ax_here.set_xlim(x_min, x_max - dx * x_scale)
                
                # Automatically adjust y-axis based on the most intense y-value in the visible region
                x_data = ax_here.get_lines()[0].get_xdata()
                y_data = ax_here.get_lines()[0].get_ydata()
                visible_indices = np.where((x_data >= x_min) & (x_data <= x_max))
                if len(visible_indices[0]) > 0:
                    if type_coordinate == 'chromatogram':
                        lines = ax_here.get_lines()
                        max_value = float('-inf')
                        for line in lines:
                            y_check_data = line.get_ydata()
                            try:
                                line_max_value = np.max(y_check_data[visible_indices])
                                max_value = max(max_value, line_max_value)
                            except:
                                pass
                        if max_value != float('-inf'):
                            max_y_value = max_value
                    else:
                        max_y_value = np.max(y_data[visible_indices])
                    ax_here.set_ylim(0, max_y_value + 0.1 * max_y_value)  # Adjust y-axis limit
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
            elif over_y:
                dy = (event.y - on_pan_right_click.last_y) * 1.2  # Adjust the panning speed here
                y_min, y_max = ax_here.get_ylim()
                y_scale = (y_max - y_min) / ax_here.bbox.height
                ax_here.set_ylim(y_min, y_max - dy * y_scale)
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
            on_pan_right_click.last_x, on_pan_right_click.last_y = event.x, event.y
global panning_enabled_right_click
panning_enabled_right_click = False
on_pan_right_click.last_x, on_pan_right_click.last_y = None, None
            
def stop_thread(t):
    '''This function attempts to finalize a function running on a non-mainthread, generated by threading package.'''
    if t.is_alive():
        t_id = ctypes.c_long(t.ident)
        res = ctypes.pythonapi.PyThreadState_SetAsyncExc(t_id, ctypes.py_object(SystemExit))
        if res == 0:
            raise ValueError("nonexistent thread id")
        elif res > 1:
            ctypes.pythonapi.PyThreadState_SetAsyncExc(t_id, 0)
            raise SystemError("PyThreadState_SetAsyncExc failed")
       
def annotate_top_y_values(ax_here, canvas_here):
    '''This function annotates the highest intensity peaks for their mz value on spectrums.'''
    # If no data is available in the axis object, ignore annotation
    if len(ax_here.get_lines()) == 0:
        return
        
    # Removes former labels
    clean_up_labels(ax_here)
    marker_coordinates = {}
    for line in ax_here.get_lines():
        if line.get_marker() == '*':
            marker_coordinates[float("%.1f" % round(line.get_data()[0][0], 1))] = line.get_label()

    # Get current limits of the plot
    x_min, x_max_here = ax_here.get_xlim()
    y_min, y_max_here = ax_here.get_ylim()
        
    # Get data within the current limits
    line = ax_here.get_lines()[0]
    x_data = line.get_xdata()
    y_data = line.get_ydata()
    mask = (x_data >= x_min) & (x_data <= x_max_here) & (y_data >= y_min) & (y_data <= y_max_here)
    visible_y_values = y_data[mask]
    visible_x_values = x_data[mask]

    if len(visible_y_values) == 0:
        return

    # Determine the number of datapoints
    min_max_here_annotations = (5, 15)
    
    number_points = int(len(visible_y_values)*0.05)
    if number_points < min_max_here_annotations[0]:
        number_points = min_max_here_annotations[0]
    elif number_points > min_max_here_annotations[1]:
        number_points = min_max_here_annotations[1]
    
    sorted_indices = np.argsort(visible_y_values)[::-1]
    
    selected_points_y = []
    selected_points_x = []
    
    annotations_density = (x_max_here-x_min)/25
    
    for i in sorted_indices:
        if len(selected_points_y) == number_points:
            break
        good = True
        for k in selected_points_x:
            if abs(visible_x_values[i]-k) < annotations_density:
                good = False
                break
        if good:
            selected_points_x.append(visible_x_values[i])
            selected_points_y.append(visible_y_values[i])
    
    # Annotate the top 5% intensities with mz
    for x_value, y_value in zip(selected_points_x, selected_points_y):
        if float("%.1f" %round(x_value, 1)) not in marker_coordinates:
            ax_here.annotate(f'{x_value:.4f}', xy=(x_value, y_value), xytext=(0, 10), textcoords='offset points', ha='center', fontsize=8, clip_on=True)
        else:
            ax_here.annotate(marker_coordinates[float("%.1f" %round(x_value, 1))], xy=(x_value, y_value), xytext=(0, 10), textcoords='offset points', ha='center', fontsize=8, clip_on=True)
    
def clean_up_labels(ax_here):
    '''Removes all annotations from plots.'''
    # Remove all text annotations
    numbers = '0123456789'
    for annotation in ax_here.texts:
        if annotation.get_fontsize() == 9:
            continue
        annotation.remove()
        
def on_key_press(event):
    """ Mark a key as pressed """
    key_states[event.keysym] = True

def on_key_release(event):
    """ Mark a key as released """
    key_states[event.keysym] = False
        
def sample_grouping_window(samples_names, parent):
    
    def on_drag_start(event):
        global dragged_widget
        
        dragged_widget = event.widget
        x = sgw.winfo_pointerx()-sgw.winfo_rootx()
        y = sgw.winfo_pointery()-sgw.winfo_rooty()
        
        widget_text = dragged_widget.tooltip.text
        dragged_widget.destroy()
            
        if len(widget_text) > 30:
            dragged_widget = ttk.Label(sgw, text=widget_text[:30]+"...")
        else:
            dragged_widget = ttk.Label(sgw, text=widget_text)
        ToolTip(dragged_widget, widget_text)
        dragged_widget.place(x=x, y=y)
    
    def on_drag(event):
        global dragged_widget
        
        if not dragged_widget:
            return
            
        x = sgw.winfo_pointerx()-sgw.winfo_rootx()
        y = sgw.winfo_pointery()-sgw.winfo_rooty()
        dragged_widget.place(x=x, y=y)
    
    def on_drop(event):
        global dragged_widget
        
        if not dragged_widget:
            return
        
        x, y = sgw.winfo_pointerx(), sgw.winfo_pointery()
        
        samples_bbox = (samples_frame.winfo_rootx(), samples_frame.winfo_rooty(), samples_frame.winfo_rootx()+samples_frame.winfo_width(), samples_frame.winfo_rooty()+samples_frame.winfo_height())
        
        # Check if it's dragged to any group
        found = False
        for group_frame in [child for child in sgw_groups_scrollable_frame.winfo_children() if isinstance(child, ttk.Frame)]:
            group_widgets = group_frame.winfo_children()
            for widget in group_widgets:
                if isinstance(widget, tk.Canvas):
                    group_samples_list_canvas = widget
                    group_samples_list = widget.winfo_children()[0]
            
            group_bbox = (group_frame.winfo_rootx(), group_frame.winfo_rooty(), group_frame.winfo_rootx()+group_frame.winfo_width(), group_frame.winfo_rooty()+group_frame.winfo_height())
            
            if is_inside_group(x, y, group_bbox):
                found = True
                widget_text = dragged_widget.tooltip.text
                dragged_widget.destroy()
                if len(widget_text) > 30:
                    label = ttk.Label(group_samples_list, text=widget_text[:30]+"...")
                else:
                    label = ttk.Label(group_samples_list, text=widget_text)
                ToolTip(label, widget_text)
                
                counter = 0
                while True:
                    free = True
                    for widget in group_samples_list.winfo_children():
                        grid_info = widget.grid_info()
                        row = grid_info.get('row', None)
                        if row == counter:
                            free = False
                            break
                    if free:
                        break
                    counter+=1
                    
                label.grid(row=counter, column=0, sticky='w')
                label.bind("<ButtonPress-1>", on_drag_start)
        
                if platform.system() == 'Windows':
                    label.bind("<MouseWheel>", lambda event: scroll_canvas(event, group_samples_list_canvas))
                else:
                    label.bind("<Button-4>", lambda event: scroll_canvas(event, group_samples_list_canvas))    # Linux Scroll Up
                    label.bind("<Button-5>", lambda event: scroll_canvas(event, group_samples_list_canvas))
                
                break
                
        if not found:
            widget_text = dragged_widget.tooltip.text
            dragged_widget.destroy()
        
            label = ttk.Label(sgw_samples_scrollable_frame, text=widget_text)
            ToolTip(label, widget_text)
            
            counter = 0
            while True:
                free = True
                for widget in sgw_samples_scrollable_frame.winfo_children():
                    grid_info = widget.grid_info()
                    row = grid_info.get('row', None)
                    if row == counter:
                        free = False
                        break
                if free:
                    break
                counter+=1
                
            label.grid(row=counter, column=0, columnspan=3, sticky='w')
            label.bind("<ButtonPress-1>", on_drag_start)
        
            if platform.system() == 'Windows':
                label.bind("<MouseWheel>", lambda event: scroll_canvas(event, sgw_samples_canvas))
            else:
                label.bind("<Button-4>", lambda event: scroll_canvas(event, sgw_samples_canvas))    # Linux Scroll Up
                label.bind("<Button-5>", lambda event: scroll_canvas(event, sgw_samples_canvas))
            
        dragged_widget = None

    def is_inside_group(x, y, group_bbox):
        """Check if the coordinates are inside the group's bounding box."""
        return group_bbox[0] <= x <= group_bbox[2] and group_bbox[1] <= y <= group_bbox[3]
    
    def add_group(parent, group_name=''):  
        def destroy_frame():
            for widget in group_scrollable_frame.winfo_children():
                counter = 0
                while True:
                    free = True
                    for widgetb in sgw_samples_scrollable_frame.winfo_children():
                        grid_info = widgetb.grid_info()
                        row = grid_info.get('row', None)
                        if row == counter:
                            free = False
                            break
                    if free:
                        break
                    counter+=1
                add_draggable(sgw_samples_scrollable_frame, widget.tooltip.text, counter, sgw_samples_canvas)
            group_frame.destroy()
            
        group_frame = ttk.Frame(parent, relief='solid')
        group_frame.pack(side=tk.LEFT, fill="y", expand=True)
    
        group_frame.grid_columnconfigure(1, weight=1)
        
        group_label = ttk.Label(group_frame, text="Group: ", font=("Segoe UI", list_font_size))
        group_label.grid(row=0, column=0, sticky='nsew', padx=5, pady=5)
        
        group_entry = ttk.Entry(group_frame, font=("Segoe UI", list_font_size))
        group_entry.grid(row=0, column=1, sticky='nsew', padx=(0, 5), pady=5)
        group_entry.insert(0, group_name)
        
        remove_button = tk.Button(group_frame, text=" X ", relief="flat", font=("Segoe UI", list_font_size), command=destroy_frame)
        remove_button.grid(row=0, column=2, sticky='nwse', padx=(0, 5), pady=5)
        
        group_canvas = tk.Canvas(group_frame, borderwidth=0, width=0, height=200)
        group_scrollbar = tk.Scrollbar(group_frame, orient="vertical", command=group_canvas.yview)
        group_scrollable_frame = tk.Frame(group_canvas)
        group_canvas.create_window((0, 0), window=group_scrollable_frame, anchor="nw")
        group_canvas.configure(yscrollcommand=group_scrollbar.set)
        group_scrollable_frame.bind("<Configure>", lambda event: group_canvas.configure(scrollregion=group_canvas.bbox("all")))
        group_canvas.grid(row=1, column=0, columnspan=3, sticky='nwse', padx=(5, 15), pady=5)
        group_scrollbar.grid(row=1, column=0, columnspan=3, sticky='nse', padx=(0, 5), pady=5)
        
        if platform.system() == 'Windows':
            group_canvas.bind("<MouseWheel>", lambda event: scroll_canvas(event, group_canvas))
        else:
            group_canvas.bind("<Button-4>", lambda event: scroll_canvas(event, group_canvas))    # Linux Scroll Up
            group_canvas.bind("<Button-5>", lambda event: scroll_canvas(event, group_canvas))
            
        return group_canvas, group_scrollable_frame
        
    def add_draggable(parent, text, row, canvas, trim = False):
        label_text = text[:30]+"..." if (trim and len(text) > 30) else text
        label = ttk.Label(parent, text=label_text)
        label.grid(row=row, column=0, sticky="w")
        ToolTip(label, text)
        
        # Bind drag events to the label
        label.bind("<ButtonPress-1>", on_drag_start)
        
        if platform.system() == 'Windows':
            label.bind("<MouseWheel>", lambda event: scroll_canvas(event, canvas))
        else:
            label.bind("<Button-4>", lambda event: scroll_canvas(event, canvas))    # Linux Scroll Up
            label.bind("<Button-5>", lambda event: scroll_canvas(event, canvas))
    
    def scroll_canvas(event, event_widget):
        """Handles mouse wheel scrolling based on OS."""
        if event.delta:  # Windows & macOS
            event_widget.yview_scroll(-1 * (event.delta // 120), "units")
        else:  # Linux (event.num is used instead of event.delta)
            if event.num == 4:  # Scroll up
                event_widget.yview_scroll(-1, "units")
            elif event.num == 5:  # Scroll down
                event_widget.yview_scroll(1, "units")
                
    def groupings_to_ds(data_type = 'dict'):
        temp_grouping = {}
        
        # Add grouped samples
        for group_frame in [child for child in sgw_groups_scrollable_frame.winfo_children() if isinstance(child, ttk.Frame)]:
            group_widgets = group_frame.winfo_children()
            for widget in group_widgets:
                if isinstance(widget, ttk.Entry):
                    if len(widget.get()) == 0:
                        error_window("Fill in the groups names before pressing confirming.")
                        return
                    temp_grouping[widget.get()] = []
                    for group_widget in group_widgets:
                        if isinstance(group_widget, tk.Canvas):
                            group_samples_list = group_widget.winfo_children()[0].winfo_children()
                            for sample in group_samples_list:
                                temp_grouping[widget.get()].append(sample.tooltip.text)
            
        if data_type == 'dict':
            return temp_grouping
            
        else:
            group_list = []
            sample_list = []
            for group, samples in temp_grouping.items():
                for sample in samples:
                    group_list.append(group)
                    sample_list.append(sample)
                    
            return group_list, sample_list
            
    def close_sgw():
        sgw.destroy()
        parent.lift()
        parent.grab_set()
                
    def sgw_ok():
        global sample_groups
        
        temp_sample_groups = groupings_to_ds()
        
        if temp_sample_groups:
            sample_groups = temp_sample_groups
            close_sgw()
            
        
    def save_to_file():
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        
        if len(file_path) != 0:
            try:
                with open(file_path, 'w') as file:
                    file.write(f"sample,group\n")
                    data = groupings_to_ds(data_type = 'other')
                    
                    for group, sample in zip(*data):
                        file.write(f"{sample},{group}\n")
                confirmation_window("File saved succesfully!")
            except Exception as e:
                error_window("The file you're trying to replace is already in use or you have no permission to write on that folder.")
            
    def load_from_file():
        global sample_groups
        file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv"), ("All files", "*.*")])
        
        if len(file_path) != 0:
            temp_grouping = {}
            try:
                with open(file_path, "r") as file:
                    for index, line in enumerate(file):
                        if index == 0:
                            continue
                        sample, group = line.strip().split(",")
                        if group in temp_grouping.keys():
                            temp_grouping[group].append(sample)
                        else:
                            temp_grouping[group] = [sample]
            except Exception as e:
                error_window("Something went wrong when loading the groups from file. Make sure it is in the correct format.")
                return
                        
            sample_groups = temp_grouping
            
            for group_frame in [child for child in sgw_groups_scrollable_frame.winfo_children() if isinstance(child, ttk.Frame)]:
                group_frame.destroy()
                
            for widget in sgw_samples_scrollable_frame.winfo_children():
                widget.destroy()
                        
            load_existing_grouping()
                        
    def load_existing_grouping():
        added_samples = []
        for group in sample_groups.keys():
            group_canvas, group_samples_list = add_group(sgw_groups_scrollable_frame, group)
            index = 0
            for sample in sample_groups[group]:
                if sample in samples_names:
                    added_samples.append(sample)
                    add_draggable(group_samples_list, sample, index, group_canvas, True)
                    index+=1
                
        index = 0
        for name in samples_names:
            if name not in added_samples:
                add_draggable(sgw_samples_scrollable_frame, name, index, sgw_samples_canvas)
                index += 1
        
    sgw = tk.Toplevel()
    sgw.geometry("510x560")
    icon = ImageTk.PhotoImage(ico_image)
    sgw.iconphoto(False, icon)
    sgw.protocol("WM_DELETE_WINDOW", close_sgw)
    sgw.grab_set()
    
    global dragged_widget
    dragged_widget = None
    
    sgw.bind("<B1-Motion>", on_drag)
    sgw.bind("<ButtonRelease-1>", on_drop)
    
    sgw.grid_rowconfigure(0, weight=0)
    sgw.grid_rowconfigure(1, weight=1)
    sgw.grid_columnconfigure(0, weight=1)
    
    groups_frame = ttk.Labelframe(sgw, text="Groups")
    groups_frame.grid(row=0, column=0, padx=10, pady=(10,0), sticky="nsew")
    
    group_add_button = ttk.Button(groups_frame, text="Add group", style="small_button_sfw_style1.TButton", command=lambda: add_group(sgw_groups_scrollable_frame))
    group_add_button.pack(anchor='w')
    ToolTip(group_add_button, "Adds a new group where you can drag your samples to. Make sure to name it.")
    
    # Make the scrollable groups window
    sgw_groups_canvas = tk.Canvas(groups_frame, borderwidth=0)
    sgw_groups_scrollbar = tk.Scrollbar(groups_frame, orient="horizontal", command=sgw_groups_canvas.xview)
    sgw_groups_scrollable_frame = tk.Frame(sgw_groups_canvas)
    sgw_groups_canvas.create_window((0, 0), window=sgw_groups_scrollable_frame, anchor="nw", height=260)
    sgw_groups_canvas.configure(xscrollcommand=sgw_groups_scrollbar.set)
    sgw_groups_scrollable_frame.bind("<Configure>", lambda event: sgw_groups_canvas.configure(scrollregion=sgw_groups_canvas.bbox("all")))
    sgw_groups_canvas.pack(fill=tk.BOTH, expand=True)
    sgw_groups_scrollbar.pack(side="bottom", fill="x")
    
    samples_frame = ttk.Labelframe(sgw, text="Samples")
    samples_frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
    
    # Make the scrollable samples window
    sgw_samples_canvas = tk.Canvas(samples_frame, borderwidth=0)
    sgw_samples_scrollbar = tk.Scrollbar(samples_frame, orient="vertical", command=sgw_samples_canvas.yview)
    sgw_samples_scrollable_frame = tk.Frame(sgw_samples_canvas)
    sgw_samples_canvas.create_window((0, 0), window=sgw_samples_scrollable_frame, anchor="nw")
    sgw_samples_canvas.configure(yscrollcommand=sgw_samples_scrollbar.set)
    sgw_samples_scrollable_frame.bind("<Configure>", lambda event: sgw_samples_canvas.configure(scrollregion=sgw_samples_canvas.bbox("all")))
    sgw_samples_scrollbar.pack(side="right", fill="y")
    sgw_samples_canvas.pack(fill=tk.BOTH, expand=True)
        
    if platform.system() == 'Windows':
        sgw_samples_canvas.bind("<MouseWheel>", lambda event: scroll_canvas(event, sgw_samples_canvas))
    else:
        sgw_samples_canvas.bind("<Button-4>", lambda event: scroll_canvas(event, sgw_samples_canvas))    # Linux Scroll Up
        sgw_samples_canvas.bind("<Button-5>", lambda event: scroll_canvas(event, sgw_samples_canvas))
    
    sgw_ok_button = ttk.Button(sgw, text="Ok", style="small_button_sfw_style1.TButton", command=sgw_ok)
    sgw_ok_button.grid(row=2, column=0, padx=(0, 100), pady=10, sticky="se")
    
    sgw_cancel_button = ttk.Button(sgw, text="Cancel", style="small_button_sfw_style1.TButton", command=close_sgw)
    sgw_cancel_button.grid(row=2, column=0, padx=10, pady=10, sticky="se")
    
    sgw_save_button = ttk.Button(sgw, text="Save to File", style="small_button_sfw_style1.TButton", command=save_to_file)
    sgw_save_button.grid(row=2, column=0, padx=10, pady=10, sticky="sw")
    ToolTip(sgw_save_button, "Saves the grouping scheme to a CSV file. It includes a header (sample,group) and each sample (file name without extension) is disposed in a line below, with its group beside it. Example for samples (sample1, sample2, sample3, sample4) in groups (group1, group2):\n    sample,group\n    sample1,group1\n    sample2,group1\n    sample3,group2\n    sample4,group2")
    
    sgw_load_button = ttk.Button(sgw, text="Load from File", style="small_button_sfw_style1.TButton", command=load_from_file)
    sgw_load_button.grid(row=2, column=0, padx=(100, 0), pady=10, sticky="sw")
    ToolTip(sgw_load_button, "Loads the grouping scheme from a CSV file. It must include a header (sample,group) and each sample (file name without extension) disposed in a line below, with its group beside it. Example for samples (sample1, sample2, sample3, sample4) in groups (group1, group2):\n    sample,group\n    sample1,group1\n    sample2,group1\n    sample3,group2\n    sample4,group2")
    
    # If no grouping was set before, just create samples in the samples frame
    load_existing_grouping()

def run_main_window():
    '''This function runs the main window loop.'''
    # Functions to be used by main_window
    
    def file_name_window(file_type = None):
        '''This function creates a window for setting the file name that will be saved when the library generation, analysis or save results finish.'''
            
        invalid_characters = [':', '\"', '/', '\\', '|', '?', '*']
        
        def fn_ok(file_type = None):
            global exp_lib_name, gg_file_name
            for i in invalid_characters:
                if i in file_name_entry.get():
                    error_window("You can use ':', '\"', '/', '\\', '|', '?', '*' in your file name")
                    return
            if file_type == 'library':
                exp_lib_name = file_name_entry.get()
                run_generate_library()
            elif file_type == 'gg':
                gg_file_name = file_name_entry.get()
                run_analysis()
            file_name_window.destroy()
        
        def fn_cancel():
            global exp_lib_name, gg_file_name
            file_name_window.destroy()
        
        # Checks to see if it's ok to proceed with analysis
        if file_type == 'gg':
            if len(samples_list) == 0:
                error_window("You must select files for analysis on 'Select Files' menu!")
                gg_file_name = ""
                return
            elif save_path == "":
                error_window("You must select a working directory in the 'Set Parameters' menu before starting an analysis!")
                gg_file_name = ""
                return
            elif library_path == "":
                error_window("You must generate or import a library before starting an analysis!")
                gg_file_name = ""
                return
            
            uncalibrated_files = False
            for file in samples_list:
                if "calibrated" not in file:
                    uncalibrated_files = True
                    break
            if uncalibrated_files:
                proceed = warning_window("Before you start an analysis, make sure your spectra files are calibrated. Calibration done after the analysis won't influence the analysis results! You can calibrate your files by clicking on the 'File Editor' icon, under the 'About' button, in the main window. If you are sure your files are calibrated enough, proceed.")
                if not proceed:
                    return
                
        elif file_type == 'library':
            if save_path == "":
                error_window("You must select a working directory in the 'Set Parameters' window before generating a library!")
                return
            if not ion_mode_selected:
                error_window("You must select an Ion Mode in the Set Parameters window before generating a library.")
                return
        
        file_name_window = tk.Toplevel()
        file_name_window.withdraw()
        file_name_window.title("File Name")
        icon = ImageTk.PhotoImage(ico_image)
        file_name_window.iconphoto(False, icon)
        file_name_window.resizable(False, False)
        file_name_window.grab_set()
        file_name_window.protocol("WM_DELETE_WINDOW", on_closing)
        
        file_name_label = ttk.Label(file_name_window, text='Type in the name for the file you want to save:', font=("Segoe UI", list_font_size))
        file_name_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        
        file_name_entry = ttk.Entry(file_name_window, width=50)
        if file_type == 'library':
            file_name_entry.insert(0, '<date>_<time>_glycans_lib')
        elif file_type == 'gg':
            file_name_entry.insert(0, '<date>_<time>_Analysis')
        file_name_entry.grid(row=1, column=0, padx=10, pady=(0, 10), sticky='w')
        
        file_name_ok_button = ttk.Button(file_name_window, text="Ok", style="small_button_spw_style1.TButton", command=lambda: fn_ok(file_type))
        file_name_ok_button.grid(row=2, column=0, padx=(10, 95), pady=(0,10), sticky="e")
        
        file_name_cancel_button = ttk.Button(file_name_window, text="Cancel", style="small_button_spw_style1.TButton", command=fn_cancel)
        file_name_cancel_button.grid(row=2, column=0, padx=(10, 10), pady=(0,10), sticky="e")
                
        file_name_window.update_idletasks()
        file_name_window.deiconify()
        window_width = file_name_window.winfo_width()
        window_height = file_name_window.winfo_height()
        screen_width = file_name_window.winfo_screenwidth()
        screen_height = file_name_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        file_name_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
    def open_select_files_window():
        run_select_files_window(samples_dropdown)
        
    def open_set_parameters_window():
        run_set_parameters_window()
            
    def open_file_dialog_import_button():
        global imp_exp_library
        file_dialog = tk.Toplevel()
        file_dialog.attributes("-topmost", True)
        file_dialog.withdraw()
        file_dialog.grab_set()
        
        file_path = filedialog.askopenfilename(filetypes=[("Glygogenius Library files", "*.ggl"), ("All files", "*.*")])
        if len(file_path) != 0:
            global library_path
            if len(library_path) > 0:
                generate_library_button_frame.config(bg=background_color)
            library_path = file_path
            import_library_info.config(state=tk.NORMAL)
            imp_exp_library = [True, False]
            import_library_button_frame.config(bg="lightgreen")
        else:
            library_path = file_path
            import_library_info.config(state=tk.DISABLED)
            imp_exp_library = [False, False]
            import_library_button_frame.config(bg=background_color)
        
        generate_library_button_frame.config(bg=background_color)
        file_dialog.destroy()
        
    def get_lib_info():
        def export_extra_files(full_library, library_metadata):
            global save_path
            
            extra_files_lactonized_ethyl_esterified = library_metadata.get('lactonized/ethyl esterified', False)
            extra_files_custom_monos = library_metadata.get('custom monosaccharides', [])
                
            extra_files_exp_lib_name = library_path.split("/")[-1].split(".")[0]
            
            if len(save_path) == 0:
                extra_files_save_path = "/".join(library_path.split("/")[:-1])
            else:
                extra_files_save_path = save_path
            
            Execution_Functions.output_extra_library_files(full_library, extra_files_lactonized_ethyl_esterified, extra_files_custom_monos, [False, False], library_path, extra_files_exp_lib_name, library_metadata, extra_files_save_path)
            
            confirmation_window("Extra files saved succesfully!")
            
        def library_browser(full_library, library_metadata):
            def autofit_columns(tree):
                for col in tree["columns"]:
                    max_width = max(len(tree.heading(col)["text"]), *[len(tree.set(item, col)) for item in tree.get_children()])
                    tree.column(col, width=max_width * 10)  # Approximate character width
                    
            def search_and_highlight(event):
                global highlighted_items, current_index
                query = event.widget.get().lower()
                
                for item in highlighted_items:
                    glycans_treeview.item(item, tags=())
                    
                highlighted_items = []
                
                if query:
                    for item in glycans_treeview.get_children():
                        values = glycans_treeview.item(item, "values")
                        if any(query in str(value).lower() for value in values):
                            highlighted_items.append(item)
                            glycans_treeview.item(item, tags=('highlight',))
                            
                current_index = 0 if highlighted_items else -1
                glycans_treeview.tag_configure('highlight', background='#f1f5b5')
                
                if highlighted_items:
                    glycans_treeview.selection_set(highlighted_items[current_index])
                    glycans_treeview.see(highlighted_items[current_index])
                    
            def navigate(direction):
                global current_index
                if not highlighted_items:
                    return
                
                current_index = (current_index + direction) % len(highlighted_items)
                glycans_treeview.selection_set(highlighted_items[current_index])
                glycans_treeview.see(highlighted_items[current_index])
        
            def glycans_treeview_sort(tv, col, reverse):
                # Get the data in the specified column
                data = [(tv.set(k, col), k) for k in tv.get_children('')]
                
                # Try to convert data to float if possible for numerical sorting
                try:
                    data = [(float(item[0]), item[1]) for item in data]
                except ValueError:
                    pass
                
                # Sort data
                data.sort(reverse=reverse)
                
                # Rearrange items in sorted positions
                for index, (val, k) in enumerate(data):
                    tv.move(k, '', index)
                
                # Reverse sort next time
                tv.heading(col, command=lambda: glycans_treeview_sort(tv, col, not reverse))
                
            def close_lib_browser():
                lib_browser.destroy()
                lib_info_window.grab_set()
                    
            global highlighted_items        
            highlighted_items = []
            
            lib_browser = tk.Toplevel()
            # lib_browser.withdraw()
            lib_browser.title("Library Browser")
            icon = ImageTk.PhotoImage(ico_image)
            lib_browser.iconphoto(False, icon)
            lib_browser.minsize(400, 700)
            lib_browser.protocol("WM_DELETE_WINDOW", close_lib_browser)
            lib_browser.grab_set()
            
            lib_browser_lactonized_ethyl_esterified = library_metadata.get('lactonized/ethyl esterified', False)
            lib_browser_custom_monos = library_metadata.get('custom monosaccharides', [])
            
            glycans_treeview_columns = ['Glycan', 'Hex', 'HexN', 'HexNAc', 'Xylose', 'dHex']
            if lib_browser_lactonized_ethyl_esterified:
                glycans_treeview_columns = glycans_treeview_columns+['a2,3-Neu5Ac', 'a2,6-Neu5Ac', 'a2,3-Neu5Gc', 'a2,6-Neu5Gc', 'UroA']
            else:
                glycans_treeview_columns = glycans_treeview_columns+['Neu5Ac', 'Neu5Gc', 'UroA']
                
            if len(lib_browser_custom_monos) > 0:
                for cm in lib_browser_custom_monos:
                    name = cm['cm_name']

                    if name in glycans_treeview_columns:
                        name+= '-custom'
                    
                    glycans_treeview_columns.append(name)
                    
            glycans_treeview_columns = glycans_treeview_columns+['Isotopic Distribution', 'Neutral Mass + Tag']
            for adduct in full_library[list(full_library.keys())[0]]['Adducts_mz']:
                glycans_treeview_columns.append(adduct)
            
            lib_browser.grid_columnconfigure(1, weight=1)
            lib_browser.grid_rowconfigure(1, weight=1)
    
            symbol_button_style1 = ttk.Style().configure("symbol_button_style1.TButton", font=("Segoe UI", list_font_size), relief="raised", padding = (0, 0), justify="center", width=5)
            
            search_field_label = ttk.Label(lib_browser, text="Search:", font=("Segoe UI", list_font_size))
            search_field_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
            
            search_field_entry = ttk.Entry(lib_browser)
            search_field_entry.grid(row=0, column=1, padx=(0, 10), pady=10, sticky="we")
            search_field_entry.bind("<KeyRelease>", search_and_highlight)
            
            search_back = ttk.Button(lib_browser, text="◄", style="symbol_button_style1.TButton", command=lambda: navigate(-1))
            search_back.grid(row=0, column=2, padx=(0, 10), pady=10, sticky="we")
            
            search_next = ttk.Button(lib_browser, text="►", style="symbol_button_style1.TButton", command=lambda: navigate(1))
            search_next.grid(row=0, column=3, padx=(0, 10), pady=10, sticky="we")
            
            glycans_treeview_scrollbar = tk.Scrollbar(lib_browser, orient=tk.VERTICAL)
            glycans_treeview = ttk.Treeview(lib_browser, columns = glycans_treeview_columns, yscrollcommand=glycans_treeview_scrollbar.set, show='headings')
        
            for col in glycans_treeview_columns:
                glycans_treeview.heading(col, text=col, command=lambda _col=col: glycans_treeview_sort(glycans_treeview, _col, False))
                
            glycans_treeview_scrollbar.config(command=glycans_treeview.yview, width=10)
            
            for glycan in full_library:
                row_values = [glycan, full_library[glycan]['Monos_Composition']['H'], full_library[glycan]['Monos_Composition']['HN'], full_library[glycan]['Monos_Composition']['N'], full_library[glycan]['Monos_Composition']['X'], full_library[glycan]['Monos_Composition']['F']]
                if lib_browser_lactonized_ethyl_esterified:
                    row_values = row_values+[full_library[glycan]['Monos_Composition']['Am'], full_library[glycan]['Monos_Composition']['E'], full_library[glycan]['Monos_Composition']['AmG'], full_library[glycan]['Monos_Composition']['EG'], full_library[glycan]['Monos_Composition']['UA']]
                else:
                    row_values = row_values+[full_library[glycan]['Monos_Composition']['S'], full_library[glycan]['Monos_Composition']['G'], full_library[glycan]['Monos_Composition']['UA']]
                if len(lib_browser_custom_monos) > 0:
                    for cm in lib_browser_custom_monos:
                        row_values.append(full_library[glycan]['Monos_Composition'].get(cm['cm_short_code'], 0))
                temp_isotopic = []
                for j in full_library[glycan]['Isotopic_Distribution']:
                    temp_isotopic.append(float("%.3f" % round(j, 3)))
                row_values.append(str(temp_isotopic)[1:-1])
                row_values.append(float("%.4f" % round(full_library[glycan]['Neutral_Mass+Tag'], 4)))
                for adduct in full_library[glycan]['Adducts_mz']:
                    row_values.append(float("%.4f" % round(full_library[glycan]['Adducts_mz'][adduct], 4)))
                    
                glycans_treeview.insert("", "end", values=row_values)
                
            glycans_treeview.grid(row=1, column=0, columnspan=4, padx=(10, 20), pady=(0, 10), sticky="nsew")
            glycans_treeview_scrollbar.grid(row=1, column=0, columnspan=4, padx=10, pady=10, sticky="nse")
            
            autofit_columns(glycans_treeview)
        
            library_length_label = ttk.Label(lib_browser, text=f"Library Length: {len(full_library)} glycans", font=("Segoe UI", list_font_size))
            library_length_label.grid(row=2, column=0, columnspan=4, padx=10, pady=10, sticky="w")
        
            close_lib_browser_button = ttk.Button(lib_browser, text="Close", style="small_button_sfw_style1.TButton", command=close_lib_browser)
            close_lib_browser_button.grid(row=2, column=0, columnspan=4, padx=10, pady=10, sticky="e")
            
        try:
            with open(library_path, 'rb') as f:
                library_data = dill.load(f)
                f.close()
        except:
            error_window("Can't load the library. The library file you are trying to load was created in a version older than 0.2.0 or is corrupted. Generate a new library or try a different file.")
            return
            
        full_library = library_data[0]
        library_metadata = library_data[1]
        
        if type(library_metadata) != dict:
            error_window("Can't load the library. The library file you are trying to load was created in a version older than 1.3.0 or is corrupted. Generate a new library or try a different file.")
            return
        
        lib_info_window = tk.Toplevel()
        # lib_info_window.attributes("-topmost", True)
        lib_info_window.withdraw()
        lib_info_window.title("Library Information")
        icon = ImageTk.PhotoImage(ico_image)
        lib_info_window.iconphoto(False, icon)
        lib_info_window.resizable(False, False)
        lib_info_window.grab_set()

        information_text = ScrolledText(lib_info_window, width=52, height=18, wrap=tk.WORD)
        information_text.grid(row=0, column=0, padx = 10, pady = 10, sticky="new")
        
        if library_metadata.get('custom glycans list', [False, []])[0]:
            information_text.insert(tk.END, f"Custom glycans list: {library_metadata.get('custom glycans list', [False, []])[1]}\n\n")
        else:
            information_text.insert(tk.END, f"Min/Max number of monosaccharides: {library_metadata.get('min/max monosaccharides', [0, 0])[0]}/{library_metadata.get('min/max monosaccharides', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of Hexoses: {library_metadata.get('min/max hexoses', [0, 0])[0]}/{library_metadata.get('min/max hexoses', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of Hexosamines: {library_metadata.get('min/max hexosamines', [0, 0])[0]}/{library_metadata.get('min/max hexosamines', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of N-Acetylhexosamines: {library_metadata.get('min/max hexnac', [0, 0])[0]}/{library_metadata.get('min/max hexnac', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of Xyloses: {library_metadata.get('min/max xyloses', [0, 0])[0]}/{library_metadata.get('min/max xyloses', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of deoxyHexoses: {library_metadata.get('min/max fucoses', [0, 0])[0]}/{library_metadata.get('min/max fucoses', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of Sialic Acids: {library_metadata.get('min/max sialic acids', [0, 0])[0]}/{library_metadata.get('min/max sialic acids', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of N-Acetylneuraminic Acids: {library_metadata.get('min/max acetyl sialic acids', [0, 0])[0]}/{library_metadata.get('min/max acetyl sialic acids', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of N-Glycolylneuraminic Acids: {library_metadata.get('min/max glycolyl sialic acids', [0, 0])[0]}/{library_metadata.get('min/max glycolyl sialic acids', [0, 0])[1]}\n")
            information_text.insert(tk.END, f"Min/Max number of Uronic Acids: {library_metadata.get('min/max uronic acids', [0, 0])[0]}/{library_metadata.get('min/max uronic acids', [0, 0])[1]}\n\n")
            
        if len(library_metadata.get('custom monosaccharides', [])) > 0:
            information_text.insert(tk.END, f"Custom monosaccharides:\n")
            for cm in library_metadata.get('custom monosaccharides', []):
                information_text.insert(tk.END, f"    - {cm['cm_name']} ({cm['cm_short_code']}): {cm['cm_chem_comp']}\n")
            information_text.insert(tk.END, f"\n")
            
        information_text.insert(tk.END, f"Glycan class: {library_metadata.get('glycan class', None)}\n")
        
        adduct_combos_dict = {General_Functions.comp_to_formula(adduct[0]): {'charges': adduct[1], 'comp': adduct[0]} for adduct in library_metadata.get('adduct combos', [])}
        adduct_combos_string = ''
        for adduct, info in adduct_combos_dict.items():
            if len(adduct_combos_string) == 0:
                adduct_combos_string += f"{adduct}: {'+' if info['charges'] > 0 else ''}{info['charges']}"
            else:
                adduct_combos_string += f"\n         {adduct}: {'+' if info['charges'] > 0 else ''}{info['charges']}"
        
        information_text.insert(tk.END, f"Adducts: {adduct_combos_string}\n")
        
        information_text.insert(tk.END, f"Maximum charges: {library_metadata.get('maximum charges', 3)}\n")
        
        if library_metadata.get('reducing end tag', 0.0) in reducing_end_tags.values():
            reducing_end_tag_text = f" ({list(reducing_end_tags.keys())[list(reducing_end_tags.values()).index(library_metadata.get('reducing end tag', 0.0))]})"
        else:
            reducing_end_tag_text = ''
            
        information_text.insert(tk.END, f"Reducing end tag: {library_metadata.get('reducing end tag', 0.0)}{reducing_end_tag_text}\n")
        information_text.insert(tk.END, f"Internal Standard: {library_metadata.get('internal standard', '0.0')}\n")
        information_text.insert(tk.END, f"Permethylated: {library_metadata.get('permethylated', False)}\n")
        information_text.insert(tk.END, f"Amidated/Ethyl-Esterified: {library_metadata.get('lactonized/ethyl esterified', False)}\n")
        information_text.insert(tk.END, f"Reduced end: {library_metadata.get('reducing end reduced', False)}\n")
        information_text.insert(tk.END, f"Min/Max number of Sulfations: {library_metadata.get('min/max sulfation', [0, 0])[0]}/{library_metadata.get('min/max sulfation', [0, 0])[1]}\n")
        information_text.insert(tk.END, f"Min/Max number of Phosphorylations: {library_metadata.get('min/max phosphorylation', [0, 0])[0]}/{library_metadata.get('min/max phosphorylation', [0, 0])[1]}\n")
        information_text.insert(tk.END, f"Fast isotopic distribution calculation: {library_metadata.get('fast isotopic pattern calculation', True)}\n")
        information_text.insert(tk.END, f"High resolution isotopic distribution: {library_metadata.get('high resolution isotopic pattern', False)}")
            
        information_text.see("end")
        
        browse_lib_button = ttk.Button(lib_info_window, text="Browse Library", style="small_button_sfw_style1.TButton", command=lambda: library_browser(full_library, library_metadata))
        browse_lib_button.grid(row=1, column=0, padx=10, pady=10, sticky="w")
        ToolTip(browse_lib_button, "Allows you to check the glycans within the library and search for specific ones.")
        
        export_extra_files_button = ttk.Button(lib_info_window, text="Export Extra Files", style="small_button_sfw_style1.TButton", command=lambda: export_extra_files(full_library, library_metadata))
        export_extra_files_button.grid(row=1, column=0, padx=(105, 10), pady=10, sticky="w")
        ToolTip(export_extra_files_button, "This button will output a human-readable excel spreadsheet containing all the glycans within the library. Additionally, it will output a Skyline-compatible transitions list.")
        
        close_lib_info_button = ttk.Button(lib_info_window, text="Close", style="small_button_sfw_style1.TButton", command=lib_info_window.destroy)
        close_lib_info_button.grid(row=1, column=0, padx=10, pady=10, sticky="e")
        
        lib_info_window.update_idletasks()
        lib_info_window.deiconify()
        window_width = lib_info_window.winfo_width()
        window_height = lib_info_window.winfo_height()
        screen_width = lib_info_window.winfo_screenwidth()
        screen_height = lib_info_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        lib_info_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
    def run_generate_library():
        def close_progress_gen_lib_window():
            sys.stdout = original_stdout
            try:
                stop_thread(t)
            except:
                pass
            progress_gen_lib.destroy()
            
        def ok_progress_gen_lib_window():
            global library_path, imp_exp_library
            library_path = os.path.join(save_path, library_name+".ggl")
            imp_exp_library = [False, False]
            generate_library_button_frame.config(bg="lightgreen")
            import_library_button_frame.config(bg=background_color)
            import_library_info.config(state=tk.NORMAL)
            close_progress_gen_lib_window()
    
        progress_gen_lib = tk.Toplevel(main_window)
        #progress_gen_lib.attributes("-topmost", True)
        progress_gen_lib.withdraw()
        progress_gen_lib.title("Generating Library")
        icon = ImageTk.PhotoImage(ico_image)
        progress_gen_lib.iconphoto(False, icon)
        progress_gen_lib.resizable(False, False)
        progress_gen_lib.grab_set()
        progress_gen_lib.protocol("WM_DELETE_WINDOW", on_closing)

        output_text = ScrolledText(progress_gen_lib, width=50, height=10, wrap=tk.WORD)
        output_text.grid(row=0, column=0, columnspan=2, padx = 10, pady = 10, sticky="new")
        
        global ok_lib_gen_button
        ok_lib_gen_button = ttk.Button(progress_gen_lib, text="Ok", style="small_button_sfw_style1.TButton", command=ok_progress_gen_lib_window, state=tk.DISABLED)
        ok_lib_gen_button.grid(row=1, column=0, padx=10, pady=10, sticky="e")
        
        cancel_lib_gen_button = ttk.Button(progress_gen_lib, text="Cancel", style="small_button_sfw_style1.TButton", command=close_progress_gen_lib_window)
        cancel_lib_gen_button.grid(row=1, column=1, padx=10, pady=10, sticky="w")
        
        progress_gen_lib.update_idletasks()
        progress_gen_lib.deiconify()
        window_width = progress_gen_lib.winfo_width()
        window_height = progress_gen_lib.winfo_height()
        screen_width = progress_gen_lib.winfo_screenwidth()
        screen_height = progress_gen_lib.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        progress_gen_lib.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        original_stdout = sys.stdout
        sys.stdout = TextRedirector_Gen_Lib(output_text)
        
        global min_max_monos, min_max_hex, min_max_hexnac, min_max_xyl, min_max_fuc, min_max_sia, min_max_ac, min_max_ac, min_max_gc, min_max_hn, min_max_ua, forced, min_max_proton, max_charges, reducing_end_tag, internal_standard, permethylated, lactonized_ethyl_esterified, reduced, min_max_sulfation, min_max_phosphorylation, fast_iso, high_res, exp_lib_name, lyase_digested
        
        output_filtered_data_args = {
                                    'curve fitting score threshold': curve_fit_score,
                                    'isotopic fitting score threshold': iso_fit_score,
                                    'signal to noise ratio threshold': s_to_n,
                                    'maximum ppm threshold': max_ppm,
                                    'minimum percentage auc threshold': percentage_auc,
                                    'reanalysis': reanalysis,
                                    'reanalysis_path': reanalysis_path,
                                    'save path': save_path,
                                    'analyze ms2': analyze_ms2[0],
                                    'unrestricted fragments': analyze_ms2[2],
                                    'reporter ions': reporter_ions,
                                    'output abundance table': plot_metaboanalyst,
                                    'output composition-separated data': compositions,
                                    'align chromatograms': align_chromatograms,
                                    'glycan class': forced,
                                    'retention time tolerance': ret_time_interval[2],
                                    'fragmentation retention time tolerance': rt_tolerance_frag,
                                    'output isotopic fitting data': iso_fittings,
                                    'output chromatogram plotting data': output_plot_data,
                                    'multithreaded': multithreaded_analysis,
                                    'number of cpu cores': number_cores,
                                    'analysis done time': 0.0,
                                    'minimum percentage of samples': min_samples,
                                    'sample groups': None,
                                    'temporary folder': None,
                                    'fill data gaps': fill_gaps,
                                    'from GUI': True
                                    }

        imp_exp_gen_library_args = {
                                    'custom glycans list': custom_glycans_list,
                                    'min/max monosaccharides': min_max_monos,
                                    'min/max hexoses': min_max_hex,
                                    'min/max hexnac': min_max_hexnac,
                                    'min/max xyloses': min_max_xyl,
                                    'min/max sialic acids': min_max_sia,
                                    'min/max fucoses': min_max_fuc,
                                    'min/max acetyl sialic acids': min_max_ac,
                                    'min/max glycolyl sialic acids': min_max_gc,
                                    'min/max hexosamines': min_max_hn,
                                    'min/max uronic acids': min_max_ua,
                                    'glycan class': forced,
                                    'min/max protons': min_max_proton,
                                    'custom adducts': custom_adducts,
                                    'maximum charges': max_charges,
                                    'reducing end tag': reducing_end_tag,
                                    'fast isotopic pattern calculation': fast_iso,
                                    'high resolution isotopic pattern': high_res,
                                    'import/export library': [False, False],
                                    'library path': library_path,
                                    'exported library name': exp_lib_name,
                                    'only generate library': True,
                                    'save path': save_path,
                                    'internal standard': internal_standard,
                                    'permethylated': permethylated,
                                    'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                                    'reducing end reduced': reduced,
                                    'min/max sulfation': min_max_sulfation,
                                    'min/max phosphorylation': min_max_phosphorylation,
                                    'lyase digested': lyase_digested,
                                    'temporary folder': None,
                                    'custom monosaccharides': custom_monosaccharides,
                                    'from GUI': True
                                    }

        list_of_data_args = {
                            'samples list': samples_list
                            }

        index_spectra_from_file_ms1_args = {
                                            'raw data': None,
                                            'ms level': 1,
                                            'multithreaded': multithreaded_analysis,
                                            'number of cpu cores': number_cores
                                            }
                                            
        index_spectra_from_file_ms2_args = {
                                            'raw data': None,
                                            'ms level': 2,
                                            'multithreaded': multithreaded_analysis,
                                            'number of cpu cores': number_cores
                                            }

        analyze_files_args = {
                            'library': None,
                            'library size': None,
                            'raw data': None,
                            'ms1 index': None,
                            'mass tolerance': tolerance,
                            'retention time interval': ret_time_interval,
                            'min isotopologue peaks': min_isotopologue_peaks,
                            'min points per peak': min_ppp,
                            'adduct combos': None,
                            'maximum charges': max_charges,
                            'custom noise level': custom_noise,
                            'only x most intense peaks': close_peaks,
                            'multithreaded': multithreaded_analysis,
                            'number of cpu cores': number_cores,
                            'analysis start time': None,
                            'temporary folder': None,
                            'from GUI': True
                            }

        analyze_ms2_args = {
                            'ms2 index': None,
                            'raw data': None,
                            'analyzed data': None,
                            'retention time interval': ret_time_interval,
                            'mass tolerance': tolerance,
                            'min/max monosaccharides': min_max_monos,
                            'min/max hexoses': min_max_hex,
                            'min/max hexnac': min_max_hexnac,
                            'min/max xyloses': min_max_xyl,
                            'min/max sialic acids': min_max_sia,
                            'min/max fucoses': min_max_fuc,
                            'min/max acetyl sialic acids': min_max_ac,
                            'min/max glycolyl sialic acids': min_max_gc,
                            'min/max hexosamines': min_max_hn,
                            'min/max uronic acids': min_max_ua,
                            'adduct combos': None,
                            'maximum charges': max_charges,
                            'reducing end tag': reducing_end_tag,
                            'glycan class': forced,
                            'permethylated': permethylated,
                            'reducing end reduced': reduced,
                            'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                            'assign fragments compatible with precursor': analyze_ms2[1],
                            'annotate all glycans in library': analyze_ms2[2],
                            'fragmentation retention time tolerance': ret_time_interval[2],
                            'multithreaded': multithreaded_analysis,
                            'number of cpu cores': number_cores,
                            'library': None,
                            'temporary folder': None,
                            'from GUI': True,
                            'custom monosaccharides': custom_monosaccharides
                            }

        arrange_raw_data_args = {
                                'analyzed data': None,
                                'sample names': samples_names,
                                'ms2 analyzed': analyze_ms2[0],
                                'save path': save_path,
                                'analysis parameters': None,
                                'adduct combos': None,
                                'library': None,
                                'temporary folder': None,
                                'results file name': None,
                                'from GUI': True,
                                'erase files': True,
                                'calculated noise levels': []
                                }
        
        t = threading.Thread(target=run_glycogenius, args=([{
            'output_filtered_data_args': output_filtered_data_args, 
            'imp_exp_gen_library_args': imp_exp_gen_library_args, 
            'list_of_data_args': list_of_data_args, 
            'index_spectra_from_file_ms1_args': index_spectra_from_file_ms1_args,
            'index_spectra_from_file_ms2_args': index_spectra_from_file_ms2_args, 
            'analyze_files_args': analyze_files_args, 
            'analyze_ms2_args': analyze_ms2_args, 
            'arrange_raw_data_args':arrange_raw_data_args
            }]))
        t.start()
            
    def run_analysis(): 
        global gg_file_name
    
        def loading_files_after_analysis():
            def close_lf():
                loading_files.destroy()
                
            def wait_thread():
                reanalysis_file.join()
                if reanalysis_file == 'Error':
                    loading_files.destroy()
                    return
                samples_dropdown.current(0)
                handle_selection(None)
                global ms1_bound, ms1_binds, ms2_bound, ms2_binds, chromatogram_bound, chromatogram_binds
                clear_plot(ax, canvas)
                clear_plot(ax_spec, canvas_spec)
                clear_plot(ax_spec_ms2, canvas_spec_ms2)
                
                # Unbind MS1 plot
                if ms1_bound:
                    for i in ms1_binds:
                        canvas_spec.mpl_disconnect(i)
                    ms1_bound = False
                    ms1_binds = []
                    
                # Unbind MS2 plot
                if ms2_bound:
                    for i in ms2_binds:
                        canvas_spec_ms2.mpl_disconnect(i)
                    ms2_bound = False
                    ms2_binds = []
                    
                # Unbind chromatogram plot
                if chromatogram_bound:
                    for i in chromatogram_binds:
                        canvas.mpl_disconnect(i)
                    chromatogram_bound = False
                    chromatogram_binds = []
                    
                close_lf()
                
            loading_files = tk.Toplevel()
            # loading_files.attributes("-topmost", True)
            loading_files.withdraw()
            loading_files.title("Loading Files")
            icon = ImageTk.PhotoImage(ico_image)
            loading_files.iconphoto(False, icon)
            loading_files.resizable(False, False)
            loading_files.grab_set()
            loading_files.protocol("WM_DELETE_WINDOW", on_closing)
            
            loading_files_label = ttk.Label(loading_files, text="Loading files, please wait.", font=("Segoe UI", list_font_size))
            loading_files_label.pack(pady=35, padx=70)
            
            loading_files.update_idletasks()
            loading_files.deiconify()
            window_width = loading_files.winfo_width()
            window_height = loading_files.winfo_height()
            screen_width = loading_files.winfo_screenwidth()
            screen_height = loading_files.winfo_screenheight()
            x_position = (screen_width - window_width) // 2
            y_position = (screen_height - window_height) // 2
            loading_files.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
            
            t = threading.Thread(target=wait_thread)
            t.start()
            
        def close_progress_run_analysis_window():
            run_analysis.destroy()
            sys.stdout = original_stdout
            try:
                stop_thread(t)
                kill_concurrent_futures()
            except:
                pass
                
        def ok_progress_run_analysis_window():
            run_analysis_button_frame.config(bg="lightgreen")
            if len(gg_file_name) != 0:
                global reanalysis_path
                reanalysis_path = os.path.join(save_path, gg_file_name)
            close_progress_run_analysis_window()
            global reanalysis_file, former_selected_sample
            reanalysis_file = threading.Thread(target=load_reanalysis, args=(reanalysis_path,))
            reanalysis_file.start()
            loading_files_after_analysis()
            check_qc_dist_button.config(state=tk.NORMAL)
            plot_graph_button.config(state=tk.NORMAL)
            
            # Clear history
            former_selected_sample = ''
        
        def get_parameters_analysis():
            def close_analysis_param():
                run_analysis.grab_set()
                analysis_info_window.destroy()
            try:
                with open(library_path, 'rb') as f:
                    library_data = dill.load(f)
                    f.close()
            except:
                error_window("Can't load the library. The library file you are trying to load was created in a version older than 0.2.0 or is corrupted. Generate a new library or try a different file.")
                return
            full_library = library_data[0]
            library_metadata = library_data[1]
            
            analysis_info_window = tk.Toplevel()
            analysis_info_window.withdraw()
            analysis_info_window.title("Analysis Parameters")
            icon = ImageTk.PhotoImage(ico_image)
            analysis_info_window.iconphoto(False, icon)
            analysis_info_window.resizable(False, False)
            analysis_info_window.grab_set()

            information_text = ScrolledText(analysis_info_window, width=55, height=30, wrap=tk.WORD)
            information_text.grid(row=0, column=0, padx = 10, pady = 10, sticky="new")
            
            information_text.insert(tk.END, "Samples analyzed:\n")
            for i in samples_list:
                i_temp_line = i.split("/")[-1]
                information_text.insert(tk.END, f"- {i_temp_line}\n")
            information_text.insert(tk.END, "\n")
            information_text.insert(tk.END, "Library properties:\n")
            
            custom_glycans = library_metadata.get('custom glycans', [False, []])
            min_max_monos = library_metadata.get('min/max monosaccharides', [0, 0])
            min_max_hex = library_metadata.get('min/max hexoses', [0, 0])
            min_max_hexnac = library_metadata.get('min/max hexnac', [0, 0])
            min_max_fuc = library_metadata.get('min/max fucoses', [0, 0])
            min_max_sia = library_metadata.get('min/max sialic acids', [0, 0])
            min_max_ac = library_metadata.get('min/max acetyl sialic acids', [0, 0])
            min_max_gc = library_metadata.get('min/max glycolyl sialic acids', [0, 0])
            min_max_xyl = library_metadata.get('min/max xyloses', [0, 0])
            min_max_hn = library_metadata.get('min/max hexosamines', [0, 0])
            min_max_ua = library_metadata.get('min/max uronic acids', [0, 0])
            forced = library_metadata.get('glycan class', None)
            min_max_proton = library_metadata.get('min/max protons', [1, 3])
            custom_adducts = library_metadata.get('custom adducts', [])
            adduct_combos = library_metadata.get('adduct combos', [])
            max_charges = library_metadata.get('maximum charges', 3)
            reducing_end_tag = library_metadata.get('reducing end tag', 0.0)
            internal_standard = library_metadata.get('internal standard', '0.0')
            permethylated = library_metadata.get('permethylated', False)
            lactonized_ethyl_esterified = library_metadata.get('lactonized/ethyl esterified', False)
            reduced = library_metadata.get('reducing end reduced', False)
            fast_iso = library_metadata.get('fast isotopic pattern calculation', True)
            high_res = library_metadata.get('high resolution isotopic pattern', False)
            min_max_sulfation = library_metadata.get('min/max sulfation', [0, 0])
            min_max_phosphorylation = library_metadata.get('min/max phosphorylation', [0, 0])
            lyase_digested = library_metadata.get('lyase digested', False)
            custom_monosaccharides = library_metadata.get('custom monosaccharides', [])
            
            if custom_glycans[0]:
                library_type = f" - Custom glycans: {str(custom_glycans[1])[1:-1]}"
            else:
                library_type = f" - Monosaccharides: {str(min_max_monos)[1:-1]}\n - Hexoses: {str(min_max_hex)[1:-1]}\n - Hexosamines: {str(min_max_hn)[1:-1]}\n - HexNAcs: {str(min_max_hexnac)[1:-1]}\n - Xyloses: {str(min_max_xyl)[1:-1]}\n - Sialic Acids: {str(min_max_sia)[1:-1]}\n - Uronic Acids: {str(min_max_ua)[1:-1]}\n - dHex: {str(min_max_fuc)[1:-1]}\n - Neu5Acs: {str(min_max_ac)[1:-1]}\n - Neu5Gcs: {str(min_max_gc)[1:-1]}"
            
            if reducing_end_tag in reducing_end_tags.values():
                reducing_end_tag_text = f" ({list(reducing_end_tags.keys())[list(reducing_end_tags.values()).index(reducing_end_tag)]})"
            else:
                reducing_end_tag_text = ''
            
            adduct_combos_dict = {General_Functions.comp_to_formula(adduct[0]): {'charges': adduct[1], 'comp': adduct[0]} for adduct in adduct_combos}
            adduct_combos_string = ''
            for adduct, info in adduct_combos_dict.items():
                if len(adduct_combos_string) == 0:
                    adduct_combos_string += f"{adduct}: {'+' if info['charges'] > 0 else ''}{info['charges']}"
                else:
                    adduct_combos_string += f"\n            {adduct}: {'+' if info['charges'] > 0 else ''}{info['charges']}"
                    
            additional_info = f" - Glycan Class: {forced}\n - Adducts: {adduct_combos_string}\n - Maximum charges: {max_charges}\n - Reducing end tag: {reducing_end_tag}{reducing_end_tag_text}\n - Permethylated: {permethylated}\n - Reduced end: {reduced}\n - Amidated/Ethyl-Esterified: {lactonized_ethyl_esterified}\n - Min/Max number of Sulfations: {str(min_max_sulfation)[1:-1]}\n - Min/Max number of Phosphorylations: {str(min_max_phosphorylation)[1:-1]}\n - Fast isotopic calculations: {fast_iso}\n - High resolution isotopic calculations: {high_res}"
                
            information_text.insert(tk.END, f"{library_type}\n")
            information_text.insert(tk.END, "\n")
            
            if len(custom_monosaccharides) > 0:
                information_text.insert(tk.END, f"Custom monosaccharides:\n")
                for cm in custom_monosaccharides:
                    information_text.insert(tk.END, f"    - {cm['cm_name']} ({cm['cm_short_code']}): {cm['cm_chem_comp']} Min/Max: {cm['cm_min']}/{cm['cm_max']}\n")
                information_text.insert(tk.END, f"\n")
                
            information_text.insert(tk.END, f"{additional_info}\n")
            if forced == 'gags':
                information_text.insert(tk.END, f"- Lyase digested: {lyase_digested}\n")
            information_text.insert(tk.END, "\n")
            information_text.insert(tk.END, "Analysis settings:\n")
            information_text.insert(tk.END, f" - Analyze MS2: {analyze_ms2[0]}\n")
            if analyze_ms2[0]:
                information_text.insert(tk.END, f" -- Limit fragments assignment to composition: {analyze_ms2[1]}\n")
                information_text.insert(tk.END, f" -- Assign MS2 of glycans not found in MS1: {analyze_ms2[2]}\n")
            information_text.insert(tk.END, f" - Tolerance unit: {tolerance[0]}, tolerance value: {tolerance[1]}\n")
            information_text.insert(tk.END, f" - Retention/Migration time interval analyzed: {ret_time_interval[0]}, {ret_time_interval[1]}\n")
            information_text.insert(tk.END, f" - Custom minimum datapoints per peaks: {min_ppp[0]}\n")
            if min_ppp[0]:
                information_text.insert(tk.END, f" -- Minimum number of datapoints per peaks: {min_ppp[1]}\n")
            information_text.insert(tk.END, f" - Limit peaks picked per chromatogram/electropherogram: {close_peaks[0]}\n")
            if close_peaks[0]:
                information_text.insert(tk.END, f" -- Number of peaks: {close_peaks[1]}\n")
            
            information_text.see("end")
        
            close_analysis_info_button = ttk.Button(analysis_info_window, text="Close", style="small_button_sfw_style1.TButton", command=close_analysis_param)
            close_analysis_info_button.grid(row=1, column=0, padx=10, pady=10)
            
            analysis_info_window.update_idletasks()
            analysis_info_window.deiconify()
            window_width = analysis_info_window.winfo_width()
            window_height = analysis_info_window.winfo_height()
            screen_width = analysis_info_window.winfo_screenwidth()
            screen_height = analysis_info_window.winfo_screenheight()
            x_position = (screen_width - window_width) // 2
            y_position = (screen_height - window_height) // 2
            analysis_info_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
        global run_analysis
        
        try:
            with open(library_path, 'rb') as f:
                library_data = dill.load(f)
                f.close()
        except:
            error_window("Can't load the library. The library file you are trying to load was created in a version older than 0.2.0 or is corrupted. Generate a new library or try a different file.")
            return
            
        full_library = library_data[0]
        library_metadata = library_data[1]
        
        if type(library_metadata) != dict:
            error_window("Can't load the library. The library file you are trying to load was created in a version older than 1.3.0 or is corrupted. Generate a new library or try a different file.")
            return
        
        run_analysis = tk.Toplevel()
        # run_analysis.attributes("-topmost", True)
        run_analysis.withdraw()
        run_analysis.title("Running Analysis")
        icon = ImageTk.PhotoImage(ico_image)
        run_analysis.iconphoto(False, icon)
        run_analysis.resizable(False, False)
        run_analysis.grab_set()
        run_analysis.protocol("WM_DELETE_WINDOW", on_closing)

        output_text = ScrolledText(run_analysis, width=50, height=10, wrap=tk.WORD)
        output_text.grid(row=0, column=0, columnspan=2, padx = 10, pady = 10, sticky="new")
        
        global progress_bar_run_analysis, progress_bar_ms1_label, progress_bar_run_analysis2, progress_bar_ms2_label
        progress_bar_ms1_label = ttk.Label(run_analysis, text="Analysis Progress: 0%", font=("Segoe UI", 10))
        progress_bar_ms1_label.grid(row=1, column=0, columnspan=2, padx=10, sticky="w")
        progress_bar_run_analysis = ttk.Progressbar(run_analysis, orient="horizontal", mode="determinate")
        progress_bar_run_analysis.grid(row=2, column=0, columnspan=2, padx=10, sticky="ew")
        if analyze_ms2[0]:
            progress_bar_ms1_label.config(text="MS1 Analysis Progress: 0%")
            progress_bar_ms2_label = ttk.Label(run_analysis, text="MS2 Analysis Progress: 0%", font=("Segoe UI", 10))
            progress_bar_ms2_label.grid(row=3, column=0, columnspan=2, padx=10, sticky="w")
            progress_bar_run_analysis2 = ttk.Progressbar(run_analysis, orient="horizontal", mode="determinate")
            progress_bar_run_analysis2.grid(row=4, column=0, columnspan=2, padx=10, sticky="ew")
        
        global ok_run_analysis_button
        ok_run_analysis_button = ttk.Button(run_analysis, text="Ok", style="small_button_sfw_style1.TButton", command=ok_progress_run_analysis_window, state=tk.DISABLED)
        ok_run_analysis_button.grid(row=5, column=0, columnspan=2, padx=10, pady=10, sticky="w")
        
        cancel_run_analysis_button = ttk.Button(run_analysis, text="Cancel", style="small_button_sfw_style1.TButton", command=close_progress_run_analysis_window)
        cancel_run_analysis_button.grid(row=5, column=0, columnspan=2, padx=(100, 10), pady=10, sticky="w")
        
        check_param_analysis_button = ttk.Button(run_analysis, text="Check Analysis Parameters", style="small_button_sfw_style1.TButton", command=get_parameters_analysis)
        check_param_analysis_button.grid(row=5, column=1, padx=10, pady=10, sticky="e")
        
        run_analysis.update_idletasks()
        run_analysis.deiconify()
        window_width = run_analysis.winfo_width()
        window_height = run_analysis.winfo_height()
        screen_width = run_analysis.winfo_screenwidth()
        screen_height = run_analysis.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        run_analysis.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        original_stdout = sys.stdout
        sys.stdout = TextRedirector_Run_Analysis(output_text)
        
        global min_max_monos, min_max_hex, min_max_hexnac, min_max_xyl, min_max_fuc, min_max_sia, min_max_ac, min_max_ac, min_max_gc, min_max_hn, min_max_ua, min_max_gc, forced, min_max_proton, custom_adducts, max_charges, reducing_end_tag, internal_standard, permethylated, lactonized_ethyl_esterified, min_max_sulfation, min_max_phosphorylation, reduced, fast_iso, high_res, lyase_digested, custom_monosaccharides
        
        min_max_monos = library_metadata.get('min/max monosaccharides', [0, 0])
        min_max_hex = library_metadata.get('min/max hexoses', [0, 0])
        min_max_hexnac = library_metadata.get('min/max hexnac', [0, 0])
        min_max_fuc = library_metadata.get('min/max fucoses', [0, 0])
        min_max_sia = library_metadata.get('min/max sialic acids', [0, 0])
        min_max_ac = library_metadata.get('min/max acetyl sialic acids', [0, 0])
        min_max_gc = library_metadata.get('min/max glycolyl sialic acids', [0, 0])
        min_max_xyl = library_metadata.get('min/max xyloses', [0, 0])
        min_max_hn = library_metadata.get('min/max hexosamines', [0, 0])
        min_max_ua = library_metadata.get('min/max uronic acids', [0, 0])
        forced = library_metadata.get('glycan class', None)
        min_max_proton = library_metadata.get('min/max protons', [1, 3])
        custom_adducts = library_metadata.get('custom adducts', [])
        max_charges = library_metadata.get('maximum charges', 3)
        reducing_end_tag = library_metadata.get('reducing end tag', 0.0)
        internal_standard = library_metadata.get('internal standard', '0.0')
        permethylated = library_metadata.get('permethylated', False)
        lactonized_ethyl_esterified = library_metadata.get('lactonized/ethyl esterified', False)
        reduced = library_metadata.get('reducing end reduced', False)
        fast_iso = library_metadata.get('fast isotopic pattern calculation', True)
        high_res = library_metadata.get('high resolution isotopic pattern', False)
        min_max_sulfation = library_metadata.get('min/max sulfation', [0, 0])
        min_max_phosphorylation = library_metadata.get('min/max phosphorylation', [0, 0])
        lyase_digested = library_metadata.get('lyase digested', False)
        custom_monosaccharides = library_metadata.get('custom monosaccharides', [])

        analysis_parameters = {
                    'library settings': {
                        'min/max monosaccharides': min_max_monos,
                        'min/max hexoses': min_max_hex,
                        'min/max hexnac': min_max_hexnac,
                        'min/max fucoses': min_max_fuc,
                        'min/max sialic acids': min_max_sia,
                        'min/max acetyl sialic acids': min_max_ac,
                        'min/max glycolyl sialic acids': min_max_gc,
                        'min/max xyloses': min_max_xyl,
                        'min/max hexosamines': min_max_hn,
                        'min/max uronic acids': min_max_ua,
                        'glycan class': forced,
                        'min/max protons': min_max_proton,
                        'custom adducts': custom_adducts,
                        'maximum charges': max_charges,
                        'import/export library': imp_exp_library,
                        'library path': library_path,
                        'exported library name': exp_lib_name,
                        'only generate library': only_gen_lib,
                        'reducing end tag': reducing_end_tag,
                        'internal standard': internal_standard,
                        'permethylated': permethylated,
                        'lactonized/ethyl esterified': lactonized_ethyl_esterified,
                        'reducing end reduced': reduced,
                        'fast isotopic pattern calculation': fast_iso,
                        'high resolution isotopic pattern': high_res,
                        'custom glycans list': custom_glycans_list,
                        'min/max sulfation': min_max_sulfation,
                        'min/max phosphorylation': min_max_phosphorylation,
                        'lyase digested': lyase_digested,
                        'custom monosaccharides': custom_monosaccharides
                        },
                    'analysis settings': {
                        'multithreaded': multithreaded_analysis,
                        'number of cpu cores': number_cores,
                        'analyze ms2': analyze_ms2,
                        'reporter ions': reporter_ions,
                        'mass tolerance': tolerance,
                        'retention time interval': ret_time_interval,
                        'fragmentation retention time tolerance': rt_tolerance_frag,
                        'min isotopologue peaks': min_isotopologue_peaks,
                        'min points per peak': min_ppp,
                        'only x most intense peaks': close_peaks,
                        'align chromatograms': align_chromatograms,
                        'minimum percentage auc threshold': percentage_auc,
                        'maximum ppm threshold': max_ppm,
                        'isotopic fitting score threshold': iso_fit_score,
                        'curve fitting score threshold': curve_fit_score,
                        'signal to noise ratio threshold': s_to_n,
                        'custom noise level': custom_noise,
                        'samples path': samples_path,
                        'save path': save_path,
                        'output abundance table': plot_metaboanalyst,
                        'output composition-separated data': compositions,
                        'output isotopic fitting data': iso_fittings,
                        'output chromatogram plotting data': output_plot_data,
                        'reanalysis': reanalysis,
                        'reanalysis_path': reanalysis_path,
                        'fill data gaps': fill_gaps
                        }
                    }
            
        output_filtered_data_args = {
                                    'curve fitting score threshold': curve_fit_score,
                                    'isotopic fitting score threshold': iso_fit_score,
                                    'signal to noise ratio threshold': s_to_n,
                                    'maximum ppm threshold': max_ppm,
                                    'minimum percentage auc threshold': percentage_auc,
                                    'reanalysis': reanalysis,
                                    'reanalysis_path': reanalysis_path,
                                    'save path': save_path,
                                    'analyze ms2': analyze_ms2[0],
                                    'unrestricted fragments': analyze_ms2[2],
                                    'reporter ions': reporter_ions,
                                    'output abundance table': plot_metaboanalyst,
                                    'output composition-separated data': compositions,
                                    'align chromatograms': align_chromatograms,
                                    'glycan class': forced,
                                    'retention time tolerance': ret_time_interval[2],
                                    'fragmentation retention time tolerance': rt_tolerance_frag,
                                    'output isotopic fitting data': iso_fittings,
                                    'output chromatogram plotting data': output_plot_data,
                                    'multithreaded': multithreaded_analysis,
                                    'number of cpu cores': number_cores,
                                    'analysis done time': 0.0,
                                    'minimum percentage of samples': min_samples,
                                    'sample groups': None,
                                    'temporary folder': None,
                                    'fill data gaps': fill_gaps,
                                    'from GUI': True
                                    }

        imp_exp_gen_library_args = {
                                    'custom glycans list': custom_glycans_list,
                                    'min/max monosaccharides': min_max_monos,
                                    'min/max hexoses': min_max_hex,
                                    'min/max hexnac': min_max_hexnac,
                                    'min/max xyloses': min_max_xyl,
                                    'min/max sialic acids': min_max_sia,
                                    'min/max fucoses': min_max_fuc,
                                    'min/max acetyl sialic acids': min_max_ac,
                                    'min/max glycolyl sialic acids': min_max_gc,
                                    'min/max hexosamines': min_max_hn,
                                    'min/max uronic acids': min_max_ua,
                                    'glycan class': forced,
                                    'min/max protons': min_max_proton,
                                    'custom adducts': custom_adducts,
                                    'maximum charges': max_charges,
                                    'reducing end tag': reducing_end_tag,
                                    'fast isotopic pattern calculation': fast_iso,
                                    'high resolution isotopic pattern': high_res,
                                    'import/export library': [True, True],
                                    'library path': library_path,
                                    'exported library name': '',
                                    'only generate library': False,
                                    'save path': save_path,
                                    'internal standard': internal_standard,
                                    'permethylated': permethylated,
                                    'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                                    'reducing end reduced': reduced,
                                    'min/max sulfation': min_max_sulfation,
                                    'min/max phosphorylation': min_max_phosphorylation,
                                    'lyase digested': lyase_digested,
                                    'temporary folder': None,
                                    'custom monosaccharides': custom_monosaccharides,
                                    'from GUI': True
                                    }

        list_of_data_args = {
                            'samples list': samples_list
                            }

        index_spectra_from_file_ms1_args = {
                                            'raw data': None,
                                            'ms level': 1,
                                            'multithreaded': multithreaded_analysis,
                                            'number of cpu cores': number_cores
                                            }
                                            
        index_spectra_from_file_ms2_args = {
                                            'raw data': None,
                                            'ms level': 2,
                                            'multithreaded': multithreaded_analysis,
                                            'number of cpu cores': number_cores
                                            }

        analyze_files_args = {
                            'library': None,
                            'library size': None,
                            'raw data': None,
                            'ms1 index': None,
                            'mass tolerance': tolerance,
                            'retention time interval': ret_time_interval,
                            'min isotopologue peaks': min_isotopologue_peaks,
                            'min points per peak': min_ppp,
                            'adduct combos': None,
                            'maximum charges': max_charges,
                            'custom noise level': custom_noise,
                            'only x most intense peaks': close_peaks,
                            'multithreaded': multithreaded_analysis,
                            'number of cpu cores': number_cores,
                            'analysis start time': None,
                            'temporary folder': None,
                            'from GUI': True
                            }

        analyze_ms2_args = {
                            'ms2 index': None,
                            'raw data': None,
                            'analyzed data': None,
                            'retention time interval': ret_time_interval,
                            'mass tolerance': tolerance,
                            'min/max monosaccharides': min_max_monos,
                            'min/max hexoses': min_max_hex,
                            'min/max hexnac': min_max_hexnac,
                            'min/max xyloses': min_max_xyl,
                            'min/max sialic acids': min_max_sia,
                            'min/max fucoses': min_max_fuc,
                            'min/max acetyl sialic acids': min_max_ac,
                            'min/max glycolyl sialic acids': min_max_gc,
                            'min/max hexosamines': min_max_hn,
                            'min/max uronic acids': min_max_ua,
                            'adduct combos': None,
                            'maximum charges': max_charges,
                            'reducing end tag': reducing_end_tag,
                            'glycan class': forced,
                            'permethylated': permethylated,
                            'reducing end reduced': reduced,
                            'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                            'assign fragments compatible with precursor': analyze_ms2[1],
                            'annotate all glycans in library': analyze_ms2[2],
                            'fragmentation retention time tolerance': ret_time_interval[2],
                            'multithreaded': multithreaded_analysis,
                            'number of cpu cores': number_cores,
                            'library': None,
                            'temporary folder': None,
                            'from GUI': True,
                            'custom monosaccharides': custom_monosaccharides
                            }

        arrange_raw_data_args = {
                                'analyzed data': None,
                                'sample names': samples_names,
                                'ms2 analyzed': analyze_ms2[0],
                                'save path': save_path,
                                'analysis parameters': analysis_parameters,
                                'adduct combos': None,
                                'library': None,
                                'temporary folder': None,
                                'results file name': gg_file_name,
                                'from GUI': True,
                                'erase files': True,
                                'calculated noise levels': []
                                }
        
        t = threading.Thread(target=run_glycogenius, args=([{
            'output_filtered_data_args': output_filtered_data_args, 
            'imp_exp_gen_library_args': imp_exp_gen_library_args, 
            'list_of_data_args': list_of_data_args, 
            'index_spectra_from_file_ms1_args': index_spectra_from_file_ms1_args,
            'index_spectra_from_file_ms2_args': index_spectra_from_file_ms2_args, 
            'analyze_files_args': analyze_files_args, 
            'analyze_ms2_args': analyze_ms2_args, 
            'arrange_raw_data_args':arrange_raw_data_args
            }]))
        t.start()
            
    def save_results_button_command():
        if reanalysis_path == "":
            error_window("No reanalysis file found! Choose an reanalysis file in the 'Select Files' menu or run an analysis on sample files before trying to save results!")
        elif save_path == "":
            error_window("You must select a working directory in the 'Set Parameters' window before saving results!")
        else:
            save_results_window()
            
    def about_button_command():
        def close_about_window():
            about_window.destroy()
            
        about_window = tk.Toplevel(main_window)
        # about_window.attributes("-topmost", True)
        about_window.withdraw()
        about_window.title("About")
        icon = ImageTk.PhotoImage(ico_image)
        about_window.iconphoto(False, icon)
        about_window.resizable(False, False)
        about_window.grab_set()
        
        banner_label_about = tk.Label(about_window, image=tk_banner)
        banner_label_about.grid(row=0, column=0, sticky='nsew')
        
        about_text = tk.Label(about_window, text=f"GlycoGenius: Glycomics Data Analysis Tool\nCopyright (C) 2023 by Hector Franco Loponte\n\nGlycoGenius version: {gg_version}    GUI version: {GUI_version}\n\nThis program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.\n\nThis program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.\n\nYou should have received a copy of the GNU General Public License along with this program. It is accessible within the program files or by typing 'license' after running it stand-alone in the terminal by typing 'glycogenius'. If not, see <https://www.gnu.org/licenses/>.", justify="center", wraplength=700)
        about_text.grid(row=1, column=0, pady=(20, 20), sticky='nsew')
        
        close_about_button = ttk.Button(about_window, text="Close", style="small_button_sfw_style1.TButton", command=close_about_window)
        close_about_button.grid(row=2, column=0, pady=(20, 20))
        
        about_window.update_idletasks()
        about_window.deiconify()
        window_width = about_window.winfo_width()
        window_height = about_window.winfo_height()
        screen_width = about_window.winfo_screenwidth()
        screen_height = about_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        about_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
    def on_resize(event):
        new_width = event.width
        new_height = event.height
        
    def handle_treeview_select(event, clear = True):
        global ax, canvas, ax_spec, canvas_spec, colors, color_number, level, two_d, compare_samples_button, samples_dropdown_options, last_xlims_chrom, last_ylims_chrom, plot_graph_button, selected_chromatograms, last_plotted
        
        compare_samples_button.config(state=tk.DISABLED)
        plot_graph_button.config(state=tk.DISABLED)
        
        if len(chromatograms_list.get_children()) == 0:
            return
        
        if len(selected_chromatograms) == 1 and last_plotted != None:
            if last_plotted not in selected_chromatograms:
                selected_chromatograms.append(last_plotted)
                
        if len(selected_chromatograms) == 1:
            clear = True
        
        if len(selected_chromatograms) > 1:
            last_plotted = None
            plot_graph_button.config(state=tk.NORMAL)
            if ax.get_xlim()[0] != 0 and ax.get_xlim()[1] != 1:
                last_xlims_chrom = ax.get_xlim()
                last_ylims_chrom = ax.get_ylim()
            color_number = 0
            clear_plot(ax, canvas)
            for chromatogram in selected_chromatograms:
                if chromatogram[0] != "Base Peak Chromatogram/Electropherogram" or chromatogram[1] != 3:
                        
                    show_graph(chromatogram[0], clear, chromatogram[1])
                    
                elif chromatogram[0] == "Base Peak Chromatogram/Electropherogram":
                    show_graph(item_text, clear)
                        
                color_number+= 1
                if color_number == len(colors):
                    color_number = 0
            ax.legend(fontsize=9)
            
        elif len(chromatograms_list.selection()) == 1:
            item_text = chromatograms_list.item(chromatograms_list.selection(), "text")  # Get the text of the selected item
            last_xlims_chrom = None
            last_ylims_chrom = None
            if item_text != 'Base Peak Chromatogram/Electropherogram':
                clear_plot(ax, canvas)
                level = 0
                parent_item = chromatograms_list.selection()
                while parent_item:
                    level += 1
                    parent_item = chromatograms_list.parent(parent_item)
                if level == 1:
                    if len(samples_dropdown_options) > 1:
                        compare_samples_button.config(state=tk.NORMAL)
                        plot_graph_button.config(state=tk.NORMAL)
                    last_plotted = [f"{item_text}", 1]
                    show_graph(f"{item_text}", clear, 1)
                if level == 2:
                    if len(samples_dropdown_options) > 1:
                        compare_samples_button.config(state=tk.NORMAL)
                        plot_graph_button.config(state=tk.NORMAL)
                    parent_item = chromatograms_list.parent(chromatograms_list.selection())
                    parent_text = chromatograms_list.item(parent_item, "text")
                    last_plotted = [f"{parent_text}+{item_text}", 2]
                    show_graph(f"{parent_text}+{item_text}", clear, 2)
                if level == 3:
                    parent_item = chromatograms_list.parent(chromatograms_list.selection())
                    parent_text = chromatograms_list.item(parent_item, "text")
                    grand_parent_item = chromatograms_list.parent(parent_item)
                    grand_parent_text = chromatograms_list.item(grand_parent_item, "text")
                    last_plotted = [f"{grand_parent_text}+{parent_text}", 2]
                    show_graph(f"{grand_parent_text}+{parent_text}", clear, 3)
            else:
                clear_plot(ax, canvas)
                # clear_plot(ax_spec, canvas_spec)
                last_plotted = [f"{item_text}", 1]
                show_graph(item_text, clear)
        
        else:
            clear_plot(ax, canvas)
            last_plotted = None
        
    def show_graph(item_text, clear = True, level = 0):
        global current_data, coordinate_label, ax, canvas, type_coordinate, colors, color_number, ms2_precursors_actual, spectra_indexes, last_xlims_chrom, last_ylims_chrom, og_x_range, og_y_range, chromatogram_bound, chromatogram_binds, rt_label
        
        if chromatogram_bound:
            for i in chromatogram_binds:
                canvas.mpl_disconnect(i)
            chromatogram_bound = False
            chromatogram_binds = []
        
        rt_label.config(text = f"")
        
        if clear:
            clear_plot(ax, canvas)
                
        if type(item_text) == list:
            x_values = item_text[0]
            y_values = item_text[1]
            label_show_graph = item_text[3]
        else:
            if item_text == "Base Peak Chromatogram/Electropherogram":
                x_values = [i/60 for i in current_data['rt_array']] if current_data['time_unit'] == 'seconds' else current_data['rt_array']
                y_values = current_data['bpc']
                label_show_graph = "Base Peak Chromatogram/Electropherogram"
            elif level == 1 or level == 2 or level == 3:
                sample_index = list(glycans_per_sample.keys()).index(selected_item)
                if level == 1:
                
                    # Load retention times from file
                    x_values = gg_file.get_rt_array(sample_index)
                        
                    y_values = []
                    counter = 0
                    for i in gg_file.list_chromatograms(sample_index):
                        if "+".join(i.split("+")[:-1]) == item_text:
                            
                            # Load intensities from file
                            y_values_temp = gg_file.get_chromatogram(sample_index, i, 'smoothed')
                            
                            if counter == 0:
                                y_values = y_values_temp
                                counter+= 1
                            else:
                                y_values = [x + y for x, y in zip(y_values, y_values_temp)]
                elif level == 2 or level == 3:
                
                    # Load retention times from file
                    x_values = gg_file.get_rt_array(sample_index)
                        
                    # Load intensities from file
                    y_values = gg_file.get_chromatogram(sample_index, item_text, 'smoothed')
                    
                if level == 1:
                    label_show_graph = f"{item_text}"
                if level == 2:
                    parent_item = chromatograms_list.parent(chromatograms_list.focus())
                    parent_text = chromatograms_list.item(parent_item, "text")
                    item_text_split_zero = item_text.split(" - ")[0]
                    label_show_graph = f"{item_text_split_zero}"
                if level == 3:
                    parent_item = chromatograms_list.parent(chromatograms_list.focus())
                    parent_text = chromatograms_list.item(parent_item, "text")
                    grand_parent_item = chromatograms_list.parent(parent_item)
                    grand_parent_text = chromatograms_list.item(grand_parent_item, "text")
                    parent_text_split_zero = parent_text.split(" ")[0]
                    label_show_graph = f"{grand_parent_text} - {parent_text_split_zero}"
            else:
                return
           
        ax.plot(x_values, y_values, linewidth=1, color=item_text[2] if type(item_text) == list else 'red', label = label_show_graph)
        
        ax.set_xlabel('Retention/Migration Time (min)')
        ax.set_ylabel('Intensity (AU)')
        
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        if type(item_text) == list and not clear:
            ax.legend(fontsize=9)
        
        vertical_line = ax.axvline(x=-10000, color='black', linestyle='--', linewidth=1)
        
        try:
            lines = ax.get_lines()

            # Initialize the maximum value to negative infinity
            max_value = float('-inf')

            # Iterate through each line
            for line in lines:
                # Get the y-data of the line
                y_data = line.get_ydata()
                
                # Find the maximum value in the y-data
                line_max_value = max(y_data)
                
                # Update the overall maximum value if needed
                max_value = max(max_value, line_max_value)
                
            x_min = min(x_values)
            x_max = max(x_values)
            og_x_range = (x_min-0.5, x_max+0.5)
            og_y_range = (0, max_value*1.1)
        except:
            og_x_range = (0, 1000)
            og_y_range = (0, 1000)
        
        if level == 3:
            parent_item = chromatograms_list.parent(chromatograms_list.selection())
            parent_text = chromatograms_list.item(parent_item, "text") #adduct
            grand_parent_item = chromatograms_list.parent(parent_item)
            grand_parent_text = chromatograms_list.item(grand_parent_item, "text") #glycan
            
            parent_text_split_zero = parent_text.split(" ")[0]
            chromatograms_list_selected_item_text = chromatograms_list.item(chromatograms_list.selection(), "text")
            if f"{grand_parent_text}+{parent_text_split_zero}_{chromatograms_list_selected_item_text}_RTs" in curve_fittings[sample_index]:
                vlines_coord = [curve_fittings[sample_index][f"{grand_parent_text}+{parent_text_split_zero}_{chromatograms_list_selected_item_text}_RTs"][0], curve_fittings[sample_index][f"{grand_parent_text}+{parent_text_split_zero}_{chromatograms_list_selected_item_text}_RTs"][-1]]
                ax.axvline(vlines_coord[0], color = "blue", linestyle='--', linewidth=1)
                ax.axvline(vlines_coord[1], color = "blue", linestyle='--', linewidth=1)
                ax.set_xlim(vlines_coord[0]-1, vlines_coord[1]+1)
                most_intense = 0
                for i_i, i in enumerate(x_values):
                    if i > vlines_coord[1]:
                        break
                    if i >= vlines_coord[0] and i <= vlines_coord[1]:
                        if y_values[i_i] > most_intense:
                            most_intense = y_values[i_i]
                ax.set_ylim(og_y_range[0], most_intense*1.5)
            else:
                vlines_coord = ''
                if last_xlims_chrom != None:
                    ax.set_xlim(last_xlims_chrom)
                    ax.set_ylim(last_ylims_chrom)
                else:
                    if og_y_range[0] == og_y_range[1]:
                        ax.set_xlim(og_x_range[0], og_x_range[1])
                        ax.set_ylim(0, 1000)
                    else:
                        ax.set_xlim(og_x_range[0], og_x_range[1])
                        ax.set_ylim(og_y_range[0], og_y_range[1])
        else:
            vlines_coord = ''
            if last_xlims_chrom != None:
                ax.set_xlim(last_xlims_chrom)
                ax.set_ylim(og_y_range[0], og_y_range[1])
            else:
                ax.set_xlim(og_x_range[0], og_x_range[1])
                ax.set_ylim(og_y_range[0], og_y_range[1])
                
        if zoom_lock_chromatogram['on']:
            ax.set_xlim(zoom_lock_chromatogram['x_axis'][0] if zoom_lock_chromatogram['x_axis'][0] != 0 else og_x_range[0], zoom_lock_chromatogram['x_axis'][1] if zoom_lock_chromatogram['x_axis'][1] != 0 else og_x_range[1])
            ax.set_ylim(zoom_lock_chromatogram['y_axis'][0] if zoom_lock_chromatogram['y_axis'][0] != 0 else og_y_range[0], zoom_lock_chromatogram['y_axis'][1] if zoom_lock_chromatogram['y_axis'][1] != 0 else og_y_range[1])
        
        canvas.draw()
        
        if not chromatogram_bound:
            zoom_selection_key_press = canvas.mpl_connect('key_press_event', lambda event: zoom_selection(event, ax, canvas, type_coordinate))
            zoom_selection_key_release = canvas.mpl_connect('key_release_event', lambda event: zoom_selection(event, ax, canvas, type_coordinate))
            zoom_selection_motion_notify = canvas.mpl_connect('motion_notify_event', lambda event: zoom_selection(event, ax, canvas, type_coordinate))
            zoom_selection_button_press = canvas.mpl_connect('button_press_event', lambda event: zoom_selection(event, ax, canvas, type_coordinate) if event.button == 1 else None)
            zoom_selection_button_release = canvas.mpl_connect('button_release_event', lambda event: zoom_selection(event, ax, canvas, type_coordinate) if event.button == 1 else None)
            on_scroll_event = canvas.mpl_connect('scroll_event', lambda event: on_scroll(event, ax, canvas, type_coordinate))
            on_double_click_event = canvas.mpl_connect('button_press_event', lambda event: on_double_click(event, ax, canvas, og_x_range, og_y_range, type_coordinate) if event.button == 1 else None)
            on_plot_hover_motion = canvas.mpl_connect('motion_notify_event', lambda event: on_plot_hover(event, ax, canvas, coordinate_label, type_coordinate, vertical_line))
            on_click_press = canvas.mpl_connect('button_press_event', lambda event: on_click(event, ax, x_values, y_values, vlines_coord) if event.button == 1 else None)
            on_click_release = canvas.mpl_connect('button_release_event', lambda event: on_click(event, ax, x_values, y_values, vlines_coord) if event.button == 1 else None)
            right_move_spectra = canvas.mpl_connect('key_press_event', lambda event: on_right_arrow(event, x_values, y_values, vlines_coord) if event.key == 'right' else None)
            left_move_spectra = canvas.mpl_connect('key_press_event', lambda event: on_left_arrow(event, x_values, y_values, vlines_coord) if event.key == 'left' else None)
            on_pan_press = canvas.mpl_connect('button_press_event', lambda event: on_pan(event, ax, canvas, type_coordinate) if event.button == 1 else None)
            on_pan_release = canvas.mpl_connect('button_release_event', lambda event: on_pan(event, ax, canvas, type_coordinate) if event.button == 1 else None)
            on_pan_motion = canvas.mpl_connect('motion_notify_event', lambda event: on_pan(event, ax, canvas, type_coordinate))
            on_pan_right_click_press = canvas.mpl_connect('button_press_event', lambda event: on_pan_right_click(event, ax, canvas, type_coordinate) if event.button == 3 else None)
            on_pan_right_click_release = canvas.mpl_connect('button_release_event', lambda event: on_pan_right_click(event, ax, canvas, type_coordinate) if event.button == 3 else None)
            on_pan_right_click_motion = canvas.mpl_connect('motion_notify_event', lambda event: on_pan_right_click(event, ax, canvas, type_coordinate))
                
            chromatogram_bound = True
            chromatogram_binds = [zoom_selection_key_press, zoom_selection_key_release, zoom_selection_motion_notify, zoom_selection_button_press, zoom_selection_button_release, on_scroll_event, on_double_click_event, on_plot_hover_motion, on_click_press, on_click_release, right_move_spectra, left_move_spectra, on_pan_press, on_pan_release, on_pan_motion, on_pan_right_click_press, on_pan_right_click_release, on_pan_right_click_motion]
            
        toggle_hide_bad_peaks()
                
    def show_graph_spectra(rt, peak_range, custom_lim = None):
        global coordinate_label_spec, ax_spec, canvas_spec, ms2_info, rt_label, ms1_bound, ms1_binds, ms2_bound, ms2_binds, precursor_label, og_x_range_spec, og_y_range_spec
            
        if ms1_bound:
            for i in ms1_binds:
                canvas_spec.mpl_disconnect(i)
            ms1_bound = False
            ms1_binds = []
            
        spectra_notebook.select(0)
            
        if len(samples_list) == 0:
            return
        
        clear_plot(ax_spec, canvas_spec)
        
        #for use with data sensitive to the formatting, such as isotopic peaks highlighting
        if current_data['time_unit'] == 'seconds':
            rt_minutes = round(rt/60, 4)
        else:
            rt_minutes = round(rt, 4)
        
        x_values_spec = current_data['access'].time[rt]['m/z array']
        y_values_spec = current_data['access'].time[rt]['intensity array']
        
        new_x_data = []
        new_y_data = []
        for index, x in enumerate(x_values_spec):
            new_x_data.append(x-0.000001)
            new_x_data.append(x)
            new_x_data.append(x+0.000001)
            new_y_data.append(0)
            new_y_data.append(y_values_spec[index])
            new_y_data.append(0)
        
        ax_spec.plot(new_x_data, new_y_data, marker='None', linewidth=1, color='black', label = 'MS1 Spectrum')
        
        scaling = scaling_dropdown.get()
        if scaling == 'Linear':
            ax_spec.set_yscale('linear')
        if scaling == 'Log':
            ax_spec.set_yscale('symlog')
        if scaling == 'Sqrt':
            ax_spec.set_yscale('function', functions=(custom_sqrt, custom_square))
        
        ax_spec.set_xlabel('m/z')
        ax_spec.set_ylabel('Intensity (AU)')
        
        ax_spec.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax_spec.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        try:
            og_x_range_spec = (min(x_values_spec)-10, max(x_values_spec)+10)
            og_y_range_spec = (0, max(y_values_spec)*1.1)
        except:
            og_x_range_spec = (0, 0)
            og_y_range_spec = (0, 0)
        if custom_lim != None:
            ax_spec.set_xlim(custom_lim[0][0], custom_lim[0][1])
            max_value = 0
            for i_i, i in enumerate(x_values_spec):
                if i > custom_lim[0][1]:
                    break
                if i > custom_lim[0][0] and i < custom_lim[0][1]:
                    if y_values_spec[i_i] > max_value:
                        max_value = y_values_spec[i_i]
            if max_value != 0:
                ax_spec.set_ylim(0, max_value*1.1)
        else:
            if og_x_range_spec[0] == og_x_range_spec[1]:
                ax_spec.set_xlim(0, 1000)
            else:    
                ax_spec.set_xlim(og_x_range_spec[0], og_x_range_spec[1])
                
            if og_y_range_spec[0] == og_y_range_spec[1]:
                ax_spec.set_ylim(0, 1000)
            else:    
                ax_spec.set_ylim(og_y_range_spec[0], og_y_range_spec[1])
        
        if len(peak_range) > 0:
            if rt_minutes >= peak_range[0] and rt_minutes <= peak_range[1]:
                sample_index = list(glycans_per_sample.keys()).index(selected_item)
                parent_item = chromatograms_list.parent(chromatograms_list.selection())
                parent_text = chromatograms_list.item(parent_item, "text") #adduct
                grand_parent_item = chromatograms_list.parent(parent_item)
                grand_parent_text = chromatograms_list.item(grand_parent_item, "text") #glycan
                temp_adduct = parent_text.split(" ")[0]
                
                mz_parent_base = parent_text.split(" ")[-1]
                if rt_minutes in isotopic_fittings[sample_index][f"{grand_parent_text}_{temp_adduct}"]:
                    peaks = isotopic_fittings[sample_index][f"{grand_parent_text}_{temp_adduct}"][rt_minutes][0]
                    
                    adduct_charge = adduct_combos_dict[temp_adduct]['charges']
                    
                    if len(peaks) != 0:
                        ax_spec.set_xlim(peaks[0]-2, peaks[-1]+2)
                        highest = 0
                        for i_i, i in enumerate(x_values_spec):
                            if i > peaks[-1]:
                                break
                            if i >= peaks[0]-0.5 and i <= peaks[-1]+0.5:
                                if y_values_spec[i_i] > highest:
                                    highest = y_values_spec[i_i]
                        if highest != 0:
                            ax_spec.set_ylim(0, highest*1.1)
                            for index, i in enumerate(peaks):
                                mz_chosen = float(mz_parent_base)+(General_Functions.h_mass/abs(adduct_charge))*index
                                calculated_width = (0.1 if gg_analysis_tolerance == None else General_Functions.tolerance_calc(gg_analysis_tolerance[0], gg_analysis_tolerance[1], mz_chosen))
                                ax_spec.add_patch(Rectangle((mz_chosen-calculated_width, ax_spec.get_ylim()[0]), (mz_chosen+calculated_width) - (mz_chosen-calculated_width), 1000000000000, color='#FEB7A1', alpha=0.3))
                                ax_spec.add_patch(Rectangle((i-0.001, ax_spec.get_ylim()[0]), (i+0.001) - (i-0.001), 1000000000000, color='#345eeb', alpha=0.3))
        
        if len(processed_data[selected_item]['ms2']) > 0:
            get_ms2(rt, ax_spec, canvas_spec, x_values_spec, y_values_spec)
                
        if zoom_lock_spectrum['on']:
            ax_spec.set_xlim(zoom_lock_spectrum['x_axis'][0] if zoom_lock_spectrum['x_axis'][0] != 0 else og_x_range_spec[0], zoom_lock_spectrum['x_axis'][1] if zoom_lock_spectrum['x_axis'][1] != 0 else og_x_range_spec[1])
            ax_spec.set_ylim(zoom_lock_spectrum['y_axis'][0] if zoom_lock_spectrum['y_axis'][0] != 0 else og_y_range_spec[0], zoom_lock_spectrum['y_axis'][1] if zoom_lock_spectrum['y_axis'][1] != 0 else og_y_range_spec[1])
            
        annotate_top_y_values(ax_spec, canvas_spec)
        
        canvas_spec.draw()
        
        rt_formatted = "%.2f" % round(rt_minutes, 2)
        rt_label.config(text = f"Retention/Migration Time: {rt_formatted}")
        
        if not ms1_bound:
            zoom_selection_key_press_spec = canvas_spec.mpl_connect('key_press_event', lambda event: zoom_selection_spec(event, ax_spec, canvas_spec, type_coordinate_spec))
            zoom_selection_key_release_spec = canvas_spec.mpl_connect('key_release_event', lambda event: zoom_selection_spec(event, ax_spec, canvas_spec, type_coordinate_spec))
            zoom_selection_motion_notify_spec = canvas_spec.mpl_connect('motion_notify_event', lambda event: zoom_selection_spec(event, ax_spec, canvas_spec, type_coordinate_spec))
            zoom_selection_button_press_spec = canvas_spec.mpl_connect('button_press_event', lambda event: zoom_selection_spec(event, ax_spec, canvas_spec, type_coordinate_spec) if event.button == 1 else None)
            zoom_selection_button_release_spec = canvas_spec.mpl_connect('button_release_event', lambda event: zoom_selection_spec(event, ax_spec, canvas_spec, type_coordinate_spec) if event.button == 1 else None)
            on_scroll_event_spec = canvas_spec.mpl_connect('scroll_event', lambda event: on_scroll(event, ax_spec, canvas_spec, type_coordinate_spec))
            on_double_click_event_spec = canvas_spec.mpl_connect('button_press_event', lambda event: on_double_click(event, ax_spec, canvas_spec, og_x_range_spec, og_y_range_spec, type_coordinate_spec) if event.button == 1 else None)
            on_plot_hover_motion_spec = canvas_spec.mpl_connect('motion_notify_event', lambda event: on_plot_hover(event, ax_spec, canvas_spec, coordinate_label_spec, type_coordinate_spec))
            pick_event_spec = canvas_spec.mpl_connect('pick_event', lambda event: on_pick_spec(event, ms2_info))
            hand_hover_spec = canvas_spec.mpl_connect('motion_notify_event', on_hover_spec)
            on_pan_press_spec = canvas_spec.mpl_connect('button_press_event', lambda event: on_pan(event, ax_spec, canvas_spec, type_coordinate_spec) if event.button == 1 else None)
            on_pan_release_spec = canvas_spec.mpl_connect('button_release_event', lambda event: on_pan(event, ax_spec, canvas_spec, type_coordinate_spec) if event.button == 1 else None)
            on_pan_motion_spec = canvas_spec.mpl_connect('motion_notify_event', lambda event: on_pan(event, ax_spec, canvas_spec, type_coordinate_spec))
            on_pan_right_click_press_spec = canvas_spec.mpl_connect('button_press_event', lambda event: on_pan_right_click(event, ax_spec, canvas_spec, type_coordinate_spec) if event.button == 3 else None)
            on_pan_right_click_release_spec = canvas_spec.mpl_connect('button_release_event', lambda event: on_pan_right_click(event, ax_spec, canvas_spec, type_coordinate_spec) if event.button == 3 else None)
            on_pan_right_click_motion_spec = canvas_spec.mpl_connect('motion_notify_event', lambda event: on_pan_right_click(event, ax_spec, canvas_spec, type_coordinate_spec))
            on_mouse_press_spectra_plot_spec = canvas_spec.mpl_connect('button_press_event', lambda event: on_mouse_press_spectra_plot(event, ax_spec))
            on_mouse_release_spectra_plot_spec = canvas_spec.mpl_connect('button_release_event', lambda event: on_mouse_release_spectra_plot(event, ax_spec))
            select_annotation_spec = canvas_spec.mpl_connect('button_press_event', select_ruler)
            
            ms1_bound = True
            ms1_binds = [zoom_selection_key_press_spec, zoom_selection_key_release_spec, zoom_selection_motion_notify_spec, zoom_selection_button_press_spec, zoom_selection_button_release_spec, on_scroll_event_spec, on_double_click_event_spec, on_plot_hover_motion_spec, pick_event_spec, hand_hover_spec, on_pan_press_spec, on_pan_release_spec, on_pan_motion_spec, on_pan_right_click_press_spec, on_pan_right_click_release_spec, on_pan_right_click_motion_spec, on_mouse_press_spectra_plot_spec, on_mouse_release_spectra_plot_spec, select_annotation_spec]
            
    def show_ms2_graph(spectra_info):
        global ms2_bound, ms2_binds, og_x_range_spec_ms2, og_y_range_spec_ms2
        
        no_spectra_loaded_label_ms2.config(text="")
        no_spectra_loaded_label_ms2.place(relx=0.01, rely=0.9)
        
        if ms2_bound:
            for i in ms2_binds:
                canvas_spec_ms2.mpl_disconnect(i)
            ms2_bound = False
            ms2_binds = []
        
        spectra_notebook.select(1)
        
        x_values_spec_ms2 = current_data['access'].time[spectra_info[1]]['m/z array']
        y_values_spec_ms2 = current_data['access'].time[spectra_info[1]]['intensity array']
        spectrum_id = current_data['access'].time[spectra_info[1]]['index'] if current_data['file_type'] == 'mzml' else int(current_data['access'].time[spectra_info[1]]['num'])-1
        
        new_x_data_spec_ms2 = []
        new_y_data_spec_ms2 = []
        for index, x in enumerate(x_values_spec_ms2):
            new_x_data_spec_ms2.append(x-0.000001)
            new_x_data_spec_ms2.append(x)
            new_x_data_spec_ms2.append(x+0.000001)
            new_y_data_spec_ms2.append(0)
            new_y_data_spec_ms2.append(y_values_spec_ms2[index])
            new_y_data_spec_ms2.append(0)
        
        if len(x_values_spec_ms2) > 0:
            og_x_range_spec_ms2 = [0, x_values_spec_ms2[-1]+50]
            og_y_range_spec_ms2 = [0, max(y_values_spec_ms2)*1.1]
        else:
            og_x_range_spec_ms2 = [0, 1000]
            og_y_range_spec_ms2 = [0, 1000]
            
        clear_plot(ax_spec_ms2, canvas_spec_ms2)
        
        ax_spec_ms2.set_xlim(og_x_range_spec_ms2[0], og_x_range_spec_ms2[1])
        ax_spec_ms2.set_ylim(og_y_range_spec_ms2[0], og_y_range_spec_ms2[1])
        
        ax_spec_ms2.plot(new_x_data_spec_ms2, new_y_data_spec_ms2, marker='', linewidth=1, color='#30034d', label = f'MS2 Spectrum - {spectrum_id}')
        
        scaling = scaling_dropdown.get()
        if scaling == 'Linear':
            ax_spec_ms2.set_yscale('linear')
        if scaling == 'Log':
            ax_spec_ms2.set_yscale('symlog')
        if scaling == 'Sqrt':
            ax_spec_ms2.set_yscale('function', functions=(custom_sqrt, custom_square))
        
        ax_spec_ms2.set_xlabel('m/z')
        ax_spec_ms2.set_ylabel('Intensity (AU)')
        
        ax_spec_ms2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax_spec_ms2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        #Check better later
        if len(x_values_spec_ms2) > 0:
            ax_spec_ms2.plot(spectra_info[0], y_values_spec_ms2[np.abs(x_values_spec_ms2-spectra_info[0]).argmin()], marker='D', color = 'red')
        
        if processed_data[selected_item]['time_unit'] == 'seconds':
            spectra_time_minutes = round(spectra_info[1]/60, 4)
        else:
            spectra_time_minutes = round(spectra_info[1], 4)
        
        precursor_mz_formatted = "%.4f" % spectra_info[0]
        precursor_label.config(text = f"RT/MT: {round(spectra_time_minutes, 2)} Precursor m/z: {precursor_mz_formatted}")
        
        if len(reanalysis_path) > 0:
            if len(fragments_library_ms2_scores[1]) > 0:
                annotate_spectrum_button.config(state=tk.NORMAL)
                files = samples_dropdown.cget('values')
                file_index = files.index(samples_dropdown.get())
                if spectrum_id in fragments_library_ms2_scores[1][file_index].keys():
                    spectrum_score = fragments_library_ms2_scores[1][file_index][spectrum_id]
                    spectrum_score_label.config(text=f"Glycan probability: {spectrum_score}/10")
                    ToolTip(spectrum_score_label, "Glycan probability is based on whether the top 10 most annotated fragments during the analysis were found within this spectrum.")
                else:                    
                    if spectra_time_minutes > parameters_gg['analysis settings']['retention time interval'][0] and spectra_time_minutes < parameters_gg['analysis settings']['retention time interval'][1]:
                        spectrum_score_label.config(text="This spectrum has been annotated.")
                        ToolTip(spectrum_score_label, "This spectrum has been annotated during the analysis, as the precursor matched a glycan in the library.")
                    else:
                        spectrum_score_label.config(text="Spectrum outside analyzed RT/MT range.")
                        ToolTip(spectrum_score_label, "This spectrum was not within the analyzed RT/MT range and thus no probability of it being the fragment spectrum of a glycan has  been calculated.")
            else:
                annotate_spectrum_button.config(state=tk.DISABLED)
                spectrum_score_label.config(text="Manual spectrum annotation unavailable.")
                ToolTip(spectrum_score_label, "This analysis was performed before the implementation of this feature or the MS2 was not analyzed.")
            
        custom_font_annotation = {'family': 'sans-serif', 'color': 'black', 'size': 9}
        
        if len(reanalysis_path) > 0:
            data_list = ms2_spectra_glycans[selected_item].get(spectra_time_minutes, None)
            
            if data_list and len(data_list) > 0:
                plotted_glycans = [line.get_label() for line in ax.get_lines() if line.get_label()[0] != "_"]
                
                for glycan in plotted_glycans:
                    glycan_splitted = glycan.split(" - ")
                    if len(glycan_splitted) == 1:
                        glycan_splitted = glycan.split("+")
                        if len(glycan_splitted) != 1:
                            glycan, adduct = glycan_splitted
                        else:
                            adduct = None
                    else:
                        glycan, adduct = glycan_splitted
                        
                    if glycan in data_list:
                        break
                    else:
                        glycan = None
                        adduct = None
                        
                if not glycan:
                    glycan = random.choice(list(data_list.keys()))
                    
                if not adduct or adduct not in data_list[glycan]:
                    adduct = random.choice(data_list[glycan])
                            
                for index, fragment_name in enumerate(glycans_per_sample[selected_item][glycan][adduct]['ms2'][spectra_time_minutes][2]):
                        frag_label = f"{fragment_name}\n{glycans_per_sample[selected_item][glycan][adduct]['ms2'][spectra_time_minutes][0][index]}"
                        ms2_marker = ax_spec_ms2.plot(
                                                      glycans_per_sample[selected_item][glycan][adduct]['ms2'][spectra_time_minutes][0][index], 
                                                      glycans_per_sample[selected_item][glycan][adduct]['ms2'][spectra_time_minutes][1][index], 
                                                      marker='*', 
                                                      markersize=4.5, 
                                                      label = frag_label, 
                                                      color="red"
                                                      )
            
                new_precursor_label = precursor_label.cget('text')+f" Composition: {glycan}"
                precursor_label.config(text = new_precursor_label)
                
        if zoom_lock_spectrum['on']:
            ax_spec_ms2.set_xlim(zoom_lock_spectrum['x_axis'][0] if zoom_lock_spectrum['x_axis'][0] != 0 else og_x_range_spec_ms2[0], zoom_lock_spectrum['x_axis'][1] if zoom_lock_spectrum['x_axis'][1] != 0 else og_x_range_spec_ms2[1])
            ax_spec_ms2.set_ylim(zoom_lock_spectrum['y_axis'][0] if zoom_lock_spectrum['y_axis'][0] != 0 else og_y_range_spec_ms2[0], zoom_lock_spectrum['y_axis'][1] if zoom_lock_spectrum['y_axis'][1] != 0 else og_y_range_spec_ms2[1])

        annotate_top_y_values(ax_spec_ms2, canvas_spec_ms2)
        
        canvas_spec_ms2.draw()
        
        if not ms2_bound:
            zoom_selection_key_press_spec_ms2 = canvas_spec_ms2.mpl_connect('key_press_event', lambda event: zoom_selection_spec_ms2(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec))
            zoom_selection_key_release_spec_ms2 = canvas_spec_ms2.mpl_connect('key_release_event', lambda event: zoom_selection_spec_ms2(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec))
            zoom_selection_motion_notify_spec_ms2 = canvas_spec_ms2.mpl_connect('motion_notify_event', lambda event: zoom_selection_spec_ms2(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec))
            zoom_selection_button_press_spec_ms2 = canvas_spec_ms2.mpl_connect('button_press_event', lambda event: zoom_selection_spec_ms2(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec) if event.button == 1 else None)
            zoom_selection_button_release_spec_ms2 = canvas_spec_ms2.mpl_connect('button_release_event', lambda event: zoom_selection_spec_ms2(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec) if event.button == 1 else None)
            on_scroll_event_spec_ms2 = canvas_spec_ms2.mpl_connect('scroll_event', lambda event: on_scroll(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec))
            on_double_click_event_spec_ms2 = canvas_spec_ms2.mpl_connect('button_press_event', lambda event: on_double_click(event, ax_spec_ms2, canvas_spec_ms2, og_x_range_spec_ms2, og_y_range_spec_ms2, type_coordinate_spec) if event.button == 1 else None)
            on_plot_hover_motion_spec_ms2 = canvas_spec_ms2.mpl_connect('motion_notify_event', lambda event: on_plot_hover(event, ax_spec_ms2, canvas_spec_ms2, coordinate_label_spec, type_coordinate_spec))
            pick_event_spec_ms2 = canvas_spec_ms2.mpl_connect('pick_event', lambda event: on_pick_spec_ms2(event, ms2_info))
            on_pan_press_spec_ms2 = canvas_spec_ms2.mpl_connect('button_press_event', lambda event: on_pan(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec) if event.button == 1 else None)
            on_pan_release_spec_ms2 = canvas_spec_ms2.mpl_connect('button_release_event', lambda event: on_pan(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec) if event.button == 1 else None)
            on_pan_motion_spec_ms2 = canvas_spec_ms2.mpl_connect('motion_notify_event', lambda event: on_pan(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate_spec))
            on_pan_right_click_press_spec_ms2 = canvas_spec_ms2.mpl_connect('button_press_event', lambda event: on_pan_right_click(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate) if event.button == 3 else None)
            on_pan_right_click_release_spec_ms2 = canvas_spec_ms2.mpl_connect('button_release_event', lambda event: on_pan_right_click(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate) if event.button == 3 else None)
            on_pan_right_click_motion_spec_ms2 = canvas_spec_ms2.mpl_connect('motion_notify_event', lambda event: on_pan_right_click(event, ax_spec_ms2, canvas_spec_ms2, type_coordinate))
            on_mouse_press_spectra_plot_spec_ms2 = canvas_spec_ms2.mpl_connect('button_press_event', lambda event: on_mouse_press_spectra_plot(event, ax_spec_ms2))
            on_mouse_release_spectra_plot_spec_ms2 = canvas_spec_ms2.mpl_connect('button_release_event', lambda event: on_mouse_release_spectra_plot(event, ax_spec_ms2))
            select_annotation_spec_ms2 = canvas_spec_ms2.mpl_connect('button_press_event', select_ruler)
            
            ms2_bound = True
            ms2_binds = [zoom_selection_key_press_spec_ms2, zoom_selection_key_release_spec_ms2, zoom_selection_motion_notify_spec_ms2, zoom_selection_button_press_spec_ms2, zoom_selection_button_release_spec_ms2, on_scroll_event_spec_ms2, on_double_click_event_spec_ms2, on_plot_hover_motion_spec_ms2, pick_event_spec_ms2, on_pan_press_spec_ms2, on_pan_release_spec_ms2, on_pan_motion_spec_ms2, on_pan_right_click_press_spec_ms2, on_pan_right_click_release_spec_ms2, on_pan_right_click_motion_spec_ms2, on_mouse_press_spectra_plot_spec_ms2, on_mouse_release_spectra_plot_spec_ms2, select_annotation_spec_ms2]
        
    def get_ms2(rt_here, ax_spec, canvas_spec, x_values_spec, y_values_spec):
        global ms2_info
        ms2_info = {}
        
        def find_closest_value(x, y, value):
            idx = np.abs(x - value).argmin()
            return x[idx], y[idx], idx
            
        ms2_precursors = []
        peaks_indexes = []
        for i in current_data['ms2'][rt_here]:
            for k in current_data['ms2'][rt_here][i]:
                ms2_info_current = find_closest_value(x_values_spec, y_values_spec, k)
                ms2_info[k] = [k, i]
                ax_spec.plot(k, ms2_info_current[1], marker='D', color="purple", picker=5)
            
    def on_pick_spec(event, ms2_info):
        global on_pan_press, on_pan_motion, on_pan_release, on_pan_press_spec, on_pan_motion_spec, on_pan_release_spec, panning_enabled
        if event.mouseevent.button == 1:
            panning_enabled = False
            # Get the artist (marker) that was picked
            artist = event.artist
            # Get the index of the picked point
            index = event.ind[0]
            # Get the x and y coordinates of the picked point
            x = artist.get_xdata()[index]
            y = artist.get_ydata()[index]
            
            show_ms2_graph(ms2_info[x])

    def on_hover_spec(event):
        # Change cursor to a hand when hovering over markers
        if event.inaxes:
            for artist in event.inaxes.get_children():
                if isinstance(artist, plt.Line2D) and artist.contains(event)[0]:
                    if artist.get_marker() == 'D':
                        event.canvas.get_tk_widget().configure(cursor="hand2")
                        return
        event.canvas.get_tk_widget().configure(cursor="")  
        
    def on_left_arrow(event, x_values_here, y_values_here, peak_range = ''):
        global marker_spectra_x_value, marker_spectra_y_value, two_d
        x_lim = ax_spec.get_xlim()
        y_lim = ax_spec.get_ylim()
        artists = ax.get_children()
        former_rt = 0
        found = False
        for artist in artists:
            if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                former_rt = artist.get_xdata()[0]
                found = True
                artist.remove()
                break
        if not found:
            return
        if x_values_here.index(former_rt) == 0:
            marker_spectra_y_value = y_values_here[x_values_here.index(former_rt)]
            marker_spectra_x_value = x_values_here[x_values_here.index(former_rt)]
            marker_spectra = ax.plot(marker_spectra_x_value, marker_spectra_y_value+((ax.get_ylim()[1]-ax.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
            canvas.draw_idle()
        else:
            marker_spectra_y_value = y_values_here[x_values_here.index(former_rt)-1]
            marker_spectra_x_value = x_values_here[x_values_here.index(former_rt)-1]
            marker_spectra = ax.plot(marker_spectra_x_value, marker_spectra_y_value+((ax.get_ylim()[1]-ax.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
            canvas.draw_idle()
            if len(samples_list) > 0:
                if selected_item in processed_data:
                    show_graph_spectra(current_data['rt_array'][x_values_here.index(former_rt)-1], peak_range, (x_lim, y_lim))
                
    def on_right_arrow(event, x_values_here, y_values_here, peak_range = ''):
        global marker_spectra_x_value, marker_spectra_y_value, two_d
        x_lim = ax_spec.get_xlim()
        y_lim = ax_spec.get_ylim()
        artists = ax.get_children()
        former_rt = 0
        found = False
        for artist in artists:
            if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                former_rt = artist.get_xdata()[0]
                found = True
                artist.remove()
                break
        if not found:
            return
        if x_values_here.index(former_rt) == len(x_values_here)-1:
            marker_spectra_y_value = y_values_here[x_values_here.index(former_rt)]
            marker_spectra_x_value = x_values_here[x_values_here.index(former_rt)]
            marker_spectra = ax.plot(marker_spectra_x_value, marker_spectra_y_value+((ax.get_ylim()[1]-ax.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
            canvas.draw_idle()
        else:
            marker_spectra_y_value = y_values_here[x_values_here.index(former_rt)+1]
            marker_spectra_x_value = x_values_here[x_values_here.index(former_rt)+1]
            marker_spectra = ax.plot(marker_spectra_x_value, marker_spectra_y_value+((ax.get_ylim()[1]-ax.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
            canvas.draw_idle()
            if len(samples_list) > 0:
                if selected_item in processed_data:
                    show_graph_spectra(current_data['rt_array'][x_values_here.index(former_rt)+1], peak_range, (x_lim, y_lim))
        
    def zoom_selection(event, ax_here, canvas_here, type_coordinate):
        global marker_spectra_x_value, marker_spectra_y_value, rect
        
        edgecolor = (0, 0, 0.9, 0.1)  # RGBA format: (red, green, blue, alpha)
        
        if not hasattr(zoom_selection, 'origin'):
            zoom_selection.origin = None
        if not hasattr(zoom_selection, 'is_pressed'):
            zoom_selection.is_pressed = False
        
        if event.name == 'button_press_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if event.xdata is not None and event.ydata is not None and event.button == 1:
                
                artists = ax_here.get_children()
                found = False
                for artist in artists:
                    if isinstance(artist, matplotlib.patches.Rectangle) and artist.get_edgecolor() == edgecolor:
                        artist.set_visible(True)
                        zoom_selection.origin = (event.xdata, event.ydata)
                        zoom_selection.is_pressed = True
                        rect.set_xy((event.xdata, 0))
                        rect.set_width(0)
                        rect.set_height(ax_here.get_ylim()[1])
                        found = True
                if not found:

                    # Create the rectangle with the semi-transparent blue background
                    rect = patches.Rectangle((0, 0), 0, 0, fill=True, edgecolor=edgecolor, facecolor=edgecolor, visible=True)
                    ax_here.add_patch(rect)
                    zoom_selection.origin = (event.xdata, event.ydata)
                    zoom_selection.is_pressed = True
                    rect.set_xy((event.xdata, 0))
                    rect.set_width(0)
                    rect.set_height(ax_here.get_ylim()[1])
        elif event.name == 'button_release_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and zoom_selection.origin is not None and event.xdata is not None and event.ydata is not None and event.button == 1:
                zoom_selection.is_pressed = False
                rect.set_visible(False)
                xlim = sorted([zoom_selection.origin[0], event.xdata])
                ax_here.set_xlim(*xlim)
                x_data = ax_here.get_lines()[0].get_xdata()
                y_data = ax_here.get_lines()[0].get_ydata()
                indices = np.where((x_data >= xlim[0]) & (x_data <= xlim[1]))
                if len(indices[0]) > 0:
                    if type_coordinate == 'chromatogram':
                        lines = ax_here.get_lines()
                        max_value = float('-inf')
                        for line in lines:
                            y_check_data = line.get_ydata()
                            try:
                                line_max_value = np.max(y_check_data[indices])
                                max_value = max(max_value, line_max_value)
                            except:
                                pass
                        if max_value != float('-inf'):
                            max_y_value = max_value
                    else:
                        max_y_value = np.max(y_data[indices])
                    ax_here.set_ylim(0, max_y_value + 0.1 * max_y_value)  # Adjust y-axis limit
                    artists = ax_here.get_children()
                    for artist in artists:
                        if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                            artist.remove()
                            marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+(ax_here.get_ylim()[1]*0.05), marker='v', markersize=5, color='black')
                            break
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
        elif event.name == 'motion_notify_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and event.xdata is not None and event.ydata is not None:
                x0, _ = zoom_selection.origin
                x1, _ = event.xdata, event.ydata
                rect.set_width(x1 - x0)
        
    def zoom_selection_spec(event, ax_here, canvas_here, type_coordinate):
        global marker_spectra_x_value, marker_spectra_y_value, rect_spec
        
        edgecolor = (0, 0, 0.9, 0.1)  # RGBA format: (red, green, blue, alpha)
        
        if not hasattr(zoom_selection, 'origin'):
            zoom_selection.origin = None
        if not hasattr(zoom_selection, 'is_pressed'):
            zoom_selection.is_pressed = False
        
        if event.name == 'button_press_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if event.xdata is not None and event.ydata is not None and event.button == 1:
                
                artists = ax_here.get_children()
                found = False
                for artist in artists:
                    if isinstance(artist, matplotlib.patches.Rectangle) and artist.get_edgecolor() == edgecolor:
                        artist.set_visible(True)
                        zoom_selection.origin = (event.xdata, event.ydata)
                        zoom_selection.is_pressed = True
                        rect_spec.set_xy((event.xdata, 0))
                        rect_spec.set_width(0)
                        rect_spec.set_height(ax_here.get_ylim()[1])
                        found = True
                if not found:

                    # Create the rectangle with the semi-transparent blue background
                    rect_spec = patches.Rectangle((0, 0), 0, 0, fill=True, edgecolor=edgecolor, facecolor=edgecolor, visible=True)
                    ax_here.add_patch(rect_spec)
                    zoom_selection.origin = (event.xdata, event.ydata)
                    zoom_selection.is_pressed = True
                    rect_spec.set_xy((event.xdata, 0))
                    rect_spec.set_width(0)
                    rect_spec.set_height(ax_here.get_ylim()[1])
        elif event.name == 'button_release_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and zoom_selection.origin is not None and event.xdata is not None and event.ydata is not None and event.button == 1:
                zoom_selection.is_pressed = False
                rect_spec.set_visible(False)
                xlim = sorted([zoom_selection.origin[0], event.xdata])
                ax_here.set_xlim(*xlim)
                x_data = ax_here.get_lines()[0].get_xdata()
                y_data = ax_here.get_lines()[0].get_ydata()
                indices = np.where((x_data >= xlim[0]) & (x_data <= xlim[1]))
                if len(indices[0]) > 0:
                    if type_coordinate == 'chromatogram':
                        lines = ax_here.get_lines()
                        max_value = float('-inf')
                        for line in lines:
                            y_check_data = line.get_ydata()
                            try:
                                line_max_value = np.max(y_check_data[indices])
                                max_value = max(max_value, line_max_value)
                            except:
                                pass
                        if max_value != float('-inf'):
                            max_y_value = max_value
                    else:
                        max_y_value = np.max(y_data[indices])
                    ax_here.set_ylim(0, max_y_value + 0.1 * max_y_value)  # Adjust y-axis limit
                    artists = ax_here.get_children()
                    for artist in artists:
                        if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                            artist.remove()
                            marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+(ax_here.get_ylim()[1]*0.05), marker='v', markersize=5, color='black')
                            break
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
        elif event.name == 'motion_notify_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and event.xdata is not None and event.ydata is not None:
                x0, _ = zoom_selection.origin
                x1, _ = event.xdata, event.ydata
                rect_spec.set_width(x1 - x0)
                
    def zoom_selection_spec_ms2(event, ax_here, canvas_here, type_coordinate):
        global marker_spectra_x_value, marker_spectra_y_value, rect_ms2
        
        edgecolor = (0, 0, 0.9, 0.1)  # RGBA format: (red, green, blue, alpha)
        
        if not hasattr(zoom_selection, 'origin'):
            zoom_selection.origin = None
        if not hasattr(zoom_selection, 'is_pressed'):
            zoom_selection.is_pressed = False
        
        if event.name == 'button_press_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if event.xdata is not None and event.ydata is not None and event.button == 1:
                
                artists = ax_here.get_children()
                found = False
                for artist in artists:
                    if isinstance(artist, matplotlib.patches.Rectangle) and artist.get_edgecolor() == edgecolor:
                        artist.set_visible(True)
                        zoom_selection.origin = (event.xdata, event.ydata)
                        zoom_selection.is_pressed = True
                        rect_ms2.set_xy((event.xdata, 0))
                        rect_ms2.set_width(0)
                        rect_ms2.set_height(ax_here.get_ylim()[1])
                        found = True
                if not found:

                    # Create the rectangle with the semi-transparent blue background
                    rect_ms2 = patches.Rectangle((0, 0), 0, 0, fill=True, edgecolor=edgecolor, facecolor=edgecolor, visible=True)
                    ax_here.add_patch(rect_ms2)
                    zoom_selection.origin = (event.xdata, event.ydata)
                    zoom_selection.is_pressed = True
                    rect_ms2.set_xy((event.xdata, 0))
                    rect_ms2.set_width(0)
                    rect_ms2.set_height(ax_here.get_ylim()[1])
        elif event.name == 'button_release_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and zoom_selection.origin is not None and event.xdata is not None and event.ydata is not None and event.button == 1:
                zoom_selection.is_pressed = False
                rect_ms2.set_visible(False)
                xlim = sorted([zoom_selection.origin[0], event.xdata])
                ax_here.set_xlim(*xlim)
                x_data = ax_here.get_lines()[0].get_xdata()
                y_data = ax_here.get_lines()[0].get_ydata()
                indices = np.where((x_data >= xlim[0]) & (x_data <= xlim[1]))
                if len(indices[0]) > 0:
                    if type_coordinate == 'chromatogram':
                        lines = ax_here.get_lines()
                        max_value = float('-inf')
                        for line in lines:
                            y_check_data = line.get_ydata()
                            try:
                                line_max_value = np.max(y_check_data[indices])
                                max_value = max(max_value, line_max_value)
                            except:
                                pass
                        if max_value != float('-inf'):
                            max_y_value = max_value
                    else:
                        max_y_value = np.max(y_data[indices])
                    ax_here.set_ylim(0, max_y_value + 0.1 * max_y_value)  # Adjust y-axis limit
                    artists = ax_here.get_children()
                    for artist in artists:
                        if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                            artist.remove()
                            marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+(ax_here.get_ylim()[1]*0.05), marker='v', markersize=5, color='black')
                            break
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
        elif event.name == 'motion_notify_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and event.xdata is not None and event.ydata is not None:
                x0, _ = zoom_selection.origin
                x1, _ = event.xdata, event.ydata
                rect_ms2.set_width(x1 - x0)
                        
    def zoom_selection_compare(event, ax_here, canvas_here, type_coordinate):
        global marker_spectra_x_value, marker_spectra_y_value, rect_compare
        
        edgecolor = (0, 0, 0.9, 0.1)  # RGBA format: (red, green, blue, alpha)
        
        if not hasattr(zoom_selection, 'origin'):
            zoom_selection.origin = None
        if not hasattr(zoom_selection, 'is_pressed'):
            zoom_selection.is_pressed = False
        
        if event.name == 'button_press_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if event.xdata is not None and event.ydata is not None and event.button == 1:
                
                artists = ax_here.get_children()
                found = False
                for artist in artists:
                    if isinstance(artist, matplotlib.patches.Rectangle) and artist.get_edgecolor() == edgecolor:
                        artist.set_visible(True)
                        zoom_selection.origin = (event.xdata, event.ydata)
                        zoom_selection.is_pressed = True
                        rect_compare.set_xy((event.xdata, 0))
                        rect_compare.set_width(0)
                        rect_compare.set_height(ax_here.get_ylim()[1])
                        found = True
                if not found:

                    # Create the rectangle with the semi-transparent blue background
                    rect_compare = patches.Rectangle((0, 0), 0, 0, fill=True, edgecolor=edgecolor, facecolor=edgecolor, visible=True)
                    ax_here.add_patch(rect_compare)
                    zoom_selection.origin = (event.xdata, event.ydata)
                    zoom_selection.is_pressed = True
                    rect_compare.set_xy((event.xdata, 0))
                    rect_compare.set_width(0)
                    rect_compare.set_height(ax_here.get_ylim()[1])
        elif event.name == 'button_release_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and zoom_selection.origin is not None and event.xdata is not None and event.ydata is not None and event.button == 1:
                zoom_selection.is_pressed = False
                rect_compare.set_visible(False)
                xlim = sorted([zoom_selection.origin[0], event.xdata])
                ax_here.set_xlim(*xlim)
                x_data = ax_here.get_lines()[0].get_xdata()
                y_data = ax_here.get_lines()[0].get_ydata()
                indices = np.where((x_data >= xlim[0]) & (x_data <= xlim[1]))
                if len(indices[0]) > 0:
                    if type_coordinate == 'chromatogram':
                        lines = ax_here.get_lines()
                        max_value = float('-inf')
                        for line in lines:
                            y_check_data = line.get_ydata()
                            try:
                                line_max_value = np.max(y_check_data[indices])
                                max_value = max(max_value, line_max_value)
                            except:
                                pass
                        if max_value != float('-inf'):
                            max_y_value = max_value
                    else:
                        max_y_value = np.max(y_data[indices])
                    ax_here.set_ylim(0, max_y_value + 0.1 * max_y_value)  # Adjust y-axis limit
                    artists = ax_here.get_children()
                    for artist in artists:
                        if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                            artist.remove()
                            marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+(ax_here.get_ylim()[1]*0.05), marker='v', markersize=5, color='black')
                            break
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
        elif event.name == 'motion_notify_event' and (key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False)):
            if zoom_selection.is_pressed and event.xdata is not None and event.ydata is not None:
                x0, _ = zoom_selection.origin
                x1, _ = event.xdata, event.ydata
                rect_compare.set_width(x1 - x0)
                
    def adjust_subplot_size(event, ax_here, canvas_here):
        # Get the current size of the graph frame
        frame_width = event.width
        frame_height = event.height
        
        # Define margin sizes (in pixels)
        left_margin = 65
        right_margin = 10
        top_margin = 21
        bottom_margin = 45
        
        # Calculate the position of the subplot relative to the frame size
        subplot_width = (frame_width - left_margin - right_margin) / frame_width
        subplot_height = (frame_height - top_margin - bottom_margin) / frame_height
        subplot_left = left_margin / frame_width
        subplot_bottom = bottom_margin / frame_height
        
        # Set the position of the subplot
        ax_here.set_position([subplot_left, subplot_bottom, subplot_width, subplot_height])

    def on_scroll(event, ax_here, canvas_here, type_coordinate):
        global marker_spectra_x_value, marker_spectra_y_value
        
        x_min, x_max = ax_here.get_xlim()
        y_min, y_max = ax_here.get_ylim()
        x_range = x_max - x_min
        y_range = y_max - y_min
        
        if event.button == 'up':
            scale_factor = 0.9
        elif event.button == 'down':
            scale_factor = 1.1
        
        # Convert pixel coordinates to data coordinates
        x_data = ax_here.transData.inverted().transform((event.x, 0))[0]
        y_data = ax_here.transData.inverted().transform((0, event.y))[1]
        
        # Check if mouse is over x-axis or y-axis
        over_x = event.y > ax_here.bbox.y1 or event.y < ax_here.bbox.y0
        over_y = event.x > ax_here.bbox.x1 or event.x < ax_here.bbox.x0
        
        if over_x:
            new_x_min = x_data - (x_data - x_min) * scale_factor
            new_x_max = x_data + (x_max - x_data) * scale_factor
            ax_here.set_xlim(new_x_min, new_x_max)
            
            # Automatically adjust y-axis based on the most intense y-value in the zoomed region
            x_data_zoomed = ax_here.get_lines()[0].get_xdata()
            y_data_zoomed = ax_here.get_lines()[0].get_ydata()
            indices = np.where((x_data_zoomed >= new_x_min) & (x_data_zoomed <= new_x_max))
            if len(indices[0]) > 0:
                if type_coordinate == 'chromatogram':
                    lines = ax_here.get_lines()
                    max_value = float('-inf')
                    for line in lines:
                        y_check_data = line.get_ydata()
                        try:
                            line_max_value = np.max(y_check_data[indices])
                            max_value = max(max_value, line_max_value)
                        except:
                            pass
                    if max_value != float('-inf'):
                        max_y_value = max_value
                else:
                    max_y_value = np.max(y_data_zoomed[indices])
                ax_here.set_ylim(0, max_y_value + 0.1 * max_y_value)  # Adjust y-axis limit
                artists = ax_here.get_children()
                for artist in artists:
                    if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                        artist.remove()
                        marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+((ax_here.get_ylim()[1]-ax_here.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
                        break
            
        elif over_y:
            new_y_min = y_min
            new_y_max = y_data + (y_max - y_data) * scale_factor
            ax_here.set_ylim(new_y_min, new_y_max)
            artists = ax_here.get_children()
            for artist in artists:
                if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                    artist.remove()
                    marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+((ax_here.get_ylim()[1]-ax_here.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
                    break
        else:
            new_x_min = x_data - (x_data - x_min) * scale_factor
            new_x_max = x_data + (x_max - x_data) * scale_factor
            ax_here.set_xlim(new_x_min, new_x_max)
            new_y_min = y_data - (y_data - y_min) * scale_factor
            new_y_max = y_data + (y_max - y_data) * scale_factor
            ax_here.set_ylim(new_y_min, new_y_max)
            artists = ax_here.get_children()
            for artist in artists:
                if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                    artist.remove()
                    marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+((ax_here.get_ylim()[1]-ax_here.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
                    break
        if type_coordinate == 'spectra':
            annotate_top_y_values(ax_here, canvas_here)
        canvas_here.draw_idle()
        
    def on_double_click(event, ax_here, canvas_here, og_x_range_here, og_y_range_here, type_coordinate):
        global marker_spectra_x_value, marker_spectra_y_value
        over_x = event.y > ax_here.bbox.y1 or event.y < ax_here.bbox.y0
        over_y = event.x > ax_here.bbox.x1 or event.x < ax_here.bbox.x0
        if event.dblclick:
            if over_x and event.button == 1:
                ax_here.set_xlim(og_x_range_here[0], og_x_range_here[1])  # Reset x-axis limits
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
                canvas_here.draw_idle()
            elif over_y and event.button == 1:
                ax_here.set_ylim(og_y_range_here[0], og_y_range_here[1])  # Reset y-axis limits
                artists = ax_here.get_children()
                for artist in artists:
                    if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                        artist.remove()
                        marker_spectra = ax_here.plot(marker_spectra_x_value, marker_spectra_y_value+((ax_here.get_ylim()[1]-ax_here.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
                        break
                if type_coordinate == 'spectra':
                    annotate_top_y_values(ax_here, canvas_here)
                canvas_here.draw_idle()
            else:
                if gg_draw_on:
                    for index, glycan in enumerate(gg_draw_list):
                        drawing = glycan[1]
                        cont, _ = drawing.contains(event)
                        if cont:
                            gg_draw_options(index)
                            break
    
    def on_plot_hover(event, ax_here, canvas_here, coordinate_label, type_coordinate, vertical_line = None, multiple_signal = False):
        global chromatograms_list
        
        graph_lines = [line for line in ax_here.get_lines() if line.get_label()[0] != "_"]
        if len(graph_lines) == 0:
            return
        x_values_here = graph_lines[0].get_xdata()
        y_values_here = graph_lines[0].get_ydata()
        if type_coordinate == 'spectra':
            mask = y_values_here != 0
            x_values_here = x_values_here[mask]
            y_values_here = y_values_here[mask]
        
        over_x = event.y > ax_here.bbox.y1 or event.y < ax_here.bbox.y0
        over_y = event.x > ax_here.bbox.x1 or event.x < ax_here.bbox.x0
        def remove_marker():
            artists = ax_here.get_children()
            for artist in artists:
                if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 4:
                    artist.remove()
                    break
                    
        if event.inaxes is None:
            if len(chromatograms_list.selection()) > 1 and vertical_line != None:
                vertical_line.set_xdata([-100000])
            elif multiple_signal:
                artists = ax_here.get_children()
                for artist in artists:
                    if isinstance(artist, matplotlib.lines.Line2D) and artist.get_linestyle() == '--':
                        artist.set_xdata([-100000])
            else:
                remove_marker()
            canvas_here.draw_idle()
            coordinate_label.config(text="")
            return

        x_data, y_data = event.xdata, event.ydata
        if not over_x and not over_y:
            if len(x_values_here) == 0:
                return
            closest_coord = np.argmin(np.abs(x_values_here - x_data))
                
            if 0 <= closest_coord < len(y_values_here):
                x_value = x_values_here[closest_coord]
                y_value = y_values_here[closest_coord]  # Get the corresponding y value from y_values_here list
                
                # Update coordinate label text
                if type_coordinate == "chromatogram":
                    if len(chromatograms_list.selection()) > 1 or multiple_signal:
                        coordinate_label.config(text=f"RT: {x_value:.2f}")
                    else:
                        coordinate_label.config(text=f"RT: {x_value:.2f}, Int: {y_value:.2f}")
                    
                elif type_coordinate == "spectra":
                    coordinate_label.config(text=f"m/z: {x_value:.2f}, Int: {y_value:.2f}")
                
                if len(chromatograms_list.selection()) > 1 and vertical_line != None:
                    vertical_line.set_xdata([x_value])
                elif multiple_signal:
                    artists = ax_here.get_children()
                    for artist in artists:
                        if isinstance(artist, matplotlib.lines.Line2D) and artist.get_linestyle() == '--':
                            artist.set_xdata([x_value])
                else:
                    remove_marker()
                    marker = ax_here.plot(x_value, y_value, 'bo', markersize=4)
                
                canvas_here.draw_idle()
            else:
                coordinate_label.config(text="")
        else:
            coordinate_label.config(text="")
            
    def on_click(event, ax, x_values, y_values, peak_range = ''):
        global marker_spectra_x_value, marker_spectra_y_value, gg_draw_on, gg_draw_list, gg_draw_selected
        if event.name == 'button_press_event':
        
            if gg_draw_on:
                if gg_draw_selected == None:
                    for index, glycan in enumerate(gg_draw_list):
                        drawing = glycan[1]
                        cont, _ = drawing.contains(event)
                        if cont:
                            gg_draw_selected = index
                            return
                    gg_draw_selected = None
                    
            # Store the initial coordinates when the button is pressed
            on_click.press_coords = (event.x, event.y)
        if on_click.press_coords == None:
            return
        elif event.name == 'button_release_event':
            if gg_draw_on:
                gg_draw_selected = None
                
            # Calculate the distance moved
            dx = abs(event.x - on_click.press_coords[0])
            dy = abs(event.y - on_click.press_coords[1])
            
            # If the mouse movement is too large, consider it a drag, not a click
            if dx > 5 or dy > 5:
                return
                        
            over_x = event.y > ax.bbox.y1 or event.y < ax.bbox.y0
            over_y = event.x > ax.bbox.x1 or event.x < ax.bbox.x0
            if not over_x and not over_y:
                x_coord = event.xdata
                distance = 999
                closest_coord = 0
                for i_i, i in enumerate(x_values):
                    if i > x_coord+0.5:
                        break
                    if abs(i-x_coord) < distance:
                        distance = abs(i-x_coord)
                        closest_coord = i_i
                # Do something with the x-coordinate
                artists = ax.get_children()
                for artist in artists:
                    if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 5:
                        artist.remove()
                        break
                marker_spectra_y_value = y_values[closest_coord]
                marker_spectra_x_value = x_values[closest_coord]
                marker_spectra = ax.plot(marker_spectra_x_value, marker_spectra_y_value+((ax.get_ylim()[1]-ax.get_ylim()[0])*0.03), marker='v', markersize=5, color='black')
                if len(samples_list) > 0:
                    if selected_item in processed_data:
                        show_graph_spectra(current_data['rt_array'][closest_coord], peak_range)
                canvas.draw_idle()
                
    # Initialize press_coords attribute
    on_click.press_coords = None
            
    def exit_window(direct_close = False):
        def ok_close_main_window():
            # Remove GG windows
            main_window.quit()
            main_window.destroy()
            
            # Remove the execution holder txt file from temp folder
            this_process_id = os.getpid()
            general_temp_folder = os.path.join(tempfile.gettempdir())
            os.remove(os.path.join(general_temp_folder, f"gg_{this_process_id}.txt"))
            try:
                os.remove(os.path.join(general_temp_folder, f"mzml_window_{this_process_id}.txt"))
            except:
                pass
            
            # Clean the temp folder
            clean_temp_folder()
            
            # Finish exitting the Python script itself
            os._exit(0)
        
        def close_exit_window():
            exit_window.destroy()
            
        exit_window = tk.Toplevel()
        exit_window.attributes("-topmost", True)
        exit_window.withdraw()
        exit_window.title("Exit")
        icon = ImageTk.PhotoImage(ico_image)
        exit_window.iconphoto(False, icon)
        exit_window.resizable(False, False)
        exit_window.grab_set()
        exit_window.protocol("WM_DELETE_WINDOW", on_closing)
        
        exit_window_label = ttk.Label(exit_window, text="Are you sure you want to exit?", font=("Segoe UI", list_font_size), justify="center")
        exit_window_label.grid(row=0, column=0, columnspan=2, padx=50, pady=20)
        
        ok_exit_window_button = ttk.Button(exit_window, text="Ok", style="small_button_sfw_style1.TButton", command=ok_close_main_window)
        ok_exit_window_button.grid(row=1, column=0, padx=10, pady=10, sticky="e")
        
        cancel_exit_window_button = ttk.Button(exit_window, text="Cancel", style="small_button_sfw_style1.TButton", command=close_exit_window)
        cancel_exit_window_button.grid(row=1, column=1, padx=10, pady=10, sticky="w")
        
        exit_window.update_idletasks()
        exit_window.deiconify()
        window_width = exit_window.winfo_width()
        window_height = exit_window.winfo_height()
        screen_width = exit_window.winfo_screenwidth()
        screen_height = exit_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        exit_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        if direct_close:
            ok_close_main_window()
    
    def handle_double_left_click(event): #peak visualizer
        selected_item_pv = chromatograms_list.selection()
        if len(selected_item_pv) == 0 or len(selected_item_pv) > 1:
            return
        else:
            level = 0
            parent_item = selected_item_pv
            while parent_item:
                level += 1
                parent_item = chromatograms_list.parent(parent_item)
            if level != 3:
                return
        
        def close_peak_visualizer():
            peak_visualizer.destroy()
            
        def on_pick_pv(event):
            # Get the artist (marker) that was picked
            artist = event.artist
            # Get the index of the picked point
            index = event.ind[0]
            # Get the x and y coordinates of the picked point
            x = artist.get_xdata()[index]
            y = artist.get_ydata()[index]
            
            artists = ax_pv.get_children()
            for artist in artists:
                if isinstance(artist, matplotlib.lines.Line2D) and artist.get_markersize() == 8:
                    artist.remove()
                    break
            marker_pv = ax_pv.plot(x, y, marker='o', markersize=8, color='black')
            canvas_pv.draw()
            
            clear_plot(ax_if, canvas_if)
            
            parent_text_split_zero = parent_text.split(" ")[0]
            if x in isotopic_fittings[sample_index][f"{grand_parent_text}_{parent_text_split_zero}"]:
                x_values_if_ideal = isotopic_fittings[sample_index][f"{grand_parent_text}_{parent_text_split_zero}"][x][0]
                y_values_if_ideal = isotopic_fittings[sample_index][f"{grand_parent_text}_{parent_text_split_zero}"][x][1]
                if len(x_values_if_ideal) > 0:
                    ax_if.set_xlim(x_values_if_ideal[0]-0.2, x_values_if_ideal[-1]+0.2)
                    x_values_if_actual = [x-(ax_if.get_xlim()[1]-ax_if.get_xlim()[0])*0.02 for x in x_values_if_ideal]
                    y_values_if_actual = isotopic_fittings[sample_index][f"{grand_parent_text}_{parent_text_split_zero}"][x][2]
                    ax_if.set_ylim(0, max(y_values_if_actual)*1.1 if max(y_values_if_actual) > max(y_values_if_ideal) else max(y_values_if_ideal)*1.1)
                    ax_if.plot(x_values_if_ideal, y_values_if_ideal, marker='', linewidth=0, label="Ideal")
                    ax_if.plot(x_values_if_actual, y_values_if_actual, marker='', linewidth=0, label="Found")
                    ax_if.vlines(x_values_if_ideal, ymin=0, ymax=y_values_if_ideal, linewidth=3, colors='red')
                    ax_if.vlines(x_values_if_actual, ymin=0, ymax=y_values_if_actual, linewidth=3, colors='blue')
                    label_rt = float("%.2f" % round(x, 2))
                    label_score = float("%.3f" % round(isotopic_fittings[sample_index][f"{grand_parent_text}_{parent_text_split_zero}"][x][3], 3))
                    info_label.config(text=f"RT: {label_rt}    Score: {label_score}")
                else:
                    label_rt = float("%.2f" % round(x, 2))
                    info_label.config(text=f"RT: {label_rt}    Score: 0.0")
                    ax_if.text(0.5, 0.5, 'No data available for this datapoint.\nThis might be caused by very poor data\nand may account for potentially low\nisotopic fitting score for the peak.', horizontalalignment='center', verticalalignment='center', transform=ax_if.transAxes, fontsize=10, wrap=True)
            else:
                label_rt = float("%.2f" % round(x, 2))
                info_label.config(text=f"RT: {label_rt}    Score: 0.0")
                ax_if.text(0.5, 0.5, 'No data available for this datapoint.\nThis might be caused by very poor data\nand may account for potentially low\nisotopic fitting score for the peak.', horizontalalignment='center', verticalalignment='center', transform=ax_if.transAxes, fontsize=10, wrap=True)
            canvas_if.draw()

        def on_hover_pv(event):
            # Change cursor to a hand when hovering over markers
            if event.inaxes:
                for artist in event.inaxes.get_children():
                    if isinstance(artist, plt.Line2D) and artist.contains(event)[0]:
                        event.canvas.get_tk_widget().configure(cursor="hand2")
                        return
            event.canvas.get_tk_widget().configure(cursor="")        
    
        peak_visualizer = tk.Toplevel()
        # peak_visualizer.attributes("-topmost", True)
        peak_visualizer.geometry("700x700")
        peak_visualizer.grid_columnconfigure(0, weight=1)
        peak_visualizer.grid_columnconfigure(1, weight=1)
        peak_visualizer.grid_rowconfigure(0, weight=1)
        peak_visualizer.grid_rowconfigure(1, weight=1)
        peak_visualizer.withdraw()
        icon = ImageTk.PhotoImage(ico_image)
        peak_visualizer.iconphoto(False, icon)
        peak_visualizer.resizable(False, False)
        
        peak_area = ttk.Labelframe(peak_visualizer, text="Peak Visualization", style="chromatogram.TLabelframe")
        peak_area.grid(row=0, rowspan=2, column=0, columnspan=2, padx=10, pady=(10, 300), sticky="nsew")
        
        isofit_area = ttk.Labelframe(peak_visualizer, text="Isotopic Peaks", style="chromatogram.TLabelframe")
        isofit_area.grid(row=1, column=0, columnspan=2, padx=(10, 280), pady=(310, 10), sticky="nsew")
        
        info_area = ttk.Labelframe(peak_visualizer, text="Peak Information", style="chromatogram.TLabelframe")
        info_area.grid(row=1, column=0, columnspan=2, padx=(430, 10), pady=(310, 40), sticky="nsew")
        info_area.propagate(0)
        
        close_peak_visualizer_button = ttk.Button(peak_visualizer, text="Close", style="small_button_sfw_style1.TButton", command=close_peak_visualizer)
        close_peak_visualizer_button.grid(row=1, column=0, columnspan=2, padx=10, pady=10, sticky="se")
        
        fig_pv = plt.figure(figsize=(0, 0))
        ax_pv = fig_pv.add_subplot(111)
        canvas_pv = FigureCanvasTkAgg(fig_pv, master=peak_area)
        canvas_pv.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        ax_pv.set_position([0.12, 0.12, 0.86, 0.82])
        ax_pv.set_xlabel('Retention/Migration Time (min)')
        ax_pv.set_ylabel('Intensity (AU)')
        
        sample_index = list(glycans_per_sample.keys()).index(selected_item)
        parent_item = chromatograms_list.parent(chromatograms_list.selection())
        parent_text = chromatograms_list.item(parent_item, "text") #adduct
        grand_parent_item = chromatograms_list.parent(parent_item)
        grand_parent_text = chromatograms_list.item(grand_parent_item, "text")
        parent_text_split_zero = parent_text.split(" ")[0]
        chromatograms_list_selected_item_text = chromatograms_list.item(chromatograms_list.selection(), "text")
        if f"{grand_parent_text}+{parent_text_split_zero}_{chromatograms_list_selected_item_text}_RTs" in curve_fittings[sample_index]:
            y_values_ideal = curve_fittings[sample_index][f"{grand_parent_text}+{parent_text_split_zero}_{chromatograms_list_selected_item_text}_Ideal_ints"]
            y_values_found = curve_fittings[sample_index][f"{grand_parent_text}+{parent_text_split_zero}_{chromatograms_list_selected_item_text}_Found_ints"]
            x_values_pv = curve_fittings[sample_index][f"{grand_parent_text}+{parent_text_split_zero}_{chromatograms_list_selected_item_text}_RTs"]
        peak_visualizer.title(f"Peak Visualizer - {grand_parent_text} - {parent_text}")
            
        ax_pv.plot(x_values_pv, y_values_ideal, marker="o", color="red", label = 'Ideal', picker=5)
        ax_pv.plot(x_values_pv, y_values_found, marker="o", color="blue", label = 'Found', picker=5)
        
        ax_pv.legend()
        
        canvas_pv.mpl_connect('pick_event', on_pick_pv)
        canvas_pv.mpl_connect('motion_notify_event', on_hover_pv)
        canvas_pv.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_pv, canvas_pv, False) if event.button == 3 else None)
        
        canvas_pv.draw()
            
        fig_if = plt.figure(figsize=(0, 0))
        ax_if = fig_if.add_subplot(111)
        canvas_if = FigureCanvasTkAgg(fig_if, master=isofit_area)
        canvas_if.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        ax_if.set_position([0.15, 0.18, 0.80, 0.74])
        ax_if.set_xlabel('m/z')
        ax_if.set_ylabel('Relative Intensity')
        
        info_label = tk.Label(isofit_area, text="", anchor="w", font=("Segoe UI", 8), bg="white")
        info_label.place(relx=1.0, rely=0, anchor='ne')
        info_label.lift()
        
        canvas_if.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_if, canvas_if, False) if event.button == 3 else None)
        
        canvas_if.draw()
        
        parent_text_split_zero = parent_text.split(" ")[0]
        chromatogram_parent_item_adduct = chromatograms_list.item(parent_item, "text").split(" ")[0]
        chromatograms_list_selected_item_text = chromatograms_list.item(chromatograms_list.selection(), "text")
        chromatograms_list_grand_parent_text = chromatograms_list.item(grand_parent_item, "text")
        
        iso_fitting_score = glycans_per_sample[selected_item][grand_parent_text][parent_text_split_zero]['iso'][glycans_per_sample[selected_item][chromatograms_list_grand_parent_text][chromatogram_parent_item_adduct]['peaks'].index(chromatograms_list_selected_item_text)]
        
        curve_fitting_score = glycans_per_sample[selected_item][grand_parent_text][parent_text_split_zero]['curve'][glycans_per_sample[selected_item][chromatograms_list_grand_parent_text][chromatogram_parent_item_adduct]['peaks'].index(chromatograms_list_selected_item_text)]
        
        s_to_n_score = glycans_per_sample[selected_item][grand_parent_text][parent_text_split_zero]['sn'][glycans_per_sample[selected_item][chromatograms_list_grand_parent_text][chromatogram_parent_item_adduct]['peaks'].index(chromatograms_list_selected_item_text)]
        
        ppm_score = glycans_per_sample[selected_item][grand_parent_text][parent_text_split_zero]['ppm'][glycans_per_sample[selected_item][chromatograms_list_grand_parent_text][chromatogram_parent_item_adduct]['peaks'].index(chromatograms_list_selected_item_text)]
        
        auc_for_label = f"{glycans_per_sample[selected_item][grand_parent_text][parent_text_split_zero]['auc'][glycans_per_sample[selected_item][chromatograms_list_grand_parent_text][chromatogram_parent_item_adduct]['peaks'].index(chromatograms_list_selected_item_text)]:.1e}"
        
        ambiguities = f"{glycans_per_sample[selected_item][grand_parent_text][parent_text_split_zero]['ambiguity']}"
        
        peak_info = [("Isotopic Fitting Score:", iso_fitting_score),
                     ("Curve Fitting Score:", curve_fitting_score),
                     ("Signal-to-Noise ratio:", s_to_n_score),
                     ("Average PPM error:", ppm_score),
                     ("Area Under Curve (AUC):", auc_for_label),
                     ("MS2 TIC explained:", f"No MS2 Data"),
                     ("Ambiguities:", ambiguities)]
        
        if 'ms2' in glycans_per_sample[selected_item][grand_parent_text][parent_text.split(" ")[0]]:
            highest_tic_explained = 0
            for i in list(glycans_per_sample[selected_item][grand_parent_text][parent_text.split(" ")[0]]['ms2'].keys()):
                if abs(float(i)-float(chromatograms_list.item(chromatograms_list.selection(), "text"))) < 1: #change number from 1?
                    current_tic_explained = (sum(glycans_per_sample[selected_item][grand_parent_text][parent_text.split(" ")[0]]['ms2'][i][1])/glycans_per_sample[selected_item][grand_parent_text][parent_text.split(" ")[0]]['ms2'][i][3])*100
                    if current_tic_explained > highest_tic_explained:
                        highest_tic_explained = current_tic_explained
            if highest_tic_explained != 0:
                ms2_tic_formatted = float("%.1f" % round(highest_tic_explained, 1))
                peak_info = [("Isotopic Fitting Score:", iso_fitting_score),
                             ("Curve Fitting Score:", curve_fitting_score),
                             ("Signal-to-Noise ratio:", s_to_n_score),
                             ("Average PPM error:", ppm_score),
                             ("Area Under Curve (AUC):", auc_for_label),
                             ("MS2 TIC explained:", f"{ms2_tic_formatted}%"),
                             ("Ambiguities:", ambiguities)]
        
        for i, (left_text, right_text) in enumerate(peak_info):
            left_label = ttk.Label(info_area, text=left_text, anchor="w", font=("Segoe UI", 10))
            right_label = ttk.Label(info_area, text=right_text, anchor="e", font=("Segoe UI", 10))
            
            left_label.grid(row=i, column=0, sticky="w")
            right_label.grid(row=i, column=1, sticky="e")
            
        info_area.columnconfigure(0, weight=1)
        info_area.columnconfigure(1, weight=1)
        
        peak_visualizer.update_idletasks()
        peak_visualizer.deiconify()
        window_width = peak_visualizer.winfo_width()
        window_height = peak_visualizer.winfo_height()
        screen_width = peak_visualizer.winfo_screenwidth()
        screen_height = peak_visualizer.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        peak_visualizer.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")

    def click_treeview(event):
        global chromatograms_list, selected_item, samples_list, selected_chromatograms, last_plotted
        
        if gg_draw_on:
            check_gg_drawings_available(chromatograms_list, False)
        
        if event.keysym == "Escape": #removes selection and clears zoom and rectangles
            selected_items = chromatograms_list.selection()
            last_plotted = None
            if len(selected_items) > 0:
                chromatograms_list.selection_remove(selected_items)
                plot_graph_button.config(state=tk.NORMAL)
                compare_samples_button.config(state=tk.DISABLED)
                clear_plot(ax, canvas)
            return 
            
        region = chromatograms_list.identify_region(event.x, event.y)
        column = chromatograms_list.identify_column(event.x)
        item = chromatograms_list.identify_row(event.y)
        
        if region == "cell" and column == "#1" and item and len(samples_list) > 0 and selected_item in processed_data:
            values = chromatograms_list.item(item, "values")
            if "MS2" in values:
                handle_treeview_select(event)
                parent_item = chromatograms_list.parent(chromatograms_list.selection())
                parent_text = chromatograms_list.item(parent_item, "text").split(" ") #adduct
                grand_parent_item = chromatograms_list.parent(parent_item)
                grand_parent_text = chromatograms_list.item(grand_parent_item, "text")
                distance = 99999
                rt_closest = 0
                rt = float(chromatograms_list.item(chromatograms_list.selection(), "text"))
                for i in glycans_per_sample[selected_item][grand_parent_text][parent_text[0]]['ms2']:
                    if abs(float(i)-rt) < distance:
                        distance = abs(float(i)-rt)
                        rt_closest = float(i)
                if rt_closest != 0:
                    if current_data['time_unit'] == 'seconds':
                        rt_seconds = float("%.4f" % round(rt_closest*60, 4))
                    else:
                        rt_seconds = float("%.4f" % round(rt_closest, 4))
                    spectra_info = (float(parent_text[-1]), rt_seconds)
                    show_ms2_graph(spectra_info)
            else:
                if len(chromatograms_list.selection()) > 1:
                    handle_treeview_select(event, False)
                    canvas.draw()
                else:
                    handle_treeview_select(event)
        else:
            
            if key_states.get("Control", False) or key_states.get("Control_L", False) or key_states.get("Control_R", False):
                glycan = chromatograms_list.item(item, "text")
                parent_item = chromatograms_list.parent(item)
                parent_text = chromatograms_list.item(parent_item, "text")
                
                if type(glycan) == float:
                    glycan = parent_text
                    parent_item = chromatograms_list.parent(parent_item)
                    parent_text = chromatograms_list.item(parent_item, "text")
                
                remove = False
                for index, element in enumerate(selected_chromatograms):
                    if glycan == element[0] or f"{parent_text}+{glycan}" == element[0]:
                        remove = True
                        del selected_chromatograms[index]
                        break
                        
                if not remove:
                    if (len(glycan.split("-")) > 1 and not f"{parent_text}+{glycan}".startswith("+")) or type(glycan) == float:
                        selected_chromatograms.append([f"{parent_text}+{glycan}", 2])
                    else:
                        selected_chromatograms.append([glycan, 1])
                            
            elif key_states.get("Shift", False) or key_states.get("Shift_L", False) or key_states.get("Shift_R", False):
                for chromatogram in chromatograms_list.selection():
                    glycan = chromatograms_list.item(chromatogram, "text")
                            
                    if type(glycan) == float or len(glycan.split("-")) > 1:
                        continue
                        
                    selected_chromatograms.append([glycan, 1])
                            
            else:
                selected_chromatograms = []
            
            if len(selected_chromatograms) > 0:
                handle_treeview_select(event, False)
                canvas.draw()
            else:
                handle_treeview_select(event)
            
    def on_treeview_motion(event):
        region = chromatograms_list.identify_region(event.x, event.y)
        column = chromatograms_list.identify_column(event.x)
        item = chromatograms_list.identify_row(event.y)
        
        if region == "cell" and column == "#1" and item and len(samples_list) > 0 and selected_item in processed_data:
            values = chromatograms_list.item(item, "values")
            # Change the cursor only if the cell value is "Value 5"
            if "MS2" in values:
                chromatograms_list.config(cursor="hand2")
            else:
                chromatograms_list.config(cursor="")
        else:
            chromatograms_list.config(cursor="")
        
    def backend_heatmap():
        global current_data, ax, ax_spec, two_d

        def prepare_data_2dplot():
            rt_list = []
            mz_list = []
            intensity_list = []
            max_int = 0
            min_int = 99999
            
            rt_list_ms1 = []
            mz_list_ms1 = []
            intensity_list_ms1 = []
            max_int_ms1 = 0
            min_int_ms1 = 99999
            
            rt_list_ms2 = []
            mz_list_ms2 = []
            intensity_list_ms2 = []
            max_int_ms2 = 0
            min_int_ms2 = 99999
                    
            if current_data['file_type'] == 'mzml':
                # Load the mzML file using pyteomics
                with mzml.read(current_data['file_path']) as reader:
                    for spectrum in reader:
                    
                        # Extract retention time from the scanList
                        scan_list = spectrum.get('scanList')
                        scan_time = scan_list['scan'][0]['scan start time']
                        
                        if current_data['time_unit'] == 'seconds':
                            rt = float("%.2f" % round(scan_time/60, 2))
                        else:
                            rt = scan_time
                        
                        mz_array = spectrum['m/z array']
                        intensity_array = spectrum['intensity array']
                        
                        if len(intensity_array) == 0:
                            continue
                        max_array = max(intensity_array)
                        min_array = min(intensity_array)
                        
                        # For all spectra
                        rt_list.extend([rt] * len(mz_array))
                        mz_list.extend(mz_array)
                        intensity_list.extend(intensity_array)
                        if max_array > max_int:
                            max_int = max_array
                        if min_array < min_int and min_array != 0:
                            min_int = min_array
                                
                        # For MS1 spectra only 
                        if spectrum.get('ms level') == 1:
                            rt_list_ms1.extend([rt] * len(mz_array))
                            mz_list_ms1.extend(mz_array)
                            intensity_list_ms1.extend(intensity_array)
                            if max_array > max_int_ms1:
                                max_int_ms1 = max_array
                            if min_array < min_int_ms1 and min_array != 0:
                                min_int_ms1 = min_array
                                    
                        # For MS2 spectra only
                        if spectrum.get('ms level') == 2:
                            rt_list_ms2.extend([rt] * len(mz_array))
                            mz_list_ms2.extend(mz_array)
                            intensity_list_ms2.extend(intensity_array)
                            if max_array > max_int_ms2:
                                max_int_ms2 = max_array
                            if min_array < min_int_ms2 and min_array != 0:
                                min_int_ms2 = min_array
            else:
                # Load the mzXML file using pyteomics
                with mzxml.read(current_data['file_path']) as reader:
                    for spectrum in reader:
                        
                        # Get retention time
                        scan_time = spectrum.get('retentionTime')
                        
                        if current_data['time_unit'] == 'seconds':
                            rt = float("%.2f" % round(scan_time/60, 2))
                        else:
                            rt = scan_time
                        
                        mz_array = spectrum['m/z array']
                        intensity_array = spectrum['intensity array']
                        
                        if len(mz_array) == 0:
                            continue
                        max_array = max(intensity_array)
                        min_array = min(intensity_array)
                        
                        rt_list.extend([rt] * len(mz_array))
                        mz_list.extend(mz_array)
                        intensity_list.extend(intensity_array)
                        if max_array > max_int:
                            max_int = max_array
                        if min_array < min_int and min_array != 0:
                            min_int = min_array
                                
                        # For MS1 spectra only 
                        if spectrum.get('msLevel') == 1:
                            rt_list_ms1.extend([rt] * len(mz_array))
                            mz_list_ms1.extend(mz_array)
                            intensity_list_ms1.extend(intensity_array)
                            if max_array > max_int_ms1:
                                max_int_ms1 = max_array
                            if min_array < min_int_ms1 and min_array != 0:
                                min_int_ms1 = min_array
                                    
                        # For MS2 spectra only
                        if spectrum.get('msLevel') == 2:
                            rt_list_ms2.extend([rt] * len(mz_array))
                            mz_list_ms2.extend(mz_array)
                            intensity_list_ms2.extend(intensity_array)
                            if max_array > max_int_ms2:
                                max_int_ms2 = max_array
                            if min_array < min_int_ms2 and min_array != 0:
                                min_int_ms2 = min_array
            
            if len(rt_list_ms2) > 0:
                return (np.array(rt_list), np.array(mz_list), np.array(intensity_list), min_int, max_int), (np.array(rt_list_ms1), np.array(mz_list_ms1), np.array(intensity_list_ms1), min_int_ms1, max_int_ms1), (np.array(rt_list_ms2), np.array(mz_list_ms2), np.array(intensity_list_ms2), min_int_ms2, max_int_ms2)
            else:
                return [(np.array(rt_list), np.array(mz_list), np.array(intensity_list), min_int, max_int)]
        
        def plot_with_scatter_density(rt_list, mz_list, intensity_list, min_int, max_int):
            fig_sd, ax_sd = plt.subplots(figsize=(1, 1), subplot_kw={'projection': 'scatter_density'})
            
            # Plot using density
            density = ax_sd.scatter_density(rt_list, mz_list, c=intensity_list, dpi=18, cmap='inferno_r', norm=LogNorm(min_int, max_int), downres_factor = 1, picker=True)

            # Add colorbar
            fig_sd.colorbar(density, label='Intensity', norm=LogNorm(min_int, max_int))

            ax_sd.set_xlabel("Retention/Migration Time (minutes)")
            ax_sd.set_ylabel("m/z")

            return fig_sd, ax_sd
            
        def plot_2d_visualization(rt_list, mz_list, intensity_list, min_int, max_int, text_notebook):

            # Create a frame for the plot and toolbar
            frame_two_d = tk.Frame(notebook_two_d)
            frame_two_d.pack(fill=tk.BOTH, expand=True)
            
            # Add the frames to notebook
            notebook_two_d.add(frame_two_d, text=text_notebook) 
            
            # Set minimum intensity to 1 if it's zero
            if min_int <= 0:
                min_int = 1
                
            # Plot the data
            fig_two_d, ax_two_d = plot_with_scatter_density(rt_list, mz_list, intensity_list, min_int, max_int)
    
            canvas_two_d = FigureCanvasTkAgg(fig_two_d, master=frame_two_d)
            canvas_two_d.get_tk_widget().pack(fill=tk.BOTH, expand=True)

            # Add the navigation toolbar
            toolbar = CustomToolbar(canvas_two_d, frame_two_d)
            toolbar.update()
            toolbar.pack(side=tk.BOTTOM, fill=tk.X)
                
            if text_notebook == 'All Spectra' or text_notebook == 'MS2 Spectra':
                lines_collection = []
                scaling_factor = 0.001
                if len(processed_data[selected_item]['ms2']) > 0:
                    for i in processed_data[selected_item]['ms2']:
                        if processed_data[selected_item]['time_unit'] == 'seconds':
                            i_mins = i/60
                        else:
                            i_mins = i
                        for k in processed_data[selected_item]['ms2'][i]:
                            if processed_data[selected_item]['time_unit'] == 'seconds':
                                k_mins = k/60
                            else:
                                k_mins = k
                            if text_notebook == 'All Spectra':
                                precursor_coordinates = [i_mins, processed_data[selected_item]['ms2'][i][k][0]]
                                
                                # Horizontal line of cross
                                horizontal_line = ax_two_d.plot([precursor_coordinates[0]-(ax_two_d.get_xlim()[1]*scaling_factor), k_mins], [precursor_coordinates[1], precursor_coordinates[1]], color='blue')
                                
                                # Vertical line of cross
                                vertical_line = ax_two_d.plot([precursor_coordinates[0], precursor_coordinates[0]], [precursor_coordinates[1]-(ax_two_d.get_ylim()[1]*scaling_factor), precursor_coordinates[1]+(ax_two_d.get_ylim()[1]*scaling_factor)], color='blue')
                                
                                cross = [horizontal_line, vertical_line]
                                
                                lines_collection.append(cross)
                            elif text_notebook == 'MS2 Spectra':
                                precursor_coordinates = [k_mins, processed_data[selected_item]['ms2'][i][k][0]]
                                
                                # Horizontal line of cross
                                horizontal_line = ax_two_d.plot([precursor_coordinates[0]-(ax_two_d.get_xlim()[1]*scaling_factor), precursor_coordinates[0]+(ax_two_d.get_xlim()[1]*scaling_factor)], [precursor_coordinates[1], precursor_coordinates[1]], color='blue')
                                
                                # Vertical line of cross
                                vertical_line = ax_two_d.plot([precursor_coordinates[0], precursor_coordinates[0]], [precursor_coordinates[1]-(ax_two_d.get_ylim()[1]*scaling_factor), precursor_coordinates[1]+(ax_two_d.get_ylim()[1]*scaling_factor)], color='blue') 
                                
                                cross = [horizontal_line, vertical_line]
                                
                                lines_collection.append(cross)
            
            canvas_two_d.draw()
            
        def process_2d_plot():
            def calculate_2d_plot():    
                # Process sample
                processed_2d_data = prepare_data_2dplot()
                
                # Process sample
                for i_i, i in enumerate(processed_2d_data):
                    rt_list, mz_list, intensity_list, min_int, max_int = i
                    if i_i == 0 and len(processed_2d_data) > 1:
                        plot_2d_visualization(rt_list, mz_list, intensity_list, min_int, max_int, 'All Spectra')
                    elif i_i == 1 or len(processed_2d_data) == 1:
                        plot_2d_visualization(rt_list, mz_list, intensity_list, min_int, max_int, 'MS1 Spectra')
                    elif i_i == 2:
                        plot_2d_visualization(rt_list, mz_list, intensity_list, min_int, max_int, 'MS2 Spectra')
                        
                processing_2d_plot.destroy()
                two_d_plot.deiconify()
            
            def wait_thread():
                calculate_2d_plot()
            
            global processing_2d_plot
            processing_2d_plot = tk.Toplevel()
            processing_2d_plot.withdraw()
            processing_2d_plot.title("Processing 2D Plot")
            icon = ImageTk.PhotoImage(ico_image)
            processing_2d_plot.iconphoto(False, icon)
            processing_2d_plot.resizable(False, False)
            processing_2d_plot.grab_set()
            processing_2d_plot.protocol("WM_DELETE_WINDOW", on_closing)
            
            processing_max_spectrum_label = ttk.Label(processing_2d_plot, text="Processing 2D plot, please wait.", font=("Segoe UI", list_font_size))
            processing_max_spectrum_label.pack(pady=35, padx=70)
            
            processing_2d_plot.update_idletasks()
            processing_2d_plot.deiconify()
            window_width = processing_2d_plot.winfo_width()
            window_height = processing_2d_plot.winfo_height()
            screen_width = processing_2d_plot.winfo_screenwidth()
            screen_height = processing_2d_plot.winfo_screenheight()
            x_position = (screen_width - window_width) // 2
            y_position = (screen_height - window_height) // 2
            processing_2d_plot.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
            
            t = threading.Thread(target=wait_thread)
            t.start()
        
        global two_d_plot
        two_d_plot = tk.Toplevel()
        icon = ImageTk.PhotoImage(ico_image)
        two_d_plot.iconphoto(False, icon)
        two_d_plot.withdraw()
        two_d_plot.minsize(450, 400)
        two_d_plot.title("2D-Plot")
        two_d_plot.resizable(True, True)
        
        # Create a notebook (tabs container)
        global notebook_two_d
        notebook_two_d = ttk.Notebook(two_d_plot)
        notebook_two_d.pack(fill=tk.BOTH, expand=True)
        
        two_d_plot.update_idletasks()
        window_width = two_d_plot.winfo_width()
        window_height = two_d_plot.winfo_height()
        screen_width = two_d_plot.winfo_screenwidth()
        screen_height = two_d_plot.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        two_d_plot.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        process_2d_plot()
        
    def compare_samples_window():
        global compare_samples, compare_samples_opened, samples_dropdown_options, level, ret_time_interval, df1, df2, aligning_samples
        
        def aligning_samples_window():
            def process_alignment():
                global former_alignments, ax_comp, canvas_comp
                
                if f"{iso_fit_score}_{curve_fit_score}_{max_ppm}_{s_to_n}" not in former_alignments:
                    total_glycans_df = make_total_glycans_df(df1, df2)
                    if len(total_glycans_df[0]['Glycan']) != 0:
                        chromatograms_delta = Execution_Functions.align_assignments(total_glycans_df, "total_glycans", multithreaded_analysis, number_cores, rt_tol = ret_time_interval[2])
                        Execution_Functions.align_assignments(df2, 'chromatograms', multithreaded_analysis, number_cores, temp_folder, reanalysis_path, chromatograms_delta[1], None, iso_fit_score, curve_fit_score, max_ppm, s_to_n)
                    former_alignments.append(f"{iso_fit_score}_{curve_fit_score}_{max_ppm}_{s_to_n}")
                
                last_xlims_chrom = ax_comp.get_xlim()
                last_ylims_chrom = ax_comp.get_ylim()
                clear_plot(ax_comp, canvas_comp)
                vertical_line = ax_comp.axvline(x=-10000, color='black', linestyle='--', linewidth=1)
                
                loops = 0
                if level == 1:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        with open(os.path.join(temp_folder, f"{i}_aligned_RTs_{iso_fit_score}_{curve_fit_score}_{max_ppm}_{s_to_n}"), 'rb') as f:
                            x_values_comp = dill.load(f)
                            f.close()
                            
                        y_values_comp = []
                        counter = 0
                        for k in gg_file.list_chromatograms(i):
                            if "+".join(k.split("+")[:-1]) == item_text_comp:
                            
                                # Load intensities from file
                                y_values_temp = gg_file.get_chromatogram(i, k, 'smoothed')
                                    
                                if counter == 0:
                                    y_values_comp = y_values_temp
                                    counter+= 1
                                else:
                                    y_values_comp = [x + y for x, y in zip(y_values_comp, y_values_temp)]
                                    
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
                if level == 2:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        with open(os.path.join(temp_folder, f"{i}_aligned_RTs_{iso_fit_score}_{curve_fit_score}_{max_ppm}_{s_to_n}"), 'rb') as f:
                            x_values_comp = dill.load(f)
                            f.close()
                            
                        # Load intensities from file
                        y_values_comp = gg_file.get_chromatogram(i, f"{grand_parent_text_comp}+{parent_text_comp}", 'smoothed')
                            
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
                    
                ax_comp.legend(fontsize=9)
                    
                lines = ax_comp.get_lines()

                # Initialize the maximum value to negative infinity
                max_value = float('-inf')

                # Iterate through each line
                for line in lines:
                    # Get the y-data of the line
                    y_data = line.get_ydata()
                    
                    # Find the maximum value in the y-data
                    line_max_value = max(y_data)
                    
                    # Update the overall maximum value if needed
                    max_value = max(max_value, line_max_value)
                
                ax_comp.set_xlim(last_xlims_chrom)
                ax_comp.set_ylim(0, max_value*1.1)
                ax_comp.set_xlabel('Retention/Migration Time (min)')
                ax_comp.set_ylabel('Intensity (AU)')
                
                canvas_comp.draw_idle()
                
                aligning_samples.destroy()
            
            global aligning_samples
            aligning_samples = tk.Toplevel()
            # aligning_samples.attributes("-topmost", True)
            aligning_samples.withdraw()
            aligning_samples.title("Aligning Samples")
            icon = ImageTk.PhotoImage(ico_image)
            aligning_samples.iconphoto(False, icon)
            aligning_samples.resizable(False, False)
            aligning_samples.grab_set()
            aligning_samples.protocol("WM_DELETE_WINDOW", on_closing)
            
            aligning_samples_label = ttk.Label(aligning_samples, text="Aligning sample files, please wait.", font=("Segoe UI", list_font_size))
            aligning_samples_label.pack(pady=35, padx=70)
            
            aligning_samples.update_idletasks()
            aligning_samples.deiconify()
            window_width = aligning_samples.winfo_width()
            window_height = aligning_samples.winfo_height()
            screen_width = aligning_samples.winfo_screenwidth()
            screen_height = aligning_samples.winfo_screenheight()
            x_position = (screen_width - window_width) // 2
            y_position = (screen_height - window_height) // 2
            aligning_samples.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
            
            t = threading.Thread(target=process_alignment)
            t.start()
        
        def on_checkbox_click(event):
            """Check or uncheck box when clicked."""            
            x, y, widget = event.x, event.y, event.widget
            elem = widget.identify("element", x, y)
            if "image" in elem:
                # a box was clicked
                item = widget.identify_row(y)
                if widget.tag_has("unchecked", item) or widget.tag_has("tristate", item):
                    widget._check_ancestor(item)
                    widget._check_descendant(item)
                else:
                    widget._uncheck_descendant(item)
                    widget._uncheck_ancestor(item)
                    
            draw_checked_samples()
                    
        def draw_checked_samples():
            
            last_xlims_chrom = ax_comp.get_xlim()
            last_ylims_chrom = ax_comp.get_ylim()
            
            clear_plot(ax_comp, canvas_comp)
            
            vertical_line = ax_comp.axvline(x=-10000, color='black', linestyle='--', linewidth=1)
            
            loops = 0
            if align_chromatograms_checkbox_state.get():
                if level == 1:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        with open(os.path.join(temp_folder, f"{i}_aligned_RTs_{iso_fit_score}_{curve_fit_score}_{max_ppm}_{s_to_n}"), 'rb') as f:
                            x_values_comp = dill.load(f)
                            f.close()
                            
                        y_values_comp = []
                        counter = 0
                        for k in gg_file.list_chromatograms(i):
                            if "+".join(k.split("+")[:-1]) == item_text_comp:
                            
                                # Load intensities from file
                                y_values_temp = gg_file.get_chromatogram(i, k, 'smoothed')
                                    
                                if counter == 0:
                                    y_values_comp = y_values_temp
                                    counter+= 1
                                else:
                                    y_values_comp = [x + y for x, y in zip(y_values_comp, y_values_temp)]
                                    
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
                if level == 2:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        with open(os.path.join(temp_folder, f"{i}_aligned_RTs_{iso_fit_score}_{curve_fit_score}_{max_ppm}_{s_to_n}"), 'rb') as f:
                            x_values_comp = dill.load(f)
                            f.close()
                            
                        # Load intensities from file
                        y_values_comp = gg_file.get_chromatogram(i, f"{grand_parent_text_comp}+{parent_text_comp}", 'smoothed')
                            
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
            else:
                if level == 1:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        x_values_comp = gg_file.get_rt_array(i)
                            
                        y_values_comp = []
                        counter = 0
                        for k in gg_file.list_chromatograms(i):
                            if "+".join(k.split("+")[:-1]) == item_text_comp:
                            
                                # Load intensities from file
                                y_values_temp = gg_file.get_chromatogram(i, k, 'smoothed')
                                    
                                if counter == 0:
                                    y_values_comp = y_values_temp
                                    counter+= 1
                                else:
                                    y_values_comp = [x + y for x, y in zip(y_values_comp, y_values_temp)]
                                    
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
                if level == 2:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        x_values_comp = gg_file.get_rt_array(i)
                            
                        # Load intensities from file
                        y_values_comp = gg_file.get_chromatogram(i, f"{grand_parent_text_comp}+{parent_text_comp}", 'smoothed')
                        
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
                    
            ax_comp.legend(fontsize=9)
            
            ax_comp.set_xlim(last_xlims_chrom)
            ax_comp.set_ylim(last_ylims_chrom)
            ax_comp.set_xlabel('Retention/Migration Time (min)')
            ax_comp.set_ylabel('Intensity (AU)')
            
            canvas_comp.draw_idle()
        
        def align_chromatograms_checkbox_state_check():
            state = align_chromatograms_checkbox_state.get()
            
            if state:
                aligning_samples_window()
            else:
                last_xlims_chrom = ax_comp.get_xlim()
                last_ylims_chrom = ax_comp.get_ylim()
                
                clear_plot(ax_comp, canvas_comp)
                vertical_line = ax_comp.axvline(x=-10000, color='black', linestyle='--', linewidth=1)
                
                loops = 0
                if level == 1:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        x_values_comp = gg_file.get_rt_array(i)
                            
                        y_values_comp = []
                        counter = 0
                        for k in gg_file.list_chromatograms(i):
                            if "+".join(k.split("+")[:-1]) == item_text_comp:
                            
                                # Load intensities from file
                                y_values_temp = gg_file.get_chromatogram(i, k, 'smoothed')
                                    
                                if counter == 0:
                                    y_values_comp = y_values_temp
                                    counter+= 1
                                else:
                                    y_values_comp = [x + y for x, y in zip(y_values_comp, y_values_temp)]
                                    
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
                if level == 2:
                    for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
                
                        # Load retention times from file
                        x_values_comp = gg_file.get_rt_array(i)
                            
                        # Load intensities from file
                        y_values_comp = gg_file.get_chromatogram(i, f"{grand_parent_text_comp}+{parent_text_comp}", 'smoothed')
                            
                        if i_i-(len(colors)*loops) == len(colors):
                            loops+= 1
                        color = colors[i_i-(len(colors)*loops)]
                        ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
                    
                ax_comp.legend(fontsize=9)
                    
                lines = ax_comp.get_lines()

                # Initialize the maximum value to negative infinity
                max_value = float('-inf')

                # Iterate through each line
                for line in lines:
                    # Get the y-data of the line
                    y_data = line.get_ydata()
                    
                    # Find the maximum value in the y-data
                    line_max_value = max(y_data)
                    
                    # Update the overall maximum value if needed
                    max_value = max(max_value, line_max_value)
                
                ax_comp.set_xlim(last_xlims_chrom)
                ax_comp.set_ylim(0, max_value*1.1)
                ax_comp.set_xlabel('Retention/Migration Time (min)')
                ax_comp.set_ylabel('Intensity (AU)')
                
                canvas_comp.draw_idle()
        
        def check_all_checkbox_state_check():
            state = check_all_checkbox_state.get()
            
            if state:
                for i_i, i in enumerate(samples_dropdown_options):
                    chromatograms_checkboxes.change_state(i_i, "checked")
            else:
                for i_i, i in enumerate(samples_dropdown_options):
                    chromatograms_checkboxes.change_state(i_i, "unchecked")
                    
            draw_checked_samples()
        
        def disable_legend_checkbox_state_check():
            state = disable_legend_checkbox_state.get()
            
            if state:
                ax_comp.legend(fontsize=9)
            else:
                ax_comp.get_legend().remove()
                
            canvas_comp.draw_idle()
            
        def exit_compare_samples():
            global compare_samples_opened
            compare_samples_opened = False
            compare_samples.destroy()
        
        if compare_samples_opened:
            compare_samples.deiconify()
            compare_samples.lift()
            compare_samples.focus_set()
            return
        
        compare_samples_opened = True
        
        global item_text_comp, parent_text_comp, grand_parent_item_comp, grand_parent_text_comp
        if level == 1:
            item_text_comp = chromatograms_list.item(chromatograms_list.selection(), "text")
        if level == 2:
            parent_text_comp = chromatograms_list.item(chromatograms_list.selection(), "text")
            grand_parent_item_comp = chromatograms_list.parent(chromatograms_list.selection())
            grand_parent_text_comp = chromatograms_list.item(grand_parent_item_comp, "text")
    
        compare_samples = tk.Toplevel()
        icon = ImageTk.PhotoImage(ico_image)
        compare_samples.iconphoto(False, icon)
        compare_samples.withdraw()
        compare_samples.minsize(720, 480)
        compare_samples.bind("<Configure>", on_resize)
        compare_samples.bind("<KeyPress>", on_key_press)
        compare_samples.bind("<KeyRelease>", on_key_release)
        if level == 1:
            compare_samples.title(f"Samples Comparison - {item_text_comp}")
        else:
            compare_samples.title(f"Samples Comparison - {grand_parent_text_comp} - {parent_text_comp}")
        compare_samples.resizable(True, True)
        compare_samples.grid_rowconfigure(0, weight=0)
        compare_samples.grid_rowconfigure(1, weight=1)
        compare_samples.grid_columnconfigure(0, weight=0)
        compare_samples.grid_columnconfigure(1, weight=1)
        compare_samples.grab_set()
        compare_samples.protocol("WM_DELETE_WINDOW", exit_compare_samples)
        
        # Check/Uncheck all checkboxes checkbox
        check_all_checkbox_state = tk.BooleanVar(value=False)
        check_all_checkbox = ttk.Checkbutton(compare_samples, text="Check/Uncheck all samples", variable=check_all_checkbox_state, command=check_all_checkbox_state_check)
        check_all_checkbox.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="nw")
        ToolTip(check_all_checkbox, "Checks or unchecks all the samples checkboxes.")
        
        # Align chromatograms checkbox, deprecated function
        align_chromatograms_checkbox_state = tk.BooleanVar(value=False)
        align_chromatograms_checkbox = ttk.Checkbutton(compare_samples, text="Align Chromatograms/Electropherograms", variable=align_chromatograms_checkbox_state, command=align_chromatograms_checkbox_state_check)
        align_chromatograms_checkbox.grid(row=0, column=1, padx=10, pady=(10, 20), sticky="ne")
        ToolTip(align_chromatograms_checkbox, "Aligns the chromatograms/electropherograms. It's very dependent on the quality of the peaks, so adjusting quality thresholds may affect the alignment quality. The first time you check this box for a given set of QC thresholds it will take a while to calculate the alignments. DEPRECATED")
        
        # Enable legend
        disable_legend_checkbox_state = tk.BooleanVar(value=True)
        disable_legend_checkbox = ttk.Checkbutton(compare_samples, text="Enable/Disable Legend", variable=disable_legend_checkbox_state, command=disable_legend_checkbox_state_check)
        disable_legend_checkbox.grid(row=0, column=1, padx=10, pady=(0,0), sticky="se")
        ToolTip(disable_legend_checkbox, "Enables or disables the legend in the graph.")
        
        # Create the panned windows
        paned_window_cs = tk.PanedWindow(compare_samples, sashwidth=5, orient=tk.HORIZONTAL)
        paned_window_cs.grid(row=1, column=0, columnspan=2, padx = 10, pady = 10, sticky='nwse')
        
        # The treeview for the chromatograms
        chromatograms_checkboxes = CheckboxTreeview(paned_window_cs)
        chromatograms_checkboxes["show"] = "tree" #removes the header
        chromatograms_checkboxes.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")
        
        for i_i, i in enumerate(samples_dropdown_options):
            chromatograms_checkboxes.insert("", "end", i_i, text=i)
            if i == samples_dropdown.get():
                chromatograms_checkboxes.change_state(i_i, "checked")
            else:
                chromatograms_checkboxes.change_state(i_i, "unchecked")
                
        TreeviewTooltip(chromatograms_checkboxes)
        
        # The frame for the plotting area
        chromatogram_plots_compare_frame = ttk.Labelframe(paned_window_cs, text="Chromatogram/Electropherogram Viewer", style="chromatogram.TLabelframe")
        chromatogram_plots_compare_frame.grid(row=1, column=1, padx=10, pady=10, sticky="nsew")

        global canvas_comp, ax_comp, coordinate_label_comp, type_coordinate_comp
        fig_comp = plt.figure(figsize=(0, 0))
        ax_comp = fig_comp.add_subplot(111)
        canvas_comp = FigureCanvasTkAgg(fig_comp, master=chromatogram_plots_compare_frame)
        canvas_comp.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        ax_comp.set_xlabel('Retention/Migration Time (min)')
        ax_comp.set_ylabel('Intensity (AU)')
    
        # Add the left side widgets frame to the paned window
        paned_window_cs.add(chromatograms_checkboxes)
        paned_window_cs.paneconfigure(chromatograms_checkboxes, minsize=140)
        
        # Add the right side widgets frame to the paned window
        paned_window_cs.add(chromatogram_plots_compare_frame)
        paned_window_cs.paneconfigure(chromatogram_plots_compare_frame, minsize=400)
            
        # The vertical dashed line that serves as the indicator of the RT/MT hovered on
        vertical_line = ax_comp.axvline(x=-10000, color='black', linestyle='--', linewidth=1)
        
        # This collects the marked checkboxes to plot the samples
        loops = 0
        if level == 1:
            for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
        
                # Load retention times from file
                x_values_comp = gg_file.get_rt_array(i)
                    
                y_values_comp = []
                counter = 0
                for k in gg_file.list_chromatograms(i):
                    if "+".join(k.split("+")[:-1]) == item_text_comp:
                    
                        # Load intensities from file
                        y_values_temp = gg_file.get_chromatogram(i, k, 'smoothed')
                            
                        if counter == 0:
                            y_values_comp = y_values_temp
                            counter+= 1
                        else:
                            y_values_comp = [x + y for x, y in zip(y_values_comp, y_values_temp)]
                if i_i-(len(colors)*loops) == len(colors):
                    loops+= 1
                color = colors[i_i-(len(colors)*loops)]
                ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))
        if level == 2:
            for i_i, i in enumerate(chromatograms_checkboxes.get_checked()):
        
                # Load retention times from file
                x_values_comp = gg_file.get_rt_array(i)
                    
                # Load intensities from file
                y_values_comp = gg_file.get_chromatogram(i, f"{grand_parent_text_comp}+{parent_text_comp}", 'smoothed')
                    
                if i_i-(len(colors)*loops) == len(colors):
                    loops+= 1
                color = colors[i_i-(len(colors)*loops)]
                ax_comp.plot(x_values_comp, y_values_comp, linewidth=1, color=color, label = chromatograms_checkboxes.item(i, "text"))

        ax_comp.legend(fontsize=9)
                
        try:
            lines = ax_comp.get_lines()

            # Initialize the maximum value to negative infinity
            max_value = float('-inf')

            # Iterate through each line
            for line in lines:
                # Get the y-data of the line
                y_data = line.get_ydata()
                
                # Find the maximum value in the y-data
                line_max_value = max(y_data)
                
                # Update the overall maximum value if needed
                max_value = max(max_value, line_max_value)
                
            x_min = min(x_values_comp)
            x_max = max(x_values_comp)
            og_x_range_comp = (x_min-0.5, x_max+0.5)
            og_y_range_comp = (0, max_value*1.1)
        except:
            og_x_range_comp = (0, 1000)
            og_y_range_comp = (0, 1000)
            
        ax_comp.set_xlim(og_x_range_comp)
        ax_comp.set_ylim(og_y_range_comp)
        
        last_xlims_chrom = og_x_range_comp
        last_ylims_chrom = og_y_range_comp
        
        canvas_comp.draw()
        
        canvas_comp.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_comp, canvas_comp, True) if event.button == 3 else None) 
        zoom_selection_key_press_comp = canvas_comp.mpl_connect('key_press_event', lambda event: zoom_selection_compare(event, ax_comp, canvas_comp, type_coordinate_comp))
        zoom_selection_key_release_comp = canvas_comp.mpl_connect('key_release_event', lambda event: zoom_selection_compare(event, ax_comp, canvas_comp, type_coordinate_comp))
        zoom_selection_motion_notify_comp = canvas_comp.mpl_connect('motion_notify_event', lambda event: zoom_selection_compare(event, ax_comp, canvas_comp, type_coordinate_comp))
        zoom_selection_button_press_comp = canvas_comp.mpl_connect('button_press_event', lambda event: zoom_selection_compare(event, ax_comp, canvas_comp, type_coordinate_comp) if event.button == 1 else None)
        zoom_selection_button_release_comp = canvas_comp.mpl_connect('button_release_event', lambda event: zoom_selection_compare(event, ax_comp, canvas_comp, type_coordinate_comp) if event.button == 1 else None)
        on_scroll_event_comp = canvas_comp.mpl_connect('scroll_event', lambda event: on_scroll(event, ax_comp, canvas_comp, type_coordinate_comp))
        on_double_click_event_comp = canvas_comp.mpl_connect('button_press_event', lambda event: on_double_click(event, ax_comp, canvas_comp, og_x_range_comp, og_y_range_comp, type_coordinate_comp) if event.button == 1 else None)
        on_plot_hover_motion_comp = canvas_comp.mpl_connect('motion_notify_event', lambda event: on_plot_hover(event, ax_comp, canvas_comp, coordinate_label_comp, type_coordinate_comp, vertical_line, True))
        on_pan_press_comp = canvas_comp.mpl_connect('button_press_event', lambda event: on_pan(event, ax_comp, canvas_comp, type_coordinate_comp) if event.button == 1 else None)
        on_pan_release_comp = canvas_comp.mpl_connect('button_release_event', lambda event: on_pan(event, ax_comp, canvas_comp, type_coordinate_comp) if event.button == 1 else None)
        on_pan_motion_comp = canvas_comp.mpl_connect('motion_notify_event', lambda event: on_pan(event, ax_comp, canvas_comp, type_coordinate_comp))
        on_pan_right_click_press_comp = canvas_comp.mpl_connect('button_press_event', lambda event: on_pan_right_click(event, ax_comp, canvas_comp, type_coordinate) if event.button == 3 else None)
        on_pan_right_click_release_comp = canvas_comp.mpl_connect('button_release_event', lambda event: on_pan_right_click(event, ax_comp, canvas_comp, type_coordinate) if event.button == 3 else None)
        on_pan_right_click_motion_comp = canvas_comp.mpl_connect('motion_notify_event', lambda event: on_pan_right_click(event, ax_comp, canvas_comp, type_coordinate))
        
        coordinate_label_comp = tk.Label(chromatogram_plots_compare_frame, text="", anchor="e", font=("Segoe UI", 8), bg="white")
        coordinate_label_comp.place(relx=1.0, rely=0, anchor='ne')
        type_coordinate_comp = "chromatogram"
        coordinate_label_comp.lift()
        
        frame_width = chromatogram_plot_frame.winfo_width()
        frame_height = chromatogram_plot_frame.winfo_height()

        left_margin = 65
        right_margin = 10
        top_margin = 20
        bottom_margin = 45

        subplot_width = (frame_width - left_margin - right_margin) / frame_width
        subplot_height = (frame_height - top_margin - bottom_margin) / frame_height
        subplot_left = left_margin / frame_width
        subplot_bottom = bottom_margin / frame_height

        ax_comp.set_position([subplot_left, subplot_bottom, subplot_width, subplot_height])
        
        chromatograms_checkboxes.bind("<Button-1>", on_checkbox_click)
        chromatogram_plots_compare_frame.bind("<Configure>", lambda event, ax=ax_comp: adjust_subplot_size(event, ax_comp, canvas_comp))
        
        compare_samples.update_idletasks()
        window_width = compare_samples.winfo_width()
        window_height = compare_samples.winfo_height()
        screen_width = compare_samples.winfo_screenwidth()
        screen_height = compare_samples.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        compare_samples.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        compare_samples.deiconify()
        
    def check_qc_dist():
        global sample_for_qc_dist, qc_dist_opened, qc_dist
            
        def adjust_subplot_size_qc(event, ax_plot, ax_kde):
            # Get the current size of the graph frame
            frame_width = event.width
            frame_height = event.height
            
            # Default top and bottom margins
            top_margin = 10
            bottom_margin = 50
            
            # Define margin sizes (in pixels)
            left_margin = 40
            right_margin = frame_width - (frame_width*0.85)
            
            # Calculate the position of the subplot relative to the frame size
            subplot_width = (frame_width - left_margin - right_margin) / frame_width
            subplot_height = (frame_height - top_margin - bottom_margin) / frame_height
            subplot_left = left_margin / frame_width
            subplot_bottom = bottom_margin / frame_height
            
            # Define margin sizes (in pixels)
            left_margin_kde = (frame_width - right_margin)
            right_margin_kde = 10
            
            # Calculate the position of the subplot relative to the frame size
            subplot_width_kde = (frame_width - left_margin_kde - right_margin_kde) / frame_width
            subplot_height_kde = (frame_height - top_margin - bottom_margin) / frame_height
            subplot_left_kde = left_margin_kde / frame_width
            subplot_bottom_kde = bottom_margin / frame_height
            
            # Set the position of the subplot
            ax_plot.set_position([subplot_left, subplot_bottom, subplot_width, subplot_height])
            
            # Set the position of the subplot
            ax_kde.set_position([subplot_left_kde, subplot_bottom_kde, subplot_width_kde, subplot_height_kde])
        
        def exit_check_qc_dist():
            global qc_dist_opened
            qc_dist_opened = False
            qc_dist.destroy()
            
        def on_hover_qc_plot(event, canvas, ax, scatter, x_values, y_values, tooltip):
            show_mz = sort_by_mz_checkbox_state.get()
            
            if event.inaxes == ax:
                cont, ind = scatter.contains(event)
                if cont:
                    ind = ind['ind'][0]
                    name = names[ind]
                    x, y = x_values[ind], y_values[ind]
                    if show_mz:
                        tooltip.set_text(f'{name}\nm/z: {x}\nY: {y}')
                    else:
                        tooltip.set_text(f'{name}\nm/z: {mzs[x]}\nY: {y}')
                    tooltip.xy = (x, y)
                    
                    # Check if the tooltip will go out of the x-axis bounds
                    xlim = ax.get_xlim()
                    difference = (((xlim[0]+((xlim[1]-xlim[0])*0.5)) - x)/(xlim[1]-xlim[0]))*250
                    if (x) > xlim[0]+((xlim[1]-xlim[0])*0.5):  # If too close to the right edge
                        tooltip.set_x(difference)  # Move tooltip to the left
                    else:
                        tooltip.set_x(10)  # Default to the right
                        
                    tooltip.set_visible(True)
                    canvas.draw_idle()
                else:
                    tooltip.set_visible(False)
                    canvas.draw_idle()
                    
        def model_electro_migrations():
            def exit_check_electro_mig():
                electro_migrations.destroy()
                
            def on_hover_electro_mig(event, canvas, ax, scatter, x_values, y_values, tooltip):
                if event.inaxes == ax:
                    cont, ind = scatter.contains(event)
                    if cont:
                        ind = ind['ind'][0]
                        name = names[ind]
                        x, y = x_values[ind], y_values[ind]
                        tooltip.set_text(f'{name}\nX: {x:.1e}\nY: {y}')
                        tooltip.xy = (x, y)
                        tooltip.set_visible(True)
                        canvas.draw_idle()
                    else:
                        tooltip.set_visible(False)
                        canvas.draw_idle()
                        
            def ph_enter(event, hover_connect, tooltip, eq_r2, ax, canvas, scatter, trendline_base, trendline_plus20, trendline_minus20):
                try:
                    float(ph_entry.get())
                except:
                    error_window("Input a float in the pH entry field.")
                    return
                
                glycans_list, quality_colors, q_list, mass_list, me_list, rt_list = calculate_electrophoretic_migrations()
                
                # Calculate the trendline (linear regression)
                coefficients = np.polyfit(me_list, rt_list, 1)
                polynomial = np.poly1d(coefficients)
                trendline = polynomial(me_list)
                    
                # Calculate the trendline equation and R squared
                trendline_equation = f"y = {coefficients[0]:2f}x + {coefficients[1]:.2f}"
                y_mean = np.mean(rt_list)
                ss_tot = np.sum((rt_list - y_mean) ** 2)
                ss_res = np.sum((rt_list - trendline) ** 2)
                r_squared = f"R² = {1 - (ss_res / ss_tot):.2f}"
                
                fig_elec_mig = plt.figure(figsize=(0, 0))
                ax_elec_mig = fig_elec_mig.add_subplot(111)
                ax_elec_mig.set_xlabel(r'$q / M^{\alpha}$')
                ax_elec_mig.set_ylabel("t'")
                
                eq_r2.set_text(f"{trendline_equation}\n{r_squared}")
                                     
                scatter.set_offsets(np.c_[me_list, rt_list])
                trendline_base[0].set_data(me_list, trendline)
                trendline_plus20[0].set_data(me_list, [i*1.2 for i in trendline])
                trendline_minus20[0].set_data(me_list, [i*0.8 for i in trendline])
            
                ax.set_ylim(min(rt_list)-0.1, max(rt_list)+0.1)
                ax.set_xlim(min(me_list)-0.0001, max(me_list)+0.0001)
                
                canvas.draw()
                
                canvas.mpl_disconnect(hover_connect)
                
                hover_tooltip = canvas.mpl_connect('motion_notify_event', lambda event: on_hover_electro_mig(event, canvas, ax, scatter, me_list, rt_list, tooltip))
            
            def adjust_subplot_size_electro_mig(event, ax1):
                # Get the current size of the graph frame
                frame_width = event.width
                frame_height = event.height
                
                # Define margin sizes (in pixels)
                left_margin = 65
                right_margin = 10
                top_margin = 20
                bottom_margin = 45
                
                # Calculate the position of the subplot relative to the frame size
                subplot_width = (frame_width - left_margin - right_margin) / frame_width
                subplot_height = (frame_height - top_margin - bottom_margin) / frame_height
                subplot_left = left_margin / frame_width
                subplot_bottom = bottom_margin / frame_height
                
                # Set the position of the subplot
                ax1.set_position([subplot_left, subplot_bottom, subplot_width, subplot_height])
            
            def calculate_electrophoretic_migrations():
                glycans_list = []
                quality_colors = []
                q_list = []
                mass_list = []
                me_list = []
                rt_list = []
                for i in glycans_per_sample[sample_for_qc_dist]:
                    for j in glycans_per_sample[sample_for_qc_dist][i]:
                        
                        adduct_charge = adduct_combos_dict[j]['charges']
                        
                        for k_k, k in enumerate(glycans_per_sample[sample_for_qc_dist][i][j]['peaks']):
                            glycans_list.append(f"{i}_{j}_{k}")
                            
                            fails = 0
                            if glycans_per_sample[sample_for_qc_dist][i][j]['ppm'][k_k] < max_ppm[0] or glycans_per_sample[sample_for_qc_dist][i][j]['ppm'][k_k] > max_ppm[1]:
                                fails += 1
                            if glycans_per_sample[sample_for_qc_dist][i][j]['iso'][k_k] < iso_fit_score:
                                fails += 1
                            if glycans_per_sample[sample_for_qc_dist][i][j]['curve'][k_k] < curve_fit_score:
                                fails += 1
                            if glycans_per_sample[sample_for_qc_dist][i][j]['sn'][k_k] < s_to_n:
                                fails += 1
                            if fails == 0:
                                quality_colors.append('green')
                            elif fails == 1:
                                quality_colors.append('#c7af12')
                            else:
                                quality_colors.append('red')
                                
                            q = 1
                            composition = General_Functions.form_to_comp_glycans(i)
                            if 'S' in composition.keys():
                                q -= (composition['S']*(1/(1+10**(2.6-float(ph_entry.get())))))
                            if 'G' in composition.keys():
                                q -= (composition['G']*(1/(1+10**(2.92-float(ph_entry.get())))))
                            q_list.append(q)
                            mass_list.append(float(glycans_per_sample[sample_for_qc_dist][i][j]['mz'])*adduct_charge)
                            me_list.append(q_list[-1]/(mass_list[-1]**1/2))
                            rt_list.append(k)
                reference_rt_id = mass_list.index(max(mass_list))
                reference_rt = rt_list[reference_rt_id]
                for i_i, i in enumerate(rt_list):
                    rt_list[i_i] = float("%.2f" % round(i/reference_rt, 2))
                    
                return glycans_list, quality_colors, q_list, mass_list, me_list, rt_list
                    
            electro_migrations = tk.Toplevel()
            icon = ImageTk.PhotoImage(ico_image)
            electro_migrations.iconphoto(False, icon)
            electro_migrations.withdraw()
            electro_migrations.minsize(500, 500)
            electro_migrations.bind("<Configure>", on_resize)
            electro_migrations.title(f"Electrophoretic Migration Modelling - {sample_for_qc_dist}")
            electro_migrations.resizable(True, True)
            electro_migrations.protocol("WM_DELETE_WINDOW", exit_check_electro_mig)
            
            electro_migrations.grid_rowconfigure(0, weight=0)
            electro_migrations.grid_rowconfigure(1, weight=1)
            electro_migrations.grid_columnconfigure(0, weight=1)
            
            ph_label = ttk.Label(electro_migrations, text='Background Electrolyte pH: ', font=("Segoe UI", list_font_size))
            ph_label.grid(row=0, column=0, padx=10, pady=(10, 10), sticky="nsw")
            
            ph_entry = ttk.Entry(electro_migrations, width=10)
            ph_entry.insert(0, "2.3")
            ph_entry.grid(row=0, column=0, padx=(180, 10), pady=(10, 10), sticky="nsw")
    
            glycans_list, quality_colors, q_list, mass_list, me_list, rt_list = calculate_electrophoretic_migrations()
            
            # Calculate the trendline (linear regression)
            coefficients = np.polyfit(me_list, rt_list, 1)
            polynomial = np.poly1d(coefficients)
            trendline = polynomial(me_list)
            
            # Calculate the trendline equation and R squared
            trendline_equation = f"y = {coefficients[0]:2f}x + {coefficients[1]:.2f}"
            y_mean = np.mean(rt_list)
            ss_tot = np.sum((rt_list - y_mean) ** 2)
            ss_res = np.sum((rt_list - trendline) ** 2)
            r_squared = f"R² = {1 - (ss_res / ss_tot):.2f}"
            
            fig_elec_mig = plt.figure(figsize=(0, 0))
            ax_elec_mig = fig_elec_mig.add_subplot(111)
            ax_elec_mig.set_xlabel(r'$q / M^{\alpha}$')
            ax_elec_mig.set_ylabel("t'")
            
            eq_r2 = ax_elec_mig.text(0.95, 0.95, f"{trendline_equation}\n{r_squared}",
                                     horizontalalignment='right',
                                     verticalalignment='top',
                                     transform=plt.gca().transAxes)
            
            global elec_mig_scatter, canvas_elec_mig
            elec_mig_scatter = ax_elec_mig.scatter(me_list, rt_list, s=1, c=quality_colors)
            trendline_base = ax_elec_mig.plot(me_list, trendline, color='red', label='Trendline')
            trendline_plus20 = ax_elec_mig.plot(me_list, [i*1.2 for i in trendline], color='red', linestyle = ":", label='Trendline')
            trendline_minus20 = ax_elec_mig.plot(me_list, [i*0.8 for i in trendline], color='red', linestyle = ":", label='Trendline')
            
            ax_elec_mig.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax_elec_mig.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            
            ax_elec_mig.set_ylim(min(rt_list)-0.1, max(rt_list)+0.1)
            ax_elec_mig.set_xlim(min(me_list)-0.0001, max(me_list)+0.0001)
            
            canvas_elec_mig = FigureCanvasTkAgg(fig_elec_mig, master=electro_migrations)
            canvas_elec_mig.draw()
            canvas_elec_mig.get_tk_widget().grid(row=1, column=0, padx=(10, 10), pady=(10, 10), sticky="nswe")
        
            tooltip_elec_mig = ax_elec_mig.annotate('', xy=(0, 0), xytext=(10, -20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.9), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'), clip_on=True)
        
            hover_tooltip = canvas_elec_mig.mpl_connect('motion_notify_event', lambda event: on_hover_electro_mig(event, canvas_elec_mig, ax_elec_mig, elec_mig_scatter, me_list, rt_list, tooltip_elec_mig))
            canvas_elec_mig.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_elec_mig, canvas_elec_mig, False) if event.button == 3 else None)
            
            ph_entry.bind("<Return>", lambda event: ph_enter(event, hover_tooltip, tooltip_elec_mig, eq_r2, ax_elec_mig, canvas_elec_mig, elec_mig_scatter, trendline_base, trendline_plus20, trendline_minus20))
        
            electro_migrations.bind("<Configure>", lambda event: adjust_subplot_size_electro_mig(event, ax_elec_mig))
        
            electro_migrations.update_idletasks()
            electro_migrations.deiconify()
            window_width = electro_migrations.winfo_width()
            window_height = electro_migrations.winfo_height()
            screen_width = electro_migrations.winfo_screenwidth()
            screen_height = electro_migrations.winfo_screenheight()
            x_position = (screen_width - window_width) // 2
            y_position = (screen_height - window_height) // 2
            electro_migrations.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
            
        def sort_by_mz_checkbox_command():
            global hover_tooltip_ppm, hover_tooltip_isofit, hover_tooltip_snplot, hover_tooltip_curvefitplot, ax_ppmplot, ax_snplot, ax_isofitplot, ax_curvefitplot
            state = sort_by_mz_checkbox_state.get()
            
            if state:
                #adjust ppm_plot
                canvas_ppmplot.mpl_disconnect(hover_tooltip_isofit)
                
                ppm_scatter.set_offsets(np.c_[mzs, ppm_list])
                ax_ppmplot.set_xlim([min(mzs)-((max(mzs)-min(mzs))*0.05), max(mzs)*1.05])
                ax_ppmplot.set_xlabel('m/z')
                
                canvas_ppmplot.draw_idle()
                
                hover_tooltip_ppm = canvas_ppmplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_ppmplot, ax_ppmplot, ppm_scatter, mzs, ppm_list, tooltip_ppmplot))
                
                #adjust isofit_plot
                canvas_isofitplot.mpl_disconnect(hover_tooltip_ppm)
                
                isofitplot_scatter.set_offsets(np.c_[mzs, iso_fit_list])
                ax_isofitplot.set_xlim([min(mzs)-((max(mzs)-min(mzs))*0.05), max(mzs)*1.05])
                ax_isofitplot.set_xlabel('m/z')
                
                canvas_isofitplot.draw_idle()
                
                hover_tooltip_isofit = canvas_isofitplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_isofitplot, ax_isofitplot, isofitplot_scatter, mzs, iso_fit_list, tooltip_isofitplot))
                
                #adjust sn_plot
                canvas_snplot.mpl_disconnect(hover_tooltip_snplot)
                
                snplot_scatter.set_offsets(np.c_[mzs, sn_list])
                ax_snplot.set_xlim([min(mzs)-((max(mzs)-min(mzs))*0.05), max(mzs)*1.05])
                ax_snplot.set_xlabel('m/z')
                
                canvas_snplot.draw_idle()
                
                hover_tooltip_snplot = canvas_snplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_snplot, ax_snplot, snplot_scatter, mzs, sn_list, tooltip_snplot))
                
                #adjust curvefit_plot
                canvas_curvefitplot.mpl_disconnect(hover_tooltip_curvefitplot)
                
                curvefitplot_scatter.set_offsets(np.c_[mzs, curve_fit_list])
                ax_curvefitplot.set_xlim([min(mzs)-((max(mzs)-min(mzs))*0.05), max(mzs)*1.05])
                ax_curvefitplot.set_xlabel('m/z')
                
                canvas_curvefitplot.draw_idle()
                
                hover_tooltip_curvefitplot = canvas_curvefitplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_curvefitplot, ax_curvefitplot, curvefitplot_scatter, mzs, curve_fit_list, tooltip_curvefitplot))
            else:
                #adjust ppm_plot
                canvas_ppmplot.mpl_disconnect(hover_tooltip_ppm)
                
                ppm_scatter.set_offsets(np.c_[range(len(ppm_list)), ppm_list])
                ax_ppmplot.set_xlim([-len(ppm_list)*0.05, len(ppm_list)*1.05])
                ax_ppmplot.set_xlabel('Feature')
                
                canvas_ppmplot.draw_idle()
                
                hover_tooltip_ppm = canvas_ppmplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_ppmplot, ax_ppmplot, ppm_scatter, range(len(ppm_list)), ppm_list, tooltip_ppmplot))
                
                #adjust isofit_plot
                canvas_isofitplot.mpl_disconnect(hover_tooltip_isofit)
                
                isofitplot_scatter.set_offsets(np.c_[range(len(iso_fit_list)), iso_fit_list])
                ax_isofitplot.set_xlim([-len(iso_fit_list)*0.05, len(iso_fit_list)*1.05])
                ax_isofitplot.set_xlabel('Feature')
                
                canvas_isofitplot.draw_idle()
                
                hover_tooltip_isofit = canvas_isofitplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_isofitplot, ax_isofitplot, isofitplot_scatter, range(len(iso_fit_list)), iso_fit_list, tooltip_isofitplot))
                
                #adjust sn_plot
                canvas_snplot.mpl_disconnect(hover_tooltip_snplot)
                
                snplot_scatter.set_offsets(np.c_[range(len(sn_list)), sn_list])
                ax_snplot.set_xlim([-len(sn_list)*0.05, len(sn_list)*1.05])
                ax_snplot.set_xlabel('Feature')
                
                canvas_snplot.draw_idle()
                
                hover_tooltip_snplot = canvas_snplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_snplot, ax_snplot, snplot_scatter, range(len(sn_list)), sn_list, tooltip_snplot))
                
                #adjust sn_plot
                canvas_curvefitplot.mpl_disconnect(hover_tooltip_curvefitplot)
                
                curvefitplot_scatter.set_offsets(np.c_[range(len(curve_fit_list)), curve_fit_list])
                ax_curvefitplot.set_xlim([-len(curve_fit_list)*0.05, len(curve_fit_list)*1.05])
                ax_curvefitplot.set_xlabel('Feature')
                
                canvas_curvefitplot.draw_idle()
                
                hover_tooltip_curvefitplot = canvas_curvefitplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_curvefitplot, ax_curvefitplot, curvefitplot_scatter, range(len(curve_fit_list)), curve_fit_list, tooltip_curvefitplot))
                
        
        if qc_dist_opened:
            qc_dist.deiconify()
            qc_dist.lift()
            qc_dist.focus_set()
            return
        
        qc_dist_opened = True
        
        sample_for_qc_dist = selected_item
            
        qc_dist = tk.Toplevel()
        icon = ImageTk.PhotoImage(ico_image)
        qc_dist.iconphoto(False, icon)
        qc_dist.withdraw()
        qc_dist.minsize(400, 400)
        qc_dist.bind("<Configure>", on_resize)
        qc_dist.title(f"QC Scores Distribution - {sample_for_qc_dist}")
        qc_dist.resizable(True, True)
        qc_dist.protocol("WM_DELETE_WINDOW", exit_check_qc_dist)
        
        qc_dist.grid_rowconfigure(0, weight=0)
        qc_dist.grid_rowconfigure(1, weight=1)
        qc_dist.grid_rowconfigure(2, weight=1)
        qc_dist.grid_rowconfigure(3, weight=0)
        qc_dist.grid_columnconfigure(0, weight=1)
        qc_dist.grid_columnconfigure(1, weight=1)
            
        names = []
        mzs = []
        quality_colors = []
        
        ppm_list = []
        iso_fit_list = []
        curve_fit_list = []
        sn_list = []
        
        for i in glycans_per_sample[sample_for_qc_dist]: #going through glycans
            for k in glycans_per_sample[sample_for_qc_dist][i]: #going through adducts
                for j_j, j in enumerate(glycans_per_sample[sample_for_qc_dist][i][k]['peaks']):
                    mzs.append(glycans_per_sample[sample_for_qc_dist][i][k]['mz'])
                    names.append(f"{i}_{k}_{j}")
                    
                    fails = check_peak_quality(ppm = glycans_per_sample[sample_for_qc_dist][i][k]['ppm'][j_j], iso = glycans_per_sample[sample_for_qc_dist][i][k]['iso'][j_j], curve = glycans_per_sample[sample_for_qc_dist][i][k]['curve'][j_j], sn = glycans_per_sample[sample_for_qc_dist][i][k]['sn'][j_j])
                    
                    if fails == 0:
                        quality_colors.append('green')
                    elif fails == 1:
                        quality_colors.append('#c7af12')
                    else:
                        quality_colors.append('red')
                        
                ppm_list+=glycans_per_sample[sample_for_qc_dist][i][k]['ppm']
                iso_fit_list+=glycans_per_sample[sample_for_qc_dist][i][k]['iso']
                curve_fit_list+=glycans_per_sample[sample_for_qc_dist][i][k]['curve']
                sn_list+=glycans_per_sample[sample_for_qc_dist][i][k]['sn']
        
        sort_by_mz_checkbox_state = tk.BooleanVar(value=False)
        sort_by_mz_checkbox = ttk.Checkbutton(qc_dist, text="Sort peaks by m/z on the x-axis", variable=sort_by_mz_checkbox_state, command=sort_by_mz_checkbox_command)
        sort_by_mz_checkbox.grid(row=0, column=1, padx=10, pady=10, sticky="ne")
        ToolTip(sort_by_mz_checkbox, "The peaks are arranged by default by feature number. This provides an even distribution of the peaks along the x-axis. If you enable this option, they'll be distributed by m/z instead.")
        sort_by_mz_checkbox_state.set(value=True)
        
        global qc_ppm_line1, qc_ppm_line2, qc_curvefit_line, qc_sn_line, qc_isofit_line, canvas_ppmplot, canvas_curvefitplot, canvas_isofitplot, canvas_snplot, ppm_scatter, isofitplot_scatter, curvefitplot_scatter, snplot_scatter, ax_ppmplot, ax_snplot, ax_curvefitplot, ax_isofitplot
        
        #PPM plot
        ppm_plot_frame = ttk.Labelframe(qc_dist, text="PPM Error:", style="qcp_frame.TLabelframe")
        ppm_plot_frame.grid(row=1, column=0, padx=10, pady=(10, 10), sticky="nsew")
        
        fig_ppmplot = plt.figure(figsize=(4.5, 3.5))
        gs_ppmplot = fig.add_gridspec(1, 2, width_ratios=[5, 1])
        
        ax_ppmplot = fig_ppmplot.add_subplot(gs_ppmplot[0])
        ax_ppmplot_kde = fig_ppmplot.add_subplot(gs_ppmplot[1], sharey=ax_ppmplot)
        ax_ppmplot.set_xlabel('Feature')
        
        ppm_scatter = ax_ppmplot.scatter(range(len(ppm_list)), ppm_list, s=1, c=quality_colors)
        qc_ppm_line1 = ax_ppmplot.axhline(max_ppm[0], linestyle='--', linewidth=1, color='blue')
        qc_ppm_line2 = ax_ppmplot.axhline(max_ppm[1], linestyle='--', linewidth=1, color='blue')
        ax_ppmplot.axhline(0, linestyle='--', linewidth=1, color='black')
        
        values_ppm = np.array(ppm_list)
        kde_ppm = gaussian_kde(values_ppm)
        y_values_ppm = np.linspace(ax_ppmplot.get_ylim()[0], ax_ppmplot.get_ylim()[1], 500)
        kde_ppm_values = kde_ppm(y_values_ppm)
        
        ax_ppmplot_kde.plot(kde_ppm_values, y_values_ppm)
        ax_ppmplot_kde.axis('off')
        
        fig_ppmplot.subplots_adjust(wspace=0)
        
        canvas_ppmplot = FigureCanvasTkAgg(fig_ppmplot, master=ppm_plot_frame)
        canvas_ppmplot.draw()
        canvas_ppmplot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        tooltip_ppmplot = ax_ppmplot.annotate('', xy=(0, 0), xytext=(10, -20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.9), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'), clip_on=True)
        
        global hover_tooltip_ppm
        hover_tooltip_ppm = canvas_ppmplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_ppmplot, ax_ppmplot, ppm_scatter, range(len(ppm_list)), ppm_list, tooltip_ppmplot))
        canvas_ppmplot.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_ppmplot, canvas_ppmplot, False) if event.button == 3 else None)
        ppm_plot_resizing = ppm_plot_frame.bind("<Configure>", lambda event: adjust_subplot_size_qc(event, ax_ppmplot, ax_ppmplot_kde))
        
        
        # Isotopic Fittings plot
        isofit_plot_frame = ttk.Labelframe(qc_dist, text="Isotopic Fittings:", style="qcp_frame.TLabelframe")
        isofit_plot_frame.grid(row=2, column=0, padx=10, pady=(10, 10), sticky="nsew")
        
        fig_isofitplot = plt.figure(figsize=(4.5, 3.5))
        gs_isofitplot = fig.add_gridspec(1, 2, width_ratios=[5, 1])
        
        ax_isofitplot = fig_isofitplot.add_subplot(gs_isofitplot[0])
        ax_isofitplot_kde = fig_isofitplot.add_subplot(gs_isofitplot[1], sharey=ax_isofitplot)
        ax_isofitplot.set_xlabel('Feature')
        
        isofitplot_scatter = ax_isofitplot.scatter(range(len(iso_fit_list)), iso_fit_list, s=1, c=quality_colors)
        qc_isofit_line = ax_isofitplot.axhline(iso_fit_score, linestyle='--', linewidth=1, color='blue')
        
        values_isofit = np.array(iso_fit_list)
        kde_isofit = gaussian_kde(values_isofit)
        y_values_isofit = np.linspace(ax_isofitplot.get_ylim()[0], ax_isofitplot.get_ylim()[1], 500)
        kde_isofit_values = kde_isofit(y_values_isofit)
        
        ax_isofitplot_kde.plot(kde_isofit_values, y_values_isofit)
        ax_isofitplot_kde.axis('off')
        
        fig_isofitplot.subplots_adjust(wspace=0)
        
        canvas_isofitplot = FigureCanvasTkAgg(fig_isofitplot, master=isofit_plot_frame)
        canvas_isofitplot.draw()
        canvas_isofitplot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        tooltip_isofitplot = ax_isofitplot.annotate('', xy=(0, 0), xytext=(10, -20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.9), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'), clip_on=True)
        
        global hover_tooltip_isofit
        hover_tooltip_isofit = canvas_isofitplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_isofitplot, ax_isofitplot, isofitplot_scatter, range(len(iso_fit_list)), iso_fit_list, tooltip_isofitplot))
        canvas_isofitplot.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_isofitplot, canvas_isofitplot, False) if event.button == 3 else None)
        isofit_plot_resizing = isofit_plot_frame.bind("<Configure>", lambda event: adjust_subplot_size_qc(event, ax_isofitplot, ax_isofitplot_kde))
        
        
        #Signal-to-Noise ratio plot
        sn_plot_frame = ttk.Labelframe(qc_dist, text="Signal-to-Noise Ratio:", style="qcp_frame.TLabelframe")
        sn_plot_frame.grid(row=1, column=1, padx=10, pady=(10, 10), sticky="nsew")
        
        fig_snplot = plt.figure(figsize=(4.5, 3.5))
        gs_snplot = fig.add_gridspec(1, 2, width_ratios=[5, 1])
        
        ax_snplot = fig_snplot.add_subplot(gs_snplot[0])
        ax_snplot_kde = fig_snplot.add_subplot(gs_snplot[1], sharey=ax_snplot)
        ax_snplot.set_xlabel('Feature')
        
        ax_snplot.set_yscale('log')
        
        snplot_scatter = ax_snplot.scatter(range(len(sn_list)), sn_list, s=1, c=quality_colors)
        qc_sn_line = ax_snplot.axhline(s_to_n, linestyle='--', linewidth=1, color='blue')
        ax_snplot.axhline(3, linestyle='--', linewidth=1, color='red')
        ax_snplot.axhline(10, linestyle='--', linewidth=1, color='green')
        
        values_sn = np.array(sn_list)
        kde_sn = gaussian_kde(values_sn)
        y_values_sn = np.linspace(ax_snplot.get_ylim()[0], ax_snplot.get_ylim()[1], 500)
        kde_sn_values = kde_sn(y_values_sn)
        
        ax_snplot_kde.plot(kde_sn_values, y_values_sn)
        ax_snplot_kde.axis('off')
        
        fig_snplot.subplots_adjust(wspace=0)
        
        canvas_snplot = FigureCanvasTkAgg(fig_snplot, master=sn_plot_frame)
        canvas_snplot.draw()
        canvas_snplot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        tooltip_snplot = ax_snplot.annotate('', xy=(0, 0), xytext=(10, -20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.9), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'), clip_on=True)
        
        global hover_tooltip_snplot
        hover_tooltip_snplot = canvas_snplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_snplot, ax_snplot, snplot_scatter, range(len(sn_list)), sn_list, tooltip_snplot))
        canvas_snplot.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_snplot, canvas_snplot, False) if event.button == 3 else None)
        sn_plot_resizing = sn_plot_frame.bind("<Configure>", lambda event: adjust_subplot_size_qc(event, ax_snplot, ax_snplot_kde))
        
        #Curve Fittings plot
        curvefit_plot_frame = ttk.Labelframe(qc_dist, text="Curve Fittings:", style="qcp_frame.TLabelframe")
        curvefit_plot_frame.grid(row=2, column=1, padx=10, pady=(10, 10), sticky="nsew")
        
        fig_curvefitplot = plt.figure(figsize=(4.5, 3.5))
        gs_curvefitplot = fig.add_gridspec(1, 2, width_ratios=[5, 1])
        
        ax_curvefitplot = fig_curvefitplot.add_subplot(gs_curvefitplot[0])
        ax_curvefitplot_kde = fig_curvefitplot.add_subplot(gs_curvefitplot[1], sharey=ax_curvefitplot)
        ax_curvefitplot.set_xlabel('Feature')
        
        curvefitplot_scatter = ax_curvefitplot.scatter(range(len(curve_fit_list)), curve_fit_list, s=1, c=quality_colors)
        qc_curvefit_line = ax_curvefitplot.axhline(curve_fit_score, linestyle='--', linewidth=1, color='blue')
        
        values_curve = np.array(curve_fit_list)
        kde_curve = gaussian_kde(values_curve)
        y_values_curve = np.linspace(ax_curvefitplot.get_ylim()[0], ax_curvefitplot.get_ylim()[1], 500)
        kde_curve_values = kde_curve(y_values_curve)
        
        ax_curvefitplot_kde.plot(kde_curve_values, y_values_curve)
        ax_curvefitplot_kde.axis('off')
        
        fig_curvefitplot.subplots_adjust(wspace=0)
        
        canvas_curvefitplot = FigureCanvasTkAgg(fig_curvefitplot, master=curvefit_plot_frame)
        canvas_curvefitplot.draw()
        canvas_curvefitplot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        tooltip_curvefitplot = ax_curvefitplot.annotate('', xy=(0, 0), xytext=(10, -20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.9), arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'), clip_on=True)
        
        global hover_tooltip_curvefitplot
        hover_tooltip_curvefitplot = canvas_curvefitplot.mpl_connect('motion_notify_event', lambda event: on_hover_qc_plot(event, canvas_curvefitplot, ax_curvefitplot, curvefitplot_scatter, range(len(curve_fit_list)), curve_fit_list, tooltip_curvefitplot))
        canvas_curvefitplot.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_curvefitplot, canvas_curvefitplot, False) if event.button == 3 else None)
        curvefit_plot_resizing = curvefit_plot_frame.bind("<Configure>", lambda event: adjust_subplot_size_qc(event, ax_curvefitplot, ax_curvefitplot_kde))
        
        
        # Important discussion necessary for this button: The sialic acids are anionic, so it should be decreasing the charge value, instead of summing... but to do that reliably, you first need to also calculate the cationic charges, otherwise values are negative. How to deal with that? I managed to model PKa values of molecules using Python, but it needs the smiles structure... expect user to input that? Make a list of tags for that? What's the default charges for neutral glycans with reduced end? And without reducing end? Currently defaulting glycans without sialic acids to q=1, independent of tag charges and such...
        
        model_electro_button = ttk.Button(qc_dist, text="Model Electrophoretic Migrations", style="small_button_style1.TButton", command=model_electro_migrations, state=tk.NORMAL)
        model_electro_button.grid(row=3, column=0, columnspan=2, padx=10, pady=(10, 10), sticky="nsw")
        ToolTip(model_electro_button, "Allows you to plot modelled electrophoretic migrations based on the Classical Polymer Model, as described in Barroso A et al., 2015 (Analytica Chimica Acta, Elsevier). EXPERIMENTAL.")
        
        sort_by_mz_checkbox_command()
        
        qc_dist.update_idletasks()
        qc_dist.deiconify()
        window_width = qc_dist.winfo_width()
        window_height = qc_dist.winfo_height()
        screen_width = qc_dist.winfo_screenwidth()
        screen_height = qc_dist.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        qc_dist.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
    def qcp_enter(event):
        global qcp_enter_delay
        try:
            main_window.after_cancel(qcp_enter_delay)
        except:
            pass
        qcp_enter_delay = main_window.after(100, qcp_enter_auxiliary_function)
        
    def qcp_enter_auxiliary_function():
        global max_ppm, iso_fit_score, curve_fit_score, s_to_n, ppm_scatter, isofitplot_scatter, curvefitplot_scatter, snplot_scatter, elec_mig_scatter, canvas_elec_mig
        # Attempt to collect the new QC thresholds from the spinboxes
        try:
            min_ppm_range = float(ppm_error_min_entry.get())
            sn_temp = float(s_n_entry.get())
            max_ppm_range = float(ppm_error_max_entry.get())
            curve_fit_temp = float(curve_fit_entry.get())
            iso_fit_temp = float(iso_fit_entry.get())
        except:
            return
        
        # Check if the PPM range is correct
        if min_ppm_range > max_ppm_range:
            return
            
        max_ppm = (min_ppm_range, max_ppm_range)
        iso_fit_score = iso_fit_temp
        curve_fit_score = curve_fit_temp
        s_to_n = sn_temp
        
        try:
            qc_ppm_line1.set_ydata([max_ppm[0]])
            qc_ppm_line2.set_ydata([max_ppm[1]])
            qc_isofit_line.set_ydata([iso_fit_score])
            qc_curvefit_line.set_ydata([curve_fit_score])
            qc_sn_line.set_ydata([s_to_n])
            
            new_colors = []
            for i in glycans_per_sample[sample_for_qc_dist]: #going through glycans
                for k in glycans_per_sample[sample_for_qc_dist][i]: #going through adducts
                    for j_j, j in enumerate(glycans_per_sample[sample_for_qc_dist][i][k]['peaks']):
                        
                        fails = check_peak_quality(ppm = glycans_per_sample[sample_for_qc_dist][i][k]['ppm'][j_j], iso = glycans_per_sample[sample_for_qc_dist][i][k]['iso'][j_j], curve = glycans_per_sample[sample_for_qc_dist][i][k]['curve'][j_j], sn = glycans_per_sample[sample_for_qc_dist][i][k]['sn'][j_j])
                        
                        if fails == 0:
                            new_colors.append('green')
                        elif fails == 1:
                            new_colors.append('#c7af12')
                        else:
                            new_colors.append('red')
            
            ppm_scatter.set_color(new_colors)
            isofitplot_scatter.set_color(new_colors)
            curvefitplot_scatter.set_color(new_colors)
            snplot_scatter.set_color(new_colors)
            
            canvas_ppmplot.draw_idle()
            canvas_curvefitplot.draw_idle()
            canvas_isofitplot.draw_idle()
            canvas_snplot.draw_idle()
        except:
            pass
            
        try:
            elec_mig_scatter.set_color(new_colors)
            canvas_elec_mig.draw_idle()
        except:
            pass
            
        populate_treeview()
        toggle_hide_bad_peaks()
        
    def on_key_release_filter(event, entry):
        # Get the current content of the entry widget
        global filter_list_update_timer
        try:
            main_window.after_cancel(filter_list_update_timer)
        except:
            pass
        filter_list_update_timer = main_window.after(500, populate_treeview)
            
    def clear_treeview_selection(event):
        global plot_graph_button
        plot_graph_button.config(state=tk.NORMAL)
        selected_items = chromatograms_list.selection()
        chromatograms_list.selection_remove(selected_items)
        
    def plot_graph_window():
        
        def exit_plot_window():
            plot_window.destroy()
            
        def copy_graph_data_clipboard():
            glycans_list = []
            abundance_list = []
            sn_list = []
            for i in glycans:
                glycans_list.append(i)
                abundance_list.append(str(glycans[i]['abundance']))
                sn_list.append(str(glycans[i]['sn']))
            plot_window.clipboard_clear()
            if mode == 'compare_samples':
                plot_window.clipboard_append(f"Sample\tAbundance\tS/N\n")
            else:
                plot_window.clipboard_append(f"Glycan\tAbundance\tS/N\n")
            for i_i, i in enumerate(glycans_list):
                plot_window.clipboard_append(f"{i}\t{abundance_list[i_i]}\t{sn_list[i_i]}\n")
            
        def on_hover_column_graph(event, canvas, ax, bars, labels, tooltip):
            if event.inaxes == ax:
                for bar, label in zip(bars, labels):
                    if bar.contains(event)[0]:
                        y = bar.get_height()
                        tooltip.set_text(f'{label}\n{y:.1e}')
                        tooltip.xy = (event.xdata, event.ydata)
                        
                        # Check if the tooltip will go out of the x-axis bounds
                        xlim = ax.get_xlim()
                        difference = (((xlim[0]+((xlim[1]-xlim[0])*0.9)) - event.xdata)/(xlim[1]-xlim[0]))*250
                        if (event.xdata) > xlim[0]+((xlim[1]-xlim[0])*0.9):  # If too close to the right edge
                            tooltip.set_x(difference)  # Move tooltip to the left
                            tooltip.set_y(10)
                        else:
                            tooltip.set_x(10)  # Default to the right
                            tooltip.set_y(10)
                            
                        tooltip.set_visible(True)
                        canvas.draw_idle()
                        break
                else:
                    tooltip.set_visible(False)
                    canvas.draw_idle()
        
        def determine_treeview_level(treeview, selected_item):
            level = 0
            parent_item = selected_item
            while parent_item:
                level += 1
                parent_item = treeview.parent(parent_item)
            return level
            
        def adjust_subplot_size_plot_window(event, ax1, graph):
            # Get the current size of the graph frame
            frame_width = event.width
            frame_height = event.height
            
            # Define margin sizes (in pixels)
            left_margin_ax1 = 85
            right_margin_ax1 = 10
            
            if graph == 'abundance':
                top_margin_ax1 = 40
                bottom_margin_ax1 = 40
            else:
                top_margin_ax1 = 10
                bottom_margin_ax1 = 70
            
            # Calculate the position of the subplot relative to the frame size
            subplot_width_ax1 = (frame_width - left_margin_ax1 - right_margin_ax1) / frame_width
            subplot_height_ax1 = (frame_height - top_margin_ax1 - bottom_margin_ax1) / frame_height
            subplot_left_ax1 = left_margin_ax1 / frame_width
            subplot_bottom_ax1 = bottom_margin_ax1 / frame_height
            
            # Set the position of the subplot
            ax1.set_position([subplot_left_ax1, subplot_bottom_ax1, subplot_width_ax1, subplot_height_ax1])
               
        if len(selected_chromatograms) > 1:
            mode = 'compare_glycans'
        elif len(chromatograms_list.selection()) == 0:
            mode = 'good'
        elif len(chromatograms_list.selection()) == 1:
            mode = 'compare_samples'
        
        if mode == 'good' or mode == 'compare_glycans':
            glycans = {}
            if mode == 'good':
                chromatogram_children = chromatograms_list.get_children()
                for i in chromatogram_children:
                    if chromatograms_list.item(i, 'text') != "Base Peak Chromatogram/Electropherogram" and "good" in chromatograms_list.item(i, 'tags'):
                        glycans[chromatograms_list.item(i, 'text')] = {'abundance': 0, 'sn': 0}
            else:
                treeview_levels = [x[1] for x in selected_chromatograms]
                for i_i, i in enumerate(selected_chromatograms):
                    if i[1] == 2:
                        glycan_text = i[0].split("+")[0]
                        adduct_text = i[0].split("+")[-1]
                        if min(treeview_levels) == 1:
                            if glycan_text not in glycans.keys():
                                glycans[glycan_text] = {'abundance': 0, 'sn': 0}
                        else:
                            if glycan_text+'_'+adduct_text.split(' ')[0] not in glycans.keys():
                                glycans[glycan_text+'_'+adduct_text.split(' ')[0]] = {'abundance': 0, 'sn': 0}
                    elif i[1] == 1:
                        glycan_text = i[0]
                        if glycan_text not in glycans.keys():
                            glycans[glycan_text] = {'abundance': 0, 'sn': 0}
            
            for i in glycans:
                if mode == 'good' or (mode == 'compare_glycans' and min(treeview_levels) == 1):
                    temp_sn = []
                    for j in glycans_per_sample[selected_item][i]:
                        for k in glycans_per_sample[selected_item][i][j]:
                            for l_l, l in enumerate(glycans_per_sample[selected_item][i][j]['auc']):
                                glycans[i]['abundance'] += l
                                temp_sn.append(glycans_per_sample[selected_item][i][j]['sn'][l_l])
                        if len(temp_sn) != 0:
                            glycans[i]['sn'] = max(temp_sn)
                elif mode == 'compare_glycans' and min(treeview_levels) == 2:
                    temp_sn = []
                    for j_j, j in enumerate(glycans_per_sample[selected_item][i.split('_')[0]][i.split('_')[1]]['auc']):
                        glycans[i]['abundance'] += j
                        temp_sn.append(glycans_per_sample[selected_item][i.split('_')[0]][i.split('_')[1]]['sn'][j_j])
                    if len(temp_sn) != 0:
                        glycans[i]['sn'] = max(temp_sn)
                else:
                    for j_j, j in enumerate(glycans_per_sample[selected_item][i.split('_')[0]][i.split('_')[1]]['peaks']):
                        if float(i.split('_')[2]) == j:
                            glycans[i]['abundance'] = glycans_per_sample[selected_item][i.split('_')[0]][i.split('_')[1]]['auc'][j_j]
                            glycans[i]['sn'] = glycans_per_sample[selected_item][i.split('_')[0]][i.split('_')[1]]['sn'][j_j]
        
        elif mode == 'compare_samples':
            glycans = {}
            glycan = ''
            treeview_level = determine_treeview_level(chromatograms_list, chromatograms_list.selection())
            if treeview_level == 3:
                rt_text = chromatograms_list.item(chromatograms_list.selection(), "text")
                adduct = chromatograms_list.parent(chromatograms_list.selection())
                adduct_text = chromatograms_list.item(adduct, "text")
                glycan = chromatograms_list.parent(adduct)
                glycan_text = chromatograms_list.item(glycan, "text")
                glycan = glycan_text+'_'+adduct_text.split(' ')[0]+'_'+str(rt_text)
            if treeview_level == 2:
                adduct_text = chromatograms_list.item(chromatograms_list.selection(), "text")
                glycan = chromatograms_list.parent(chromatograms_list.selection())
                glycan_text = chromatograms_list.item(glycan, "text")
                glycan = glycan_text+'_'+adduct_text.split(' ')[0]
            if treeview_level == 1:
                glycan_text = chromatograms_list.item(chromatograms_list.selection(), "text")
                glycan = glycan_text
            
            for i in glycans_per_sample:
                if treeview_level == 1:
                    temp_abundance = 0
                    temp_sn = []
                    if glycan in glycans_per_sample[i].keys():
                        for j in glycans_per_sample[i][glycan]:
                            for k_k, k in enumerate(glycans_per_sample[i][glycan][j]['auc']):
                                temp_abundance += k
                                temp_sn.append(glycans_per_sample[i][glycan][j]['sn'][k_k])
                    glycans[i] = {'abundance': temp_abundance, 'sn': max(temp_sn) if len(temp_sn) > 0 else 0}
                if treeview_level == 2:
                    temp_abundance = 0
                    temp_sn = []
                    if glycan.split('_')[0] in glycans_per_sample[i].keys():
                        if glycan.split('_')[1] in glycans_per_sample[i][glycan.split('_')[0]].keys():
                            for j_j, j in enumerate(glycans_per_sample[i][glycan.split('_')[0]][glycan.split('_')[1]]['auc']):
                                temp_abundance += j
                                temp_sn.append(glycans_per_sample[i][glycan.split('_')[0]][glycan.split('_')[1]]['sn'][j_j])
                    glycans[i] = {'abundance': temp_abundance, 'sn': max(temp_sn) if len(temp_sn) > 0 else 0}
                if treeview_level == 3:
                    error_window('Select a glycan or adduct. Per peak comparison between samples graph plotting implementation is still work in progress.')
                    return
                    
        plot_window = tk.Toplevel()
        icon = ImageTk.PhotoImage(ico_image)
        plot_window.iconphoto(False, icon)
        plot_window.minsize(900, 900)
        plot_window.withdraw()
        plot_window.bind("<Configure>", on_resize)
        plot_window_title = glycan if mode == "compare_samples" else selected_item
        plot_window.title(f"Plot Graph - {plot_window_title}")
        plot_window.resizable(True, True)
        plot_window.protocol("WM_DELETE_WINDOW", exit_plot_window)
        
        #abundance plot
        abundance_graph_frame = tk.Frame(plot_window)
        abundance_graph_frame.pack(fill=tk.BOTH, expand=True)
        
        fig_plot_window = plt.figure(figsize=(0, 0))
        ax_plot_window = fig_plot_window.add_subplot(111)
        canvas_plot_window = FigureCanvasTkAgg(fig_plot_window, master=abundance_graph_frame)
        canvas_plot_window.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        fontsize = 100/len(glycans) if len(glycans) > 10 else 10
        
        glycans_keys = list(glycans.keys())
        for i_i, i in enumerate(glycans_keys):
            if len(i) > 10:
                glycans_keys[i_i] = i[:len(i)//2]+"\n"+i[len(i)//2:]
        
        diagonal = False
        if fontsize < 5:
            fontsize = 5
            diagonal = True
            
        zipped_plot_data = zip([glycans[i]['abundance'] for i in glycans], glycans_keys)
        zipped_plot_data = sorted(zipped_plot_data, reverse = True)
        plot_data_y_abundance, plot_data_x_abundance = zip(*zipped_plot_data)
        
        bar_graph = ax_plot_window.bar(plot_data_x_abundance, plot_data_y_abundance, color = 'black')
        ax_plot_window.set_ylabel('Abundance (AUC)')
        ax_plot_window.set_yscale('log')
        ax_plot_window.set_xticks(range(len(glycans)))
        ax_plot_window.set_xticklabels(plot_data_x_abundance, fontdict={'fontsize': fontsize})
        ax_plot_window_title = glycan if mode == "compare_samples" else selected_item
        ax_plot_window.set_title(f"{ax_plot_window_title}")
        
        if diagonal:
            plt.setp(ax_plot_window.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        tooltip_plot_window = ax_plot_window.annotate('', xy=(0, 0), xytext=(10, -20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.9), clip_on=True)
        
        canvas_plot_window.mpl_connect('motion_notify_event', lambda event: on_hover_column_graph(event, canvas_plot_window, ax_plot_window, bar_graph, plot_data_x_abundance, tooltip_plot_window))
        canvas_plot_window.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_plot_window, canvas_plot_window, False) if event.button == 3 else None)
        abundance_graph_frame.bind("<Configure>", lambda event: adjust_subplot_size_plot_window(event, ax_plot_window, 'abundance'))
        
        #signal-to-noise plot
        sn_graph_frame = tk.Frame(plot_window)
        sn_graph_frame.pack(fill=tk.BOTH, expand=True)
        
        fig_plot_window1 = plt.figure(figsize=(0, 0))
        ax_plot_window1 = fig_plot_window1.add_subplot(111)
        canvas_plot_window1 = FigureCanvasTkAgg(fig_plot_window1, master=sn_graph_frame)
        canvas_plot_window1.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            
        zipped_plot_data = zip([glycans[i]['sn'] for i in glycans], glycans_keys)
        zipped_plot_data = sorted(zipped_plot_data, reverse = True)
        plot_data_y_sn, plot_data_x_sn = zip(*zipped_plot_data)
        
        bar_graph1 = ax_plot_window1.bar(plot_data_x_sn, plot_data_y_sn, color = 'blue')
        ax_plot_window1_xlabel = "Samples" if mode == "compare_samples" else "Glycans"
        ax_plot_window1.set_xlabel(f"{ax_plot_window1_xlabel}")
        ax_plot_window1.set_ylabel('Signal-to-Noise Ratio')
        ax_plot_window1.set_yscale('log')
        ax_plot_window1.set_xticks(range(len(glycans)))
        ax_plot_window1.set_xticklabels(plot_data_x_sn, fontdict={'fontsize': fontsize})
        
        if diagonal:
            plt.setp(ax_plot_window1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        
        tooltip_plot_window1 = ax_plot_window1.annotate('', xy=(0, 0), xytext=(10, -20), textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.9), clip_on=True)
        
        ax_plot_window.set_position([0.0944, 0.1209, 0.8909, 0.8252])
        ax_plot_window1.set_position([0.0944, 0.1209, 0.8909, 0.8252])
        
        canvas_plot_window1.mpl_connect('motion_notify_event', lambda event: on_hover_column_graph(event, canvas_plot_window1, ax_plot_window1, bar_graph1, plot_data_x_sn, tooltip_plot_window1))
        canvas_plot_window1.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_plot_window1, canvas_plot_window1, False) if event.button == 3 else None)
        sn_graph_frame.bind("<Configure>", lambda event: adjust_subplot_size_plot_window(event, ax_plot_window1, 'sn'))
        
        
        canvas_plot_window.draw()
        canvas_plot_window1.draw()
        
        copy_data_to_clipboard = ttk.Button(plot_window, text="Copy Data to Clipboard", style="small_button_style1.TButton", command=copy_graph_data_clipboard, state=tk.NORMAL)
        copy_data_to_clipboard.pack()
        ToolTip(copy_data_to_clipboard, "Copies the values used in this graph to the clipboard. You can paste somewhere else (for example, on Excel) to use your data as you wish.")
        
        plot_window.update_idletasks()
        plot_window.deiconify()
        window_width = plot_window.winfo_width()
        window_height = plot_window.winfo_height()
        screen_width = plot_window.winfo_screenwidth()
        screen_height = plot_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        plot_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
    def quick_check_window():
        global max_spectrum_window, max_spectrum_opened, maximum_spectra, current_data, selected_item, samples_names
            
        def process_maximum_spectrum():
            def calculate_maximum_spectrum():
                global maximum_spectra, ax_mis, canvas_mis, coordinate_label_mis, type_coordinate_mis, selected_item
                
                if selected_item in maximum_spectra.keys():
                    maximum_spectrum = maximum_spectra[selected_item]
                else:
                    maximum = defaultdict(list)
                    
                    cpu_number = (os.cpu_count())-2 if os.cpu_count() < 60 else 60
                    if cpu_number <= 0:
                        cpu_number = 1
                        
                    interval_len = float(current_data['access'][-1]['scanList']['scan'][0]['scan start time'])//cpu_number if current_data['file_type'] == 'mzml' else float(current_data['access'][-1]['retentionTime'])//cpu_number
                    indexes = []

                    for i in range(cpu_number):
                        indexes.append(current_data['access'].time[i*interval_len]['index'] if current_data['file_type'] == 'mzml' else int(current_data['access'].time[i*interval_len]['num'])-1)
                    indexes.append(current_data['access'][-1]['index']+1 if current_data['file_type'] == 'mzml' else int(current_data['access'][-1]['num']))

                    results = []
                    with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_number) as executor:
                        for i_i, i in enumerate(indexes):
                            if i_i == len(indexes)-1:
                                break
                            min_max_index = [i, indexes[i_i+1]]
                            result = executor.submit(analyze_fraction, current_data['access'], 'ms level' if current_data['file_type'] == 'mzml' else 'msLevel', min_max_index)
                            results.append(result)
                            
                        for index, i in enumerate(results):
                            result = i.result()
                            for i in result:
                                maximum[i].extend(result[i])
                            results[index] = ''

                    for i in maximum:
                        maximum[i] = max(maximum[i])
                        
                    sorted_maximum = {key: maximum[key] for key in sorted(maximum)}
                    maximum = sorted_maximum
                    
                    maximum_spectrum = maximum
                    maximum_spectra[selected_item] = maximum_spectrum
                    
                x_data_mis = []
                y_data_mis = []
                for x in maximum_spectrum:
                    x_data_mis.extend([x - 0.000001, x, x + 0.000001])
                    y_data_mis.extend([0, maximum_spectrum[x], 0])
                    
                ax_mis.plot(x_data_mis, y_data_mis, marker='None', linewidth=1, color='black', label = 'MIS')
                
                annotate_top_y_values(ax_mis, canvas_mis)
                
                global og_x_range_mis, og_y_range_mis
                og_x_range_mis = [x_data_mis[0], x_data_mis[-1]]
                og_y_range_mis = [-1, max(y_data_mis)*1.1]
                ax_mis.set_xlim(og_x_range_mis)
                ax_mis.set_ylim(og_y_range_mis)
                
                canvas_mis.draw()
                
                on_plot_hover_motion_mis = canvas_mis.mpl_connect('motion_notify_event', lambda event: on_plot_hover(event, ax_mis, canvas_mis, coordinate_label_mis, type_coordinate_mis))
                
                processing_max_spectrum.destroy()
                max_spectrum_window.deiconify()
            
            def wait_thread():
                calculate_maximum_spectrum()
            
            global processing_max_spectrum
            processing_max_spectrum = tk.Toplevel()
            # processing_max_spectrum.attributes("-topmost", True)
            processing_max_spectrum.withdraw()
            processing_max_spectrum.title("Processing Maximum Intensity Spectrum")
            icon = ImageTk.PhotoImage(ico_image)
            processing_max_spectrum.iconphoto(False, icon)
            processing_max_spectrum.resizable(False, False)
            processing_max_spectrum.grab_set()
            processing_max_spectrum.protocol("WM_DELETE_WINDOW", on_closing)
            
            processing_max_spectrum_label = ttk.Label(processing_max_spectrum, text="Processing the Maximum Intensity Spectrum\nfor this sample, please wait.", font=("Segoe UI", list_font_size))
            processing_max_spectrum_label.pack(pady=35, padx=70)
            
            processing_max_spectrum.update_idletasks()
            processing_max_spectrum.deiconify()
            window_width = processing_max_spectrum.winfo_width()
            window_height = processing_max_spectrum.winfo_height()
            screen_width = processing_max_spectrum.winfo_screenwidth()
            screen_height = processing_max_spectrum.winfo_screenheight()
            x_position = (screen_width - window_width) // 2
            y_position = (screen_height - window_height) // 2
            processing_max_spectrum.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
            
            t = threading.Thread(target=wait_thread)
            t.start()
            
        def exit_window_qcw():
            global max_spectrum_opened, glycans_list_quickcheck_save
            values_list = []
            for item_id in glycans_list.get_children():
                values = glycans_list.item(item_id, "values")
                values_list.append(values)
            library_name_from_path = library_path.split("/")[-1]
            glycans_list_quickcheck_save[f"{selected_item}_{tolerance}_{library_name_from_path}"] = values_list
            max_spectrum_opened = False
            max_spectrum_window.destroy()
        
        def glycans_list_sort(tv, col, reverse):
            # Get the data in the specified column
            data = [(tv.set(k, col), k) for k in tv.get_children('')]
            
            # Try to convert data to float if possible for numerical sorting
            try:
                data = [(float(item[0]), item[1]) for item in data]
            except ValueError:
                pass
            
            # Sort data
            data.sort(reverse=reverse)
            
            # Rearrange items in sorted positions
            for index, (val, k) in enumerate(data):
                tv.move(k, '', index)
            
            # Reverse sort next time
            tv.heading(col, command=lambda: glycans_list_sort(tv, col, not reverse))
            
        def copy_selected_rows(event):
            selected_items = glycans_list.selection()
            if len(selected_items) == 0:
                return
                
            rows_data = []

            # Retrieve the data from selected rows
            for item in selected_items:
                row_values = glycans_list.item(item, "values")
                rows_data.append('\t'.join(row_values))
            
            # Join the rows with newline character
            clipboard_data = '\n'.join(rows_data)
            
            # Copy to clipboard
            main_window.clipboard_clear()
            main_window.clipboard_append(clipboard_data)
        
        def click_glycans_list(event):
            global glycans_list_quickcheck, rectangles
            key = event.keysym
            
            if key == "Escape": #removes selection and clears zoom and rectangles
                selected_items = glycans_list.selection()
                glycans_list.selection_remove(selected_items)
                ax_mis.set_xlim(og_x_range_mis)
                ax_mis.set_ylim(og_y_range_mis)
                selected_glycan_list = ''
                canvas_mis.draw_idle()
            else:
                selected_glycan_list = glycans_list.focus()
            
            for rect in rectangles:
                rect.remove()
            rectangles.clear()
            
            if selected_glycan_list == '':
                return
                
            selected_glycan_list_content = glycans_list.item(selected_glycan_list, "values")
            selected_glycan_list_content = (selected_glycan_list_content[0], selected_glycan_list_content[1], float(selected_glycan_list_content[2]), float(selected_glycan_list_content[3]), float(selected_glycan_list_content[4]))
            
            mz_list = []
            library_name_from_path = library_path.split("/")[-1]
            for i in glycans_list_quickcheck[f"{selected_item}_{tolerance}_{library_name_from_path}"]:
                if i[0] == selected_glycan_list_content:
                    mz_list = i[1]
                    max_int = i[2]
                    break
            
            mz_array = list(maximum_spectra[selected_item].keys())
            int_array = list(maximum_spectra[selected_item].values())
            
            ax_mis.set_xlim(mz_list[0]-2, mz_list[-1]+2)
            ax_mis.set_ylim(0, max_int*1.1)
            
            adduct_charge = General_Functions.determine_adduct_charge('isotopic_spacing', mz_list)
            for index, i in enumerate(mz_list):
                current_peak_found_mz = selected_glycan_list_content[2]+(General_Functions.h_mass/abs(adduct_charge))*index
                calculated_width = General_Functions.tolerance_calc(tolerance[0], tolerance[1], current_peak_found_mz)
                rectangles.append(ax_mis.add_patch(Rectangle((current_peak_found_mz-calculated_width, ax_mis.get_ylim()[0]), (current_peak_found_mz+calculated_width) - (current_peak_found_mz-calculated_width), 1000000000000, color='#FEB7A1', alpha=0.3)))
                rectangles.append(ax_mis.add_patch(Rectangle((i-0.001, ax_mis.get_ylim()[0]), (i+0.001) - (i-0.001), 1000000000000, color='#345eeb', alpha=0.3)))
                
            canvas_mis.draw_idle()
            
        def analyze_glycan(mz_array, int_array, glycan_info, target_mz, tolerance, max_charges, adduct_charge):
            mz_array_len = len(mz_array)-1
            if target_mz > mz_array[-1]:
                mz_id = -1
            else:
                mz_id = General_Functions.binary_search_with_tolerance(mz_array, target_mz, 0, mz_array_len, General_Functions.tolerance_calc(tolerance[0], tolerance[1], target_mz), int_array)
                
            if mz_id != -1:
                found_mz = mz_array[mz_id]
                charge_range = range(1, abs(max_charges)*2)
                mono_int = int_array[mz_id]
                ppm_error = General_Functions.calculate_ppm_diff(mz_array[mz_id], target_mz)
                iso_actual = [1]
                bad = False
                margin = 0.0
                max_int = int_array[mz_id]
                temp_id = General_Functions.binary_search_with_tolerance(mz_array, found_mz+(General_Functions.h_mass/abs(adduct_charge)), mz_id, mz_array_len, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz+(General_Functions.h_mass/abs(adduct_charge))), int_array)
                
                if temp_id == -1:
                    bad = True
                
                if not bad:
                    for i in charge_range: #check if it's monoisotopic and correct charge
                        temp_id = General_Functions.binary_search_with_tolerance(mz_array, found_mz-(General_Functions.h_mass/i), 0, mz_id, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz-(General_Functions.h_mass/i)), int_array) #check monoisotopic
                        if temp_id != -1 and int_array[temp_id] > 0:
                            expected_value = (mz_array[temp_id]*i*0.0006)+0.1401
                            if (mono_int/int_array[temp_id] < expected_value*(1+(margin*2))):
                                bad = True
                                break                   
                        if i == 1 or i == abs(adduct_charge) or (i == 2 and abs(adduct_charge) == 4) or (i == 3 and abs(adduct_charge) == 6): #ignores charge 1 due to the fact that any charge distribution will find a hit on that one
                            continue
                        temp_id = General_Functions.binary_search_with_tolerance(mz_array, found_mz+(General_Functions.h_mass/i), mz_id, mz_array_len, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz+(General_Functions.h_mass/i)), int_array) #check for correct charge
                        if temp_id != -1:
                            expected_value = (target_mz*i*0.0006)+0.1401
                            if (int_array[temp_id]/mono_int > expected_value*(1-margin)):
                                bad = True
                                break
                if not bad:
                    isos_found = 0
                    mz_isos = []
                    for i_i, i in enumerate(glycan_info['Isotopic_Distribution_Masses']): #check isotopic peaks and add to the intensity
                        if i_i == 0: #ignores monoisotopic this time around
                            continue
                        temp_id = General_Functions.binary_search_with_tolerance(mz_array, found_mz+(i_i*(General_Functions.h_mass/abs(adduct_charge))), mz_id, mz_array_len, General_Functions.tolerance_calc(tolerance[0], tolerance[1], found_mz+(i_i*(General_Functions.h_mass/abs(adduct_charge)))), int_array)
                        if temp_id != -1 and int_array[temp_id] > 0:
                            isos_found += 1
                            mz_isos.append(mz_array[temp_id])
                            iso_actual.append(int_array[temp_id]/mono_int)
                            if isos_found < 3:
                                max_int = max(max_int, int_array[temp_id])
                        else:
                            if isos_found == 0: #a compound needs at least 2 identifiable peaks (monoisotopic + 1 from isotopic envelope)
                                bad = True
                                break
                            else:
                                break                
                if not bad and (iso_actual[1] < 0.2 or iso_actual[1] > 5): #this should avoid situations where it's obvious that it's picking the wrong charge because the second peak is almost invisible compared to the third and first, which when z=2 means that it's very likely actually a singly charge compound, for example
                    if len(iso_actual) > 2:
                        if iso_actual[2] > iso_actual[1]*10 or iso_actual[1] > 5:
                            bad = True
                    else:
                        if iso_actual[1] < 0.2 or iso_actual[1] > 5: #smallest glycan should have the second iso_actual somewhere around 0.5, so a cutoff lower than that is fine
                            bad = True
                
                if bad:
                    return "bad"
                else:
                    iso_target = glycan_info['Isotopic_Distribution'][:isos_found+1]
                    
                    ratios = []
                    weights = []
                    for i_i, i in enumerate(iso_target):
                        if i_i == 0:
                            continue
                        intensities = [i, iso_actual[i_i]]
                        ratio = min(intensities)/max(intensities)
                        
                        #scales the score in a sigmoid, with steepness determine by k_value
                        k_value = 10
                        corrected_ratio = 1 / (1 + np.exp(-k_value * (ratio - 0.5)))
                        
                        ratios.append(corrected_ratio)
                        weights.append(1/(np.exp(1.25*i_i)))
                    
                    iso_quali = np.average(ratios, weights = weights)
                
                    #reduces score if fewer isotopic peaks are found: punishing for only 1 peaks, normal score from 2 and over (besides the monoisotopic)
                    if len(iso_actual) == 2:
                        iso_quali = (iso_quali*0.8)
                            
                    return [found_mz]+mz_isos, iso_quali, ppm_error, max_int
            else:
                return "bad"
                
        def quick_check_glycans():
            global min_max_monos, min_max_hex, min_max_hexnac, min_max_xyl, min_max_fuc, min_max_sia, min_max_ac, min_max_ac, min_max_gc, min_max_hn, min_max_ua, forced, min_max_proton, max_charges, reducing_end_tag, internal_standard, permethylated, lactonized_ethyl_esterified, min_max_sulfation, min_max_phosphorylation, reduced, fast_iso, high_res, glycans_list_quickcheck
            
            if len(library_path) == 0:
                error_window("You must first generate or import a library\nto do a quick check.")
                max_spectrum_window.grab_set()
                return
            library_name_from_path = library_path.split("/")[-1]
            if f"{selected_item}_{tolerance}_{library_name_from_path}" in glycans_list_quickcheck.keys():
                glycans_list_temp = glycans_list_quickcheck[f"{selected_item}_{tolerance}_{library_name_from_path}"]
            else:
                try:
                    with open(library_path, 'rb') as f:
                        library_data = dill.load(f)
                        f.close()
                except:
                    error_window("Can't load the library. The library file you are trying to load was created in a version older than 0.2.0 or is corrupted. Generate a new library or try a different file.")
                    return
                    
                full_library = library_data[0]
                library_metadata = library_data[1]
                
                min_max_monos = library_metadata[0]
                min_max_hex = library_metadata[1]
                min_max_hexnac = library_metadata[2]
                min_max_fuc = library_metadata[3]
                min_max_sia = library_metadata[4]
                min_max_ac = library_metadata[5]
                min_max_gc = library_metadata[6]
                forced = library_metadata[7]
                min_max_proton = library_metadata[8]
                max_charges = library_metadata[9]
                reducing_end_tag = library_metadata[10]
                internal_standard = library_metadata[11]
                permethylated = library_metadata[12]
                lactonized_ethyl_esterified = library_metadata[13]
                reduced = library_metadata[14]
                fast_iso = library_metadata[15]
                high_res = library_metadata[16]
                min_max_xyl = library_metadata[18]
                min_max_hn = library_metadata[19]
                min_max_ua = library_metadata[20]
                min_max_sulfation = library_metadata[21]
                min_max_phosphorylation = library_metadata[22]
                
                glycans_list_temp = []
                
                mz_array = list(maximum_spectra[selected_item].keys())
                int_array = list(maximum_spectra[selected_item].values())
                
                for i in full_library:
                    for j in full_library[i]['Adducts_mz']:
                        adduct_charge = General_Functions.determine_adduct_charge('neutral_mass', full_library[i]['Neutral_Mass'], full_library[i]['Adducts_mz'][j], max_charges)
                        result = analyze_glycan(mz_array, int_array, full_library[i], full_library[i]['Adducts_mz'][j], tolerance, max_charges, adduct_charge)
                        if result != 'bad':
                            glycans_list_temp.append([(i, j, float("%.4f" % round(full_library[i]['Adducts_mz'][j], 4)), float("%.2f" % round(result[1], 2)), float("%.1f" % round(result[2], 1))), result[0], result[3]])
                            
                library_name_from_path = library_path.split("/")[-1]
                glycans_list_quickcheck[f"{selected_item}_{tolerance}_{library_name_from_path}"] = glycans_list_temp
                    
            for item in glycans_list.get_children():
                glycans_list.delete(item)        
            for i in glycans_list_temp:
                glycans_list.insert("", "end", values=i[0])
                
        if selected_item not in samples_names:
            error_window("Raw spectra file is not loaded for this sample. Load it in the 'Select Files' menu and try again.")
            return
        
        if max_spectrum_opened:
            max_spectrum_window.deiconify()
            max_spectrum_window.lift()
            max_spectrum_window.focus_set()
            return
        
        max_spectrum_opened = True
        
        max_spectrum_window = tk.Toplevel()
        icon = ImageTk.PhotoImage(ico_image)
        max_spectrum_window.iconphoto(False, icon)
        max_spectrum_window.withdraw()
        max_spectrum_window.minsize(900, 480)
        max_spectrum_window.bind("<Configure>", on_resize)
        max_spectrum_window.bind("<KeyPress>", on_key_press)
        max_spectrum_window.bind("<KeyRelease>", on_key_release)
        max_spectrum_window.title(f"Maximum Intensity Spectrum - {selected_item}")
        max_spectrum_window.resizable(True, True)
        max_spectrum_window.grid_rowconfigure(0, weight=0)
        max_spectrum_window.grid_rowconfigure(1, weight=1)
        max_spectrum_window.grid_columnconfigure(0, weight=0)
        max_spectrum_window.grid_columnconfigure(1, weight=1)
        max_spectrum_window.grab_set()
        max_spectrum_window.protocol("WM_DELETE_WINDOW", exit_window_qcw)
        
        global quick_analysis_button
        quick_analysis_button = ttk.Button(max_spectrum_window, text="Quick Analysis", style="small_button_sfw_style1.TButton", command=quick_check_glycans)
        quick_analysis_button.grid(row=0, column=1, padx=10, pady=10, sticky="ne")
        ToolTip(quick_analysis_button, "If you have a library loaded, you can click here to get a quick glimpse of the presence of glycans in your sample.")
        
        glycans_list_scrollbar = tk.Scrollbar(max_spectrum_window, orient=tk.VERTICAL)
        glycans_list = ttk.Treeview(max_spectrum_window, columns=("Glycan", "Adduct", "m/z", "Iso. Fit. Score", "PPM error"), height=25, yscrollcommand=glycans_list_scrollbar.set, show='headings')
        
        glycans_list_columns = ["Glycan", "Adduct", "m/z", "Iso. Fit. Score", "PPM error"]
        for col in glycans_list_columns:
            glycans_list.heading(col, text=col, command=lambda _col=col: glycans_list_sort(glycans_list, _col, False))
        
        glycans_list.column("Glycan", width=75)
        glycans_list.column("Adduct", width=50)
        glycans_list.column("m/z", width=70)
        glycans_list.column("Iso. Fit. Score", width=80)
        glycans_list.column("PPM error", width=70)
        glycans_list_scrollbar.config(command=glycans_list.yview, width=10)
        glycans_list.grid(row=1, column=0, padx=(10, 20), pady=10, sticky="nsew")
        glycans_list_scrollbar.grid(row=1, column=0, padx=10, pady=10, sticky="nse")
        
        glycans_list.bind("<ButtonRelease-1>", click_glycans_list)
        glycans_list.bind("<KeyRelease-Up>", click_glycans_list)
        glycans_list.bind("<KeyRelease-Down>", click_glycans_list)
        max_spectrum_window.bind("<Control-c>", copy_selected_rows)
        max_spectrum_window.bind("<KeyRelease-Escape>", click_glycans_list)
        
        library_name_from_path = library_path.split("/")[-1]
        if f"{selected_item}_{tolerance}_{library_name_from_path}" in glycans_list_quickcheck_save.keys():
            for i_i, i in enumerate(glycans_list_quickcheck_save[f"{selected_item}_{tolerance}_{library_name_from_path}"]):
                glycans_list.insert("", "end", values=i)
    
        max_spectrum_plot_frame = ttk.Labelframe(max_spectrum_window, text="Maximum Intensity Spectrum", style="chromatogram.TLabelframe")
        max_spectrum_plot_frame.grid(row=1, column=1, padx=10, pady=10, sticky="nsew")

        global canvas_mis, ax_mis, coordinate_label_mis, type_coordinate_mis
        fig_mis = plt.figure(figsize=(0, 0))
        ax_mis = fig_mis.add_subplot(111)
        canvas_mis = FigureCanvasTkAgg(fig_mis, master=max_spectrum_plot_frame)
        canvas_mis.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        ax_mis.set_xlabel('m/z')
        ax_mis.set_ylabel('Intensity (AU)')
        
        global rectangles
        rectangles = []
        
        coordinate_label_mis = tk.Label(max_spectrum_plot_frame, text="", anchor="e", font=("Segoe UI", 8), bg="white")
        coordinate_label_mis.place(relx=1.0, rely=0, anchor='ne')
        type_coordinate_mis = "spectra"
        coordinate_label_mis.lift()
        
        frame_width = max_spectrum_plot_frame.winfo_width()
        frame_height = max_spectrum_plot_frame.winfo_height()

        left_margin = 65
        right_margin = 10
        top_margin = 20
        bottom_margin = 45

        subplot_width = (frame_width - left_margin - right_margin) / frame_width
        subplot_height = (frame_height - top_margin - bottom_margin) / frame_height
        subplot_left = left_margin / frame_width
        subplot_bottom = bottom_margin / frame_height

        ax_mis.set_position([subplot_left, subplot_bottom, subplot_width, subplot_height])
        
        max_spectrum_plot_frame.bind("<Configure>", lambda event, ax=ax_mis: adjust_subplot_size(event, ax_mis, canvas_mis))
        canvas_mis.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_mis, canvas_mis, True) if event.button == 3 else None) 
        zoom_selection_key_press_mis = canvas_mis.mpl_connect('key_press_event', lambda event: zoom_selection_compare(event, ax_mis, canvas_mis, type_coordinate_mis))
        zoom_selection_key_release_mis = canvas_mis.mpl_connect('key_release_event', lambda event: zoom_selection_compare(event, ax_mis, canvas_mis, type_coordinate_mis))
        zoom_selection_motion_notify_mis = canvas_mis.mpl_connect('motion_notify_event', lambda event: zoom_selection_compare(event, ax_mis, canvas_mis, type_coordinate_mis))
        zoom_selection_button_press_mis = canvas_mis.mpl_connect('button_press_event', lambda event: zoom_selection_compare(event, ax_mis, canvas_mis, type_coordinate_mis) if event.button == 1 else None)
        zoom_selection_button_release_mis = canvas_mis.mpl_connect('button_release_event', lambda event: zoom_selection_compare(event, ax_mis, canvas_mis, type_coordinate_mis) if event.button == 1 else None)
        on_scroll_event_mis = canvas_mis.mpl_connect('scroll_event', lambda event: on_scroll(event, ax_mis, canvas_mis, type_coordinate_mis))
        on_double_click_event_mis = canvas_mis.mpl_connect('button_press_event', lambda event: on_double_click(event, ax_mis, canvas_mis, og_x_range_mis, og_y_range_mis, type_coordinate_mis) if event.button == 1 else None)
        on_pan_press_mis = canvas_mis.mpl_connect('button_press_event', lambda event: on_pan(event, ax_mis, canvas_mis, type_coordinate_mis) if event.button == 1 else None)
        on_pan_release_mis = canvas_mis.mpl_connect('button_release_event', lambda event: on_pan(event, ax_mis, canvas_mis, type_coordinate_mis) if event.button == 1 else None)
        on_pan_motion_mis = canvas_mis.mpl_connect('motion_notify_event', lambda event: on_pan(event, ax_mis, canvas_mis, type_coordinate_mis))
        on_pan_right_click_press_mis = canvas_mis.mpl_connect('button_press_event', lambda event: on_pan_right_click(event, ax_mis, canvas_mis, type_coordinate_mis) if event.button == 3 else None)
        on_pan_right_click_release_mis = canvas_mis.mpl_connect('button_release_event', lambda event: on_pan_right_click(event, ax_mis, canvas_mis, type_coordinate_mis) if event.button == 3 else None)
        on_pan_right_click_motion_mis = canvas_mis.mpl_connect('motion_notify_event', lambda event: on_pan_right_click(event, ax_mis, canvas_mis, type_coordinate_mis))
        
        max_spectrum_window.update_idletasks()
        window_width = max_spectrum_window.winfo_width()
        window_height = max_spectrum_window.winfo_height()
        screen_width = max_spectrum_window.winfo_screenwidth()
        screen_height = max_spectrum_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        max_spectrum_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        process_maximum_spectrum()
    
    def run_quick_trace_window():
        global quick_trace_opened, quick_trace_window, colors, current_data, selected_item, quick_traces_all, quick_traces_list_save, samples_list
                    
        def trace_loading():
            global loading_eic_window
            loading_eic_window = tk.Toplevel()
            loading_eic_window.withdraw()
            loading_eic_window.title("Tracing the EIC/EIE")
            icon = ImageTk.PhotoImage(ico_image)
            loading_eic_window.iconphoto(False, icon)
            loading_eic_window.resizable(False, False)
            loading_eic_window.grab_set()
            loading_eic_window.protocol("WM_DELETE_WINDOW", on_closing)
            
            loading_eic_window_label = ttk.Label(loading_eic_window, text="Tracing EIC/EIE for the specified m/z in this sample.\nPlease wait.", font=("Segoe UI", list_font_size))
            loading_eic_window_label.pack(pady=35, padx=70)
            
            loading_eic_window.update_idletasks()
            loading_eic_window.deiconify()
            window_width = loading_eic_window.winfo_width()
            window_height = loading_eic_window.winfo_height()
            screen_width = loading_eic_window.winfo_screenwidth()
            screen_height = loading_eic_window.winfo_screenheight()
            x_position = (screen_width - window_width) // 2
            y_position = (screen_height - window_height) // 2
            loading_eic_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        def trace_mz():
            global loading_eic_window, quick_traces_all
            selected_items_qtl = quick_trace_list.selection()
            if len(samples_list) == 0 or len(selected_items_qtl) == 0:
                loading_eic_window.destroy()
                quick_trace_window.lift()
                quick_trace_list.focus_set()
                return
                
            if len(selected_items_qtl) > 1:
                clear_plot(ax, canvas)
                
            for i_i, i in enumerate(selected_items_qtl):
                values_qtl = quick_trace_list.item(i, "values")
                target_mz = values_qtl[0]
                tolerance = float(values_qtl[1])
                color_to_trace = quick_trace_list.item(i, 'tags')[0]
                
                rt_array = current_data['rt_array']
                if current_data['time_unit'] == 'seconds':
                    rt_array = [x/60 for x in rt_array]
                if f"{selected_item}_{target_mz}_{tolerance}" in quick_traces_all.keys():
                    int_array = quick_traces_all[f"{selected_item}_{target_mz}_{tolerance}"]
                else:
                    targets = target_mz.split(";")
                    int_array = []
                    if current_data['file_type'] == 'mzml':
                        for i in current_data['ms1_array']:
                            if current_data['access'][i]['ms level'] == 1:
                                int_array.append(0)
                                for j in targets:
                                    temp_id = General_Functions.binary_search_with_tolerance(current_data['access'][i]['m/z array'], float(j), 0, len(current_data['access'][i]['m/z array'])-1, tolerance, current_data['access'][i]['intensity array'])
                                    if temp_id != -1:
                                        int_array[-1] += current_data['access'][i]['intensity array'][temp_id]
                    elif current_data['file_type'] == 'mzxml':
                        for i in current_data['ms1_array']:
                            if current_data['access'][i]['msLevel'] == 1:
                                int_array.append(0)
                                for j in targets:
                                    temp_id = General_Functions.binary_search_with_tolerance(current_data['access'][i]['m/z array'], float(j), 0, len(current_data['access'][i]['m/z array'])-1, tolerance, current_data['access'][i]['intensity array'])
                                    if temp_id != -1:
                                        int_array[-1] += current_data['access'][i]['intensity array'][temp_id]
                    quick_traces_all[f"{selected_item}_{target_mz}_{tolerance}"] = int_array
                
                if len(quick_trace_list.selection()) > 1:
                    show_graph([rt_array, int_array, color_to_trace, f"{target_mz}±{tolerance}"], clear = False)
                    if i_i >= len(selected_items_qtl)-1:
                        loading_eic_window.destroy()
                        quick_trace_window.lift()
                        quick_trace_list.focus_set()
                    
                else:
                    show_graph([rt_array, int_array, color_to_trace, f"{target_mz}±{tolerance}"])
                    if i_i >= len(selected_items_qtl)-1:
                        loading_eic_window.destroy()
                        quick_trace_window.lift()
                        quick_trace_list.focus_set()
    
        def on_quick_trace_list_motion(event):
            region = quick_trace_list.identify_region(event.x, event.y)
            column = quick_trace_list.identify_column(event.x)
            item = quick_trace_list.identify_row(event.y)
            
            if region == "tree" or (region == "cell" and column == "#3"):
                values = quick_trace_list.item(item, "values")
                quick_trace_list.config(cursor="hand2")
            else:
                quick_trace_list.config(cursor="")
        
        def exit_window_qt():
            global quick_trace_opened, quick_traces_list_save
            values_list = []
            tags_list = []
            for item_id in quick_trace_list.get_children():
                values = quick_trace_list.item(item_id, "values")
                tags = quick_trace_list.item(item_id, "tags")
                values_list.append(values)
                tags_list.append(tags)
            quick_traces_list_save = [values_list, tags_list]
            quick_trace_opened = False
            quick_trace_window.destroy()
            
        def draw_rectangle(color):
            # Create a small green rectangle image (16x16 pixels)
            size = (32, 16)
            rect_image = Image.new("RGBA", size, (255, 255, 255, 0))  # Transparent background
            draw = ImageDraw.Draw(rect_image)
            draw.rectangle((0, 0, 32, 16), fill=color)  # Draw green rectangle
            photo = ImageTk.PhotoImage(rect_image)
            
            return photo

        def click_quick_trace_list(event):
            region_qtl = quick_trace_list.identify_region(event.x, event.y)
            column_qtl = quick_trace_list.identify_column(event.x)
            item_qtl = quick_trace_list.identify_row(event.y)
            values_qtl = quick_trace_list.item(item_qtl, "values")
            
            if region_qtl == "tree":
                color_code = colorchooser.askcolor(title="Choose a color")[1]
                quick_trace_window.lift()
                quick_trace_window.focus_set()
                if color_code:
                    # quick_trace_list.tag_configure(color_code, foreground=color_code)
                    row_image = draw_rectangle(color_code)
                    quick_trace_list.item(item_qtl, tags=(color_code,), image=row_image)
                    quick_trace_list_image_refs.append(row_image)
                    
            elif region_qtl == "cell" and column_qtl == "#3":
                quick_trace_list.delete(item_qtl)
            
            elif region_qtl == "cell":
                trace_loading()
                threading.Thread(target = trace_mz).start()
                    
        
        def add_qtl():
            if mz_entry.get().lower() == 'mis':
                library_name_from_path = library_path.split("/")[-1]
                if f"{selected_item}_{tolerance}_{library_name_from_path}" in glycans_list_quickcheck_save.keys():
                    for i_i, i in enumerate(glycans_list_quickcheck_save[f"{selected_item}_{tolerance}_{library_name_from_path}"]):
                        random_color = random.choice(colors)
                        # quick_trace_list.tag_configure(random_color, foreground=random_color)
                        row_image = draw_rectangle(random_color)
                        quick_trace_list.insert("", "end", text='', image=row_image, values=(f"{i[2]}", f"{round(General_Functions.tolerance_calc(tolerance[0], tolerance[1], mz = 1000.0), 4)}", "❌"), tags=(random_color,))
                        quick_trace_list_image_refs.append(row_image)
                else:
                    error_window(f"No MIS analyzed for this sample.")
                    quick_trace_window.lift()
                    quick_trace_window.focus_set()
                    return
            else:
                try:
                    mz_list = mz_entry.get().split(";")
                    for i in mz_list:
                        if len(i) > 0:
                            float(i)
                except:
                    error_window(f"Invalid value on m/z field. Correct it and try again.")
                    quick_trace_window.lift()
                    quick_trace_window.focus_set()
                    return
                try:
                    float(tol_entry.get())
                except:
                    error_window(f"Invalid value on tolerance field. Correct it and try again.")
                    quick_trace_window.lift()
                    quick_trace_window.focus_set()
                    return
                random_color = random.choice(colors)
                # quick_trace_list.tag_configure(random_color, foreground=random_color)
                row_image = draw_rectangle(random_color)
                quick_trace_list.insert("", "end", text='', image=row_image, values=(f"{mz_entry.get()}", f"{tol_entry.get()}", "❌"), tags=(random_color,))
                quick_trace_list_image_refs.append(row_image)
                
        def save_qtraces_csv():
            # Create the save file dialog
            csv_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
            
            # If no file is selected, cancel the saving
            if len(csv_path) == 0:
                return
                
            # Initialize list of lines
            lines = []
            
            # Create header
            header = ['color_hex', 'mz', 'tolerance']
            
            # Add header to lines
            lines.append(",".join(header)+"\n")
            
            # Gather the rows from the list and add to lines
            for row in quick_trace_list.get_children():
                # Get the values from the row
                row_values = quick_trace_list.item(row, "values")
                mz_qtr_save = row_values[0]
                tol_qtr_save = row_values[1]
                
                # Get the color from the row
                row_color = quick_trace_list.item(row, "tags")[0]
                
                # Create a new line
                new_line = [row_color, mz_qtr_save, tol_qtr_save]
            
                # Add new_line to lines
                lines.append(",".join(new_line)+"\n")
                
            with open(csv_path, "w") as csv_file:
                for line in lines:
                    csv_file.write(line)
                
        def load_qtraces_csv():
            # Open the open file dialog
            csv_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
            
            # If no file is selected, return None
            if len(csv_path) == 0:
                return
                
            # Initialize the new quick traces list
            new_qtraces = []
            
            try:
                # Open the csv file
                with open(csv_path, "r") as csv_file:
                    for index, line in enumerate(csv_file):
                        # Ignore the header
                        if index == 0:
                            continue
                            
                        # Get the data from the row
                        color_load_qtr, mz_load_qtr, tol_load_qtr = line.strip().split(",")
                        
                        # Check if there's a color specified, if not, pick random color
                        if len(color_load_qtr) == 0:
                            color_load_qtr = random.choice(colors)
                            
                        # Test the mz and tol_entry
                        float(mz_load_qtr)
                        float(tol_load_qtr)
                            
                        # Add the line to the new_qtraces
                        new_qtraces.append(([mz_load_qtr, tol_load_qtr, "❌"], color_load_qtr))
            except:
                error_window("Couldn't load the quick traces from the file. Make sure it is formatted correctly.")
                quick_trace_window.deiconify()
                quick_trace_window.lift()
                quick_trace_window.focus_set()
                return
                
            for row_id in quick_trace_list.get_children():
                quick_trace_list.delete(row_id)
                
            for row_values, color_load_qtr in new_qtraces:
                # Create the tag for color
                # quick_trace_list.tag_configure(color_load_qtr, foreground=color_load_qtr)
                row_image = draw_rectangle(color_load_qtr)
                quick_trace_list.insert("", "end", text='', image=row_image, values=row_values, tags=(color_load_qtr,))
                quick_trace_list_image_refs.append(row_image)
                
            quick_trace_window.deiconify()
            quick_trace_window.lift()
            quick_trace_window.focus_set()
        
        if quick_trace_opened:
            quick_trace_window.deiconify()
            quick_trace_window.lift()
            quick_trace_window.focus_set()
            return
        
        quick_trace_opened = True
        
        global quick_trace_list_image_refs
        quick_trace_list_image_refs = []
        
        small_button_qtw_style1 = ttk.Style().configure("small_button_qtw_style1.TButton", font=("Segoe UI", list_font_size), relief="raised", padding = (0, 0), justify="center")
        
        quick_trace_window = tk.Toplevel()
        icon = ImageTk.PhotoImage(ico_image)
        quick_trace_window.iconphoto(False, icon)
        quick_trace_window.withdraw()
        quick_trace_window.minsize(380, 480)
        quick_trace_window.bind("<Configure>", on_resize)
        quick_trace_window.title(f"Quick Traces")
        quick_trace_window.resizable(False, True)
        quick_trace_window.grid_rowconfigure(0, weight=0)
        quick_trace_window.grid_rowconfigure(1, weight=1)
        quick_trace_window.grid_columnconfigure(0, weight=0)
        quick_trace_window.protocol("WM_DELETE_WINDOW", exit_window_qt)
    
        mz_label = ttk.Label(quick_trace_window, text='m/z:', font=("Segoe UI", list_font_size_smaller))
        mz_label.grid(row=0, column=0, padx=(10, 10), pady=(10, 0), sticky="wns")
        ToolTip(mz_label, "Type in the m/z value you want to trace or 'MIS' to import the Quick Analysis results from the MIS window.")
    
        mz_entry = ttk.Entry(quick_trace_window, width=20)
        mz_entry.grid(row=0, column=0, padx=(40, 10), pady=(10, 0), sticky='wns')
        ToolTip(mz_entry, "Type in the m/z value you want to trace or 'MIS' to import the Quick Analysis results from the MIS window..")
        mz_entry.bind("<KeyRelease-Return>", lambda event: (add_qtl(), mz_entry.focus_set(), mz_entry.selection_range(0, tk.END)))
    
        tol_label = ttk.Label(quick_trace_window, text='Tolerance (m/z):', font=("Segoe UI", list_font_size_smaller))
        tol_label.grid(row=0, column=0, padx=(170, 135), pady=(10, 0), sticky="ens")
        ToolTip(tol_label, "Type in the tolerance value for the traced m/z value.")
    
        tol_entry = ttk.Entry(quick_trace_window, width=5)
        tol_entry.grid(row=0, column=0, padx=(10, 100), pady=(10, 0), sticky='ens')
        ToolTip(tol_entry, "Type in the tolerance value for the traced m/z value.")
        tol_entry.bind("<KeyRelease-Return>", lambda event: (add_qtl(), mz_entry.focus_set(), mz_entry.selection_range(0, tk.END)))
        
        add_button = ttk.Button(quick_trace_window, text="Add Trace", style="small_button_qtw_style1.TButton", command=add_qtl)
        add_button.grid(row=0, column=0, padx=(10, 10), pady=(10,0), sticky="ens")
        
        quick_trace_list_columns = ["m/z", "Tolerance", "Remove"]
        quick_trace_list_scrollbar = tk.Scrollbar(quick_trace_window, orient=tk.VERTICAL)
        quick_trace_list = ttk.Treeview(quick_trace_window, columns=quick_trace_list_columns, height=25, yscrollcommand=quick_trace_list_scrollbar.set, show='tree headings')
        
        for col in quick_trace_list_columns:
            quick_trace_list.heading(col, text=col)
            
        quick_trace_list.heading("#0", text="Color")
        quick_trace_list.column("#0", width=10, anchor='w')
        quick_trace_list.column("m/z", width=80)
        quick_trace_list.column("Tolerance", width=10)
        quick_trace_list.column("Remove", width=5, anchor='center')
        quick_trace_list_scrollbar.config(command=quick_trace_list.yview, width=10)
        quick_trace_list.grid(row=1, column=0, padx=(10, 20), pady=10, sticky="nsew")
        quick_trace_list_scrollbar.grid(row=1, column=0, padx=10, pady=10, sticky="nse")
        
        if len(quick_traces_list_save) > 0:
            for i_i, i in enumerate(quick_traces_list_save[0]):
                quick_trace_list.tag_configure(quick_traces_list_save[1][i_i][0], foreground=quick_traces_list_save[1][i_i][0])
                quick_trace_list.insert("", "end", values=i, tags=quick_traces_list_save[1][i_i])
        
        save_qtraces_button = ttk.Button(quick_trace_window, text="Save Traces", style="small_button_qtw_style1.TButton", command=save_qtraces_csv)
        save_qtraces_button.grid(row=2, column=0, padx=(10, 10), pady=(0,10), sticky="sw")
        ToolTip(save_qtraces_button, "Saves this list to a .csv file.")
        
        load_qtraces_button = ttk.Button(quick_trace_window, text="Load Traces", style="small_button_qtw_style1.TButton", command=load_qtraces_csv)
        load_qtraces_button.grid(row=2, column=0, padx=(10, 10), pady=(0,10), sticky="se")
        ToolTip(load_qtraces_button, "Loads traces from a .csv file. It must contain a header (color,mz,tol), and each following line will contain the color hexadecimal code (if none is provided, a random color is picked), mz value and tolerance (in mz value), like this:\n        color,mz,tolerance\n        #FFA500,785.2,0.5\n        #808000,800.5,0.5")
        
        quick_trace_list.bind("<ButtonRelease-1>", click_quick_trace_list)
        quick_trace_list.bind("<KeyRelease-Up>", click_quick_trace_list)
        quick_trace_list.bind("<KeyRelease-Down>", click_quick_trace_list)
        quick_trace_list.bind("<Motion>", on_quick_trace_list_motion)
        
        quick_trace_window.update_idletasks()
        quick_trace_window.deiconify()
        window_width = quick_trace_window.winfo_width()
        window_height = quick_trace_window.winfo_height()
        screen_width = quick_trace_window.winfo_screenwidth()
        screen_height = quick_trace_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        quick_trace_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
    def on_focus_in_filter_list(event):
        if filter_list.get() == "Filter the list of glycans...":
            filter_list.configure(fg='black')
            filter_list.delete(0, tk.END)
            
    def on_focus_out_filter_list(event):
        if filter_list.get() == '':
            filter_list.configure(fg='grey')
            filter_list.insert(0, "Filter the list of glycans...")
            
    def on_scale_selected(event):
        global plot_scale
        selected_scale = scaling_dropdown.get()
        if selected_scale == 'Linear':
            ax_spec.set_yscale('linear')
            ax_spec_ms2.set_yscale('linear')
        if selected_scale == 'Log':
            ax_spec.set_yscale('symlog')
            ax_spec_ms2.set_yscale('symlog')
        if selected_scale == 'Sqrt':
            ax_spec.set_yscale('function', functions=(custom_sqrt, custom_square))
            ax_spec_ms2.set_yscale('function', functions=(custom_sqrt, custom_square))
        
        ax_spec.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax_spec.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        ax_spec_ms2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax_spec_ms2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        canvas_spec.draw_idle()
        canvas_spec_ms2.draw_idle()

    def toggle_spectrum_ruler():
        """Toggle the ruler tool on or off and ensure exclusivity with pan/zoom."""
        global spectrum_ruler_on
        
        spectrum_ruler_on = not spectrum_ruler_on
        if spectrum_ruler_on:
            spectrum_ruler_frame.config(bg="blue")
        else:
            spectrum_ruler_frame.config(bg=background_color)
            
    def remove_spectrum_ruler():
        global current_line_ruler, distance_ruler_text, ruler_lines_list, selected_ruler
        
        if selected_ruler != None:
            for i in selected_ruler[0]:
                i[0].remove()
            del ruler_lines_list[selected_ruler[1]]
            selected_ruler = None
        else:
            for i in ruler_lines_list:
                for j in i:
                    j[0].remove()
            ruler_lines_list = []
        
        canvas_spec.draw_idle()
        canvas_spec_ms2.draw_idle()
    
    def check_mzml_window(data):
        '''
        '''     
        check_folder = os.path.join(tempfile.gettempdir())
        
        this_process_id = os.getpid()
        
        if f"mzml_window_{this_process_id}.txt" not in os.listdir(check_folder):
            with open(os.path.join(check_folder, f"mzml_window_{this_process_id}.txt"), 'w') as f:
                f.write("mzML editing window is open")
                f.close()
                
            handle_selection_vigilance()
            
            mzml_window_start(data)
                
        else:
            mzml_window_start(to_lift = True)
            
    def toggle_gg_draw(treeview):
        '''
        '''
        global gg_draw_on, gg_draw_list, photo_ggdraw, photo_ggdraw_off, gg_draw_zoom_scale, gg_draw_zoom_scale_value, gg_draw_zoom_scale_value_float
        if gg_draw_on:
            gg_draw_on = False
            gg_draw.config(image=photo_ggdraw)
            gg_draw_zoom_scale_value_float = gg_draw_zoom_scale_value.get()
            
            gg_draw_zoom_scale.destroy()
            
            for glycan in gg_draw_list:
                drawing = glycan[1]
                drawing.remove()
            gg_draw_list = []
            canvas.draw_idle()
        else:
            gg_draw_on = True
            gg_draw.config(image=photo_ggdraw_off)
            
            gg_draw_zoom_scale_value = tk.DoubleVar()
            gg_draw_zoom_scale = ttk.Scale(main_window, from_=0.1, to=0.5, orient='horizontal', variable=gg_draw_zoom_scale_value, command=gg_draw_zoom_scale_moved)
            gg_draw_zoom_scale.grid(row=0, column=2, padx=(10,200), sticky='se')
            gg_draw_zoom_scale_value.set(gg_draw_zoom_scale_value_float)
            ToolTip(gg_draw_zoom_scale, "Set the size of the glycans drawings on the chromatogram/electropherogram plot.")
            
            if reanalysis_path != "":
                check_gg_drawings_available(treeview)
    
    def gg_draw_zoom_scale_moved(value):
        '''
        '''
        global gg_draw_list
        for index, glycan in enumerate(gg_draw_list):
            recreate_gg_draw(index, glycan[-1], zoom_level = float(value))
    
    global check_gg_drawings_available
    def check_gg_drawings_available(treeview, draw = True):
        '''
        '''
        directory_content = os.listdir(gg_draw_glycans_path)
        glycans_to_draw = []
        selected_glycans = []
        for element in treeview.selection():
        
            # Get the level of the treeview
            glycan = element
            while True:
                test_parent = treeview.parent(glycan)
                if test_parent != '':
                    glycan = treeview.parent(glycan)
                else:
                    break
            selected_glycans.append(treeview.item(glycan, "text"))
            
        for glycan in selected_glycans:
            comp = General_Functions.form_to_comp_glycans(glycan)
            possibilities = [glycan]
            
            if 'S' in comp.keys() or 'G' in comp.keys():
                if 'S' in comp.keys():
                    new_combos = der_sialics_combo('S', comp['S'])
                    for i in new_combos:
                        new_comp = copy.deepcopy(comp)
                        new_comp = replace_dict_keys(new_comp, 'S', i)
                        possibilities.append(General_Functions.comp_to_formula(new_comp))
                        
                if 'G' in comp.keys():
                    new_combos = der_sialics_combo('G', comp['G'])
                    for i in new_combos:
                        new_comp = copy.deepcopy(comp)
                        new_comp = replace_dict_keys(new_comp, 'G', i)
                        possibilities.append(General_Functions.comp_to_formula(new_comp))
                        
            for possibility in possibilities:
                if possibility not in directory_content:
                    glycans_to_draw.append(possibility)
                    
        if len(glycans_to_draw) > 0:
            start_glycans_drawing(glycans_to_draw, treeview)
            
        if len(glycans_to_draw) == 0 and draw:
            draw_glycans_on_chromatogram(treeview)
                
    def start_glycans_drawing(glycans_list, treeview):
        '''
        '''
        def close_gd():
            background_glycan_drawing.destroy()
            
        def draw_glycans_list(glycans_list, treeview):
            cpu_number = (os.cpu_count())-2 if os.cpu_count() < 60 else 60
            if cpu_number <= 0:
                cpu_number = 1
            results = []
            with concurrent.futures.ProcessPoolExecutor(max_workers = cpu_number) as executor:
                for glycan in glycans_list:
                    result = executor.submit(draw_glycan, 
                                             gg_draw_glycans_path,
                                             30,
                                             glycan,
                                             'none')
                    results.append(result)
                        
            for index, i in enumerate(results):
                result_data = i.result()
                results[index] = None
                result_data = None
                
            close_gd()
            draw_glycans_on_chromatogram(treeview)
            
        background_glycan_drawing = tk.Toplevel()
        # background_glycan_drawing.attributes("-topmost", True)
        background_glycan_drawing.withdraw()
        background_glycan_drawing.title("Drawing Glycans")
        icon = ImageTk.PhotoImage(ico_image)
        background_glycan_drawing.iconphoto(False, icon)
        background_glycan_drawing.resizable(False, False)
        background_glycan_drawing.grab_set()
        background_glycan_drawing.protocol("WM_DELETE_WINDOW", on_closing)
        
        generating_glycans_label = ttk.Label(background_glycan_drawing, text=f"Generating glycans drawings...\nThis may take up to a few minutes.\nPlease wait.", font=("Segoe UI", list_font_size))
        generating_glycans_label.pack(pady=20, padx=50)
        
        background_glycan_drawing.update_idletasks()
        background_glycan_drawing.deiconify()
        window_width = background_glycan_drawing.winfo_width()
        window_height = background_glycan_drawing.winfo_height()
        screen_width = background_glycan_drawing.winfo_screenwidth()
        screen_height = background_glycan_drawing.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        background_glycan_drawing.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
        t = threading.Thread(target=draw_glycans_list, args=(glycans_list, treeview))
        t.start()
        
    def draw_glycans_on_chromatogram(treeview):
        '''
        '''
        global gg_draw_list, forced, gg_draw_zoom_scale_value
            
        for glycan in gg_draw_list:
            drawing = glycan[1]
            try:
                drawing.remove()
            except:
                pass
                
        gg_draw_list = []
        
        to_plot = []
        
        for element in ax.get_lines():
            label = element.get_label()
            
            if label[0] == "_":
                continue
            
            if len(label.split("+")) > 1:
                glycan = label.split("+")[0]
                adduct = label.split("+")[-1]
                level = 2
            elif len(label.split(" - ")) > 1:
                glycan = label.split(" - ")[0]
                adduct = label.split(" - ")[-1]
                level = 3
            else:
                glycan = label
                adduct = None
                level = 1
            
            rts = []
            
            if level == 1:
                for adduct in glycans_per_sample[samples_dropdown.get()][glycan]:
                    for index, peak in enumerate(glycans_per_sample[samples_dropdown.get()][glycan][adduct]['peaks']):
                        
                        if hide_bad_peaks_checkbox_state.get():
                            peak_quality_fails = check_peak_quality(ppm = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['ppm'][index], iso = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['iso'][index], curve = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['curve'][index], sn = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['sn'][index])
                            if peak_quality_fails == 0 and peak not in rts:
                                rts.append(peak)
                        else:
                            if peak not in rts:
                                rts.append(peak)
                                
            elif level == 2 or level == 3:
                for index, peak in enumerate(glycans_per_sample[samples_dropdown.get()][glycan][adduct]['peaks']):
                    if peak not in rts:
                        rts.append(peak)
            
            if len(rts) == 0:
                continue
            
            intensities = []
            y_check_data = element.get_ydata()
            x_check_data = element.get_xdata()
            max_int = np.max(y_check_data)
            for rt in rts:
                x_coord = np.abs(x_check_data - rt).argmin()
                intensities.append(y_check_data[x_coord])
            
            # Find the possibilities of drawings for the given glycan
            possibilities = []
            if glycan in os.listdir(gg_draw_glycans_path):
                directory_files = os.listdir(os.path.join(gg_draw_glycans_path, f"{glycan}"))
                for file in directory_files:
                    if parameters_gg['library settings']['glycan class'] == 'n_glycans':
                        if file.split(".")[0].split("_")[2] == "N":
                            possibilities.append(file)
                    elif parameters_gg['library settings']['glycan class'] == 'o_glycans':
                        if file.split(".")[0].split("_")[2] == "O":
                            possibilities.append(file)
                    else:
                        possibilities.append(file)
            
            if level == 1:
                sorted_pairs = sorted(zip(rts, intensities))
                rts, intensities = zip(*sorted_pairs)
                
                last_rt = 0
                new_rts = []
                new_ints = []
                for index, rt in enumerate(rts):
                    if abs(rt-last_rt) > 0.2:
                        new_rts.append(rt)
                        new_ints.append(intensities[index])
                    else:
                        new_ints[-1] = max([new_ints[-1], intensities[index]])
                    last_rt = rt
                
                rts = new_rts
                intensities = new_ints
            
            for x, y in zip(rts, intensities):
                if len(possibilities) > 0:
                    to_plot.append([glycan, x, y, f"{glycan}/{random.choice(possibilities)}"])
                
        to_plot = sorted(to_plot, key=lambda item: item[1])
        for index, glycan in enumerate(to_plot):
            x_offset = -1 if index-math.floor(len(to_plot)/2) < 0 else 1
            x = glycan[1]
            y = glycan[2]
            max_int = ax.get_ylim()[1]
            xbox = x+(x_offset*((ax.get_xlim()[1]-ax.get_xlim()[0])*0.1))
            if y < max_int*0.2:
                ybox = y+(max_int*0.4)
            elif y > max_int*0.8:
                ybox = y-(max_int*0.2)
            else:
                ybox = y+(max_int*0.2)
            
            image = mpimg.imread(os.path.join(gg_draw_glycans_path, f"{glycan[3]}"))
            imagebox = OffsetImage(image, zoom=gg_draw_zoom_scale_value.get())
            ab = AnnotationBbox(imagebox, xy=[x, y], xybox=[xbox, ybox], frameon=True, bboxprops=dict(boxstyle='square', fc='none', ec='none', lw=1, linestyle=':'), xycoords='data', arrowprops=dict(arrowstyle="->", color='black', connectionstyle="angle3,angleA=90,angleB=0"))
            ax.add_artist(ab)
            
            gg_draw_list.append([glycan[0], ab, [x,y], [xbox, ybox], glycan[3]])
            
        canvas.draw_idle()
        
    def recreate_gg_draw(gg_draw_list_index, new_image, zoom_level = '', mirrored = False):
        '''
        '''
        global gg_draw_zoom_scale_value
        if zoom_level == '':
            zoom_level = gg_draw_zoom_scale_value.get()
        gg_draw_list[gg_draw_list_index][1].remove()
        image = mpimg.imread(os.path.join(gg_draw_glycans_path, f"{new_image}"))
        if mirrored:
            image = np.flip(image, axis=1)
        imagebox = OffsetImage(image, zoom=zoom_level)
        ab = AnnotationBbox(imagebox, xy=gg_draw_list[gg_draw_list_index][2], xybox=gg_draw_list[gg_draw_list_index][3], frameon=True, bboxprops=dict(boxstyle='square', fc='none', ec='none', lw=1, linestyle=':'), xycoords='data', arrowprops=dict(arrowstyle="->", color='black', connectionstyle="angle3,angleA=90,angleB=0"))
        ax.add_artist(ab)
        
        gg_draw_list[gg_draw_list_index] = [gg_draw_list[gg_draw_list_index][0], ab, gg_draw_list[gg_draw_list_index][2], gg_draw_list[gg_draw_list_index][3], new_image]
        
        canvas.draw_idle()
        
    def glycan_browser(gg_draw_list_index):
        '''
        '''
        def on_closing_glycan_browser():
            global gg_draw_selected, panning_enabled
            browse_glycans_window.destroy()
            gg_draw_selected = None
            panning_enabled = False
            
        def bissecting_glcnac_checkbox_variable_check():
            state = bissecting_glcnac_checkbox_variable.get()
            update_pictures(None)
        
        def core_fuc_checkbox_variable_check():
            state = core_fuc_checkbox_variable.get()
            update_pictures(None)
        
        def mirror_checkbox_variable_check():
            state = mirror_checkbox_variable.get()
            
        def change_glycan():
            if gallery.get_selected():
                recreate_gg_draw(gg_draw_list_index, "/".join(gallery.get_selected().split("/")[-2:]).split("\\")[-1], mirrored=mirror_checkbox_variable.get())
            
        def remove_glycan():
            gg_draw_list[gg_draw_list_index][1].remove()
            del gg_draw_list[gg_draw_list_index]
            canvas.draw_idle()
            on_closing_glycan_browser()
            
        def update_pictures(event):
            global gallery_warning_label
            
            glycan_class = glycan_class_combobox.get()
            
            if glycan_class == 'All':
                glycan_class = ''
            elif glycan_class == 'N-Glycans':
                glycan_class = 'N'
            elif glycan_class == 'O-Glycans':
                glycan_class = 'O'
            
            antennas_no = antennas_no_combobox.get()
            if antennas_no == 'All':
                antennas_no = ''
            elif antennas_no == 'Generic':
                antennas_no = 'G'
                
            core_fuc = core_fuc_checkbox_variable.get()
            bissecting_glcnac = bissecting_glcnac_checkbox_variable.get()
            
            new_possibilities = []
            
            for picture in possibilities:
                picture_parts = picture.split("/")[-1].split(".")[0].split("_")
                
                if glycan_class in picture_parts[2]:
                    if antennas_no in picture_parts[1]:
                        if core_fuc and "F" in picture_parts:
                            if bissecting_glcnac and "B" in picture_parts:
                                new_possibilities.append(picture)
                            elif not bissecting_glcnac:
                                new_possibilities.append(picture)
                        elif not core_fuc:
                            if bissecting_glcnac and "B" in picture_parts:
                                new_possibilities.append(picture)
                            elif not bissecting_glcnac:
                                new_possibilities.append(picture)
                                
            if len(new_possibilities) > 750:
                gallery.update_images(new_possibilities[:750])
                try:
                    gallery_warning_label.destroy()
                except:
                    pass
                gallery.grid_configure(pady=(30, 30))
                gallery_warning_label = ttk.Label(browse_glycans_window, text='Displaying 750 options. Use the filters to explore all structures.', font=("Segoe UI", list_font_size))
                gallery_warning_label.grid(row=3, column=0, padx=(10, 10), pady=(5, 5), sticky='sew')
            else:
                gallery.update_images(new_possibilities)
                try:
                    gallery_warning_label.destroy()
                except:
                    pass
                gallery.grid_configure(pady=(30, 10))
        
        comp = General_Functions.form_to_comp_glycans(gg_draw_list[gg_draw_list_index][0])
        possibilities_glycans = [gg_draw_list[gg_draw_list_index][0]]
        
        if 'S' in comp.keys() or 'G' in comp.keys():
            if 'S' in comp.keys():
                new_combos = der_sialics_combo('S', comp['S'])
                for i in new_combos:
                    new_comp = copy.deepcopy(comp)
                    new_comp = replace_dict_keys(new_comp, 'S', i)
                    possibilities_glycans.append(General_Functions.comp_to_formula(new_comp))
                    
            if 'G' in comp.keys():
                new_combos = der_sialics_combo('G', comp['G'])
                for i in new_combos:
                    new_comp = copy.deepcopy(comp)
                    new_comp = replace_dict_keys(new_comp, 'G', i)
                    possibilities_glycans.append(General_Functions.comp_to_formula(new_comp))
                # to_add = []
                # for new_glycan_wip in possibilities_glycans[1:-1]:
                    # comp = General_Functions.form_to_comp(new_glycan_wip)
                    # new_combos = der_sialics_combo('G', comp['G'])
                    # for i in new_combos:
                        # new_comp = copy.deepcopy(comp)
                        # new_comp = replace_dict_keys(new_comp, 'G', i)
                        # to_add.append(General_Functions.comp_to_formula(new_comp))
                # for item in to_add:
                    # possibilities_glycans.append(to_add)
        
        possibilities = []
        generic_n_added = False
        generic_o_added = False
        directory_files = os.listdir(gg_draw_glycans_path)
        for file in directory_files:
            if file in possibilities_glycans:
                for glycan in os.listdir(os.path.join(gg_draw_glycans_path, f"{file}")):
                    if "G" in glycan.split("_") and "N" in glycan.split(".")[0].split("_"):
                        if generic_n_added:
                            continue
                        else:
                            generic_n_added = True
                    if "G" in glycan.split("_") and "O" in glycan.split(".")[0].split("_"):
                        if generic_o_added:
                            continue
                        else:
                            generic_o_added = True
                    possibilities.append(os.path.join(gg_draw_glycans_path, f"{file}/{glycan}"))
        
        browse_glycans_window = tk.Toplevel()
        icon = ImageTk.PhotoImage(ico_image)
        browse_glycans_window.iconphoto(False, icon)
        browse_glycans_window.withdraw()
        browse_glycans_window.grab_set()
        browse_glycans_window.bind("<Configure>", on_resize)
        browse_glycans_window.title(f"Glycan Browser - {gg_draw_list[gg_draw_list_index][0]}")
        browse_glycans_window.resizable(False, False)
        browse_glycans_window.geometry("400x480")
        browse_glycans_window.grid_rowconfigure(3, weight=1)
        browse_glycans_window.protocol("WM_DELETE_WINDOW", on_closing_glycan_browser)
        
        # Glycan class
        glycan_class_options = ['All', 'N-Glycans', 'O-Glycans']
        glycan_class_label = ttk.Label(browse_glycans_window, text='Glycan Class:', font=("Segoe UI", list_font_size))
        glycan_class_label.grid(row=0, column=0, padx=(10, 10), pady=(0, 10), sticky='nw')
        glycan_class_combobox = ttk.Combobox(browse_glycans_window, state="readonly", values=glycan_class_options, width=10)
        glycan_class_combobox.grid(row=0, column=0, padx=(95, 10), pady=(0, 10), sticky='nw')
        glycan_class_combobox.set('All')
        glycan_class_combobox.bind("<<ComboboxSelected>>", update_pictures)
        ToolTip(glycan_class_label, "Filter the glycans by their class (e.g. N-glycans, O-Glycans, etc.)")
        ToolTip(glycan_class_combobox, "Filter the glycans by their class (e.g. N-glycans, O-Glycans, etc.)")
        
        # Antennas No.
        antennas_no_options = ['All', 'Generic', '1', '2', '3', '4']
        antennas_no_label = ttk.Label(browse_glycans_window, text='Number of Antennas:', font=("Segoe UI", list_font_size))
        antennas_no_label.grid(row=0, column=0, padx=(188, 10), pady=(0, 10), sticky='nw')
        antennas_no_combobox = ttk.Combobox(browse_glycans_window, state="readonly", values=antennas_no_options, width=7)
        antennas_no_combobox.grid(row=0, column=0, padx=(320, 10), pady=(0, 10), sticky='nw')
        antennas_no_combobox.set('All')
        antennas_no_combobox.bind("<<ComboboxSelected>>", update_pictures)
        ToolTip(antennas_no_label, "Filter the glycans by the amount of antennas they have.")
        ToolTip(antennas_no_combobox, "Filter the glycans by the amount of antennas they have.")
        
        # Bissecting GlcNAc
        bissecting_glcnac_checkbox_variable = tk.BooleanVar(value=False)
        bissecting_glcnac_checkbox = ttk.Checkbutton(browse_glycans_window, text="Only with Bissecting GlcNAc", variable=bissecting_glcnac_checkbox_variable, command=bissecting_glcnac_checkbox_variable_check)
        bissecting_glcnac_checkbox.grid(row=1, column=0, padx=(10, 10), pady=(0, 10), sticky='nw')
        ToolTip(bissecting_glcnac_checkbox, "Display only glycans with bissecting GlcNAc.")
        
        # Change Glycan button
        change_glycan_button = ttk.Button(browse_glycans_window, text="    Change Glycan   ", style="small_button_sfw_style1.TButton", command=change_glycan)
        change_glycan_button.grid(row=2, column=0, padx=(10, 20), pady=(0, 0), sticky='ne')
        
        # Core Fucosylation
        core_fuc_checkbox_variable = tk.BooleanVar(value=False)
        core_fuc_checkbox = ttk.Checkbutton(browse_glycans_window, text="Only with Core Fucosylation", variable=core_fuc_checkbox_variable, command=core_fuc_checkbox_variable_check)
        core_fuc_checkbox.grid(row=2, column=0, padx=(10, 10), pady=(0, 10), sticky='nw')
        ToolTip(core_fuc_checkbox, "Display only glycans with core fucosylation.")
        
        # Mirrored
        mirror_checkbox_variable = tk.BooleanVar(value=False)
        mirror_checkbox = ttk.Checkbutton(browse_glycans_window, text="Mirror the glycan", variable=mirror_checkbox_variable, command=mirror_checkbox_variable_check)
        mirror_checkbox.grid(row=3, column=0, padx=(10, 10), pady=(0, 10), sticky='nw')
        ToolTip(mirror_checkbox, "Mirrors the glycan on the y axis (e.g. left elements will be on the right side and vice-versa.)")
        
        # Remove Glycan button
        remove_glycan_button = ttk.Button(browse_glycans_window, text="   Remove Glycan   ", style="small_button_sfw_style1.TButton", command=remove_glycan)
        remove_glycan_button.grid(row=3, column=0, padx=(10, 20), pady=(0, 10), sticky='ne')
        
        if len(possibilities) > 750:
            gallery = ImageGallery(browse_glycans_window, possibilities[:750], columns=3, cell_size=(125, 125))
            gallery.grid(row=3, column=0, padx=(0, 10), pady=(30, 30), sticky="nsew")
            
            global gallery_warning_label
            gallery_warning_label = ttk.Label(browse_glycans_window, text='Displaying 750 options. Use the filters to explore all structures.', font=("Segoe UI", list_font_size))
            gallery_warning_label.grid(row=3, column=0, padx=(10, 10), pady=(5, 5), sticky='sew')
        else:
            gallery = ImageGallery(browse_glycans_window, possibilities, columns=3, cell_size=(125, 125))
            gallery.grid(row=3, column=0, padx=(0, 10), pady=(30, 10), sticky="nsew")
        
        browse_glycans_window.update_idletasks()
        browse_glycans_window.deiconify()
        window_width = browse_glycans_window.winfo_width()
        window_height = browse_glycans_window.winfo_height()
        screen_width = browse_glycans_window.winfo_screenwidth()
        screen_height = browse_glycans_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        browse_glycans_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        global gathering_glycan_options    
        gathering_glycan_options.destroy()
                
    def gg_draw_options(gg_draw_list_index):
        '''
        '''
        
        global gathering_glycan_options
        gathering_glycan_options = tk.Toplevel()
        # gathering_glycan_options.attributes("-topmost", True)
        gathering_glycan_options.withdraw()
        gathering_glycan_options.title("Gathering Glycans Options")
        icon = ImageTk.PhotoImage(ico_image)
        gathering_glycan_options.iconphoto(False, icon)
        gathering_glycan_options.resizable(False, False)
        gathering_glycan_options.grab_set()
        gathering_glycan_options.protocol("WM_DELETE_WINDOW", on_closing)
        
        gathering_glycan_options_label = ttk.Label(gathering_glycan_options, text=f"Gathering the structure options for this glycan.\nPlease wait.", font=("Segoe UI", list_font_size))
        gathering_glycan_options_label.pack(pady=20, padx=50)
        
        gathering_glycan_options.update_idletasks()
        gathering_glycan_options.deiconify()
        window_width = gathering_glycan_options.winfo_width()
        window_height = gathering_glycan_options.winfo_height()
        screen_width = gathering_glycan_options.winfo_screenwidth()
        screen_height = gathering_glycan_options.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        gathering_glycan_options.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
        t = threading.Thread(target=glycan_browser, args=(gg_draw_list_index,))
        t.start()
        
    def loading_files_from_explorer():
        def close_lf():
            loading_files.destroy()
            
        def wait_thread():
            if len(reanalysis_path) > 0:
                reanalysis_file.join()
            if len(samples_dropdown['values']) != 0:
                samples_dropdown.current(0)
                handle_selection(None)
            close_lf()
        
        global loading_files
        loading_files = tk.Toplevel()
        # loading_files.attributes("-topmost", True)
        loading_files.withdraw()
        loading_files.title("Loading Files")
        icon = ImageTk.PhotoImage(ico_image)
        loading_files.iconphoto(False, icon)
        loading_files.resizable(False, False)
        loading_files.grab_set()
        loading_files.protocol("WM_DELETE_WINDOW", on_closing)
        
        loading_files_label = ttk.Label(loading_files, text="Loading files, please wait.", font=("Segoe UI", list_font_size))
        loading_files_label.pack(pady=35, padx=70)
        
        loading_files.update_idletasks()
        loading_files.deiconify()
        window_width = loading_files.winfo_width()
        window_height = loading_files.winfo_height()
        screen_width = loading_files.winfo_screenwidth()
        screen_height = loading_files.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        loading_files.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        t = threading.Thread(target=wait_thread)
        t.start()
        
    def create_zoomlock_window(plot_type):
        def on_type(event, plot_type):
            lims = {'x_min' : x_label_min_entry.get(), 'x_max' : x_label_max_entry.get(), 'y_min' : y_label_min_entry.get(), 'y_max' : y_label_max_entry.get()}
            
            lims_save_status = [True, True, True, True]
            
            for axis_limit, value in lims.items():
                try:
                    lims[axis_limit] = float(value)
                except Exception:
                    pass
                
                if axis_limit == 'x_min' and (type(lims[axis_limit]) == str or lims[axis_limit] == 0):
                    lims[axis_limit] = original_xlims[0]
                    lims_save_status[0] = False
                if axis_limit == 'x_max' and (type(lims[axis_limit]) == str or lims[axis_limit] == 0):
                    lims[axis_limit] = original_xlims[1]
                    lims_save_status[1] = False
                if axis_limit == 'y_min' and (type(lims[axis_limit]) == str or lims[axis_limit] == 0):
                    lims[axis_limit] = original_ylims[0]
                    lims_save_status[2] = False
                if axis_limit == 'y_max' and (type(lims[axis_limit]) == str or lims[axis_limit] == 0):
                    lims[axis_limit] = original_ylims[1]
                    lims_save_status[3] = False
            
            global zoom_lock_chromatogram, zoom_lock_spectrum
            if plot_type == 'chromatogram':
                ax.set_xlim(lims['x_min'], lims['x_max'])
                ax.set_ylim(lims['y_min'], lims['y_max'])
                canvas.draw_idle()
                
                zoom_lock_chromatogram = {'on': True, 'x_axis': (lims['x_min'] if lims_save_status[0] else 0, lims['x_max'] if lims_save_status[1] else 0), 'y_axis': (lims['y_min'] if lims_save_status[2] else 0, lims['y_max'] if lims_save_status[3] else 0)}
            else:
                ax_spec.set_xlim(lims['x_min'], lims['x_max'])
                ax_spec.set_ylim(lims['y_min'], lims['y_max'])
                canvas_spec.draw_idle()
                
                ax_spec_ms2.set_xlim(lims['x_min'], lims['x_max'])
                ax_spec_ms2.set_ylim(lims['y_min'], lims['y_max'])
                canvas_spec_ms2.draw_idle()
                
                zoom_lock_spectrum = {'on': True, 'x_axis': (lims['x_min'] if lims_save_status[0] else 0, lims['x_max'] if lims_save_status[1] else 0), 'y_axis': (lims['y_min'] if lims_save_status[2] else 0, lims['y_max'] if lims_save_status[3] else 0)}
            
        def reset_axis(plot_type):
            global zoom_lock_chromatogram, zoom_lock_spectrum
            
            x_label_min_entry.delete(0, tk.END)
            x_label_min_entry.insert(0, "0")
            
            x_label_max_entry.delete(0, tk.END)
            x_label_max_entry.insert(0, "0")
            
            y_label_min_entry.delete(0, tk.END)
            y_label_min_entry.insert(0, "0")
            
            y_label_max_entry.delete(0, tk.END)
            y_label_max_entry.insert(0, "0")
            
            if plot_type == 'chromatogram':
                zoom_lock_chromatogram = {'on': False, 'x_axis': (0, 0), 'y_axis': (0, 0)}
                try:
                    ax.set_xlim(og_x_range[0], og_x_range[1])
                    ax.set_ylim(og_y_range[0], og_y_range[1])
                    canvas.draw_idle()
                except Exception:
                    pass
            else:
                zoom_lock_spectrum = {'on': False, 'x_axis': (0, 0), 'y_axis': (0, 0)}
                try:
                    ax_spec.set_xlim(og_x_range_spec[0], og_x_range_spec[1])
                    ax_spec.set_ylim(og_y_range_spec[0], og_y_range_spec[1])
                    canvas_spec.draw_idle()
                    
                    ax_spec_ms2.set_xlim(og_x_range_spec_ms2[0], og_x_range_spec_ms2[1])
                    ax_spec_ms2.set_ylim(og_y_range_spec_ms2[0], og_y_range_spec_ms2[1])
                    canvas_spec_ms2.draw_idle()
                except Exception:
                    pass
            
        def on_focus_loss(event):
            if not any(widget.winfo_exists() and widget.focus_displayof() for widget in zoomlock_window.winfo_children()):
                zoomlock_window.destroy()
            
        x = main_window.winfo_pointerx()+15
        y = main_window.winfo_pointery()+15
        screen_width = main_window.winfo_screenwidth()
        screen_height = main_window.winfo_screenheight()
        
        widget_width = 150
        widget_height = 80
        
        if x + widget_width > screen_width:  # Tooltip goes off the screen
            x = main_window.winfo_pointerx()-15 - widget_width  # Move to the left side
        if y + widget_height > screen_height:  # Tooltip goes off the screen
            y = main_window.winfo_pointery()-15 - widget_height  # Move to the left side
        
        if plot_type == 'chromatogram':
            original_xlims = ax.get_xlim()
            original_ylims = ax.get_ylim()
        else:
            original_xlims = ax_spec.get_xlim()
            original_ylims = ax_spec.get_ylim()
        
        zoomlock_window = tk.Toplevel()
        zoomlock_window.overrideredirect(True)
        zoomlock_window.geometry(f"+{x}+{y}")
        zoomlock_window.attributes('-topmost', True)
        zoomlock_window.bind("<FocusOut>", on_focus_loss)
        
        zoomlock_window.focus_set()

        zoom_lock_window_frame = tk.Frame(zoomlock_window, highlightbackground="grey", highlightthickness=1)
        zoom_lock_window_frame.pack(fill=tk.BOTH, expand=True)
        
        if plot_type == 'chromatogram':
            x_label = 'RT/MT'
        else:
            x_label = 'm/z'
        y_label = 'Int'
        
        # x-label widgets
        x_label_widget = ttk.Label(zoom_lock_window_frame, text=x_label, font=("Segoe UI", list_font_size_smaller))
        x_label_widget.grid(row=0, column=0, padx=(5, 0), pady=(5, 0))
        
        x_label_min_entry = ttk.Entry(zoom_lock_window_frame, width = 6)
        x_label_min_entry.grid(row=0, column=1, pady=(5, 0))
        x_label_min_entry.insert(0, f"{zoom_lock_chromatogram['x_axis'][0] if plot_type == 'chromatogram' else zoom_lock_spectrum['x_axis'][0]}")
        x_label_min_entry.bind("<KeyRelease>", lambda event: on_type(event, plot_type))
        
        x_hyphen_label = ttk.Label(zoom_lock_window_frame, text="-", font=("Segoe UI", list_font_size_smaller))
        x_hyphen_label.grid(row=0, column=2, pady=(5, 0))
        
        x_label_max_entry = ttk.Entry(zoom_lock_window_frame, width = 6)
        x_label_max_entry.grid(row=0, column=3, padx=(0, 5), pady=(5, 0))
        x_label_max_entry.insert(0, f"{zoom_lock_chromatogram['x_axis'][1] if plot_type == 'chromatogram' else zoom_lock_spectrum['x_axis'][1]}")
        x_label_max_entry.bind("<KeyRelease>", lambda event: on_type(event, plot_type))
        
        # y-label widgets
        y_label_widget = ttk.Label(zoom_lock_window_frame, text=y_label, font=("Segoe UI", list_font_size_smaller))
        y_label_widget.grid(row=1, column=0, pady=10, padx=(5, 0))
        
        y_label_min_entry = ttk.Entry(zoom_lock_window_frame, width = 6)
        y_label_min_entry.grid(row=1, column=1, pady=10)
        y_label_min_entry.insert(0, f"{zoom_lock_chromatogram['y_axis'][0] if plot_type == 'chromatogram' else zoom_lock_spectrum['y_axis'][0]}")
        y_label_min_entry.bind("<KeyRelease>", lambda event: on_type(event, plot_type))
        
        x_hyphen_label = ttk.Label(zoom_lock_window_frame, text="-", font=("Segoe UI", list_font_size_smaller))
        x_hyphen_label.grid(row=1, column=2, pady=10)
        
        y_label_max_entry = ttk.Entry(zoom_lock_window_frame, width = 6)
        y_label_max_entry.grid(row=1, column=3, pady=10, padx=(0, 5))
        y_label_max_entry.insert(0, f"{zoom_lock_chromatogram['y_axis'][1] if plot_type == 'chromatogram' else zoom_lock_spectrum['y_axis'][1]}")
        y_label_max_entry.bind("<KeyRelease>", lambda event: on_type(event, plot_type))
        
        # reset button
        reset_lock_button = ttk.Button(zoom_lock_window_frame, text="Reset", command=lambda:reset_axis(plot_type))
        reset_lock_button.grid(row=2, column=0, columnspan=4, pady=(0, 5))
        
    def toggle_hide_bad_peaks():
        global original_lines
        state = hide_bad_peaks_checkbox_state.get()
        
        graph_lines = [line for line in ax.get_lines() if line.get_label()[0] != "_"]
        
        if len(graph_lines) > 1:
            fill_patches = [obj for obj in ax.collections]
            for patch in fill_patches:
                patch.remove()
        
        if state:
            # Access the plotted lines
            for element in graph_lines:
                label = element.get_label()
                
                if label[0] == "_" or label == 'Base Peak Chromatogram/Electropherogram' or '±' in label:
                    continue
                
                if len(label.split("+")) > 1:
                    label_split = label.split("+")
                    glycan = label_split[0]
                    
                    for index, part in enumerate(label_split):
                        if "(s)" in part:
                            glycan+= "+"+label_split[index]
                        if "(p)" in part:
                            glycan+= "+"+label_split[index]
                        if "1Na" in part or "2Na" in part:
                            glycan+= "+"+label_split[index]
                            
                    if len(label.split(" - ")) > 1:
                        adduct = glycan.split(" - ")[-1]
                        glycan = glycan.split(" - ")[0]
                        level = 3
                    elif label_split[-1] in glycan:
                        adduct = None
                        level = 1
                    else:
                        adduct = label_split[-1]
                        level = 2
                elif len(label.split(" - ")) > 1:
                    glycan = label.split(" - ")[0]
                    adduct = label.split(" - ")[-1]
                    level = 3
                else:
                    glycan = label
                    adduct = None
                    level = 1
                
                bad_peaks = []
                
                if level == 1:
                    good_rts = []
                    for adduct in glycans_per_sample[samples_dropdown.get()][glycan]:
                        for index, peak in enumerate(glycans_per_sample[samples_dropdown.get()][glycan][adduct]['peaks']):
                            # To filter by quality
                            current_ppm = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['ppm'][index]
                            current_iso = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['iso'][index]
                            current_curve = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['curve'][index]
                            snratio = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['sn'][index]
                            
                            if current_ppm > max_ppm[1] or current_ppm < max_ppm[0] or current_iso < iso_fit_score or current_curve < curve_fit_score or snratio < s_to_n:
                                found = False
                                for peak_good in good_rts:
                                    peak_rt = float(peak_good.split("_")[-1])
                                    if abs(peak_rt - peak) < 0.2:
                                        found = True
                                        break
                                if not found:    
                                    bad_peaks.append(f"{glycan}_{adduct}_{peak}")
                            else:
                                for index_bad, peak_bad in enumerate(bad_peaks):
                                    peak_rt = float(peak_bad.split("_")[-1])
                                    if abs(peak_rt - peak) < 0.2:
                                        del bad_peaks[index_bad]
                                        break
                                good_rts.append(f"{glycan}_{adduct}_{peak}")
                                
                elif level == 2:
                    for index, peak in enumerate(glycans_per_sample[samples_dropdown.get()][glycan][adduct]['peaks']):
                        # To filter by quality
                        current_ppm = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['ppm'][index]
                        current_iso = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['iso'][index]
                        current_curve = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['curve'][index]
                        snratio = glycans_per_sample[samples_dropdown.get()][glycan][adduct]['sn'][index]
                
                        # count = 0
                        if current_ppm > max_ppm[1] or current_ppm < max_ppm[0] or current_iso < iso_fit_score or current_curve < curve_fit_score or snratio < s_to_n:
                            bad_peaks.append(f"{glycan}_{adduct}_{peak}")
            
                # Find the boundaries of the bad peaks
                sample_index = list(glycans_per_sample.keys()).index(samples_dropdown.get())
                
                zeroing_intervals = []
                for peak in bad_peaks:
                    glycan, adduct, rt = peak.split("_")
                    if f"{glycan}+{adduct}_{rt}_RTs" in curve_fittings[sample_index]:
                        peak_limits = [curve_fittings[sample_index][f"{glycan}+{adduct}_{rt}_RTs"][0], curve_fittings[sample_index][f"{glycan}+{adduct}_{rt}_RTs"][-1]]
                        zeroing_intervals.append(peak_limits)
                        
                x_data = element.get_xdata()
                if f"{samples_dropdown.get()}_{label}" in original_lines.keys():
                    y_data = copy.deepcopy(original_lines[f"{samples_dropdown.get()}_{label}"])
                else:
                    y_data = element.get_ydata()
                    original_lines[f"{samples_dropdown.get()}_{label}"] = copy.deepcopy(y_data)
                
                for start, end in zeroing_intervals:
                    mask = (x_data >= start) & (x_data <= end)
                    y_data[mask] = 0  # Set values to zero
                    
                element.set_ydata(y_data)
        else:
            for element in graph_lines:
                label = element.get_label()
                
                if label[0] == "_" or "±" in label:
                    continue
                    
                if f"{samples_dropdown.get()}_{label}" in original_lines:
                    element.set_ydata(original_lines[f"{samples_dropdown.get()}_{label}"])
        
        for index, line in enumerate(graph_lines):
            if "±" in line.get_label():
                continue
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            while index > len(colors)-1:
                index -= len(colors)
            color = colors[index]
            line.set_color(color)
            if len(graph_lines) > 1:
                ax.fill_between(x_data, y_data, color=color, alpha=0.25)
                
        if gg_draw_on:
            toggle_gg_draw(chromatograms_list)
            toggle_gg_draw(chromatograms_list)
                
        canvas.draw_idle()
        
    def on_tab_selected(event):
        notebook = event.widget  # the Notebook widget
        tab_id = notebook.select()  # ID of the selected tab
        tab_name = notebook.tab(tab_id, "text")
        if tab_name == 'MS2':
            previous_spectrum_button.place(relx=1, x=-70, y=0, anchor='ne')
            next_spectrum_button.place(relx=1, x=-35, y=0, anchor='ne')
            
            if len([line.get_label() for line in ax_spec_ms2.get_lines() if line.get_label()[0] != "_"]) == 0:
                previous_spectrum_button.config(state=tk.DISABLED)
                next_spectrum_button.config(state=tk.DISABLED)
                annotate_spectrum_button.config(state=tk.DISABLED)
            else:
                previous_spectrum_button.config(state=tk.NORMAL)
                next_spectrum_button.config(state=tk.NORMAL)
            
            if len(reanalysis_path) > 0 and len(fragments_library_ms2_scores) > 0:
                annotate_spectrum_button.place(relx=1, x=-105, y=0, anchor='ne')
                spectrum_score_label.place(relx=1, x=-190, y=0, anchor='ne')
        else:
            try:
                previous_spectrum_button.place_forget()
                next_spectrum_button.place_forget()
                annotate_spectrum_button.place_forget()
                spectrum_score_label.place_forget()
            except:
                pass
                
    def annotating_spectrum(target_glycan = None):
        def close_anw():
            annotating_spectrum_window.destroy()
            
        def wait_thread():
            
            if target_glycan:
                target_glycan_comp = General_Functions.form_to_comp_glycans(target_glycan)
            else:
                target_glycan_comp = None
                
            annotated_fragments = annotate_ms2_spectrum(current_data['access'], int([line.get_label() for line in ax_spec_ms2.get_lines() if line.get_label()[0] != "_"][0].split(" - ")[-1]), fragments_library_ms2_scores[0], tolerance, target_glycan_comp)
            
            markers_list = [line for line in ax_spec_ms2.lines if isinstance(line, matplotlib.lines.Line2D) and line.get_marker() == '*']
            for marker in markers_list:
                marker.remove()
            
            for index, fragment in enumerate(annotated_fragments):
                frag_label = fragment[2]
                ms2_marker = ax_spec_ms2.plot(
                                              fragment[3], 
                                              fragment[4], 
                                              marker='*', 
                                              markersize=4.5, 
                                              label = f"{frag_label}\n{fragment[3]:.4f}", 
                                              color="red"
                                              )
            
            if len(annotated_fragments) > 0:
                percent_tic_annotated = (sum([fragment[4] for fragment in annotated_fragments])/annotated_fragments[0][7])*100
            else:
                percent_tic_annotated = 0
            
            new_precursor_label = precursor_label.cget('text').split(" TIC Annotated:")[0].split(" Composition:")[0]
            new_precursor_label += f" TIC Annotated: {percent_tic_annotated:.1f}%"
            if target_glycan:
                new_precursor_label += f" Composition: {target_glycan}"
            precursor_label.config(text = new_precursor_label)
            
            annotate_top_y_values(ax_spec_ms2, canvas_spec_ms2)
            canvas_spec_ms2.draw_idle()
            
            close_anw()
            
        annotating_spectrum_window = tk.Toplevel()
        # calibrating_spectra_window.attributes("-topmost", True)
        annotating_spectrum_window.withdraw()
        annotating_spectrum_window.title("Annotating Spectrum...")
        icon = ImageTk.PhotoImage(ico_image)
        annotating_spectrum_window.iconphoto(False, icon)
        annotating_spectrum_window.resizable(False, False)
        annotating_spectrum_window.grab_set()
        annotating_spectrum_window.protocol("WM_DELETE_WINDOW", on_closing)
        
        annotating_spectrum_window_label = ttk.Label(annotating_spectrum_window, text="Annotating the peaks of this spectrum. Please wait.", font=("Segoe UI", list_font_size))
        annotating_spectrum_window_label.pack(pady=35, padx=70)
        
        annotating_spectrum_window.update_idletasks()
        annotating_spectrum_window.deiconify()
        window_width = annotating_spectrum_window.winfo_width()
        window_height = annotating_spectrum_window.winfo_height()
        screen_width = annotating_spectrum_window.winfo_screenwidth()
        screen_height = annotating_spectrum_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        annotating_spectrum_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        t = threading.Thread(target=wait_thread)
        t.start()
        
    def move_ms2_spectrum(direction):
        current_spectrum_number = int([line.get_label() for line in ax_spec_ms2.get_lines() if line.get_label()[0] != "_"][0].split(" - ")[-1])
        
        if direction == 'back':
            direction = -1
        else:
            direction = 1
            
        current_spectrum_number += direction
        if current_spectrum_number < 0:
            current_spectrum_number = len(current_data['access'])-1
        elif current_spectrum_number > len(current_data['access'])-1:
            current_spectrum_number = 0
            
        while True:
            current_spectrum = current_data['access'][current_spectrum_number]
            
            try:
                ms_level = current_spectrum['msLevel']
            except:
                ms_level = current_spectrum['ms level']
                
            if ms_level != 2:
                current_spectrum_number += direction
                if current_spectrum_number < 0:
                    current_spectrum_number = len(current_data['access'])-1
                elif current_spectrum_number > len(current_data['access'])-1:
                    current_spectrum_number = 0
            else:
                try:
                    spectrum_time = current_spectrum['retentionTime']
                    precursor_mz = current_spectrum['precursorMz'][0]['precursorMz']
                except:
                    try:
                        spectrum_time = current_spectrum['scanList']['scan'][0]['scan start time']
                        precursor_mz = current_spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']
                    except:
                        spectrum_time = current_spectrum['scanList']['scan'][0]['scan start time']
                        precursor_mz = current_spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                break
                
        show_ms2_graph([precursor_mz, spectrum_time])
        
    def open_target_glycan_window():
        ''''''
        def tg_ok():
            entry_text = target_glycan_entry.get()
            if entry_text != 'All':
                glycan = target_glycan_entry.get()
                annotating_spectrum(glycan)
            else:
                annotating_spectrum()
            target_glycan_window.destroy()
        
        def tg_cancel():
            target_glycan_window.destroy()
        
        target_glycan_window = tk.Toplevel()
        target_glycan_window.withdraw()
        target_glycan_window.title("Target Glycan")
        icon = ImageTk.PhotoImage(ico_image)
        target_glycan_window.iconphoto(False, icon)
        target_glycan_window.resizable(False, False)
        target_glycan_window.grab_set()
        target_glycan_window.protocol("WM_DELETE_WINDOW", tg_cancel)
        
        target_glycan_label = ttk.Label(target_glycan_window, text="If you'd like to force the fragments \ninto a specific glycan composition, \nenter the target glycan name:", font=("Segoe UI", list_font_size))
        target_glycan_label.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        
        target_glycan_entry = ttk.Entry(target_glycan_window)
        target_glycan_entry.insert(0, 'All')
        target_glycan_entry.grid(row=1, column=0, padx=10, pady=(0, 10), sticky='we')
        
        target_glycan_ok_button = ttk.Button(target_glycan_window, text="Ok", style="small_button_spw_style1.TButton", command=tg_ok)
        target_glycan_ok_button.grid(row=2, column=0, padx=(10, 95), pady=(0,10), sticky="e")
        
        target_glycan_cancel_button = ttk.Button(target_glycan_window, text="Cancel", style="small_button_spw_style1.TButton", command=tg_cancel)
        target_glycan_cancel_button.grid(row=2, column=0, padx=(10, 10), pady=(0,10), sticky="e")
                
        target_glycan_window.update_idletasks()
        target_glycan_window.deiconify()
        window_width = target_glycan_window.winfo_width()
        window_height = target_glycan_window.winfo_height()
        screen_width = target_glycan_window.winfo_screenwidth()
        screen_height = target_glycan_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        target_glycan_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
            
    # Create the main window
    global main_window
    main_window = tk.Tk()
    main_window.attributes("-topmost", True)
    main_window.withdraw()
    main_window.grab_set()

    # Configure window properties
    main_window.title("GlycoGenius")
    icon = ImageTk.PhotoImage(ico_image)
    main_window.iconphoto(False, icon)
    main_window.minsize(1025, 720)
    main_window.bind("<Configure>", on_resize)
    main_window.grid_columnconfigure(0, weight=0)
    main_window.grid_columnconfigure(1, weight=1)
    main_window.grid_columnconfigure(2, weight=0)
    main_window.grid_rowconfigure(0, weight=0)
    main_window.grid_rowconfigure(1, weight=1)
    main_window.grid_rowconfigure(2, weight=1)
    main_window.bind("<KeyRelease-Escape>", click_treeview)
    main_window.bind("<KeyPress>", on_key_press)
    main_window.bind("<KeyRelease>", on_key_release)
    main_window.protocol("WM_DELETE_WINDOW", exit_window)
    main_window.geometry("1025x720")
    
    global background_color
    background_color = main_window.cget("background")

    # Load assets
    logo = Image.open(os.path.join(current_dir, "Assets/logo.png"))
    logo_size = logo.size
    logo = logo.resize((int(logo_size[0]/4), int(logo_size[1]/4)))
    tk_logo = ImageTk.PhotoImage(logo)
    
    image_two_d = Image.open(os.path.join(current_dir, "Assets/heatmap_small.png"))
    photo_two_d = ImageTk.PhotoImage(image_two_d)
    
    image_mis = Image.open(os.path.join(current_dir, "Assets/mis.png"))
    photo_mis = ImageTk.PhotoImage(image_mis)
    
    image_eic = Image.open(os.path.join(current_dir, "Assets/eic.png"))
    photo_eic = ImageTk.PhotoImage(image_eic)
    
    global photo_ggdraw, photo_ggdraw_off
    image_ggdraw = Image.open(os.path.join(current_dir, "Assets/gg_draw.png"))
    photo_ggdraw = ImageTk.PhotoImage(image_ggdraw)
    
    image_ggdraw_off = Image.open(os.path.join(current_dir, "Assets/gg_draw_off.png"))
    photo_ggdraw_off = ImageTk.PhotoImage(image_ggdraw_off)
    
    image_zoom_lock = Image.open(os.path.join(current_dir, "Assets/zoom_lock.png"))
    photo_zoom_lock = ImageTk.PhotoImage(image_zoom_lock)
    
    image_previous_arrow = Image.open(os.path.join(current_dir, "Assets/back_arrow.png"))
    photo_previous_arrow = ImageTk.PhotoImage(image_previous_arrow)
    
    image_next_arrow = Image.open(os.path.join(current_dir, "Assets/next_arrow.png"))
    photo_next_arrow = ImageTk.PhotoImage(image_next_arrow)
    
    image_spectra_edit = Image.open(os.path.join(current_dir, "Assets/spectra_edit.png"))
    photo_spectra_edit = ImageTk.PhotoImage(image_spectra_edit)
    
    global photo_ruler
    image_ruler = Image.open(os.path.join(current_dir, "Assets/ruler.png"))
    photo_ruler = ImageTk.PhotoImage(image_ruler)
    
    image_ruler_small = Image.open(os.path.join(current_dir, "Assets/ruler_small.png"))
    photo_ruler_small = ImageTk.PhotoImage(image_ruler_small)
    
    image_remove_ruler_small = Image.open(os.path.join(current_dir, "Assets/remove_ruler_small.png"))
    photo_remove_ruler_small = ImageTk.PhotoImage(image_remove_ruler_small)
    
    banner = Image.open(os.path.join(current_dir, "Assets/banner.png"))
    banner_size = banner.size
    banner = banner.resize((int(banner_size[0]/4), int(banner_size[1]/4)))
    tk_banner = ImageTk.PhotoImage(banner)

    right_arrow = Image.open(os.path.join(current_dir, "Assets/right_arrow.png"))
    right_arrow_size = right_arrow.size
    right_arrow = right_arrow.resize((int(right_arrow_size[0]/2), int(right_arrow_size[1]/2)))
    tk_right_arrow = ImageTk.PhotoImage(right_arrow)

    # Add ttk styles
    chromatograms_list_style = ttk.Style().configure("chromatograms_list.Treeview", bg="white", fg="black", font=("Segoe UI", list_font_size))

    big_button_style = ttk.Style().configure("big_button_style.TButton", font=("Segoe UI", button_font_size), relief="raised", padding = (big_button_size[0], big_button_size[1]), justify="center")

    small_button_style1 = ttk.Style().configure("small_button_style1.TButton", font=("Segoe UI", button_font_size), relief="raised", padding = (int(big_button_size[0]), int(big_button_size[1]/3)), justify="center")

    small_button_style2 = ttk.Style().configure("small_button_style2.TButton", font=("Segoe UI", button_font_size), relief="raised", padding = (int(big_button_size[0]*3.7), int(big_button_size[1]/3)), justify="center")
    
    small_button_style3 = ttk.Style().configure("small_button_style3.TButton", font=("Segoe UI", button_font_size), relief="raised", justify="center")
    
    small_button_style4 = ttk.Style().configure("small_button_style4.TButton", font=("Segoe UI", button_font_size), relief="raised", justify="center")

    about_button_style = ttk.Style().configure("about_button_style.TButton", font=("Segoe UI", button_font_size), relief="raised", padding = (0, int(big_button_size[1]*1.45)), justify="center")
    
    annotate_spectrum_button_style = ttk.Style().configure("annotate_spectrum_button_style.TButton", font=("Segoe UI", button_font_size-4), relief="raised", justify="center")
    
    two_d_button_style = ttk.Style().configure("two_d_button_style.TButton", font=("Segoe UI", button_font_size), relief="raised", padding = (0, 0), justify="center")
    
    spectrum_ruler_button_style = ttk.Style().configure("spectrum_ruler_button_style.TButton", font=("Segoe UI", button_font_size), relief="raised", padding = (0, 0), justify="center", background = background_color)

    chromatogram_plot_frame_style = ttk.Style().configure("chromatogram.TLabelframe", font=("Segoe UI", list_font_size))
    
    qcp_frame_style = ttk.Style().configure("qcp_frame.TLabelframe", font=("Segoe UI", list_font_size_smaller))

    # First row of widgets
    banner_label = tk.Label(main_window, image=tk_logo)
    banner_label.grid(row=0, column=0, sticky='w')
    
    global select_files_frame
    select_files_frame = tk.Frame(main_window, bd=3, relief="flat")
    select_files_frame.grid(row=0, column=0, columnspan=3, padx=(160, 0), sticky='w')
    select_files = ttk.Button(select_files_frame, text="Select\nFiles", style="big_button_style.TButton", command=open_select_files_window)
    select_files.pack(padx=0, pady=0)
    ToolTip(select_files, "Select the MzML or MzXML files to analyze, a .gg file to check the results of a completed analysis, or both of them to access all the data together.")

    right_arrow_label1 = tk.Label(main_window, image=tk_right_arrow)
    right_arrow_label1.grid(row=0, column=0, columnspan=3, padx=(275, 0), sticky='w')

    global set_parameters_frame
    set_parameters_frame = tk.Frame(main_window, bd=3, relief="flat")
    set_parameters_frame.grid(row=0, column=0, columnspan=3, padx=(305, 0), sticky='w')
    set_parameters = ttk.Button(set_parameters_frame, text="Set\nParameters", style="big_button_style.TButton", command=open_set_parameters_window)
    set_parameters.pack(padx=0, pady=0)
    ToolTip(set_parameters, "Set the parameters for the library building and the analysis here.")

    right_arrow_label2 = tk.Label(main_window, image=tk_right_arrow)
    right_arrow_label2.grid(row=0, column=0, columnspan=3, padx=(420, 0), pady=42, sticky='nw')

    right_arrow_label3 = tk.Label(main_window, image=tk_right_arrow)
    right_arrow_label3.grid(row=0, column=0, columnspan=3, padx=(420, 0), pady=42, sticky='sw')

    global generate_library
    global generate_library_button_frame
    generate_library_button_frame = tk.Frame(main_window, bd=3, relief="flat")
    generate_library_button_frame.grid(row=0, column=0, columnspan=3, padx=(450, 0), sticky='nw', pady=29)
    generate_library = ttk.Button(generate_library_button_frame, text="Generate Library", style="small_button_style2.TButton", command=lambda: file_name_window(file_type = 'library'))
    generate_library.pack(padx=0, pady=0)
    ToolTip(generate_library, "Generates a new library based on the parameters set on Set Parameters menu.")

    global import_library
    global import_library_button_frame
    import_library_button_frame = tk.Frame(main_window, bd=3, relief="flat")
    import_library_button_frame.grid(row=0, column=0, columnspan=3, padx=(450, 0), sticky='sw', pady=29)
    import_library = ttk.Button(import_library_button_frame, text="Import\nLibrary", style="small_button_style3.TButton", command=open_file_dialog_import_button)
    import_library.pack(padx=0, pady=0)
    ToolTip(import_library, "Import a previously generated library in .ggl format.")
    
    global import_library_info
    global import_library_info_button_frame
    import_library_info_button_frame = tk.Frame(main_window, bd=3, relief="flat")
    import_library_info_button_frame.grid(row=0, column=0, columnspan=3, padx=(542, 0), sticky='sw', pady=29)
    import_library_info = ttk.Button(import_library_info_button_frame, text="Check\nInfo", style="small_button_style4.TButton", command=get_lib_info, state=tk.DISABLED)
    import_library_info.pack(padx=0, pady=0)
    ToolTip(import_library_info, "Click to access information regarding how the imported library was built.")

    right_arrow_label4 = tk.Label(main_window, image=tk_right_arrow)
    right_arrow_label4.grid(row=0, column=0, columnspan=3, padx=(640, 0), sticky='w')
    
    global run_analysis_button_frame
    global run_analysis_button
    run_analysis_button_frame = tk.Frame(main_window, bd=3, relief="flat")
    run_analysis_button_frame.grid(row=0, column=0, columnspan=3, padx=(670, 0), sticky='w')
    run_analysis_button = ttk.Button(run_analysis_button_frame, text="Run\nAnalysis", style="big_button_style.TButton", command=lambda: file_name_window(file_type = 'gg'))
    run_analysis_button.pack(padx=0, pady=0)
    ToolTip(run_analysis_button, "Run an analysis on the selected MzML and/or MzXML files you selected in the Select Files menu, with the parameters set on the Set Parameters menu.")

    right_arrow_label5 = tk.Label(main_window, image=tk_right_arrow)
    right_arrow_label5.grid(row=0, column=0, columnspan=3, padx=(785, 0), sticky='w')

    global save_results_button_frame
    save_results_button_frame = tk.Frame(main_window, bd=3, relief="flat")
    save_results_button_frame.grid(row=0, column=0, columnspan=3, padx=(815, 0), sticky='w')
    save_results = ttk.Button(save_results_button_frame, text="Save\nResults", style="big_button_style.TButton", command=save_results_button_command)
    save_results.pack(padx=0, pady=0)
    ToolTip(save_results, "Exports the results to files that can be opened on excel, in a comprehensible format.")

    about = ttk.Button(main_window, text="About", style="about_button_style.TButton", command=about_button_command)
    about.grid(row=0, column=2, padx=(10,10), sticky='e')
    ToolTip(about, "More information about GlycoGenius.")
    
    global two_d
    two_d = ttk.Button(main_window, image=photo_two_d, style="two_d_button_style.TButton", command=backend_heatmap, state=tk.DISABLED)
    two_d.grid(row=0, column=2, padx=(10,55), sticky='se')
    ToolTip(two_d, "Creates a 2D-Plot (Retention/Migration Time x m/z) of the current axis ranges of the displayed chromatogram/electropherogram and spectrum in a new window. Wider ranges will take longer to load.")
    
    global quick_check
    quick_check = ttk.Button(main_window, image=photo_mis, style="two_d_button_style.TButton", command=quick_check_window, state=tk.DISABLED)
    quick_check.grid(row=0, column=2, padx=(10,90), sticky='se')
    ToolTip(quick_check, "Calculates an aggregated spectra called Maximum Intensity Spectrum (MIS) of the whole chromatographic/electropherographic run based on the maximum intensity values each m/z achieves. The sample can be quickly checked for the presence of glycans on the MIS window, if you have a library loaded.\nThe first time you click on this button for a given sample it will take up to a few minutes to calculate the MIS. Subsequent times will load instantly.")
    
    global quick_trace
    quick_trace = ttk.Button(main_window, image=photo_eic, style="two_d_button_style.TButton", command=run_quick_trace_window, state=tk.NORMAL)
    quick_trace.grid(row=0, column=2, padx=(10,125), sticky='se')
    ToolTip(quick_trace, "Opens a small window that allows you to create custom Extracted Ion Chromatograms/Electropherograms for your samples.")
    
    global spectra_file_editing
    spectra_file_editing = ttk.Button(main_window, image=photo_spectra_edit, style="two_d_button_style.TButton", command=lambda: check_mzml_window([samples_dropdown, samples_list, library_path, [ax_spec, canvas_spec, scaling_dropdown, rt_label], processed_data]), state=tk.DISABLED)
    spectra_file_editing.grid(row=0, column=2, padx=(10,20), sticky='se')
    ToolTip(spectra_file_editing, "Opens a window that allows you to calibrate, trim or align your spectra files, creating new files.")
    
    global gg_draw
    gg_draw = ttk.Button(main_window, image=photo_ggdraw, style="two_d_button_style.TButton", command=lambda: toggle_gg_draw(chromatograms_list), state=tk.NORMAL)
    gg_draw.grid(row=0, column=2, padx=(10,160), sticky='se')
    ToolTip(gg_draw, "Display on the chromatogram/electropherogram viewer the glycans cartoons. The initial glycan structures are chosen randomly. Double click a structure to change or remove the displayed structure. Click and drag on a structure to reposition it in the plot.")
    
    # Create panned window for bottom side widgets
    paned_window = tk.PanedWindow(main_window, sashwidth=5, orient=tk.HORIZONTAL)
    paned_window.grid(row=1, column=0, columnspan=3, padx = 10, pady = 10, sticky='nwse')
    
    # Create the left side widgets frame in the main window paned window
    left_side_widgets_frame = ttk.Frame(paned_window)
    left_side_widgets_frame.grid_columnconfigure(0, weight=1)
    left_side_widgets_frame.grid_columnconfigure(1, weight=1)
    left_side_widgets_frame.grid_rowconfigure(0, weight=0)
    left_side_widgets_frame.grid_rowconfigure(1, weight=0)
    left_side_widgets_frame.grid_rowconfigure(2, weight=1)
    left_side_widgets_frame.grid_rowconfigure(3, weight=0)
    left_side_widgets_frame.grid_rowconfigure(4, weight=0)
    
    # Second row of widgets
    global samples_dropdown_options, samples_dropdown, chromatograms_list, selected_item
    samples_dropdown = ttk.Combobox(left_side_widgets_frame, state="readonly", values=samples_dropdown_options)
    samples_dropdown.grid(row=0, column=0, columnspan=2, padx = 10, sticky='new')
    samples_dropdown.bind("<<ComboboxSelected>>", handle_selection)
    ToolTip(samples_dropdown, "")
    
    global filter_list
    filter_list = tk.Entry(left_side_widgets_frame, fg='grey', width=6)
    filter_list.grid(row=1, column=0, columnspan=2, padx=10, pady=(0, 0), sticky='new')
    filter_list.bind("<KeyRelease>", lambda event: on_key_release_filter(event, filter_list))
    filter_list.insert(0, "Filter the list of glycans...")
    ToolTip(filter_list, "Type here to filter the glycans list. Can filter by the glycan name, by their QC threshold status (e.g.: 'good', 'bad', 'average'), by whether or not they have MS2 and a combination of these, by connecting them with a '+' sign (e.g.: 'good+ms2+S2' will filter the list with all glycans that meet all thresholds of quality, have MS2 annotated and contain two sialic acids).")
    filter_list.bind("<FocusIn>", on_focus_in_filter_list)
    filter_list.bind("<FocusOut>", on_focus_out_filter_list)

    chromatograms_list_scrollbar = tk.Scrollbar(left_side_widgets_frame, orient=tk.VERTICAL)
    chromatograms_list = ttk.Treeview(left_side_widgets_frame, height=100, style="chromatograms_list.Treeview", yscrollcommand=chromatograms_list_scrollbar.set)
        
    chromatograms_list.tag_configure('bad', background='#FED5CD')
    chromatograms_list.tag_configure('warning', background='#FEFACD')
    chromatograms_list.tag_configure('good', background='#E3FECD')
    
    chromatograms_list["show"] = "tree" #removes the header
    chromatograms_list["columns"] = ("#1")
    
    chromatograms_list.column("#0", width=230)
    chromatograms_list.column("#1", width=35, stretch=False) #this column is for showing ambiguities
    
    chromatograms_list_scrollbar.config(command=chromatograms_list.yview, width=10)
    chromatograms_list.grid(row=2, column=0, columnspan=2, padx=10, pady=(0, 0), sticky="nsew")
    chromatograms_list_scrollbar.grid(row=2, column=0, columnspan=2, pady=(0, 0), sticky="nse")
    
    chromatograms_list.bind("<KeyRelease-Up>", handle_treeview_select)
    chromatograms_list.bind("<KeyRelease-Down>", handle_treeview_select)
    chromatograms_list.bind("<ButtonRelease-1>", click_treeview)
    chromatograms_list.bind("<Double-Button-1>", handle_double_left_click)
    chromatograms_list.bind("<Motion>", on_treeview_motion)
    
    global compare_samples_button
    compare_samples_button = ttk.Button(left_side_widgets_frame, text="Compare samples", style="small_button_style1.TButton", command=compare_samples_window, state=tk.DISABLED)
    compare_samples_button.grid(row=3, column=0, padx=(10, 0), pady=(0, 0), sticky="sew")
    ToolTip(compare_samples_button, "Opens a window for comparing the chromatograms/electropherograms for the selected compound on different samples. It features an option for alignment of the chromatograms/electropherograms that can be toggled on or off.")
    
    global plot_graph_button
    plot_graph_button = ttk.Button(left_side_widgets_frame, text="Plot Graph", style="small_button_style1.TButton", command=plot_graph_window, state=tk.DISABLED)
    plot_graph_button.grid(row=3, column=1, padx=(0, 10), pady=(0, 0), sticky="sew")
    ToolTip(plot_graph_button, "Plots graphs of the selected glycans' abundance. If one glycan is selected and more than one sample is loaded, plots comparison between samples. If more than one glycan is selected, plots comparison between selected glycans within the same sample. If no glycans are selected, plots all the 'good' glycans for the current sample.")
    
    # Quality criteria parameters frame
    global s_n_entry, curve_fit_entry, ppm_error_min_entry, ppm_error_max_entry, iso_fit_entry
    qcp_frame = ttk.Labelframe(left_side_widgets_frame, text="Quality Scores Thresholds:", style="qcp_frame.TLabelframe")
    qcp_frame.grid(row=4, column=0, columnspan=2, padx=10, pady=(0, 0), sticky="sew")
    
    # Adjust column resizing settings
    qcp_frame.grid_columnconfigure(0, weight=1)
    qcp_frame.grid_columnconfigure(1, weight=0)
    qcp_frame.grid_rowconfigure(0, weight=0)
    qcp_frame.grid_rowconfigure(1, weight=0)
    qcp_frame.grid_rowconfigure(2, weight=0)
    qcp_frame.grid_rowconfigure(3, weight=0)
    qcp_frame.grid_rowconfigure(4, weight=0)
    
    global check_qc_dist_button
    check_qc_dist_button = ttk.Button(qcp_frame, text="Check Scores Distribution", style="small_button_style1.TButton", command=check_qc_dist, state=tk.DISABLED)
    check_qc_dist_button.grid(row=0, column=0, columnspan=2, padx=10, pady=(5, 0), sticky="new")
    ToolTip(check_qc_dist_button, "Plots all the quality scores of all the peaks identified for all the glycans in the analysis. This can work as a base for setting the Quality Scores Criteria below.")
    
    iso_fit_label = ttk.Label(qcp_frame, text='Minimum Isotopic Fitting Score:', font=("Segoe UI", list_font_size_smaller))
    iso_fit_label.grid(row=1, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(iso_fit_label, "Insert here the minimum isotopic fitting score for peaks. Values allowed from 0.0 to 1.0.")
    
    iso_fit_entry = ttk.Spinbox(qcp_frame, width=5, from_=0, to=1.0, increment=0.01)
    iso_fit_entry.insert(0, iso_fit_score)
    iso_fit_entry.grid(row=1, column=1, padx=(5, 10), pady=(5, 0), sticky='e')
    ToolTip(iso_fit_entry, "Insert here the minimum isotopic fitting score for peaks. Values allowed from 0.0 to 1.0.")
    iso_fit_entry.bind("<KeyRelease>", qcp_enter)
    iso_fit_entry.bind("<<Increment>>", qcp_enter)
    iso_fit_entry.bind("<<Decrement>>", qcp_enter)
    
    curve_fit_label = ttk.Label(qcp_frame, text='Minimum Curve Fitting Score:', font=("Segoe UI", list_font_size_smaller))
    curve_fit_label.grid(row=2, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(curve_fit_label, "Insert here the minimum curve-fitting score for peaks. Values allowed from 0.0 to 1.0.")
    
    curve_fit_entry = ttk.Spinbox(qcp_frame, width=5, from_=0, to=1.0, increment=0.01)
    curve_fit_entry.insert(0, curve_fit_score)
    curve_fit_entry.grid(row=2, column=1, padx=(5, 10), pady=(5, 0), sticky='e')
    ToolTip(curve_fit_entry, "Insert here the minimum curve-fitting score for peaks. Values allowed from 0.0 to 1.0.")
    curve_fit_entry.bind("<KeyRelease>", qcp_enter)
    curve_fit_entry.bind("<<Increment>>", qcp_enter)
    curve_fit_entry.bind("<<Decrement>>", qcp_enter)
    
    s_n_label = ttk.Label(qcp_frame, text='Minimum Signal-to-Noise Ratio:', font=("Segoe UI", list_font_size_smaller))
    s_n_label.grid(row=3, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(s_n_label, "Insert here the minimum amount of signal-to-noise ratio. Values under 1.0 won't make a difference.")
    
    s_n_entry = ttk.Spinbox(qcp_frame, width=5, from_=0, to=9999, increment=1)
    s_n_entry.insert(0, s_to_n)
    s_n_entry.grid(row=3, column=1, padx=(5, 10), pady=(5, 0), sticky='e')
    ToolTip(s_n_entry, "Insert here the minimum amount of signal-to-noise ratio. Values under 1.0 won't make a difference.")
    s_n_entry.bind("<KeyRelease>", qcp_enter)
    s_n_entry.bind("<<Increment>>", qcp_enter)
    s_n_entry.bind("<<Decrement>>", qcp_enter)
    
    ppm_error_label = ttk.Label(qcp_frame, text='Min/Max PPM Error:', font=("Segoe UI", list_font_size_smaller))
    ppm_error_label.grid(row=4, column=0, padx=(10, 10), pady=(5, 10), sticky="w")
    ToolTip(ppm_error_label, "Insert here the PPM error.")
    
    ppm_error_min_entry = ttk.Spinbox(qcp_frame, width=5, from_=-999, to=999, increment=1)
    ppm_error_min_entry.insert(0, max_ppm[0])
    ppm_error_min_entry.grid(row=4, column=0, padx=(10, 0), pady=(5, 10), sticky="e")
    ToolTip(ppm_error_min_entry, "Insert here the minimum PPM error.")
    ppm_error_min_entry.bind("<KeyRelease>", qcp_enter)
    ppm_error_min_entry.bind("<<Increment>>", qcp_enter)
    ppm_error_min_entry.bind("<<Decrement>>", qcp_enter)
    
    ppm_error_hyphen_label = ttk.Label(qcp_frame, text='-', font=("Segoe UI", list_font_size_smaller))
    ppm_error_hyphen_label.grid(row=4, column=1, padx=(0, 60), pady=(5, 10), sticky="e")
    ToolTip(ppm_error_hyphen_label, "Insert here the PPM error.")
    
    ppm_error_max_entry = ttk.Spinbox(qcp_frame, width=5, from_=-999, to=999, increment=1)
    ppm_error_max_entry.insert(0, max_ppm[1])
    ppm_error_max_entry.grid(row=4, column=1, padx=(5, 10), pady=(5, 10), sticky='e')
    ToolTip(ppm_error_max_entry, "Insert here the maximum PPM error.")
    ppm_error_max_entry.bind("<KeyRelease>", qcp_enter)
    ppm_error_max_entry.bind("<<Increment>>", qcp_enter)
    ppm_error_max_entry.bind("<<Decrement>>", qcp_enter)
    
    global chromatograms_qc_numbers
    chromatograms_qc_numbers = ttk.Label(left_side_widgets_frame, text=f"Compositions Quality:\n        Good: {0}    Average: {0}    Bad: {0}\n        Ambiguities: {0}", font=("Segoe UI", list_font_size_smaller))
    chromatograms_qc_numbers.grid(row=5, column=0, columnspan=2, padx=10, pady=(0, 0), sticky="sew")
    ToolTip(chromatograms_qc_numbers, "Good compositions have at least one peak that matches all quality scores criteria set above; Average have at least one peak that fails only one criteria; Bad have all peaks failing at least two criterias.")
    
    # Add a new paned window to the right_side_widgets_frame
    paned_window_plots = tk.PanedWindow(paned_window, sashwidth=5, orient=tk.VERTICAL)

    chromatogram_plot_frame = ttk.Labelframe(paned_window_plots, text="Chromatogram/Electropherogram Viewer", style="chromatogram.TLabelframe", height = 400)
    chromatogram_plot_frame.pack(fill=tk.BOTH, expand=True)
    
    hide_bad_peaks_checkbox_state = tk.BooleanVar(value=True)
    hide_bad_peaks_checkbox = ttk.Checkbutton(chromatogram_plot_frame, text="Hide peaks below QC thresholds", variable=hide_bad_peaks_checkbox_state, command=toggle_hide_bad_peaks)
    hide_bad_peaks_checkbox.place(relx = 0, rely = 0)
    ToolTip(hide_bad_peaks_checkbox, "Hides the peaks that don't fit the quality scores thresholds you've set.")
    
    zoom_lock_button = ttk.Button(chromatogram_plot_frame, image=photo_zoom_lock, style="spectrum_ruler_button_style.TButton", command=lambda: create_zoomlock_window('chromatogram'), state=tk.NORMAL)
    zoom_lock_button.pack(side="top", anchor="ne")
    ToolTip(zoom_lock_button, "Locks the zoom of the chromatogram viewer to selected ranges.")

    global canvas, ax, coordinate_label, type_coordinate, scaling_dropdown
    fig = plt.figure(figsize=(0, 4))
    ax = fig.add_subplot(111)
    canvas = FigureCanvasTkAgg(fig, master=chromatogram_plot_frame)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    coordinate_label = tk.Label(chromatogram_plot_frame, text="", anchor="e", font=("Segoe UI", 8), bg="white")
    coordinate_label.place(relx=1.0, y=20, anchor='ne')
    type_coordinate = "chromatogram"
    coordinate_label.lift()
    
    ax.set_position([0.0944, 0.1209, 0.8909, 0.8252])
    chromatogram_plot_frame.bind("<Configure>", lambda event, ax=ax: adjust_subplot_size(event, ax, canvas))
    canvas.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax, canvas, True) if event.button == 3 else None)
    
    global spectra_plot_frame
    spectra_plot_frame = ttk.Labelframe(paned_window_plots, text="Spectra Viewer", style="chromatogram.TLabelframe")
    spectra_plot_frame.pack(fill=tk.BOTH, expand=True)
    
    # Create a notebook for the spectra
    spectra_notebook = ttk.Notebook(spectra_plot_frame)
    spectra_notebook.pack(fill=tk.BOTH, expand=True)
    spectra_notebook.bind("<<NotebookTabChanged>>", on_tab_selected)
    
    # Add dropdown menu for spectra plot scale
    scaling_dropdown_options = ['Linear', 'Sqrt', 'Log']
    scaling_dropdown = ttk.Combobox(spectra_plot_frame, state="readonly", values=scaling_dropdown_options, width=7)
    scaling_dropdown.set('Linear')
    scaling_dropdown.place(x=70, y=0, anchor='nw')
    scaling_dropdown.bind("<<ComboboxSelected>>", on_scale_selected)
    ToolTip(scaling_dropdown, "Allows you to select different scalings for the y-axis of the spectra.")
    
    # Ruler button
    global spectrum_ruler_on, spectrum_ruler_frame, spectrum_rulers
    spectrum_ruler_on = False
    spectrum_rulers = []
    spectrum_ruler_frame = tk.Frame(spectra_plot_frame, bd=2, relief="flat")
    spectrum_ruler_frame.place(x=140, y=0, anchor='nw')
    spectrum_ruler_button = ttk.Button(spectrum_ruler_frame, image=photo_ruler_small, style="spectrum_ruler_button_style.TButton", command=toggle_spectrum_ruler, state=tk.NORMAL)
    spectrum_ruler_button.pack(padx=0, pady=0)
    ToolTip(spectrum_ruler_button, "Activate to draw rulers on the spectrum plot and measure the distance between different peaks or to select existing ones.")
    
    clear_spectrum_ruler_button = ttk.Button(spectra_plot_frame, image=photo_remove_ruler_small, style="spectrum_ruler_button_style.TButton", command=remove_spectrum_ruler, state=tk.NORMAL)
    clear_spectrum_ruler_button.place(x=175, y=2, anchor='nw')
    ToolTip(clear_spectrum_ruler_button, "Delete the currently selected ruler or clear them all if none is selected.")
    
    # MS1 plotting field on the notebook
    ms1_plot_frame = tk.Frame(spectra_notebook)
    ms1_plot_frame.pack(fill=tk.BOTH, expand=True)
    spectra_notebook.add(ms1_plot_frame, text="MS1")
    
    # MS2 plotting field on the notebook
    ms2_plot_frame = tk.Frame(spectra_notebook)
    ms2_plot_frame.pack(fill=tk.BOTH, expand=True)
    spectra_notebook.add(ms2_plot_frame, text="MS2")
    
    # MS1 plotting elements
    spectra_plot_frame_canvas = tk.Canvas(ms1_plot_frame, bg="white")
    spectra_plot_frame_canvas.place(relwidth=1, relheight=1)
    
    global canvas_spec, ax_spec, coordinate_label_spec, type_coordinate_spec, no_spectra_loaded_label, rt_label
    fig_spec = plt.figure(figsize=(0, 4))
    ax_spec = fig_spec.add_subplot(111)
    canvas_spec = FigureCanvasTkAgg(fig_spec, master=ms1_plot_frame)
    canvas_spec.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    spectrum_score_label = tk.Label(spectra_plot_frame, text="", anchor="e", font=("Segoe UI", 8))
    
    global annotate_spectrum_button, previous_spectrum_button, next_spectrum_button
    annotate_spectrum_button = ttk.Button(spectra_plot_frame, text="Annotate Spectrum", style="annotate_spectrum_button_style.TButton", command=open_target_glycan_window)
    ToolTip(annotate_spectrum_button, "Annotates the peaks in this spectrum with the fragments library generated during the analysis.")
    
    previous_spectrum_button = ttk.Button(spectra_plot_frame, image=photo_previous_arrow, style="spectrum_ruler_button_style.TButton", command=lambda:move_ms2_spectrum('back'))
    ToolTip(previous_spectrum_button, "Go back one MS2 spectrum.")
    
    next_spectrum_button = ttk.Button(spectra_plot_frame, image=photo_next_arrow, style="spectrum_ruler_button_style.TButton", command=lambda:move_ms2_spectrum('forward'))
    ToolTip(next_spectrum_button, "Advance one MS2 spectrum.")
    
    zoom_lock_spec_button = ttk.Button(spectra_plot_frame, image=photo_zoom_lock, style="spectrum_ruler_button_style.TButton", command=lambda: create_zoomlock_window('spectrum'), state=tk.NORMAL)
    zoom_lock_spec_button.place(relx=1, y=0, anchor='ne')
    ToolTip(zoom_lock_spec_button, "Locks the zoom of the spectra viewer to selected ranges.")
    
    coordinate_label_spec = tk.Label(spectra_plot_frame, text="", anchor="e", font=("Segoe UI", 8), bg="white")
    coordinate_label_spec.place(relx=0.99, y=24, anchor='ne')
    type_coordinate_spec = "spectra"
    
    rt_label = tk.Label(ms1_plot_frame, text=f"", font=("Segoe UI", 8), bg="white")
    rt_label.place(relx=0.5, y=0, anchor='n')
    
    no_spectra_loaded_label = tk.Label(ms1_plot_frame, text=f"", anchor="center", font=("Segoe UI", 14), bg="white")
    no_spectra_loaded_label.place(relx=0.01, rely=0.9, anchor='center')
    
    coordinate_label_spec.lift()
    rt_label.lift()
        
    ax_spec.set_position([0.0944, 0.1209, 0.8909, 0.8052])
    ms1_plot_frame.bind("<Configure>", lambda event, ax_spec=ax_spec: adjust_subplot_size(event, ax_spec, canvas_spec))
    canvas_spec.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_spec, canvas_spec, True) if event.button == 3 else None)
    
    # MS2 plotting elements
    spectra_plot_ms2_frame_canvas = tk.Canvas(ms2_plot_frame, bg="white")
    spectra_plot_ms2_frame_canvas.place(relwidth=1, relheight=1)
    
    global canvas_spec_ms2, ax_spec_ms2, no_spectra_loaded_label_ms2, precursor_label
    fig_spec_ms2 = plt.figure(figsize=(0, 4))
    ax_spec_ms2 = fig_spec_ms2.add_subplot(111)
    canvas_spec_ms2 = FigureCanvasTkAgg(fig_spec_ms2, master=ms2_plot_frame)
    canvas_spec_ms2.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    precursor_label = tk.Label(ms2_plot_frame, text=f"", font=("Segoe UI", 8), bg="white")
    precursor_label.place(relx=0.5, y=0, anchor='n')
    
    no_spectra_loaded_label_ms2 = tk.Label(ms2_plot_frame, text=f"", anchor="center", font=("Segoe UI", 14), bg="white")
    no_spectra_loaded_label_ms2.place(relx=0.01, rely=0.9, anchor='center')
        
    ax_spec_ms2.set_position([0.0944, 0.1209, 0.8909, 0.8052])
    ms2_plot_frame.bind("<Configure>", lambda event, ax_spec_ms2=ax_spec_ms2: adjust_subplot_size(event, ax_spec_ms2, canvas_spec_ms2))
    canvas_spec_ms2.mpl_connect('button_press_event', lambda event: on_right_click_plot(event, ax_spec_ms2, canvas_spec_ms2, True) if event.button == 3 else None)
    
    # Add the chromatogram plot frame to the paned_window_plots
    paned_window_plots.add(chromatogram_plot_frame)
    paned_window_plots.paneconfigure(chromatogram_plot_frame, minsize=300)
    
    # Add the spectra plot frame to the paned_window_plots
    paned_window_plots.add(spectra_plot_frame)
    paned_window_plots.paneconfigure(spectra_plot_frame, minsize=300)
    
    # Add the left side widgets frame to the paned window
    paned_window.add(left_side_widgets_frame)
    paned_window.paneconfigure(left_side_widgets_frame, minsize=280)
    
    # Add the right side widgets frame to the paned window
    paned_window.add(paned_window_plots)
    paned_window.paneconfigure(paned_window_plots, minsize=600)
    
    main_window.deiconify()
    main_window.attributes("-topmost", False)
    
    # Check file associations
    if platform.system() == "Windows":
        association_exists = check_file_association(".gg", "ggfile")
        if not association_exists[0]:
            if association_exists[1] == 2:
                create_new_association = warning_window("Former association of analysis result files (.gg) created by GlycoGenius is broken. Do you wish to fix it?")
            elif association_exists[1] == 3:
                create_new_association = warning_window("Analysis result files (.gg) created by GlycoGenius are associated with another installation of GlycoGenius. Do you wish to change it to this installation?")
            else:
                create_new_association = warning_window("Analysis result files (.gg) created by GlycoGenius are currently not associated to open with GlycoGenius. Do you wish to create this association?")
            if create_new_association:
                if association_exists[1] == 2 or association_exists[1] == 3:
                    remove_file_association(".gg", "ggfile")
                create_file_association(".gg", "ggfile", "GlycoGenius Analysis Results File")
                
        # Check if running GG from .gg file
        if len(sys.argv) > 1:
            if sys.argv[1].split(".")[-1] == 'gg':
                global reanalysis_path
                try:
                    test_access = gg_archive(sys.argv[1])
                    if str(test_access) == 'File is unloaded.':
                        error_window("Invalid reanalysis file. Check it and try again.")
                        exit_window(True)
                    else:
                        global reanalysis_file
                        reanalysis_path = sys.argv[1]
                        reanalysis_file = threading.Thread(target=load_reanalysis, args=(reanalysis_path,))
                        reanalysis_file.start()
                        select_files_frame.config(bg='lightgreen')
                        loading_files_from_explorer()
                except:
                    error_window("Invalid reanalysis file. Check it and try again.")
    
    # Start the event loop
    main_window.mainloop()
    
def run_select_files_window(samples_dropdown):
    '''
    '''
    # Functions used by this window
    def close_sf_window():
        select_files_window.destroy()
        
    def open_file_dialog_sfw():
        file_dialog = tk.Toplevel()
        file_dialog.attributes("-topmost", True)
        file_dialog.withdraw()
        file_dialog.grab_set()
        
        file_path = filedialog.askopenfilenames(filetypes=[("MzXML or MzML files", "*.mzML *.mzXML"), ("All files", "*.*")])
            
        for i in file_path:
            files_list.insert('', 'end', text=i)
        longest_len = 0
        for i in files_list.get_children():
            if len(files_list.item(i, "text")) > longest_len:
                longest_len = len(files_list.item(i, "text"))
        if longest_len*7 > 448:
            files_list.column("#0", width=longest_len*7)
        
        if len(files_list.get_children()) > 0:
            if len(files_list.get_children()) == 1:
                file_amount_label.config(text=f'{len(files_list.get_children())} file selected')
            else:
                file_amount_label.config(text=f'{len(files_list.get_children())} files selected')
        else:
            file_amount_label.config(text='')
            
        file_dialog.destroy()
        select_files_window.grab_set()
            
    def open_file_dialog_sfw_reanalysis():
        file_dialog = tk.Toplevel()
        file_dialog.attributes("-topmost", True)
        file_dialog.withdraw()
        file_dialog.grab_set()
        
        file_path = filedialog.askopenfilename(filetypes=[("Glycogenius files", "*.gg"), ("All files", "*.*")])
        if len(file_path) > 0:
            file_info_button.config(state=tk.NORMAL)
        else:
            file_info_button.config(state=tk.DISABLED)
        gg_file_label.config(text=file_path)
        
        file_dialog.destroy()
        select_files_window.grab_set()
            
    def remove_selected_item_sfw():
        selected_items = files_list.selection()  # Get the selected item(s)
        for item in selected_items:
            files_list.delete(item)
        
        if len(files_list.get_children()) > 0:
            if len(files_list.get_children()) == 1:
                file_amount_label.config(text=f'{len(files_list.get_children())} file selected')
            else:
                file_amount_label.config(text=f'{len(files_list.get_children())} files selected')
        else:
            file_amount_label.config(text='')
            
    def loading_files():
        global quick_check
        def close_lf():
            loading_files.destroy()
            
        def wait_thread():
            if len(samples_list) > 0:
                processed_samples.join()
            if len(reanalysis_path) > 0:
                reanalysis_file.join()
            if len(samples_dropdown['values']) != 0:
                samples_dropdown.current(0)
                handle_selection(None)
            else:
                global ms1_bound, ms1_binds, ms2_bound, ms2_binds, chromatograms_list, chromatogram_bound, chromatogram_binds
                samples_dropdown.set('')
                chromatograms_list.delete(*chromatograms_list.get_children())
                clear_plot(ax, canvas)
                clear_plot(ax_spec, canvas_spec)
                clear_plot(ax_spec_ms2, canvas_spec_ms2)
                
                # Unbind MS1 plot
                if ms1_bound:
                    for i in ms1_binds:
                        canvas_spec.mpl_disconnect(i)
                    ms1_bound = False
                    ms1_binds = []
                    
                # Unbind MS2 plot
                if ms2_bound:
                    for i in ms2_binds:
                        canvas_spec_ms2.mpl_disconnect(i)
                    ms2_bound = False
                    ms2_binds = []
                    
                # Unbind chromatogram plot
                if chromatogram_bound:
                    for i in chromatogram_binds:
                        canvas.mpl_disconnect(i)
                    chromatogram_bound = False
                    chromatogram_binds = []
                    
            close_lf()
        
        global loading_files
        loading_files = tk.Toplevel()
        # loading_files.attributes("-topmost", True)
        loading_files.withdraw()
        loading_files.title("Loading Files")
        icon = ImageTk.PhotoImage(ico_image)
        loading_files.iconphoto(False, icon)
        loading_files.resizable(False, False)
        loading_files.grab_set()
        loading_files.protocol("WM_DELETE_WINDOW", on_closing)
        
        loading_files_label = ttk.Label(loading_files, text="Loading files, please wait.", font=("Segoe UI", list_font_size))
        loading_files_label.pack(pady=35, padx=70)
        
        loading_files.update_idletasks()
        loading_files.deiconify()
        window_width = loading_files.winfo_width()
        window_height = loading_files.winfo_height()
        screen_width = loading_files.winfo_screenwidth()
        screen_height = loading_files.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        loading_files.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        t = threading.Thread(target=wait_thread)
        t.start()
            
    def ok_button_sfw():
        global reanalysis_path, samples_list, samples_names, samples_dropdown_options, two_d, quick_check, spectra_file_editing, compare_samples_button, check_qc_dist_button, filter_list, former_alignments, plot_graph_button, processed_samples, reanalysis_file, former_selected_sample
        
        # Get the list of item identifiers (IDs) in the Treeview
        item_ids = files_list.get_children()
        
        # Check integrity of gg file
        if len(gg_file_label.cget("text")) > 0:
            try:
                test_access = gg_archive(gg_file_label.cget("text"))
                if str(test_access) == 'File is unloaded.':
                    error_window("Invalid reanalysis file. Check it and try again.")
                    return
            except:
                error_window("Invalid reanalysis file. Check it and try again.")
                return
        
        # Checks whether the sample files loaded are within gg file
        if len(gg_file_label.cget("text")) > 0:
            gg_files_list = list(test_access.list_samples().values())
            samples_list = [files_list.item(item_id, 'text') for item_id in item_ids]
            samples_names = Execution_Functions.sample_names(samples_list)
            
            for name in samples_names:
                # Get the text of the item (i.e., its content)
                if name not in gg_files_list:
                    error_window("You can't select a reanalysis file and samples files at once, unless the sample files are present in the reanalysis file. Remove the reanalysis file or the sample files or select sample files that are present in the reanalysis file.")
                    test_access.close_gg()
                    return
                    
        # Closes the test access and clears the temporary folder
        if len(gg_file_label.cget("text")) > 0 and str(test_access) != 'File is unloaded.':
            test_access.close_gg()
            
        # Activates and deactivates the proper buttons and fields
        filter_list.delete(0, tk.END)
        filter_list.configure(fg='grey')
        filter_list.insert(0, "Filter the list of glycans...")
        compare_samples_button.config(state=tk.DISABLED)
        if len(gg_file_label.cget("text")) > 0:
            run_analysis_button.config(state=tk.DISABLED)
            check_qc_dist_button.config(state=tk.NORMAL)
            plot_graph_button.config(state=tk.NORMAL)
        else:
            check_qc_dist_button.config(state=tk.DISABLED)
            plot_graph_button.config(state=tk.DISABLED)
                
        # Clear the QC numbers
        chromatograms_qc_numbers.config(text=f"Compositions Quality:\n        Good: {0}    Average: {0}    Bad: {0}\n        Ambiguities: {0}")
            
        # Empty the samples_list and samples_names if nothing is listed
        if len(item_ids) == 0 and len(gg_file_label.cget("text")) == 0:
            samples_list = []
            samples_names = []
            reanalysis_path = ''
        
            # Clear labels
            rt_label.config(text="")
            precursor_label.config(text="")
            coordinate_label_spec.config(text='')
            
        elif len(item_ids) > 0 and len(gg_file_label.cget("text")) == 0:
            samples_list = [files_list.item(item_id, 'text') for item_id in item_ids]
            samples_names = Execution_Functions.sample_names(samples_list)
        
        if reanalysis_path != gg_file_label.cget("text"):
            # Try to unload the formerly loaded gg file
            try:
                gg_file.close_gg()
            except:
                pass
        
            # Now that everything is confirmed good, change reanalysis path.
            reanalysis_path = gg_file_label.cget("text")
            
            # Process the gg file
            if len(reanalysis_path) > 0:
                reanalysis_file = threading.Thread(target=load_reanalysis, args=(reanalysis_path,))
                reanalysis_file.start()
                
        # If nothing is loaded, make sure to empty the samples dropdown
        if len(samples_list) == 0 and len(reanalysis_path) == 0:
            samples_dropdown['values'] = []
            ToolTip(samples_dropdown, "")
            
        quick_check.config(state=tk.DISABLED)
        two_d.config(state=tk.DISABLED)
        spectra_file_editing.config(state=tk.DISABLED)
            
        # Process the sample files and activate the proper buttons
        if len(samples_list) > 0:
            if len(reanalysis_path) == 0:
                samples_dropdown_options = samples_names
                samples_dropdown['values'] = samples_dropdown_options
            processed_samples = threading.Thread(target=pre_process, args=(samples_list,))
            processed_samples.start()
            
            # Only activate the buttons if the sample passes the processing.
            quick_check.config(state=tk.NORMAL)
            two_d.config(state=tk.NORMAL)
            spectra_file_editing.config(state=tk.NORMAL)
            run_analysis_button.config(state=tk.NORMAL)
        
        # Creates the loading window
        loading_files()
                    
        # If passes tests, color the button border green
        if len(reanalysis_path) > 0 or len(samples_list) > 0:
            select_files_frame.config(bg='lightgreen')
        else:
            select_files_frame.config(bg=background_color)
        
        # Clear the former_alignments
        former_alignments = []
        
        # Erase last file loaded history
        former_selected_sample = ''
        
        # Close select file window
        close_sf_window()
        
    def get_gg_parameters():        
        test_access = gg_archive(gg_file_label.cget("text"))
        if str(test_access) == 'File is unloaded.':
            error_window("This GG file is corrupted or incompatible with this version of GlycoGenius.")
            return
        parameters_gg = test_access.get_metadata()
        samples_list = list(test_access.list_samples().values())
        version_gg = test_access.get_version()
        
        if type(parameters_gg) != dict:
            error_window(f"This analysis file was created in a version of GlycoGenius before 1.3.0, and thus can only be opened in that or older versions.")
            return
            
        analysis_info_window = tk.Toplevel()
        analysis_info_window.attributes("-topmost", True)
        analysis_info_window.withdraw()
        analysis_info_window.title("Analysis Information")
        icon = ImageTk.PhotoImage(ico_image)
        analysis_info_window.iconphoto(False, icon)
        analysis_info_window.resizable(False, False)
        analysis_info_window.grab_set()

        information_text = ScrolledText(analysis_info_window, width=55, height=30, wrap=tk.WORD)
        information_text.grid(row=0, column=0, padx = 10, pady = 10, sticky="new")
        
        data_formatted = f"{parameters_gg['analysis time'][4:6]}/{parameters_gg['analysis time'][2:4]}/20{parameters_gg['analysis time'][:2]} - {parameters_gg['analysis time'][7:9]}:{parameters_gg['analysis time'][9:11]}.{parameters_gg['analysis time'][11:13]}"
        information_text.insert(tk.END, f"Date and time of analysis: {data_formatted}\n")
        information_text.insert(tk.END, f"Glycogenius version used: {version_gg}\n")
        information_text.insert(tk.END, "\n")
        information_text.insert(tk.END, "Samples analyzed:\n")
        for i in samples_list:
            information_text.insert(tk.END, f"- {i}\n")
        information_text.insert(tk.END, "\n")
        information_text.insert(tk.END, "Library properties:\n")
        library_type = ''
        additional_info = ''
        if parameters_gg['library settings']['import/export library'][0]: #if imported
            library_type = 'Imported Library'
            additional_info = f" - Library file name: {parameters_gg['library settings']['library path'].split('/')[-1]}\n\n"

            if parameters_gg['library settings']['custom glycans list'][0]: #if custom and imported
                additional_info+= f" - Glycans list/path: {str(parameters_gg['library settings']['custom glycans list'][1])[1:-1]}"
                
            else: #if generated and imported
                additional_info+= f"- Monosaccharides: {str(parameters_gg['library settings']['min/max monosaccharides'])[1:-1]}\n - Hexoses: {str(parameters_gg['library settings']['min/max hexoses'])[1:-1]}\n - Hexosamines: {str(parameters_gg['library settings']['min/max hexosamines'])[1:-1]}\n - HexNAcs: {str(parameters_gg['library settings']['min/max hexnac'])[1:-1]}\n - Xyloses: {str(parameters_gg['library settings']['min/max xyloses'])[1:-1]}\n - Sialic Acids: {str(parameters_gg['library settings']['min/max sialic acids'])[1:-1]}\n - dHex: {str(parameters_gg['library settings']['min/max fucoses'])[1:-1]}\n - Neu5Acs: {str(parameters_gg['library settings']['min/max acetyl sialic acids'])[1:-1]}\n - Neu5Gcs: {str(parameters_gg['library settings']['min/max glycolyl sialic acids'])[1:-1]}\n - Uronic Acids: {str(parameters_gg['library settings']['min/max uronic acids'])[1:-1]}"
            
        elif parameters_gg['library settings']['custom glycans list'][0]: #if custom and not imported
            library_type = 'Custom glycans list'
            additional_info = f" - Glycans list/path: {str(parameters_gg['library settings']['custom glycans list'][1])[1:-1]}"
            
        else: #if generated and not imported
            library_type = 'Generated Library'
            additional_info+= f"- Monosaccharides: {str(parameters_gg['library settings']['min/max monosaccharides'])[1:-1]}\n - Hexoses: {str(parameters_gg['library settings']['min/max hexoses'])[1:-1]}\n - Hexosamines: {str(parameters_gg['library settings']['min/max hexosamines'])[1:-1]}\n - HexNAcs: {str(parameters_gg['library settings']['min/max hexnac'])[1:-1]}\n - Xyloses: {str(parameters_gg['library settings']['min/max xyloses'])[1:-1]}\n - Sialic Acids: {str(parameters_gg['library settings']['min/max sialic acids'])[1:-1]}\n - dHex: {str(parameters_gg['library settings']['min/max fucoses'])[1:-1]}\n - Neu5Acs: {str(parameters_gg['library settings']['min/max acetyl sialic acids'])[1:-1]}\n - Neu5Gcs: {str(parameters_gg['library settings']['min/max glycolyl sialic acids'])[1:-1]}\n - Uronic Acids: {str(parameters_gg['library settings']['min/max uronic acids'])[1:-1]}"
        
        if len(parameters_gg['library settings']['custom monosaccharides']) > 0:
            additional_info+= f"\n\n - Custom monosaccharides:\n"
            for cm in parameters_gg['library settings']['custom monosaccharides']:
                additional_info+= f"    - {cm['cm_name']} ({cm['cm_short_code']}): {cm['cm_chem_comp']} "
                if not parameters_gg['library settings']['custom glycans list'][0]:
                    additional_info+= f"Min/Max: {cm['cm_min']}/{cm['cm_max']}"
            
        if parameters_gg['library settings']['reducing end tag'] in reducing_end_tags.values():
            reducing_end_tag_text = f" ({list(reducing_end_tags.keys())[list(reducing_end_tags.values()).index(parameters_gg['library settings']['reducing end tag'])]})"
        else:
            reducing_end_tag_text = ''
        
        adduct_combos_dict = {General_Functions.comp_to_formula(adduct[0]): {'charges': adduct[1], 'comp': adduct[0]} for adduct in parameters_gg['adduct combos']}
        adduct_combos_string = ''
        for adduct, info in adduct_combos_dict.items():
            if len(adduct_combos_string) == 0:
                adduct_combos_string += f"{adduct}: {'+' if info['charges'] > 0 else ''}{info['charges']}"
            else:
                adduct_combos_string += f"\n            {adduct}: {'+' if info['charges'] > 0 else ''}{info['charges']}"
                
        additional_info+= f"\n\n - Forced glycan class: {parameters_gg['library settings']['glycan class']}\n - Adducts: {adduct_combos_string}\n - Maximum charges: {parameters_gg['library settings']['maximum charges']}\n - Reducing end tag: {parameters_gg['library settings']['reducing end tag'] if parameters_gg['library settings']['reducing end tag'] != 0.0 else False}{reducing_end_tag_text}\n - Permethylated: {parameters_gg['library settings']['permethylated']}\n - Reduced end: {parameters_gg['library settings']['reducing end reduced']}\n - Amidated/Ethyl-Esterified: {parameters_gg['library settings']['lactonized/ethyl esterified']}"
        
        additional_info+= f"\n - Min/Max number of Sulfations: {str(parameters_gg['library settings']['min/max sulfation'])[1:-1]}"
        additional_info+= f"\n - Min/Max number of Phosphorylations: {str(parameters_gg['library settings']['min/max phosphorylation'])[1:-1]}"
            
        additional_info+= f"\n - Fast isotopic calculations: {parameters_gg['library settings']['fast isotopic pattern calculation']}\n - High resolution isotopic calculations: {parameters_gg['library settings']['high resolution isotopic pattern']}\n - Internal standard: {parameters_gg['library settings']['internal standard'] if parameters_gg['library settings']['internal standard'] != 0.0 else False}"
        information_text.insert(tk.END, f" - Library type: {library_type}\n")
        information_text.insert(tk.END, "\n")
        information_text.insert(tk.END, f"{additional_info}\n")
        information_text.insert(tk.END, "\n")
        information_text.insert(tk.END, "Analysis settings:\n")
        information_text.insert(tk.END, f" - Analyze MS2: {parameters_gg['analysis settings']['analyze ms2'][0]}\n")
        if parameters_gg['analysis settings']['analyze ms2'][0]:
            information_text.insert(tk.END, f" -- Limit fragments assignment to composition: {parameters_gg['analysis settings']['analyze ms2'][1]}\n")
            information_text.insert(tk.END, f" -- Assign MS2 of glycans not found in MS1: {parameters_gg['analysis settings']['analyze ms2'][2]}\n")
        information_text.insert(tk.END, f" - Tolerance unit: {parameters_gg['analysis settings']['mass tolerance'][0]}, tolerance value: {parameters_gg['analysis settings']['mass tolerance'][1]}\n")
        information_text.insert(tk.END, f" - Retention/Migration time interval analyzed: {parameters_gg['analysis settings']['retention time interval'][0]}, {parameters_gg['analysis settings']['retention time interval'][1]}\n")
        information_text.insert(tk.END, f" - Custom minimum datapoints per peaks: {parameters_gg['analysis settings']['min points per peak'][0]}\n")
        if parameters_gg['analysis settings']['min points per peak'][0]:
            information_text.insert(tk.END, f" -- Minimum number of datapoints per peaks: {parameters_gg['analysis settings']['min points per peak'][1]}\n")
        information_text.insert(tk.END, f" - Limit peaks picked per chromatogram/electropherogram: {parameters_gg['analysis settings']['only x most intense peaks'][0]}\n")
        if parameters_gg['analysis settings']['only x most intense peaks'][0]:
            information_text.insert(tk.END, f" -- Number of peaks: {parameters_gg['analysis settings']['only x most intense peaks'][1]}\n")
            
        information_text.see("end")
        
        close_analysis_info_button = ttk.Button(analysis_info_window, text="Close", style="small_button_sfw_style1.TButton", command=analysis_info_window.destroy)
        close_analysis_info_button.grid(row=1, column=0, padx=10, pady=10)
        
        analysis_info_window.update_idletasks()
        analysis_info_window.deiconify()
        window_width = analysis_info_window.winfo_width()
        window_height = analysis_info_window.winfo_height()
        screen_width = analysis_info_window.winfo_screenwidth()
        screen_height = analysis_info_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        analysis_info_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        test_access.close_gg()
        
    def prepare_sample_groups_window():
        if len(gg_file_label.cget("text")) == 0:
            item_ids = files_list.get_children()
            samples_list = [files_list.item(item_id, 'text') for item_id in item_ids]
            samples_names = Execution_Functions.sample_names(samples_list)
        else:
            test_access = gg_archive(gg_file_label.cget("text"))
            if str(test_access) == 'File is unloaded.':
                error_window("This GG file is corrupted or incompatible with this version of GlycoGenius.")
                return
            samples_names = list(test_access.list_samples().values())
            
        sample_grouping_window(samples_names, select_files_window)
        
    # Create a new top-level window
    select_files_window = tk.Toplevel()
    #select_files_window.attributes("-topmost", True)
    select_files_window.withdraw()
    select_files_window.title("Select Files")
    icon = ImageTk.PhotoImage(ico_image)
    select_files_window.iconphoto(False, icon)
    select_files_window.resizable(False, False)
    select_files_window.grab_set()
    
    # Add Ttk Styles
    files_list_style = ttk.Style().configure("files_list.Treeview", bg="white", fg="black", font=("Segoe UI", list_font_size))
    
    small_button_sfw_style1 = ttk.Style().configure("small_button_sfw_style1.TButton", font=("Segoe UI", list_font_size), relief="raised", padding = (0, 0), justify="center")
    
    small_button_sfw_style2 = ttk.Style().configure("small_button_sfw_style2.TButton", font=("Segoe UI", list_font_size), relief="raised", padding = (10, 0), justify="center")

    # Add widgets to the new window
    main_label = ttk.Label(select_files_window, text="Select files for analysis:", font=("Segoe UI", list_font_size))
    main_label.grid(row=0, column=0, padx=(10, 10), pady=(10, 10), sticky="nw")
    
    files_list_frame = tk.Frame(select_files_window, bd=0, relief="flat", width=450, height=400)
    files_list_frame.pack_propagate(0)
    files_list_frame.grid(row=0, column=0, padx=(10,20), pady=(30, 10), sticky="nsew")
    
    files_list_scrollbar = tk.Scrollbar(select_files_window, orient=tk.HORIZONTAL)
    files_list_scrollbar_v = tk.Scrollbar(select_files_window, orient=tk.VERTICAL)
    files_list = ttk.Treeview(files_list_frame, height=20, style="files_list.Treeview", show="tree", xscrollcommand=files_list_scrollbar.set, yscrollcommand=files_list_scrollbar_v.set)
    files_list.column("#0", stretch=False, width=448)
    files_list_scrollbar.config(command=files_list.xview, width=10)
    files_list_scrollbar_v.config(command=files_list.yview, width=10)
    files_list.pack(fill="both", expand=True, padx=0, pady=0)
    files_list_scrollbar.grid(row=0, column=0, padx=(10, 20), pady=0, sticky="sew")
    files_list_scrollbar_v.grid(row=0, column=0, padx=(0, 10), pady=(30, 10), sticky="nse")
        # Load the samples_list from global variable
    global samples_list
    for item in samples_list:
        files_list.insert('', 'end', text=item)
    longest_len = 0
    for i in files_list.get_children():
        if len(files_list.item(i, "text")) > longest_len:
            longest_len = len(files_list.item(i, "text"))
    if longest_len*7 > 448:
        files_list.column("#0", width=longest_len*7)
    
    add_file_button = ttk.Button(select_files_window, text="Add File", style="small_button_sfw_style1.TButton", command=open_file_dialog_sfw)
    add_file_button.grid(row=0, column=1, padx=(0,10), pady=(30,10), sticky="new")
    ToolTip(add_file_button, "Select one or more MzML or MzXML files to add to the list of files to open and analyze with GlycoGenius.")
    
    remove_selected_button = ttk.Button(select_files_window, text="Remove Selected", style="small_button_sfw_style1.TButton", command=remove_selected_item_sfw)
    remove_selected_button.grid(row=0, column=1, padx=(0,10), pady=(60,10), sticky="new")
    ToolTip(remove_selected_button, "Remove the currently selected file from the list of files to load.")
    
    arrange_groups_button = ttk.Button(select_files_window, text="Arrange Groups", style="small_button_sfw_style1.TButton", command=prepare_sample_groups_window)
    arrange_groups_button.grid(row=0, column=1, padx=(0,10), pady=(90,10), sticky="new")
    ToolTip(arrange_groups_button, "Allows you to separate the samples into groups through an easy-to-use interface. Samples left ungrouped will be excluded from the analysis, unless all samples are ungrouped. Relevant for results saving settings.")
    
    file_amount_label = ttk.Label(select_files_window, text="", font=("Segoe UI", list_font_size), wraplength=310)
    file_amount_label.grid(row=0, column=1, padx=(0,10), pady=(120,10), sticky="new")
    if len(files_list.get_children()) > 0:
        if len(files_list.get_children()) == 1:
            file_amount_label.config(text=f'{len(files_list.get_children())} file selected')
        else:
            file_amount_label.config(text=f'{len(files_list.get_children())} files selected')
    else:
        file_amount_label.config(text='')
    ToolTip(file_amount_label, "The number of files loaded.")
    
    reanalyze_gg_button = ttk.Button(select_files_window, text="    Load .gg file   ", style="small_button_sfw_style2.TButton", command=open_file_dialog_sfw_reanalysis)
    reanalyze_gg_button.grid(row=1, column=0, padx=(10,10), pady=(10,40), sticky="wns")
    ToolTip(reanalyze_gg_button, "Select a GlycoGenius analysis file (.gg) to check the results and save new files of them.")
    
    gg_file_label = ttk.Label(select_files_window, text="", font=("Segoe UI", list_font_size), wraplength=310)
    gg_file_label.grid(row=1, column=0, columnspan=2, padx=(150, 10), pady=(5, 10), sticky="w")
    global reanalysis_path
    gg_file_label.config(text=reanalysis_path)
    ToolTip(gg_file_label, "The currently selected .gg file.")
    
    file_info_button_state = tk.DISABLED
    if len(gg_file_label.cget("text")) > 0:
        file_info_button_state = tk.NORMAL
    file_info_button = ttk.Button(select_files_window, text="    File Information    ", style="small_button_sfw_style1.TButton", command=get_gg_parameters, state=file_info_button_state)
    file_info_button.grid(row=1, column=0, columnspan=2, padx=(10, 10), pady=(40,10), sticky="nsw")
    ToolTip(file_info_button, "Check the information of how the analysis on the .gg file was performed.")
    
    ok_button = ttk.Button(select_files_window, text="Ok", style="small_button_sfw_style1.TButton", command=ok_button_sfw)
    ok_button.grid(row=1, column=1, padx=(0,100), pady=(40,15), sticky="nse")
    
    cancel_button = ttk.Button(select_files_window, text="Cancel", style="small_button_sfw_style1.TButton", command=close_sf_window)
    cancel_button.grid(row=1, column=1, padx=(0,10), pady=(40,15), sticky="nse")
    
    select_files_window.update_idletasks()
    select_files_window.deiconify()
    window_width = select_files_window.winfo_width()
    window_height = select_files_window.winfo_height()
    screen_width = select_files_window.winfo_screenwidth()
    screen_height = select_files_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    select_files_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
def run_set_parameters_window():
    '''
    '''
    def load_library_data():
        try:
            with open(library_path, 'rb') as f:
                library_data = dill.load(f)
                f.close()
        except:
            return
            
        full_library = library_data[0]
        library_metadata = library_data[1]
        
        if type(library_metadata) != dict:
            error_window("Can't load the library. The library file you are trying to load was created in a version older than 1.3.0 or is corrupted. Generate a new library or try a different file.")
            return "Error"
        
        global min_max_monos, min_max_hex, min_max_hexnac, min_max_xyl, min_max_fuc, min_max_sia, min_max_ac, min_max_gc, min_max_hn, min_max_ua, forced, min_max_proton, custom_adducts, max_charges, reducing_end_tag, internal_standard, permethylated, lactonized_ethyl_esterified, min_max_sulfation, min_max_phosphorylation, reduced, fast_iso, high_res, lyase_digested, custom_monosaccharides, custom_glycans_list, ion_mode_selected
        
        if library_metadata.get('custom glycans list', [False, []])[0]:
            is_custom = True
            custom_glycans_list[1] = library_metadata.get('custom glycans list', [False, []])[1]
        min_max_monos = library_metadata.get('min/max monosaccharides', [0, 0])
        min_max_hex = library_metadata.get('min/max hexoses', [0, 0])
        min_max_hexnac = library_metadata.get('min/max hexnac', [0, 0])
        min_max_fuc = library_metadata.get('min/max fucoses', [0, 0])
        min_max_sia = library_metadata.get('min/max sialic acids', [0, 0])
        min_max_ac = library_metadata.get('min/max acetyl sialic acids', [0, 0])
        min_max_gc = library_metadata.get('min/max glycolyl sialic acids', [0, 0])
        min_max_xyl = library_metadata.get('min/max xyloses', [0, 0])
        min_max_hn = library_metadata.get('min/max hexosamines', [0, 0])
        min_max_ua = library_metadata.get('min/max uronic acids', [0, 0])
        min_max_sulfation = library_metadata.get('min/max sulfation', [0, 0])
        min_max_phosphorylation = library_metadata.get('min/max phosphorylation', [0, 0])
        forced = library_metadata.get('glycan class', None)
        min_max_proton = library_metadata.get('min/max protons', [0, 0])
        custom_adducts = library_metadata.get('custom adducts', [])
        max_charges = library_metadata.get('maximum charges', 3)
        reducing_end_tag = library_metadata.get('reducing end tag', 0.0)
        internal_standard = library_metadata.get('internal standard', '0.0')
        permethylated = library_metadata.get('permethylated', False)
        lactonized_ethyl_esterified = library_metadata.get('lactonized/ethyl esterified', False)
        reduced = library_metadata.get('reducing end reduced', False)
        fast_iso = library_metadata.get('fast isotopic pattern calculation', True)
        high_res = library_metadata.get('high resolution isotopic pattern', False)
        lyase_digested = library_metadata.get('lyase digested', False)
        custom_monosaccharides = library_metadata.get('custom monosaccharides', [])
            
        ion_mode_selected = True
            
    # Fetch global variables, whenever possible
    global custom_glycans_list, analyze_ms2, custom_monosaccharides
    
    # If library path is set, load the library parameters
    if len(library_path) > 0:
        lib_check = load_library_data()
        if lib_check == 'Error':
            return
    
    local_analyze_ms2 = analyze_ms2
    local_custom_glycans_list = copy.deepcopy(custom_glycans_list)
    local_custom_monos_list = copy.deepcopy(custom_monosaccharides)
    local_custom_adducts_list = copy.deepcopy(custom_adducts)
    
    # Functions used by this window
    def close_sp_window():
        set_parameters_window.destroy()  
        
    def select_working_dir_button():
        file_dialog = tk.Toplevel()
        file_dialog.attributes("-topmost", True)
        file_dialog.withdraw()
        file_dialog.grab_set()
        
        dir_path = filedialog.askdirectory()
        
        if dir_path != '':
            if dir_path[-1] != "/":
                dir_path += "/"
            
        working_dir_label.config(text=dir_path)
        
        file_dialog.destroy()
        set_parameters_window.grab_set()
        
    def use_custom_library_checkbox_state_check():
        state = use_custom_library_checkbox_state.get()
        if int(state) == 1:
            local_custom_glycans_list[0] = True
            from_file_button.config(state=tk.NORMAL)
            min_monosaccharides_entry.config(state=tk.DISABLED)
            max_monosaccharides_entry.config(state=tk.DISABLED)
            min_hex_entry.config(state=tk.DISABLED)
            max_hex_entry.config(state=tk.DISABLED)
            min_hn_entry.config(state=tk.DISABLED)
            max_hn_entry.config(state=tk.DISABLED)
            min_hexnac_entry.config(state=tk.DISABLED)
            max_hexnac_entry.config(state=tk.DISABLED)
            min_xyl_entry.config(state=tk.DISABLED)
            max_xyl_entry.config(state=tk.DISABLED)
            min_dhex_entry.config(state=tk.DISABLED)
            max_dhex_entry.config(state=tk.DISABLED)
            min_sia_entry.config(state=tk.DISABLED)
            max_sia_entry.config(state=tk.DISABLED)
            min_ua_entry.config(state=tk.DISABLED)
            max_ua_entry.config(state=tk.DISABLED)
            min_neu5ac_entry.config(state=tk.DISABLED)
            max_neu5ac_entry.config(state=tk.DISABLED)
            min_neu5gc_entry.config(state=tk.DISABLED)
            max_neu5gc_entry.config(state=tk.DISABLED)
        else:
            local_custom_glycans_list[0] = False
            from_file_button.config(state=tk.DISABLED)
            min_monosaccharides_entry.config(state=tk.NORMAL)
            max_monosaccharides_entry.config(state=tk.NORMAL)
            min_hex_entry.config(state=tk.NORMAL)
            max_hex_entry.config(state=tk.NORMAL)
            min_hn_entry.config(state=tk.NORMAL)
            max_hn_entry.config(state=tk.NORMAL)
            min_hexnac_entry.config(state=tk.NORMAL)
            max_hexnac_entry.config(state=tk.NORMAL)
            min_xyl_entry.config(state=tk.NORMAL)
            max_xyl_entry.config(state=tk.NORMAL)
            min_dhex_entry.config(state=tk.NORMAL)
            max_dhex_entry.config(state=tk.NORMAL)
            min_sia_entry.config(state=tk.NORMAL)
            max_sia_entry.config(state=tk.NORMAL)
            min_ua_entry.config(state=tk.NORMAL)
            max_ua_entry.config(state=tk.NORMAL)
            min_neu5ac_entry.config(state=tk.NORMAL)
            max_neu5ac_entry.config(state=tk.NORMAL)
            min_neu5gc_entry.config(state=tk.NORMAL)
            max_neu5gc_entry.config(state=tk.NORMAL)
            
    def reducing_end_tag_checkbox_state_check():
        state = reducing_end_tag_checkbox_state.get()
        if int(state) == 1:
            reducing_end_tag_dropdown.state(['!disabled'])
            reduced_checkbox.config(state=tk.DISABLED)
            if reducing_end_tag_dropdown.get() == 'Custom':
                reducing_end_tag_entry.config(state=tk.NORMAL)
            else:
                reducing_end_tag_entry.config(state=tk.DISABLED)
        else:
            reducing_end_tag_entry.delete(0, tk.END)
            reducing_end_tag_entry.insert(0, '0.0')
            reducing_end_tag_entry.config(state=tk.DISABLED) 
            reducing_end_tag_dropdown.set('Custom')
            reducing_end_tag_dropdown.state(['disabled'])
            reduced_checkbox.config(state=tk.NORMAL)
            
    def on_tag_selected(event):
        selected_tag = reducing_end_tag_dropdown.get()
        if selected_tag == 'Custom' or selected_tag == 'Peptide':
            reducing_end_tag_entry.config(state=tk.NORMAL)
        else:
            reducing_end_tag_entry.config(state=tk.NORMAL)
            reducing_end_tag_entry.delete(0, tk.END)
            reducing_end_tag_entry.insert(0, reducing_end_tags[selected_tag])
            reducing_end_tag_entry.config(state=tk.DISABLED)
            
    def permethylated_checkbox_state_check():
        state = permethylated_checkbox_state.get()
        
    def reduced_checkbox_state_check():
        state = reduced_checkbox_state.get()
        if int(state) == 1:
            reducing_end_tag_entry.config(state=tk.DISABLED)   
            reducing_end_tag_checkbox.config(state=tk.DISABLED)
        else:
            if reducing_end_tag_checkbox_state.get():
                reducing_end_tag_entry.config(state=tk.NORMAL) 
            reducing_end_tag_checkbox.config(state=tk.NORMAL)
        
    def lac_ee_checkbox_state_check():
        state = lac_ee_checkbox_state.get()
        
    def fast_iso_checkbox_state_check():
        state = fast_iso_checkbox_state.get()
        if int(state) == 1:
            hires_iso_checkbox.config(state=tk.DISABLED)
        else:
            hires_iso_checkbox.config(state=tk.NORMAL)  
        
    def hires_iso_checkbox_state_check():
        state = hires_iso_checkbox_state.get()
        
    def multithreaded_checkbox_state_check():
        state = multithreaded_checkbox_state.get()
        if int(state) == 1:
            number_cores_entry.config(state=tk.NORMAL)   
        else:
            number_cores_entry.config(state=tk.DISABLED)
            
    def analyze_ms2_checkbox_state_check():
        state = analyze_ms2_checkbox_state.get()
        if int(state) == 1:
            local_analyze_ms2[0] = True
            force_ms2_comp_checkbox.config(state=tk.NORMAL)
            unrestricted_frags_checkbox.config(state=tk.NORMAL)
        else:
            local_analyze_ms2[0] = False
            force_ms2_comp_checkbox.config(state=tk.DISABLED)
            unrestricted_frags_checkbox.config(state=tk.DISABLED)
            
    def force_ms2_comp_checkbox_state_check():
        state = force_ms2_comp_checkbox_state.get()
        if int(state) == 1:
            local_analyze_ms2[1] = True
        else:
            local_analyze_ms2[1] = False
        
    def unrestricted_frags_checkbox_state_check():
        state = unrestricted_frags_checkbox_state.get()
        if int(state) == 1:
            local_analyze_ms2[2] = True
        else:
            local_analyze_ms2[2] = False
                    
    def custom_ppp_checkbox_state_check():
        state = custom_ppp_checkbox_state.get()
        if int(state) == 1:
            custom_ppp_entry.config(state=tk.NORMAL)
        else:
            custom_ppp_entry.config(state=tk.DISABLED)
            
    def close_peaks_checkbox_state_check():
        state = close_peaks_checkbox_state.get()
        if int(state) == 1:
            close_peaks_entry.config(state=tk.NORMAL)
        else:
            close_peaks_entry.config(state=tk.DISABLED)
        
    def lyase_digested_checkbox_state_check():
        state = lyase_digested_checkbox_state.get()
            
    def on_forced_class_selected(event):
        global gallery_warning_label
        
        forced_class = forced_class_dropdown.get()
        
        if forced_class == 'GAGs':
            lyase_digested_checkbox.config(state=tk.NORMAL)
        else:
            lyase_digested_checkbox.config(state=tk.DISABLED)
            
    def ok_sp_window():
        global save_path, min_max_monos, min_max_hex, min_max_hexnac, min_max_xyl, min_max_fuc,  min_max_sia, min_max_ac, min_max_gc, min_max_hn, min_max_ua, internal_standard, reducing_end_boolean, reducing_end_tag, permethylated, reduced, lactonized_ethyl_esterified,  min_max_sulfation,  min_max_phosphorylation,  forced, fast_iso, high_res, multithreaded_analysis, number_cores, min_ppp,  close_peaks, iso_fit_score, curve_fit_score, s_to_n, max_ppm, max_charges, min_max_proton, custom_adducts, custom_glycans_list, lyase_digested, custom_monosaccharides, ion_mode_selected
        
        parameters_dict = {'Min. Monosaccharides':min_monosaccharides_entry.get(), 
                           'Max. Monosaccharides':max_monosaccharides_entry.get(), 
                           'Min. Hexoses':min_hex_entry.get(), 
                           'Max. Hexoses':max_hex_entry.get(), 
                           'Min. Hexosamines':min_hn_entry.get(), 
                           'Max. Hexosamines':max_hn_entry.get(), 
                           'Min. HexNAcs':min_hexnac_entry.get(), 
                           'Max. HexNAcs':max_hexnac_entry.get(),
                           'Min. Fucoses':min_dhex_entry.get(),
                           'Max. Fucoses':max_dhex_entry.get(),
                           'Min. Sialic Acids':min_sia_entry.get(),
                           'Max. Sialic Acids':max_sia_entry.get(),
                           'Min. Uronic Acids':min_ua_entry.get(),
                           'Max. Uronic Acids':max_ua_entry.get(),
                           'Min. Neu5Ac':min_neu5ac_entry.get(),
                           'Max. Neu5Ac':max_neu5ac_entry.get(),
                           'Min. Neu5Gc':min_neu5gc_entry.get(),
                           'Max. Neu5Gc':max_neu5gc_entry.get(), 
                           'Internal Standard Mass':intstandard_entry.get(), 
                           'Reducing End Tag':reducing_end_tag_entry.get(), 
                           'Number of Cores':number_cores_entry.get(), 
                           'Min. Datapoints per Peak':custom_ppp_entry.get(), 
                           'Max. Number of Peaks Picked':close_peaks_entry.get(),
                           'Min. H Adduct':hydrogen_min_entry.get(), 
                           'Max. H Adduct':hydrogen_max_entry.get(), 
                           'Max. Charges':max_charges_entry.get(),
                           'Min. Sulfation':min_sulfation_entry.get(),
                           'Max. Sulfation':max_sulfation_entry.get(),
                           'Min. Phosphorylation':min_phosphorylation_entry.get(),
                           'Max. Phosphorylation':max_phosphorylation_entry.get(),
                           'Min. Retention/Migration Time':rt_int_min_entry.get(),
                           'Max. Retention/Migration Time':rt_int_max_entry.get(),
                           'Accuracy Value':acc_value_entry.get()}
                           
        for adduct, widgets in custom_adducts_widgets.items():
            parameters_dict[f'Min. {adduct} Adduct'] = widgets[0].get()
            parameters_dict[f'Max. {adduct} Adduct'] = widgets[1].get()
            
        cores_good = True
        temp_number_cores = number_cores_entry.get()
        try:
            temp_number_cores = int(number_cores_entry.get())
        except:
            if temp_number_cores.lower() == 'all':
                pass
            else:
                cores_good=False
                error_window("Invalid input on the field:\n\nNumber of Cores\n\nCorrect it and try again.")
        if not ion_mode_choice.get() and len(samples_list) > 0:
            error_window("You must set an ion mode before confirming settings!")
            return
        if cores_good:
            try:
                if int(hydrogen_max_entry.get()) < int(hydrogen_min_entry.get()):
                    error_window("Max. H value must be greater than min. H value.")
                    return
                min_max_proton = [int(hydrogen_min_entry.get()), int(hydrogen_max_entry.get())]
                for adduct, widgets in custom_adducts_widgets.items():
                    if int(widgets[0].get()) > int(widgets[1].get()):
                        error_window(f"Max. custom adduct {adduct} value must be greater than the min.")
                        return
                for index, adduct_list in enumerate(local_custom_adducts_list):
                    local_custom_adducts_list[index][1] = custom_adducts_widgets[adduct_list[0]][0].get()
                    local_custom_adducts_list[index][2] = custom_adducts_widgets[adduct_list[0]][1].get()
                custom_adducts = local_custom_adducts_list
                
                if float(rt_int_max_entry.get()) < float(rt_int_min_entry.get()):
                    error_window("Max. Retention/Migration Time value must be greater than min. Retention/Migration Time value.")
                    return
                ret_time_interval[0] = float(rt_int_min_entry.get())
                ret_time_interval[1] = float(rt_int_max_entry.get())
                tolerance[1] = float(acc_value_entry.get())
                close_peaks = [close_peaks_checkbox_state.get(), int(close_peaks_entry.get())]
                
                if int(max_monosaccharides_entry.get()) < int(min_monosaccharides_entry.get()):
                    error_window("Max. Monosaccharides value must be greater than min. Monosaccharides value.")
                    return
                min_max_monos = [int(min_monosaccharides_entry.get()), int(max_monosaccharides_entry.get())]
                
                if int(max_hex_entry.get()) < int(min_hex_entry.get()):
                    error_window("Max. Hexoses value must be greater than min. Hexoses value.")
                    return
                min_max_hex = [int(min_hex_entry.get()), int(max_hex_entry.get())]
                
                if int(max_hexnac_entry.get()) < int(min_hexnac_entry.get()):
                    error_window("Max. HexNAc value must be greater than min. HexNAc value.")
                    return
                min_max_hexnac = [int(min_hexnac_entry.get()), int(max_hexnac_entry.get())]
                
                if int(max_xyl_entry.get()) < int(min_xyl_entry.get()):
                    error_window("Max. Xyl value must be greater than min. Xyl value.")
                    return
                min_max_xyl = [int(min_xyl_entry.get()), int(max_xyl_entry.get())]
                
                if int(max_dhex_entry.get()) < int(min_dhex_entry.get()):
                    error_window("Max. dHex value must be greater than min. dHex value.")
                    return
                min_max_fuc = [int(min_dhex_entry.get()), int(max_dhex_entry.get())]
                
                if int(max_sia_entry.get()) < int(min_sia_entry.get()):
                    error_window("Max. Sialic Acids value must be greater than min. Sialic Acids value.")
                    return
                min_max_sia = [int(min_sia_entry.get()), int(max_sia_entry.get())]
                
                if int(max_neu5ac_entry.get()) < int(min_neu5ac_entry.get()):
                    error_window("Max. Neu5Ac value must be greater than min. Neu5Ac value.")
                    return
                min_max_ac = [int(min_neu5ac_entry.get()), int(max_neu5ac_entry.get())]
                
                if int(max_neu5gc_entry.get()) < int(min_neu5gc_entry.get()):
                    error_window("Max. Neu5Gc value must be greater than min. Neu5Gc value.")
                    return
                min_max_gc = [int(min_neu5gc_entry.get()), int(max_neu5gc_entry.get())]
                
                if int(max_hn_entry.get()) < int(min_hn_entry.get()):
                    error_window("Max. Hexosamines value must be greater than min. Hexosamines value.")
                    return
                min_max_hn = [int(min_hn_entry.get()), int(max_hn_entry.get())]
                
                if int(max_ua_entry.get()) < int(min_ua_entry.get()):
                    error_window("Max. Uronic Acids value must be greater than min. Uronic Acids value.")
                    return
                min_max_ua = [int(min_ua_entry.get()), int(max_ua_entry.get())]
                
                # Check the internal standard
                if intstandard_entry.get().split("-")[0] == 'glycan':
                    test_is = General_Functions.form_to_comp_glycans(intstandard_entry.get().split("-")[1])
                    for i in test_is.keys():
                        if i not in list(General_Functions.monosaccharides.keys())+[k['cm_short_code'] for k in custom_monosaccharides]:
                            error_window("Invalid input on Internal Standard.\nYou must input either a mass, a chemical formula or a glycan formula.")
                            return
                        
                        elif i in ['Am', 'E', 'AmG', 'EG'] and not lac_ee_checkbox_state.get():
                            error_window("You must toggle 'Amidated/Ethyl-Esterified' checkbox in order to use Amidated or Ethyl-Esterified sialic acids in Internal Standard composition.")
                            return
                            
                        elif (i == 'S' or i == 'G') and lac_ee_checkbox_state.get():
                            error_window("You can't have a non-derivatized sialic acid in your internal standard when you have the 'Amidated/Ethyl-Esterified' checkbox toggled.")
                            return
                
                else:
                    try:
                        float(intstandard_entry.get())
                    except:
                        test_is = General_Functions.form_to_comp_atoms(intstandard_entry.get())
                        for i in test_is.keys():
                            if i not in ['H', 'C', 'O', 'N', 'S', 'P', 'Cl', 'Na', 'K']:
                                error_window("Invalid input on Internal Standard.\nYou must input either a mass, a chemical formula or a glycan formula.")
                                return
                        
                internal_standard = intstandard_entry.get()
                
                min_ppp = [custom_ppp_checkbox_state.get(), int(custom_ppp_entry.get())]
                if acc_unit_dropdown.get() == 'PPM':
                    tolerance[0] = 'ppm'
                elif acc_unit_dropdown.get() == 'mz':
                    tolerance[0] = 'mz'
                if len(working_dir_label.cget('text')) != 0:
                    save_path = working_dir_label.cget('text')
                        
                    try:
                        pathlib.Path(save_path).mkdir(exist_ok = True, parents = True)
                    except Exception:
                        error_window("There was an issue while creating the working directory.")
                        return
                else:
                    save_path = ""
                custom_glycans_list = local_custom_glycans_list
                custom_monosaccharides = local_custom_monos_list
                analyze_ms2 = local_analyze_ms2
                
                # Test reducing end tag
                reducing_end_boolean = reducing_end_tag_checkbox_state.get()
                if reducing_end_tag_dropdown.get() == 'Peptide' or reducing_end_tag_entry.get().split("-")[0] == 'pep':
                    if len(reducing_end_tag_entry.get().split("-")) > 1:
                        to_check = reducing_end_tag_entry.get().split("-")[1]
                    else:
                        to_check = reducing_end_tag_entry.get()
                    if any(ch.isdigit() for ch in to_check):
                        error_window("A peptide sequence must be only built of letters.")
                        return
                reducing_end_tag = f"{f'pep-{reducing_end_tag_entry.get()}' if reducing_end_tag_dropdown.get() == 'Peptide' else reducing_end_tag_entry.get()}"
                permethylated = permethylated_checkbox_state.get()
                reduced = reduced_checkbox_state.get()
                lactonized_ethyl_esterified = lac_ee_checkbox_state.get()
                
                if int(max_sulfation_entry.get()) < int(min_sulfation_entry.get()):
                    error_window("Max. Sulfation value must be greater than min. Sulfation value.")
                    return
                min_max_sulfation = [int(min_sulfation_entry.get()), int(max_sulfation_entry.get())]
                
                if int(max_phosphorylation_entry.get()) < int(min_phosphorylation_entry.get()):
                    error_window("Max. Phosphorylation value must be greater than min. Phosphorylation value.")
                    return
                min_max_phosphorylation = [int(min_phosphorylation_entry.get()), int(max_phosphorylation_entry.get())]
                
                forced = forced_classes[forced_class_dropdown.get()]
                lyase_digested = lyase_digested_checkbox_state.get()
                fast_iso = fast_iso_checkbox_state.get()
                high_res = hires_iso_checkbox_state.get()
                multithreaded_analysis = multithreaded_checkbox_state.get()
                number_cores = number_cores_entry.get()
                if len(save_path) > 0 and ion_mode_choice.get():
                    set_parameters_frame.config(bg="lightgreen")
                else:
                    set_parameters_frame.config(bg="red")
                if ion_mode_choice.get():
                    ion_mode_selected = True
                    max_charges = -int(max_charges_entry.get()) if ion_mode_choice.get() == 'neg' else int(max_charges_entry.get())
                close_sp_window()
                
            except Exception as e:
                print(e)
                last_line = traceback.format_exc().split(" ")[-1][1:-2]
                for i_i, i in enumerate(list(parameters_dict.values())):
                    if i == last_line:
                        error_window("Invalid input on the field:\n\n"+list(parameters_dict.keys())[i_i]+"\n\nCorrect it and try again.")
                        break
    
    def save_param_to_file():
        global s_n_entry, curve_fit_entry, ppm_error_min_entry, ppm_error_max_entry, iso_fit_entry
        file_dialog = tk.Toplevel()
        file_dialog.attributes("-topmost", True)
        file_dialog.withdraw()
        file_dialog.grab_set()
        
        file_path = filedialog.asksaveasfilename(defaultextension=".ini", filetypes=[("INI files", "*.ini"), ("All files", "*.*")])
        if len(file_path) == 0:
            file_dialog.destroy()
            set_parameters_window.grab_set()
            return
        if file_path.split(".")[-1] != "ini":
            file_path = file_path.split(".")[0]+".ini"
        with open(file_path, "w") as f:
            f.write(f"[running_modes]\n")
            f.write(f"mode = analysis\n")
            f.write(f"use_multiple_CPU_cores = {multithreaded_checkbox_state.get()}\n")
            f.write(f"number_cores = {number_cores_entry.get()}\n")
            
            export_save_path = working_dir_label.cget('text')
            if len(export_save_path) != 0:
                if export_save_path[-1] != "/":
                    export_save_path = export_save_path+"/"
            if len(export_save_path) != 0:
                f.write(f"working_directory = {export_save_path}\n")
                f.write(f"samples_directory = {export_save_path}\n")
            else:
                f.write(f"working_directory = C:/Glycogenius\n")
                f.write(f"samples_directory = C:/Glycogenius\n")
            f.write(f"exported_library_name = \n")
            f.write(f"file_for_reanalysis = {reanalysis_path}\n")
            f.write("\n")
            f.write(f"[library_building_modes]\n")
            if not use_custom_library_checkbox_state.get():
                f.write(f"mode = generate_library\n")
                f.write(f"export_library = False\n")
                f.write(f"exported_library_name = \n")
                f.write(f"import_library_path = \n")
                f.write(f"custom_glycans_list = \n")
                f.write(f"total_monosaccharides = {int(min_monosaccharides_entry.get())}, {int(max_monosaccharides_entry.get())}\n")
                f.write(f"hexoses = {int(min_hex_entry.get())}, {int(max_hex_entry.get())}\n")
                f.write(f"hexosamines = {int(min_hn_entry.get())}, {int(max_hn_entry.get())}\n")
                f.write(f"hexnacs = {int(min_hexnac_entry.get())}, {int(max_hexnac_entry.get())}\n")
                f.write(f"xyloses = {int(min_xyl_entry.get())}, {int(max_xyl_entry.get())}\n")
                f.write(f"sialic_acids = {int(min_sia_entry.get())}, {int(max_sia_entry.get())}\n")
                f.write(f"uronic_acids = {int(min_ua_entry.get())}, {int(max_ua_entry.get())}\n")
                f.write(f"fucoses = {int(min_dhex_entry.get())}, {int(max_dhex_entry.get())}\n")
                f.write(f"neu5ac = {int(min_neu5ac_entry.get())}, {int(max_neu5ac_entry.get())}\n")
                f.write(f"neu5gc = {int(min_neu5gc_entry.get())}, {int(max_neu5gc_entry.get())}\n")
            else:
                f.write(f"mode = custom_library\n")
                f.write(f"export_library = False\n")
                f.write(f"exported_library_name = \n")
                f.write(f"import_library_path = \n")
                custom_glycans_save_param = ''
                for i in local_custom_glycans_list[1]:
                    custom_glycans_save_param+=f"{i}, " 
                f.write(f"custom_glycans_list = {custom_glycans_save_param[:-2]}\n")
                f.write(f"total_monosaccharides = {int(min_monosaccharides_entry.get())}, {int(max_monosaccharides_entry.get())}\n")
                f.write(f"hexoses = {int(min_hex_entry.get())}, {int(max_hex_entry.get())}\n")
                f.write(f"hexosamines = {int(min_hn_entry.get())}, {int(max_hn_entry.get())}\n")
                f.write(f"hexnacs = {int(min_hexnac_entry.get())}, {int(max_hexnac_entry.get())}\n")
                f.write(f"xyloses = {int(min_xyl_entry.get())}, {int(max_xyl_entry.get())}\n")
                f.write(f"sialic_acids = {int(min_sia_entry.get())}, {int(max_sia_entry.get())}\n")
                f.write(f"uronic_acids = {int(min_ua_entry.get())}, {int(max_ua_entry.get())}\n")
                f.write(f"fucoses = {int(min_dhex_entry.get())}, {int(max_dhex_entry.get())}\n")
                f.write(f"neu5ac = {int(min_neu5ac_entry.get())}, {int(max_neu5ac_entry.get())}\n")
                f.write(f"neu5gc = {int(min_neu5gc_entry.get())}, {int(max_neu5gc_entry.get())}\n")
            
            custom_monosaccharides_string = f""
            for index, cm in enumerate(local_custom_monos_list):
                custom_monosaccharides_string += f"({cm['cm_name']}, {cm['cm_short_code']}, {cm['cm_chem_comp']}, {cm['cm_min']}, {cm['cm_max']}, {'yes' if cm['sialic'] else 'no'})"
                if index != len(local_custom_monos_list)-1:
                    custom_monosaccharides_string += f", "
                    
            f.write(f"custom_monosaccharides = {custom_monosaccharides_string}\n")
                
            f.write("\n")
            f.write(f"[common_library_building_settings]\n")
            f.write(f"force_class_structure = {forced_classes[forced_class_dropdown.get()]}\n")
            
            f.write(f"min_max_proton_adducts = {hydrogen_min_entry.get()}, {hydrogen_max_entry.get()}\n")
            
            custom_adducts_line_content = []
            for adduct in local_custom_adducts_list:
                custom_adducts_line_content.append(f"({adduct[0]}, {custom_adducts_widgets[adduct[0]][0].get()}, {custom_adducts_widgets[adduct[0]][1].get()}, {adduct[3]})")
            f.write(f"custom_adducts = {', '.join(custom_adducts_line_content)}\n")
            
            f.write(f"max_charges = {'-'+max_charges_entry.get() if ion_mode_choice.get() == 'neg' else max_charges_entry.get()}\n")
            f.write(f"reducing_end_tag = {reducing_end_tag_entry.get()}\n")
            f.write(f"permethylated = {permethylated_checkbox_state.get()}\n")
            f.write(f"reduced = {reduced_checkbox_state.get()}\n")
            f.write(f"min_max_sulfation_per_glycan = {min_sulfation_entry.get()}, {max_sulfation_entry.get()}\n")
            f.write(f"min_max_phosphorylation_per_glycan = {min_phosphorylation_entry.get()}, {max_phosphorylation_entry.get()}\n")
            f.write(f"lyase_digested = {lyase_digested_checkbox_state.get()}\n")
            f.write(f"aminated_ethyl_esterified = {lac_ee_checkbox_state.get()}\n")
            f.write(f"fast_iso = {fast_iso_checkbox_state.get()}\n")
            f.write(f"high_resolution_isotopic_dist = {hires_iso_checkbox_state.get()}\n")
            f.write(f"internal_standard_mass = {intstandard_entry.get()}\n")
            f.write("\n")
            f.write(f"[analysis_parameters]\n")
            f.write(f"analyze_ms2 = {local_analyze_ms2[0]}\n")
            f.write(f"force_fragments_to_glycans = {local_analyze_ms2[1]}\n")
            f.write(f"unrestricted_fragments = {local_analyze_ms2[2]}\n")
            
            export_acc_unit = ''
            if acc_unit_dropdown.get() == 'PPM':
                export_acc_unit = 'ppm'
            elif acc_unit_dropdown.get() == 'mz':
                export_acc_unit = 'mz'
                
            f.write(f"accuracy_unit = {export_acc_unit}\n")
            f.write(f"accuracy_value = {acc_value_entry.get()}\n")
            f.write(f"ret_time_interval = {rt_int_min_entry.get()}, {rt_int_max_entry.get()}\n")
            f.write(f"custom_min_points_per_peak = {custom_ppp_checkbox_state.get()}\n")
            f.write(f"number_points_per_peak = {int(custom_ppp_entry.get())}\n")
            f.write(f"limit_peaks_picked = {close_peaks_checkbox_state.get()}\n")
            f.write(f"max_number_peaks = {int(close_peaks_entry.get())}\n")
            f.write("\n")
            f.write(f"[post-analysis/reanalysis]\n")
            f.write(f"filter_ms2_by_reporter_ions = \n")
            f.write(f"align_chromatograms = False\n")
            f.write(f"auc_percentage_threshold = 1\n")
            f.write(f"minimum_samples = 0\n")
            f.write(f"max_ppm_threshold = {float(ppm_error_min_entry.get())}, {float(ppm_error_max_entry.get())}\n")
            f.write(f"isotopic_fitting_score_threshold = {float(iso_fit_entry.get())}\n")
            f.write(f"curve_fitting_score_threshold = {float(curve_fit_entry.get())}\n")
            f.write(f"signal_to_noise_threshold = {float(s_n_entry.get())}\n")
            f.write(f"output_compositions_analysis = True\n")
            f.write(f"output_metaboanalyst_file = False\n")
            f.write(f"metaboanalyst_groups = \n")
            f.write(f"output_fittings_data = False\n")
            f.write(f"output_plot_data = False")
            f.close()
        file_dialog.destroy()
        set_parameters_window.grab_set()
        
    def load_param_from_file():
        file_dialog = tk.Toplevel()
        file_dialog.attributes("-topmost", True)
        file_dialog.withdraw()
        file_dialog.grab_set()
        
        file_path = filedialog.askopenfilename(filetypes=[("INI files", "*.ini"), ("All files", "*.*")])
        if len(file_path) == 0:
            file_dialog.destroy()
            set_parameters_window.grab_set()
            return
        try:
            parameters = Config_Handler.config_handler(True, file_path)
        except Exception as e:
            file_dialog.destroy()
            error_window("Unable to read the configurations from the file you selected. Try again with a different file.")
            set_parameters_window.grab_set()
            return
                    
        # Insert adducts in the entries
        hydrogen_max_entry.delete(0, tk.END)
        hydrogen_min_entry.delete(0, tk.END)
        hydrogen_max_entry.insert(0, parameters['library settings']['min/max protons'][1])
        hydrogen_min_entry.insert(0, parameters['library settings']['min/max protons'][0])
        
        # Add the custom adducts
        for custom_adduct in list(custom_adducts_widgets.keys()):
            remove_custom_adduct(custom_adduct)
        for custom_adduct in parameters['library settings']['custom adducts']:
            create_custom_adduct_window(custom_adduct[0], custom_adduct[3])
            custom_adducts_widgets[custom_adduct[0]][0].delete(0, tk.END)
            custom_adducts_widgets[custom_adduct[0]][1].delete(0, tk.END)
            custom_adducts_widgets[custom_adduct[0]][0].insert(0, custom_adduct[1])
            custom_adducts_widgets[custom_adduct[0]][1].insert(0, custom_adduct[2])
            
        # Insert max charge and positive/negative mode
        max_charges_entry.delete(0, tk.END)
        max_charges_entry.insert(0, abs(parameters['library settings']['maximum charges']))
        if parameters['library settings']['maximum charges'] < 0:
            ion_mode_choice.set('neg')
        else:
            ion_mode_choice.set('pos')
            
        # Insert score values in the entry boxes
        ppm_error_min_entry.delete(0, tk.END)
        ppm_error_max_entry.delete(0, tk.END)
        if type(parameters['analysis settings']['maximum ppm threshold']) == float:
            ppm_error_min_entry.insert(0, 0-parameters['analysis settings']['maximum ppm threshold'])
            ppm_error_max_entry.insert(0, parameters['analysis settings']['maximum ppm threshold'])
        else:
            ppm_error_min_entry.insert(0, parameters['analysis settings']['maximum ppm threshold'][0])
            ppm_error_max_entry.insert(0, parameters['analysis settings']['maximum ppm threshold'][1])
        iso_fit_entry.delete(0, tk.END)
        iso_fit_entry.insert(0, parameters['analysis settings']['isotopic fitting score threshold'])
        curve_fit_entry.delete(0, tk.END)
        curve_fit_entry.insert(0, parameters['analysis settings']['curve fitting score threshold'])
        s_n_entry.delete(0, tk.END)
        s_n_entry.insert(0, parameters['analysis settings']['signal to noise ratio threshold'])
        
        # RT/MT range entry
        rt_int_min_entry.delete(0, tk.END)
        rt_int_min_entry.insert(0, parameters['analysis settings']['retention time interval'][0])
        rt_int_max_entry.delete(0, tk.END)
        rt_int_max_entry.insert(0, parameters['analysis settings']['retention time interval'][1])
        
        # Insert the accuracy value into its entry
        acc_value_entry.delete(0, tk.END)
        acc_value_entry.insert(0, parameters['analysis settings']['mass tolerance'][1])
        
        # Close peaks entry
        if parameters['analysis settings']['only x most intense peaks'][0]:
            close_peaks_checkbox_state.set(True)
            close_peaks_entry_state = tk.NORMAL
            close_peaks_entry.config(state = close_peaks_entry_state)
        else:
            close_peaks_checkbox_state.set(False)
            close_peaks_entry_state = tk.DISABLED
            close_peaks_entry.config(state = close_peaks_entry_state)
            
        # Monosaccharides entries
        min_monosaccharides_entry.delete(0, tk.END)
        min_monosaccharides_entry.insert(0, parameters['library settings']['min/max monosaccharides'][0])
        max_monosaccharides_entry.delete(0, tk.END)
        max_monosaccharides_entry.insert(0, parameters['library settings']['min/max monosaccharides'][1])
        min_hex_entry.delete(0, tk.END)
        min_hex_entry.insert(0, parameters['library settings']['min/max hexoses'][0])
        max_hex_entry.delete(0, tk.END)
        max_hex_entry.insert(0, parameters['library settings']['min/max hexoses'][1])
        min_hn_entry.delete(0, tk.END)
        min_hn_entry.insert(0, parameters['library settings']['min/max hexosamines'][0])
        max_hn_entry.delete(0, tk.END)
        max_hn_entry.insert(0, parameters['library settings']['min/max hexosamines'][1])
        min_hexnac_entry.delete(0, tk.END)
        min_hexnac_entry.insert(0, parameters['library settings']['min/max hexnac'][0])
        max_hexnac_entry.delete(0, tk.END)
        max_hexnac_entry.insert(0, parameters['library settings']['min/max hexnac'][1])
        min_xyl_entry.delete(0, tk.END)
        min_xyl_entry.insert(0, parameters['library settings']['min/max xyloses'][0])
        max_xyl_entry.delete(0, tk.END)
        max_xyl_entry.insert(0, parameters['library settings']['min/max xyloses'][1])
        min_dhex_entry.delete(0, tk.END)
        min_dhex_entry.insert(0, parameters['library settings']['min/max fucoses'][0])
        max_dhex_entry.delete(0, tk.END)
        max_dhex_entry.insert(0, parameters['library settings']['min/max fucoses'][1])
        min_sia_entry.delete(0, tk.END)
        min_sia_entry.insert(0, parameters['library settings']['min/max sialic acids'][0])
        max_sia_entry.delete(0, tk.END)
        max_sia_entry.insert(0, parameters['library settings']['min/max sialic acids'][1])
        min_ua_entry.delete(0, tk.END)
        min_ua_entry.insert(0, parameters['library settings']['min/max uronic acids'][0])
        max_ua_entry.delete(0, tk.END)
        max_ua_entry.insert(0, parameters['library settings']['min/max uronic acids'][1])
        min_neu5ac_entry.delete(0, tk.END)
        min_neu5ac_entry.insert(0, parameters['library settings']['min/max acetyl sialic acids'][0])
        max_neu5ac_entry.delete(0, tk.END)
        max_neu5ac_entry.insert(0, parameters['library settings']['min/max acetyl sialic acids'][1])
        min_neu5gc_entry.delete(0, tk.END)
        min_neu5gc_entry.insert(0, parameters['library settings']['min/max glycolyl sialic acids'][0])
        max_neu5gc_entry.delete(0, tk.END)
        max_neu5gc_entry.insert(0, parameters['library settings']['min/max glycolyl sialic acids'][1])
        
        # Internal standard entry
        intstandard_entry.delete(0, tk.END)
        intstandard_entry.insert(0, parameters['library settings']['internal standard'])
        
        # Custom points per peak entry
        if parameters['analysis settings']['min points per peak'][0]:
            custom_ppp_checkbox_state.set(True)
            custom_ppp_entry_state = tk.NORMAL
            custom_ppp_entry.config(state = custom_ppp_entry_state)
            custom_ppp_entry.delete(0, tk.END)
            custom_ppp_entry.insert(0, parameters['analysis settings']['min points per peak'][1])
        else:
            custom_ppp_checkbox_state.set(False)
            custom_ppp_entry_state = tk.DISABLED
            custom_ppp_entry.config(state = custom_ppp_entry_state)
            
        # Accuracy unit dropdown menu
        if parameters['analysis settings']['mass tolerance'][0] == 'ppm':
            acc_unit_dropdown.set(acc_unit_dropdown_options[0])
        elif parameters['analysis settings']['mass tolerance'][0] == 'mz':
            acc_unit_dropdown.set(acc_unit_dropdown_options[1])
            
        # Custom monosaccharides
        local_custom_monos_list[:] = parameters['library settings']['custom monosaccharides']
        
        # Working directory
        working_dir_label.config(text=parameters['analysis settings']['save path'])
        
        # If custom glycans
        if parameters['library settings']['custom glycans list'][0]:
            use_custom_library_checkbox_state.set(True)
            from_file_button_state = tk.NORMAL
            from_file_button.config(state = from_file_button_state)
            local_custom_glycans_list[0] = True
            local_custom_glycans_list[1] = parameters['library settings']['custom glycans list'][1]
            
            min_monosaccharides_entry.config(state=tk.DISABLED)
            max_monosaccharides_entry.config(state=tk.DISABLED)
            min_hex_entry.config(state=tk.DISABLED)
            max_hex_entry.config(state=tk.DISABLED)
            min_hn_entry.config(state=tk.DISABLED)
            max_hn_entry.config(state=tk.DISABLED)
            min_hexnac_entry.config(state=tk.DISABLED)
            max_hexnac_entry.config(state=tk.DISABLED)
            min_xyl_entry.config(state=tk.DISABLED)
            max_xyl_entry.config(state=tk.DISABLED)
            min_dhex_entry.config(state=tk.DISABLED)
            max_dhex_entry.config(state=tk.DISABLED)
            min_sia_entry.config(state=tk.DISABLED)
            max_sia_entry.config(state=tk.DISABLED)
            min_ua_entry.config(state=tk.DISABLED)
            max_ua_entry.config(state=tk.DISABLED)
            min_neu5ac_entry.config(state=tk.DISABLED)
            max_neu5ac_entry.config(state=tk.DISABLED)
            min_neu5gc_entry.config(state=tk.DISABLED)
            max_neu5gc_entry.config(state=tk.DISABLED)
        else:
            use_custom_library_checkbox_state.set(False)
            from_file_button_state = tk.DISABLED
            from_file_button.config(state = from_file_button_state)
            local_custom_glycans_list[0] = False
            
            min_monosaccharides_entry.config(state=tk.NORMAL)
            max_monosaccharides_entry.config(state=tk.NORMAL)
            min_hex_entry.config(state=tk.NORMAL)
            max_hex_entry.config(state=tk.NORMAL)
            min_hn_entry.config(state=tk.NORMAL)
            max_hn_entry.config(state=tk.NORMAL)
            min_hexnac_entry.config(state=tk.NORMAL)
            max_hexnac_entry.config(state=tk.NORMAL)
            min_xyl_entry.config(state=tk.NORMAL)
            max_xyl_entry.config(state=tk.NORMAL)
            min_dhex_entry.config(state=tk.NORMAL)
            max_dhex_entry.config(state=tk.NORMAL)
            min_sia_entry.config(state=tk.NORMAL)
            max_sia_entry.config(state=tk.NORMAL)
            min_ua_entry.config(state=tk.NORMAL)
            max_ua_entry.config(state=tk.NORMAL)
            min_neu5ac_entry.config(state=tk.NORMAL)
            max_neu5ac_entry.config(state=tk.NORMAL)
            min_neu5gc_entry.config(state=tk.NORMAL)
            max_neu5gc_entry.config(state=tk.NORMAL)
            
        # MS2 analysis
        if parameters['analysis settings']['analyze ms2'][0]:
            analyze_ms2_checkbox_state.set(True)
            force_ms2_comp_checkbox_active_state = tk.NORMAL
            force_ms2_comp_checkbox.config(state = force_ms2_comp_checkbox_active_state)
            if parameters['analysis settings']['analyze ms2'][1]:
                force_ms2_comp_checkbox_state.set(True)
            else:
                force_ms2_comp_checkbox_state.set(False)
            unrestricted_frags_checkbox_active_state = tk.NORMAL
            unrestricted_frags_checkbox.config(state = force_ms2_comp_checkbox_active_state)
            if parameters['analysis settings']['analyze ms2'][2]:
                unrestricted_frags_checkbox_state.set(True)
            else:
                unrestricted_frags_checkbox_state.set(False)
            local_analyze_ms2[0] = True
            local_analyze_ms2[1] = parameters['analysis settings']['analyze ms2'][1]
            local_analyze_ms2[2] = parameters['analysis settings']['analyze ms2'][2]
        else:
            analyze_ms2_checkbox_state.set(False)
            force_ms2_comp_checkbox_active_state = tk.DISABLED
            force_ms2_comp_checkbox.config(state = force_ms2_comp_checkbox_active_state)
            unrestricted_frags_checkbox_active_state = tk.DISABLED
            unrestricted_frags_checkbox.config(state = force_ms2_comp_checkbox_active_state)
            local_analyze_ms2[0] = False
            
        # Reducing end tag
        if parameters['library settings']['reducing end tag'] == 0.0 or parameters['library settings']['reducing end tag'] == '0.0':
            reducing_end_tag_dropdown.set('Custom')
            reducing_end_tag_dropdown.state(['disabled'])
            reducing_end_tag_checkbox_state.set(False)
            reducing_end_tag_entry_state = tk.DISABLED
            reducing_end_tag_entry.config(state = reducing_end_tag_entry_state)
            reduced_checkbox_activation_state = tk.NORMAL
            reduced_checkbox.config(state = reduced_checkbox_activation_state)
        elif parameters['library settings']['reducing end tag'] in reducing_end_tags.values():
            reducing_end_tag_checkbox_state.set(True)
            swapped_reducing_end_tags = {v: k for k, v in reducing_end_tags.items()}
            reducing_end_tag_dropdown.set(swapped_reducing_end_tags[parameters['library settings']['reducing end tag']])
            reducing_end_tag_dropdown.state(['!disabled'])
            reducing_end_tag_entry_state = tk.NORMAL
            reducing_end_tag_entry.config(state = reducing_end_tag_entry_state)
            reducing_end_tag_entry.delete(0, tk.END)
            reducing_end_tag_entry.insert(0, parameters['library settings']['reducing end tag'])
            reducing_end_tag_entry_state = tk.DISABLED
            reducing_end_tag_entry.config(state = reducing_end_tag_entry_state)
            reduced_checkbox_activation_state = tk.DISABLED
            reduced_checkbox.config(state = reduced_checkbox_activation_state)
        else:
            reducing_end_tag_checkbox_state.set(True)
            reducing_end_tag_dropdown.set('Custom')
            reducing_end_tag_entry_state = tk.NORMAL
            reducing_end_tag_entry.config(state = reducing_end_tag_entry_state)
            reducing_end_tag_entry.delete(0, tk.END)
            reducing_end_tag_entry.insert(0, parameters['library settings']['reducing end tag'])
            reduced_checkbox_activation_state = tk.DISABLED
            reduced_checkbox.config(state = reduced_checkbox_activation_state)
            
        # Permethylated
        if parameters['library settings']['permethylated']:
            permethylated_checkbox_state.set(True)
        else:
            permethylated_checkbox_state.set(False)
            
        # Reduced end
        if parameters['library settings']['reducing end reduced']:
            reduced_checkbox_state.set(True)
        else:
            reduced_checkbox_state.set(False)
            
        # Amidated and Ethyl-esterified
        if parameters['library settings']['lactonized/ethyl esterified']:
            lac_ee_checkbox_state.set(True)
        else:
            lac_ee_checkbox_state.set(False)
            
        # Glycan class
        if type(parameters['library settings']['glycan class']) == bool:
            if parameters['library settings']['glycan class']:
                forced_class_dropdown.set('N-Glycans')
            else:
                forced_class_dropdown.set('None')
        else:
            if parameters['library settings']['glycan class'] == 'n_glycans':
                forced_class_dropdown.set('N-Glycans')
            elif parameters['library settings']['glycan class'] == 'o_glycans':
                forced_class_dropdown.set('O-Glycans')
            elif parameters['library settings']['glycan class'] == 'gags':
                forced_class_dropdown.set('GAGs')
            else:
                forced_class_dropdown.set('None')
        
        # Fast isotopic envelope calculation
        if parameters['library settings']['fast isotopic pattern calculation']:
            fast_iso_checkbox_state.set(True)
        else:
            fast_iso_checkbox_state.set(False)
            
        # High-res isotopic envelope
        if parameters['library settings']['high resolution isotopic pattern']:
            hires_iso_checkbox_state.set(True)
        else:
            hires_iso_checkbox_state.set(False)
            
        # Multithreaded
        if parameters['analysis settings']['multithreaded']:
            multithreaded_checkbox_state.set(True)
            number_cores_entry_state = tk.NORMAL
            number_cores_entry.config(state = number_cores_entry_state)
            number_cores_entry.delete(0, tk.END)
            number_cores_entry.insert(0, parameters['analysis settings']['number of cpu cores'])
        else:
            multithreaded_checkbox_state.set(False)
            number_cores_entry_state = tk.DISABLED
            number_cores_entry.config(state = number_cores_entry_state)
        
        # Lyase digested
        if parameters['library settings']['lyase digested'] and parameters['library settings']['glycan class'] == 'gags':
            lyase_digested_checkbox.config(state=tk.NORMAL)
            lyase_digested_checkbox_state.set(True)
        else:
            lyase_digested_checkbox.config(state=tk.DISABLED)
            lyase_digested_checkbox_state.set(False)
            
        # Sulfation/phosphorylation
        min_sulfation_entry.delete(0, tk.END)
        min_sulfation_entry.insert(0, parameters['library settings']['min/max sulfation'][0])
        max_sulfation_entry.delete(0, tk.END)
        max_sulfation_entry.insert(0, parameters['library settings']['min/max sulfation'][1])
        min_phosphorylation_entry.delete(0, tk.END)
        min_phosphorylation_entry.insert(0, parameters['library settings']['min/max phosphorylation'][0])
        max_phosphorylation_entry.delete(0, tk.END)
        max_phosphorylation_entry.insert(0, parameters['library settings']['min/max phosphorylation'][1])
        
        file_dialog.destroy()
        set_parameters_window.grab_set()
        
    def open_file_dialog_spw_custom_glycan():
        global custom_glycans_text, custom_glycans_window
        file_dialog = tk.Toplevel()
        file_dialog.attributes("-topmost", True)
        file_dialog.withdraw()
        file_dialog.grab_set()
        
        file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        temp_custom_glycans_list = ""
        if file_path != "":
            with open(file_path, 'r') as f:
                for i in f:
                    temp_custom_glycans_list += i
                f.close()
                
            if len(temp_custom_glycans_list.strip()) == 0:
                error_window("Custom glycans file is incorrect. Make sure the glycans are comma-separated or line-separated.")
                file_dialog.destroy()
                custom_glycans_window.grab_set()
                return
                
            # Attempt to split by commas
            ocgw_custom_glycans = temp_custom_glycans_list.split(",")
            
            # Attempt to split by new lines
            if len(ocgw_custom_glycans) == 1:
                ocgw_custom_glycans = ocgw_custom_glycans[0].split("\n")
            
            # Cleanup
            to_remove = []
            
            # Remove empty entries
            for i_i, i in enumerate(ocgw_custom_glycans):
                ocgw_custom_glycans[i_i] = i.strip()
                if len(i) == 0:
                    to_remove.append(i_i)
            for i in sorted(to_remove, reverse = True):
                del ocgw_custom_glycans[i]
                
            # Delete the current text in the entry
            custom_glycans_text.delete(1.0, tk.END)
            
            # Insert the new glycans
            for i_i, i in enumerate(ocgw_custom_glycans):
                custom_glycans_text.insert(tk.END, i)
                if i_i != len(ocgw_custom_glycans)-1:
                    custom_glycans_text.insert(tk.END, ", ")
                    
        file_dialog.destroy()
        custom_glycans_window.grab_set()
    
    def open_custom_glycans_window():
        global custom_glycans_text, custom_glycans_window
        def ok_ocg_window():
            ocgw_custom_glycans = custom_glycans_text.get("1.0", "end-1c").split(",")
            if len(ocgw_custom_glycans) == 1:
                ocgw_custom_glycans = ocgw_custom_glycans[0].split("\n")
            to_remove = []
            to_add = []
            if len(ocgw_custom_glycans) >= 1:
                for i_i, i in enumerate(ocgw_custom_glycans):
                    ocgw_custom_glycans[i_i] = i.strip().strip("'").strip('"')
                    if len(ocgw_custom_glycans[i_i]) == 0:
                        to_remove.append(i_i)
                        continue
                    if len(ocgw_custom_glycans[i_i].split("/")) > 1:
                        splitted_glycan = ocgw_custom_glycans[i_i].split("/")
                        ocgw_custom_glycans[i_i] = splitted_glycan[0]
                        to_add.append(splitted_glycan[1])
            for i in sorted(to_remove, reverse = True):
                del ocgw_custom_glycans[i]
            for i in to_add:
                ocgw_custom_glycans.append(i)
            for i in ocgw_custom_glycans:
                glycan_comp = General_Functions.form_to_comp_glycans(i)
                for i in glycan_comp:
                    if i == 'Am' or i == 'E' or i == 'AmG' or i == 'EG':
                        lac_ee_checkbox_state.set(True)
                    if i not in General_Functions.monosaccharides and i not in [cm['cm_short_code'] for cm in local_custom_monos_list]:
                        error_window(f"Unrecognized monosaccharide in glycan list: {i}\nCheck your custom glycans list.")
                        return
            local_custom_glycans_list[1] = ocgw_custom_glycans
            custom_glycans_window.destroy()
            set_parameters_window.grab_set()
        
        custom_glycans_window = tk.Toplevel()
        #custom_glycans_window.attributes("-topmost", True)
        custom_glycans_window.withdraw()
        custom_glycans_window.title("Custom Glycans List")
        icon = ImageTk.PhotoImage(ico_image)
        custom_glycans_window.iconphoto(False, icon)
        custom_glycans_window.resizable(False, False)
        custom_glycans_window.grab_set()
        custom_glycans_window.geometry("600x625")
        custom_glycans_window.columnconfigure(0, weight=1)
        
        description_custom_glycans_label = ttk.Label(custom_glycans_window, text=f"Here you can type in the glycans you wish to insert in the custom glycans list. Use the following nomenclature:\n\n    H: Hexose, HN: Hexosamine, N: HexNAc, X: Xylose, F: Deoxyhexose, S: Neu5Ac, G: Neu5Gc, UA: Uronic Acid {''.join([', '+cm['cm_short_code']+': '+cm['cm_name'] for cm in local_custom_monos_list])}\n\nand in case of amidation/ethyl-esterification of sialic acids:\n\n    Am: Amidated Neu5Ac (alpha2,3), E: Ethyl-Esterified Neu5Ac (alpha2,6)\n    AmG: Amidated Neu5Gc (alpha2,3), EG: Ethyl-Esterified Neu5Gc (alpha2,6)\n\nExample: H5N4S2F1 refers to a glycan with 5 Hexoses, 4 HexNacs, 2 Neu5Ac and 1 Deoxyhexose", font=("Segoe UI", list_font_size), wraplength=600, justify = "left")
        description_custom_glycans_label.grid(row=0, column=0, padx=(10, 10), pady=(10, 10), sticky="nsew")
        
        custom_glycans_text = ScrolledText(custom_glycans_window, height=22, width=100)
        custom_glycans_text.grid(row=1, column=0, padx=(10, 10), pady=(0, 10), sticky="nsew")
        for i_i, i in enumerate(local_custom_glycans_list[1]):
            custom_glycans_text.insert(tk.END, i)
            if i_i != len(local_custom_glycans_list[1])-1:
                custom_glycans_text.insert(tk.END, ", ")
        
        ff_cg_window_button = ttk.Button(custom_glycans_window, text="From File", style="small_button_spw_style1.TButton", command=open_file_dialog_spw_custom_glycan)
        ff_cg_window_button.grid(row=2, column=0, padx=(10, 10), pady=(0,10), sticky="w")
        ToolTip(ff_cg_window_button, "Allows the use of a list of glycans from a text file, comma- or line-separated.")
        
        ok_cg_window_button = ttk.Button(custom_glycans_window, text="Ok", style="small_button_spw_style1.TButton", command=ok_ocg_window)
        ok_cg_window_button.grid(row=2, column=0, padx=(10, 10), pady=(0,10), sticky="e")
        
        custom_glycans_window.update_idletasks()
        custom_glycans_window.deiconify()
        window_width = custom_glycans_window.winfo_width()
        window_height = custom_glycans_window.winfo_height()
        screen_width = custom_glycans_window.winfo_screenwidth()
        screen_height = custom_glycans_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        custom_glycans_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
    def open_custom_monosaccharides_window():
        def ok_ocm_window(local_custom_monos_list):
            temp_custom_mono = []
            for index, row in enumerate(rows_widgets):
                if len(row['entry_full_name'].get()) > 0:
                    
                    single_letters_in_use = set([x[-1] for x in General_Functions.monosaccharides.values()]+[x['cm_single_letter_code'] for x in temp_custom_mono]+[x['cm_single_letter_code'] for x in custom_monosaccharides]+['T'])
                    
                    short_codes_in_use = set([x.upper() for x in General_Functions.monosaccharides.keys()]+[x['cm_short_code'] for x in temp_custom_mono]+[x['cm_short_code'] for x in custom_monosaccharides]+['T'])
                    
                    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
                    for letter in alphabet:
                        if letter not in single_letters_in_use:
                            single_letter_code = letter
                            break
                            
                    if row['entry_short_code'].get().upper().strip() in short_codes_in_use:
                        error_window(f"Short code for glycan {row['entry_full_name'].get()} already used by other monosaccharides. Please, input a different one.")
                        return
                        
                    if len(row['entry_short_code'].get().strip()) > 3:
                        error_window(f"Short code for glycan {row['entry_full_name'].get()} is too long. Please, input a short code with up to three letters.")
                        return
                        
                    if int(row['entry_min'].get()) > int(row['entry_max'].get()):
                        error_window(f"Error in glycan {row['entry_full_name'].get()}. Maximum amount of the monosaccharide can't be smaller than the minimum amount.")
                        return
                        
                    try:
                        test_chem_comp = General_Functions.form_to_comp_atoms(row['entry_chemical_comp'].get().strip())
                        for atom in test_chem_comp.keys():
                            if atom not in mass.nist_mass.keys():
                                error_window(f"Something is wrong with the chemical formula of the monosaccharide {row['entry_full_name'].get()}. Atom {atom} not identified.")
                                return
                    except:
                        error_window(f"Something is wrong with the chemical formula of the monosaccharide {row['entry_full_name'].get()}.")
                        return
                        
                    for number in range(10):
                        if str(number) in row['entry_short_code'].get().strip():
                            error_window(f"Short code or single letter code for glycan {row['entry_full_name'].get()} has a number in it. Only letters are accepted.")
                            return
                    temp_custom_mono.append({
                        'cm_name' : row['entry_full_name'].get().strip(), 
                        'cm_short_code' : row['entry_short_code'].get().strip(), 
                        'cm_single_letter_code' : single_letter_code, 
                        'cm_chem_comp' : row['entry_chemical_comp'].get().strip(), 
                        'cm_min' : int(row['entry_min'].get()), 
                        'cm_max': int(row['entry_max'].get()), 
                        'sialic' : checkbox_states[index].get()
                        })
            
            # Update custom glycans list so no conflict occurs, in case custom monosaccharide used in custom glycans is removed
            removed_monos = []
            for mono in local_custom_monos_list:
                if mono not in temp_custom_mono:
                    removed_monos.append(mono['cm_short_code'])
            to_remove = []        
            for index, glycan in enumerate(local_custom_glycans_list[1]):
                glycan_comp = General_Functions.form_to_comp_glycans(glycan)
                for mono in glycan_comp:
                    if mono in removed_monos:
                        to_remove.append(index)
            for index in sorted(to_remove, reverse=True):
                del local_custom_glycans_list[1][index]
            
            local_custom_monos_list[:] = temp_custom_mono
            
            custom_monosaccharides_window.destroy()
            set_parameters_window.grab_set()
            
        def add_new_row(rows_widgets, checkbox_states):
            global row_number
            
            checkbox_states.append(tk.BooleanVar(value=False))
            default_row = {
                'entry_full_name' : ttk.Entry(cm_scrollable_frame, width=15),
                'entry_short_code' : ttk.Entry(cm_scrollable_frame, width=7),
                'entry_chemical_comp' : ttk.Entry(cm_scrollable_frame, width=20),
                'entry_min' : ttk.Spinbox(cm_scrollable_frame, width=3, from_=0, to=99),
                'entry_max' : ttk.Spinbox(cm_scrollable_frame, width=3, from_=0, to=99),
                'checkbox' : tk.Checkbutton(cm_scrollable_frame, variable=checkbox_states[-1])
                }
                
            default_row['entry_min'].insert(0, '0')
            default_row['entry_max'].insert(0, '1')
            
            for column, widget in enumerate(default_row.values()):
                widget.grid(row=row_number, column=column, padx=10, pady=5)
                
                if widget.winfo_class() != "TSpinbox":
                    if current_os == "Windows":
                        widget.bind("<MouseWheel>", cm_scroll_off_bar)
                    else:
                        widget.bind("<Button-4>", cm_scroll_off_bar)    # Linux Scroll Up
                        widget.bind("<Button-5>", cm_scroll_off_bar)
                
            rows_widgets.append(default_row)
            row_number += 1
        
        def cm_scroll_off_bar(event):
            """Handles mouse wheel scrolling based on OS."""
            if event.delta:  # Windows & macOS
                canvas_cm_window.yview_scroll(-1 * (event.delta // 120), "units")
            else:  # Linux (event.num is used instead of event.delta)
                if event.num == 4:  # Scroll up
                    canvas_cm_window.yview_scroll(-1, "units")
                elif event.num == 5:  # Scroll down
                    canvas_cm_window.yview_scroll(1, "units")
        
        custom_monosaccharides_window = tk.Toplevel(main_window)
        #custom_glycans_window.attributes("-topmost", True)
        custom_monosaccharides_window.withdraw()
        custom_monosaccharides_window.title("Custom Monosaccharides")
        icon = ImageTk.PhotoImage(ico_image)
        custom_monosaccharides_window.iconphoto(False, icon)
        custom_monosaccharides_window.resizable(False, False)
        custom_monosaccharides_window.grab_set()
        custom_monosaccharides_window.geometry("580x450")
        custom_monosaccharides_window.columnconfigure(0, weight=1)
        
        # Create a main frame to contain everything
        cm_main_frame = ttk.Frame(custom_monosaccharides_window)
        cm_main_frame.pack(fill="both", expand=True)
        
        # Canvas and scrollbar for the scrollable frame
        canvas_cm_window = tk.Canvas(cm_main_frame, borderwidth=0, highlightthickness=0)
        cm_scrollbar = ttk.Scrollbar(cm_main_frame, orient="vertical", command=canvas_cm_window.yview)
        cm_scrollable_frame = ttk.Frame(custom_monosaccharides_window, borderwidth=0)
        
        current_os = platform.system()
        if current_os == "Windows":
            cm_scrollable_frame.bind("<MouseWheel>", cm_scroll_off_bar)
        else:
            cm_scrollable_frame.bind("<Button-4>", cm_scroll_off_bar)    # Linux Scroll Up
            cm_scrollable_frame.bind("<Button-5>", cm_scroll_off_bar)
        
        canvas_cm_window.configure(yscrollcommand=cm_scrollbar.set)
        cm_scrollbar.pack(side="right", fill="y")
        canvas_cm_window.pack(side="left", fill="both", expand=True)
        
        canvas_cm_window.create_window((0, 0), window=cm_scrollable_frame, anchor="nw")
        cm_scrollable_frame.bind("<Configure>", lambda e: canvas_cm_window.configure(scrollregion=canvas_cm_window.bbox("all")))
        
        # Title labels
        column1_title_label = tk.Label(cm_scrollable_frame, text="Full Name", font=("Segoe UI", 10, "bold"))
        column1_title_label.grid(row=0, column=0, sticky="e", padx=20)
        ToolTip(column1_title_label, "The full name of the monosaccharide. It will be displayed in the XLSX library file.")
        
        column2_title_label = tk.Label(cm_scrollable_frame, text="Short-Code", font=("Segoe UI", 10, "bold"))
        column2_title_label.grid(row=0, column=1, sticky="w", padx=5)
        ToolTip(column2_title_label, "The short code for the monosaccharide. It will be used in the glycan formula (e.g. in H5N4Am1 the monosaccharides 'H', 'N' and 'Am' are in short code).")
        
        column4_title_label = tk.Label(cm_scrollable_frame, text="Chemical Composition", font=("Segoe UI", 10, "bold"))
        column4_title_label.grid(row=0, column=2, sticky="w", padx=10)
        ToolTip(column4_title_label, "The residue chemical composition of the glycan. If only one of a given atom is present, the number must still be indicated (e.g. C8O5NH13 must be inputted as C8O5N1H13).")
        
        column5_title_label = tk.Label(cm_scrollable_frame, text="Min", font=("Segoe UI", 10, "bold"))
        column5_title_label.grid(row=0, column=3, sticky="w", padx=5)
        ToolTip(column5_title_label, "The minimum amount of this monosaccharide in glycans in your library.")
        
        column6_title_label = tk.Label(cm_scrollable_frame, text="Max", font=("Segoe UI", 10, "bold"))
        column6_title_label.grid(row=0, column=4, sticky="w", padx=5)
        ToolTip(column6_title_label, "The maximum amount of this monosaccharide in glycans in your library.")
        
        column7_title_label = tk.Label(cm_scrollable_frame, text="Sialic Acid?", font=("Segoe UI", 10, "bold"))
        column7_title_label.grid(row=0, column=5, sticky="w", padx=5)
        ToolTip(column7_title_label, "Whether or not the monosaccharide is a sialic acid. This will impact internal calculations when building the library.")
        
        # Create labels and entry widgets and pack them into the scrollable frame
        rows_widgets = []
        checkbox_states = []
        
        global row_number
        row_number = 1
        
        for custom_monosaccharide in local_custom_monos_list:
            checkbox_states.append(tk.BooleanVar(value=custom_monosaccharide['sialic']))
            default_row = {
                'entry_full_name' : ttk.Entry(cm_scrollable_frame, width=15),
                'entry_short_code' : ttk.Entry(cm_scrollable_frame, width=7),
                'entry_chemical_comp' : ttk.Entry(cm_scrollable_frame, width=20),
                'entry_min' : ttk.Spinbox(cm_scrollable_frame, width=3, from_=0, to=99),
                'entry_max' : ttk.Spinbox(cm_scrollable_frame, width=3, from_=0, to=99),
                'checkbox' : tk.Checkbutton(cm_scrollable_frame, variable=checkbox_states[-1])
                }
                
            default_row['entry_full_name'].insert(0, custom_monosaccharide['cm_name'])
            default_row['entry_short_code'].insert(0, custom_monosaccharide['cm_short_code'])
            default_row['entry_chemical_comp'].insert(0, custom_monosaccharide['cm_chem_comp'])
            default_row['entry_min'].insert(0, custom_monosaccharide['cm_min'])
            default_row['entry_max'].insert(0, custom_monosaccharide['cm_max'])
            
            for column, widget in enumerate(default_row.values()):
                widget.grid(row=row_number, column=column, padx=10, pady=5)
                
                if widget.winfo_class() != "TSpinbox":
                    if current_os == "Windows":
                        widget.bind("<MouseWheel>", cm_scroll_off_bar)
                    else:
                        widget.bind("<Button-4>", cm_scroll_off_bar)    # Linux Scroll Up
                        widget.bind("<Button-5>", cm_scroll_off_bar)
                
            rows_widgets.append(default_row)
                
            row_number += 1
            
        checkbox_states.append(tk.BooleanVar(value=False))
        default_row = {
            'entry_full_name' : ttk.Entry(cm_scrollable_frame, width=15),
            'entry_short_code' : ttk.Entry(cm_scrollable_frame, width=7),
            'entry_chemical_comp' : ttk.Entry(cm_scrollable_frame, width=20),
            'entry_min' : ttk.Spinbox(cm_scrollable_frame, width=3, from_=0, to=99),
            'entry_max' : ttk.Spinbox(cm_scrollable_frame, width=3, from_=0, to=99),
            'checkbox' : tk.Checkbutton(cm_scrollable_frame, variable=checkbox_states[-1])
            }
        default_row['entry_min'].insert(0, '0')
        default_row['entry_max'].insert(0, '1')
        
        for column, widget in enumerate(default_row.values()):
            widget.grid(row=row_number, column=column, padx=10, pady=5)
            
            if widget.winfo_class() != "TSpinbox":
                if current_os == "Windows":
                    widget.bind("<MouseWheel>", cm_scroll_off_bar)
                else:
                    widget.bind("<Button-4>", cm_scroll_off_bar)    # Linux Scroll Up
                    widget.bind("<Button-5>", cm_scroll_off_bar)
        
        rows_widgets.append(default_row)
            
        row_number += 1
        
        # Add a new frame for the OK button and pack it at the bottom
        buttons_frame = ttk.Frame(custom_monosaccharides_window)
        buttons_frame.pack(side="bottom", fill="x")
        
        ok_cm_window_button = ttk.Button(buttons_frame, text="Ok", style="small_button_spw_style1.TButton", command=lambda:ok_ocm_window(local_custom_monos_list))
        ok_cm_window_button.pack(side="right", anchor="se", padx=(10, 30), pady=(10, 15))
        
        add_new_row_button = ttk.Button(buttons_frame, text="Add New Monosaccharide", style="small_button_spw_style1.TButton", command=lambda: add_new_row(rows_widgets, checkbox_states))
        add_new_row_button.pack(side="left", anchor="sw", padx=(30, 10), pady=(10, 15))
        
        custom_monosaccharides_window.update_idletasks()
        custom_monosaccharides_window.deiconify()
        window_width = custom_monosaccharides_window.winfo_width()
        window_height = custom_monosaccharides_window.winfo_height()
        screen_width = custom_monosaccharides_window.winfo_screenwidth()
        screen_height = custom_monosaccharides_window.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        custom_monosaccharides_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
    def edit_acc_value(event):
        previous_value = float(acc_value_entry.get())
        if acc_unit_dropdown.get() == 'mz':
            acc_value_entry.configure(increment=0.01)
            acc_value_entry.delete(0, tk.END)
            new_value = (previous_value*1000)/(10**6)
            acc_value_entry.insert(0, new_value)
        else:
            acc_value_entry.configure(increment=1)
            acc_value_entry.delete(0, tk.END)
            new_value = int((previous_value/1000)*(10**6))
            acc_value_entry.insert(0, new_value)
            
    def create_standard_library_param_window():    
        def on_focus_loss(event):
            if not any(widget.winfo_exists() and widget.focus_displayof() for widget in std_param_window.winfo_children()):
                std_param_window.destroy()
                
        def set_parameters(params):
            # Custom library settings
            use_custom_library_checkbox_state.set(False)
            from_file_button_state = tk.DISABLED
            from_file_button.config(state = from_file_button_state)
            local_custom_glycans_list[0] = False
            
            # Activate monosaccharide number entries
            min_monosaccharides_entry.config(state=tk.NORMAL)
            max_monosaccharides_entry.config(state=tk.NORMAL)
            min_hex_entry.config(state=tk.NORMAL)
            max_hex_entry.config(state=tk.NORMAL)
            min_hn_entry.config(state=tk.NORMAL)
            max_hn_entry.config(state=tk.NORMAL)
            min_hexnac_entry.config(state=tk.NORMAL)
            max_hexnac_entry.config(state=tk.NORMAL)
            min_xyl_entry.config(state=tk.NORMAL)
            max_xyl_entry.config(state=tk.NORMAL)
            min_dhex_entry.config(state=tk.NORMAL)
            max_dhex_entry.config(state=tk.NORMAL)
            min_sia_entry.config(state=tk.NORMAL)
            max_sia_entry.config(state=tk.NORMAL)
            min_ua_entry.config(state=tk.NORMAL)
            max_ua_entry.config(state=tk.NORMAL)
            min_neu5ac_entry.config(state=tk.NORMAL)
            max_neu5ac_entry.config(state=tk.NORMAL)
            min_neu5gc_entry.config(state=tk.NORMAL)
            max_neu5gc_entry.config(state=tk.NORMAL)
            
            # Set monosaccharide numbers
            min_monosaccharides_entry.delete(0, tk.END)
            min_monosaccharides_entry.insert(0, params['monos'][0])
            max_monosaccharides_entry.delete(0, tk.END)
            max_monosaccharides_entry.insert(0, params['monos'][1])
            min_hex_entry.delete(0, tk.END)
            min_hex_entry.insert(0, params['h'][0])
            max_hex_entry.delete(0, tk.END)
            max_hex_entry.insert(0, params['h'][1])
            min_hn_entry.delete(0, tk.END)
            min_hn_entry.insert(0, params['hn'][0])
            max_hn_entry.delete(0, tk.END)
            max_hn_entry.insert(0, params['hn'][1])
            min_hexnac_entry.delete(0, tk.END)
            min_hexnac_entry.insert(0, params['n'][0])
            max_hexnac_entry.delete(0, tk.END)
            max_hexnac_entry.insert(0, params['n'][1])
            min_xyl_entry.delete(0, tk.END)
            min_xyl_entry.insert(0, params['x'][0])
            max_xyl_entry.delete(0, tk.END)
            max_xyl_entry.insert(0, params['x'][1])
            min_dhex_entry.delete(0, tk.END)
            min_dhex_entry.insert(0, params['f'][0])
            max_dhex_entry.delete(0, tk.END)
            max_dhex_entry.insert(0, params['f'][1])
            min_sia_entry.delete(0, tk.END)
            min_sia_entry.insert(0, params['sia'][0])
            max_sia_entry.delete(0, tk.END)
            max_sia_entry.insert(0, params['sia'][1])
            min_ua_entry.delete(0, tk.END)
            min_ua_entry.insert(0, params['ua'][0])
            max_ua_entry.delete(0, tk.END)
            max_ua_entry.insert(0, params['ua'][1])
            min_neu5ac_entry.delete(0, tk.END)
            min_neu5ac_entry.insert(0, params['s'][0])
            max_neu5ac_entry.delete(0, tk.END)
            max_neu5ac_entry.insert(0, params['s'][1])
            min_neu5gc_entry.delete(0, tk.END)
            min_neu5gc_entry.insert(0, params['g'][0])
            max_neu5gc_entry.delete(0, tk.END)
            max_neu5gc_entry.insert(0, params['g'][1])
            
            # Set glycan class
            lyase_digested_checkbox.config(state=tk.NORMAL)
            lyase_digested_checkbox_state.set(False)
            lyase_digested_checkbox.config(state=tk.DISABLED)
            
            if params['class'] == 'n_glycans':
                forced_class_dropdown.set('N-Glycans')
            elif params['class'] == 'o_glycans':
                forced_class_dropdown.set('O-Glycans')
            elif params['class'] == 'gags':
                forced_class_dropdown.set('GAGs')
                lyase_digested_checkbox.config(state=tk.NORMAL)
                lyase_digested_checkbox_state.set(params['lyase'])
            else:
                forced_class_dropdown.set('None')
        
        # Mouse cursor location
        x = main_window.winfo_pointerx()+15
        y = main_window.winfo_pointery()+15
        
        # Screen size
        screen_width = main_window.winfo_screenwidth()
        screen_height = main_window.winfo_screenheight()
        
        # Widget size
        widget_width = 150
        widget_height = 80
        
        # Adjustment for tooltip position on screen
        if x + widget_width > screen_width:  # Tooltip goes off the screen
            x = main_window.winfo_pointerx()-15 - widget_width  # Move to the left side
        if y + widget_height > screen_height:  # Tooltip goes off the screen
            y = main_window.winfo_pointery()-15 - widget_height  # Move to the left side
        
        # Window creation
        std_param_window = tk.Toplevel(set_parameters_window)
        std_param_window.overrideredirect(True)
        std_param_window.geometry(f"+{x}+{y}")
        std_param_window.attributes('-topmost', True)
        std_param_window.bind("<FocusOut>", on_focus_loss)
        
        std_param_window.focus_set()

        std_param_window_frame = tk.Frame(std_param_window, highlightbackground="grey", highlightthickness=1)
        std_param_window_frame.pack(fill=tk.BOTH, expand=True)
        
        # Buttons creation
        buttons = {}
        options_count = 0
        for option, params in standard_library_parameters.items():
            buttons[option] = ttk.Button(std_param_window_frame, text=option, command=lambda p=params: set_parameters(p))
            buttons[option].grid(row=options_count//3, column=options_count%3, pady=5, padx=5)
            options_count += 1
            
        # Warning label
        warning_label = ttk.Label(std_param_window_frame, text="Adjust the parameters manually based on\nyour sample for optimal performance.", foreground="red", font=("Segoe UI", list_font_size, "bold"))
        warning_label.grid(row=(options_count//3)+1, column=0, columnspan=3, pady=(0, 5), padx=5)
            
    def create_custom_adduct_window(chem_form = None, charges = None):    
        def on_focus_loss(event):
            if not any(widget.winfo_exists() and widget.focus_displayof() for widget in custom_adduct_window.winfo_children()):
                custom_adduct_window.destroy()
                
        def add_adduct(chem_form, charges):
            if type(chem_form) != str:
                chem_form = chem_form.get()
                charges = charges.get()
            
            row_number = 18
            for adduct, widgets in custom_adducts_widgets.items():
                if widgets[0].grid_info()['row'] > row_number:
                    row_number = widgets[0].grid_info()['row']
                    
            if chem_form in custom_adducts_widgets.keys():
                error_window("You can't add the same adduct twice!")
                return
            
            if len(chem_form) == 0:
                error_window("Please, insert a chemical formula for the custom adduct before confirming or click outside the window to cancel the addition")
                return
                
            if len(charges) == 0 or charges == '0' or charges == 0:
                error_window("Please, insert a valid value for the adduct charge.")
                return
            
            custom_adducts_widgets[chem_form] = []
            custom_adducts_widgets[chem_form].append(ttk.Spinbox(library_building_frame, width=4, from_=0, to=99))
            custom_adducts_widgets[chem_form][-1].grid(row=row_number+1, column=1, padx=(25, 0), pady=(0, 0), sticky='w')
            custom_adducts_widgets[chem_form][-1].insert(0, 0)
            ToolTip(custom_adducts_widgets[chem_form][-1], "Set the minimum number of this adduct.")
            
            custom_adducts_widgets[chem_form].append(ttk.Spinbox(library_building_frame, width=4, from_=0, to=99))
            custom_adducts_widgets[chem_form][-1].grid(row=row_number+1, column=1, padx=(75, 0), pady=(0, 0), sticky='w')
            custom_adducts_widgets[chem_form][-1].insert(0, 0)
            ToolTip(custom_adducts_widgets[chem_form][-1], "Set the maximum number of this adduct.")
            
            custom_adducts_widgets[chem_form].append(ttk.Label(library_building_frame, text=chem_form, font=("Segoe UI", list_font_size)))
            custom_adducts_widgets[chem_form][-1].grid(row=row_number+1, column=1, padx=(125, 10), pady=(0, 0), sticky="w")
            
            custom_adducts_widgets[chem_form].append(tk.Button(library_building_frame, text='x', command=lambda adduct=chem_form: remove_custom_adduct(adduct), relief='flat', padx=0, pady=0))
            custom_adducts_widgets[chem_form][-1].grid(row=row_number+1, column=1, padx=(10, 0), pady=(0, 0), sticky='w')
            
            if add_custom_adduct_button.grid_info()['row'] < 21:
                add_custom_adduct_button.grid_configure(row=add_custom_adduct_button.grid_info()['row']+1)
            else:
                add_custom_adduct_button.grid_remove()
                
            local_custom_adducts_list.append([chem_form, 0, 0, charges])
                
            try:
                custom_adduct_window.destroy()
            except:
                pass
                
        if chem_form and charges:
            add_adduct(chem_form, charges)
            return
        
        # Mouse cursor location
        x = main_window.winfo_pointerx()+15
        y = main_window.winfo_pointery()+15
        
        # Screen size
        screen_width = main_window.winfo_screenwidth()
        screen_height = main_window.winfo_screenheight()
        
        # Widget size
        widget_width = 150
        widget_height = 80
        
        # Adjustment for tooltip position on screen
        if x + widget_width > screen_width:  # Tooltip goes off the screen
            x = main_window.winfo_pointerx()-15 - widget_width  # Move to the left side
        if y + widget_height > screen_height:  # Tooltip goes off the screen
            y = main_window.winfo_pointery()-15 - widget_height  # Move to the left side
        
        # Window creation
        custom_adduct_window = tk.Toplevel(set_parameters_window)
        custom_adduct_window.overrideredirect(True)
        custom_adduct_window.geometry(f"+{x}+{y}")
        custom_adduct_window.attributes('-topmost', True)
        custom_adduct_window.bind("<FocusOut>", on_focus_loss)
        
        custom_adduct_window.focus_set()

        custom_adduct_window_frame = tk.Frame(custom_adduct_window, highlightbackground="grey", highlightthickness=1)
        custom_adduct_window_frame.pack(fill=tk.BOTH, expand=True)
        
        label_chem_form = ttk.Label(custom_adduct_window_frame, text="Chem. Formula:")
        label_chem_form.grid(row=0, column=0, padx=(10,10), pady=(5,5), sticky="nsew")
        ToolTip(label_chem_form, "Insert the chemical formula for the adduct.")
        
        label_charges = ttk.Label(custom_adduct_window_frame, text="Charges:")
        label_charges.grid(row=0, column=1, padx=(10,10), pady=(5,5), sticky="nsew")
        ToolTip(label_charges, "Insert the charge of each adduct molecule. Can be a positive or negative value.")
        
        label_chem_form_entry = ttk.Entry(custom_adduct_window_frame)
        label_chem_form_entry.grid(row=1, column=0, padx=(10,10), pady=(0,5), sticky="nsew")
        ToolTip(label_chem_form_entry, "Insert the chemical formula for the adduct.")
        
        label_charges_entry = ttk.Spinbox(custom_adduct_window_frame, width=4, state=tk.NORMAL, from_=-99, to=99)
        label_charges_entry.grid(row=1, column=1, padx=(10,10), pady=(0,5), sticky="nsew")
        label_charges_entry.insert(0, '0')
        ToolTip(label_charges_entry, "Insert the charge of each adduct molecule. Can be a positive or negative value.")
        
        confirm_button = ttk.Button(custom_adduct_window_frame, text="Add Adduct", command=lambda: add_adduct(label_chem_form_entry, label_charges_entry))
        confirm_button.grid(row=2, column=0, columnspan=2, padx=(10,10), pady=(0,5))
        
    def remove_custom_adduct(adduct):
        row = custom_adducts_widgets[adduct][0].grid_info()['row']
        
        for widget in custom_adducts_widgets[adduct]:
            widget.destroy()
        del custom_adducts_widgets[adduct]
        
        for adduct_move, widgets in custom_adducts_widgets.items():
            if widgets[0].grid_info()['row'] > row:
                for widget in widgets:
                    widget.grid_configure(row=widget.grid_info()['row']-1)
                    
        if add_custom_adduct_button.winfo_ismapped():
            add_custom_adduct_button.grid_configure(row=add_custom_adduct_button.grid_info()['row']-1)
        else:
            add_custom_adduct_button.grid()
            add_custom_adduct_button.grid_configure(row=21)
        
        for index, adduct_list in enumerate(local_custom_adducts_list):
            if adduct_list[0] == adduct:
                del local_custom_adducts_list[index]
                break
        
    # Create a new top-level window
    set_parameters_window = tk.Toplevel()
    #set_parameters_window.attributes("-topmost", True)
    set_parameters_window.withdraw()
    set_parameters_window.title("Set Parameters")
    icon = ImageTk.PhotoImage(ico_image)
    set_parameters_window.iconphoto(False, icon)
    set_parameters_window.resizable(False, False)
    set_parameters_window.grab_set()
    
    # Add Ttk Styles
    library_building_frame_style = ttk.Style().configure("library_building.TLabelframe", font=("Segoe UI", list_font_size))
    
    small_button_spw_style1 = ttk.Style().configure("small_button_spw_style1.TButton", font=("Segoe UI", list_font_size), relief="raised", padding = (0, 0), justify="center")
    
    small_button_spw_style2 = ttk.Style().configure("small_button_spw_style2.TButton", font=("Segoe UI", list_font_size), relief="raised", padding = (10, 0), justify="center")
    
    # Add widgets to the new window
    library_building_frame = ttk.Labelframe(set_parameters_window, text="Library Building", style="library_building.TLabelframe")
    library_building_frame.grid(row=0, column=0, padx=(10, 5), pady=(0, 0), sticky="nsew")
    
    analysis_frame = ttk.Labelframe(set_parameters_window, text="Analysis", style="library_building.TLabelframe")
    analysis_frame.grid(row=0, column=1, padx=(5, 10), pady=(0, 0), sticky="new")
    
    set_working_dir_button = ttk.Button(set_parameters_window, text="Set Working Directory", style="small_button_spw_style2.TButton", command=select_working_dir_button)
    set_working_dir_button.grid(row=1, column=0, padx=(10,10), pady=(15,15), sticky="wns")
    
    working_dir_label = ttk.Label(set_parameters_window, text="", font=("Segoe UI", list_font_size), wraplength=300)
    working_dir_label.grid(row=1, column=0, columnspan=2, padx=(180, 10), pady=(5, 10), sticky="w")
    working_dir_label.config(text=save_path)
    
    save_param_button = ttk.Button(set_parameters_window, text="    Save Parameters to File    ", style="small_button_spw_style1.TButton", command=save_param_to_file)
    save_param_button.grid(row=0, column=1, padx=(10,10), pady=(15,45), sticky="se")
    
    load_param_button = ttk.Button(set_parameters_window, text="  Load Parameters from File  ", style="small_button_spw_style1.TButton", command=load_param_from_file)
    load_param_button.grid(row=0, column=1, padx=(10, 10), pady=(15,15), sticky="se")
    
    ok_button = ttk.Button(set_parameters_window, text="Ok", style="small_button_spw_style1.TButton", command=ok_sp_window)
    ok_button.grid(row=1, column=1, padx=(10, 100), pady=(15,15), sticky="nse")
    
    cancel_button = ttk.Button(set_parameters_window, text="Cancel", style="small_button_spw_style1.TButton", command=close_sp_window)
    cancel_button.grid(row=1, column=1, padx=(10,10), pady=(15,15), sticky="nse")
    
    # Widgets to library_building_frame
    use_custom_library_checkbox_state = tk.BooleanVar(value=custom_glycans_list[0])
    use_custom_library_checkbox = ttk.Checkbutton(library_building_frame, text="Use Custom Library", variable=use_custom_library_checkbox_state, command=use_custom_library_checkbox_state_check)
    use_custom_library_checkbox.grid(row=0, column=0, padx=10, sticky="nw")
    ToolTip(use_custom_library_checkbox, "Allows to use a list containing specific glycans.")
    
    from_file_button_state = tk.DISABLED
    if use_custom_library_checkbox_state.get():
        from_file_button_state = tk.NORMAL
    from_file_button = ttk.Button(library_building_frame, text=" Set custom glycans ", style="small_button_spw_style1.TButton", command=open_custom_glycans_window, state=from_file_button_state)
    from_file_button.grid(row=0, column=1, padx=(10,10), sticky="e")
    ToolTip(from_file_button, "Allows to use a list containing specific glycans.")
    
    # Generate library part
    generate_library_label = ttk.Label(library_building_frame, text='Library Generation Settings:', font=("Segoe UI", list_font_size))
    generate_library_label.grid(row=2, column=0, columnspan=2, padx=(10, 10), pady=(10, 0), sticky="w")
    
    total_monosaccharides_label = ttk.Label(library_building_frame, text='Total Monosaccharides:', font=("Segoe UI", list_font_size))
    total_monosaccharides_label.grid(row=3, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_monosaccharides_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_monosaccharides_label.grid(row=4, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_monosaccharides_label, "Insert the minimum number of monosaccharides in your glycans.")
    
    min_monosaccharides_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_monosaccharides_entry_state = tk.DISABLED
    min_monosaccharides_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_monosaccharides_entry.insert(0, min_max_monos[0])
    min_monosaccharides_entry.config(state=min_monosaccharides_entry_state)
    min_monosaccharides_entry.grid(row=4, column=0, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_monosaccharides_entry, "Insert the minimum number of monosaccharides in your glycans.")
    
    max_monosaccharides_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_monosaccharides_label.grid(row=4, column=0, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_monosaccharides_label, "Insert the maximum number of monosaccharides in your glycans.")
    
    max_monosaccharides_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_monosaccharides_entry_state = tk.DISABLED
    max_monosaccharides_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_monosaccharides_entry.insert(0, min_max_monos[1])
    max_monosaccharides_entry.config(state=max_monosaccharides_entry_state)
    max_monosaccharides_entry.grid(row=4, column=0, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_monosaccharides_entry, "Insert the maximum number of monosaccharides in your glycans.")
    
    hexoses_label = ttk.Label(library_building_frame, text='Hex:', font=("Segoe UI", list_font_size))
    hexoses_label.grid(row=3, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_hex_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_hex_label.grid(row=4, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_hex_label, "Insert the minimum number of hexoses in your glycans.")
    
    min_hex_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_hex_entry_state = tk.DISABLED
    min_hex_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_hex_entry.insert(0, min_max_hex[0])
    min_hex_entry.config(state=min_hex_entry_state)
    min_hex_entry.grid(row=4, column=1, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_hex_entry, "Insert the minimum number of hexoses in your glycans.")
    
    max_hex_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_hex_label.grid(row=4, column=1, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_hex_label, "Insert the maximum number of hexoses in your glycans.")
    
    max_hex_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_hex_entry_state = tk.DISABLED
    max_hex_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_hex_entry.insert(0, min_max_hex[1])
    max_hex_entry.config(state=max_hex_entry_state)
    max_hex_entry.grid(row=4, column=1, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_hex_entry, "Insert the maximum number of hexoses in your glycans.")
    
    dhex_label = ttk.Label(library_building_frame, text='dHex:', font=("Segoe UI", list_font_size))
    dhex_label.grid(row=5, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_dhex_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_dhex_label.grid(row=6, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_dhex_label, "Insert the minimum number of deoxyhexoses in your glycans.")
    
    min_dhex_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_dhex_entry_state = tk.DISABLED
    min_dhex_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_dhex_entry.insert(0, min_max_fuc[0])
    min_dhex_entry.config(state=min_dhex_entry_state)
    min_dhex_entry.grid(row=6, column=0, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_dhex_entry, "Insert the minimum number of deoxyhexoses in your glycans.")
    
    max_dhex_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_dhex_label.grid(row=6, column=0, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_dhex_label, "Insert the maximum number of deoxyhexoses in your glycans.")
    
    max_dhex_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_dhex_entry_state = tk.DISABLED
    max_dhex_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_dhex_entry.insert(0, min_max_fuc[1])
    max_dhex_entry.config(state=max_dhex_entry_state)
    max_dhex_entry.grid(row=6, column=0, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_dhex_entry, "Insert the maximum number of deoxyhexoses in your glycans.")
    
    hexnac_label = ttk.Label(library_building_frame, text='HexNAc:', font=("Segoe UI", list_font_size))
    hexnac_label.grid(row=5, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_hexnac_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_hexnac_label.grid(row=6, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_hexnac_label, "Insert the minimum number of N-Acetylhexosamines in your glycans.")
    
    min_hexnac_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_hexnac_entry_state = tk.DISABLED
    min_hexnac_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_hexnac_entry.insert(0, min_max_hexnac[0])
    min_hexnac_entry.config(state=min_hexnac_entry_state)
    min_hexnac_entry.grid(row=6, column=1, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_hexnac_entry, "Insert the minimum number of N-Acetylhexosamines in your glycans.")
    
    max_hexnac_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_hexnac_label.grid(row=6, column=1, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_hexnac_label, "Insert the maximum number of N-Acetylhexosamines in your glycans.")
    
    max_hexnac_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_hexnac_entry_state = tk.DISABLED
    max_hexnac_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_hexnac_entry.insert(0, min_max_hexnac[1])
    max_hexnac_entry.config(state=max_hexnac_entry_state)
    max_hexnac_entry.grid(row=6, column=1, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_hexnac_entry, "Insert the maximum number of N-Acetylhexosamines in your glycans.")
    
    sia_label = ttk.Label(library_building_frame, text='Sialic Acids:', font=("Segoe UI", list_font_size))
    sia_label.grid(row=7, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_sia_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_sia_label.grid(row=8, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_sia_label, "Insert the minimum number of sialic acids in your glycans.")
    
    min_sia_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_sia_entry_state = tk.DISABLED
    min_sia_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_sia_entry.insert(0, min_max_sia[0])
    min_sia_entry.config(state=min_sia_entry_state)
    min_sia_entry.grid(row=8, column=0, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_sia_entry, "Insert the minimum number of sialic acids in your glycans.")
    
    max_sia_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_sia_label.grid(row=8, column=0, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_sia_label, "Insert the maximum number of sialic acids in your glycans.")
    
    max_sia_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_sia_entry_state = tk.DISABLED
    max_sia_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_sia_entry.insert(0, min_max_sia[1])
    max_sia_entry.config(state=max_sia_entry_state)
    max_sia_entry.grid(row=8, column=0, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_sia_entry, "Insert the maximum number of sialic acids in your glycans.")
    
    neu5ac_label = ttk.Label(library_building_frame, text='Neu5Ac:', font=("Segoe UI", list_font_size))
    neu5ac_label.grid(row=7, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_neu5ac_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_neu5ac_label.grid(row=8, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_neu5ac_label, "Insert the minimum number of Neu5Ac in your glycans.")
    
    min_neu5ac_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_neu5ac_entry_state = tk.DISABLED
    min_neu5ac_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_neu5ac_entry.insert(0, min_max_ac[0])
    min_neu5ac_entry.config(state=min_neu5ac_entry_state)
    min_neu5ac_entry.grid(row=8, column=1, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_neu5ac_entry, "Insert the minimum number of Neu5Ac in your glycans.")
    
    max_neu5ac_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_neu5ac_label.grid(row=8, column=1, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_neu5ac_label, "Insert the maximum number of Neu5Ac in your glycans.")
    
    max_neu5ac_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_neu5ac_entry_state = tk.DISABLED
    max_neu5ac_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_neu5ac_entry.insert(0, min_max_ac[1])
    max_neu5ac_entry.config(state=max_neu5ac_entry_state)
    max_neu5ac_entry.grid(row=8, column=1, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_neu5ac_entry, "Insert the maximum number of Neu5Ac in your glycans.")
    
    neu5gc_label = ttk.Label(library_building_frame, text='Neu5Gc:', font=("Segoe UI", list_font_size))
    neu5gc_label.grid(row=9, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_neu5gc_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_neu5gc_label.grid(row=10, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_neu5gc_label, "Insert the minimum number of Neu5Gc in your glycans.")
    
    min_neu5gc_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_neu5gc_entry_state = tk.DISABLED
    min_neu5gc_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_neu5gc_entry.insert(0, min_max_gc[0])
    min_neu5gc_entry.config(state=min_neu5gc_entry_state)
    min_neu5gc_entry.grid(row=10, column=0, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_neu5gc_entry, "Insert the minimum number of Neu5Gc in your glycans.")
    
    max_neu5gc_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_neu5gc_label.grid(row=10, column=0, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_neu5gc_label, "Insert the maximum number of Neu5Gc in your glycans.")
    
    max_neu5gc_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_neu5gc_entry_state = tk.DISABLED
    max_neu5gc_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_neu5gc_entry.insert(0, min_max_gc[1])
    max_neu5gc_entry.config(state=max_neu5gc_entry_state)
    max_neu5gc_entry.grid(row=10, column=0, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_neu5gc_entry, "Insert the maximum number of Neu5Gc in your glycans.")
    
    xyl_label = ttk.Label(library_building_frame, text='Xylose:', font=("Segoe UI", list_font_size))
    xyl_label.grid(row=9, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_xyl_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_xyl_label.grid(row=10, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_xyl_label, "Insert the minimum number of Xylose in your glycans.")
    
    min_xyl_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_xyl_entry_state = tk.DISABLED
    min_xyl_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_xyl_entry.insert(0, min_max_xyl[0])
    min_xyl_entry.config(state=min_xyl_entry_state)
    min_xyl_entry.grid(row=10, column=1, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_xyl_entry, "Insert the minimum number of Xylose in your glycans.")
    
    max_xyl_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_xyl_label.grid(row=10, column=1, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_xyl_label, "Insert the maximum number of Xylose in your glycans.")
    
    max_xyl_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_xyl_entry_state = tk.DISABLED
    max_xyl_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_xyl_entry.insert(0, min_max_xyl[1])
    max_xyl_entry.config(state=max_xyl_entry_state)
    max_xyl_entry.grid(row=10, column=1, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_xyl_entry, "Insert the maximum number of Xylose in your glycans.")
    
    hn_label = ttk.Label(library_building_frame, text='Hexosamines:', font=("Segoe UI", list_font_size))
    hn_label.grid(row=11, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_hn_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_hn_label.grid(row=12, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_hn_label, "Insert the minimum number of Hexosamines in your glycans.")
    
    min_hn_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_hn_entry_state = tk.DISABLED
    min_hn_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_hn_entry.insert(0, min_max_hn[0])
    min_hn_entry.config(state=min_hn_entry_state)
    min_hn_entry.grid(row=12, column=0, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_hn_entry, "Insert the minimum number of Hexosamines in your glycans.")
    
    max_hn_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_hn_label.grid(row=12, column=0, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_hn_label, "Insert the maximum number of Hexosamines in your glycans.")
    
    max_hn_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_hn_entry_state = tk.DISABLED
    max_hn_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_hn_entry.insert(0, min_max_hn[1])
    max_hn_entry.config(state=max_hn_entry_state)
    max_hn_entry.grid(row=12, column=0, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_hn_entry, "Insert the maximum number of Hexosamines in your glycans.")
    
    ua_label = ttk.Label(library_building_frame, text='Uronic Acids:', font=("Segoe UI", list_font_size))
    ua_label.grid(row=11, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    
    min_ua_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", list_font_size))
    min_ua_label.grid(row=12, column=1, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(min_ua_label, "Insert the minimum number of Uronic Acids in your glycans.")
    
    min_ua_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        min_ua_entry_state = tk.DISABLED
    min_ua_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    min_ua_entry.insert(0, min_max_ua[0])
    min_ua_entry.config(state=min_ua_entry_state)
    min_ua_entry.grid(row=12, column=1, padx=(40, 0), pady=0, sticky='w')
    ToolTip(min_ua_entry, "Insert the minimum number of Uronic Acids in your glycans.")
    
    max_ua_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", list_font_size))
    max_ua_label.grid(row=12, column=1, padx=(85, 10), pady=(0, 0), sticky="w")
    ToolTip(max_ua_label, "Insert the maximum number of Uronic Acids in your glycans.")
    
    max_ua_entry_state = tk.NORMAL
    if use_custom_library_checkbox_state.get():
        max_ua_entry_state = tk.DISABLED
    max_ua_entry = ttk.Spinbox(library_building_frame, width=4, state=tk.NORMAL, from_=0, to=99)
    max_ua_entry.insert(0, min_max_ua[1])
    max_ua_entry.config(state=max_ua_entry_state)
    max_ua_entry.grid(row=12, column=1, padx=(125, 0), pady=0, sticky='w')
    ToolTip(max_ua_entry, "Insert the maximum number of Uronic Acids in your glycans.")
    
    intstandard_label = ttk.Label(library_building_frame, text='Internal Standard:', font=("Segoe UI", list_font_size))
    intstandard_label.grid(row=13, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(intstandard_label, "Insert the internal standard used. You can input a mass, a chemical formula or a glycan formula. If you use a glycan formula, add 'glycan-' before the formula and the reducing end modification (tag or reduced) or permethylation will be applied to it, if you selected these modifications.\nIn case of Amidation/ Ethyl-Esterification, input the modified sialic acid (e.g.: Am, E, AmG or EG).")
    
    intstandard_entry = ttk.Entry(library_building_frame, width=5, state=tk.NORMAL)
    intstandard_entry.insert(0, internal_standard)
    intstandard_entry.grid(row=14, column=0, padx=(10, 10), pady=0, sticky='we')
    ToolTip(intstandard_entry, "Insert the internal standard used. You can input a mass, a chemical formula or a glycan formula. If you use a glycan formula, add 'glycan-' before the formula and the reducing end modification (tag or reduced) or permethylation will be applied to it, if you selected these modifications.\nIn case of Amidation/ Ethyl-Esterification, input the modified sialic acid (e.g.: Am, E, AmG or EG).")
    
    custom_monosaccharide_button = ttk.Button(library_building_frame, text="Custom Monosaccharides", style="small_button_spw_style1.TButton", command=open_custom_monosaccharides_window)
    custom_monosaccharide_button.grid(row=13, column=1, padx=(10, 10), pady=(5, 0), sticky='we')
    ToolTip(custom_monosaccharide_button, "Opens a window that let's you input custom monosaccharides information, allowing you to build search spaces with different derivatizations or exotic monosaccharides.")
    
    standard_lib_param_button = ttk.Button(library_building_frame, text="Standard Library\nParameters", style="small_button_spw_style1.TButton", command=create_standard_library_param_window)
    standard_lib_param_button.grid(row=14, rowspan=2, column=1, padx=(10, 10), pady=0, sticky='we')
    ToolTip(standard_lib_param_button, "Provides standard parameters to be used as a starting point for your data analysis.\nWarning: As these default parameters are based of mammalian samples, you may need to adjust them to better fit samples from other organisms.")
    
    modifications_label = ttk.Label(library_building_frame, text='Modifications:', font=("Segoe UI", list_font_size))
    modifications_label.grid(row=15, column=0, columnspan=2, padx=(10, 10), pady=(10, 0), sticky="w")
    
    reducing_end_tag_checkbox_activation_state = tk.NORMAL
    if reduced:
        reducing_end_tag_checkbox_activation_state = tk.DISABLED
    reducing_end_tag_checkbox_state = tk.BooleanVar(value=reducing_end_boolean)
    reducing_end_tag_checkbox = ttk.Checkbutton(library_building_frame, text="Reducing End Tag", variable=reducing_end_tag_checkbox_state, command=reducing_end_tag_checkbox_state_check, state=reducing_end_tag_checkbox_activation_state)
    reducing_end_tag_checkbox.grid(row=16, column=0, padx=10, sticky="nw")
    ToolTip(reducing_end_tag_checkbox, "Select if you added a reducing end tag to your samples' glycans or if you're analyzing a glycopeptide.")
    
    reducing_end_tag_label = ttk.Label(library_building_frame, text='Tag:', font=("Segoe UI", list_font_size))
    reducing_end_tag_label.grid(row=17, column=0, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(reducing_end_tag_label, "Select the reducing end tag or type in the ADDED mass of the reducing end tag (ie. 133.0644 for Girard Reagent P (GirP) tag or 219.1735 for Procainamide (ProA) tag), or type in the chemical formula of the reducing end tag (ie. C7H7N3 for GirP or C13H21N3 for ProA) or type in 'pep-' followed by a peptide sequence for analyzing a glycopeptide (ie. pep-NK for the dipeptide made out of an asparagine and a lysine residue).")
    
    reducing_end_tag_options = list(reducing_end_tags.keys())+['Peptide', 'Custom']
    reducing_end_tag_dropdown = ttk.Combobox(library_building_frame, state="readonly", values=reducing_end_tag_options, width=7)
    if reducing_end_tag == 0.0 or reducing_end_tag == '0.0':
        reducing_end_tag_dropdown.set('Custom')
        reducing_end_tag_dropdown.state(['disabled'])
        reducing_end_tag_checkbox_state.set(False)
    elif reducing_end_tag in reducing_end_tags.values():
        swapped_reducing_end_tags = {v: k for k, v in reducing_end_tags.items()}
        reducing_end_tag_dropdown.set(swapped_reducing_end_tags[reducing_end_tag])
        reducing_end_tag_checkbox_state.set(True)
    elif reducing_end_tag.startswith('pep-'):
        reducing_end_tag_dropdown.set('Peptide')
        reducing_end_tag_checkbox_state.set(True)
    else:
        reducing_end_tag_dropdown.set('Custom')
        reducing_end_tag_checkbox_state.set(True)
    reducing_end_tag_dropdown.bind("<<ComboboxSelected>>", on_tag_selected)
    reducing_end_tag_dropdown.grid(row=17, column=0, padx=(40, 0), pady=0, sticky='w')
    ToolTip(reducing_end_tag_dropdown, "Select the reducing end tag or type in the ADDED mass of the reducing end tag (ie. 133.0644 for Girard Reagent P (GirP) tag or 219.1735 for Procainamide (ProA) tag), or type in the chemical formula of the reducing end tag (ie. C7H7N3 for GirP or C13H21N3 for ProA) or type in 'pep-' followed by a peptide sequence for analyzing a glycopeptide (ie. pep-NK for the dipeptide made out of an asparagine and a lysine residue).")
    
    reducing_end_tag_entry_state = tk.DISABLED
    if reducing_end_tag_checkbox_state.get() and reducing_end_tag_dropdown.get() not in reducing_end_tags:
        reducing_end_tag_entry_state = tk.NORMAL
    reducing_end_tag_entry = ttk.Entry(library_building_frame, width=14, state = tk.NORMAL)
    reducing_end_tag_entry.insert(0, f"{reducing_end_tag.split('-')[-1] if type(reducing_end_tag) == str and reducing_end_tag.startswith('pep-') else reducing_end_tag}") 
    reducing_end_tag_entry.config(state = reducing_end_tag_entry_state)
    reducing_end_tag_entry.grid(row=17, column=0, columnspan=2, padx=(110, 10), pady=0, sticky='w')
    ToolTip(reducing_end_tag_entry, "Select the reducing end tag or type in the ADDED mass of the reducing end tag (ie. 133.0644 for Girard Reagent P (GirP) tag or 219.1735 for Procainamide (ProA) tag), or type in the chemical formula of the reducing end tag (ie. C7H7N3 for GirP or C13H21N3 for ProA) or type in 'pep-' followed by a peptide sequence for analyzing a glycopeptide (ie. pep-NK for the dipeptide made out of an asparagine and a lysine residue).")
    
    permethylated_checkbox_state = tk.BooleanVar(value=permethylated)
    permethylated_checkbox = ttk.Checkbutton(library_building_frame, text="Permethylated", variable=permethylated_checkbox_state, command=permethylated_checkbox_state_check)
    permethylated_checkbox.grid(row=18, column=0, padx=10, sticky="w")
    ToolTip(permethylated_checkbox, "Select if the sample was permethylated.")
    
    reduced_checkbox_activation_state = tk.NORMAL
    if reducing_end_tag_checkbox_state.get():
        reduced_checkbox_activation_state = tk.DISABLED
    reduced_checkbox_state = tk.BooleanVar(value=reduced)
    reduced_checkbox = ttk.Checkbutton(library_building_frame, text="Reduced End", variable=reduced_checkbox_state, command=reduced_checkbox_state_check, state = reduced_checkbox_activation_state)
    reduced_checkbox.grid(row=19, column=0, padx=10, sticky="w")
    ToolTip(reduced_checkbox, "Select if the sample doesn't have a reducing end tag and had it's reducing end reduced.")
    
    lac_ee_checkbox_state = tk.BooleanVar(value=lactonized_ethyl_esterified)
    lac_ee_checkbox = ttk.Checkbutton(library_building_frame, text="Amidated/Ethyl-Esterified", variable=lac_ee_checkbox_state, command=lac_ee_checkbox_state_check)
    lac_ee_checkbox.grid(row=20, column=0, padx=10, sticky="w")
    ToolTip(lac_ee_checkbox, "Select if the sample has been amidated and ethyl-esterified, allowing to distinguish between alpha-2,3 and alpha-2,6 N-Acetilneuraminic Acids.")
    
    sulfation_label = ttk.Label(library_building_frame, text='Sulfations:', font=("Segoe UI", list_font_size))
    sulfation_label.grid(row=21, column=0, padx=10, sticky="w")
    ToolTip(sulfation_label, "Input the minimum and maximum number of sulfations per glycan.")
    
    min_sulfation_entry = ttk.Spinbox(library_building_frame, width=3, from_=0, to=99)
    min_sulfation_entry.grid(row=21, column=0, padx=(115, 0), sticky='w')
    min_sulfation_entry.insert(0, min_max_sulfation[0])
    ToolTip(min_sulfation_entry, "Input the minimum and maximum number of sulfations per glycan.")
    
    sulfation_hyphen_label = ttk.Label(library_building_frame, text='-', font=("Segoe UI", list_font_size))
    sulfation_hyphen_label.grid(row=21, column=0, padx=(150, 0), sticky="w")
    ToolTip(sulfation_hyphen_label, "Input the minimum and maximum number of sulfations per glycan.")
    
    max_sulfation_entry = ttk.Spinbox(library_building_frame, width=3, from_=0, to=99)
    max_sulfation_entry.grid(row=21, column=0, padx=(160, 0), sticky='w')
    max_sulfation_entry.insert(0, min_max_sulfation[1])
    ToolTip(max_sulfation_entry, "Input the minimum and maximum number of sulfations per glycan.")
    
    phosphorylation_label = ttk.Label(library_building_frame, text='Phosphorylations:', font=("Segoe UI", list_font_size))
    phosphorylation_label.grid(row=22, column=0, padx=10, sticky="w")
    ToolTip(phosphorylation_label, "Input the minimum and maximum number of phosphorylations per glycan.")
    
    min_phosphorylation_entry = ttk.Spinbox(library_building_frame, width=3, from_=0, to=99)
    min_phosphorylation_entry.grid(row=22, column=0, padx=(115, 0), sticky='w')
    min_phosphorylation_entry.insert(0, min_max_phosphorylation[0])
    ToolTip(min_phosphorylation_entry, "Input the minimum and maximum number of phosphorylations per glycan.")
    
    phosphorylation_hyphen_label = ttk.Label(library_building_frame, text='-', font=("Segoe UI", list_font_size))
    phosphorylation_hyphen_label.grid(row=22, column=0, padx=(150, 0), sticky="w")
    ToolTip(phosphorylation_hyphen_label, "Input the minimum and maximum number of phosphorylations per glycan.")
    
    max_phosphorylation_entry = ttk.Spinbox(library_building_frame, width=3, from_=0, to=99)
    max_phosphorylation_entry.grid(row=22, column=0, padx=(160, 0), sticky='w')
    max_phosphorylation_entry.insert(0, min_max_phosphorylation[1])
    ToolTip(max_phosphorylation_entry, "Input the minimum and maximum number of phosphorylations per glycan.")
    
    ion_mode_label = ttk.Label(library_building_frame, text='Ion mode:', font=("Segoe UI", list_font_size, 'bold'))
    ion_mode_label.grid(row=23, column=0, columnspan=2, padx=10, pady=(0, 10), sticky="w")
    ToolTip(ion_mode_label, "Select the ionization mode used during the data acquisition.")
    
    ion_mode_choice = tk.StringVar(value=None if not ion_mode_selected else ('neg' if max_charges<0 else "pos"))
    
    radio_button_style = ttk.Style().configure("Bold.TRadiobutton", font=("Segoe UI", list_font_size, 'bold'))

    positive_mode_radio = ttk.Radiobutton(library_building_frame, text="Positive", value='pos', variable=ion_mode_choice, style="Bold.TRadiobutton")
    positive_mode_radio.grid(row=23, column=0, columnspan=2, padx=(200, 0), pady=(0, 10), sticky="w")
    ToolTip(positive_mode_radio, "Select this option if data was acquired in positive mode.")
    
    negative_mode_radio = ttk.Radiobutton(library_building_frame, text="Negative", value='neg', variable=ion_mode_choice, style="Bold.TRadiobutton")
    negative_mode_radio.grid(row=23, column=0, columnspan=2, padx=(280, 0), pady=(0, 10), sticky="w")
    ToolTip(negative_mode_radio, "Select this option if data was acquired in negative mode.")
    
    add_settings_label = ttk.Label(library_building_frame, text='Additional Settings:', font=("Segoe UI", list_font_size))
    add_settings_label.grid(row=24, column=0, columnspan=2, padx=(10, 10), pady=(0, 0), sticky="w")
    
    forced_class_label = ttk.Label(library_building_frame, text='Force Glycan Class:', font=("Segoe UI", list_font_size))
    forced_class_label.grid(row=25, column=0, columnspan=2, padx=10, sticky="w")
    ToolTip(forced_class_label, "Forces the compositions to match a known glycan class biological features, reducing the search space, allowing for faster analysis and reducing false positives. Details on the constraints imposed can be found on the User Guide.")
    
    forced_class_options = list(forced_classes.keys())
    forced_class_dropdown = ttk.Combobox(library_building_frame, state="readonly", values=forced_class_options, width=10)
    swapped_forced_classes = {v: k for k, v in forced_classes.items()}
    forced_class_dropdown.set(swapped_forced_classes[forced])
    forced_class_dropdown.grid(row=25, column=0, columnspan=2, padx=(140, 10), sticky="w")
    forced_class_dropdown.bind("<<ComboboxSelected>>", on_forced_class_selected)
    ToolTip(forced_class_dropdown, "Forces the compositions to match a known glycan class biological features, reducing the search space, allowing for faster analysis and reducing false positives. Details on the constraints imposed can be found on the User Guide.")
    
    fast_iso_checkbox_state = tk.BooleanVar(value=fast_iso)
    fast_iso_checkbox = ttk.Checkbutton(library_building_frame, text="Fast Isotopic Distribution Calculation", variable=fast_iso_checkbox_state, command=fast_iso_checkbox_state_check)
    fast_iso_checkbox.grid(row=26, column=0, columnspan=2, padx=10, sticky="w")
    ToolTip(fast_iso_checkbox, "If enabled, calculates the isotopic distribution only based on carbon atoms and corrects it to improve accuracy. If disabled, all atoms are considered in the isotopic distribution calculation and the library building process will take SIGNIFICANTLY longer.")
    
    hires_iso_checkbox_active_state = tk.NORMAL
    if fast_iso_checkbox_state.get():
        hires_iso_checkbox_active_state = tk.DISABLED
    hires_iso_checkbox_state = tk.BooleanVar(value=high_res)
    hires_iso_checkbox = ttk.Checkbutton(library_building_frame, text="High Resolution Distribution Calculation", variable=hires_iso_checkbox_state, command=hires_iso_checkbox_state_check, state=hires_iso_checkbox_active_state)
    hires_iso_checkbox.grid(row=27, column=0, columnspan=2, padx=10, pady=(0, 10), sticky="w")
    ToolTip(hires_iso_checkbox, "Only available if fast isotopic distribution calculation is disabled. Calculates the isotopic distribution peaks in very high resolution. Useful for data acquired in very high resolution MS instruments, such as FT analyzers.")
    
    adducts_label = ttk.Label(library_building_frame, text='Adducts:', font=("Segoe UI", list_font_size))
    adducts_label.grid(row=16, column=1, padx=(25, 10), pady=(0, 0), sticky="w")
    
    adducts_min_label = ttk.Label(library_building_frame, text='Min:', font=("Segoe UI", 9))
    adducts_min_label.grid(row=17, column=1, padx=(35, 10), pady=(0, 0), sticky="w")
    
    adducts_max_label = ttk.Label(library_building_frame, text='Max:', font=("Segoe UI", 9))
    adducts_max_label.grid(row=17, column=1, padx=(80, 10), pady=(0, 0), sticky="w")
    
    hydrogen_min_entry = ttk.Spinbox(library_building_frame, width=4, from_=0, to=99)
    hydrogen_min_entry.grid(row=18, column=1, padx=(25, 0), pady=(0, 0), sticky='w')
    hydrogen_min_entry.insert(0, min_max_proton[0])
    ToolTip(hydrogen_min_entry, "Set the minimum number of proton adducts.")
    
    hydrogen_max_entry = ttk.Spinbox(library_building_frame, width=4, from_=0, to=99)
    hydrogen_max_entry.grid(row=18, column=1, padx=(75, 0), pady=(0, 0), sticky='w')
    hydrogen_max_entry.insert(0, min_max_proton[1])
    ToolTip(hydrogen_max_entry, "Set the maximum number of proton adducts.")
    
    hydrogen_label = ttk.Label(library_building_frame, text='H', font=("Segoe UI", list_font_size))
    hydrogen_label.grid(row=18, column=1, padx=(125, 10), pady=(0, 0), sticky="w")
    ToolTip(hydrogen_label, "Set the range of proton adducts.")
    
    global custom_adducts_widgets
    custom_adducts_widgets = {}
    for index, custom_adduct in enumerate(custom_adducts):
        custom_adducts_widgets[custom_adduct[0]] = []
        custom_adducts_widgets[custom_adduct[0]].append(ttk.Spinbox(library_building_frame, width=4, from_=0, to=99))
        custom_adducts_widgets[custom_adduct[0]][-1].grid(row=19+index, column=1, padx=(25, 0), pady=(0, 0), sticky='w')
        custom_adducts_widgets[custom_adduct[0]][-1].insert(0, custom_adduct[1])
        ToolTip(custom_adducts_widgets[custom_adduct[0]][-1], "Set the minimum number of this adduct.")
        
        custom_adducts_widgets[custom_adduct[0]].append(ttk.Spinbox(library_building_frame, width=4, from_=0, to=99))
        custom_adducts_widgets[custom_adduct[0]][-1].grid(row=19+index, column=1, padx=(75, 0), pady=(0, 0), sticky='w')
        custom_adducts_widgets[custom_adduct[0]][-1].insert(0, custom_adduct[2])
        ToolTip(custom_adducts_widgets[custom_adduct[0]][-1], "Set the maximum number of this adduct.")
        
        custom_adducts_widgets[custom_adduct[0]].append(ttk.Label(library_building_frame, text=custom_adduct[0], font=("Segoe UI", list_font_size)))
        custom_adducts_widgets[custom_adduct[0]][-1].grid(row=19+index, column=1, padx=(125, 10), pady=(0, 0), sticky="w")
        
        custom_adducts_widgets[custom_adduct[0]].append(tk.Button(library_building_frame, text='x', command=lambda adduct=custom_adduct[0]: remove_custom_adduct(adduct), relief='flat', padx=0, pady=0))
        custom_adducts_widgets[custom_adduct[0]][-1].grid(row=19+index, column=1, padx=(10, 0), pady=(0, 0), sticky='w')
            
    if len(custom_adducts) < 3:
        add_custom_adduct_button = ttk.Button(library_building_frame, text="Add Adduct", style="small_button_spw_style1.TButton", command=create_custom_adduct_window)
        add_custom_adduct_button.grid(row=19+len(custom_adducts), column=1, padx=(10, 10), pady=(0, 0), sticky="we")
        ToolTip(add_custom_adduct_button, "Add a custom adduct to be included in the library.")
    
    max_charges_label = ttk.Label(library_building_frame, text='Max Charges:', font=("Segoe UI", list_font_size))
    max_charges_label.grid(row=22, column=1, padx=(20, 10), pady=(5, 10), sticky="w")
    ToolTip(max_charges_label, "Set the maximum amount of charges for any given combination of the adducts selected. Adducts combination that exceed this number won't be used for analysis.")
    
    max_charges_entry = ttk.Spinbox(library_building_frame, width=4, from_=0, to=99)
    max_charges_entry.grid(row=22, column=1, padx=(110, 0), pady=(5, 10), sticky='w')
    max_charges_entry.insert(0, abs(max_charges))
    ToolTip(max_charges_entry, "Set the maximum amount of charges for any given combination of the adducts selected. Adducts combination that exceed this number won't be used for analysis.")
    
    lyase_digested_checkbox_activation_state = tk.NORMAL
    if forced_class_dropdown.get() != 'GAGs':
        lyase_digested_checkbox_activation_state = tk.DISABLED
    lyase_digested_checkbox_state = tk.BooleanVar(value=lyase_digested)
    lyase_digested_checkbox = ttk.Checkbutton(library_building_frame, text="Lyase Digested", variable=lyase_digested_checkbox_state, command=lyase_digested_checkbox_state_check, state=lyase_digested_checkbox_activation_state)
    lyase_digested_checkbox.grid(row=25, column=1, padx=(30, 10), sticky="w")
    ToolTip(lyase_digested_checkbox, "Check this if your GAGs were digested using Lyase enzyme.")
    
    # Widgets to analysis_frame
    multithreaded_checkbox_state = tk.BooleanVar(value=multithreaded_analysis)
    multithreaded_checkbox = ttk.Checkbutton(analysis_frame, text="Multithreading", variable=multithreaded_checkbox_state, command=multithreaded_checkbox_state_check)
    multithreaded_checkbox.grid(row=0, column=0, padx=10, sticky="nw")
    ToolTip(multithreaded_checkbox, "If on, allows the use of multiple CPU cores on the analysis, reducing the processing time.")
    
    number_cores_label = ttk.Label(analysis_frame, text='Number of CPU Cores:', font=("Segoe UI", list_font_size))
    number_cores_label.grid(row=0, column=0, columnspan=2, padx=(10, 60), pady=(0, 0), sticky="e")
    ToolTip(number_cores_label, "Set the number of CPU cores used in the analysis, if multithreading is enabled. Will always leave 2 free cores. Value can be 'all' or a integer.")
    
    number_cores_entry_state = tk.DISABLED
    if multithreaded_checkbox_state.get():
        number_cores_entry_state = tk.NORMAL
    number_cores_entry = ttk.Entry(analysis_frame, width=6, state = tk.NORMAL)
    number_cores_entry.insert(0, number_cores) 
    number_cores_entry.config(state = number_cores_entry_state)
    number_cores_entry.grid(row=0, column=0, columnspan=2, padx=(10, 10), pady=0, sticky='e')
    ToolTip(number_cores_entry, "Set the number of CPU cores used in the analysis, if multithreading is enabled. Will always leave 2 free cores. Value can be 'all' or a integer.")
    
    analyze_ms2_checkbox_state = tk.BooleanVar(value=analyze_ms2[0])
    analyze_ms2_checkbox = ttk.Checkbutton(analysis_frame, text="Analyze MS2", variable=analyze_ms2_checkbox_state, command=analyze_ms2_checkbox_state_check)
    analyze_ms2_checkbox.grid(row=1, column=0, padx=10, pady=(10,5), sticky="w")
    ToolTip(analyze_ms2_checkbox, "Allows for automatic annotation of MS2 spectra. This will increase the processing time significantly if the data contains MS2 spectra.")
    
    ms2_analysis_frame = ttk.Labelframe(analysis_frame, text="MS2 Analysis Settings", style="library_building.TLabelframe")
    ms2_analysis_frame.grid(row=2, column=0, columnspan=2, padx=(10, 10), pady=(0, 10), sticky="new")
    
    force_ms2_comp_checkbox_active_state = tk.DISABLED
    if analyze_ms2_checkbox_state.get():
        force_ms2_comp_checkbox_active_state = tk.NORMAL
    force_ms2_comp_checkbox_state = tk.BooleanVar(value=analyze_ms2[1])
    force_ms2_comp_checkbox = ttk.Checkbutton(ms2_analysis_frame, text="Only assign fragments compatible with\nthe precursor composition", variable=force_ms2_comp_checkbox_state, command=force_ms2_comp_checkbox_state_check, state=force_ms2_comp_checkbox_active_state)
    force_ms2_comp_checkbox.grid(row=0, column=0, padx=10, sticky="w")
    ToolTip(force_ms2_comp_checkbox, "If on, limits the assignments of MS2 peaks to match fragments of which the composition is compatible with the precursor glycan. If off, annotates peaks that match any possible fragment, include those which composition are different than the precursor glycan. Enabling this decreases processing time and might reduce the number of false-positive annotations. Disabling this might help with identifying if a glycan really has the composition identified.")
    
    unrestricted_frags_checkbox_active_state = tk.DISABLED
    if analyze_ms2_checkbox_state.get():
        unrestricted_frags_checkbox_active_state = tk.NORMAL
    unrestricted_frags_checkbox_state = tk.BooleanVar(value=analyze_ms2[2])
    unrestricted_frags_checkbox = ttk.Checkbutton(ms2_analysis_frame, text="Look for fragments of glycans not found on MS1", variable=unrestricted_frags_checkbox_state, command=unrestricted_frags_checkbox_state_check, state = unrestricted_frags_checkbox_active_state)
    unrestricted_frags_checkbox.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="w")
    ToolTip(unrestricted_frags_checkbox, "If on, annotates every MS2 spectra of which the precursor mz matches ANY glycan in the library, regardless if it was found in MS1 or not. If off, only annotates MS2 spectra with precursor mz matching glycans found in MS1 analysis.")
    
    acc_unit_label = ttk.Label(analysis_frame, text='Accuracy Unit:', font=("Segoe UI", list_font_size))
    acc_unit_label.grid(row=3, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(acc_unit_label, "Choose the accuracy unit to use during the analysis: Particles-per-Million (PPM), or absolute mz.")
    
    acc_unit_dropdown_options = ['PPM', 'mz']
    acc_unit_dropdown = ttk.Combobox(analysis_frame, state="readonly", values=acc_unit_dropdown_options, width=6)
    if tolerance[0] == 'ppm':
        acc_unit_dropdown.set('PPM')
    elif tolerance[0] == 'mz':
        acc_unit_dropdown.set('mz')
    acc_unit_dropdown.grid(row=3, column=0, padx=(100, 0), pady=(5, 0), sticky='w')
    ToolTip(acc_unit_dropdown, "Choose the accuracy unit to use during the analysis: Particles-per-Million (PPM), or absolute mz.")
    acc_unit_dropdown.bind("<<ComboboxSelected>>", edit_acc_value)
    
    acc_value_label = ttk.Label(analysis_frame, text='Accuracy Value:', font=("Segoe UI", list_font_size))
    acc_value_label.grid(row=3, column=1, padx=(10, 65), pady=(5, 0), sticky="e")
    ToolTip(acc_value_label, "Enter the value for the accuracy: If using PPM, recommended a value > 1; if using mz, recommended a value between 0 and 1.")
    
    acc_value_entry = ttk.Spinbox(analysis_frame, width=6, from_=0, to=99, increment=1)
    acc_value_entry.configure(increment=0.01 if acc_unit_dropdown.get() == 'mz' else 1)
    acc_value_entry.insert(0, tolerance[1])
    acc_value_entry.grid(row=3, column=1, padx=(10, 10), pady=(5, 0), sticky='e')
    ToolTip(acc_value_entry, "Enter the value for the accuracy: If using PPM, recommended a value > 1; if using mz, recommended a value between 0 and 1.")
    
    rt_int_label = ttk.Label(analysis_frame, text='RT/MT range to analyze (min):', font=("Segoe UI", list_font_size))
    rt_int_label.grid(row=4, column=0, columnspan=2, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(rt_int_label, "Crops the chromatogram/electropherogram based on retention/migration time. This should be done, whenever it's possible, since it has a big impact in the analysis speed.")
    
    rt_int_min_entry = ttk.Spinbox(analysis_frame, width=5, from_=0.0, to=999.0, increment=0.1)
    rt_int_min_entry.insert(0, ret_time_interval[0])
    rt_int_min_entry.grid(row=4, column=0, columnspan=2, padx=(0, 66), pady=(5, 0), sticky='e')
    ToolTip(rt_int_min_entry, "Set the beginning of the retention/migration time interval to analyze.")
    
    rt_int_dash_label = ttk.Label(analysis_frame, text='-', font=("Segoe UI", list_font_size))
    rt_int_dash_label.grid(row=4, column=0, columnspan=2, padx=(0, 57), pady=(5, 0), sticky="e")
    ToolTip(rt_int_dash_label, "Set the end of the retention/migration time interval to analyze.")
    
    rt_int_max_entry = ttk.Spinbox(analysis_frame, width=5, from_=0.0, to=999.0, increment=0.1)
    rt_int_max_entry.insert(0, ret_time_interval[1])
    rt_int_max_entry.grid(row=4, column=0, columnspan=2, padx=(0, 10), pady=(5, 0), sticky='e')
    ToolTip(rt_int_max_entry, "Set the end of the retention/migration time interval to analyze.")
    
    custom_ppp_checkbox_state = tk.BooleanVar(value=min_ppp[0])
    custom_ppp_checkbox = ttk.Checkbutton(analysis_frame, text="Custom minimum datapoints\nper peak", variable=custom_ppp_checkbox_state, command=custom_ppp_checkbox_state_check)
    custom_ppp_checkbox.grid(row=5, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(custom_ppp_checkbox, "Imposes a minimum amount of datapoints per peak to consider any given peak viable. If off, the program automatically determines an amount based on the number of datapoints in 0.2 minutes of chromatogram/electropherogram, implying that peaks are at least 0.2 minutes wide.")
    
    custom_ppp_label = ttk.Label(analysis_frame, text='Number of\nDatapoints:', font=("Segoe UI", list_font_size))
    custom_ppp_label.grid(row=5, column=1, columnspan=2, padx=(0, 55), pady=(5, 0), sticky="e")
    ToolTip(custom_ppp_label, "Defines the minimum amount of datapoints to consider a peak viable.")
    
    custom_ppp_entry_state = tk.DISABLED
    if custom_ppp_checkbox_state.get():
        custom_ppp_entry_state = tk.NORMAL
    custom_ppp_entry = ttk.Spinbox(analysis_frame, width=4, state=tk.NORMAL, from_=0, to=99, increment=1)
    custom_ppp_entry.insert(0, min_ppp[1])
    custom_ppp_entry.config(state=custom_ppp_entry_state)
    custom_ppp_entry.grid(row=5, column=0, columnspan=2, padx=(0, 10), pady=(5, 0), sticky='e')
    ToolTip(custom_ppp_entry, "Defines the minimum amount of datapoints to consider a peak viable.")
    
    close_peaks_checkbox_state = tk.BooleanVar(value=close_peaks[0])
    close_peaks_checkbox = ttk.Checkbutton(analysis_frame, text="Limit peaks picked per\nEIC/EIE", variable=close_peaks_checkbox_state, command=close_peaks_checkbox_state_check)
    close_peaks_checkbox.grid(row=6, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(close_peaks_checkbox, "Imposes a limit to the number of peaks picked in each extracted ion chromatogram/electropherogram. Peaks picked will always be selected based on the most intense peak, followed by the peaks closest to it.")
    
    close_peaks_label = ttk.Label(analysis_frame, text='Number of\nPeaks:', font=("Segoe UI", list_font_size))
    close_peaks_label.grid(row=6, column=1, columnspan=2, padx=(0, 55), pady=(5, 0), sticky="e")
    ToolTip(close_peaks_label, "Defines the maximum number of peaks picked.")
    
    close_peaks_entry_state = tk.DISABLED
    if close_peaks_checkbox_state.get():
        close_peaks_entry_state = tk.NORMAL
    close_peaks_entry = ttk.Spinbox(analysis_frame, width=4, state=tk.NORMAL, from_=0, to=99, increment=1)
    close_peaks_entry.insert(0, close_peaks[1])
    close_peaks_entry.config(state=close_peaks_entry_state)
    close_peaks_entry.grid(row=6, column=0, columnspan=2, padx=(0, 10), pady=(5, 0), sticky='e')
    ToolTip(close_peaks_entry, "Defines the maximum number of peaks picked.")
    
    set_parameters_window.update_idletasks()
    set_parameters_window.deiconify()
    window_width = set_parameters_window.winfo_width()
    window_height = set_parameters_window.winfo_height()
    screen_width = set_parameters_window.winfo_screenwidth()
    screen_height = set_parameters_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    set_parameters_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
def save_results_window():
    '''
    '''
    def close_sr_window():
        sr_window.destroy()  
    
    def save_composition_checkbox_state_check():
        state = save_composition_checkbox_state.get()
        
    def n_glycans_class_checkbox_state_check():
        state = n_glycans_class_checkbox_state.get()
        
    def iso_fits_checkbox_state_check():
        state = iso_fits_checkbox_state.get()
        
    def fill_gaps_checkbox_state_check():
        state = fill_gaps_checkbox_state.get()
        
        if state:
            fill_gaps_perc_samples_spinbox.config(state=tk.NORMAL)
            fill_gaps_rt_tol_spinbox.config(state=tk.NORMAL)
            if metaboanalyst_checkbox_state.get():
                fill_gaps_noise_checkbox.config(state=tk.NORMAL)
            else:
                fill_gaps_noise_checkbox.config(state=tk.DISABLED)
        else:
            fill_gaps_perc_samples_spinbox.config(state=tk.DISABLED)
            fill_gaps_rt_tol_spinbox.config(state=tk.DISABLED)
            fill_gaps_noise_checkbox.config(state=tk.DISABLED)
            
    def fill_gaps_noise_checkbox_state_check():
        state = fill_gaps_noise_checkbox_state.get()
        
    def plot_data_checkbox_state_check():
        state = plot_data_checkbox_state.get()
        
    def metaboanalyst_checkbox_state_check():
        state = metaboanalyst_checkbox_state.get()
        if int(state) == 1:
            select_groups_button.config(state=tk.NORMAL)
            if fill_gaps_checkbox_state.get():
                fill_gaps_noise_checkbox.config(state=tk.NORMAL)
            else:
                fill_gaps_noise_checkbox.config(state=tk.DISABLED)
        else:
            select_groups_button.config(state=tk.DISABLED)
            fill_gaps_noise_checkbox.config(state=tk.DISABLED)
        
    def align_chromatograms_sr_checkbox_state_check():
        state = align_chromatograms_sr_checkbox_state.get()
            
    def ok_sr_window():
        global sample_groups
        def close_progress_save_result_window():
            progress_save_result.destroy()
            sys.stdout = original_stdout
            try:
                stop_thread(t)
                kill_concurrent_futures()
                folders = [folder for folder in os.listdir(save_path) if os.path.isdir(os.path.join(save_path, folder))]
                for i in folders:
                    if "Temp" in i:
                        shutil.rmtree(os.path.join(save_path, i))
            except:
                pass
            
        def ok_progress_save_result_window():
            save_results_button_frame.config(bg="lightgreen")
            close_progress_save_result_window()
            
        try:
            float(curve_fit_sr_entry.get())
        except:
            error_window("Invalid value used in the Curve Fitting Score Threshold field.")
            return
        try:
            float(iso_fit_sr_entry.get())
        except:
            error_window("Invalid value used in the Isotopic Fitting Score Threshold field.")
            return
        try:
            float(s_n_sr_entry.get())
        except:
            error_window("Invalid value used in the Signal-to-Noise Ratio Threshold field.")
            return
        try:
            float(ppm_error_min_sr_entry.get())
        except:
            error_window("Invalid value used in the Minimum PPM Error Threshold field.")
            return
        try:
            float(ppm_error_max_sr_entry.get())
        except:
            error_window("Invalid value used in the Maximum PPM Error Threshold field.")
            return
        try:
            float(auc_percentage_sr_entry.get())/100
        except:
            error_window("Invalid value used in the AUC Percentage Threshold field.")
            return
        try:
            int(min_samples_sr_entry.get())
        except:
            error_window("Invalid value used in the Minimum Samples field.")
            return
            
        progress_save_result = tk.Toplevel()
        # progress_save_result.attributes("-topmost", True)
        progress_save_result.withdraw()
        progress_save_result.title("Saving Results")
        icon = ImageTk.PhotoImage(ico_image)
        progress_save_result.iconphoto(False, icon)
        progress_save_result.resizable(False, False)
        progress_save_result.grab_set()
        progress_save_result.protocol("WM_DELETE_WINDOW", on_closing)

        output_text = ScrolledText(progress_save_result, width=50, height=10, wrap=tk.WORD)
        output_text.grid(row=0, column=0, columnspan=2, padx = 10, pady = 10, sticky="new")
        
        global ok_save_result_button
        ok_save_result_button = ttk.Button(progress_save_result, text="Ok", style="small_button_sfw_style1.TButton", command=ok_progress_save_result_window, state=tk.DISABLED)
        ok_save_result_button.grid(row=1, column=0, padx=10, pady=10, sticky="e")
        
        cancel_save_result_button = ttk.Button(progress_save_result, text="Cancel", style="small_button_sfw_style1.TButton", command=close_progress_save_result_window)
        cancel_save_result_button.grid(row=1, column=1, padx=10, pady=10, sticky="w")
        
        progress_save_result.update_idletasks()
        progress_save_result.deiconify()
        window_width = progress_save_result.winfo_width()
        window_height = progress_save_result.winfo_height()
        screen_width = progress_save_result.winfo_screenwidth()
        screen_height = progress_save_result.winfo_screenheight()
        x_position = (screen_width - window_width) // 2
        y_position = (screen_height - window_height) // 2
        progress_save_result.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
        original_stdout = sys.stdout
        sys.stdout = TextRedirector_save_result(output_text)
        
        global reporter_ions
        reporter_ions = reporter_ions_entry.get().split(",")
        for i_i, i in enumerate(reporter_ions):
            reporter_ions[i_i] = i.strip()
            if len(i) == 0:
                reporter_ions = reporter_ions[:i_i]+reporter_ions[i_i+1:]
        
        global min_max_monos, min_max_hex, min_max_hexnac, min_max_fuc, min_max_sia, min_max_ac, min_max_ac, min_max_gc, forced, min_max_proton, max_charges, reducing_end_tag, internal_standard, permethylated, lactonized_ethyl_esterified, reduced, fast_iso, high_res
        
        output_filtered_data_args = {
                                    'curve fitting score threshold': float(curve_fit_sr_entry.get()),
                                    'isotopic fitting score threshold': float(iso_fit_sr_entry.get()),
                                    'signal to noise ratio threshold': float(s_n_sr_entry.get()),
                                    'maximum ppm threshold': (float(ppm_error_min_sr_entry.get()), float(ppm_error_max_sr_entry.get())),
                                    'minimum percentage auc threshold': float(auc_percentage_sr_entry.get())/100,
                                    'reanalysis': True,
                                    'reanalysis_path': reanalysis_path,
                                    'save path': save_path,
                                    'analyze ms2': analyze_ms2[0],
                                    'unrestricted fragments': analyze_ms2[2],
                                    'reporter ions': reporter_ions,
                                    'output abundance table': [metaboanalyst_checkbox_state.get(), ''],
                                    'output composition-separated data': save_composition_checkbox_state.get(),
                                    'align chromatograms': align_chromatograms_sr_checkbox_state.get(),
                                    'glycan class': 'n_glycans' if n_glycans_class_checkbox_state.get() else 'none',
                                    'retention time tolerance': ret_time_interval[2],
                                    'fragmentation retention time tolerance': rt_tolerance_frag,
                                    'output isotopic fitting data': False,
                                    'output chromatogram plotting data': plot_data_checkbox_state.get(),
                                    'multithreaded': multithreaded_analysis,
                                    'number of cpu cores': number_cores,
                                    'analysis done time': 0.0,
                                    'minimum percentage of samples': int(min_samples_sr_entry.get()),
                                    'sample groups': sample_groups,
                                    'temporary folder': None,
                                    'fill data gaps': (fill_gaps_checkbox_state.get(), float(fill_gaps_perc_samples_spinbox.get()), float(fill_gaps_rt_tol_spinbox.get()), fill_gaps_noise_checkbox_state.get()),
                                    'from GUI': True
                                    }
        
        imp_exp_gen_library_args = {
                                    'custom glycans list': custom_glycans_list,
                                    'min/max monosaccharides': min_max_monos,
                                    'min/max hexoses': min_max_hex,
                                    'min/max hexnac': min_max_hexnac,
                                    'min/max xyloses': min_max_xyl,
                                    'min/max sialic acids': min_max_sia,
                                    'min/max fucoses': min_max_fuc,
                                    'min/max acetyl sialic acids': min_max_ac,
                                    'min/max glycolyl sialic acids': min_max_gc,
                                    'min/max hexosamines': min_max_hn,
                                    'min/max uronic acids': min_max_ua,
                                    'glycan class': forced,
                                    'min/max protons': min_max_proton,
                                    'custom adducts': custom_adducts,
                                    'maximum charges': max_charges,
                                    'reducing end tag': reducing_end_tag,
                                    'fast isotopic pattern calculation': fast_iso,
                                    'high resolution isotopic pattern': high_res,
                                    'import/export library': imp_exp_library,
                                    'library path': library_path,
                                    'exported library name': exp_lib_name,
                                    'only generate library': False,
                                    'save path': save_path,
                                    'internal standard': internal_standard,
                                    'permethylated': permethylated,
                                    'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                                    'reducing end reduced': reduced,
                                    'min/max sulfation': min_max_sulfation,
                                    'min/max phosphorylation': min_max_phosphorylation,
                                    'lyase digested': lyase_digested,
                                    'temporary folder': None,
                                    'custom monosaccharides': custom_monosaccharides,
                                    'from GUI': True
                                    }

        list_of_data_args = {
                            'samples list': samples_list
                            }

        index_spectra_from_file_ms1_args = {
                                            'raw data': None,
                                            'ms level': 1,
                                            'multithreaded': multithreaded_analysis,
                                            'number of cpu cores': number_cores
                                            }
                                            
        index_spectra_from_file_ms2_args = {
                                            'raw data': None,
                                            'ms level': 2,
                                            'multithreaded': multithreaded_analysis,
                                            'number of cpu cores': number_cores
                                            }

        analyze_files_args = {
                            'library': None,
                            'library size': None,
                            'raw data': None,
                            'ms1 index': None,
                            'mass tolerance': tolerance,
                            'retention time interval': ret_time_interval,
                            'min isotopologue peaks': min_isotopologue_peaks,
                            'min points per peak': min_ppp,
                            'adduct combos': None,
                            'maximum charges': max_charges,
                            'custom noise level': custom_noise,
                            'only x most intense peaks': close_peaks,
                            'multithreaded': multithreaded_analysis,
                            'number of cpu cores': number_cores,
                            'analysis start time': None,
                            'temporary folder': None,
                            'from GUI': True
                            }

        analyze_ms2_args = {
                            'ms2 index': None,
                            'raw data': None,
                            'analyzed data': None,
                            'retention time interval': ret_time_interval,
                            'mass tolerance': tolerance,
                            'min/max monosaccharides': min_max_monos,
                            'min/max hexoses': min_max_hex,
                            'min/max hexnac': min_max_hexnac,
                            'min/max xyloses': min_max_xyl,
                            'min/max sialic acids': min_max_sia,
                            'min/max fucoses': min_max_fuc,
                            'min/max acetyl sialic acids': min_max_ac,
                            'min/max glycolyl sialic acids': min_max_gc,
                            'min/max hexosamines': min_max_hn,
                            'min/max uronic acids': min_max_ua,
                            'adduct combos': None,
                            'maximum charges': max_charges,
                            'reducing end tag': reducing_end_tag,
                            'glycan class': forced,
                            'permethylated': permethylated,
                            'reducing end reduced': reduced,
                            'lactonized/ethyl-esterified': lactonized_ethyl_esterified,
                            'assign fragments compatible with precursor': analyze_ms2[1],
                            'annotate all glycans in library': analyze_ms2[2],
                            'fragmentation retention time tolerance': ret_time_interval[2],
                            'multithreaded': multithreaded_analysis,
                            'number of cpu cores': number_cores,
                            'library': None,
                            'temporary folder': None,
                            'from GUI': True,
                            'custom monosaccharides': custom_monosaccharides
                            }

        arrange_raw_data_args = {
                                'analyzed data': None,
                                'sample names': samples_names,
                                'ms2 analyzed': analyze_ms2[0],
                                'save path': save_path,
                                'analysis parameters': [],
                                'adduct combos': None,
                                'library': None,
                                'temporary folder': None,
                                'results file name': None,
                                'from GUI': False,
                                'erase files': True,
                                'calculated noise levels': []
                                }
        
        close_sr_window()
        
        t = threading.Thread(target=run_glycogenius, args=([{
            'output_filtered_data_args': output_filtered_data_args, 
            'imp_exp_gen_library_args': imp_exp_gen_library_args, 
            'list_of_data_args': list_of_data_args, 
            'index_spectra_from_file_ms1_args': index_spectra_from_file_ms1_args,
            'index_spectra_from_file_ms2_args': index_spectra_from_file_ms2_args, 
            'analyze_files_args': analyze_files_args, 
            'analyze_ms2_args': analyze_ms2_args, 
            'arrange_raw_data_args':arrange_raw_data_args
            }]))
        t.start()
            
    def set_groups_window():
        
        samples = gg_file.list_samples().values()
        
        sample_grouping_window(samples, sr_window)
        
    sr_window = tk.Toplevel()
    #sr_window.attributes("-topmost", True)
    sr_window.withdraw()
    sr_window.title("Save Results")
    icon = ImageTk.PhotoImage(ico_image)
    sr_window.iconphoto(False, icon)
    sr_window.resizable(False, False)
    sr_window.grab_set()
    
    save_composition_checkbox_state = tk.BooleanVar(value=compositions)
    save_composition_checkbox = ttk.Checkbutton(sr_window, text="Include whole composition information, in addition to peak-separated\ninformation", variable=save_composition_checkbox_state, command=save_composition_checkbox_state_check)
    save_composition_checkbox.grid(row=0, column=0, columnspan=2, padx=(10, 10), pady=(10, 0), sticky="we")
    ToolTip(save_composition_checkbox, "Normally the data will be peak-separated, meaning that the same composition will result in multiple results, one for each peak in the chromatogram/electropherogram. If you enable this option, you'll also have a sheet where all the peaks area under curve value of a given composition are combined.")
    
    n_glycans_class_checkbox_state = tk.BooleanVar(value=True if forced == 'n_glycans' else False)
    n_glycans_class_checkbox = ttk.Checkbutton(sr_window, text="Determine N-Glycans class", variable=n_glycans_class_checkbox_state, command=n_glycans_class_checkbox_state_check)
    n_glycans_class_checkbox.grid(row=1, column=0, columnspan=2, padx=(10, 10), pady=(5, 0), sticky="we")
    ToolTip(n_glycans_class_checkbox, "If used, outputs extra columns informing if a glycan's composition fits any of the N-Glycans classes, such as oligomannose, complex, hybrid or paucimannose. Can be enabled even if 'Force N-Glycans Compositions' is off in parameters' window, but it's meant to be used for N-Glycans analysis.")
    
    align_chromatograms_sr_checkbox_state = tk.BooleanVar(value=align_chromatograms)
    align_chromatograms_sr_checkbox = ttk.Checkbutton(sr_window, text="Align Results and Chromatograms/Electropherograms by Retention/\nMigration Time", variable=align_chromatograms_sr_checkbox_state, command=align_chromatograms_sr_checkbox_state_check)
    align_chromatograms_sr_checkbox.grid(row=2, column=0, columnspan=2, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(align_chromatograms_sr_checkbox, "Aligns the results and the chromatograms/electropherograms based on peaks assignment. DEPRECATED, will still work, but might take VERY LONG to finish.")
    
    select_groups_button = ttk.Button(sr_window, text="Separate Samples into Groups", style="small_button_spw_style1.TButton", command=set_groups_window)
    select_groups_button.grid(row=3, column=0, columnspan=2, padx=(10, 10), pady=(5, 0), sticky="we")
    ToolTip(select_groups_button, "Allows you to separate the samples into groups through an easy-to-use interface. Samples left ungrouped will be excluded from the analysis, unless all samples are ungrouped. This will change the way the % of samples is calculated on certain parameters.")
    
    additional_files_frame = ttk.Labelframe(sr_window, text="Additional Files:", style="library_building.TLabelframe")
    additional_files_frame.grid(row=4, column=0, columnspan=2, padx=(10, 10), pady=(10, 0), sticky="nsew")
    
    metaboanalyst_checkbox_state = tk.BooleanVar(value=plot_metaboanalyst[0])
    metaboanalyst_checkbox = ttk.Checkbutton(additional_files_frame, text="Save Glycan Abundance Table (.csv)", variable=metaboanalyst_checkbox_state, command=metaboanalyst_checkbox_state_check)
    metaboanalyst_checkbox.grid(row=0, column=0, columnspan=2, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(metaboanalyst_checkbox, "Saves a .csv file that contains abundance of the identified glycans in your dataset on all the samples (only those that meet the QC thresholds you set). It is compatible with Metaboanalyst, which allows for automated statistical analysis of the dataset. You need to set at least 2 groups with 3 samples each or Metaboanalyst website won't accept the file.")
    
    plot_data_checkbox_state = tk.BooleanVar(value=output_plot_data)
    plot_data_checkbox = ttk.Checkbutton(additional_files_frame, text="Save Plot Data Table (.xlsx)", variable=plot_data_checkbox_state, command=plot_data_checkbox_state_check)
    plot_data_checkbox.grid(row=1, column=0, columnspan=2, padx=(10, 10), pady=(10, 10), sticky="w")
    ToolTip(plot_data_checkbox, "Saves a .xlsx file that contains the EIC/EIE-plotting data of all the glycans in the analysis, enable this option.")
    
    # iso_fits_checkbox_state = tk.BooleanVar(value=iso_fittings)
    # iso_fits_checkbox = ttk.Checkbutton(additional_files_frame, text="Save Isotopic- and Curve-Fitting Data Table (.xlsx)", variable=iso_fits_checkbox_state, command=iso_fits_checkbox_state_check)
    # iso_fits_checkbox.grid(row=2, column=0, columnspan=2, padx=(10, 10), pady=(0, 10), sticky="w")
    # ToolTip(iso_fits_checkbox, "Creates files containing the data used for isotopic-fitting and curve-fitting information. Useful for diagnosing possible issues with the analysis. NOT RECOMMENDED: outputting this file will take a LONG time, depending on the size of your original data/dataset.")
    
    qcp_sr_frame = ttk.Labelframe(sr_window, text="Thresholds for Saved Files:", style="library_building.TLabelframe")
    qcp_sr_frame.grid(row=5, column=0, columnspan=2, padx=(10, 10), pady=(10, 10), sticky="nsew")
    
    iso_fit_sr_label = ttk.Label(qcp_sr_frame, text='Minimum Isotopic Fitting Score:', font=("Segoe UI", list_font_size))
    iso_fit_sr_label.grid(row=0, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(iso_fit_sr_label, "Insert here the minimum isotopic fitting score for peaks to be saved to the excel file. Values allowed from 0.0 to 1.0.")
    
    iso_fit_sr_entry = ttk.Spinbox(qcp_sr_frame, width=5, from_=0, to=1.0, increment=0.01)
    iso_fit_sr_entry.insert(0, iso_fit_score)
    iso_fit_sr_entry.grid(row=0, column=1, padx=(130, 10), pady=(5, 0), sticky='e')
    ToolTip(iso_fit_sr_entry, "Insert here the minimum isotopic fitting score for peaks to be saved to the excel file. Values allowed from 0.0 to 1.0.")
    
    curve_fit_sr_label = ttk.Label(qcp_sr_frame, text='Minimum Curve Fitting Score:', font=("Segoe UI", list_font_size))
    curve_fit_sr_label.grid(row=1, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(curve_fit_sr_label, "Insert here the minimum curve-fitting score for peaks to be saved to the excel file. Values allowed from 0.0 to 1.0.")
    
    curve_fit_sr_entry = ttk.Spinbox(qcp_sr_frame, width=5, from_=0, to=1.0, increment=0.01)
    curve_fit_sr_entry.insert(0, curve_fit_score)
    curve_fit_sr_entry.grid(row=1, column=1, padx=(130, 10), pady=(5, 0), sticky='e')
    ToolTip(curve_fit_sr_entry, "Insert here the minimum curve-fitting score for peaks to be saved to the excel file. Values allowed from 0.0 to 1.0.")
    
    s_n_sr_label = ttk.Label(qcp_sr_frame, text='Minimum Signal-to-Noise Ratio:', font=("Segoe UI", list_font_size))
    s_n_sr_label.grid(row=2, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(s_n_sr_label, "Insert here the minimum amount of signal-to-noise ratio necessary to save a result to the excel file. Values under 1.0 won't make a difference.")
    
    s_n_sr_entry = ttk.Spinbox(qcp_sr_frame, width=5, from_=0, to=9999, increment=1)
    s_n_sr_entry.insert(0, s_to_n)
    s_n_sr_entry.grid(row=2, column=1, padx=(130, 10), pady=(5, 0), sticky='e')
    ToolTip(s_n_sr_entry, "Insert here the minimum amount of signal-to-noise ratio necessary to save a result to the excel file. Values under 1.0 won't make a difference.")
    
    ppm_error_sr_label = ttk.Label(qcp_sr_frame, text='Min/Max PPM Error:', font=("Segoe UI", list_font_size))
    ppm_error_sr_label.grid(row=3, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(ppm_error_sr_label, "Insert here the PPM error range allowed to save a peak to the result excel file.")
    
    ppm_error_min_sr_entry = ttk.Spinbox(qcp_sr_frame, width=5, from_=-999, to=999, increment=1)
    ppm_error_min_sr_entry.insert(0, max_ppm[0])
    ppm_error_min_sr_entry.grid(row=3, column=1, padx=(5, 70), pady=(5, 0), sticky='e')
    ToolTip(ppm_error_min_sr_entry, "Insert here the minimum PPM error allowed to save a peak to the result excel file.")
    
    ppm_error_hyphen_sr_label = ttk.Label(qcp_sr_frame, text='-', font=("Segoe UI", list_font_size))
    ppm_error_hyphen_sr_label.grid(row=3, column=1, padx=(5, 60), pady=(5, 0), sticky="e")
    ToolTip(ppm_error_hyphen_sr_label, "Insert here the PPM error range allowed to save a peak to the result excel file.")
    
    ppm_error_max_sr_entry = ttk.Spinbox(qcp_sr_frame, width=5, from_=-999, to=999, increment=1)
    ppm_error_max_sr_entry.insert(0, max_ppm[1])
    ppm_error_max_sr_entry.grid(row=3, column=1, padx=(130, 10), pady=(5, 0), sticky='e')
    ToolTip(ppm_error_max_sr_entry, "Insert here the maximum PPM error allowed to save a peak to the result excel file.")
    
    auc_percentage_sr_label = ttk.Label(qcp_sr_frame, text='Minimum AUC % of maximum intensity:', font=("Segoe UI", list_font_size))
    auc_percentage_sr_label.grid(row=4, column=0, padx=(10, 10), pady=(5, 0), sticky="w")
    ToolTip(auc_percentage_sr_label, "Filters the peaks by a % of the Area-Under-Curve (AUC) of the highest intense peak of that composition. If highest intensity peak is 100 AUC and this setting is set to 5, then only peaks with an AUC of 5 or greater will be saved.")
    
    auc_percentage_sr_entry = ttk.Spinbox(qcp_sr_frame, width=5, from_=0, to=100, increment=1)
    auc_percentage_sr_entry.insert(0, percentage_auc*100)
    auc_percentage_sr_entry.grid(row=4, column=1, padx=(130, 10), pady=(5, 0), sticky='e')
    ToolTip(auc_percentage_sr_entry, "Filters the peaks by a % of the Area-Under-Curve (AUC) of the highest intense peak of that composition. If highest intensity peak is 100 AUC and this setting is set to 5, then only peaks with an AUC of 5 or greater will be saved.")
    
    fill_gaps_checkbox_state = tk.BooleanVar(value=fill_gaps[0])
    fill_gaps_checkbox = ttk.Checkbutton(qcp_sr_frame, text="Fill missing values", variable=fill_gaps_checkbox_state, command=fill_gaps_checkbox_state_check)
    fill_gaps_checkbox.grid(row=5, column=0, padx=(10, 10), pady=(0, 5), sticky='w')
    ToolTip(fill_gaps_checkbox, "Fills missing values with peaks that were detected, but didn't meet the quality thresholds set. Will only fill the gap if the peak matches RT/MT of samples where it is found to be good, and it must be good in at least half the samples. Beware: May produce data containing peaks with quality below the thresholds set.")
    
    fill_gaps_perc_samples_label = ttk.Label(qcp_sr_frame, text='% Samples:', font=("Segoe UI", list_font_size))
    fill_gaps_perc_samples_label.grid(row=5, column=0, columnspan=2, padx=(0, 200), pady=(0, 5), sticky='e')
    ToolTip(fill_gaps_perc_samples_label, "The minimum % of samples at which the glycan peak must be found in a good state to fill the gaps. If samples are set into groups, takes the % of samples within each group instead.")
    
    fill_gaps_perc_samples_spinbox = ttk.Spinbox(qcp_sr_frame, width=5, from_=0, to=100, increment=1)
    fill_gaps_perc_samples_spinbox.insert(0, fill_gaps[1])
    fill_gaps_perc_samples_spinbox.config(state=tk.DISABLED)
    fill_gaps_perc_samples_spinbox.grid(row=5, column=0, columnspan=2, padx=(0, 150), pady=(5, 5), sticky='e')
    ToolTip(fill_gaps_perc_samples_spinbox, "The minimum % of samples at which the glycan peak must be found in a good state to fill the gaps. If samples are set into groups, takes the % of samples within each group instead.")
    
    fill_gaps_rt_tol_label = ttk.Label(qcp_sr_frame, text='RT/MT Tol.:', font=("Segoe UI", list_font_size))
    fill_gaps_rt_tol_label.grid(row=5, column=0, columnspan=2, padx=(0, 60), pady=(0, 5), sticky='e')
    ToolTip(fill_gaps_rt_tol_label, "The retention time tolerance for considering two the same elution time when filling gaps.")
    
    fill_gaps_rt_tol_spinbox = ttk.Spinbox(qcp_sr_frame, width=5, from_=0, to=999, increment=0.1)
    fill_gaps_rt_tol_spinbox.insert(0, fill_gaps[2])
    fill_gaps_rt_tol_spinbox.config(state=tk.DISABLED)
    fill_gaps_rt_tol_spinbox.grid(row=5, column=0, columnspan=2, padx=(0, 10), pady=(5, 5), sticky='e')
    ToolTip(fill_gaps_rt_tol_spinbox, "The retention time tolerance for considering two the same elution time when filling gaps.")
    
    min_samples_sr_label = ttk.Label(qcp_sr_frame, text='Glycan in minimum % of samples:', font=("Segoe UI", list_font_size))
    min_samples_sr_label.grid(row=6, column=0, padx=(10, 10), pady=(0, 5), sticky="w")
    ToolTip(min_samples_sr_label, "Filters the glycans by the percentage of samples they are found in (e.g.: if set to 50 and a glycan is only present in 30% of the samples, that glycan will be removed from the results). If samples are set into groups, takes the % of samples within each group instead. This is applied after the gap filling option above and before filling the gaps with noise (below).")
    
    min_samples_sr_entry = ttk.Spinbox(qcp_sr_frame, width=5, from_=0, to=999, increment=1)
    min_samples_sr_entry.insert(0, min_samples)
    min_samples_sr_entry.grid(row=6, column=1, padx=(130, 10), pady=(0, 5), sticky='e')
    ToolTip(min_samples_sr_entry, "Filters the glycans by the percentage of samples they are found in (e.g.: if set to 50 and a glycan is only present in 30% of the samples, that glycan will be removed from the results). If samples are set into groups, takes the % of samples within each group instead. This is applied after the gap filling option above and before filling the gaps with noise (below).")
    
    fill_gaps_noise_checkbox_state = tk.BooleanVar(value=fill_gaps[3])
    fill_gaps_noise_checkbox = ttk.Checkbutton(qcp_sr_frame, text="Fill remaining gaps with local noise level in glycans abundance table", variable=fill_gaps_noise_checkbox_state, command=fill_gaps_noise_checkbox_state_check)
    fill_gaps_noise_checkbox.config(state=tk.DISABLED)
    fill_gaps_noise_checkbox.grid(row=7, column=0, columnspan=2, padx=(10, 10), pady=(0, 10), sticky='w')
    ToolTip(fill_gaps_noise_checkbox, "Fills the remaining gaps with local noise level in Glycan Abundance Table only. Requires GlycoGenius version 1.1.13 or later (GUI version 1.0.12 or later), otherwise will not fill the gaps.")
    
    reporter_ions_label = ttk.Label(qcp_sr_frame, text='Filter by reporter ions:', font=("Segoe UI", list_font_size))
    reporter_ions_label.grid(row=8, column=0, columnspan=2, padx=(10, 10), pady=(0, 0), sticky="w")
    ToolTip(reporter_ions_label, "Fragments inputted here will be used to filter MS2 data on output. Only MS2 spectra containing these ions will be considered. You can input fragments using the glycans formula (ie. H1N1T1, where T refers to the reducing end) or the mz of the fragment. Can be reapplied on reanalysis.")
    
    reporter_ions_entry = ttk.Entry(qcp_sr_frame, width=52)
    reporter_ions_entry.grid(row=9, column=0, columnspan=2, padx=(10, 10), pady=(0, 10), sticky='we')
    ToolTip(reporter_ions_entry, "Fragments inputted here will be used to filter MS2 data on output. Only MS2 spectra containing these ions will be considered. You can input fragments using the glycans formula (ie. H1N1T1, where T refers to the reducing end) or the mz of the fragment. Can be reapplied on reanalysis.")
    
    ok_button = ttk.Button(sr_window, text="Ok", style="small_button_spw_style1.TButton", command=ok_sr_window)
    ok_button.grid(row=6, column=1, padx=(10, 100), pady=(15,15), sticky="nse")
    
    cancel_button = ttk.Button(sr_window, text="Cancel", style="small_button_spw_style1.TButton", command=close_sr_window)
    cancel_button.grid(row=6, column=1, padx=(10,10), pady=(15,15), sticky="nse")
    
    sr_window.update_idletasks()
    sr_window.deiconify()
    window_width = sr_window.winfo_width()
    window_height = sr_window.winfo_height()
    screen_width = sr_window.winfo_screenwidth()
    screen_height = sr_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    sr_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")

if __name__ == "__main__":
    from File_Editor import mzml_window_start
    from GG_Draw import draw_glycan
    
    # Avoid creating multiple windows when creating multiple processes
    multiprocessing.freeze_support()
    
    # Determine global variables
    global splash_screen, date, begin_time, temp_folder, mzml_window_start, draw_glycan
    
    # Load matplotlib backends
    matplotlib.use("TkAgg")
    matplotlib.use("svg")
    
    # Load the current data and time for using for the temp folder
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    
    # Determine temp folder name and create it
    temp_folder = os.path.join(tempfile.gettempdir(), "gg_GUI_"+begin_time)
    os.makedirs(temp_folder, exist_ok=True)
    
    # If splash screen is running, kill it
    try:
        splash_screen.destroy()
    except:
        pass
        
    # Run main window loop
    run_main_window()
    
def main():
    from .File_Editor import mzml_window_start
    from .GG_Draw import draw_glycan
    
    # Avoid creating multiple windows when creating multiple processes
    multiprocessing.freeze_support()
    
    # Determine global variables
    global splash_screen, date, begin_time, temp_folder, mzml_window_start, draw_glycan
    
    # Load matplotlib backends
    matplotlib.use("TkAgg")
    matplotlib.use("svg")
    
    # Load the current data and time for using for the temp folder
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    
    # Determine temp folder name and create it
    temp_folder = os.path.join(tempfile.gettempdir(), "gg_GUI_"+begin_time)
    os.makedirs(temp_folder, exist_ok=True)
    
    # If splash screen is running, kill it
    try:
        splash_screen.destroy()
    except:
        pass
        
    # Run main window loop
    run_main_window()