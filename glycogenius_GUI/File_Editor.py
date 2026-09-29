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

from pyteomics import mzml, mzxml
from lxml import etree
from PIL import Image, ImageTk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import ScalarFormatter
from matplotlib.patches import Rectangle, Ellipse
from reportlab.lib.pagesizes import letter, A4
from reportlab.pdfgen import canvas as pdf_canvas
from reportlab.lib.utils import ImageReader
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
from io import BytesIO
import pastaq as pq
import threading
import matplotlib.pyplot as plt
import tkinter as tk
import pathlib
import base64
import struct
import os
import dill
import datetime
import tempfile
import numpy
import psutil
import webbrowser

standard_calibrants = {
    '': '',
    'Na Formate Pos': [90.976644, 158.964068, 226.951491, 294.938915, 362.926338, 430.913762, 498.901186, 566.888609, 634.876033, 702.863456, 770.85088, 838.838303, 906.825727, 974.81315, 1042.800574, 1110.787998, 1178.775421, 1246.762845, 1314.750268, 1382.737692, 1450.725115, 1518.712539],
    'Na Formate Neg': [112.985626, 180.97305, 248.960474, 316.947897, 384.935321, 452.922744, 520.910168, 588.897591, 656.885015, 724.872438, 792.859862, 860.847285, 928.834709, 996.822133, 1064.809556, 1132.79698, 1200.784403, 1268.771827, 1336.75925, 1404.746674, 1472.734097],
    'Na Acetate Pos': [104.992294, 186.995368, 268.998442, 351.001515, 433.004589, 515.007662, 597.010736, 679.01381, 761.016883, 843.019957, 925.02303, 1007.026104, 1089.029178, 1171.032251, 1253.035325, 1335.038399, 1417.041472, 1499.044546, 1581.047619, 1663.050693],
    'Na Acetate Neg': [141.016927, 223.02, 305.023074, 387.026147, 469.029221, 551.032295, 633.035368, 715.038442, 797.041515, 879.044589, 961.047663, 1043.050736, 1125.05381, 1207.056884, 1289.059957, 1371.063031, 1453.066104, 1535.069178, 1617.072252, 1699.075325, 1336.75925, 1404.746674, 1472.734097],
    'Na TFA Pos': [158.964029, 294.938837, 430.913645, 566.888453, 702.863262, 838.83807, 974.812878, 1110.787686, 1246.762494, 1382.737303, 1518.712111, 1654.686919, 1790.661727, 1926.636535, 2062.611343],
    'Na TFA Neg': [112.985587, 248.960396, 384.935204, 520.910012, 656.88482, 792.859628, 928.834437, 1064.809245, 1200.784053, 1336.758861, 1472.733669, 1608.708477, 1744.683286, 1880.658094, 2016.632902, 2152.60771, 2288.582518, 2424.557327, 2560.532135, 2696.506943, 2832.481751, 2968.456559, 1535.069178, 1617.072252, 1699.075325, 1336.75925, 1404.746674, 1472.734097],
    'Tune Mix ES Pos': [118.086255, 322.048122, 622.02896, 922.009799, 1321.98425, 1521.971476, 2121.933153, 2721.894831],
    'Tune Mix ES Neg': [112.985587, 431.98233, 601.978978, 1033.98811, 1433.962561, 1633.949787, 2233.911465, 2833.873142]
}

global file_path
file_path = ''

global data_storage
data_storage = {}

global plotted_datapoints
plotted_datapoints = [[], []]

global equation
equation = []

global trimming_range
trimming_range = []

global last_selection
last_selection = {}

global rt_bpc
rt_bpc = {}

global stop_run_event
stop_run_event = threading.Event()

# Background task running in the alignment tab ('alignment', 'auto parameters' or None), and whether the window should close when the alignment stops
global background_task, close_after_alignment
background_task = None
close_after_alignment = False

global alignment_settings
alignment_settings = {
    'instrument_type': 'TOF', # Seems to have no impact up to resampling.
    'resolution_ms1': 70000, # BIG impact on m/z values resampling. Could be extrapolated from data. BIG impact on resampling time (higher res = longer to resample).
    'reference_mz': 200, # Reference mz for the specified resolution. Probably has an impact between different equipments, like TOF and ORBITRAP.
    'fwhm_rt': 9, # BIG impact on RT time resampling. Could be extrapolated from data. Small impact on resampling time. Smaller values takes longer. BIG impact on alignment.
    'polarity': '+',
    'num_samples_mz': 5, # The number of sub-samples per mz step in the mesh
    'num_samples_rt': 5, # The number of sub-samples per rt step in the mesh
    'smoothing_coefficient_mz': 0.8, # Smoothing coefficient in the mz direction
    'smoothing_coefficient_rt': 0.4, # Smoothing coefficient in the rt direction
    'slack': 5, # How far the end of each segment can shift during alignment (in virtual “points”)
    'segment_length': 10, # Number of virtual “points” per segment used for alignment
    'num_points': 500, # Total number of virtual “points” in the chromatogram; divided by segment_length gives the number of segments
    'expand_factor': 0.2, # Allows the chromatogram to expand to accommodate the alignment at the beginning and end of the chromatogram
    'peaks_per_window': 20 # Number of most abundant peaks in each segment to use for alignment
}

# Expected type of each numeric alignment setting, used to convert the values typed in the settings window
alignment_settings_types = {
    'resolution_ms1': float,
    'reference_mz': float,
    'fwhm_rt': float,
    'num_samples_mz': int,
    'num_samples_rt': int,
    'smoothing_coefficient_mz': float,
    'smoothing_coefficient_rt': float,
    'slack': int,
    'segment_length': int,
    'num_points': int,
    'expand_factor': float,
    'peaks_per_window': int
}

alignment_settings_descriptions = {
    'instrument_type': 'Type of mass spectrometer used for data acquisition. This can affect the characteristics of the data and may influence the choice of alignment parameters.',
    'resolution_ms1': 'The resolving power of the mass spectrometer at the specified reference m/z for MS1 spectra. Higher resolution can provide more accurate m/z values but may increase processing time.',
    'reference_mz': 'The m/z value at which the specified resolution is measured. This can help in extrapolating the expected peak widths across the m/z range and can influence the resampling of the chromatogram.',
    'fwhm_rt': 'The full width at half maximum (FWHM) of chromatographic peaks in the retention time dimension. This parameter can influence the resampling of the chromatogram and the subsequent alignment.',
    'polarity': 'The ionization mode (positive or negative) used during data acquisition. This can affect the types of ions detected and may influence the alignment process.',
    'num_samples_mz': 'The number of sub-samples to take in the m/z dimension when creating the mesh for alignment. More samples can improve alignment accuracy but will increase processing time.',
    'num_samples_rt': 'The number of sub-samples to take in the retention time dimension when creating the mesh for alignment. More samples can improve alignment accuracy but will increase processing time.',
    'smoothing_coefficient_mz': 'The coefficient used for smoothing the data in the m/z dimension. Higher values can lead to smoother data but may obscure fine details.',
    'smoothing_coefficient_rt': 'The coefficient used for smoothing the data in the retention time dimension. Higher values can lead to smoother data but may obscure fine details.',
    'slack': 'The maximum allowed shift at the end of each segment during alignment, measured in virtual “points”. This parameter can help accommodate variations in retention time across segments.',
    'segment_length': 'The number of virtual “points” in each segment used for alignment. This parameter can influence the granularity of the alignment process.',
    'num_points': 'The total number of virtual “points” in the chromatogram. This determines the resolution of the chromatogram mesh and can affect both alignment accuracy and processing time.',
    'expand_factor': 'The factor by which the chromatogram can expand to accommodate alignment at the beginning and end of the chromatogram. This can help ensure that peaks near the edges are properly aligned.',
    'peaks_per_window': 'The number of most abundant peaks in each segment to use for alignment. This can help focus the alignment on the most relevant features in the data.'
    }

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
                x, y, _, _ = self.widget.bbox("insert")
                x += self.widget.winfo_pointerx()+15
                y += self.widget.winfo_pointery()+15

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

def binary_search_with_tolerance(arr, target, low, high, tolerance, int_arr = [], black_list = []):
    '''A function to quickly find a target in an array by splitting the array recurssively in two and looking for the mid point, finding out if value is bigger or smaller than  target, then splitting again. It also checks if found target in target array is within a tolerance, and picks the most intense one within the tolerance if intensity array is available, else picks the closest one to target.
    
    Parameters
    ----------
    arr : list/np.array
        Target array to search for target.
        
    target : float
        Float to find in target array.
        
    low : int
        Index of the first element in target array.
        
    high : int
        Index of the last element in target array.
        
    tolerance : float
        Tolerance to check for target in target array.
        
    int_arr : list/np.array
        List of intensities, synchronized with arr, to help choose the best target within tolerance.
        
    Uses
    ----
    numpy.argmax : int
        Outputs the index of the highest value in a given array.
        
    numpy.argmin : int
        Outputs the index of the value with smallest difference to a given target, in an array.
        
    Returns
    -------
    selected_id : index
        The index of the selected target.
    '''
    # If array is empty, skip
    if len(arr) == 0:
        return -1
    # Base case: if the range is invalid, the target is not in the array
    if low > high:
        return -1  # Target not found
        
    if target < arr[0] or target > arr[-1]:
        return -1
    
    # Find the middle index
    mid = (low + high) // 2
    
    # Check if the target is within the tolerance range of the middle element
    if abs(arr[mid] - target) <= tolerance:
    
        range_width = 5
        range_search = [mid, mid+1]
        for i in range(mid, mid-range_width, -1):
            if i == 0 or arr[i] < target-tolerance:
                break
            range_search[0] = i
        for i in range(mid+1, mid+range_width+1):
            range_search[1] = i
            if i > high or i >= len(arr)-1 or arr[i] > target+tolerance:
                break
                
        if len(int_arr) != 0:
            array_slice = int_arr[range_search[0]:range_search[1]]
        else:
            array_slice = arr[range_search[0]:range_search[1]]
            
        if len(array_slice) == 0:
            return -1

        if len(int_arr) != 0:
            relative_id = numpy.argmax(array_slice)
        else:
            relative_id = (numpy.abs(array_slice - target).argmin())
        selected_id = range_search[0]+relative_id
            
        # This avoids picking the same peak twice
        forbidden_ids = []
        while arr[selected_id] in black_list:
            forbidden_ids.append(relative_id)
            if len(int_arr) != 0:
                array_slice[relative_id] = 0
            else:
                array_slice[relative_id] += 1000
            if len(int_arr) != 0:
                relative_id = numpy.argmax(array_slice)
            else:
                relative_id = (numpy.abs(array_slice - target).argmin())
            if relative_id in forbidden_ids:
                return -1
            selected_id = range_search[0]+relative_id
            
        return selected_id
    elif arr[mid] < target:
        # If target is greater, ignore the left half
        return binary_search_with_tolerance(arr, target, mid + 1, high, tolerance, int_arr, black_list)
    else:
        # If target is smaller, ignore the right half
        return binary_search_with_tolerance(arr, target, low, mid - 1, tolerance, int_arr, black_list)
        
def calculate_ppm_diff(mz, target):
    '''Calculates the PPM difference between a mz and a target mz.
    
    Parameters
    ----------
    mz : float
        A float of the mz you want to check for the difference.
        
    target : float
        A float of the target mz you're comparing your mz to.
        
    Returns
    -------
    float
        A float of the PPM difference between the mz and target mz.
    '''
    return ((target-mz)/target)*(10**6)

def sample_name_from_path(path):
    '''Returns the sample name of a spectra file path (the file name without its extension), the same way GlycoGenius names samples.'''
    file_name = path.split("/")[-1].split("\\")[-1]
    return ".".join(file_name.split(".")[:-1]) if "." in file_name else file_name

def find_sample_path(sample_name, samples_list):
    '''Finds the path of the spectra file corresponding to a sample name.'''
    for file in samples_list:
        if sample_name_from_path(file) == sample_name:
            return file
    return sample_name

def error_window(text):
    '''This function makes an error window with the text.'''
    messagebox.showerror("Error", text)

def convert_to_base64(float_array):
    '''
    '''
    # Convert the float array to binary data
    binary_data = struct.pack('d' * len(float_array), *float_array)
    
    # Encode the compressed data to base64
    return base64.b64encode(binary_data).decode('utf-8')
    
def make_calibration_report(output_path, file_name, calibrants_list, labels, plot_area, original_file_path = None):
    '''
    '''
    global equation
    
    def change_transparency_to_white_in_memory(input_image_path):
        # Open the image
        img = Image.open(input_image_path).convert("RGBA")  # Ensure it's in RGBA mode

        # Create a new image with a white background
        white_background = Image.new("RGBA", img.size, (255, 255, 255, 255))  # White background

        # Paste the original image onto the white background
        white_background.paste(img, (0, 0), img)  # Use img as a mask to retain visible parts

        # Convert to RGB if needed for further usage
        white_background_rgb = white_background.convert("RGB")  # Change to RGB for non-PNG use

        return white_background_rgb
    
    date = datetime.datetime.now()
    begin_time = str(date)[2:4]+str(date)[5:7]+str(date)[8:10]+"_"+str(date)[11:13]+str(date)[14:16]+str(date)[17:19]
    
    if original_file_path != None:
        original_file_name = original_file_path.split("/")[-1].split("\\")[-1]
    else:
        original_file_name = "_".join(file_name.split("_")[:-1])+"."+file_name.split(".")[-1]
    
    table = [['Calibrant', 'Ref. m/z', 'Found m/z', 'PPM error', 'm/z error', 'Peak Time']]
    for target in calibrants_list.get_children():
        # Get the values from the treeview
        target_values = list(calibrants_list.item(target, "values"))
        
        if target_values[2] != '':
            table.append(target_values)
        
    table = Table(table)
    
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # Header background
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # Header text color
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),  # Center alignment
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),  # Header font
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),  # Row background
        ('GRID', (0, 0), (-1, -1), 1, colors.black),  # Grid lines
    ])
    
    table.setStyle(style)
    
        
    # Start creating a PDF
    c = pdf_canvas.Canvas(os.path.join(output_path, f"{begin_time}_{file_name}_calibration_report.pdf"), pagesize=A4)
    
    # Paper size and default line spacing
    width, height = A4
    space_between_lines = 16
    
    # Load the image
    gg_logo_image = ImageReader(change_transparency_to_white_in_memory(os.path.join(current_dir, "Assets/logo.png")))
    
    # Draw the image on the PDF at coordinates (x, y) and set the width and height
    c.drawImage(gg_logo_image, 50, height-170, width=120, height=120)
    
    # Title of the file
    c.setFont("Helvetica", 20)
    c.drawString(180, height-110, "GlycoGenius Calibration Report")
    
    # Starting line of the text
    text_line = 12
    
    # Document body font
    c.setFont("Helvetica", 12)
    
    c.drawString(50, height-(text_line*space_between_lines), f"Calibration Date and Time: {date}")
    text_line += 1
    
    if len(original_file_name) > 50:
        c.drawString(50, height-(text_line*space_between_lines), f"Calibrated File: {original_file_name[:50]}")
        text_line += 1
        c.drawString(50, height-(text_line*space_between_lines), f"                       {original_file_name[50:]}")
        text_line += 1
    else:
        c.drawString(50, height-(text_line*space_between_lines), f"Calibrated File: {original_file_name}")
        text_line += 1
        
    if len(file_name) > 50:
        c.drawString(50, height-(text_line*space_between_lines), f"New File Name: {file_name[:50]}")
        text_line += 1
        c.drawString(50, height-(text_line*space_between_lines), f"                    {file_name[50:]}")
        text_line += 1
    else:
        c.drawString(50, height-(text_line*space_between_lines), f"New File Name: {file_name}")
        text_line += 1
        
    c.drawString(50, height-(text_line*space_between_lines), f"Calibration Fit: {equation[0]}")
    text_line += 1
    
    c.drawString(50, height-(text_line*space_between_lines), f"Calibration Stats:")
    text_line += 1
    
    text = labels[0].cget("text")
    c.drawString(75, height-(text_line*space_between_lines), f"{text}")
    text_line += 1
    
    text = labels[1].cget("text")
    c.drawString(75, height-(text_line*space_between_lines), f"{text}")
    text_line += 1
    
    text = labels[2].cget("text")
    c.drawString(75, height-(text_line*space_between_lines), f"{text}")
    text_line += 2
    
    c.drawString(50, height-(text_line*space_between_lines), f"Calibration Plot:")
    image_stream = BytesIO()
    plot_area[2].savefig(image_stream, format='png', bbox_inches='tight', transparent=True, dpi=600)
    image_stream.seek(0)
    plot_image = Image.open(image_stream).convert("RGBA")
    plot_img_reader = ImageReader(plot_image)
    c.drawImage(plot_img_reader, (width-360)/2, height-(text_line*space_between_lines)-150, width=360, height=150)
    text_line += 10
    
    c.drawString(50, height-(text_line*space_between_lines), f"Calibrants Table:")
    text_line += 1
    
    table_width, table_height = table.wrapOn(c, 0, 0)
    table.drawOn(c, (width-table_width)/2, height-(text_line*space_between_lines)-table_height)
    
    c.save()

def create_mzml(spectra_file, mode, samples_list, calibration_parameters = [], trimming_parameters = [], alignment_parameters = []):
    '''
    '''
    global equation, file_path
    
    if mode not in ["calibration", "trimming", "aligning"]:
        print("Wrong mode used to create mzml file.")
        return
    
    if mode == 'calibration':
        if len(equation) == 0:
            return
    elif mode == 'trimming':
        if len(trimming_range) == 0:
            return

    if mode == 'aligning':
        spectra_file_path = alignment_parameters[0]
    elif len(samples_list) != 0:
        spectra_file_path = find_sample_path(spectra_file.get(), samples_list)
    else:
        spectra_file_path = file_path
        
    output_file_path = "/".join(spectra_file_path.split("/")[:-1])
    
    if mode == 'calibration':
        new_file_name = sample_name_from_path(spectra_file_path)+"_calibrated.mzML"
    elif mode == 'trimming':
        new_file_name = sample_name_from_path(spectra_file_path)+"_trimmed.mzML"
    elif mode == 'aligning':
        new_file_name = sample_name_from_path(spectra_file_path)+"_aligned.mzML"
        
    new_file_base = new_file_name[:-len(".mzML")]
    counter = 1
    while new_file_name in os.listdir(output_file_path):
        new_file_name = f"{new_file_base}({counter}).mzML"
        counter += 1

    # Read existing mzML file
    if spectra_file_path.split(".")[-1].lower() == 'mzml':
        spectra = mzml.MzML(spectra_file_path)
        file_type = 'mzml'
    elif spectra_file_path.split(".")[-1].lower() == 'mzxml':
        spectra = mzxml.MzXML(spectra_file_path)
        file_type = 'mzxml'

    if file_type == 'mzml' and float(spectra[-1]['scanList']['scan'][0]['scan start time']) > 300:
        time_unit = 'second'
        uo_time_unit = "UO:0000010"
    else:
        time_unit = 'minute'
        uo_time_unit = "UO:0000031"

    # Create the root element for the new mzML file
    mzML = etree.Element("indexedmzML")

    mzML_subelement = etree.SubElement(mzML, 
                                       "mzML")
                                       
    softwareList = etree.SubElement(mzML_subelement,
                                    "softwareList",
                                    count="1")
    software = etree.SubElement(softwareList,
                                "software",
                                id="GlycoGenius_Calibration",
                                version="")
    etree.SubElement(software,
                     'cvParam',
                     cvRef="MS",
                     accession="MS:1000531",
                     name="software",
                     value="")
                                       
    instrumentConfigurationList = etree.SubElement(mzML_subelement,
                                                   "instrumentConfigurationList",
                                                   count="1")
    instrumentConfiguration = etree.SubElement(instrumentConfigurationList,
                                               "instrumentConfiguration",
                                               id="instrument")
    etree.SubElement(instrumentConfiguration,
                     "cvParam",
                     cvRef="MS",
                     accession="MS:1000529",
                     name="instrument serial number",
                     value="")
    componentList = etree.SubElement(instrumentConfiguration,
                                     "componentList",
                                     count="3")
    source = etree.SubElement(componentList,
                              "source",
                              order="1")
    etree.SubElement(source,
                     "cvParam",
                     cvRef="MS",
                     accession="MS:1000008",
                     name="ionization type")
    analyzer = etree.SubElement(componentList,
                                "analyzer",
                                order="2")
    etree.SubElement(analyzer,
                     "cvParam",
                     cvRef="MS",
                     accession="MS:1000443",
                     name="mass analyzer type")
    detector = etree.SubElement(componentList,
                                "detector",
                                order="2")
    etree.SubElement(detector,
                     "cvParam",
                     cvRef="MS",
                     accession="MS:1000253",
                     name="electron multiplier")
                     
    dataProcessingList = etree.SubElement(mzML_subelement,
                                          "dataProcessingList",
                                          count="1")
    dataProcessing = etree.SubElement(dataProcessingList,
                                      "dataProcessing",
                                      id="gg_calibration")
    processingMethod = etree.SubElement(dataProcessing,
                                        "processingMethod",
                                        order="1",
                                        softwareRef="GlycoGenius_Calibration")
    etree.SubElement(processingMethod,
                     "cvParam",
                     cvRef="MS",
                     accession="MS:1000544",
                     name="Conversion to mzML",
                     value="")                 
                                       
    run = etree.SubElement(mzML_subelement,
                           "run",
                           defaultInstrumentConfigurationRef="instrument",
                           id=new_file_name)
                           
    # Create a spectrumList element
    spectrumList = etree.SubElement(run, 
                                    "spectrumList",
                                    count=f"{len(spectra)}",
                                    defaultDataProcessingRef="gg_calibration")
    
    if mode == 'trimming' or mode == 'aligning':
        first_spectrum_num = None

    # Counter of the spectra written to the new file, used to number them sequentially
    written_spectra = 0

    # Loop through each spectrum in the existing file
    for spectrum in spectra:
        if file_type == 'mzxml':
            scan_num = spectrum['num']
            temp_spectrum = {'id':f'scan={scan_num}', 'ms level':spectrum['msLevel'], 'scanList':{'scan': [{'scan start time':spectrum['retentionTime']}]}}
            
            # Polarity (optional in mzXML)
            if spectrum.get('polarity') == '+':
                temp_spectrum['positive scan'] = ''
            elif spectrum.get('polarity') == '-':
                temp_spectrum['negative scan'] = ''
                
            # BPC
            if 'basePeakIntensity' in spectrum.keys():
                temp_spectrum['base peak intensity'] = spectrum['basePeakIntensity']
                
            # TIC
            if 'totIonCurrent' in spectrum.keys():
                temp_spectrum['total ion current'] = spectrum['totIonCurrent']
                
            if spectrum['msLevel'] == 2:
                temp_spectrum['precursorList'] = {}
                temp_spectrum['precursorList']['precursor'] = []
                temp_spectrum['precursorList']['precursor'].append({})
                temp_spectrum['precursorList']['precursor'][0]['selectedIonList'] = {}
                temp_spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'] = []
                temp_spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'].append({})
                temp_spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'] = spectrum['precursorMz'][0]['precursorMz']
                
            temp_spectrum['m/z array'] = spectrum['m/z array']
            temp_spectrum['intensity array'] = spectrum['intensity array']
            
            spectrum = temp_spectrum
        
        if mode == 'trimming':
            if spectrum['scanList']['scan'][0]['scan start time'] > (trimming_range[1] if time_unit == 'minute' else trimming_range[1]*60):
                break
            if first_spectrum_num is None and spectrum['scanList']['scan'][0]['scan start time'] < (trimming_range[0] if time_unit == 'minute' else trimming_range[0]*60):
                continue
            elif first_spectrum_num is None:
                if int(spectrum['ms level']) == 1:
                    first_spectrum_num = int(spectrum['id'].split("=")[-1])
                else:
                    continue
            spectrum['scanList']['scan'][0]['scan start time'] = f"{(float(spectrum['scanList']['scan'][0]['scan start time']) - trimming_range[0]) if time_unit == 'minute' else (float(spectrum['scanList']['scan'][0]['scan start time']) - trimming_range[0]*60)}"
            
        if mode == 'aligning':
            below = False
            segment_found = False
            spectrum_rt_seconds = spectrum['scanList']['scan'][0]['scan start time'] if time_unit == 'second' else spectrum['scanList']['scan'][0]['scan start time']*60
            for segment in alignment_parameters[1]:
                if (spectrum['scanList']['scan'][0]['scan start time'] if time_unit == 'second' else spectrum['scanList']['scan'][0]['scan start time']*60) >= segment[2] and (spectrum['scanList']['scan'][0]['scan start time'] if time_unit == 'second' else spectrum['scanList']['scan'][0]['scan start time']*60) <= segment[3]:
                    if segment[0] < 0 and segment[1] <  0:
                        below = True
                        break

                    spectrum['scanList']['scan'][0]['scan start time'] = (segment[0] if time_unit == 'second' else segment[0]/60) + (spectrum['scanList']['scan'][0]['scan start time'] - (segment[2] if time_unit == 'second' else segment[2]/60)) * ((segment[1]-segment[0])/(segment[3]-segment[2]))
                    segment_found = True

                    break

            # Spectra outside the time map are mapped with the closest segment at the edge, so the retention times stay continuous
            if not below and not segment_found:
                edge_segment = alignment_parameters[1][0] if spectrum_rt_seconds < alignment_parameters[1][0][2] else alignment_parameters[1][-1]
                if edge_segment[3] != edge_segment[2]:
                    spectrum['scanList']['scan'][0]['scan start time'] = (edge_segment[0] if time_unit == 'second' else edge_segment[0]/60) + (spectrum['scanList']['scan'][0]['scan start time'] - (edge_segment[2] if time_unit == 'second' else edge_segment[2]/60)) * ((edge_segment[1]-edge_segment[0])/(edge_segment[3]-edge_segment[2]))

            if below:
                continue
                       
            if first_spectrum_num is None and spectrum['scanList']['scan'][0]['scan start time'] < 0:
                continue
            elif first_spectrum_num is None:
                if int(spectrum['ms level']) == 1:
                    first_spectrum_num = int(spectrum['id'].split("=")[-1])
                else:
                    continue

        # Create a spectrum element, numbered sequentially
        written_spectra += 1
        spectrum['id'] = f"scan={written_spectra}"
        spectrum_id = spectrum['id']
        spectrum_element = etree.SubElement(spectrumList,
                                            "spectrum",
                                            index=f"{written_spectra-1}",
                                            id=f"{spectrum_id}", 
                                            defaultArrayLength=str(len(spectrum['m/z array'])))

        # Add ms level
        etree.SubElement(spectrum_element, 
                         'cvParam',
                         cvRef = 'MS',
                         accession = "MS:1000511",
                         name = 'ms level',
                         value = str(spectrum['ms level']))
                         
        # Optional (for GG) base information
        
        # Scan level name
        if spectrum['ms level'] == 1:
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000579",
                             name="MS1 spectrum",
                             value="")
        else:
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000579",
                             name="MSn spectrum",
                             value="")
        
        # Scan mode name
        if 'positive scan' in spectrum.keys():
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000130",
                             name="positive scan",
                             value="")
        if 'negative scan' in spectrum.keys():
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000130",
                             name="negative scan",
                             value="")
                             
        # Base peak information
        if 'base peak intensity' in spectrum.keys():
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000505",
                             name="base peak intensity",
                             value=str(spectrum['base peak intensity']))
                             
        # TIC                  
        if 'total ion current' in spectrum.keys():
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000285",
                             name="total ion current",
                             value=str(spectrum['total ion current']))
                             
        # Centroided or not
        if 'centroid spectrum' in spectrum.keys():
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000127",
                             name="centroid spectrum",
                             value='')
        elif 'profile spectrum' in spectrum.keys(): # have to check for this one
            None
            
        # Spectrum title
        if 'spectrum title' in spectrum.keys():
            etree.SubElement(spectrum_element, 
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000796",
                             name="spectrum title",
                             value=str(spectrum['spectrum title']))
                         
        # Add retention time
        scanlist = etree.SubElement(spectrum_element, 
                                    'scanList',
                                    count="1")
        etree.SubElement(scanlist,
                         'cvParam',
                         cvRef='MS',
                         accession='MS:1000795',
                         name='no combination',
                         value='')
                         
        scan = etree.SubElement(scanlist, 
                                'scan')
        etree.SubElement(scan,
                         'cvParam',
                         cvRef = 'MS',
                         accession = "MS:1000016",
                         name = 'scan start time',
                         value = str(spectrum['scanList']['scan'][0]['scan start time']),
                         unitCvRef="UO",
                         unitAccession=uo_time_unit,
                         unitName=time_unit)
                         
        if 'scanWindowList' in spectrum['scanList']['scan'][0].keys():
            scanWindowList = etree.SubElement(scan,
                                              'scanWindowList',
                                              count="1")
            scanWindow = etree.SubElement(scanWindowList,
                                          'scanWindow')
            etree.SubElement(scanWindow,
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000501",
                             name="scan window lower limit",
                             value=str(spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']),
                             unitCvRef="MS",
                             unitAccession="MS:1000040",
                             unitName="m/z")
            etree.SubElement(scanWindow,
                             'cvParam',
                             cvRef="MS",
                             accession="MS:1000500",
                             name="scan window upper limit",
                             value=str(spectrum['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']),
                             unitCvRef="MS",
                             unitAccession="MS:1000040",
                             unitName="m/z")
        
        # MS2 data
        if spectrum['ms level'] == 2:
            precursorList = etree.SubElement(spectrum_element,
                                             "precursorList",
                                             count="1")
            precursor = etree.SubElement(precursorList,
                                         "precursor")
                                         
            isolationWindow = etree.SubElement(precursor,
                                               'isolationWindow')
            if 'isolationWindow' in spectrum['precursorList']['precursor'][0].keys():
                etree.SubElement(isolationWindow,
                                 'cvParam',
                                 cvRef="MS", 
                                 accession="MS:1000827",
                                 name="isolation window target m/z",
                                 value= str(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z']),
                                 unitCvRef="MS",
                                 unitAccession="MS:1000040",
                                 unitName="m/z")
                etree.SubElement(isolationWindow,
                                 'cvParam',
                                 cvRef="MS", 
                                 accession="MS:1000828",
                                 name="isolation window lower offset",
                                 value= str(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window lower offset']),
                                 unitCvRef="MS",
                                 unitAccession="MS:1000040",
                                 unitName="m/z")
                etree.SubElement(isolationWindow,
                                 'cvParam',
                                 cvRef="MS", 
                                 accession="MS:1000829",
                                 name="isolation window upper offset",
                                 value= str(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window upper offset']),
                                 unitCvRef="MS",
                                 unitAccession="MS:1000040",
                                 unitName="m/z")
                             
            selectedIonList = etree.SubElement(precursor,
                                               "selectedIonList",
                                               count="1")
            selectedIon = etree.SubElement(selectedIonList,
                                           "selectedIon")
            
            if mode == 'calibration':
                mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                if equation[0] == 'linear':
                    corrected_mz = mz + ((equation[1]*mz)+equation[2])
                elif equation[0] == 'quadratic':
                    corrected_mz = mz+((equation[1]*(mz**2))+(equation[2]*mz)+equation[3])
                spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'] = corrected_mz
                
            etree.SubElement(selectedIon,
                             "cvParam",
                             cvRef = "MS",
                             accession = "MS:1000744",
                             name = "selected ion m/z",
                             value = str(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']),
                             unitCvRef="MS",
                             unitAccession="MS:1000040",
                             unitName="m/z")
            if 'charge state' in spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0].keys():
                etree.SubElement(selectedIon,
                                 "cvParam",
                                 cvRef = "MS",
                                 accession = "MS:1000041",
                                 name = "charge state",
                                 value = str(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']))
                             
            activation = etree.SubElement(precursor,
                                          'activation')
                         
        # Add m/z array
        mz_array_binary_string = ''
        if len(spectrum['m/z array']) > 0:
            if mode == 'calibration':
                corrected_mz_array = correct_mz_array(spectrum['m/z array'])
                spectrum['m/z array'] = corrected_mz_array
            mz_array_binary_string = convert_to_base64(spectrum['m/z array'])
            
        binaryDataArrayList = etree.SubElement(spectrum_element,
                                               'binaryDataArrayList',
                                               count="2")
        mz_array_binary = etree.SubElement(binaryDataArrayList,
                                           'binaryDataArray',
                                           encodedLength=f"{len(mz_array_binary_string)}")
        etree.SubElement(mz_array_binary,
                         "cvParam",
                         cvRef='MS',
                         accession='MS:1000523',
                         name='64-bit float')
        etree.SubElement(mz_array_binary,
                         "cvParam",
                         cvRef="MS",
                         accession="MS:1000576",
                         name="no compression",
                         value="")
        etree.SubElement(mz_array_binary,
                         "cvParam",
                         cvRef='MS',
                         accession='MS:1000514',
                         name='m/z array',
                         unitAccession='MS:1000040',
                         unitCvRef='MS',
                         unitName='m/z',
                         value='')
        mz_array_binary_element = etree.SubElement(mz_array_binary,
                                                   'binary')
        if len(spectrum['m/z array']) > 0:
            mz_array_binary_element.text = mz_array_binary_string

        # Add intensity array
        int_array_binary_string = ''
        if len(spectrum['intensity array']) > 0:
            int_array_binary_string = convert_to_base64(spectrum['intensity array'])
            
        int_array_binary = etree.SubElement(binaryDataArrayList,
                                           'binaryDataArray',
                                           encodedLength=f"{len(int_array_binary_string)}")
        etree.SubElement(int_array_binary,
                         "cvParam",
                         cvRef='MS',
                         accession='MS:1000523',
                         name='64-bit float')
        etree.SubElement(int_array_binary,
                         "cvParam",
                         cvRef="MS",
                         accession="MS:1000576",
                         name="no compression",
                         value="")
        etree.SubElement(int_array_binary,
                         "cvParam",
                         cvRef='MS',
                         accession='MS:1000515',
                         name='intensity array',
                         unitAccession='MS:1000131',
                         unitCvRef='MS',
                         unitName='number of detector counts',
                         value='')
        int_array_binary_element = etree.SubElement(int_array_binary,
                                                   'binary')
        if len(spectrum['intensity array']) > 0:
            int_array_binary_element.text = int_array_binary_string
    
    # Update the number of spectra in the file
    spectrumList.set('count', str(written_spectra))

    # Convert to a pretty XML string, decoded for indexing
    xml_str = etree.tostring(mzML, pretty_print=True, xml_declaration=True, encoding='UTF-8').decode('utf-8')

    # Find all positions of <spectrum>
    tag = "<spectrum "
    positions = {}
    counter = 1
    index_list_offset = len(xml_str)-13

    for index, char in enumerate(xml_str):
        try:
            if char == '<' and xml_str[index:index+len(tag)] == tag:
                positions[f'scan={counter}'] = index
                counter+= 1
        except:
            break

    indexList = etree.SubElement(mzML,
                                 'indexList',
                                 count='1')
    index = etree.SubElement(indexList,
                             'index',
                             name='spectrum')
    for scan in positions:
        current_index = etree.SubElement(index,
                                         'offset',
                                         idRef=scan)
        current_index.text = str(positions[scan])
        
    indexListOffset = etree.SubElement(mzML,
                                       'indexListOffset')
    indexListOffset.text = str(index_list_offset)

    # Convert to a pretty XML string, encoded, with the indexed part
    xml_str = etree.tostring(mzML, pretty_print=True, xml_declaration=True, encoding='UTF-8')

    # Save to a new mzML file
    with open(os.path.join(output_file_path, new_file_name), "wb") as f:
        f.write(xml_str)
    
    if mode == 'calibration':
        make_calibration_report(output_file_path, new_file_name, calibration_parameters[0], calibration_parameters[1], calibration_parameters[2], spectra_file_path)
        
    return os.path.join(output_file_path, new_file_name)
    
# Custom square root function that handles negative values
def custom_sqrt(x):
    return numpy.where(x >= 0, numpy.sqrt(x), -numpy.sqrt(-x))

# Custom square function (inverse of sqrt)
def custom_square(x):
    return numpy.where(x >= 0, numpy.square(x), -numpy.square(-x))
    
def add_calibrant(treeview, file_name, mz, standard, glycan, library):
    '''
    '''
    global glycans, data_storage, last_selection
    
    if type(file_name) == list:
        spectra_file = find_sample_path(file_name[0].get(), file_name[1])
    else:
        spectra_file = file_path
    
    # Get file name
    file_name = sample_name_from_path(spectra_file)
    
    targets = [mz.get(), standard.get(), glycan.get()]
    count = 0
    good_target = ''
    good_target_type = ''
    for index, target in enumerate(targets):
        if len(target) > 0:
            count += 1
            good_target = target
            if index == 0:
                good_target_type = 'mz'
            if index == 1:
                good_target_type = 'standard'
            if index == 2:
                good_target_type = 'glycan'
    if count > 1:
        error_window("You must selected only one of the three options: m/z, Standard or From Library.")
    elif count == 1:
        if good_target_type == 'mz':
            try:
                float(targets[0])
            except:
                error_window("Invalid input on the m/z field.")
                return
            item_to_add = [targets[0], targets[0], '', '', '', '']
            treeview.insert('', 'end', values=item_to_add)
        elif good_target_type == 'standard':
            for number, value in enumerate(standard_calibrants[good_target]):
                item_to_add = [f"{good_target}({number+1})", value, '', '', '', '']
                treeview.insert('', 'end', values=item_to_add)
        elif good_target_type == 'glycan':
            item_to_add = [good_target, glycans[good_target], '', '', '', '']
            treeview.insert('', 'end', values=item_to_add)
            calibrants_list_sort(treeview, 'Ref. m/z', False)
        
        # Clear the entries
        mz.delete(0, tk.END)
        standard.set('')
        glycan.set('')
        
    last_selection[file_name] = [treeview.item(calibrant, "values") for calibrant in treeview.get_children()]
        
def remove_selected_item(treeview, file_name, plot_area, labels, mode_combobox, search_tolerance, intensity_threshold):
    '''
    '''
    global file_path, last_selection
    selected = treeview.selection()
    for item in selected:
        treeview.delete(item)
    
    if type(file_name) == list:
        spectra_file = find_sample_path(file_name[0].get(), file_name[1])
    else:
        spectra_file = file_path
    
    # Get file name
    file_name = sample_name_from_path(spectra_file)
    file_name_plus_data = f"{file_name}_{search_tolerance}_{intensity_threshold}"
    
    edit_calibrant_list(treeview, file_name_plus_data, plot_area, labels, mode_combobox, search_tolerance)
    
    if len(treeview.get_children()) == 0:
        plot_area[0].cla()
        plot_area[0].set_xlabel('m/z', fontsize=8)
        plot_area[0].set_title('m/z error', fontsize=8)
        plot_area[1].draw()
    
        labels[0].config(text="")
        labels[1].config(text="")
        labels[2].config(text="")
        
    last_selection[file_name] = [treeview.item(calibrant, "values") for calibrant in treeview.get_children()]
        
def remove_all_items(treeview, file_name, plot_area, labels):
    '''
    '''
    global last_selection
    
    if type(file_name) == list:
        spectra_file = find_sample_path(file_name[0].get(), file_name[1])
    else:
        spectra_file = file_path
    
    # Get file name
    file_name = sample_name_from_path(spectra_file)
    
    selected = treeview.get_children()
    for item in selected:
        treeview.delete(item)
        
    plot_area[0].cla()
    plot_area[0].set_xlabel('m/z', fontsize=8)
    plot_area[0].set_title('m/z error', fontsize=8)
    plot_area[1].draw()
    
    labels[0].config(text="")
    labels[1].config(text="")
    labels[2].config(text="")
        
    last_selection[file_name] = [treeview.item(calibrant, "values") for calibrant in treeview.get_children()]
        
def load_file(labels, mzml_window, trimming_args=[]):
    '''
    '''
    global file_path
    file_path = filedialog.askopenfilename(filetypes=[("MzXML or MzML files", "*.mzML *.mzXML"), ("All files", "*.*")])
    file_name = file_path.split("/")[-1]
    if len(file_name) > 40:
        file_name = file_name[:40]+'...'
    for label in labels:
        ToolTip(label, f"{file_path}")
        label.config(text=f"File: {file_name}")
    if len(file_path) > 0 and len(trimming_args) > 0:
        loading_bpc_data(file_path, trimming_args[0], trimming_args[1], trimming_args[2])
    mzml_window.lift()
    
    
def loading_bpc_data(file_path, trimmer_fig, trimmer_ax, trimmer_canvas):
    '''
    '''

    def close_cs():
        loading_file_window.destroy()

    def on_closing():
        '''This function is used to remove the function from the Close Window button (x button).'''
        return
        
    def wait_thread():
        bpc_data = load_bpc_data(file_path)
        if len(bpc_data[0]) > 0:
            plot_sample_bpc(trimmer_fig, trimmer_ax, trimmer_canvas, bpc_data[0], bpc_data[1])
        
        close_cs()
        
    loading_file_window = tk.Toplevel()
    loading_file_window.withdraw()
    loading_file_window.title("Calibrating Spectra")
    icon = ImageTk.PhotoImage(ico_image)
    loading_file_window.iconphoto(False, icon)
    loading_file_window.resizable(False, False)
    loading_file_window.grab_set()
    loading_file_window.protocol("WM_DELETE_WINDOW", on_closing)
    
    calibrating_spectra_window_label = ttk.Label(loading_file_window, text="Loading spectra file...\nPlease wait.", font=("Segoe UI", fontsize))
    calibrating_spectra_window_label.pack(pady=35, padx=70)
    
    loading_file_window.update_idletasks()
    loading_file_window.deiconify()
    window_width = loading_file_window.winfo_width()
    window_height = loading_file_window.winfo_height()
    screen_width = loading_file_window.winfo_screenwidth()
    screen_height = loading_file_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    loading_file_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
    t = threading.Thread(target=wait_thread)
    t.start()
    
def load_bpc_data(sample):
    file_type = sample.split('.')[-1].lower()
    
    rt_array = []
    bpc = []
    
    if file_type == 'mzml':
        access = mzml.MzML(sample)
        rt_key = 'scanList'
        rt_subkey = 'scan start time'
        ms_key = 'ms level'
        intensity_key = 'base peak intensity'
        
    elif file_type == 'mzxml':
        access = mzxml.MzXML(sample)
        rt_key = 'retentionTime'
        rt_subkey = None
        ms_key = 'msLevel'
        intensity_key = 'basePeakIntensity'
        
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

            if ms_level == 1:
                rt_array.append(current_rt)

                # Get base peak intensity
                intensity_value = i.get(intensity_key, numpy.max(i['intensity array']) if len(i['intensity array']) > 0 else 0.0)
                bpc.append(float(intensity_value))
                
        return [rt_array, bpc], time_unit
    except:
        return [], None
    
def quit_mzml_window(mzml_window, from_GG):
    '''
    '''
    global file_path, close_after_alignment
    
    # While an alignment runs, cancel it and close the window once the alignment thread stops
    if background_task == 'alignment':
        close_after_alignment = True
        stop_run_event.set()
        return
    # The automatic parameters calculation can't be cancelled, so the window can only be closed after it finishes
    if background_task == 'auto parameters':
        return
    
    if from_GG != False:
        check_folder = os.path.join(tempfile.gettempdir())
        
        this_process_id = os.getpid()
        this_process = psutil.Process(this_process_id)
        this_process_ppid = this_process.ppid()
        
        general_temp_folder = os.path.join(tempfile.gettempdir())
        os.remove(os.path.join(general_temp_folder, f"mzml_window_{this_process_id}.txt"))
        
    file_path = ''
    mzml_window.destroy()
    
def load_library_data(library):
    '''
    '''
    with open(library, 'rb') as f:
        library_data = dill.load(f)
        f.close()
    full_library = library_data[0]
    glycans_list = {}
    for glycan in full_library:
        for adduct in full_library[glycan]['Adducts_mz']:
            glycans_list[f"{glycan}({adduct})"] = full_library[glycan]['Adducts_mz'][adduct]
    return glycans_list
    
def on_fit_selected(event, treeview, file_name, plot_area, labels, mode_combobox, search_tolerance, intensity_threshold):
    '''
    '''
    global file_path
        
    if type(file_name) == list:
        spectra_file = find_sample_path(file_name[0].get(), file_name[1])
    else:
        spectra_file = file_path
    
    # Get file name
    file_name = sample_name_from_path(spectra_file)
    file_name = f"{file_name}_{search_tolerance}_{intensity_threshold}"
    
    edit_calibrant_list(treeview, file_name, plot_area, labels, mode_combobox, search_tolerance)
        
def calibrants_list_sort(tv, col, reverse):
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
    tv.heading(col, command=lambda: calibrants_list_sort(tv, col, not reverse))
    
def search_all_samples(spectra_file, treeview, search_tolerance, intensity_threshold, plot_area, labels, mode_combobox, samples_list = [], draw = False):
    '''
    '''
    for sample in samples_list:
        trace(sample_name_from_path(sample), treeview, search_tolerance, intensity_threshold, plot_area, labels, mode_combobox, samples_list, draw = False)
        
    spectra_file = find_sample_path(spectra_file.get(), samples_list)
    file_name = sample_name_from_path(spectra_file)
    file_name_plus_data = f"{file_name}_{search_tolerance}_{intensity_threshold}"
    edit_calibrant_list(treeview, file_name_plus_data, plot_area, labels, mode_combobox, search_tolerance)
    
    
def trace(spectra_file, treeview, search_tolerance, intensity_threshold, plot_area, labels, mode_combobox, samples_list = [], draw=True):
    '''
    '''
    global data_storage, file_path, last_selection
    
    if len(treeview.get_children()) == 0:
        return
    
    if len(samples_list) != 0:
        if draw:
            spectra_file = spectra_file.get()
        else:
            spectra_file = spectra_file
        spectra_file = find_sample_path(spectra_file, samples_list)
    else:
        spectra_file = file_path
    
    if len(spectra_file) == 0 or len(spectra_file.split("/")) == 1:
        return
    
    # Get file name
    file_name = sample_name_from_path(spectra_file)
    file_name_plus_data = f"{file_name}_{search_tolerance}_{intensity_threshold}"
    
    if file_name_plus_data not in data_storage.keys():
        # Add file name to data storage dictionary
        data_storage[file_name_plus_data] = {}
    
    # Start rt_array
    rt_array = []
    
    for target_index, target in enumerate(treeview.get_children()):
    
        # Get the values from the treeview
        target_values = treeview.item(target, "values")
        
        # skip if already searched
        if file_name_plus_data in data_storage.keys():
            if target_values[0] in data_storage[file_name_plus_data].keys():
                continue
        
        # Add target to data storage
        data_storage[file_name_plus_data][target_values[0]] = {}
        
        # Load the file
        if spectra_file.split(".")[-1].lower() == 'mzml':
            spectra = mzml.MzML(spectra_file)
            file_type = 'mzml'
        elif spectra_file.split(".")[-1].lower() == 'mzxml':
            spectra = mzxml.MzXML(spectra_file)
            file_type = 'mzxml'
            
        if file_type == 'mzml' and float(spectra[-1]['scanList']['scan'][0]['scan start time']) > 300:
            time_unit = 'second'
        else:
            time_unit = 'minute'
        
        int_array = []
        found_mz_array = []
            
        for spectrum in spectra:
            if file_type == 'mzxml':
                spectrum['scanList'] = {'scan': [{'scan start time':spectrum['retentionTime']}]}
                spectrum['ms level'] = spectrum['msLevel']
                
            if spectrum['ms level'] != 1:
                continue
                
            if target_index == 0:
                if time_unit == 'second':
                    rt_array.append(spectrum['scanList']['scan'][0]['scan start time']/60)
                else:
                    rt_array.append(spectrum['scanList']['scan'][0]['scan start time'])
            
            found_index = binary_search_with_tolerance(spectrum['m/z array'], float(target_values[1]), 0, len(spectrum['m/z array']), search_tolerance, spectrum['intensity array'])
            
            if found_index != -1:
                if spectrum['intensity array'][found_index] < intensity_threshold:
                    int_array.append(0)
                    found_mz_array.append(0)
                else:
                    int_array.append(spectrum['intensity array'][found_index])
                    found_mz_array.append(spectrum['m/z array'][found_index])
            if found_index == -1:
                int_array.append(0)
                found_mz_array.append(0)
        
        # Add trace to the data storage
        data_storage[file_name_plus_data][target_values[0]]['found_mz_array'] = found_mz_array
        data_storage[file_name_plus_data][target_values[0]]['int_array'] = int_array
    
    if len(rt_array) > 0:
        data_storage[file_name_plus_data]['rt_array'] = rt_array
        
    last_selection[file_name] = [treeview.item(calibrant, "values") for calibrant in treeview.get_children()]
    
    if draw:
        edit_calibrant_list(treeview, file_name_plus_data, plot_area, labels, mode_combobox, search_tolerance)
    
def edit_calibrant_list(treeview, file_name, plot_area, labels, mode_combobox, search_tolerance, start=False):
    '''
    '''
    global data_storage, plotted_datapoints, equation
    
    if file_name not in data_storage:
        return
    
    plot_area[0].cla()
    plotted_datapoints = [[], []]
    r_square_datapoints = [[],[]]
    
    for target in treeview.get_children():
        # Get the values from the treeview
        target_values = list(treeview.item(target, "values"))
        
        if target_values[0] not in data_storage[file_name]:
            continue
            
        if len(data_storage[file_name][target_values[0]]['found_mz_array']) == 0:
            continue
        
        max_id = numpy.argmax(data_storage[file_name][target_values[0]]['int_array'])
        
        if data_storage[file_name][target_values[0]]['found_mz_array'][max_id] == 0:
            continue
        
        max_id_found_mz_array = data_storage[file_name][target_values[0]]['found_mz_array'][max_id]
        max_id_rt_array = data_storage[file_name]['rt_array'][max_id]
        target_values[2] = f"{max_id_found_mz_array:.4f}"
        target_values[3] = f"{calculate_ppm_diff(float(target_values[2]), float(target_values[1])):.2f}"
        target_values[4] = f"{float(target_values[1])-float(target_values[2]):.4f}"
        target_values[5] = f"{max_id_rt_array:.2f}"
        
        plotted_datapoints[0].append(float(target_values[1]))
        plotted_datapoints[1].append(float(target_values[4]))
        
        treeview.item(target, values=target_values)
        
    if len(plotted_datapoints[0]) < 3:
        error_window("You need at least three datapoints to calculate the calibration.")
        return
        
    for index, item in enumerate(plotted_datapoints[0]):
        plot_area[0].scatter([plotted_datapoints[0][index]], [plotted_datapoints[1][index]], marker='o', s=50, c='blue')
        
    if mode_combobox.get() == 'Linear':
        # Calculate the standard deviation of the datapoints
        y_std = numpy.std(plotted_datapoints[1])
        
        # Calculate the trendline
        coefficients = numpy.polyfit(plotted_datapoints[0], plotted_datapoints[1], 1)  # 1 means a linear fit (degree 1)
        slope, intercept = coefficients
        equation = ['linear', slope, intercept]
        trendline = numpy.polyval(coefficients, plotted_datapoints[0])
        
        # Add the trendline equation to the interface
        if intercept >= 0:
            trendline_equation = f"y = {slope:.2f}x + {intercept:.2f}"
        else:
            trendline_equation = f"y = {slope:.2f}x {intercept:.2f}"
    elif mode_combobox.get() == 'Quadratic':
        # Calculate the standard deviation of the datapoints
        y_std = numpy.std(plotted_datapoints[1])
        
        # Calculate the trendline
        coefficients = numpy.polyfit(plotted_datapoints[0], plotted_datapoints[1], 2)  # 1 means a linear fit (degree 1)
        a, b, c = coefficients
        equation = ['quadratic', a, b, c]
        trendline = numpy.polyval(coefficients, plotted_datapoints[0])
        
        # Add the trendline equation to the interface
        if b >= 0 and c < 0:
            trendline_equation = f"y = {a:.2f}x² + {b:.2f}x {c:.2f}"
        elif b < 0 and c >= 0:
            trendline_equation = f"y = {a:.2f}x² {b:.2f}x + {c:.2f}"
        elif b < 0 and c < 0:
            trendline_equation = f"y = {a:.2f}x² {b:.2f}x {c:.2f}"
        else:
            trendline_equation = f"y = {a:.2f}x² + {b:.2f}x + {c:.2f}"
        
    # Calculate the R² and add it to the interface
    corrected_y_data = numpy.array(plotted_datapoints[1])
    corrected_trendline = numpy.array(trendline)
    
    ss_total = numpy.sum((corrected_y_data - numpy.mean(corrected_y_data)) ** 2)  # Total sum of squares
    ss_residual = numpy.sum((corrected_y_data - corrected_trendline) ** 2)  # Residual sum of squares
    r_squared = 1 - (ss_residual / ss_total)

    # Plot the trendline
    plot_area[0].plot(plotted_datapoints[0], trendline, color='black')
    
    labels[0].config(text=f"m/z error SD: {y_std:.4f}")
    labels[1].config(text=f"Fit R²: {r_squared:.4f}")
    labels[2].config(text=f"Fit equation: {trendline_equation}")
    
    plot_area[0].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plot_area[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plot_area[0].set_xlabel('m/z', fontsize=8)
    plot_area[0].set_title('m/z error', fontsize=8)
    
    plot_area[0].set_ylim([-search_tolerance, search_tolerance])
    if not start:
        plot_area[1].draw_idle()
    
def correct_mz_array(mz_array):
    '''
    '''
    global equation
    
    new_array = []
    if equation[0] == 'linear':
        for mz in mz_array:
            new_array.append(mz+((equation[1]*mz)+equation[2]))
    elif equation[0] == 'quadratic':
        for mz in mz_array:
            new_array.append(mz+((equation[1]*(mz**2))+(equation[2]*mz)+equation[3]))
    
    return numpy.array(new_array)
    
def finding_calibrants(spectra_file, treeview, search_tolerance, intensity_threshold, plot_area, labels, mode_combobox, samples_list = [], draw=True):

    def close_fc():
        finding_calibrant_window.destroy()

    def on_closing():
        '''This function is used to remove the function from the Close Window button (x button).'''
        return
        
    def wait_thread():
        if draw:
            trace(spectra_file, treeview, search_tolerance, intensity_threshold, plot_area, labels, mode_combobox, samples_list, draw)
        else:
            search_all_samples(spectra_file, treeview, search_tolerance, intensity_threshold, plot_area, labels, mode_combobox, samples_list, draw)
        
        close_fc()
        
    finding_calibrant_window = tk.Toplevel()
    # finding_calibrant_window.attributes("-topmost", True)
    finding_calibrant_window.withdraw()
    finding_calibrant_window.title("Finding Calibrants")
    icon = ImageTk.PhotoImage(ico_image)
    finding_calibrant_window.iconphoto(False, icon)
    finding_calibrant_window.resizable(False, False)
    finding_calibrant_window.grab_set()
    finding_calibrant_window.protocol("WM_DELETE_WINDOW", on_closing)
    
    if draw:
        finding_calibrant_window_label = ttk.Label(finding_calibrant_window, text="Finding calibrants in the spectra file...", font=("Segoe UI", fontsize))
    else:
        finding_calibrant_window_label = ttk.Label(finding_calibrant_window, text="Finding calibrants in all spectra files...\nThis may take a while, please wait.", font=("Segoe UI", fontsize))
    finding_calibrant_window_label.pack(pady=35, padx=70)
    
    finding_calibrant_window.update_idletasks()
    finding_calibrant_window.deiconify()
    window_width = finding_calibrant_window.winfo_width()
    window_height = finding_calibrant_window.winfo_height()
    screen_width = finding_calibrant_window.winfo_screenwidth()
    screen_height = finding_calibrant_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    finding_calibrant_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
    t = threading.Thread(target=wait_thread)
    t.start()
    
def calibrating_spectra(spectra_file_path, mode, calibrants_list, labels, plot_area, file_list = []):
    '''
    '''

    def close_cs():
        calibrating_spectra_window.destroy()

    def on_closing():
        '''This function is used to remove the function from the Close Window button (x button).'''
        return
        
    def wait_thread():
        new_file_name = create_mzml(spectra_file_path, mode, file_list, calibration_parameters = [calibrants_list, labels, plot_area])
        
        close_cs()
        
        if len(file_list) > 0:
            calibration_done(new_file_name, True)
        else:
            calibration_done(new_file_name, False)
        
    global file_path
    
    if len(file_list) != 0:
        test_spectra_file_path = find_sample_path(spectra_file_path.get(), file_list)
    else:
        test_spectra_file_path = file_path
        
    calibrating_spectra_window = tk.Toplevel()
    # calibrating_spectra_window.attributes("-topmost", True)
    calibrating_spectra_window.withdraw()
    calibrating_spectra_window.title("Calibrating Spectra")
    icon = ImageTk.PhotoImage(ico_image)
    calibrating_spectra_window.iconphoto(False, icon)
    calibrating_spectra_window.resizable(False, False)
    calibrating_spectra_window.grab_set()
    calibrating_spectra_window.protocol("WM_DELETE_WINDOW", on_closing)
    
    calibrating_spectra_window_label = ttk.Label(calibrating_spectra_window, text="Calibrating spectra file...\nPlease wait.", font=("Segoe UI", fontsize))
    calibrating_spectra_window_label.pack(pady=35, padx=70)
    
    calibrating_spectra_window.update_idletasks()
    calibrating_spectra_window.deiconify()
    window_width = calibrating_spectra_window.winfo_width()
    window_height = calibrating_spectra_window.winfo_height()
    screen_width = calibrating_spectra_window.winfo_screenwidth()
    screen_height = calibrating_spectra_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    calibrating_spectra_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
    t = threading.Thread(target=wait_thread)
    t.start()
    
def trimming_spectra(spectra_file_path, mode, trimming_range, file_list = []):
    '''
    '''

    def close_cs():
        calibrating_spectra_window.destroy()

    def on_closing():
        '''This function is used to remove the function from the Close Window button (x button).'''
        return
        
    def wait_thread():
        new_file_name = create_mzml(spectra_file_path, mode, file_list, trimming_range)
        
        close_cs()
        
        if len(file_list) > 0:
            trimming_done(new_file_name, True)
        else:
            trimming_done(new_file_name, False)
        
    global file_path
    
    if len(file_list) != 0:
        test_spectra_file_path = find_sample_path(spectra_file_path.get(), file_list)
    else:
        test_spectra_file_path = file_path
        
    calibrating_spectra_window = tk.Toplevel()
    # calibrating_spectra_window.attributes("-topmost", True)
    calibrating_spectra_window.withdraw()
    calibrating_spectra_window.title("Calibrating Spectra")
    icon = ImageTk.PhotoImage(ico_image)
    calibrating_spectra_window.iconphoto(False, icon)
    calibrating_spectra_window.resizable(False, False)
    calibrating_spectra_window.grab_set()
    calibrating_spectra_window.protocol("WM_DELETE_WINDOW", on_closing)
    
    calibrating_spectra_window_label = ttk.Label(calibrating_spectra_window, text="Trimming the spectra file...\nPlease wait.", font=("Segoe UI", fontsize))
    calibrating_spectra_window_label.pack(pady=35, padx=70)
    
    calibrating_spectra_window.update_idletasks()
    calibrating_spectra_window.deiconify()
    window_width = calibrating_spectra_window.winfo_width()
    window_height = calibrating_spectra_window.winfo_height()
    screen_width = calibrating_spectra_window.winfo_screenwidth()
    screen_height = calibrating_spectra_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    calibrating_spectra_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
    t = threading.Thread(target=wait_thread)
    t.start()
    
def trimming_done(file_name, from_GG):
    '''
    '''
    def close_cd():
        calibration_done_window.destroy()

    def on_closing():
        '''This function is used to remove the function from the Close Window button (x button).'''
        return
        
    calibration_done_window = tk.Toplevel()
    # calibration_done_window.attributes("-topmost", True)
    calibration_done_window.withdraw()
    calibration_done_window.title("Trimming Done!")
    icon = ImageTk.PhotoImage(ico_image)
    calibration_done_window.iconphoto(False, icon)
    calibration_done_window.resizable(False, False)
    calibration_done_window.grab_set()
    calibration_done_window.protocol("WM_DELETE_WINDOW", on_closing)
    
    if not from_GG:
        file_name_split = file_name.split("/")[-1].split("\\")[-1]
        calibrating_spectra_window_label = ttk.Label(calibration_done_window, text=f"\nTrimmed file name:\n\n{file_name_split}", font=("Segoe UI", fontsize))
        calibrating_spectra_window_label.pack(pady=10, padx=70)
    else:
        file_name_split = file_name.split("/")[-1].split("\\")[-1]
        calibrating_spectra_window_label = ttk.Label(calibration_done_window, text=f"\nTrimmed file name:\n\n{file_name_split}\n\nYou must load it in GlycoGenius to use it.", font=("Segoe UI", fontsize))
        calibrating_spectra_window_label.pack(pady=10, padx=70)
    
    # ok button
    ok_button = ttk.Button(calibration_done_window, text="Ok", style="small_button_sfw_style1.TButton", command=close_cd)
    ok_button.pack(pady=10, padx=70)
    
    calibration_done_window.update_idletasks()
    calibration_done_window.deiconify()
    window_width = calibration_done_window.winfo_width()
    window_height = calibration_done_window.winfo_height()
    screen_width = calibration_done_window.winfo_screenwidth()
    screen_height = calibration_done_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    calibration_done_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
def calibration_done(file_name, from_GG):
    '''
    '''
    def close_cd():
        calibration_done_window.destroy()

    def on_closing():
        '''This function is used to remove the function from the Close Window button (x button).'''
        return
        
    calibration_done_window = tk.Toplevel()
    # calibration_done_window.attributes("-topmost", True)
    calibration_done_window.withdraw()
    calibration_done_window.title("Calibration Done!")
    icon = ImageTk.PhotoImage(ico_image)
    calibration_done_window.iconphoto(False, icon)
    calibration_done_window.resizable(False, False)
    calibration_done_window.grab_set()
    calibration_done_window.protocol("WM_DELETE_WINDOW", on_closing)
    
    if not from_GG:
        file_name_split = file_name.split("/")[-1].split("\\")[-1]
        calibrating_spectra_window_label = ttk.Label(calibration_done_window, text=f"\nCalibrated file name:\n\n{file_name_split}", font=("Segoe UI", fontsize))
        calibrating_spectra_window_label.pack(pady=10, padx=70)
    else:
        file_name_split = file_name.split("/")[-1].split("\\")[-1]
        calibrating_spectra_window_label = ttk.Label(calibration_done_window, text=f"\nCalibrated file name:\n\n{file_name_split}\n\nYou must load it in GlycoGenius to use it.", font=("Segoe UI", fontsize))
        calibrating_spectra_window_label.pack(pady=10, padx=70)
    
    # ok button
    ok_button = ttk.Button(calibration_done_window, text="Ok", style="small_button_sfw_style1.TButton", command=close_cd)
    ok_button.pack(pady=10, padx=70)
    
    calibration_done_window.update_idletasks()
    calibration_done_window.deiconify()
    window_width = calibration_done_window.winfo_width()
    window_height = calibration_done_window.winfo_height()
    screen_width = calibration_done_window.winfo_screenwidth()
    screen_height = calibration_done_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    calibration_done_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
    
def click_treeview(event, treeview, file_data, spectra_viewer, search_tolerance, intensity_threshold):
    '''
    '''
    global data_storage
    selected_items = treeview.selection()
    
    if len(selected_items) > 1:
        return
    
    spectra_file = find_sample_path(file_data[0].get(), file_data[1])
                
    file_name = sample_name_from_path(spectra_file)
    file_name_plus_data = f"{file_name}_{search_tolerance}_{intensity_threshold}"
    
    for row in selected_items:
        row_data = list(treeview.item(row, "values"))
        
        if row_data[2] == '':
            return
        
        spectra_viewer[0].cla()
        
        # Load the file
        if spectra_file.split(".")[-1].lower() == 'mzml':
            spectra = mzml.MzML(spectra_file)
            file_type = 'mzml'
        elif spectra_file.split(".")[-1].lower() == 'mzxml':
            spectra = mzxml.MzXML(spectra_file)
            file_type = 'mzxml'
            
        if file_type == 'mzml' and float(spectra[-1]['scanList']['scan'][0]['scan start time']) > 300:
            time_unit = 'second'
        else:
            time_unit = 'minute'
            
        rt = float(row_data[5])
        if time_unit == 'second':
            rt = rt*60
        
        x_values = spectra.time[rt]['m/z array']
        y_values = spectra.time[rt]['intensity array']
        
        new_x_data = []
        new_y_data = []
        for index, x in enumerate(x_values):
            new_x_data.append(x-0.000001)
            new_x_data.append(x)
            new_x_data.append(x+0.000001)
            new_y_data.append(0)
            new_y_data.append(y_values[index])
            new_y_data.append(0)
            
        spectra_viewer[0].plot(new_x_data, new_y_data, marker='None', linewidth=1, color='black')
        
        scaling = spectra_viewer[2].get()
        if scaling == 'Linear':
            spectra_viewer[0].set_yscale('linear')
        if scaling == 'Log':
            spectra_viewer[0].set_yscale('symlog')
        if scaling == 'Sqrt':
            spectra_viewer[0].set_yscale('function', functions=(custom_sqrt, custom_square))
        
        spectra_viewer[0].set_xlabel('m/z')
        spectra_viewer[0].set_ylabel('Intensity (AU)')
        
        spectra_viewer[0].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        spectra_viewer[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        
        spectra_viewer[0].set_xlim([float(row_data[1])-5, float(row_data[1])+5])
        spectra_viewer[0].set_ylim([0, max(data_storage[file_name_plus_data][row_data[0]]['int_array'])*1.2])
        
        spectra_viewer[0].add_patch(Rectangle((float(row_data[1])-search_tolerance, spectra_viewer[0].get_ylim()[0]), (float(row_data[1])+search_tolerance) - (float(row_data[1])-search_tolerance), 1000000000000, color='#FEB7A1', alpha=0.3))
        spectra_viewer[0].add_patch(Rectangle((float(row_data[2])-0.001, spectra_viewer[0].get_ylim()[0]), (float(row_data[2])+0.001) - (float(row_data[2])-0.001), 1000000000000, color='#345eeb', alpha=0.3))
        
        rt_formatted = "%.2f" % round(rt if time_unit == 'minute' else rt/60, 2)
        spectra_viewer[3].config(text = f"Retention/Migration Time: {rt_formatted}")
        
        # Calculate the absolute differences
        differences = numpy.abs(x_values - float(row_data[2]))
        
        # Find the index of the minimum difference
        closest_index = numpy.argmin(differences)
        
        spectra_viewer[0].annotate(f'{float(row_data[2]):.4f}', xy=(float(row_data[2]), y_values[closest_index]), xytext=(0, 10), textcoords='offset points', ha='center', fontsize=8, clip_on=True)
        
        spectra_viewer[1].draw()
        
def switch_sample_button(direction, combobox, samples_list, mode, params = []):
    '''
    '''
    sample_name = combobox.get()
    sample_path = find_sample_path(sample_name, samples_list)
            
    index = samples_list.index(sample_path)
    
    if direction == 'back':
        index -= 1
        if index < 0:
            index = 0
    
    elif direction == 'forward':
        index += 1
        if index > len(samples_list)-1:
            index = len(samples_list)-1
    file_name = sample_name_from_path(samples_list[index])
            
    # Set combobox here
    combobox.set(file_name)
            
    mzml_window.title(f"Spectra File Editor - {file_name}")
    
    if mode == 'calibration':
        # load parameters from args
        treeview, plot_area, labels, mode_combobox, search_tolerance, intensity_threshold = params
        
        # load data for given file
        file_name_plus_data = f"{file_name}_{search_tolerance}_{intensity_threshold}"
        edit_calibrant_list(treeview, file_name_plus_data, plot_area, labels, mode_combobox, search_tolerance)
    elif mode == 'trimming':
        # load parameters from args
        trimmer_fig, trimmer_ax, trimmer_canvas, samples_data = params
        
        # plot the BPC
        plot_sample_bpc(trimmer_fig, trimmer_ax, trimmer_canvas, (samples_data[combobox.get()]['rt_array'], samples_data[combobox.get()]['bpc']), samples_data[combobox.get()]['time_unit'])

def plot_sample_bpc(trimmer_fig, trimmer_ax, trimmer_canvas, data, time_unit):
    global selected_rt_range, trimming_range, slider_min, slider_max
    
    trimmer_ax.clear()
    
    if time_unit == 'seconds':
        plot_data_x = [x/60 for x in data[0]]
    else:
        plot_data_x = data[0]
    trimmer_ax.plot(plot_data_x, data[1], linewidth=1, color='red')
    
    trimmer_ax.set_xlim(0, plot_data_x[-1])
        
    trimmer_ax.set_xlabel('Retention/Migration Time (min)')
    trimmer_ax.set_ylabel('Intensity (AU)')
        
    trimmer_ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    trimmer_ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    
    if len(trimming_range) != 0:
        selected_rt_range = trimmer_ax.add_patch(Rectangle((trimming_range[0], -1000000000000), trimming_range[1]-trimming_range[0], 2000000000000, color='#03ecfc', alpha=0.3))
        
        slider_min.config(to=plot_data_x[-1])
        slider_max.config(to=plot_data_x[-1])
        slider_max.set(trimming_range[1])
    else:
        selected_rt_range = trimmer_ax.add_patch(Rectangle((0, -1000000000000), trimmer_ax.get_xlim()[1], 2000000000000, color='#03ecfc', alpha=0.3))

        trimming_range = [float(selected_rt_range.get_xy()[0]), float(selected_rt_range.get_xy()[0])+float(selected_rt_range.get_width())]
        
        slider_min.config(to=plot_data_x[-1])
        slider_max.config(to=plot_data_x[-1])
        slider_max.set(trimming_range[1])
    
    trimmer_canvas.draw()
    
def update_min_rt(value, slider_value, other_slider_value, rectangle, trimmer_canvas, value_entry, from_spinbox = False, typed = False):
    global trimming_range
    
    if value == '':
        return
        
    if from_spinbox:
        value = value.get()
        slider_value.set(value)
        
    if float(value) > float(other_slider_value.get()):
        slider_value.set(float(other_slider_value.get()))
        value = other_slider_value.get()
        if from_spinbox and not typed:
            value_entry.delete(0, tk.END)
            value_entry.insert(0, f"{float(other_slider_value.get())-0.1:.1f}")
    
    if not from_spinbox:
        value_entry.delete(0, tk.END)
        value_entry.insert(0, f"{float(value):.1f}")
    
    x, y = rectangle.get_xy()
    width = rectangle.get_width()
    height = rectangle.get_height()
    
    new_width = float(width)-(float(value)-float(x))
    
    rectangle.set_xy((float(value), -1000000000000))
    rectangle.set_width(new_width)
    trimmer_canvas.draw_idle()
    
    trimming_range = [float(rectangle.get_xy()[0]), float(rectangle.get_xy()[0])+float(rectangle.get_width())]
    
def update_max_rt(value, slider_value, other_slider_value, rectangle, trimmer_canvas, value_entry, from_spinbox = False, typed = False):
    global trimming_range
    
    if value == '':
        return
        
    if from_spinbox:
        value = value.get()
        slider_value.set(value)
        
    if float(value) < float(other_slider_value.get()):
        slider_value.set(float(other_slider_value.get()))
        value = other_slider_value.get()
        if from_spinbox and not typed:
            value_entry.delete(0, tk.END)
            value_entry.insert(0, f"{float(other_slider_value.get())+0.1:.1f}")
    
    if not from_spinbox:
        value_entry.delete(0, tk.END)
        value_entry.insert(0, f"{float(value):.1f}")
    
    x, y = rectangle.get_xy()
    width = rectangle.get_width()
    height = rectangle.get_height()
    
    new_width = float(value)-float(x)
    
    rectangle.set_xy((float(x), -1000000000000))
    rectangle.set_width(new_width)
    trimmer_canvas.draw_idle()
    
    trimming_range = [float(rectangle.get_xy()[0]), float(rectangle.get_xy()[0])+float(rectangle.get_width())]

def set_reference(files_list):
    selected_items = files_list.selection()
    
    if len(selected_items) > 1 or len(selected_items) == 0:
        return
    
    if "reference_sample" not in files_list.item(selected_items[0], "tags"):
        for other_item in files_list.get_children():
            if "reference_sample" in files_list.item(other_item, "tags"):
                files_list.item(other_item, tags=())
        files_list.item(selected_items[0], tags=("reference_sample",))
    else:
        files_list.item(selected_items[0], tags=())
        
    files_list.selection_remove(selected_items)

def update_alignment_settings(settings_scrolltext):
    settings_scrolltext.config(state="normal")
    settings_scrolltext.delete(1.0, tk.END)
    for setting, value in alignment_settings.items():
        settings_scrolltext.insert(tk.END, f"{setting}: {value}\n")
    settings_scrolltext.config(state="disabled")

def open_alignment_settings_window(settings_scrolltext):
    global alignment_settings

    alignment_settings_window = tk.Toplevel()
    alignment_settings_window.withdraw()
    alignment_settings_window.title("Alignment Settings")
    alignment_settings_window.iconphoto(False, icon)

    setting_entries = {}

    def update_numeric_setting(setting):
        # Only store the value if it can be converted to the expected type, otherwise keep the last valid value
        try:
            alignment_settings[setting] = alignment_settings_types[setting](setting_entries[setting].get())
        except ValueError:
            return
        update_alignment_settings(settings_scrolltext)

    for setting, value in alignment_settings.items():
        # Crete a label for the setting
        ttk.Label(alignment_settings_window, text=setting, font=("Segoe UI", fontsize)).grid(row=len(setting_entries), column=0, padx=10, pady=5, sticky="w")

        if setting == 'instrument_type':
            setting_entries[setting] = ttk.Combobox(alignment_settings_window, values=['ORBITRAP', 'TOF', 'QUAD', 'FTICR'], state="readonly")
            setting_entries[setting].set(value)
            setting_entries[setting].grid(row=len(setting_entries)-1, column=1, padx=10, pady=5)
            setting_entries[setting].bind("<<ComboboxSelected>>", lambda event, s=setting: (alignment_settings.update({s: setting_entries[s].get()}), update_alignment_settings(settings_scrolltext)))

        elif setting == 'polarity':
            setting_entries[setting] = ttk.Combobox(alignment_settings_window, values=['+', '-'], state="readonly")
            setting_entries[setting].set(value)
            setting_entries[setting].grid(row=len(setting_entries)-1, column=1, padx=10, pady=5)
            setting_entries[setting].bind("<<ComboboxSelected>>", lambda event, s=setting: (alignment_settings.update({s: setting_entries[s].get()}), update_alignment_settings(settings_scrolltext)))

        else:
            # Create the entry widget
            setting_entries[setting] = ttk.Entry(alignment_settings_window)
            setting_entries[setting].insert(0, str(value))
            setting_entries[setting].grid(row=len(setting_entries)-1, column=1, padx=10, pady=5)
            setting_entries[setting].bind("<KeyRelease>", lambda event, s=setting: update_numeric_setting(s))

        ToolTip(setting_entries[setting], alignment_settings_descriptions[setting])

    ttk.Button(alignment_settings_window, text="Close", style="small_button_sfw_style1.TButton", command=alignment_settings_window.destroy).grid(row=len(setting_entries), column=0, columnspan=2, pady=10)

    # Parent window position and size
    parent_x = mzml_window.winfo_x()
    parent_y = mzml_window.winfo_y()
    parent_width = mzml_window.winfo_width()

    # Position to the right of parent
    new_x = parent_x + parent_width
    new_y = parent_y

    # Place to the right of parent
    alignment_settings_window.geometry(f"+{new_x}+{new_y}")
    alignment_settings_window.resizable(False, False)
    alignment_settings_window.deiconify()

def infer_res(x, y, max_rel_overlap=0.01, min_sigma=1e-12):
    '''This functions plots non-overlapping gaussians over the data to infer the resolution.'''

    x = numpy.asarray(x, dtype=float)
    y = numpy.asarray(y, dtype=float)
    if x.size < 10:
        return numpy.array([])

    if not (0 < max_rel_overlap < 1):
        raise ValueError("max_rel_overlap must be in (0,1)")

    n = x.size

    # distances between consecutive points in the inumpyut order (absolute)
    d_left = numpy.empty(n)
    d_right = numpy.empty(n)
    d_left[0] = numpy.inf
    d_left[1:] = numpy.abs(x[1:] - x[:-1])
    d_right[-1] = numpy.inf
    d_right[:-1] = numpy.abs(x[1:] - x[:-1])

    # sigma limit so that contribution at neighbor center <= max_rel_overlap
    factor = numpy.sqrt(-2.0 * numpy.log(max_rel_overlap))
    sigma_lim_left = numpy.where(numpy.isfinite(d_left), d_left / factor, numpy.inf)
    sigma_lim_right = numpy.where(numpy.isfinite(d_right), d_right / factor, numpy.inf)

    sigmas = numpy.minimum(sigma_lim_left, sigma_lim_right)

    if numpy.isinf(sigmas).all():
        sigmas = numpy.full(n, max(min_sigma, 1e-6))
    else:
        finite_min = numpy.min(sigmas[numpy.isfinite(sigmas)])
        sigmas[numpy.isinf(sigmas)] = finite_min

    sigmas = numpy.maximum(sigmas, min_sigma)

    # vectorized FWHM with safe fallback and cap
    fwhms = 2.3548 * sigmas
    fallback_fwhm = 1.0
    max_fwhm = 1e6
    fwhms = numpy.where(numpy.isfinite(fwhms) & (fwhms > 0), fwhms, fallback_fwhm)
    fwhms = numpy.minimum(fwhms, max_fwhm)

    # vectorized resolution (safe divide)
    res = numpy.divide(x, fwhms, out=numpy.zeros_like(x, dtype=float), where=(fwhms != 0))

    return res

def auto_parameters(samples_list_widget, settings_scrolltext, buttons, from_GG=False):
    global alignment_settings, background_task
    
    background_task = 'auto parameters'
    for name, button in buttons.items():
        button.config(state="disabled")
    block_tabs(True)

    try:
        def get_ms_level(spectrum):
            return spectrum.get('ms level', spectrum.get('msLevel'))

        if len(samples_list_widget.get_children()) == 0:
            error_window("No samples found. Please add at least one sample to the list.")
            return
    
        reference_sample = None
        samples_list = []
        for item_id in samples_list_widget.get_children():
            if "reference_sample" in samples_list_widget.item(item_id, "tags"):
                reference_sample = samples_list_widget.item(item_id, "text")
            else:
                samples_list.append(samples_list_widget.item(item_id, "text"))

        if reference_sample == None:
            reference_sample = samples_list[0]

        # Check whether it's a mzML or mzXML file and load accordingly
        if reference_sample.split(".")[-1].lower() == 'mzml':
            file_type = 'mzml'
        elif reference_sample.split(".")[-1].lower() == 'mzxml':
            file_type = 'mzxml'
    
        # Infer the resolution from the reference sample
        if file_type == 'mzml':
            spectra = mzml.MzML(reference_sample)
        elif file_type == 'mzxml':
            spectra = mzxml.MzXML(reference_sample)

        # Figure out time unit
        if file_type == 'mzml' and float(spectra[-1]['scanList']['scan'][0]['scan start time']) > 300:
            time_unit = 'second'
        else:
            time_unit = 'minute'

        perc_array_of_resolutions = []

        for spectrum in spectra:
            if get_ms_level(spectrum) != 1:
                continue
            centers = spectrum['m/z array']
            amps = spectrum['intensity array']
            resolutions = infer_res(centers, amps, max_rel_overlap=0.01)
            if len(resolutions) > 0:
                perc_array_of_resolutions.append(numpy.percentile(resolutions, 95))
    
        inferred_resolution = numpy.percentile(perc_array_of_resolutions, 95)

        # Infer the FWHM_rt from the reference sample
        if file_type == 'mzml':
            spectra = mzml.MzML(reference_sample)
        elif file_type == 'mzxml':
            spectra = mzxml.MzXML(reference_sample)

        # Find the optimal target
        highest_intensity = 0
        highest_intensity_mz = None
        for spectrum in spectra:
            if get_ms_level(spectrum) != 1:
                continue
            if len(spectrum['intensity array']) > 0 and max(spectrum['intensity array']) > highest_intensity:
                highest_intensity = max(spectrum['intensity array'])
                highest_intensity_mz = spectrum['m/z array'][spectrum['intensity array'].argmax()]

        # Trace EIC
        tolerance = 0.05
        eic = []
        eic_rts = []
        if file_type == 'mzml':
            spectra = mzml.MzML(reference_sample)
        elif file_type == 'mzxml':
            spectra = mzxml.MzXML(reference_sample)
        for spectrum in spectra:
            if get_ms_level(spectrum) != 1:
                continue
            id = binary_search_with_tolerance(spectrum['m/z array'], highest_intensity_mz, 0, len(spectrum['m/z array'])-1, tolerance, spectrum['intensity array'])
            eic.append(spectrum['intensity array'][id] if id != -1 else 0)
            if file_type == 'mzml':
                eic_rts.append(float(spectrum['scanList']['scan'][0]['scan start time'])*(1 if time_unit == 'second' else 60))
            else:
                eic_rts.append(float(spectrum['retentionTime'])*60)

        # Find FWHM
        max_intensity = max(eic)
        max_intensity_id = eic.index(max_intensity)
        half_max = max_intensity / 2

        # Seek on the left side
        for i in range(max_intensity_id, -1, -1):
            if eic[i] <= half_max:
                left_id = i
                break

        # Seek on the right side
        for i in range(max_intensity_id, len(eic)):
            if eic[i] <= half_max:
                right_id = i
                break
    
        fwhm = eic_rts[right_id] - eic_rts[left_id]

        # To finalize, get the polarity from the sample
        if file_type == 'mzml':
            polarity = '+' if 'positive scan' in spectra[-1] else '-'
        elif file_type == 'mzxml':
            polarity = spectra[-1].get('polarity', alignment_settings['polarity'])

        # Update the alignment settings
        alignment_settings['resolution_ms1'] = int(inferred_resolution)
        alignment_settings['fwhm_rt'] = round(fwhm, 1)
        alignment_settings['polarity'] = polarity
        update_alignment_settings(settings_scrolltext)
        
    except Exception as e:
        error_window(f"Could not set the alignment parameters automatically.\nError message: {str(e)}")
        mzml_window.lift()
        mzml_window.focus_set()
    finally:
        # Always re-enable the interface, even if the parameters calculation failed
        background_task = None
        block_tabs(False)
        for name, button in buttons.items():
            if name == 'align' or name == 'add' or name == 'remove' or name == 'ref' or name == 'auto' or name == 'param' or name == 'close':
                button.config(state="normal")
            if name == 'cancel':
                button.config(state="disabled")

def align_sample_pair(ref, sample, progress_bar_info=None):
    try:
        # Update progress bar text
        if stop_run_event.is_set():
            stop_run_event.clear()
            raise Exception("Process stopped by user.")
        progress_bar_info[1].config(text=f"File: {progress_bar_info[3]}/{progress_bar_info[4]}. Loading data...")

        # Check whether it's a mzML or mzXML file and load accordingly
        if ref.split(".")[-1].lower() == 'mzml':
            file_type = 'mzml'
        elif ref.split(".")[-1].lower() == 'mzxml':
            file_type = 'mzxml'

        # Check if ref and sample file are same extension
        if sample.split(".")[-1].lower() != ref.split(".")[-1].lower():
            raise Exception("The files you are trying to align do not have the same extension. Please make sure both files are either mzML or mzXML.")

        # Load data into PASTAQ. ~5% of the progress bar...
        if file_type == 'mzml':
            ref_data = pq.read_mzml(ref, instrument_type = alignment_settings['instrument_type'], resolution_ms1 = float(alignment_settings['resolution_ms1']), resolution_msn=17500, reference_mz = float(alignment_settings['reference_mz']), fwhm_rt = float(alignment_settings['fwhm_rt']), polarity = alignment_settings['polarity'], ms_level = 1)
            align_data = pq.read_mzml(sample, instrument_type = alignment_settings['instrument_type'], resolution_ms1 = float(alignment_settings['resolution_ms1']), resolution_msn=17500, reference_mz = float(alignment_settings['reference_mz']), fwhm_rt = float(alignment_settings['fwhm_rt']), polarity = alignment_settings['polarity'], ms_level = 1)
        elif file_type == 'mzxml':
            ref_data = pq.read_mzxml(ref, instrument_type = alignment_settings['instrument_type'], resolution_ms1 = float(alignment_settings['resolution_ms1']), resolution_msn=17500, reference_mz = float(alignment_settings['reference_mz']), fwhm_rt = float(alignment_settings['fwhm_rt']), polarity = alignment_settings['polarity'], ms_level = 1)
            align_data = pq.read_mzxml(sample, instrument_type = alignment_settings['instrument_type'], resolution_ms1 = float(alignment_settings['resolution_ms1']), resolution_msn=17500, reference_mz = float(alignment_settings['reference_mz']), fwhm_rt = float(alignment_settings['fwhm_rt']), polarity = alignment_settings['polarity'], ms_level = 1)

        # Update progress bar
        if stop_run_event.is_set():
            stop_run_event.clear()
            raise Exception("Process stopped by user.")
        progress_bar_info[0]['value'] += (100/progress_bar_info[4])*0.05
        progress_bar_info[1].config(text=f"File: {progress_bar_info[3]}/{progress_bar_info[4]}. Resampling data...")
        progress_bar_info[2].config(text=f"{int(progress_bar_info[0]['value'])}% Done")

        # Resample the data. ~60% of the progress bar...
        grid1 = pq.resample(ref_data, alignment_settings['num_samples_mz'], alignment_settings['num_samples_rt'], alignment_settings['smoothing_coefficient_mz'], alignment_settings['smoothing_coefficient_rt'])
        grid2 = pq.resample(align_data, alignment_settings['num_samples_mz'], alignment_settings['num_samples_rt'], alignment_settings['smoothing_coefficient_mz'], alignment_settings['smoothing_coefficient_rt'])

        # Update progress bar
        if stop_run_event.is_set():
            stop_run_event.clear()
            raise Exception("Process stopped by user.")
        progress_bar_info[0]['value'] += (100/progress_bar_info[4])*0.6
        progress_bar_info[1].config(text=f"File: {progress_bar_info[3]}/{progress_bar_info[4]}. Identifying peaks...")
        progress_bar_info[2].config(text=f"{int(progress_bar_info[0]['value'])}% Done")

        # Identify peaks. ~5% of the progress bar...
        peaks1 = pq.find_peaks(ref_data, grid1, 1000000)
        peaks2 = pq.find_peaks(align_data, grid2, 1000000)

        # Update progress bar
        if stop_run_event.is_set():
            stop_run_event.clear()
            raise Exception("Process stopped by user.")
        progress_bar_info[0]['value'] += (100/progress_bar_info[4])*0.05
        progress_bar_info[1].config(text=f"File: {progress_bar_info[3]}/{progress_bar_info[4]}. Calculating time map...")
        progress_bar_info[2].config(text=f"{int(progress_bar_info[0]['value'])}% Done")

        # Generate time alignment map. ~5% of the progress bar...
        time_map_1r_2s = pq.calculate_time_map(peaks1, peaks2, alignment_settings['slack'], alignment_settings['segment_length'], alignment_settings['num_points'], alignment_settings['expand_factor'], alignment_settings['peaks_per_window'])

        # Update progress bar
        if stop_run_event.is_set():
            stop_run_event.clear()
            raise Exception("Process stopped by user.")
        progress_bar_info[0]['value'] += (100/progress_bar_info[4])*0.05
        progress_bar_info[1].config(text=f"File: {progress_bar_info[3]}/{progress_bar_info[4]}. Creating aligned MzML file...")
        progress_bar_info[2].config(text=f"{int(progress_bar_info[0]['value'])}% Done")

        # Convert time alignment map to useable format. ~5% of the progress bar...
        index = 0
        alignment_parameters_sample = [sample, []]
        for time_ref, time_target in zip(time_map_1r_2s.rt_start, time_map_1r_2s.sample_rt_start):
            alignment_parameters_sample[1].append([time_ref, time_map_1r_2s.rt_end[index], time_target, time_map_1r_2s.sample_rt_end[index]])
            index += 1

        # Create mzML file with aligned retention times
        create_mzml(None, 'aligning', [], calibration_parameters = [], trimming_parameters = [], alignment_parameters = alignment_parameters_sample)

        # Update progress bar
        if stop_run_event.is_set():
            stop_run_event.clear()
            raise Exception("Process stopped by user.")
        progress_bar_info[0]['value'] += (100/progress_bar_info[4])*0.20
        progress_bar_info[2].config(text=f"{int(progress_bar_info[0]['value'])}% Done")

        return "success"
    except Exception as e:
        if str(e) == "Process stopped by user.":
            if not close_after_alignment:
                messagebox.showinfo("Alignment Cancelled", "The alignment was cancelled.")
        else:
            error_window(f"An error occurred while aligning {sample.split("/")[-1]} to {ref.split("/")[-1]}.\nError message: {str(e)}")
        mzml_window.lift()
        mzml_window.focus_set()
        return "error"

def align_samples(samples_list_widget, buttons, progress_bar_info, from_GG=False):
    global background_task, close_after_alignment
    
    reference_sample = None
    samples_list = []
    for item_id in samples_list_widget.get_children():
        if "reference_sample" in samples_list_widget.item(item_id, "tags"):
            reference_sample = samples_list_widget.item(item_id, "text")
        else:
            samples_list.append(samples_list_widget.item(item_id, "text"))

    if reference_sample == None:
        error_window("You must select a reference sample to align to.")
        mzml_window.lift()
        mzml_window.focus_set()
        return
    
    for name, button in buttons.items():
        if name == 'align' or name == 'add' or name == 'remove' or name == 'ref' or name == 'auto' or name == 'param' or name == 'close':
            button.config(state="disabled")
        if name == 'cancel':
            button.config(state="normal")

    progress_bar_info[0]['value'] = 0
    progress_bar_info[1].config(text=f"")
    progress_bar_info[2].config(text=f"{0}% Done")

    # Clear any cancel request left over from a previous run
    stop_run_event.clear()
    background_task = 'alignment'
    close_after_alignment = False
    block_tabs(True)

    for index, sample in enumerate(samples_list):
        if align_sample_pair(reference_sample, sample, [progress_bar_info[0], progress_bar_info[1], progress_bar_info[2], index+1, len(samples_list)]) == "error":
            background_task = None
            
            # The window was closed during the alignment, so close it now that the alignment stopped
            if close_after_alignment:
                mzml_window.after(0, lambda: quit_mzml_window(mzml_window, from_GG))
                return
            
            progress_bar_info[0]['value'] = 0
            progress_bar_info[1].config(text=f"")
            progress_bar_info[2].config(text=f"{0}% Done")
            block_tabs(False)
            
            for name, button in buttons.items():
                if name == 'align' or name == 'add' or name == 'remove' or name == 'ref' or name == 'auto' or name == 'param' or name == 'close':
                    button.config(state="normal")
                if name == 'cancel':
                    button.config(state="disabled")

            return

    background_task = None
    
    # The window was closed right as the last sample finished aligning
    if close_after_alignment:
        mzml_window.after(0, lambda: quit_mzml_window(mzml_window, from_GG))
        return

    for name, button in buttons.items():
        if name == 'align' or name == 'add' or name == 'remove' or name == 'ref' or name == 'auto' or name == 'param' or name == 'close':
            button.config(state="normal")
        if name == 'cancel':
            button.config(state="disabled")

    progress_bar_info[0]['value'] = 100
    progress_bar_info[1].config(text=f"")
    progress_bar_info[2].config(text=f"{100}% Done")
    block_tabs(False)
    alignment_done(from_GG=from_GG)
    
def alignment_done(from_GG=False):
    '''
    '''
    def close_ad():
        alignment_done_window.destroy()

    def on_closing():
        '''This function is used to remove the function from the Close Window button (x button).'''
        return
        
    alignment_done_window = tk.Toplevel()
    # alignment_done_window.attributes("-topmost", True)
    alignment_done_window.withdraw()
    alignment_done_window.title("Alignment Done!")
    alignment_done_window.iconphoto(False, icon)
    alignment_done_window.resizable(False, False)
    alignment_done_window.grab_set()
    alignment_done_window.protocol("WM_DELETE_WINDOW", on_closing)
    
    calibrating_spectra_window_label = ttk.Label(alignment_done_window, text=f"Alignment Done!{"\nLoad the new files into GlycoGenius to analyze them." if from_GG else ""}", font=("Segoe UI", fontsize))
    calibrating_spectra_window_label.pack(pady=10, padx=70)
    
    # ok button
    ok_button = ttk.Button(alignment_done_window, text="Ok", style="small_button_sfw_style1.TButton", command=close_ad)
    ok_button.pack(pady=10, padx=70)
    
    alignment_done_window.update_idletasks()
    alignment_done_window.deiconify()
    window_width = alignment_done_window.winfo_width()
    window_height = alignment_done_window.winfo_height()
    screen_width = alignment_done_window.winfo_screenwidth()
    screen_height = alignment_done_window.winfo_screenheight()
    x_position = (screen_width - window_width) // 2
    y_position = (screen_height - window_height) // 2
    alignment_done_window.geometry(f"{window_width}x{window_height}+{x_position}+{y_position}")
        
def open_file_dialog_add_alignment_file(files_list):
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
    if longest_len*7 > 428:
        files_list.column("#0", width=longest_len*7)
        
    file_dialog.destroy()
    mzml_window.grab_set()
            
def remove_selected_item_alignment_file(files_list):
    selected_items = files_list.selection()  # Get the selected item(s)
    for item in selected_items:
        files_list.delete(item)

def block_tabs(state = False):
    if state:
        function_tabs.bind("<Button-1>", lambda e: "break")
    else:
        function_tabs.unbind("<Button-1>")

def mzml_window_start(from_GG=False, to_lift=False, change_sample=None):
    ''' from_GG: 
        0.  combobox (contains file name), 
        1.  samples path list, 
        2.  library, 
        3.  [spectra ax, spectra canvas, scaling_dropdown, rt_label], 
        4.  processed_data (with BPC and such)
    '''
    global current_dir, ico_image, mzml_window, last_selection, calibrants_list, mztol_entry, min_int_entry, ax_calibrate, canvas_calibrate, fig_calibrate, stats_label_SD, stats_label_fitr2, stats_label_eq, fit_combobox, trimmer_fig, trimmer_ax, trimmer_canvas, rt_bpc, function_tabs
    
    current_dir = pathlib.Path(__file__).parent.resolve()
    
    if change_sample != None:
        last_calibrants = last_selection.get(change_sample)
        for item in calibrants_list.get_children():
            calibrants_list.delete(item)
        
        ax_calibrate.cla()
        ax_calibrate.set_xlabel('m/z', fontsize=8)
        ax_calibrate.set_title('m/z error', fontsize=8)
        canvas_calibrate.draw()
        
        stats_label_SD.config(text="")
        stats_label_fitr2.config(text="")
        stats_label_eq.config(text="")
        
        if last_calibrants != None:
            for row in last_calibrants:
                calibrants_list.insert("", "end", values=row)
            if from_GG == 'start':
                edit_calibrant_list(calibrants_list, f"{change_sample}_{float(mztol_entry.get())}_{float(min_int_entry.get())}", [ax_calibrate, canvas_calibrate, fig_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()), True)
            else:
                edit_calibrant_list(calibrants_list, f"{change_sample}_{float(mztol_entry.get())}_{float(min_int_entry.get())}", [ax_calibrate, canvas_calibrate, fig_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()))
              
        plot_sample_bpc(trimmer_fig, trimmer_ax, trimmer_canvas, (rt_bpc[change_sample][0], rt_bpc[change_sample][1]), rt_bpc[change_sample][2])
        return
    
    if to_lift:
        mzml_window.deiconify()
        mzml_window.lift()
        return
        
    if not from_GG:
        library = ''
    else:
        library = from_GG[2]
    
    # Load icon to be used in the program
    ico_path = os.path.join(current_dir, "Assets/gg_icon.ico")
    ico_image = Image.open(ico_path)
    
    if not from_GG:
        mzml_window = tk.Tk()
    else:
        mzml_window = tk.Toplevel()
    mzml_window.withdraw()
    mzml_window.title("Spectra File Editor")
    if from_GG != False:
        mzml_window.title(f"Spectra File Editor - {from_GG[0].get()}")
    global icon
    icon = ImageTk.PhotoImage(ico_image)
    mzml_window.resizable(False, False)
    mzml_window.iconphoto(False, icon)
    mzml_window.grid_columnconfigure(0, weight=1)
    mzml_window.grid_columnconfigure(1, weight=1)
    mzml_window.grid_rowconfigure(0, weight=1)
    mzml_window.grid_rowconfigure(1, weight=1)
    mzml_window.grid_rowconfigure(2, weight=1)
    mzml_window.grid_rowconfigure(3, weight=1)
    mzml_window.grid_rowconfigure(4, weight=1)
    mzml_window.grid_rowconfigure(5, weight=1)
    mzml_window.protocol("WM_DELETE_WINDOW", lambda: quit_mzml_window(mzml_window, from_GG))
    mzml_window.geometry("580x620")
    
    # Default font size
    global fontsize
    fontsize = 10
    
    # Button style
    small_button_sfw_style1 = ttk.Style().configure("small_button_sfw_style1.TButton", font=("Segoe UI", fontsize), relief="raised", padding = (0, 0), justify="center")
    
    symbol_button_style1 = ttk.Style().configure("symbol_button_style1.TButton", font=("Segoe UI", fontsize), relief="raised", padding = (0, 0), justify="center", width=5)
    
    # Editing tabs framework
    function_tabs = ttk.Notebook(mzml_window)
    function_tabs.grid(row=0, column=0, sticky='nsew')
    
    # close button
    close_button = ttk.Button(mzml_window, text="Close", style="small_button_sfw_style1.TButton", command=lambda: quit_mzml_window(mzml_window, from_GG))
    close_button.grid(row=1, column=0, padx=(10,10), pady=(0,10), sticky="nse")
    
    # Calibration frame
    calibration_frame = tk.Frame(function_tabs, bd=0, relief="flat")
    calibration_frame.pack(fill=tk.BOTH, expand=True)
    function_tabs.add(calibration_frame, text="Calibration")
    
    # Calibration tab components
    
    # raw mz add
    mz_label = ttk.Label(calibration_frame, text='m/z:', font=("Segoe UI", fontsize))
    mz_label.grid(row=0, column=0, padx=(10, 10), pady=(10, 10), sticky='nw')
    mz_entry = ttk.Entry(calibration_frame, width=10)
    mz_entry.grid(row=0, column=0, padx=(40, 10), pady=(10, 10), sticky='nw')
    ToolTip(mz_label, "Add a single specific m/z value to the calibrants' list.")
    ToolTip(mz_entry, "Add a single specific m/z value to the calibrants' list.")
    
    # standard add
    standard_combobox_options = list(standard_calibrants.keys())
    standard_label = ttk.Label(calibration_frame, text='Standard:', font=("Segoe UI", fontsize))
    standard_label.grid(row=0, column=0, padx=(110, 10), pady=(10, 10), sticky='nw')
    standard_combobox = ttk.Combobox(calibration_frame, state="readonly", values=standard_combobox_options, width=15)
    if from_GG:
        standard_combobox.grid(row=0, column=0, padx=(170, 10), pady=(10, 10), sticky='nw')
    else:
        standard_combobox.grid(row=0, column=0, padx=(170, 180), pady=(10, 10), sticky='nw')
    ToolTip(standard_label, "Choose from a list of standard calibrants.")
    ToolTip(standard_combobox, "Choose from a list of standard calibrants.")
    
    # from library add
    fromlibrary_combobox_options = ['']
    if library != '':
        global glycans
        glycans = load_library_data(library)
        fromlibrary_combobox_options = fromlibrary_combobox_options+list(glycans.keys())
    fromlibrary_label = ttk.Label(calibration_frame, text='From Library:', font=("Segoe UI", fontsize))
    fromlibrary_combobox = ttk.Combobox(calibration_frame, state="readonly", values=fromlibrary_combobox_options, width=10)
    if from_GG:
        fromlibrary_label.grid(row=0, column=0, padx=(285, 10), pady=(10, 10), sticky='nw')
        fromlibrary_combobox.grid(row=0, column=0, padx=(365, 30), pady=(10, 10), sticky='nw')
    ToolTip(fromlibrary_label, "Choose a glycan from your glycans' library to use as a calibrant.")
    ToolTip(fromlibrary_combobox, "Choose a glycan from your glycans' library to use as a calibrant.")
    
    # add button
    if not from_GG:
        add_button = ttk.Button(calibration_frame, text="    Add Calibrant    ", style="small_button_sfw_style1.TButton", command=lambda: add_calibrant(calibrants_list, file_label, mz_entry, standard_combobox, fromlibrary_combobox, library))
    else:
        add_button = ttk.Button(calibration_frame, text="    Add Calibrant    ", style="small_button_sfw_style1.TButton", command=lambda: add_calibrant(calibrants_list, [from_GG[0], from_GG[1]], mz_entry, standard_combobox, fromlibrary_combobox, library))
    add_button.grid(row=0, column=0, columnspan=2, padx=(10, 10), pady=(8, 10), sticky="ne")
    ToolTip(add_button, "Add the selected calibrant to the list.")
    
    # calibrants list
    calibrants_list_scrollbar = tk.Scrollbar(calibration_frame, orient=tk.VERTICAL)
    calibrants_list = ttk.Treeview(calibration_frame, columns=("Calibrant", "Ref. m/z", "Found m/z", "PPM error", "m/z error", "Peak Time"), height=10, yscrollcommand=calibrants_list_scrollbar.set, show='headings')
    
    calibrants_list_columns = ["Calibrant", "Ref. m/z", "Found m/z", "PPM error", "m/z error", "Peak Time"]
    for col in calibrants_list_columns:
        calibrants_list.heading(col, text=col, command=lambda _col=col: calibrants_list_sort(calibrants_list, _col, False))
        
    calibrants_list.column("Calibrant", width=10)
    calibrants_list.column("Ref. m/z", width=5)
    calibrants_list.column("Found m/z", width=5)
    calibrants_list.column("PPM error", width=1)
    calibrants_list.column("m/z error", width=1)
    calibrants_list.column("Peak Time", width=1)
            
    calibrants_list_scrollbar.config(command=calibrants_list.yview, width=10)
    calibrants_list.grid(row=1, column=0, padx=(10, 10), pady=(0, 0), sticky="nsew")
    calibrants_list_scrollbar.grid(row=1, column=0, pady=(0, 0), sticky="nse")
    
    if from_GG != False:
        calibrants_list.bind("<ButtonRelease-1>", lambda event: click_treeview(event, calibrants_list, [from_GG[0], from_GG[1]], from_GG[3], float(mztol_entry.get()), float(min_int_entry.get())))
    
    # remove selected button
    if not from_GG:
        remove_button = ttk.Button(calibration_frame, text="Remove\nSelected", style="small_button_sfw_style1.TButton", command=lambda: remove_selected_item(calibrants_list, file_label, [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()), float(min_int_entry.get())))
        remove_button.grid(row=1, column=1, padx=(10, 10), pady=(0, 10), sticky="new")
    else:
        remove_button = ttk.Button(calibration_frame, text="Remove\nSelected", style="small_button_sfw_style1.TButton", command=lambda: remove_selected_item(calibrants_list, [from_GG[0], from_GG[1]], [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()), float(min_int_entry.get())))
        remove_button.grid(row=1, column=1, padx=(10, 10), pady=(0, 10), sticky="new")
    ToolTip(remove_button, "Remove the selected calibrant(s) from the list.")
    
    # remove all button
    if not from_GG:
        remove_all_button = ttk.Button(calibration_frame, text="Remove\nAll", style="small_button_sfw_style1.TButton", command=lambda: remove_all_items(calibrants_list, file_label, [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq]))
    else:
        remove_all_button = ttk.Button(calibration_frame, text="Remove\nAll", style="small_button_sfw_style1.TButton", command=lambda: remove_all_items(calibrants_list, [from_GG[0], from_GG[1]], [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq]))
    remove_all_button.grid(row=1, column=1, padx=(10, 10), pady=(45, 10), sticky="new")
    ToolTip(remove_all_button, "Remove all the calibrants from the list.")
    
    # sample navigation buttons
    if from_GG != False:
        sample_back_button = ttk.Button(calibration_frame, text="◄", style="symbol_button_style1.TButton", command=lambda: switch_sample_button('back', from_GG[0], from_GG[1], 'calibration', [calibrants_list, [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()), float(min_int_entry.get())]))
        sample_back_button.grid(row=1, column=1, padx=(10, 10), pady=(10, 67), sticky="sw")
        ToolTip(sample_back_button, "Go back one sample on your loaded samples list.")
        
        sample_forward_button = ttk.Button(calibration_frame, text="►", style="symbol_button_style1.TButton", command=lambda: switch_sample_button('forward', from_GG[0], from_GG[1], 'calibration', [calibrants_list, [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()), float(min_int_entry.get())]))
        sample_forward_button.grid(row=1, column=1, padx=(10, 10), pady=(10, 67), sticky="se")
        ToolTip(sample_forward_button, "Move to the next sample on your loaded samples list.")
        
        search_all_samples_button = ttk.Button(calibration_frame, text="Search\nAll Samples", style="small_button_sfw_style1.TButton", command=lambda: finding_calibrants(from_GG[0], calibrants_list, float(mztol_entry.get()), float(min_int_entry.get()), [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, from_GG[1], draw = False))
        search_all_samples_button.grid(row=1, column=1, padx=(10, 10), pady=(10, 25), sticky="sew")
        ToolTip(search_all_samples_button, "Search for the listed calibrants in all loaded sample files.")
    
    # search button
    if not from_GG:
        search_button = ttk.Button(calibration_frame, text="Search", style="small_button_sfw_style1.TButton", command=lambda: finding_calibrants(file_label, calibrants_list, float(mztol_entry.get()), float(min_int_entry.get()), [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox))
        search_button.grid(row=1, column=1, padx=(10, 10), pady=(10, 0), sticky="sew")
    else:
        search_button = ttk.Button(calibration_frame, text="Search", style="small_button_sfw_style1.TButton", command=lambda: finding_calibrants(from_GG[0], calibrants_list, float(mztol_entry.get()), float(min_int_entry.get()), [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, from_GG[1]))
        search_button.grid(row=1, column=1, padx=(10, 10), pady=(10, 0), sticky="sew")
    ToolTip(search_button, "Search for the listed calibrants in the sample.")
    
    # calibration status
    stats_label_SD = ttk.Label(calibration_frame, text='', font=("Segoe UI", fontsize))
    stats_label_SD.grid(row=2, column=0, columnspan=2, padx=(10, 10), pady=(0, 0), sticky='nw')
    ToolTip(stats_label_SD, "Standard deviation of the m/z error.")
    
    stats_label_fitr2 = ttk.Label(calibration_frame, text='', font=("Segoe UI", fontsize))
    stats_label_fitr2.grid(row=2, column=0, columnspan=2, padx=(180, 10), pady=(0, 0), sticky='nw')
    ToolTip(stats_label_fitr2, "R² value of the equation fitted to the datapoints.")
    
    stats_label_eq = ttk.Label(calibration_frame, text='', font=("Segoe UI", fontsize))
    stats_label_eq.grid(row=2, column=0, columnspan=2, padx=(300, 10), pady=(0, 0), sticky='nw')
    ToolTip(stats_label_eq, "Equation of the curve fitted to the datapoints, used in the calibration.")
    
    # plot area
    plot_frame = tk.Frame(calibration_frame)
    plot_frame.grid(row=3, column=0, columnspan=2, padx=(10, 10), pady=(0, 10), sticky='nwse')
    fig_calibrate = plt.figure(figsize=(0, 2))
    ax_calibrate = fig_calibrate.add_subplot(111)
    canvas_calibrate = FigureCanvasTkAgg(fig_calibrate, master=plot_frame)
    canvas_calibrate.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    ax_calibrate.set_position([0.07, 0.17, 0.91, 0.73])
    ax_calibrate.tick_params(axis='both', which='major', labelsize=6)
    ax_calibrate.set_xlabel('m/z', fontsize=8)
    ax_calibrate.set_title('m/z error', fontsize=8)
    
    # fit selection
    fit_combobox_options = ['Linear', 'Quadratic']
    fit_label = ttk.Label(calibration_frame, text='Fit:', font=("Segoe UI", fontsize))
    fit_label.grid(row=4, column=0, padx=(10, 10), pady=(0, 10), sticky='nw')
    fit_combobox = ttk.Combobox(calibration_frame, state="readonly", values=fit_combobox_options, width=7)
    fit_combobox.grid(row=4, column=0, padx=(35, 10), pady=(0, 10), sticky='nw')
    fit_combobox.set('Linear')
    if not from_GG:
        fit_combobox.bind("<<ComboboxSelected>>", lambda event: on_fit_selected(event, calibrants_list, file_label, [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()), float(min_int_entry.get())))
    else:
        fit_combobox.bind("<<ComboboxSelected>>", lambda event: on_fit_selected(event, calibrants_list, [from_GG[0], from_GG[1]], [ax_calibrate, canvas_calibrate], [stats_label_SD, stats_label_fitr2, stats_label_eq], fit_combobox, float(mztol_entry.get()), float(min_int_entry.get())))
    ToolTip(fit_label, "The polynomial degree of the equation fitted to the datapoints.")
    ToolTip(fit_combobox, "The polynomial degree of the equation fitted to the datapoints.")
    
    # search m/z tolerance
    mztol_label = ttk.Label(calibration_frame, text='Search m/z Tolerance:', font=("Segoe UI", fontsize))
    mztol_label.grid(row=4, column=0, columnspan=2, padx=(115, 10), pady=(0, 10), sticky='nw')
    mztol_entry = ttk.Entry(calibration_frame, width=10)
    mztol_entry.grid(row=4, column=0, columnspan=2, padx=(250, 10), pady=(0, 10), sticky='nw')
    mztol_entry.insert(0, "0.1")
    ToolTip(mztol_label, "The tolerance in m/z used for searching for calibrants in the spectra file.")
    ToolTip(mztol_entry, "The tolerance in m/z used for searching for calibrants in the spectra file.")
    
    # intensity threshold
    min_int_label = ttk.Label(calibration_frame, text='Search Intensity Threshold:', font=("Segoe UI", fontsize))
    min_int_label.grid(row=4, column=0, columnspan=2, padx=(330, 10), pady=(0, 10), sticky='nw')
    min_int_entry = ttk.Entry(calibration_frame, width=10)
    min_int_entry.grid(row=4, column=0, columnspan=2, padx=(490, 10), pady=(0, 10), sticky='nw')
    min_int_entry.insert(0, "500")
    ToolTip(min_int_label, "The intensity threshold for peaks that are found within the spectra file that corresponds to calibrants on the list.")
    ToolTip(min_int_entry, "The intensity threshold for peaks that are found within the spectra file that corresponds to calibrants on the list.")
    
    # load file button
    load_button = ttk.Button(calibration_frame, text="Load File", style="small_button_sfw_style1.TButton", command=lambda: load_file([file_label, file_label_trimming], mzml_window, [trimmer_fig, trimmer_ax, trimmer_canvas]))
    file_label = ttk.Label(calibration_frame, text=f'File: {file_path}', font=("Segoe UI", fontsize-1))
    if not from_GG:
        load_button.grid(row=5, column=0, padx=(10, 10), pady=(0, 0), sticky="sw")
        file_label.grid(row=5, column=0, padx=(100, 10), pady=(0, 2), sticky="sw")
        file_name = file_path.split("/")[-1]
        ToolTip(load_button, "Select a spectra file to load.")
        ToolTip(file_label, f"{file_name}")
    
    # calibrate button
    if not from_GG:
        calibrate_button = ttk.Button(calibration_frame, text="    Calibrate File    ", style="small_button_sfw_style1.TButton", command=lambda: calibrating_spectra(file_label, 'calibration', calibrants_list, [stats_label_SD, stats_label_fitr2, stats_label_eq], [ax_calibrate, canvas_calibrate, fig_calibrate]))
        calibrate_button.grid(row=5, column=0, columnspan=2, padx=(10, 10), pady=(0, 0), sticky="se")
    else:
        calibrate_button = ttk.Button(calibration_frame, text="    Calibrate File    ", style="small_button_sfw_style1.TButton", command=lambda: calibrating_spectra(from_GG[0], 'calibration', calibrants_list, [stats_label_SD, stats_label_fitr2, stats_label_eq], [ax_calibrate, canvas_calibrate, fig_calibrate], from_GG[1]))
        calibrate_button.grid(row=5, column=0, columnspan=2, padx=(10, 10), pady=(0, 0), sticky="se")
    ToolTip(calibrate_button, "Start the calibration using the selected calibrants derived equation.")
    
    ########################################################################################################
    #####################################          Trimming Tab          ###################################
    ########################################################################################################
    
    if from_GG != False:
        # Trimming frame
        rt_bpc = {key:[value['rt_array'], value['bpc'], value['time_unit']] for key, value in from_GG[4].items()}
    
    trimming_frame = tk.Frame(function_tabs, bd=0, relief="flat")
    trimming_frame.pack(fill=tk.BOTH, expand=True)
    function_tabs.add(trimming_frame, text="Trimming")
    
    trimming_frame.grid_rowconfigure(0, weight=2)
    trimming_frame.grid_rowconfigure(4, weight=1)
    trimming_frame.grid_columnconfigure(1, weight=1)
    
    trimmer_plot_frame = tk.Frame(trimming_frame)
    trimmer_plot_frame.grid(row=0, column=0, columnspan=3, padx=10, pady=10, sticky="nsew")
    
    trimmer_fig = plt.figure(figsize=(0,0))
    trimmer_ax = trimmer_fig.add_subplot(111)
    trimmer_canvas = FigureCanvasTkAgg(trimmer_fig, master=trimmer_plot_frame)
    trimmer_canvas.get_tk_widget().pack(fill = tk.BOTH, expand=True)
    
    sliders_label = ttk.Label(trimming_frame, text=f'RT/MT Range:', font=("Segoe UI", fontsize))
    sliders_label.grid(row=1, column=0, columnspan=3, padx=10, pady=(0, 10), sticky="ew")
    
    min_slider_label = ttk.Label(trimming_frame, text=f'Min:', font=("Segoe UI", fontsize))
    min_slider_label.grid(row=2, column=0, padx=10, pady=(0, 10), sticky="w")
    
    global slider_min
    slider_min_value = tk.DoubleVar()
    slider_min = ttk.Scale(trimming_frame, from_=0, to=trimmer_ax.get_xlim()[1],  orient='horizontal', variable=slider_min_value, command=lambda value:update_min_rt(value, slider_min_value, slider_max_value, selected_rt_range, trimmer_canvas, min_value_entry))
    slider_min.grid(row=2, column=0, columnspan=2, padx=(50, 10), pady=(0, 10), sticky="ew")
    
    min_value_entry = ttk.Spinbox(trimming_frame, width = 5, from_=0, to=trimmer_ax.get_xlim()[1], increment=0.1)
    min_value_entry.grid(row=2, column=2, padx=(10, 45), pady=(0, 10), sticky="e")
    min_value_entry.insert(0, 0.0)
    min_value_entry.bind("<KeyRelease>", lambda event:update_min_rt(min_value_entry, slider_min_value, slider_max_value, selected_rt_range, trimmer_canvas, min_value_entry, True, True))
    min_value_entry.bind("<<Increment>>", lambda event:update_min_rt(min_value_entry, slider_min_value, slider_max_value, selected_rt_range, trimmer_canvas, min_value_entry, True))
    min_value_entry.bind("<<Decrement>>", lambda event:update_min_rt(min_value_entry, slider_min_value, slider_max_value, selected_rt_range, trimmer_canvas, min_value_entry, True))
    
    min_value_label = ttk.Label(trimming_frame, text=f'min', font=("Segoe UI", fontsize))
    min_value_label.grid(row=2, column=2, padx=10, pady=(0, 10), sticky="e")
    
    max_sliders_label = ttk.Label(trimming_frame, text=f'Max:', font=("Segoe UI", fontsize))
    max_sliders_label.grid(row=3, column=0, padx=10, pady=(0, 10), sticky="w")
    
    global slider_max
    slider_max_value = tk.DoubleVar()
    slider_max = ttk.Scale(trimming_frame, from_=0, to=trimmer_ax.get_xlim()[1],  orient='horizontal', variable=slider_max_value, command=lambda value:update_max_rt(value, slider_max_value, slider_min_value, selected_rt_range, trimmer_canvas, max_value_entry))
    slider_max.grid(row=3, column=0, columnspan=2, padx=(50, 10), pady=(0, 10), sticky="ew")
    slider_max_value.set(trimmer_ax.get_xlim()[1])
    
    max_value_entry = ttk.Spinbox(trimming_frame, width = 5, from_=0, to=trimmer_ax.get_xlim()[1], increment=0.1)
    max_value_entry.grid(row=3, column=2, padx=(10, 45), pady=(0, 10), sticky="e")
    max_value_entry.insert(0, f"{trimmer_ax.get_xlim()[1]:.1f}")
    max_value_entry.bind("<KeyRelease>", lambda event:update_max_rt(max_value_entry, slider_max_value, slider_min_value, selected_rt_range, trimmer_canvas, max_value_entry, True, True))
    max_value_entry.bind("<<Increment>>", lambda event:update_max_rt(max_value_entry, slider_max_value, slider_min_value, selected_rt_range, trimmer_canvas, max_value_entry, True))
    max_value_entry.bind("<<Decrement>>", lambda event:update_max_rt(max_value_entry, slider_max_value, slider_min_value, selected_rt_range, trimmer_canvas, max_value_entry, True))
    
    max_value_label = ttk.Label(trimming_frame, text=f'min', font=("Segoe UI", fontsize))
    max_value_label.grid(row=3, column=2, padx=10, pady=(0, 10), sticky="e")
    
    if not from_GG:
        pass
    else:
        plot_sample_bpc(trimmer_fig, trimmer_ax, trimmer_canvas, (from_GG[4][from_GG[0].get()]['rt_array'], from_GG[4][from_GG[0].get()]['bpc']), from_GG[4][from_GG[0].get()]['time_unit'])
    
    # load file button
    load_button_trimming = ttk.Button(trimming_frame, text="Load File", style="small_button_sfw_style1.TButton", command=lambda: load_file([file_label, file_label_trimming], mzml_window, [trimmer_fig, trimmer_ax, trimmer_canvas]))
    file_label_trimming = ttk.Label(trimming_frame, text=f'File: {file_path}', font=("Segoe UI", fontsize-1))
    if from_GG != False:
        sample_back_button_trimming = ttk.Button(trimming_frame, text="◄", style="symbol_button_style1.TButton", command=lambda: switch_sample_button('back', from_GG[0], from_GG[1], "trimming", [trimmer_fig, trimmer_ax, trimmer_canvas, from_GG[4]]))
        sample_back_button_trimming.grid(row=5, column=0, padx=(10, 10), pady=(0, 2), sticky="sw")
        ToolTip(sample_back_button_trimming, "Go back one sample on your loaded samples list.")
        
        sample_forward_button_trimming = ttk.Button(trimming_frame, text="►", style="symbol_button_style1.TButton", command=lambda: switch_sample_button('forward', from_GG[0], from_GG[1], "trimming", [trimmer_fig, trimmer_ax, trimmer_canvas, from_GG[4]]))
        sample_forward_button_trimming.grid(row=5, column=0, padx=(55, 10), pady=(0, 2), sticky="sw")
        ToolTip(sample_forward_button_trimming, "Move to the next sample on your loaded samples list.")
    else:
        load_button_trimming.grid(row=5, column=0, padx=(10, 10), pady=(0, 2), sticky="sw")
        ToolTip(load_button_trimming, "Select a spectra file to load.")
        file_label_trimming.grid(row=5, column=0, padx=(100, 10), pady=(0, 2), sticky="sw")
        file_name_trimming = file_path.split("/")[-1]
        ToolTip(file_label_trimming, f"{file_name_trimming}")
        
    if not from_GG:
        trim_sample_button = ttk.Button(trimming_frame, text="Trim Sample", style="small_button_sfw_style1.TButton", command=lambda: trimming_spectra(file_label_trimming, 'trimming', trimming_range))
        trim_sample_button.grid(row=5, column=2, padx=10, pady=(0, 2), sticky="se")
    else:
        trim_sample_button = ttk.Button(trimming_frame, text="Trim Sample", style="small_button_sfw_style1.TButton", command=lambda: trimming_spectra(from_GG[0], 'trimming', trimming_range, from_GG[1]))
        trim_sample_button.grid(row=5, column=2, padx=10, pady=(0, 2), sticky="se")
    
    ########################################################################################################
    #####################################          Alignment Tab         ###################################
    ########################################################################################################
    
    align_frame = tk.Frame(function_tabs, bd=0, relief="flat")
    align_frame.pack(fill=tk.BOTH, expand=True)
    function_tabs.add(align_frame, text="Alignment")
    
    # Add Ttk Styles
    files_list_style = ttk.Style().configure("files_list.Treeview", bg="white", fg="black", font=("Segoe UI", fontsize-1))
    
    # Alignment window widgets
    main_label = ttk.Label(align_frame, text="Select files to align:", font=("Segoe UI", fontsize))
    main_label.grid(row=0, column=0, padx=(10, 10), pady=(10, 10), sticky="nw")
    
    files_list_frame = tk.Frame(align_frame, bd=0, relief="flat", width=430, height=250)
    files_list_frame.pack_propagate(0)
    files_list_frame.grid(row=0, column=0, padx=(10,20), pady=(30, 10), sticky="nsew")
    
    files_list_scrollbar = tk.Scrollbar(align_frame, orient=tk.HORIZONTAL)
    files_list_scrollbar_v = tk.Scrollbar(align_frame, orient=tk.VERTICAL)
    files_list = ttk.Treeview(files_list_frame, height=20, style="files_list.Treeview", show="tree", xscrollcommand=files_list_scrollbar.set, yscrollcommand=files_list_scrollbar_v.set)
    files_list.column("#0", stretch=False, width=428)
    files_list_scrollbar.config(command=files_list.xview, width=10)
    files_list_scrollbar_v.config(command=files_list.yview, width=10)
    files_list.pack(fill="both", expand=True, padx=0, pady=0)
    files_list_scrollbar.grid(row=0, column=0, padx=(10, 20), pady=0, sticky="sew")
    files_list_scrollbar_v.grid(row=0, column=0, padx=(0, 10), pady=(30, 10), sticky="nse")
    files_list.tag_configure("reference_sample", background="lightblue")
    
    if from_GG:
        # Load the samples_list from GG if possible
        for item in from_GG[1]:
            files_list.insert('', 'end', text=item)
        longest_len = 0
        for i in files_list.get_children():
            if len(files_list.item(i, "text")) > longest_len:
                longest_len = len(files_list.item(i, "text"))
        if longest_len*7 > 428:
            files_list.column("#0", width=longest_len*7)
    
    add_file_button = ttk.Button(align_frame, text="Add File", style="small_button_sfw_style1.TButton", command=lambda: open_file_dialog_add_alignment_file(files_list))
    
    remove_selected_button = ttk.Button(align_frame, text="Remove Selected", style="small_button_sfw_style1.TButton", command=lambda: remove_selected_item_alignment_file(files_list))
    
    if not from_GG:
        add_file_button.grid(row=0, column=1, padx=(0,10), pady=(30,10), sticky="new")
        ToolTip(add_file_button, "Select one or more MzML or MzXML files to add to the list of files to align.")
        
        remove_selected_button.grid(row=0, column=1, padx=(0,10), pady=(60,10), sticky="new")
        ToolTip(remove_selected_button, "Remove the currently selected file from the list of files.")
    
    select_reference_button = ttk.Button(align_frame, text="Select as\nReference", style="small_button_sfw_style1.TButton", command=lambda: set_reference(files_list))
    select_reference_button.grid(row=0, column=1, padx=(0,10), pady=(90,10), sticky="new")
    ToolTip(select_reference_button, "Select file as reference for alignment.")
    
    automatic_button = ttk.Button(align_frame, text="Automatically Set Parameters", style="small_button_sfw_style1.TButton", command=lambda: threading.Thread(target=auto_parameters, args=(files_list, selected_parameters, {'cancel': cancel_align_button, 'align': align_sample_button, 'auto': automatic_button, 'param': check_parameters_button, 'close': close_button, 'add': add_file_button, 'remove': remove_selected_button, 'ref': select_reference_button}, from_GG)).start())
    automatic_button.grid(row=1, rowspan=2, column=0, padx=(10,10), pady=(15,0), sticky="nw")
    
    check_parameters_button = ttk.Button(align_frame, text="   Manually Set Parameters   ", style="small_button_sfw_style1.TButton", command=lambda: open_alignment_settings_window(selected_parameters))
    check_parameters_button.grid(row=1, rowspan=2, column=0, padx=(10,10), pady=(50,10), sticky="nw")
    
    selected_parameters_label = ttk.Label(align_frame, text="Selected Parameters:", font=("Segoe UI", fontsize))
    selected_parameters_label.grid(row=1, column=0, columnspan=2, padx=(230, 10), pady=(5, 0), sticky="nw")
    selected_parameters = ScrolledText(align_frame, width=35, height=10, wrap=tk.WORD)
    selected_parameters.grid(row=2, column=0, columnspan=2, padx=(10,10), pady=(0,10), sticky="e")
    selected_parameters.config(state="disabled")
    update_alignment_settings(selected_parameters)
    
    progress_bar_label = ttk.Label(align_frame, text="", font=("Segoe UI", fontsize))
    progress_bar_label.grid(row=2, column=0, padx=10, pady=(0, 2), sticky="sw")
    
    progress_bar_align = ttk.Progressbar(align_frame, orient="horizontal", mode="determinate")
    progress_bar_align.grid(row=3, column=0, columnspan=2, padx=10, pady=(0,0), sticky="ew")
    
    progress_bar_label2 = ttk.Label(align_frame, text="0% Done", font=("Segoe UI", fontsize))
    progress_bar_label2.grid(row=4, column=0, columnspan=2, padx=10, pady=(0,10))
    
    powered_by_label = ttk.Label(align_frame, text="Powered by PASTAQ", font=("Segoe UI", fontsize), cursor="hand2")
    powered_by_label.grid(row=4, column=0, padx=10, pady=(5, 2), sticky="sw")
    powered_by_label.bind("<Button-1>", lambda event: webbrowser.open("https://github.com/PASTAQ-MS/PASTAQ"))
    ToolTip(powered_by_label, "Go to PASTAQ's GitHub page.")
        
    cancel_align_button = ttk.Button(align_frame, text="Cancel Alignment", style="small_button_sfw_style1.TButton", command=lambda: (stop_run_event.set(), progress_bar_label.config(text="Cancelling...")))
    cancel_align_button.grid(row=4, column=0, columnspan=2, padx=(10, 100), pady=(5, 2), sticky="se")
    cancel_align_button.config(state="disabled")
        
    align_sample_button = ttk.Button(align_frame, text="Align Samples", style="small_button_sfw_style1.TButton", command=lambda: threading.Thread(target=align_samples, args=(files_list, {'cancel': cancel_align_button, 'align': align_sample_button, 'auto': automatic_button, 'param': check_parameters_button, 'close': close_button, 'add': add_file_button, 'remove': remove_selected_button, 'ref': select_reference_button}, [progress_bar_align, progress_bar_label, progress_bar_label2], from_GG)).start())
    align_sample_button.grid(row=4, column=1, padx=10, pady=(5, 2), sticky="se")
    
    mzml_window.deiconify()
    
    if from_GG != False and change_sample == None:
        mzml_window_start(from_GG = 'start', change_sample=from_GG[0].get())
            
    mzml_window.mainloop()
    
if __name__ == "__main__":
    mzml_window_start()