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

"""
```
The coordinates system to be used is like this:
   0   1    2    3    4    5    6    7   (y)
0
1                     X
2                     |
3                     |
4                     |
5                     V
6   
7
8
9
10
11
12

(x)
```

"X" indicates the starting drawing point. These coordinates draws the glycan from top to bottom on the canvas.

The image is later rotated 180 degrees to have the glycan reducing end at the bottom.

By default coordinates within the same line will have 2 sizes difference between each other.
"""

from PIL import Image, ImageDraw, ImageFont
from itertools import combinations_with_replacement, permutations
from re import split
import datetime
import pathlib
import dill
import math
import copy
import os

# The glycan colors following SNFG parameters.
glycans_colors = {
    'green' : '#00a651',
    'yellow' : '#ffd400',
    'blue' : '#0072bc',
    'red' : '#ed1c24',
    'purple' : '#a54399',
    'light blue' : '#8fcce9',
    'orange' : '#f47920'
}

# This registers the highest and lowest x and y coordinates to crop the image accordingly.
global highest_coordinate, lowest_coordinate
highest_coordinate = [0, 0]
lowest_coordinate = [99999, 99999]

global image_width, image_height
# Define image size and background color
image_width, image_height = 9999,9999
background_color = (0, 0, 0, 0)  # Transparent

"""Functions necessary for the figure generation:"""

def calculate_coordinates(x, y, size):
    '''Receives the coordinates as described in the header and translate it into canvas coordinates.

    Parameters
    ----------
    x : int
        The horizontal coordinate.

    y : int
        The vertical coordinate.

    size : int
        The size to scale the coordinates.

    Returns
    -------
    int, int
        Coordinates to draw on canvas.
    '''
    return (2*size)+(x*2*size), y*2*size

def form_to_comp(string):
    '''Separates a molecular formula or monosaccharides formula of glycans into a
    dictionary with each atom/monosaccharide as a key and its amount as value.

    Parameters
    ----------
    string : str
        A string in the form of C6O6N0H12 or H5N4S1F1G1.

    Returns
    -------
    counts : dict
        A dictionary with keys containing the monosaccharides/atoms letter(s) and values
        containing the amounts of each. ie. {"H": 5, "N": 4, "S": 1, "F": 1, "G": 1}.
    '''
    string = string.split("+")[0] #avoids getting phosphorylation and sulfation symbols
    counts = {}
    split_str = split('(\\d+)', string)
    negative = False
    for i_i, i in enumerate(split_str):
        if i != '' and i[-1] == '-':
            split_str[i_i] = i[:-1]
            negative = True
        if i_i%2 != 0 and i != '' and negative:
            split_str[i_i] = '-'+i
    if len(split_str)%2 != 0:
        split_str.append('1')
    for i in range(len(split_str)-1):
        if i%2 == 0:
            counts[split_str[i]] = int(split_str[i+1])
    if '' in counts:
        del counts['']
    return counts

def comp_to_formula(composition):
    '''Transforms a composition dictionary into string formula.

    Parameters
    ----------
    composition : dict
        Dictionary containing the composition of the molecule or glycan.

    Returns
    -------
    formula : string
        Formula of the atomic or monosaccharides composition in string form.
    '''
    formula = ''
    for i in composition:
        if composition[i] != 0:
            formula+=i+str(composition[i])
    return formula

def permute_comp_keys(input_comp):
    '''Permutes the composition dictionary key:value positions.

    Parameters
    ----------
    input_comp : dict
        A dictionary containing the compositions.

    Returns
    -------
    permuted_dicts : list
        A list of the possible dictionary permutations.
    '''
    # Extract items from the dictionary as a list of tuples
    items = list(input_comp.items())

    # Generate all permutations of the items
    permutations_dict = list(permutations(items))

    # Convert each permutation back into a dictionary
    permuted_dicts = [dict(permutation) for permutation in permutations_dict]

    return permuted_dicts

def unique_permutations(s):
    '''Calculates unique permutations of letters in a string.

    Parameters
    ----------
    s : string
        String to permutate letters.

    Returns
    -------
    perm_strings : list
        List containing permutated strings.
    '''
    # Generate all permutations
    perms = set(permutations(s))

    # Join tuples to form strings
    perm_strings = [''.join(p) for p in perms]

    return perm_strings

def replace_s_and_g_with_chars(str_list, s_char_string, g_char_string):
    '''Replaces 'S' with characters from s_char_string and 'G' with characters from g_char_string in the list of strings.'''

    # Create two iterators: one for replacing 'S' and one for replacing 'G'
    s_char_iter = iter(s_char_string)
    g_char_iter = iter(g_char_string)

    # Process each string in the list
    result = []
    for s in str_list:
        new_string = ''
        for char in s:
            if char == 'S':
                # Replace 'S' with the next character from s_char_iter
                new_string += next(s_char_iter)
            elif char == 'G':
                # Replace 'G' with the next character from g_char_iter
                new_string += next(g_char_iter)
            else:
                new_string += char
        result.append(new_string)

    return result

def draw_upside_down_text(image, text, coordinates, size, font, fill="black"):
    '''Writes upside down on an image so that it looks upside up when rotating.

    Parameters
    ----------
    image
        The target canvas to draw the text.

    text : string
        The text to write.

    coordinates : list
        The x and y coordinates for the text.

    size : int
        The size of the text.

    font : string
        A font or path to the font to be used in the writing.

    fill : str
        A string containing the color or the hex code for the color of the text.
    '''
    # Canvas width and height
    width = image.width
    height = image.height

    # Create a temporary image for the text
    temp_image = Image.new('RGBA', (width, height), (255, 255, 255, 0))  # Transparent background
    temp_draw = ImageDraw.Draw(temp_image)

    # Calculate upside down position
    coordinates = [width-coordinates[0]-size*0.75, height-coordinates[1]-size]

    # Draw the text on the temporary image
    temp_draw.text(coordinates, text, font=font, fill=fill)

    # Rotate the temporary image by 180 degrees
    rotated_text_image = temp_image.rotate(180, expand=True)

    # Paste the rotated text onto the main canvas
    image.paste(rotated_text_image, (0, 0), rotated_text_image)

def list_available_fonts():
    '''Lists the available fonts to be used by the script.

    Returns
    -------
    available_fonts
        A list of the available fonts.
    '''
    # Common font directories on different platforms
    font_dirs = [
        "/usr/share/fonts",            # Linux
        "/Library/Fonts",              # macOS
        "C:/Windows/Fonts"             # Windows
    ]

    available_fonts = []

    for font_dir in font_dirs:
        if os.path.exists(font_dir):
            for root, dirs, files in os.walk(font_dir):
                for file in files:
                    if file.endswith(".ttf") or file.endswith(".otf"):
                        available_fonts.append(os.path.join(root, file))

    return available_fonts

def diamond(target_canvas, center_x, center_y, size, fill_color="blue"):
    """Draws a diamond shape at the set coordinates, set size and fill color directly on the target canvas.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    center_x : int
        The horizontal coordinate for the drawing.

    center_y : int
        The vertical coordinate for the drawing.

    size : int
        The size of the geometric shape.

    fill_color : str
        A string containing the color or the hex code for the color of the shape.
    """
    # Calculate the four points of the diamond
    points = [
        center_x, center_y - size*0.6,  # Top point
        center_x + size*0.6, center_y,  # Right point
        center_x, center_y + size*0.6,  # Bottom point
        center_x - size*0.6, center_y   # Left point
    ]

    # Draw the diamond
    target_canvas.polygon(points, fill=fill_color, outline='black')

def triangle(target_canvas, center_x, center_y, size, fill_color, orientation=0):
    """Draws a triangle shape at the set coordinates, set size and fill color directly on the target canvas.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    center_x : int
        The horizontal coordinate for the drawing.

    center_y : int
        The vertical coordinate for the drawing.

    size : int
        The size of the geometric shape.

    fill_color : str
        A string containing the color or the hex code for the color of the shape.

    orientation : int
        The degree of rotation of the triangle.
    """
    # Calculate the height of the equilateral triangle
    height = (math.sqrt(3) / 2) * size*1.3

    # Calculate the vertices of the triangle relative to the center
    point1 = (center_x, center_y - 2/3 * height)  # Top vertex
    point2 = (center_x - size * 0.65, center_y + 1/3 * height)  # Bottom left vertex
    point3 = (center_x + size * 0.65, center_y + 1/3 * height)  # Bottom right vertex

    # Function to rotate a point around the center
    def rotate_point(x, y, angle, cx, cy):
        radians = math.radians(angle)
        cos_val = math.cos(radians)
        sin_val = math.sin(radians)
        x -= cx
        y -= cy
        new_x = x * cos_val - y * sin_val + cx
        new_y = x * sin_val + y * cos_val + cy
        return new_x, new_y

    # Rotate each point around the center by the specified orientation angle
    point1 = rotate_point(*point1, orientation, center_x, center_y-size*0.05)
    point2 = rotate_point(*point2, orientation, center_x, center_y-size*0.05)
    point3 = rotate_point(*point3, orientation, center_x, center_y-size*0.05)

    # Draw the triangle using polygon
    target_canvas.polygon((point1, point2, point3), fill=fill_color, outline='black')

def draw_glcnac(target_canvas, x_y_coordinates, size):
    '''Draws a GlcNAc on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0], x_y_coordinates[1], x_y_coordinates[0]+size, x_y_coordinates[1]+size]
    target_canvas.rectangle(shape_coordinates, fill=glycans_colors['blue'], outline='black')

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.05 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.05
    if x_y_coordinates[1]+size*1.05 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.05
    if x_y_coordinates[0]-size*0.05 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.05
    if x_y_coordinates[1]-size*0.05 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.05

def draw_galnac(target_canvas, x_y_coordinates, size):
    '''Draws a GalNAc on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0], x_y_coordinates[1], x_y_coordinates[0]+size, x_y_coordinates[1]+size]
    target_canvas.rectangle(shape_coordinates, fill=glycans_colors['yellow'], outline='black')

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.05 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.05
    if x_y_coordinates[1]+size*1.05 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.05
    if x_y_coordinates[0]-size*0.05 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.05
    if x_y_coordinates[1]-size*0.05 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.05

def draw_hexnac(target_canvas, x_y_coordinates, size):
    '''Draws a generic HexNAc on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0], x_y_coordinates[1], x_y_coordinates[0]+size, x_y_coordinates[1]+size]
    target_canvas.rectangle(shape_coordinates, fill='white', outline='black')

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.05 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.05
    if x_y_coordinates[1]+size*1.05 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.05
    if x_y_coordinates[0]-size*0.05 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.05
    if x_y_coordinates[1]-size*0.05 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.05

def draw_mannose(target_canvas, x_y_coordinates, size):
    '''Draws a Mannose on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0], x_y_coordinates[1], x_y_coordinates[0]+size, x_y_coordinates[1]+size]
    target_canvas.ellipse(shape_coordinates, fill=glycans_colors['green'], outline='black')

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.05 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.05
    if x_y_coordinates[1]+size*1.05 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.05
    if x_y_coordinates[0]-size*0.05 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.05
    if x_y_coordinates[1]-size*0.05 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.05

def draw_galactose(target_canvas, x_y_coordinates, size):
    '''Draws a Galactose on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0], x_y_coordinates[1], x_y_coordinates[0]+size, x_y_coordinates[1]+size]
    target_canvas.ellipse(shape_coordinates, fill=glycans_colors['yellow'], outline='black')

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.05 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.05
    if x_y_coordinates[1]+size*1.05 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.05
    if x_y_coordinates[0]-size*0.05 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.05
    if x_y_coordinates[1]-size*0.05 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.05

def draw_glucose(target_canvas, x_y_coordinates, size):
    '''Draws a Glucose on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0], x_y_coordinates[1], x_y_coordinates[0]+size, x_y_coordinates[1]+size]
    target_canvas.ellipse(shape_coordinates, fill=glycans_colors['blue'], outline='black')

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.05 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.05
    if x_y_coordinates[1]+size*1.05 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.05
    if x_y_coordinates[0]-size*0.05 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.05
    if x_y_coordinates[1]-size*0.05 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.05

def draw_hexose(target_canvas, x_y_coordinates, size):
    '''Draws a generic hexose on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0], x_y_coordinates[1], x_y_coordinates[0]+size, x_y_coordinates[1]+size]
    target_canvas.ellipse(shape_coordinates, fill='white', outline='black')

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.05 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.05
    if x_y_coordinates[1]+size*1.05 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.05
    if x_y_coordinates[0]-size*0.05 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.05
    if x_y_coordinates[1]-size*0.05 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.05

def draw_neu5ac(target_canvas, x_y_coordinates, size):
    '''Draws a Neu5Ac sialic acid on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0]+size*0.5, x_y_coordinates[1]+size*0.5]
    diamond(target_canvas, shape_coordinates[0], shape_coordinates[1], size, fill_color=glycans_colors['purple'])

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.15 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.15
    if x_y_coordinates[1]+size*1.15 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.15
    if x_y_coordinates[0]-size*0.15 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.15
    if x_y_coordinates[1]-size*0.15 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.15

def draw_neu5gc(target_canvas, x_y_coordinates, size):
    '''Draws a Neu5Gc sialic acid on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0]+size*0.5, x_y_coordinates[1]+size*0.5]
    diamond(target_canvas, shape_coordinates[0], shape_coordinates[1], size, fill_color=glycans_colors['light blue'])

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.15 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.15
    if x_y_coordinates[1]+size*1.15 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.15
    if x_y_coordinates[0]-size*0.15 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.15
    if x_y_coordinates[1]-size*0.15 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.15

def draw_fucose(target_canvas, x_y_coordinates, size):
    '''Draws a Fucose on the target_canvas, at the target coordinates, with the given size.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    x_y_coordinates : list
        A list containing the x and y coordinates to draw the monosaccharide.

    size : int
        The size of the geometric shape.
    '''
    shape_coordinates = [x_y_coordinates[0]+size*0.5, x_y_coordinates[1]+size*0.5]
    triangle(target_canvas, shape_coordinates[0], shape_coordinates[1], size, fill_color=glycans_colors['red'], orientation = 180)

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if x_y_coordinates[0]+size*1.2 > highest_coordinate[0]:
        highest_coordinate[0] = x_y_coordinates[0]+size*1.2
    if x_y_coordinates[1]+size*1.2 > highest_coordinate[1]:
        highest_coordinate[1] = x_y_coordinates[1]+size*1.2
    if x_y_coordinates[0]-size*0.2 < lowest_coordinate[0]:
        lowest_coordinate[0] = x_y_coordinates[0]-size*0.2
    if x_y_coordinates[1]-size*0.2 < lowest_coordinate[1]:
        lowest_coordinate[1] = x_y_coordinates[1]-size*0.2

def draw_connecting_line(target_canvas, former_coord, new_coord, size, angled = False, side = False, o_glycan = False):
    '''Draws a connecting line between two pair of coordinates that might be more angled or not.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    former_coord : list
        A list containing the x and y coordinates of the former monosaccharide.

    new_coord : list
        A list containing the x and y coordinates of the new monosaccharide.

    angled : boolean, string
        Whether or not the connecting line is in an angle or not. Can be 'False', 'right' or 'left'.
    '''
    if angled == 'right':
        if o_glycan:
            coord = (former_coord[0]+size*0.1, former_coord[1]+size), (new_coord[0]+size*0.8, new_coord[1]+size*0.2)
            target_canvas.line(coord, fill="black", width=round(size/10))
        else:
            coord = (former_coord[0]+size*0.2, former_coord[1]+size*0.9), (new_coord[0]+size*0.8, new_coord[1]+size*0.2)
            target_canvas.line(coord, fill="black", width=round(size/10))
    elif angled =='left':
        if o_glycan:
            coord = (former_coord[0]+size*0.9, former_coord[1]+size), (new_coord[0]+size*0.2, new_coord[1]+size*0.2)
            target_canvas.line(coord, fill="black", width=round(size/10))
        else:
            coord = (former_coord[0]+size*0.8, former_coord[1]+size*0.9), (new_coord[0]+size*0.2, new_coord[1]+size*0.2)
            target_canvas.line(coord, fill="black", width=round(size/10))
    elif side == 'right':
        coord = (former_coord[0], former_coord[1]+size*0.5), (new_coord[0]+size*0.5, new_coord[1]+size*0.5)
        target_canvas.line(coord, fill="black", width=round(size/10))
    elif side =='left':
        coord = (former_coord[0]+size, former_coord[1]+size*0.5), (new_coord[0]+size*0.5, new_coord[1]+size*0.5)
        target_canvas.line(coord, fill="black", width=round(size/10))
    else:
        coord = (former_coord[0]+size/2, former_coord[1]+size), (new_coord[0]+size/2, new_coord[1])
        target_canvas.line(coord, fill="black", width=round(size/10))

def draw_bracket(target_canvas, center, half_width, size):
    '''Draw a bracket for generic antennas.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    center : list
        List of x and y coordinates of the center of the bracket.

    half_width : int
        Half width of the bracket, in GG coordinates.

    size : int
        The size of the geometric shape.
    '''
    # Scale the width from GG coordinates
    half_width = half_width*size*2

    # Calculate coordinates of sides
    left_side = (center[0]-half_width, center[1]+size*1.5)
    right_side = (center[0]+half_width+size, center[1]+size*1.5)

    # Draw horizontal line
    target_canvas.line([left_side, right_side], fill="black", width = round(size/15))

    # Draw left and right vertical lines
    target_canvas.line([left_side, (left_side[0], left_side[1]-size/5)], fill="black", width = round(size/15))
    target_canvas.line([right_side, (right_side[0], right_side[1]-size/5)], fill="black", width = round(size/15))

    # Updates the glycan coordinate in the canvas
    global highest_coordinate, lowest_coordinate
    if right_side[0]+size*0.2 > highest_coordinate[0]:
        highest_coordinate[0] = right_side[0]+size*0.2
    if left_side[0]-size*0.2 < lowest_coordinate[0]:
        lowest_coordinate[0] = left_side[0]-size*0.2

def draw_n_glycan_core(target_canvas, size, wide = [False, False], complete = True, starting_coordinates = [4, 1]):
    '''Draws an N-glycan chitobiose core at fixed coordinates on the target canvas.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    size : int
        The size of the geometric shape.

    complete : boolean
        Whether or not to draw the complete core or a core lacking one mannose.
    '''
    monosaccharides_coordinates = []

    # First GlcNAc
    gg_coordinates = starting_coordinates
    x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
    monosaccharides_coordinates.append(x_y_coordinates)
    draw_glcnac(target_canvas, x_y_coordinates, size)

    # Second GlcNac
    gg_coordinates[1] += 1
    x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
    draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size)
    monosaccharides_coordinates.append(x_y_coordinates)
    draw_glcnac(target_canvas, x_y_coordinates, size)

    # First Mannose, bifurcation
    gg_coordinates[1] += 1
    x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
    draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size)
    monosaccharides_coordinates.append(x_y_coordinates)
    draw_mannose(target_canvas, x_y_coordinates, size)

    # Second Mannose, angled to the right
    antenna_gg_coordinates = [gg_coordinates[0]-1, gg_coordinates[1]+1]
    if wide[1]:
        antenna_gg_coordinates[0] -= 1
    x_y_coordinates = calculate_coordinates(*antenna_gg_coordinates, size)
    draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled = 'right')
    draw_mannose(target_canvas, x_y_coordinates, size)

    if complete:
        # Third Mannose, angled to the left
        antenna_gg_coordinates = [gg_coordinates[0]+1, gg_coordinates[1]+1]
        if wide[0]:
            antenna_gg_coordinates[0] += 1
        x_y_coordinates = calculate_coordinates(*antenna_gg_coordinates, size)
        draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled = 'left')
        draw_mannose(target_canvas, x_y_coordinates, size)

def draw_o_glycan_core(target_canvas, size, wide = [False, False], starting_coordinates = [4, 1]):
    '''Draws an O-glycan core at fixed coordinates on the target canvas.

    Parameters
    ----------
    target_canvas
        The drawing target. Can be a tkinter canvas or a PIL imageDraw object.

    size : int
        The size of the geometric shape.

    core_number : int
        The o-glycan core.
    '''
    monosaccharides_coordinates = []

    gg_coordinates = starting_coordinates

    # First GalNAc
    x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
    monosaccharides_coordinates.append(x_y_coordinates)
    draw_galnac(target_canvas, x_y_coordinates, size)

def draw_antenna(target_canvas, size, monosaccharides, wide = False, angled = False, starting_coordinates = [4, 1], antenna_number = 0, galnac_start = False, o_glycan = False):
    '''Work in Progress. Need to implement Am and E sialic acids.
    '''
    monosaccharides_coordinates = []
    gg_coordinates = starting_coordinates
    monosaccharides_coordinates.append(calculate_coordinates(*gg_coordinates, size))
    original_angled = angled

    for index, mono in enumerate(monosaccharides):
        from_glcnac = False
        if index != 0:
            angled = False
        if mono == 'N':
            if index == 0 or monosaccharides[index-1] != 'N':
                gg_coordinates[1] += 1
                if angled == 'left' and index == 0:
                    gg_coordinates[0] += 1
                if angled == 'right' and index == 0:
                    gg_coordinates[0] -= 1
                x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
                draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, o_glycan = o_glycan)
                monosaccharides_coordinates.append(x_y_coordinates)
                if galnac_start and index == 0:
                    draw_galnac(target_canvas, x_y_coordinates, size)
                else:
                    draw_glcnac(target_canvas, x_y_coordinates, size)
            elif index != 0 and monosaccharides[index-1] == 'N':
                gg_coordinates[1] += 1
                if angled == 'left' and index == 0:
                    gg_coordinates[0] += 1
                if angled == 'right' and index == 0:
                    gg_coordinates[0] -= 1
                x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
                draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled)
                monosaccharides_coordinates.append(x_y_coordinates)
                draw_galnac(target_canvas, x_y_coordinates, size)
        if mono == 'H':
            if monosaccharides.count('H') == len(monosaccharides) and not o_glycan: # High Mannoses drawing
                if index < 2:
                    gg_coordinates[1] += 1
                    if angled == 'left' and index == 0:
                        gg_coordinates[0] += 1
                    if angled == 'right' and index == 0:
                        gg_coordinates[0] -= 1
                    x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
                    draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled)
                    monosaccharides_coordinates.append(x_y_coordinates)
                    draw_mannose(target_canvas, x_y_coordinates, size)
                else:
                    gg_coordinates[1] += 1
                    if angled == 'left' and index == 0:
                        gg_coordinates[0] += 1
                    if angled == 'right' and index == 0:
                        gg_coordinates[0] -= 1
                    x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
                    draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled)
                    monosaccharides_coordinates.append(x_y_coordinates)
                    draw_glucose(target_canvas, x_y_coordinates, size)
            else: # If not High Mannose
                if len(monosaccharides) >= 3 and monosaccharides[index-1] in ['S', 'G', 'E', 'L', 'R', 'A'] and (monosaccharides[index-2] == 'N' or (monosaccharides[index-2] == 'F' and (len(monosaccharides) >= 4 and monosaccharides[index-3] == 'N'))): # Draw angled
                    if (monosaccharides[index-1] in ['S', 'G'] and antenna_number%2 == 0) or monosaccharides[index-1] in ['L', 'A']:
                        angled = 'right'
                        gg_coordinates[0] -= 0.70
                        gg_coordinates[1] += 0.75
                        from_glcnac = True
                    else:
                        angled = 'left'
                        gg_coordinates[0] += 0.70
                        gg_coordinates[1] += 0.75
                        from_glcnac = True
                elif index == 0:
                    if angled == 'left':
                        gg_coordinates[0] += 0.70
                        gg_coordinates[1] += 0.75
                    else:
                        gg_coordinates[0] -= 0.70
                        gg_coordinates[1] += 0.75
                else: # Draw linear
                    gg_coordinates[1] += 1
                x_y_coordinates = calculate_coordinates(*gg_coordinates, size)
                if index == 0:
                    draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, o_glycan = o_glycan)
                else:
                    draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, o_glycan = from_glcnac)
                monosaccharides_coordinates.append(x_y_coordinates)
                draw_galactose(target_canvas, x_y_coordinates, size)
        if mono == 'S':
            special_gg_coordinates = copy.deepcopy(gg_coordinates)
            if index != len(monosaccharides)-1: # Bisect
                from_glcnac = True
                if antenna_number%2 == 0:
                    special_gg_coordinates[0] += 0.70
                    special_gg_coordinates[1] += 0.75
                    angled = 'left'
                else:
                    special_gg_coordinates[0] -= 0.70
                    special_gg_coordinates[1] += 0.75
                    angled = 'right'
            else: # Linear
                from_glcnac = False
                if index == 0 and starting_coordinates == [4, 1]:
                    angled = original_angled
                    if angled == 'left' and index == 0:
                        special_gg_coordinates[0] += 0.70
                        special_gg_coordinates[1] += 0.75
                        from_glcnac = True
                    if angled == 'right' and index == 0:
                        special_gg_coordinates[0] -= 0.70
                        special_gg_coordinates[1] += 0.75
                        from_glcnac = True
                else:
                    special_gg_coordinates[1] += 1
            x_y_coordinates = calculate_coordinates(*special_gg_coordinates, size)
            draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, False, from_glcnac)
            draw_neu5ac(target_canvas, x_y_coordinates, size)
        if mono == 'L':
            special_gg_coordinates = copy.deepcopy(gg_coordinates)
            special_gg_coordinates[0] += 0.70
            special_gg_coordinates[1] += 0.75
            x_y_coordinates = calculate_coordinates(*special_gg_coordinates, size)
            angled = 'left'
            from_glcnac = False
            if monosaccharides[index-1] == 'N':
                from_glcnac = True
            draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, False, o_glycan = from_glcnac)
            draw_neu5ac(target_canvas, x_y_coordinates, size)
        if mono == 'E':
            special_gg_coordinates = copy.deepcopy(gg_coordinates)
            special_gg_coordinates[0] -= 0.70
            special_gg_coordinates[1] += 0.75
            x_y_coordinates = calculate_coordinates(*special_gg_coordinates, size)
            angled = 'right'
            from_glcnac = False
            if monosaccharides[index-1] == 'N':
                from_glcnac = True
            draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, False, o_glycan = from_glcnac)
            draw_neu5ac(target_canvas, x_y_coordinates, size)
        if mono == 'G':
            special_gg_coordinates = copy.deepcopy(gg_coordinates)
            if index != len(monosaccharides)-1: # Bisect
                from_glcnac = True
                if antenna_number%2 == 0:
                    special_gg_coordinates[0] += 0.70
                    special_gg_coordinates[1] += 0.75
                    angled = 'left'
                else:
                    special_gg_coordinates[0] -= 0.70
                    special_gg_coordinates[1] += 0.75
                    angled = 'right'
            else: # Linear
                from_glcnac = False
                if index == 0 and starting_coordinates == [4, 1]:
                    angled = original_angled
                    if angled == 'left' and index == 0:
                        special_gg_coordinates[0] += 0.70
                        special_gg_coordinates[1] += 0.75
                        from_glcnac = True
                    if angled == 'right' and index == 0:
                        special_gg_coordinates[0] -= 0.70
                        special_gg_coordinates[1] += 0.75
                        from_glcnac = True
                else:
                    special_gg_coordinates[1] += 1
            x_y_coordinates = calculate_coordinates(*special_gg_coordinates, size)
            draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, False, from_glcnac)
            draw_neu5gc(target_canvas, x_y_coordinates, size)
        if mono == 'A': #amidated acetyl sialic acid
            special_gg_coordinates = copy.deepcopy(gg_coordinates)
            special_gg_coordinates[0] += 0.70
            special_gg_coordinates[1] += 0.75
            x_y_coordinates = calculate_coordinates(*special_gg_coordinates, size)
            angled = 'left'
            from_glcnac = False
            if monosaccharides[index-1] == 'N':
                from_glcnac = True
            draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, False, o_glycan = from_glcnac)
            draw_neu5gc(target_canvas, x_y_coordinates, size)
        if mono == 'R':
            special_gg_coordinates = copy.deepcopy(gg_coordinates)
            special_gg_coordinates[0] -= 0.70
            special_gg_coordinates[1] += 0.75
            x_y_coordinates = calculate_coordinates(*special_gg_coordinates, size)
            angled = 'right'
            from_glcnac = False
            if monosaccharides[index-1] == 'N':
                from_glcnac = True
            draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, False, o_glycan = from_glcnac)
            draw_neu5gc(target_canvas, x_y_coordinates, size)
        if mono == 'F':
            special_gg_coordinates = copy.deepcopy(gg_coordinates)
            if antenna_number%2 == 0:
                special_gg_coordinates[0] += 0.70
            else:
                special_gg_coordinates[0] -= 0.70
            x_y_coordinates = calculate_coordinates(*special_gg_coordinates, size)
            inclination = 'left' if antenna_number%2 == 0 else 'right'
            draw_connecting_line(target_canvas, monosaccharides_coordinates[-1], x_y_coordinates, size, angled, side = inclination)
            draw_fucose(target_canvas, x_y_coordinates, size)

def draw_generic_antenna(target_canvas, image, size, composition, glycan_class):
    '''Work in Progress.
    '''
    global font
    
    composition = copy.deepcopy(composition)
    
    to_remove = []
    for i in composition:
        if composition[i] <= 0:
            to_remove.append(i)
    for i in to_remove:
        del composition[i]

    if glycan_class == 'n_glycan':
        gg_center_coord = (4, 4)
        center_coordinate = calculate_coordinates(*gg_center_coord, size)
        draw_bracket(target_canvas, center_coordinate, 1, size)

        width = (len(composition)-1)/2
        starting_coords = [4+width, 5]

    if glycan_class == 'o_glycan':
        gg_center_coord = (4, 1)
        center_coordinate = calculate_coordinates(*gg_center_coord, size)
        draw_bracket(target_canvas, center_coordinate, 1, size)

        width = (len(composition)-1)/2
        starting_coords = [4+width, 2]

    ac_sialic_ok = False
    gc_sialic_ok = False
    for index, i in enumerate(composition):
        if index == 0:
            gg_current_coord = starting_coords
        else:
            gg_current_coord[0] -= 1
        current_coords = calculate_coordinates(*gg_current_coord, size)
        if (i == 'S' or i == 'Am' or i == 'E') and composition[i] > 0 and not ac_sialic_ok:
            draw_neu5ac(target_canvas, current_coords, size)
            total_ac = 0
            for j in composition:
                if j == 'S' or j == 'Am' or j == 'E':
                    total_ac += composition[j]
            draw_upside_down_text(image, str(total_ac), current_coords, size, font)
            ac_sialic_ok = True
        if (i == 'G' or i == 'AmG' or i == 'EG') and composition[i] > 0 and not gc_sialic_ok:
            draw_neu5gc(target_canvas, current_coords, size)
            total_gc = 0
            for j in composition:
                if j == 'G' or j == 'AmG' or j == 'EG':
                    total_gc += composition[j]
            draw_upside_down_text(image, str(total_gc), current_coords, size, font)
            gc_sialic_ok = True
        if i == 'H' and composition[i] > 0:
            draw_hexose(target_canvas, current_coords, size)
            draw_upside_down_text(image, str(composition[i]), current_coords, size, font)
        if i == 'N' and composition[i] > 0:
            draw_hexnac(target_canvas, current_coords, size)
            draw_upside_down_text(image, str(composition[i]), current_coords, size, font)
        if i == 'F' and composition[i] > 0:
            draw_fucose(target_canvas, current_coords, size)
            draw_upside_down_text(image, str(composition[i]), current_coords, size, font)

def draw_full_n_glycan(i, index_combo, derivatized_sialic_acids, core_fucose, bissecting_glcnac, figures_folder, glycan_input, size):
    '''
    '''
    global lowest_coordinate, highest_coordinate

    # This avoids high mannoses with more than 3 antennas
    high_mannose = True
    monos_number_antenna = []
    high_mannose_antennas = []
    for index, j in enumerate(i):
        if j.count('H') != len(j):
            high_mannose = False
        else:
            monos_number_antenna.append(len(j))
            high_mannose_antennas.append(index)

    if high_mannose and len(i) > 3:
        return index_combo

    # This avoids more than one antenna in high mannoses having glucoses
    skip = False
    glucose_antenna = False
    glucose_antenna_id = None
    if high_mannose:
        for index_j, j in enumerate(monos_number_antenna):
            if j > 2:
                if not glucose_antenna:
                    glucose_antenna = True
                    glucose_antenna_id = index_j
                else:
                    skip = True
                    break
    if skip:
        return index_combo

    # This changes the glucose-containing antenna in high mannose N-glycans to the leftmost one
    if high_mannose and glucose_antenna:
        new_antenna_combo = list(copy.deepcopy(i))
        glucose_antenna_sequence = new_antenna_combo.pop(glucose_antenna_id)
        new_antenna_combo.insert(1, glucose_antenna_sequence)
        i = new_antenna_combo
    
    # Avoids weird hybrids
    if not high_mannose and len(high_mannose_antennas) > 2:
        return index_combo
    if not high_mannose and len(high_mannose_antennas) == 1 and len(i) >= 3:
        return index_combo
    
    # Put high_mannose antennas together on the right if hybrid
    if not high_mannose and len(high_mannose_antennas) > 0:
        new_antenna_combo = list(copy.deepcopy(i))
        for index, hm_antenna in enumerate(high_mannose_antennas):
            temp_antenna_grab = new_antenna_combo.pop(hm_antenna)
            new_antenna_combo.insert(2*index, temp_antenna_grab)
        i = new_antenna_combo

    # Skip special case of 'HS' or 'S' antenna used specially on o-glycans
    if 'HS' in i or 'S' in i:
        return index_combo

    # Calculate derivatized sialic acids permutations
    derivatized = False
    for g in derivatized_sialic_acids:
        if len(g) > 0:
            derivatized = True
            break
    if derivatized:
        acetyl_permutations = unique_permutations(derivatized_sialic_acids[0])
        glycolyl_permutations = unique_permutations(derivatized_sialic_acids[1])

        sia_possibilities = []
        for perm_a in acetyl_permutations:
            for perm_g in glycolyl_permutations:
                sia_possibilities.append(replace_s_and_g_with_chars(i, perm_a, perm_g))

        # Goes through each permutation of derivatized sialic acids and fixes antennas
        for new_i in sia_possibilities:

            # Create a new image with transparent background
            image = Image.new("RGBA", (image_width, image_height), background_color)

            # Create an object to draw on the image
            draw = ImageDraw.Draw(image)

            # Checks if any one antenna has bifurcation to widen the N-glycan core mannoses and correct starting coordinates for antennas
            wide = [False, False]
            starting_coords = [[5,4], [3,4]]
            if len(new_i) > 2:
                wide[1] = True
                starting_coords[1][0] -= 1
                if len(new_i) == 4:
                    wide[0] = True
                    starting_coords[0][0] += 1

            # Draw the N-glycan core
            draw_n_glycan_core(draw, size, wide, starting_coordinates = [4,1])

            # Draw the N-glycan antennas
            for index, j in enumerate(new_i):
                angled_antenna = False
                if len(new_i) > 3:
                    if index < 2:
                        angled_antenna = 'left'
                    else:
                        angled_antenna = 'right'
                elif len(new_i) == 3:
                    if index == 0:
                        angled_antenna = 'left'
                    if index == 2:
                        angled_antenna = 'right'
                draw_antenna(draw, size, j, starting_coordinates = copy.deepcopy(starting_coords[(index+1)%2]), angled = angled_antenna, antenna_number = index+1)

            if core_fucose:
                # Draw the core fucose
                x_y_coordinates = calculate_coordinates(3, 1, size)
                inclination = 'right'
                draw_connecting_line(draw, calculate_coordinates(4, 1, size), x_y_coordinates, size, side = inclination)
                draw_fucose(draw, x_y_coordinates, size)

            if bissecting_glcnac:
                # Draw the bissecting GlcNAc
                x_y_coordinates = calculate_coordinates(4, 4, size)
                draw_connecting_line(draw, calculate_coordinates(4, 3, size), x_y_coordinates, size)
                draw_glcnac(draw, x_y_coordinates, size)

            # Define the crop box (left, upper, right, lower)
            crop_box = (lowest_coordinate[0], lowest_coordinate[1], highest_coordinate[0], highest_coordinate[1])

            # Crop the image
            image = image.crop(crop_box)

            # Rotate the image by 180 degrees
            image = image.rotate(180, expand=True)

            # Save the image to a file
            core_f = ''
            if core_fucose:
                core_f = '_F'
            bissecting = ''
            if bissecting_glcnac:
                bissecting = '_B'
            
            os.makedirs(os.path.join(figures_folder, str(glycan_input)), exist_ok=True)
            image.save(os.path.join(os.path.join(figures_folder, f"{glycan_input}"), str(index_combo+1)+f"_{len(i)}"+f"_N"+core_f+bissecting+'.png'))

            # Reset the cropping coordinates
            highest_coordinate = [0, 0]
            lowest_coordinate = [99999, 99999]
            index_combo+= 1

    else:
        # Create a new image with transparent background
        image = Image.new("RGBA", (image_width, image_height), background_color)

        # Create an object to draw on the image
        draw = ImageDraw.Draw(image)

        # Checks if any one antenna has bifurcation to widen the N-glycan core mannoses and correct starting coordinates for antennas
        wide = [False, False]
        starting_coords = [[5,4], [3,4]]
        if len(i) > 2:
            wide[1] = True
            starting_coords[1][0] -= 1
            if len(i) == 4:
                wide[0] = True
                starting_coords[0][0] += 1

        # Draw the N-glycan core
        draw_n_glycan_core(draw, size, wide, starting_coordinates = [4,1])

        # Draw the N-glycan antennas
        for index, j in enumerate(i):
            angled_antenna = False
            if len(i) > 3:
                if index < 2:
                    angled_antenna = 'left'
                else:
                    angled_antenna = 'right'
            elif len(i) == 3:
                if index == 0:
                    angled_antenna = 'left'
                if index == 2:
                    angled_antenna = 'right'
            draw_antenna(draw, size, j, starting_coordinates = copy.deepcopy(starting_coords[(index+1)%2]), angled = angled_antenna, antenna_number = index+1)

        if core_fucose:
            # Draw the core fucose
            x_y_coordinates = calculate_coordinates(3, 1, size)
            inclination = 'right'
            draw_connecting_line(draw, calculate_coordinates(4, 1, size), x_y_coordinates, size, side = inclination)
            draw_fucose(draw, x_y_coordinates, size)

        if bissecting_glcnac:
            # Draw the bissecting GlcNAc
            x_y_coordinates = calculate_coordinates(4, 4, size)
            draw_connecting_line(draw, calculate_coordinates(4, 3, size), x_y_coordinates, size)
            draw_glcnac(draw, x_y_coordinates, size)

        # Define the crop box (left, upper, right, lower)
        crop_box = (lowest_coordinate[0], lowest_coordinate[1], highest_coordinate[0], highest_coordinate[1])

        # Crop the image
        image = image.crop(crop_box)

        # Rotate the image by 180 degrees
        image = image.rotate(180, expand=True)

        # Save the image to a file
        core_f = ''
        if core_fucose:
            core_f = '_F'
        bissecting = ''
        if bissecting_glcnac:
            bissecting = '_B'
            
        os.makedirs(os.path.join(figures_folder, str(glycan_input)), exist_ok=True)
        image.save(os.path.join(os.path.join(figures_folder, f"{glycan_input}"), str(index_combo+1)+f"_{len(i)}"+f"_N"+core_f+bissecting+'.png'))

        # Reset the cropping coordinates
        highest_coordinate = [0, 0]
        lowest_coordinate = [99999, 99999]
        index_combo+= 1
    return index_combo

def draw_full_o_glycan(i, index_combo, derivatized_sialic_acids, figures_folder, glycan_input, size):
    '''
    '''
    global lowest_coordinate, highest_coordinate

    # This loop allows to make options with antennas starting with GalNAc and GlcNAc
    hexnac_antennas_positions = [False, False]
    hexnac_sorting_range = 0
    for index, l in enumerate(i):
        if l[0] == 'N':
            hexnac_antennas_positions[index] = True
            hexnac_sorting_range += 2
    if hexnac_sorting_range == 0:
        hexnac_sorting_range = 1


    # Calculate derivatized sialic acids permutations
    derivatized = False
    for g in derivatized_sialic_acids:
        if len(g) > 0:
            derivatized = True
            break
    if derivatized:
        acetyl_permutations = unique_permutations(derivatized_sialic_acids[0])
        glycolyl_permutations = unique_permutations(derivatized_sialic_acids[1])

        sia_possibilities = []
        for perm_a in acetyl_permutations:
            for perm_g in glycolyl_permutations:
                sia_possibilities.append(replace_s_and_g_with_chars(i, perm_a, perm_g))

        # Goes through each permutation of derivatized sialic acids and fixes antennas
        for new_i in sia_possibilities:

            for l in range(hexnac_sorting_range):

                # Create a new image with transparent background
                image = Image.new("RGBA", (image_width, image_height), background_color)

                # Create an object to draw on the image
                draw = ImageDraw.Draw(image)

                # Change starting point for antennas depending on the core
                starting_coords = [[4,1], [4,1]]

                # Draw the O-glycan core
                draw_o_glycan_core(draw, size, starting_coordinates = [4,1])

                # Draw the O-glycan antennas
                angled_antenna = False
                for index, j in enumerate(new_i):
                    if index%2 == 0:
                        angled_antenna = 'left'
                    elif index%2 != 0:
                        angled_antenna = 'right'

                    # The possibilities of antennas starting with GalNAc or GlcNAc
                    if l == 0:
                        antennas_galnac_start = [False, False]
                    elif l == 1 and hexnac_sorting_range == 2:
                        antennas_galnac_start = hexnac_antennas_positions
                    elif l == 1 and hexnac_sorting_range > 2:
                        antennas_galnac_start = [True, False]
                    elif l == 2:
                        antennas_galnac_start = [False, True]
                    elif l == 3:
                        antennas_galnac_start = [True, True]

                    draw_antenna(draw, size, j, starting_coordinates = copy.deepcopy(starting_coords[index%2]), angled = angled_antenna, antenna_number = index, o_glycan = True, galnac_start = antennas_galnac_start[index%2])

                # Define the crop box (left, upper, right, lower)
                crop_box = (lowest_coordinate[0], lowest_coordinate[1], highest_coordinate[0], highest_coordinate[1])

                # Crop the image
                image = image.crop(crop_box)

                # Rotate the image by 180 degrees
                image = image.rotate(180, expand=True)

                # Save the image to a file (Optional)
            
                os.makedirs(os.path.join(figures_folder, str(glycan_input)), exist_ok=True)
                image.save(os.path.join(os.path.join(figures_folder, f"{glycan_input}"), str(index_combo+1)+f"_{len(i)}"+f"_O"+'.png'))

                # Reset the cropping coordinates
                highest_coordinate = [0, 0]
                lowest_coordinate = [99999, 99999]
                index_combo += 1
    else:
        for l in range(hexnac_sorting_range):

            # Create a new image with transparent background
            image = Image.new("RGBA", (image_width, image_height), background_color)

            # Create an object to draw on the image
            draw = ImageDraw.Draw(image)

            # Change starting point for antennas depending on the core
            starting_coords = [[4,1], [4,1]]

            # Draw the O-glycan core
            draw_o_glycan_core(draw, size, starting_coordinates = [4,1])

            # Draw the O-glycan antennas
            angled_antenna = False
            for index, j in enumerate(i):
                if index%2 == 0:
                    angled_antenna = 'left'
                elif index%2 != 0:
                    angled_antenna = 'right'

                # The possibilities of antennas starting with GalNAc or GlcNAc
                if l == 0:
                    antennas_galnac_start = [False, False]
                elif l == 1 and hexnac_sorting_range == 2:
                    antennas_galnac_start = hexnac_antennas_positions
                elif l == 1 and hexnac_sorting_range > 2:
                    antennas_galnac_start = [True, False]
                elif l == 2:
                    antennas_galnac_start = [False, True]
                elif l == 3:
                    antennas_galnac_start = [True, True]

                draw_antenna(draw, size, j, starting_coordinates = copy.deepcopy(starting_coords[index%2]), angled = angled_antenna, antenna_number = index, o_glycan = True, galnac_start = antennas_galnac_start[index%2])

            # Define the crop box (left, upper, right, lower)
            crop_box = (lowest_coordinate[0], lowest_coordinate[1], highest_coordinate[0], highest_coordinate[1])

            # Crop the image
            image = image.crop(crop_box)

            # Rotate the image by 180 degrees
            image = image.rotate(180, expand=True)

            # Save the image to a file (Optional)
            os.makedirs(os.path.join(figures_folder, str(glycan_input)), exist_ok=True)
            image.save(os.path.join(os.path.join(figures_folder, f"{glycan_input}"), str(index_combo+1)+f"_{len(i)}"+f"_O"+'.png'))

            # Reset the cropping coordinates
            highest_coordinate = [0, 0]
            lowest_coordinate = [99999, 99999]
            index_combo += 1

    return index_combo
    
def draw_glycan(figures_folder, size, glycan_input, glycan_class):
    '''
    '''
    global font, lowest_coordinate, highest_coordinate

    current_dir = pathlib.Path(__file__).parent.resolve()
    
    # Load the font (use default font if custom font is not available)
    try:
        font = ImageFont.truetype("Arial.ttf", round(size*0.8))  # Specify font and size
    except IOError:
        try:
            font = ImageFont.truetype(list_available_fonts()[0], round(size*0.8))  # Specify font and size
        except IOError:
            font = ImageFont.load_default()
            
    # Change 'Am', 'EG' and 'AmG' to 'S' and count their amounts for later usage
    derivatization_test = form_to_comp(glycan_input)
    am_count = derivatization_test.get('Am')
    amg_count = derivatization_test.get('AmG')
    e_count = derivatization_test.get('E')
    eg_count = derivatization_test.get('EG')
    derivatized_acetyl_acids = ''
    derivatized_glycolyl_acids = ''
    for i in range(am_count if am_count != None else 0):
        derivatized_acetyl_acids += 'L'
    for i in range(e_count if e_count != None else 0):
        derivatized_acetyl_acids += 'E'
    for i in range(amg_count if amg_count != None else 0):
        derivatized_glycolyl_acids += 'A'
    for i in range(eg_count if eg_count != None else 0):
        derivatized_glycolyl_acids += 'R'
    derivatized_sialic_acids = [derivatized_acetyl_acids, derivatized_glycolyl_acids]

    if am_count != None or e_count != None:
        derivatization_test['S'] = 0
        if am_count != None:
            derivatization_test['S'] += derivatization_test.pop('Am')
        if e_count != None:
            derivatization_test['S'] += derivatization_test.pop('E')

    if amg_count != None or eg_count != None:
        derivatization_test['G'] = 0
        if amg_count != None:
            derivatization_test['G'] += derivatization_test.pop('AmG')
        if eg_count != None:
            derivatization_test['G'] += derivatization_test.pop('EG')

    glycan = comp_to_formula(derivatization_test)

    index_combo = 0

    if glycan_class == 'n_glycans' or glycan_class == 'none':

        # Loads antennas from files
        with open(os.path.join(current_dir, "Assets/antennas_n.gga"), "rb") as f:
            combinations_dict = dill.load(f)
            f.close()

        # Turn the glycan into composition
        glycan_comp = form_to_comp(glycan)

        # Check if it could be a N-glycan
        if ('N' in glycan_comp and glycan_comp['N'] >=2) and ('H' in glycan_comp and glycan_comp['H'] >= 3):
            
            # Remove core from glycan
            glycan_comp['N'] -= 2
            glycan_comp['H'] -= 3

            if glycan_comp['N'] >= 0 and glycan_comp['H'] >= 0:
                
                # Draw generic formula
                image = Image.new("RGBA", (image_width, image_height), background_color)
                draw = ImageDraw.Draw(image)
                draw_n_glycan_core(draw, size, [False, False], starting_coordinates = [4,1])
                draw_generic_antenna(draw, image, size, glycan_comp, 'n_glycan')
                crop_box = (lowest_coordinate[0], lowest_coordinate[1], highest_coordinate[0], highest_coordinate[1])
                image = image.crop(crop_box)
                image = image.rotate(180, expand=True)
            
                os.makedirs(os.path.join(figures_folder, str(glycan_input)), exist_ok=True)
                image.save(os.path.join(os.path.join(figures_folder, f"{glycan_input}"), f'{index_combo+1}_G_N.png'))
                
                index_combo += 1
                
                # Reset the cropping coordinates
                highest_coordinate = [0, 0]
                lowest_coordinate = [99999, 99999]
                
                # Find matching antenna composition
                for i in permute_comp_keys(glycan_comp):
                    possible_antennas_glycan = combinations_dict.get(comp_to_formula(i))
                    if possible_antennas_glycan != None:
                        break

                # If antenna composition found, start drawing the possibilities
                if possible_antennas_glycan != None:

                    for i in possible_antennas_glycan:
                        index_combo = draw_full_n_glycan(i, index_combo, derivatized_sialic_acids, False, False, figures_folder, glycan_input, size)

                # Core fucosylation check
                if 'F' in glycan_comp.keys():

                    current_comp = copy.deepcopy(glycan_comp)

                    current_comp['F'] -= 1

                    # Find matching antenna composition
                    for i in permute_comp_keys(current_comp):
                        possible_antennas_glycan = combinations_dict.get(comp_to_formula(i))
                        if possible_antennas_glycan != None:
                            break

                    # If antenna composition found, start drawing the possibilities
                    if possible_antennas_glycan != None:
                        for i in possible_antennas_glycan:
                            index_combo = draw_full_n_glycan(i, index_combo, derivatized_sialic_acids, True, False, figures_folder, glycan_input, size)
                
                # Bissecting GlcNAc check
                if glycan_comp['N'] >= 1:

                    current_comp = copy.deepcopy(glycan_comp)

                    current_comp['N'] -= 1

                    # Find matching antenna composition
                    for i in permute_comp_keys(current_comp):
                        possible_antennas_glycan = combinations_dict.get(comp_to_formula(i))
                        if possible_antennas_glycan != None:
                            break

                    # If antenna composition found, start drawing the possibilities
                    if possible_antennas_glycan != None:
                        for i in possible_antennas_glycan:
                            index_combo = draw_full_n_glycan(i, index_combo, derivatized_sialic_acids, False, True, figures_folder, glycan_input, size)

                    # This checks for both core fucosylation and bissecting GlcNAc
                    if 'F' in current_comp:

                        current_comp['F'] -= 1

                        # Find matching antenna composition
                        for i in permute_comp_keys(current_comp):
                            possible_antennas_glycan = combinations_dict.get(comp_to_formula(i))
                            if possible_antennas_glycan != None:
                                break

                        # If antenna composition found, start drawing the possibilities
                        if possible_antennas_glycan != None:
                            for i in possible_antennas_glycan:
                                index_combo = draw_full_n_glycan(i, index_combo, derivatized_sialic_acids, True, True, figures_folder, glycan_input, size)

    if glycan_class == 'o_glycans' or glycan_class == 'none':

        # Loads antennas from files
        with open(os.path.join(current_dir, "Assets/antennas_o.gga"), "rb") as f:
            combinations_dict = dill.load(f)
            f.close()

        # Turn the glycan into composition
        glycan_comp = form_to_comp(glycan)

        # Check if 'N' is in the composition
        if 'N' in glycan_comp and glycan_comp['N'] > 0:

            # Remove core from glycan
            glycan_comp['N'] -= 1

            # Draw generic formula
            image = Image.new("RGBA", (image_width, image_height), background_color)
            draw = ImageDraw.Draw(image)
            draw_o_glycan_core(draw, size, starting_coordinates = [4,1])
            draw_generic_antenna(draw, image, size, glycan_comp, 'o_glycan')
            crop_box = (lowest_coordinate[0], lowest_coordinate[1], highest_coordinate[0], highest_coordinate[1])
            image = image.crop(crop_box)
            image = image.rotate(180, expand=True)
            
            os.makedirs(os.path.join(figures_folder, str(glycan_input)), exist_ok=True)
            image.save(os.path.join(os.path.join(figures_folder, f"{glycan_input}"), f'{index_combo+1}_G_O.png'))
            index_combo += 1

            # Reset the cropping coordinates
            highest_coordinate = [0, 0]
            lowest_coordinate = [99999, 99999]

            # Find matching antenna composition
            for i in permute_comp_keys(glycan_comp):
                possible_antennas_glycan = combinations_dict.get(comp_to_formula(i))
                if possible_antennas_glycan != None:
                    break

            # If antenna composition found, start drawing the possibilities
            if possible_antennas_glycan != None:
                for i in possible_antennas_glycan:
                    index_combo = draw_full_o_glycan(i, index_combo, derivatized_sialic_acids, figures_folder, glycan_input, size)
                    
if __name__ == '__main__':
    draw_glycan("D:/Arquivos/Desktop/test/gg_draw_tests", 30, "H1N1E1", 'none')