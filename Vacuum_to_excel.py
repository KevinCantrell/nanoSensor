#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 12:48:27 2021

@author: zionirving-singh
"""

import os
import xlsxwriter


#imports the Tkinter module in a way that works for python 2 or 3
try: #works on Python 2.7
    import Tkinter as tk
    from tkFileDialog import askopenfilename
except: #works on Python 3.7
    import tkinter as tk 
    from tkinter.filedialog import askopenfilename
    
#gets rid of tkinter blank window
root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)
 
#the following line lets the user select a file that can be processed later
file_path = askopenfilename()
hooper = 8
#creates an alias "startpath" for whenever you want os to get a directory name
startpath = os.path.dirname(file_path)

#creates a new excel workbook named pressure_summary
workbook = xlsxwriter.Workbook(os.path.join(startpath, 'pressure_summary.xlsx'))

#creates the two subsheets in the workbook
rawData = workbook.add_worksheet('Raw Data')
press = workbook.add_worksheet('Pressure Graph')

#adds a blank chart into the pressure subsheet
plot = workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})

#converts the scientific notation of the pressure data into regular numbers
num_format = workbook.add_format({'num_format': '@'})

#opens the file you select and then picks out only the pressure data
filename = file_path                                 #
vac_data = []                                         # Declare an empty list named mylines.
with open (filename , 'rt') as myfile:               # Open pump data text for reading.
    for myline in myfile:                            # For each line in the file 
        vac_data.append(myline[4:12])    # isolates the pressure readout

row = 0
col = 0
rowSum = 0
colSum = 0
filetoprocess = [file_path]
#writes the pressure values to the raw data worksheet in number format
for pressure in vac_data:
    rawData.write(row, col, pressure, num_format)
    row += 1
lengthvacdata = len(vac_data)
time = list(range(0, lengthvacdata, int(0.5))
            

.add_series({
    'categories': ['Time', time],
    'values': ['Pressure', vac_data]})
 
.set_y_axis({
    'crossing': '-2',
    'major_gridlines': {
    'visible': False},
    'name': 'Pressure (torr)',
    'min': 0.2, 
    'max': 5,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})
.set_x_axis({
    'crossing': '-2',
    'major_gridlines': {
    'visible': False},
    'name': 'Time (mins)',
    'min': 0, 
    'max': 120,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})

plot.insert_chart('B1', press)
