#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:25:54 2019

@author: hallw and cantrellk

edited by irvingsinghz
"""
# this fits spectral data to two Gaussians
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import xlsxwriter

import os
import statistics as st
import scipy.signal
from scipy.optimize import curve_fit


# defines a single Gaussian function
def gaussian(x, mu, sigma,
             amp):  # this defines a function with the paramaters x=wavelength, mu=peak center, sigma=peak width, amp=peak height
    # Note that x data MUST be the first argument when defining a function that will be used with CurveFit
    y = (amp / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp((-(x - mu) ** 2) / (2 * sigma ** 2))
    return y


import tkinter as tk
from tkinter.filedialog import askopenfilename

# gets rid of tkinter blank window
root = tk.Tk()
root.withdraw()
root.wm_attributes('-topmost', 1)

# the following line lets the user select a file that can be processed later
file_path = askopenfilename()

# creates an alias "startpath" for whenever you want os to get a directory name
startpath = os.path.dirname(file_path)

# the following block of code lets the user choose whether to process ONLY the file selected during the askopenfilename (Option 1)
# or whether to also process all other files in that directory (option 2)
# or whether to also process all other files in that directory AND all subdirectories (option 3)
mode = input('root is' + os.path.dirname(
    file_path) + '\n' + '\n' + 'Do you want to process one file (1), all files in this directory (2), or all files in this directory and all subdirectories (3)?')  # asks user which files to process
filestoprocess = []  # initializes "files to process" as an empty list. We can then add files to process to this list based on option selected by user
if mode == "1":
    filestoprocess.append(file_path)  # only file processed will be the file selected during askopenfilename
if mode == "2":
    for filename in os.listdir(startpath):  # gets all the filenames in the directory
        if filename[-3:] == "csv":  # if last 3 characters in filename are csv, then do the next line
            filestoprocess.append(os.path.join(startpath,
                                               filename))  # joins the directory to the filename to create a filepath, then adds that filepath to the "files to process" list
if mode == "3":
    for (dirpath, dirnames, filenames) in os.walk(
            startpath):  # os.walk walks through the directory AND all subdirectories where a given file is stored; dirpath dirnames filenames lists all files in the directory and subdirectories
        for filename in filenames:
            if filename[-3:] == "csv":  # if last 3 characters in filename are csv, then do the next line
                filestoprocess.append(os.path.join(os.path.normpath(dirpath), filename))

# creating workbook
workbook = xlsxwriter.Workbook(os.path.join(startpath, 'Data_Summary.xlsx'))

# adding worksheets with names
rawData = workbook.add_worksheet('Raw Data')
smoothedDataSheet = workbook.add_worksheet('Smoothed Data')
GaussDataSheet = workbook.add_worksheet('Gaussian Data')
summaryStats = workbook.add_worksheet('Summary')
averages = workbook.add_worksheet('Avg and Stdev')
# refractiveIndexSheet = workbook.add_worksheet('Refractive Index')
refractPlotSheet = workbook.add_worksheet('Refractive Index Plot')
rawDataPlotSheet = workbook.add_worksheet('Raw Data Plot')
smoothedDataPlotSheet = workbook.add_worksheet('Smoothed Data Plot')
GaussDataPlotSheet = workbook.add_worksheet('Gaussian Data Plot')

# creating charts
rawDataPlot = workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
smoothedDataPlot = workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
GaussDataPlot = workbook.add_chart({'type': 'scatter', 'subtype': 'smooth'})
refractPlot = workbook.add_chart({'type': 'scatter'})
row = 1
col = 0
rowSum = 0
colSum = 0

# creates a figure to plot
# fig,ax=plt.subplots()

# creates lists to store the values from each csv file, to be used at the end
gauss_lam_max = []
gauss_intensity = []
raw_lam_max = []
raw_intensity = []
smooth_lam_max = []
smooth_intensity = []
# for each file in filesToProcess performs data analysis and writes to excel file
for filename in filestoprocess:
    baseName = os.path.basename(os.path.basename(filename))
    # reads data from array containing list of files
    dfHead = pd.read_csv(filename, skiprows=1, skipfooter=513, delimiter=',', index_col=0, header=None, engine='python')
    dfType = pd.read_csv(filename, skiprows=64, skipfooter=523, delimiter=',', header=None, engine='python')
    if dfType.iloc[0, 0] == 'sel_pixel_start':
        skiprows = 89
    else:
        skiprows = 80
    dfSpectrum = pd.read_csv(filename, skiprows=skiprows, delimiter=',', index_col=0, header=None,
                             names=['Waves', 'Wavenumber', 'Raman_Shift', 'Dark', 'Reference', 'Raw', 'Dark_Subtracted',
                                    'Transmission', 'Absorbance', 'Irradiance', 'blk'], engine='python')
    dfSpectrum['Wavelength'] = float(dfHead.loc['coefs_a0'].values[0]) + float(
        dfHead.loc['coefs_a1'].values[0]) * dfSpectrum.index + float(
        dfHead.loc['coefs_a2'].values[0]) * dfSpectrum.index ** 2 + float(
        dfHead.loc['coefs_a3'].values[0]) * dfSpectrum.index ** 3

    dfSpectrum['Absorbance'] = dfSpectrum.iloc[:, 8]

    # truncates data to region of interest (the plasmonic peak region) based on wavelength
    rangeBool = (dfSpectrum['Wavelength'] >= 500) & (dfSpectrum['Wavelength'] <= 600)
    # This takes the x data from dfSpectrum and makes a Boolean array with values "true" between 500 and 600 and "false" everywhere else

    # returns full row of data for row where absorbance is a maximum within the wavelength range defined in rangeBool
    maxVals = dfSpectrum.loc[dfSpectrum['Absorbance'][rangeBool].idxmax()]

    # truncates absorbance array to only the upper 25% of absorbance values
    rangeBoolIntens = (dfSpectrum['Absorbance'] >= 0.75 * maxVals['Absorbance']) & (dfSpectrum['Wavelength'] >= 500) & (
                dfSpectrum['Wavelength'] <= 600)
    # y = dfSpectrum['Absorbance'][rangeBoolIntens]
    # x = np.array(dfSpectrum['Wavelength'][rangeBoolIntens])

    # calculates Gaussian Fit to upper 25% of spectrum popt does parameter optimization for the gaussian function and
    # specifies the fit will be for the x and y data from rangeBoolIntens with initial parameter guesses for mu,
    # sigma, amp of 520, 30, and 87
    popt, pcov = curve_fit(gaussian, dfSpectrum['Wavelength'][rangeBoolIntens],
                           dfSpectrum['Absorbance'][rangeBoolIntens], p0=[520, 30, 87])
    # this creates a gaussian curve (y data) using the optimized parameters determined using popt above
    GaussData = gaussian(dfSpectrum['Wavelength'][rangeBoolIntens], popt[0], popt[1], popt[2])

    # smoothing data, 2nd argument is the pixels to smooth, 3rd argument is polynomial order)
    smoothedData = scipy.signal.savgol_filter(dfSpectrum['Absorbance'], 51, 5)

    # calculating smoothed max
    y = smoothedData[rangeBool]  # get smoothed absorbance spectrum values
    x = np.array(dfSpectrum['Wavelength'][rangeBool].values)  # gets wavelengths for absorbance values
    smoothedMaxAbsIndex = np.argmax(y)  # finds index where absorbance is a maximum
    smoothedDataLambda = x[smoothedMaxAbsIndex]  # gives wavelength at that index

    # Creates labels for columns and writes Wavelength values once in "raw data", "smoothed data" and "Gaussian data"
    # worksheet in Excel
    if col == 0:
        rawData.write(row - 1, col, 'Wavelength (nm)')
        smoothedDataSheet.write(row - 1, col, 'Wavelength (nm)')
        GaussDataSheet.write(row - 1, col, 'Wavelength (nm)')
        rawData.write_column(row, col, dfSpectrum['Wavelength'])
        smoothedDataSheet.write_column(row, col, dfSpectrum['Wavelength'])
        GaussDataSheet.write_column(row, col, dfSpectrum['Wavelength'][rangeBoolIntens])
        col += 1

    # writes Absorbance values for each file along with filename in data worksheets
    rawData.write(row - 1, col, os.path.basename(filename))
    rawData.write_column(row, col, dfSpectrum['Absorbance'])
    smoothedDataSheet.write(row - 1, col, os.path.basename(filename))
    smoothedDataSheet.write_column(row, col, smoothedData)
    GaussDataSheet.write(row - 1, col, os.path.basename(filename))
    GaussDataSheet.write_column(row, col, GaussData)
    col += 1

    # calculates length of wavelength column for raw and smoothed data
    lengthWaveCol = len(dfSpectrum['Wavelength'])

    # adding data to be plotted to rawData
    rawDataPlot.add_series({
        'name': ['Raw Data', row - 1, col - 1],
        'categories': ['Raw Data', row, 0, lengthWaveCol, 0],
        'values': ['Raw Data', row, col - 1, lengthWaveCol, col - 1]})

    # adding data to be plotted to smoothedData
    smoothedDataPlot.add_series({
        'name': ['Smoothed Data', row - 1, col - 1],
        'categories': ['Smoothed Data', row, 0, lengthWaveCol, 0],
        'values': ['Smoothed Data', row, col - 1, lengthWaveCol, col - 1]})

    # calculates length of wavelength column for Gaussian Data
    lengthWaveColGauss = len(dfSpectrum['Wavelength'][rangeBoolIntens])

    # adding data to be plotted to GaussData
    GaussDataPlot.add_series({
        'name': ['Gaussian Data', row - 1, col - 1],
        'categories': ['Gaussian Data', row, 0, lengthWaveColGauss, 0],
        'values': ['Gaussian Data', row, col - 1, lengthWaveColGauss, col - 1]})

    # adjusting column width for rawData and smoothedData and Gaussian Data
    rawData.set_column(0, col, 30)
    smoothedDataSheet.set_column(0, col, 30)
    GaussDataSheet.set_column(0, col, 30)

    # labels columns in data sheets
    if rowSum == 0:
        summaryStats.write(rowSum, colSum, 'Filename')
        summaryStats.write(rowSum, colSum + 1, 'Raw Lambda Max (nm)')
        summaryStats.write(rowSum, colSum + 2, 'Raw Absorbance')
        summaryStats.write(rowSum, colSum + 3, 'Smoothed Lambda Max (nm)')
        summaryStats.write(rowSum, colSum + 4, 'Smoothed Absorbance')
        summaryStats.write(rowSum, colSum + 5, 'Gaussian Lambda Max (nm)')
        summaryStats.write(rowSum, colSum + 6, 'Gaussian Intensity (nm)')
        averages.write(rowSum, colSum, 'Type')
        averages.write(rowSum, colSum + 1, 'Avg Gauss Lambda max (nm)')
        averages.write(rowSum, colSum + 2, 'Avg Gauss Intensity (nm)')
        averages.write(rowSum, colSum + 3, 'Stdev of Gauss Lambda Max')
        averages.write(rowSum, colSum + 4, 'Stdev of Gauss Intensity')
        averages.write(rowSum, colSum + 5, 'Variance of Gauss Lambda Max')
        averages.write(rowSum, colSum + 6, 'Variance of Gauss Intensity')

        averages.write(rowSum, colSum + 7, 'Avg Raw Lambda max (nm)')
        averages.write(rowSum, colSum + 8, 'Avg raw Intensity (nm)')
        averages.write(rowSum, colSum + 9, 'Stdev of raw Lambda Max')
        averages.write(rowSum, colSum + 10, 'Stdev of raw Intensity')
        averages.write(rowSum, colSum + 11, 'Variance of raw Lambda Max')
        averages.write(rowSum, colSum + 12, 'Variance of raw Intensity')

        averages.write(rowSum, colSum + 13, 'Avg smooth Lambda max (nm)')
        averages.write(rowSum, colSum + 14, 'Avg smooth Intensity (nm)')
        averages.write(rowSum, colSum + 15, 'Stdev of smooth Lambda Max')
        averages.write(rowSum, colSum + 16, 'Stdev of smooth Intensity')
        averages.write(rowSum, colSum + 17, 'Variance of smooth Lambda Max')
        averages.write(rowSum, colSum + 18, 'Variance of smooth Intensity')
        #    refractiveIndexSheet.write(rowSum, colSum, 'Medium')
        #   refractiveIndexSheet.write(rowSum, colSum + 1, 'Refractive Index')
        #    refractiveIndexSheet.write(rowSum, colSum + 2, 'Smoothed Lambda Max (nm)')
        #     refractiveIndexSheet.write(rowSum, colSum + 3, 'Gaussian Lambda Max (nm)')
        rowSum += 1

    # little blurb of code that outputs these values form each individual csv to the overall code
    if popt[0] > 0:
        gauss_lam_max.append(popt[0])
        gauss_intensity.append(popt[2])
        raw_lam_max.append(maxVals['Wavelength'])
        raw_intensity.append(maxVals['Absorbance'])
        smooth_lam_max.append(smoothedDataLambda)
        smooth_intensity.append(smoothedData[smoothedMaxAbsIndex])

    # avg_lam_max = st.fmean(gauss_lam_max)
    # avg_intensity = st.fmeanlist(gauss_intensity)
    # lam_dev = st.stdev(gauss_lam_max)
    # intensity_dev = st.stdev(gauss_intensity)

    # dfResults=pd.DataFrame( data = (popt[0], popt[2]),  index= 'Filename' , columns=['lambda max', 'intensity'])
    # dfResults.iloc[index,0]= basename
    # dfResults.iloc[index,1]=popt[0]
    # dfResults.iloc[index,2]=popt[2]
    # dfResults.iloc[index,3]=mean(dfResults['lambda max'])
    # dfResults.iloc[index,4]=mean(dfResults['intensity'])
    # dfResults.iloc[index,5]=statistics.stdev(dfResults['lambda max'])
    # dfResults.iloc[index,6]=statistics.stdev(dfResults['intensity'])

    # determining refractive index from the 7th - 10th characters in filename. Makes lowercase so that name is case-insensitive
    #   refractiveIndexNum = 0
    #  sensorMedium = baseName[7:10].lower()
    #  if sensorMedium  == 'air':
    #      refractiveIndexNum = 1.0000
    #  if sensorMedium  == 'h20':
    #      refractiveIndexNum = 1.3333
    #   if sensorMedium  == '10%':
    #       refractiveIndexNum = 1.3478
    #   if sensorMedium  == '20%':
    #        refractiveIndexNum = 1.3639
    #   if sensorMedium  == '30%':
    #       refractiveIndexNum = 1.3812
    #    if sensorMedium  == '40%':
    #       refractiveIndexNum = 1.3999
    #   if sensorMedium  == '50%':
    #     refractiveIndexNum = 1.4201

    # writes neccessary data into correct column
    summaryStats.write(rowSum, colSum, baseName)
    summaryStats.write(rowSum, colSum + 1, maxVals['Wavelength'])
    summaryStats.write(rowSum, colSum + 2, maxVals['Absorbance'])
    summaryStats.write(rowSum, colSum + 3, smoothedDataLambda)
    summaryStats.write(rowSum, colSum + 4, smoothedData[smoothedMaxAbsIndex])
    summaryStats.write(rowSum, colSum + 5, popt[0])
    summaryStats.write(rowSum, colSum + 6, popt[2])
    # summaryStats.write(rowSum, colSum + 7, avg_lam_max)
    # summaryStats.write(rowSum, colSum + 8, avg_intensity)
    # summaryStats.write(rowSum, colSum + 9, lam_dev)
    # summaryStats.write(rowSum, colSum + 10, intensity_dev)
    #  refractiveIndexSheet.write(rowSum, colSum + 0, sensorMedium)
    #   refractiveIndexSheet.write(rowSum, colSum +1, refractiveIndexNum)
    #   refractiveIndexSheet.write(rowSum, colSum + 2, smoothedDataLambda)
    #   refractiveIndexSheet.write(rowSum, colSum + 3, popt[0])
    rowSum += 1
    # adjusting column width for summaryStats, hardcoded RIP
    summaryStats.set_column(0, 6, 30)
    # refractiveIndexSheet.set_column(0, 6, 25)

    # ax.plot(dfSpectrum['Wavelength'], smoothedData, label = 'Smoothed Data')
    # ax.plot(dfSpectrum['Wavelength'], dfSpectrum['Absorbance'], label = 'Raw Data')
    # plots upper 25% of spectrum only
    # fig,ax=plt.subplots()
    # ax.plot(dfSpectrum['Wavelength'][rangeBoolIntens],dfSpectrum['Absorbance'][rangeBoolIntens])

# calculates the mean and standard deviation for gaussian lambda max and intensity values
gauss_avg_lam_max = st.fmean(gauss_lam_max)
gauss_avg_intensity = st.fmean(gauss_intensity)
gauss_lam_dev = st.stdev(gauss_lam_max)
gauss_intensity_dev = st.stdev(gauss_intensity)
gauss_lam_variance = gauss_lam_dev ** 2
gauss_intensity_variance = gauss_intensity_dev ** 2

raw_avg_lam_max = st.fmean(raw_lam_max)
raw_avg_intensity = st.fmean(raw_intensity)
raw_lam_dev = st.stdev(raw_lam_max)
raw_intensity_dev = st.stdev(raw_intensity)
raw_lam_variance = raw_lam_dev ** 2
raw_intensity_variance = raw_intensity_dev ** 2

smooth_avg_lam_max = st.fmean(smooth_lam_max)
smooth_avg_intensity = st.fmean(smooth_intensity)
smooth_lam_dev = st.stdev(smooth_lam_max)
smooth_intensity_dev = st.stdev(smooth_intensity)
smooth_lam_variance = smooth_lam_dev ** 2
smooth_intensity_variance = smooth_intensity_dev ** 2
# writres the values to the new gaussian calculated values

row = 0

averages.write(rowSum, colSum, 'Name')
averages.write(rowSum, colSum + 1, gauss_avg_lam_max)
averages.write(rowSum, colSum + 2, gauss_avg_intensity)
averages.write(rowSum, colSum + 3, gauss_lam_dev)
averages.write(rowSum, colSum + 4, gauss_intensity_dev)
averages.write(rowSum, colSum + 5, gauss_lam_variance)
averages.write(rowSum, colSum + 6, gauss_intensity_variance)

averages.write(rowSum, colSum + 7, raw_avg_lam_max)
averages.write(rowSum, colSum + 8, raw_avg_intensity)
averages.write(rowSum, colSum + 9, raw_lam_dev)
averages.write(rowSum, colSum + 10, raw_intensity_dev)
averages.write(rowSum, colSum + 11, raw_lam_variance)
averages.write(rowSum, colSum + 12, raw_intensity_variance)

averages.write(rowSum, colSum + 13, smooth_avg_lam_max)
averages.write(rowSum, colSum + 14, smooth_avg_intensity)
averages.write(rowSum, colSum + 15, smooth_lam_dev)
averages.write(rowSum, colSum + 16, smooth_intensity_dev)
averages.write(rowSum, colSum + 17, smooth_lam_variance)
averages.write(rowSum, colSum + 18, smooth_intensity_variance)
averages.set_column(0, 7, 30)
rowSum += 1

# plotting lambdamax vs refractive index
# refractPlot.add_series({
#          'categories': ['Refractive Index', 2, 1, len(filestoprocess), 1],
#          'values': ['Refractive Index', 2, 2, len(filestoprocess), 2]
# })

# adjusting axis for refract plot
# refractPlot.set_y_axis({
#   'crossing': '-2',
#    'major_gridlines': {
#            'visible': False},
#    'name': 'Lambda Max (nm)',
#    'min': 500,
#    'max': 580,
#    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
# })
# refractPlot.set_x_axis({
#    'crossing': '-2',
#    'major_gridlines': {
#            'visible': False},
#    'name': 'Refractive Index',
#   'min': 0.90,
#  'max': 1.5,
# 'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
# })

# adjusting axis for rawDataPlot
rawDataPlot.set_y_axis({
    'crossing': '-2',
    'major_gridlines': {
        'visible': False},
    'name': 'Extinction',
    'min': 0,
    'max': 0.2,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})
rawDataPlot.set_x_axis({
    'crossing': '-2',
    'major_gridlines': {
        'visible': False},
    'name': 'Wavelength (nm)',
    'min': 400,
    'max': 1000,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})

# adjusting axis for smoothedDataPlot
smoothedDataPlot.set_y_axis({
    'crossing': '-2',
    'major_gridlines': {
        'visible': False},
    'name': 'Extinction',
    'min': 0,
    'max': 0.2,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})
smoothedDataPlot.set_x_axis({
    'crossing': '-2',
    'major_gridlines': {
        'visible': False},
    'name': 'Wavelength (nm)',
    'min': 400,
    'max': 1000,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})

# adjusting axis for Gaussian Fit plot
GaussDataPlot.set_y_axis({
    'crossing': '-2',
    'major_gridlines': {
        'visible': False},
    'name': 'Extinction',
    'min': 0,
    'max': 0.2,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})
GaussDataPlot.set_x_axis({
    'crossing': '-2',
    'major_gridlines': {
        'visible': False},
    'name': 'Wavelength (nm)',
    'min': 400,
    'max': 1000,
    'num_font': {'name': 'Calibri', 'size': 10, 'bold': True}
})

# inserting charts
rawDataPlotSheet.insert_chart('B1', rawDataPlot)
smoothedDataPlotSheet.insert_chart('B1', smoothedDataPlot)
# refractPlotSheet.insert_chart('B1', refractPlot)
GaussDataPlotSheet.insert_chart('B1', GaussDataPlot)

# closing workbook
workbook.close()

# making legend
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# ax.set_xlabel('Wavelength (nm)')
# ax.set_ylabel('Absorbance')
# ax.set_xlim([400,1000])
# ax.set_ylim([0,0.1])
