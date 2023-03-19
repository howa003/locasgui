# -*- coding: utf-8 -*-
from numpy import delete, s_, matrix, sqrt, power, array, append, zeros, ones, dot, reshape, amax, linalg, absolute, log10, amin, where, arange, shape, save, argmax, argmin, asarray, savetxt
import matplotlib.pyplot as plt
import gif
import datetime
from decimal import Decimal
from pandas import read_excel, DataFrame
from pathlib import Path
from scipy import interpolate

import os.path

plt.rcParams.update({'figure.max_open_warning': 0})

def hlavni_vypocet(Tstavba, Ti0, Te, duration, concrThick, steelThick, polomerVnitrni, sigmaP0, plochaPredp, density0, waterCont, thermExpan, modulusConc, modulusSteel, emiss, Lair, pe, Le, dt1, dt2, dt3, dt4, dt5):

    ####################################################################
    ############################## FUNKCE ##############################
    ####################################################################

    ############################### PLOT ###############################

    # Crop funkce potřebné pro sprváné vykreslování čar pro ocel a beton
    def cropVectSteelIn(fullVect):
        if (steelThick != 0):
            cropedVect = fullVect[:indSl]
        else:
            cropedVect = []
        return cropedVect

    def cropVectConcr(fullVect):
        if ((steelThick != 0) and (steelThickOut != 0)):
            cropedVect = fullVect[indSl:indSl2]
        elif (steelThick != 0):
            cropedVect = fullVect[indSl:]
        elif (steelThickOut != 0):
            cropedVect = fullVect[:indSl2]
        else:
            cropedVect = fullVect
        return cropedVect

    def cropVectSteelOut(fullVect):
        if (steelThickOut != 0):
            cropedVect = fullVect[indSl2:]
        else:
            cropedVect = []
        return cropedVect

    # Frame pro GIF soubor
    @gif.frame
    def plotGifFrame(x, temp, airTempVal, strnT, strnR, strnD, strssSE, strssFE, strssD, ovrprssStrn, prestrssStrn, time):
        plt.figure(figsize=(20,10))

        plt.subplot(2, 2, 1)
        plt.title('Temperature in '+str(datetime.timedelta(seconds=time)))
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(x, temp, 'b-', label = 'Temperature: <' + numToStr1dec(amin(temp)) + ',' + numToStr1dec(amax(temp)) + '> °C')
        plt.plot(0, airTempVal, 'ro', label = 'Air temperature ' + numToStr1dec(airTempVal) + ' °C')
        plt.ylabel('Temperature')
        plt.xlabel('Thickness [mm]')
        plt.xlim((0, length*1000))
        plt.ylim((0, 100))
        plt.xticks(arange(0, int(length*1000), 100))
        plt.legend(loc='lower right')

        plt.subplot(2, 2, 2)
        plt.title('Strain in '+str(datetime.timedelta(seconds=time)))
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        yAxisMin = int((10000*strnMin/1)-1)/10000
        yAxisMax = int((10000*strnMax/1)+1)/10000
        plt.plot(x, strnT, 'k--', label = 'Theor. temperature strain: <' + numToStr3dec(0) + ' %,' + numToStr3dec(100*amax(strnT)) + ' %>')
        if (pretlak == 1):
            plt.plot(x, ovrprssStrn, linestyle=(0, (3, 3, 1, 3)), color='k', label = 'Theor. overpressure strain: <' + numToStr3dec(100*amin(ovrprssStrn)) + ' %,' + numToStr3dec(100*amax(ovrprssStrn)) + ' %>')
        if (predpeti == 1):
            plt.plot(x, prestrssStrn, linestyle=(0, (3, 3, 1, 3, 1, 3)), color='k', label = 'Theor. prestress strain: <' + numToStr3dec(100*amin(prestrssStrn)) + ' %,' + numToStr3dec(100*amax(prestrssStrn)) + ' %>')
        plt.plot(x, 0*x, 'r-', label = 'Strain for fixed-elong: 0 %,')
        plt.plot(x, strnD, 'm-', label = 'Strain for free-elong: ' + numToStr3dec(100*amin(strnD)) + ' %,')
        plt.plot(x, strnR, 'b-', label = 'Strain for free-end: <' + numToStr3dec(100*amin(strnR)) + ' %,' + numToStr3dec(100*amax(strnR)) + ' %>')
        plt.ylabel('Strain')
        plt.xlabel('Thickness [mm]')
        plt.xlim((0, length*1000))
        plt.ylim((yAxisMin, yAxisMax))
        plt.xticks(arange(0, int(length*1000), 100))
        plt.legend(loc='lower right')

        plt.subplot(2, 2, 3)
        plt.title('Stress in '+str(datetime.timedelta(seconds=time)))
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        yAxisMin = 5*int((strssMin/5)-1)
        yAxisMax = 5*int((strssMax/5)+1)
        plt.plot(x, strssFE, 'r-', label = 'Stress for fixed-end: <' + numToStr1dec(amin(strssFE)) + ' ,' + numToStr1dec(amax(strssFE)) + '>')
        plt.plot(x, strssD, 'm-', label = 'Stress for free-elong: <' + numToStr1dec(amin(strssD)) + ' ,' + numToStr1dec(amax(strssD)) + '>')
        plt.plot(x, strssSE, 'b-', label = 'Stress for free-end: <' + numToStr1dec(amin(strssSE)) + ' ,' + numToStr1dec(amax(strssD)) + '>')
        plt.ylabel('Stress [MPa]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((0, length*1000))
        plt.ylim((yAxisMin,yAxisMax))
        plt.xticks(arange(0, int(length*1000), 100))
        plt.yticks(arange(yAxisMin, yAxisMax, 10))
        plt.legend(loc='lower right')

        plt.subplot(2, 2, 4)
        plt.title('Stress in concrete in '+str(datetime.timedelta(seconds=time)))
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        x1 = cropVectConcr(x)
        y1 = cropVectConcr(strssFE)
        y2 = cropVectConcr(strssD)
        y3 = cropVectConcr(strssSE)
        yAxisMin = 1*int((strssMinC/1)-1)
        yAxisMax = 1*int((strssMaxC/1)+1)
        plt.plot(x1, y1, 'r-', label = 'Stress for fixed-end: <' + numToStr3dec(amin(y1)) + ' ,' + numToStr3dec(amax(y1)) + '>')
        plt.plot(x1, y2, 'm-', label = 'Stress for free-elong: <' + numToStr3dec(amin(y2)) + ' ,' + numToStr3dec(amax(y2)) + '>')
        plt.plot(x1, y3, 'b-', label = 'Stress for free-end: <' + numToStr3dec(amin(y3)) + ' ,' + numToStr3dec(amax(y3)) + '>')
        plt.ylabel('Stress [MPa]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((0, length*1000))
        plt.ylim((yAxisMin,yAxisMax))
        plt.xticks(arange(0, int(length*1000), 100))
        plt.yticks(arange(yAxisMin, yAxisMax, 1))
        plt.legend(loc='lower right')

    # Plot průběhu teplot v konstrukci v daných časech
    def plotTemp3fig():
        figLbl = 'Temperature distribution when maximal temperature in air/concrete/steel is reached'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        # >>>>> x-axis declaration <<<<<
        xAxisVal = xAxisThick
        xAxisMax = max(xAxisVal)
        xAxisVal1 = cropVectSteelIn(xAxisVal)
        xAxisVal2 = cropVectConcr(xAxisVal)
        xAxisVal3 = cropVectSteelOut(xAxisVal)

        # >>>>> y-axis global maximum <<<<<
        yAxisMax = 100*int((max(max(vectTmaxTa),max(vectTmaxC),max(vectTmaxS))/100)+1)

        # Plot temperature distribution for time when maximal air temperature is reached
        plt.subplot(1, 3, 1)
        timeStep = maxTempAirStep
        plotTime = xAxisTime[timeStep]
        plt.title('Maximal temp. in air (in '+str(datetime.timedelta(seconds=plotTime))+')')

        yAxisVal = vectTmaxTa
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(0, airTempVect[timeStep], 'ko', label = 'Air temp. ' + str(numToStr1dec(airTempVect[timeStep])) + ' °C')
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Concrete temp.: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> °C')
        if (steelThick != 0):
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Inner steel temp.: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> °C')
        if (steelThickOut != 0):
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Outer steel temp.: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> °C')
        plt.ylabel('Temperature [°C]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((0, yAxisMax))
        plt.legend(loc='upper right')

        # Plot temperature distribution for time when maximal concrete temperature is reached
        plt.subplot(1, 3, 2)
        timeStep = maxTempCstep
        plotTime = xAxisTime[timeStep]
        plt.title('Maximal temp. in concrete (in '+str(datetime.timedelta(seconds=plotTime))+')')

        yAxisVal = vectTmaxC
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(0, airTempVect[timeStep], 'ko', label = 'Air temp. ' + str(numToStr1dec(airTempVect[timeStep])) + ' °C')
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Concrete temp.: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> °C')
        if (steelThick != 0):
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Inner steel temp.: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> °C')
        if (steelThickOut != 0):
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Outer steel temp.: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> °C')
        plt.ylabel('Temperature [°C]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((0, yAxisMax))
        plt.legend(loc='upper right')

        # Plot temperature distribution for time when maximal steel temperature is reached
        plt.subplot(1, 3, 3)
        timeStep = maxTempSstep
        plotTime = xAxisTime[timeStep]
        plt.title('Maximal temp. in steel (in '+str(datetime.timedelta(seconds=plotTime))+')')

        yAxisVal = vectTmaxS
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(0, airTempVect[timeStep], 'ko', label = 'Air temp. ' + str(numToStr1dec(airTempVect[timeStep])) + ' °C')
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Concrete temp.: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> °C')
        if (steelThick != 0):
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Inner steel temp.: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> °C')
        if (steelThickOut != 0):
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Outer steel temp.: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> °C')
        plt.ylabel('Temperature [°C]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((0, yAxisMax))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_Temp'+'.png')



    # Plot průběhu napětí v konstrukci v daných časech
    def plotStrss3fig(label,strssFxE,stepFxE,strssFrL,stepFrL,strssFrE,stepFrE,strssType):
        # >>>>> Plot declaration <<<<<
        print('Plotting: ' + label)
        plt.figure(figsize=(16,7), dpi=300)
        plt.suptitle(label, fontsize=16)


        # >>>>> Label selection and min/max value calculation <<<<<
        if (strssType == 'T'):
            nameLabel = 'StrssCncrTnsn'
        elif (strssType == 'C'):
            nameLabel = 'StrssCncrCmprsn'
        elif (strssType == 'ST'):
            nameLabel = 'StrssStlTnsn'
        elif (strssType == 'SC'):
            nameLabel = 'StrssStlCmprsn'
        elif (strssType == 'PRESURE'):
            nameLabel = 'StrssMaxPreStrss'
        elif (strssType == 'TIME0'):
            nameLabel = 'StrssInT0'
        elif (strssType == 'TIME0temp'):
            nameLabel = 'StrssInT0temp'
        elif (strssType == 'FIXED'):
            nameLabel = 'StrssAtTime'

        yAxisMin = int((min(min(strssFxE),min(strssFrL),min(strssFrE))/1)-1)
        yAxisMax = int((max(max(strssFxE),max(strssFrL),max(strssFrE))/1)+1)

        # >>>>> x-axis declaration <<<<<
        xAxisVal = xAxisThick
        xAxisMax = max(xAxisVal)
        xAxisVal1 = cropVectSteelIn(xAxisVal)
        xAxisVal2 = cropVectConcr(xAxisVal)
        xAxisVal3 = cropVectSteelOut(xAxisVal)

        # >>>>> Subplot FxE for time step with max stress <<<<<
        plt.subplot(1, 3, 1)
        timeStep = stepFxE
        plotTime = xAxisTime[timeStep]
        plt.title('Stress in Fixed-end structure in '+str(datetime.timedelta(seconds=plotTime)))

        yAxisVal = strssFxE
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Stress in concrete: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> MPa')
        if (steelThick != 0):
            plt.plot(array([xAxisVal1[-1],xAxisVal2[0]]), array([yAxisVal1[-1],yAxisVal2[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Stress in inner steel: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> MPa')
        if (steelThickOut != 0):
            plt.plot(array([xAxisVal2[-1],xAxisVal3[0]]), array([yAxisVal2[-1],yAxisVal3[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Stress in outer steel: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> MPa')
        plt.ylabel('Stress [MPa]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.legend(loc='lower right')

        # >>>>> Subplot FrL for time step with max stress <<<<<
        plt.subplot(1, 3, 2)
        timeStep = stepFrL
        plotTime = xAxisTime[timeStep]
        plt.title('Stress in Free-elonganation structure in '+str(datetime.timedelta(seconds=plotTime)))

        yAxisVal = strssFrL
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Stress in concrete: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> MPa')
        if (steelThick != 0):
            plt.plot(array([xAxisVal1[-1],xAxisVal2[0]]), array([yAxisVal1[-1],yAxisVal2[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Stress in inner steel: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> MPa')
        if (steelThickOut != 0):
            plt.plot(array([xAxisVal2[-1],xAxisVal3[0]]), array([yAxisVal2[-1],yAxisVal3[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Stress in outer steel: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> MPa')
        plt.ylabel('Stress [MPa]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.legend(loc='lower right')

        # >>>>> Subplot FrE for time step with max stress <<<<<
        plt.subplot(1, 3, 3)
        timeStep = stepFrE
        plotTime = xAxisTime[timeStep]
        plt.title('Stress in Free-end structure in '+str(datetime.timedelta(seconds=plotTime)))

        yAxisVal = strssFrE
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Stress in concrete: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> MPa')
        if (steelThick != 0):
            plt.plot(array([xAxisVal1[-1],xAxisVal2[0]]), array([yAxisVal1[-1],yAxisVal2[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Stress in inner steel: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> MPa')
        if (steelThickOut != 0):
            plt.plot(array([xAxisVal2[-1],xAxisVal3[0]]), array([yAxisVal2[-1],yAxisVal3[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Stress in outer steel: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> MPa')
        plt.ylabel('Stress [MPa]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.legend(loc='lower right')

        # >>>>> Save figure <<<<<
        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+nameLabel+'.png')

    # Plot časový průběh maximálních napětí v konstrukci
    def plotStrssTime(xLimit):
        figLbl = 'Evolution of maximal stresses'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef

        plt.subplot(1, 3, 1)
        plt.title('Maximal tensile stress in concrete')
        y1 = strssCtTimeFxE
        y2 = strssCtTimeFrL
        y3 = strssCtTimeFrE
        yAxisMin = 1*int((min(min(y1),min(y2),min(y3))/1)-1)
        yAxisMax = 1*int((max(max(y1),max(y2),max(y3))/1)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='dashed', color=[1, 0, 0], label = 'Fixed-end')
        plt.plot(xAxisPlot, y2, linestyle=(0, (3, 3, 1, 3)), color='m', label = 'Free-elongation')
        plt.plot(xAxisPlot, y3, linestyle=(0, (3, 3, 1, 3, 1, 3)), color=[0, 0, 1], label = 'Free-end')
        plt.ylabel('Stress [MPa]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.yticks(arange(yAxisMin, yAxisMax, 0.5))
        plt.legend(loc='upper right')

        plt.subplot(1, 3, 2)
        plt.title('Maximal compressive stress in concrete')
        y1 = strssCcTimeFxE
        y2 = strssCcTimeFrL
        y3 = strssCcTimeFrE
        yAxisMin = 1*int((min(min(y1),min(y2),min(y3))/1)-1)
        yAxisMax = 1*int((max(max(y1),max(y2),max(y3))/1)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='dashed', color=[1, 0, 0], label = 'Fixed-end')
        plt.plot(xAxisPlot, y2, linestyle=(0, (3, 3, 1, 3)), color='m', label = 'Free-elongation')
        plt.plot(xAxisPlot, y3, linestyle=(0, (3, 3, 1, 3, 1, 3)), color=[0, 0, 1], label = 'Free-end')
        plt.ylabel('Stress [MPa]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.yticks(arange(yAxisMin, yAxisMax, 1))
        plt.legend(loc='upper right')

        plt.subplot(1, 3, 3)
        plt.title('Maximal stress in steel')
        y1 = strssScTimeFxE
        y2 = strssScTimeFrL
        y3 = strssSTimeFrE
        yAxisMin = 5*int((min(min(y1),min(y2),min(y3))/5)-1)
        yAxisMax = 5*int((max(max(y1),max(y2),max(y3))/5)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='dashed', color=[1, 0, 0], label = 'Fixed-end')
        plt.plot(xAxisPlot, y2, linestyle=(0, (3, 3, 1, 3)), color='m', label = 'Free-elongation')
        plt.plot(xAxisPlot, y3, linestyle=(0, (3, 3, 1, 3, 1, 3)), color=[0, 0, 1], label = 'Free-end')
        plt.ylabel('Stress [MPa]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.yticks(arange(yAxisMin, yAxisMax, 5))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'strssMaxTime'+'_xLim'+str(xLimit)+'.png')

    # Plot časový průběh maximálních teplot v air/konstrukci
    def plotTempTime(xLimit):
        figLbl = 'Evolution of maximal temperature'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef

        y1 = airTempVect
        y2 = tempConcrMaxTime
        y3 = tempSteelMaxTime
        yAxisMax = 10*int((max(max(y1),max(y2),max(y3))/10)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='dashed', color=[1, 0, 0], label = 'Maximal temperature in inner atmosphere: <' + numToStr1dec(amin(y1)) + ',' + numToStr1dec(amax(y1)) + '> °C')
        plt.plot(xAxisPlot, y2, linestyle=(0, (3, 3, 1, 3)), color='m', label = 'Maximal temperature in concrete: <' + numToStr1dec(amin(y2)) + ',' + numToStr1dec(amax(y2)) + '> °C')
        plt.plot(xAxisPlot, y3, linestyle=(0, (3, 3, 1, 3, 1, 3)), color=[0, 0, 1], label = 'Maximal temperature in steel: <' + numToStr1dec(amin(y3)) + ',' + numToStr1dec(amax(y3)) + '> °C')
        plt.ylabel('Temperature [°C]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((0, yAxisMax))
        if (xLimit <= 600):
            plt.xticks(arange(0, xAxisMax, 25))
        plt.yticks(arange(0, yAxisMax, 25))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'tempMaxTime'+'_xLim'+str(xLimit)+'.png')

    # Plot časový průběh napětí od předpětí a přetlaku
    def plotNonTempTime(xLimit):
        figLbl = 'Evolution of non-temperature (i.e. overpressure and prestress) stresses in concrete'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef

        y1 = vectMaxStrssOvrprss
        y2 = vectMaxStrssPrstrss
        y3 = vectMaxStrssNontemp
        yAxisMin = 1*int((min(min(y1),min(y2),min(y3))/1)-1)
        yAxisMax = 1*int((max(max(y1),max(y2),max(y3))/1)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='dashed', color=[1, 0, 0], label = 'Overpressure induced stresses: <' + numToStr1dec(amin(y1)) + ',' + numToStr1dec(amax(y1)) + '> MPa')
        plt.plot(xAxisPlot, y2, linestyle=(0, (3, 3, 1, 3)), color=[0, 0, 1], label = 'Prestress induced stresses: <' + numToStr1dec(amin(y2)) + ',' + numToStr1dec(amax(y2)) + '> MPa')
        plt.plot(xAxisPlot, y3, linestyle='solid', color='m', label = 'Overall non-temperature induced stresses: <' + numToStr1dec(amin(y3)) + ',' + numToStr1dec(amax(y3)) + '> MPa')
        plt.ylabel('Stress [MPa]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        if (xLimit <= 600):
            plt.xticks(arange(0, xAxisMax, 25))
        plt.yticks(arange(yAxisMin, yAxisMax, 1))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'nonTempTime'+'_xLim'+str(xLimit)+'.png')

    # Plot časový průběh tlaku v air
    def plotLOCA(xLimit):
        figLbl = 'Evolution of temperature in inner atmosphere of the guard vessel'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef

        y1 = airTempVect
        yAxisMax = int(max(y1)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='solid', color=[1, 0, 0], label = 'Temperature in the inner atmosphere: <' + numToStr1dec(amin(y1)) + ',' + numToStr1dec(amax(y1)) + '> °C')
        plt.ylabel('Temperature [°C]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((0, yAxisMax))
        if (xLimit <= 600):
            plt.xticks(arange(0, xAxisMax, 25))
        plt.yticks(arange(0, yAxisMax, 25))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'tempLOCA'+'_xLim'+str(xLimit)+'.png')

    # Plot průběhu teplot ve vnitřní atmosféře
    def plotPressTime(xLimit):
        figLbl = 'Evolution of inner overpressure'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef

        y1 = airPressVect
        yAxisMax = int(max(y1)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='solid', color=[1, 0, 0], label = 'Pressure in the inner atmosphere: <' + numToStr3dec(amin(y1)) + ',' + numToStr3dec(amax(y1)) + '> MPa')
        plt.ylabel('Pressure [MPa]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((0, yAxisMax))
        if (xLimit <= 600):
            plt.xticks(arange(0, xAxisMax, 25))
        plt.yticks(arange(0, yAxisMax, 0.25))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'pressLOCA'+'_xLim'+str(xLimit)+'.png')

    # Plot průběhu napětí od předpětí v konstrukci
    def plotPresstressStress():
        figLbl = 'Stress distribution from concrete presstress'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        # >>>>> x-axis declaration <<<<<
        xAxisVal = xAxisThick
        xAxisMax = max(xAxisVal)
        xAxisVal1 = cropVectSteelIn(xAxisVal)
        xAxisVal2 = cropVectConcr(xAxisVal)
        xAxisVal3 = cropVectSteelOut(xAxisVal)

        # >>>>> y-axis declaration <<<<<
        yAxisVal = vectStrssPrstrss
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)
        yAxisMin = 5*int((min(yAxisVal)/5)-1)
        yAxisMax = 5*int((max(yAxisVal)/5)+1)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Stress in concrete: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> MPa')
        if (steelThick != 0):
            plt.plot(array([xAxisVal1[-1],xAxisVal2[0]]), array([yAxisVal1[-1],yAxisVal2[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Stress in inner steel: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> MPa')
        if (steelThickOut != 0):
            plt.plot(array([xAxisVal2[-1],xAxisVal3[0]]), array([yAxisVal2[-1],yAxisVal3[0]]), linestyle=(0, (1, 3)), color=[0.5, 0.5, 0.5])
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Stress in outer steel: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> MPa')
        plt.ylabel('Stress [MPa]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_PresstressStress'+'.png')


    def plotTempTstart():
        figLbl = 'Temperature distribution before LOCA (in t=0)'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        # >>>>> x-axis declaration <<<<<
        xAxisVal = xAxisThick
        xAxisMax = max(xAxisVal)
        xAxisVal1 = cropVectSteelIn(xAxisVal)
        xAxisVal2 = cropVectConcr(xAxisVal)
        xAxisVal3 = cropVectSteelOut(xAxisVal)

        plt.title(figLbl)

        yAxisVal = Tstart
        yAxisVal1 = cropVectSteelIn(yAxisVal)
        yAxisVal2 = cropVectConcr(yAxisVal)
        yAxisVal3 = cropVectSteelOut(yAxisVal)
        yAxisMax = 10*int((Ti0/10)+1)

        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(0, Ti0, 'ko', label = 'Inner air temp. ' + str(Ti0) + ' °C')
        plt.plot(xAxisMax, Te, 'ks', label = 'Outer air temp. ' + str(Te) + ' °C')
        plt.plot(xAxisVal2, yAxisVal2, 'b--', label = 'Concrete temp.: <' + numToStr1dec(amin(yAxisVal2)) + ',' + numToStr1dec(amax(yAxisVal2)) + '> °C')
        if (steelThick != 0):
            plt.plot(xAxisVal1, yAxisVal1, 'r-', label = 'Inner steel temp.: <' + numToStr1dec(amin(yAxisVal1)) + ',' + numToStr1dec(amax(yAxisVal1)) + '> °C')
        if (steelThickOut != 0):
            plt.plot(xAxisVal3, yAxisVal3, 'g-', label = 'Outer steel temp.: <' + numToStr1dec(amin(yAxisVal3)) + ',' + numToStr1dec(amax(yAxisVal3)) + '> °C')
        plt.ylabel('Temperature [°C]')
        plt.xlabel('Thickness [mm]')
        plt.xlim((0, xAxisMax))
        plt.ylim((20, yAxisMax))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_Tstart'+'.png')

    def plotSoucinitele():
        plt.figure(figsize=(12,6), dpi=300)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()

        xAxisMax = max(xAxisTime)
        xAxisMin = -xAxisMax/100
        yAxisMax = 5*int((max(max(htcFvect),max(htcCvect))/5)+1)
        yAxisMin = 5*int((min(min(htcFvect),min(htcCvect))/5)-1)

        plt.plot(xAxisTime, htcFvect, linestyle='dashed', color='r', label = 'Inner : <' + numToStr1dec(amin(htcFvect)) + ',' + numToStr1dec(amax(htcFvect)) + '> W/m2K')
        plt.plot(xAxisTime, htcCvect, linestyle='dashdot', color='b', label = 'Outer : <' + numToStr1dec(amin(htcCvect)) + ',' + numToStr1dec(amax(htcCvect)) + '> W/m2K')
        plt.title('Heat transfer coefficients')
        plt.ylabel('W/m2K')
        plt.xlabel('Time')
        plt.xlim((xAxisMin, xAxisMax))
        plt.ylim((yAxisMin, yAxisMax))
        plt.legend(loc='upper right')
        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_HeatTransCoef'+'.png')

    # Plot časový průběh maximálních teplot v air/konstrukci
    def plotGradient(xLimit):
        figLbl = 'Evolution of temperature gradient'
        print('Plotting: ' + figLbl)
        plt.figure(figsize=(15,6), dpi=300)
        plt.suptitle(figLbl, fontsize=16)

        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef

        y1 = tempGradVect
        yAxisMax = 5*int((max(y1)/5)+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='solid', color=[0, 128/255, 0], label = 'Temperature gradient: <' + numToStr1dec(amin(y1)) + ',' + numToStr1dec(amax(y1)) + '> °C')
        plt.ylabel('Temperature [°C]')
        plt.xlabel(xAxisLbl)
        plt.xlim((-xAxisMax/100, xAxisMax))
        plt.ylim((0, yAxisMax))
        if (xLimit <= 600):
            plt.xticks(arange(0, xAxisMax, 25))
        plt.yticks(arange(0, yAxisMax, 5))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'Gradient'+'_xLim'+str(xLimit)+'.png')


    # Plot průběhu normálové síly v čase
    def plotNormalTime(xLimit,y1,ySteel,yConcr,yType,yLbl):
        plt.figure(figsize=(15,6), dpi=300)
        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef
        figLbl = 'Evolution of the normal force in the cross-section for '+yLbl

        print('Plotting: ' + figLbl)
        plt.suptitle(figLbl, fontsize=16)

        yAxisMin = int(min(min(y1),min(ySteel),min(yConcr))-1)
        yAxisMax = int(max(max(y1),max(ySteel),max(yConcr))+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='solid', color=[1, 0, 1], label = 'Normal force total: <' + numToStr1dec(amin(y1)) + ',' + numToStr1dec(amax(y1)) + '> kN')
        plt.plot(xAxisPlot, ySteel, linestyle='dashed', color=[1, 0, 0], label = 'Normal force in steel: <' + numToStr1dec(amin(ySteel)) + ',' + numToStr1dec(amax(ySteel)) + '> kN')
        plt.plot(xAxisPlot, yConcr, linestyle='dashed', color=[0, 0, 1], label = 'Normal force in concrete: <' + numToStr1dec(amin(yConcr)) + ',' + numToStr1dec(amax(yConcr)) + '> kN')
        plt.ylabel('Normal force [kN]')
        plt.xlabel(xAxisLbl)
        plt.xlim((0, xAxisMax))
        if (abs(yAxisMax)>abs(yAxisMin)):
            plt.ylim((yAxisMin, yAxisMax))
        else:
            plt.ylim((yAxisMax, yAxisMin))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'normal'+yType+'_xLim'+str(xLimit)+'.png')

    # Plot průběhu ohybového momentu v čase
    def plotBendingTime(xLimit,y1,ySteel,yConcr,yType,yLbl):
        plt.figure(figsize=(15,6), dpi=300)
        if (xLimit <= 600):
            timeRedCoef = 1
            xAxisLbl = 'Time [s]'
        else:
            timeRedCoef = 60
            xAxisLbl = 'Time [min]'
        xAxisPlot = xAxisTime/timeRedCoef
        xAxisMax = xLimit/timeRedCoef
        figLbl = 'Evolution of the bending moment in the cross-section for '+yLbl

        print('Plotting: ' + figLbl)
        plt.suptitle(figLbl, fontsize=16)

        yAxisMin = int(min(min(y1),min(ySteel),min(yConcr))-1)
        yAxisMax = int(max(max(y1),max(ySteel),max(yConcr))+1)
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.grid()
        plt.plot(xAxisPlot, y1, linestyle='solid', color=[1, 0, 1], label = 'Bending moment total: <' + numToStr1dec(amin(y1)) + ',' + numToStr1dec(amax(y1)) + '> kNm')
        plt.plot(xAxisPlot, ySteel, linestyle='dashed', color=[1, 0, 0], label = 'Bending moment from steel: <' + numToStr1dec(amin(ySteel)) + ',' + numToStr1dec(amax(ySteel)) + '> kNm')
        plt.plot(xAxisPlot, yConcr, linestyle='dashed', color=[0, 0, 1], label = 'Bending moment from concrete: <' + numToStr1dec(amin(yConcr)) + ',' + numToStr1dec(amax(yConcr)) + '> kNm')
        plt.ylabel('Bending moment [kNm]')
        plt.xlabel(xAxisLbl)
        plt.xlim((0, xAxisMax))
        if (abs(yAxisMax)>abs(yAxisMin)):
            plt.ylim((yAxisMin, yAxisMax))
        else:
            plt.ylim((yAxisMax, yAxisMin))
        plt.legend(loc='upper right')

        plt.savefig('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'bending'+yType+'_xLim'+str(xLimit)+'.png')


    ############################ OBECNÉ FUNKCE ############################

    def rowWithMaxValue(matrix):
        rowLength = len(matrix[0]) # Délka řádku
        positFlat = argmax(matrix) # Pozice maximální hodnoty ve vektoru vzniklém zploštění matice
        return int(positFlat/rowLength)
        # return argmax(amax(matrix, axis=1))

    def rowWithMinValue(matrix):
        rowLength = len(matrix[0]) # Délka řádku
        positFlat = argmin(matrix) # Pozice maximální hodnoty ve vektoru vzniklém zploštění matice
        return int(positFlat/rowLength)
        # return argmin(amin(matrix, axis=1))

    def numToStr1dec(number):
        return str(Decimal(Decimal(number).quantize(Decimal('.1'))))

    def numToStr3dec(number):
        return str(Decimal(Decimal(number).quantize(Decimal('.001'))))

    def numToStr4dec(number):
        return str(Decimal(Decimal(number).quantize(Decimal('.0001'))))

    def kelvin(C):
        return (C+273.15)

    def prumer(a,b):
        return (a+b)/2

    def chyba(a,b):
        return abs((a/b)-1)


    #################### FUNKCE MATERIÁLOVÉ VLASTNOSTI ####################

    def conduc(T):
        if (Ti0 == 20):
            return (1.36 - 0.136*(T/100) + 0.0057*power((T/100),2)) # dolní mezní hodnota
        else:
            return (2 - 0.2451*(T/100) + 0.0107*power((T/100),2)) # horní mezní hodnota

    def density(T):
        if T < 115:
            return density0
        elif T < 200:
            return density0*(1 - (0.02*(T - 115)/85))
        elif T < 400:
            return density0*(0.98 - (0.03*(T - 200)/200))
        else:
            return density0*(0.95 - (0.07*(T - 400)/800))

    def heatCapac(T):
        cpeak = 900+(waterCont/3)*(2020-900)
        if T < 100:
            return 900
        elif T < 115:
            return cpeak
        elif T < 200:
            return (cpeak - ((T-115)/(200-115))*(cpeak-1000))
        elif T < 400:
            return (1000+(T-200)/2)
        else:
            return 1100

    def volHeatCap(T):
        return (density(T)*heatCapac(T))

    def conducSteel(T):
        if T<= 800:
            return (54 - 3.33*T/100)
        else:
            return (27.3)

    def volHeatCapSteel(T):
        rhoSt = 7850
        if T <= 600:
            return rhoSt*(425 + 7.73*T/10 - 1.69*T*T/1000 + 2.22*T*T*T/1000000)
        elif T <= 735:
            return rhoSt*(666 - (13002/(T-738)))
        elif T <= 900:
            return rhoSt*(545 + 17820/(T-731))
        else:
            return rhoSt*650

    #################### FUNKCE TOK, VODIVOST A KAPACITA ####################
    def fluxFire(Tair, Tnode):
        Tair = kelvin(Tair)
        Tnode = kelvin(Tnode)
        return clcConvection(Lair, Tnode, Tair)*(Tnode-Tair) + esbFcoef*esb*(power(Tnode,4)-power(Tair,4))

    def fluxCool(Tair, Tnode):
        Tair = kelvin(Tair)
        Tnode = kelvin(Tnode)
        return clcConvection(Lair, Tnode, Tair)*(Tnode-Tair) + esbCcoef*esb*(power(Tnode,4)-power(Tair,4))

    def fluxFireDx(Tair, Tnode):
        Tair = kelvin(Tair)
        Tnode = kelvin(Tnode)
        return clcConvection(Lair, Tnode, Tair) + esbFcoef*4*esb*(power(Tnode,3))

    def fluxCoolDx(Tair, Tnode):
        Tair = kelvin(Tair)
        Tnode = kelvin(Tnode)
        return clcConvection(Lair, Tnode, Tair) + esbFcoef*4*esb*(power(Tnode,3))

    # Vektor tepelných toků a jejich derivace
    def vectorF(tempVec,valAirTemp):
        vector = zeros((nNode))
        vector[0] = fluxFire(valAirTemp,tempVec[0])
        vector[-1] = fluxCool(Te,tempVec[-1])
        if (radialni == 1):
            vector[0] = vectElemCentAbs[0]*vector[0]
            vector[-1] = vectElemCentAbs[-1]*vector[-1]
        return vector

    def vectordF(tempVec,valAirTemp):
        vector = zeros((nNode))
        vector[0] = fluxFireDx(valAirTemp, tempVec[0])
        vector[-1] = fluxCoolDx(Te, tempVec[-1])
        if (radialni == 1):
            vector[0] = vectElemCentAbs[0]*vector[0]
            vector[-1] = vectElemCentAbs[-1]*vector[-1]
        return vector

    # Matice vodivosti K a kapacity C
    def matrixK(vectorT):
        struct = zeros((nNode,nNode))                          # matice konstrukce
        for i in elemID:                                       # Výpočet matice prvku a její lokalizace od matice konstrukce
            tempElem = (vectorT[i]+vectorT[i+1])/2
            matConduc = conduc(tempElem)                       # vodivost betonu
            if ((steelThick != 0) and (i <= (indSl-1))):
                matConduc = conducSteel(tempElem)                   # vodivost oceli
            elif ((steelThickOut != 0) and (i >= (indSl2))):
                matConduc = conducSteel(tempElem)                   # vodivost oceli
            element = matConduc * matrix([[1, -1],[-1, 1]]) / Le  # matice prvku; W/Km2
            if (radialni == 1):
                element = vectElemCentAbs[i] * element
            struct[i,i] = struct[i,i] + element[0,0]
            struct[i,i+1] = struct[i,i+1]+element[0,1]
            struct[i+1,i] = struct[i+1,i]+element[1,0]
            struct[i+1,i+1] = struct[i+1,i+1]+element[1,1]
        return struct

    def matrixC(vectorT):
        struct = zeros((nNode,nNode))                         # matice konstrukce
        for i in elemID:                                      # výpočet matice prvku a její lokalizace od matice konstrukce
            tempElem = (vectorT[i]+vectorT[i+1])/2
            matVHC = volHeatCap(tempElem)                     # objemová teplotní kapacita betonu
            if ((steelThick != 0) and (i <= (indSl-1))):
                matVHC = volHeatCapSteel(tempElem)                 # objemová teplotní kapacita oceli
            elif ((steelThickOut != 0) and (i >= (indSl2))):
                matVHC = volHeatCapSteel(tempElem)                 # objemová teplotní kapacita oceli
            element = matVHC * matrix([[1/3, 1/6],[1/6, 1/3]]) * Le # matice elementu; Ws/Km3
            if (radialni == 1):
                element = vectElemCentAbs[i] * element
            struct[i,i] = struct[i,i] + element[0,0]
            struct[i,i+1] = struct[i,i+1]+element[0,1]
            struct[i+1,i] = struct[i+1,i]+element[1,0]
            struct[i+1,i+1] = struct[i+1,i+1]+element[1,1]
        return struct

    def conducAir(Ta):
        Ta = Ta - 273.15        # °C!
        return 418.4*5.75*(1 + 0.00317*Ta - 0.0000021*Ta*Ta)/power(10,5)

    def dynVisAir(Ta):
        return 1.716*power(Ta/273.11,3/2)*((273.11 + 110.56)/(Ta + 110.56))/power(10,5)

    def densityAir(Ta):
        p = 101325      # Pa
        Rs = 287.058    # J/kgK
        return p/(Rs*Ta)

    def heatCapacAir(Ta):
        return 717.8 + 0.07075*(Ta - 300) + 0.00026125*power(Ta - 300,2)

    def clcConvection(L, Ts, Ta):
        # TEPLOTY MUSÍ PŘIJÍT V KELVINECH!
        k = conducAir(Ta)       # W/mK
        u = dynVisAir(Ta)       # Ns/m2
        p = densityAir(Ta)      # kg/m3
        cv = heatCapacAir(Ta)   # J/kgK
        g = 9.81	            # m/s2
        B = (1/Ta)              # 1/K
        v = u/p                 # m2/s
        a = k/(p*cv)            # m2/s
        Ra = (g*B/(v*a))*abs(Ts-Ta)*power(L,3)
        Pr = v/a
        hc = (k/L) * (0.68+0.67*power(Ra,1/4)/(power(1+power(0.492/Pr,9/16),4/9))) # W/m2K
        return hc

    def clcRprestup(L,emiss,sbc,Ts,Ta):
        # TEPLOTY MUSÍ PŘIJÍT V °C!
        Ts = kelvin(Ts)
        Ta = kelvin(Ta)
        hc = clcConvection(L, Ts, Ta)
        hr = emiss*sbc*(Ts+Ta)*(Ts*Ts+Ta*Ta)
        # print('   hc = ' + str(hc))
        # print('   hr = ' + str(hr))
        h = hc+hr
        return (1/h)

    def statHeatTrans():    # Stanovení rozložení teplot v konstrukci v t=0 (stacionární vedení tepla)
        RsIn = steelThick/conducSteel(Ti0)               # Odpor oceli v t=0
        RsOut = steelThickOut/conducSteel(Te)               # Odpor oceli v t=0
        Rc = concrThick/conduc(((Ti0+Te)/2))           # Odpor betonu v t=0
        mtrxStat=array([[-Ri, -1, 0, 0, 0],[-RsIn, 1, -1, 0, 0],[-Rc, 0, 1, -1, 0],[-RsOut, 0, 0, 1, -1],[-Re, 0, 0, 0, 1]])
        vctrStat = array([-Ti0, 0, 0, 0, Te])
        rsltStat = linalg.solve(mtrxStat,vctrStat)                                          # Vektor toku a teplot konstrukce
        T1 = rsltStat[1]
        T2 = rsltStat[2]
        T3 = rsltStat[3]
        T4 = rsltStat[4]
        tempInit = zeros(nNode)
        for i in nodesRange:
            if ((steelThick != 0) and (i < indSl)):
                tempInit[i] = T1-(i*Le/steelThick)*(T1-T2)
            elif ((steelThickOut != 0) and (i >= indSl2)):
                tempInit[i] = T3-((i*Le-steelThick-concrThick)/steelThickOut)*(T3-T4)
            else:
                tempInit[i] = T2-((i*Le-steelThick)/concrThick)*(T2-T3)
        return tempInit

    ####################################################################
    ############################## SKRIPT ##############################
    ####################################################################

    print('===== Program started. =====')

    # >>>>>>>>>>>>>>>>>>>>>>>>>> INPUTS <<<<<<<<<<<<<<<<<<<<<<<<<<

    # Ti0 = 41 # Teplota v interiéru před LOCA; 41 or 44
    # Te = 20 # Teplota v exteriéru před LOCA
    # Tstavba = 20            # Teplota při betonáži (před uvedením do provozu)
    # duration = 60*60*1*1   # Doba trvání LOCA [s]

    # density0 = 2500         # Objemová hmotnost betonu pro T=20°C [kg/m3]
    # waterCont = 1.5         # Obsah vody [hm. %]
    # thermExpan = 12/1000000
    # modulusConc = 35000     # MPa
    # modulusSteel = 210000   # MPa

    # concrThick = 1.5        # m; Tloušťka betonu
    # steelThick = 0.008      # m; Tloušťka ocelové vystýlky
    steelThickOut = 0  # m; Tloušťka ocelové vnější vrstvy
    # polomerVnitrni = 19/2   # [m]; Vnitřní poloměr válce
    # sigmaP0 = 1280              # MPa; Napětí v předpínací výztuži v čase t=0
    # plochaPredp = 5*2850/1000000 # m2; Plocha výztuže na 1 metr běžný stěny (výška)
    # pe = 0.1                    # MPa; Tlak vzduchu vně válce

    # Lair = 0.05                 # Šířka vrstv vzduchu, ve které dochází k výměně tepla prouděním
    # emiss = 0.7                  # emisivita

    # Le = 0.001                        # m; Prostorový krok

    # dt1 = 1                 # s
    # dt2 = 15                # s
    # dt3 = 60                # s
    # dt4 = 300               # s
    # dt5 = 1800              # s



    # >>>>>>>>>>>>>>>>>>>>>>>>>> SETTINGS <<<<<<<<<<<<<<<<<<<<<<<<<<
    lblVersion = 'locas2gui'      # Pojmenování výstupních souborů

    provozniNapeti = 1      # Má být uvažován vznik napětí od provozních teplot? 1 (ano), 0 (ne) ... pureTemp = 0; teplotaTlakPredpeti = 1
    predpeti = 1            # Uvažovat předpětí? 1 (ano) or 0 (ne) ... pureTemp = 0; teplotaTlakPredpeti = 1
    pretlak = 1             # Uvažovat přetlak? 1 (ano) or 0 (ne) ... pureTemp = 0; teplotaTlakPredpeti = 1

    radialni = 1            # Má se použít vedení tepla pro radiální směr? (1-ano; 0-ne).. 0 značí klasické 1D vedení

    # Koeficienty přestupu tepla
    sbc = 5.670374419/power(10,8)  # Stefan-Boltzmannova konstanta
    stefBolt = 5.67/power(10,8)  # Stefan-Boltzmannova konstanta
    esb = emiss*stefBolt    # emisivita x Stefan-Boltzmannova konstanta [W/(m2K4)]
    esbFcoef = 1            # Uvažovat sálání požáru? 0 ne; 1 ano;
    esbCcoef = 1            # Uvažovat sálání nevystaveného povrchu? 0 ne; 1 ano;

    printPercent = 1       # Hlášení progresu výpočtu (menší číslo = častěji)




    # Geometrie
    length = concrThick + steelThick + steelThickOut  # m; Celková tloušťka (ocel+beton)
    modulusTotal = (modulusConc*concrThick+modulusSteel*(steelThick+steelThickOut))/length

    # Předpětí a Přetlak
    predpinaciNapeti = sigmaP0*predpeti
    prssSftCoeff = 1.5*pretlak
    lblType = ''
    prestressStrss = (-predpinaciNapeti*plochaPredp)/(1*length-plochaPredp) # MPa; Napětí v betonu od předpětí

    lblVersion = lblVersion + '-T' + str(Tstavba)
    lblVersion = lblVersion + '-P' + str(prssSftCoeff)
    lblVersion = lblVersion + '-' + str(provozniNapeti) + str(predpeti) + str(pretlak)
    lblVersion = lblVersion + '-' + str(Ti0)

    # Prostorová diskretizace
    nElem = int(length/Le)                      # Počet elementů
    nNode = nElem+1                             # Počet uzlů
    elemID = array([i for i in range(nElem)])   # Čísla elementů
    nodeID = append(elemID, elemID[-1]+1)       # Čísla uzlů
    xAxisThick = 1000*nodeID*Le

    nodesRange = range(nNode)
    nodesPosit = array([(i*Le) for i in nodesRange])
    elemRange = range(nElem)
    elemCenters = array([(i*Le+Le/2) for i in elemRange])
    vectElemCentAbs = elemCenters+polomerVnitrni

    indSl = int(((steelThick/Le)+1))                # Slice index
    indSl2 = int((((length-steelThickOut)/Le)))     # Slice index 2

    # Časová diskretizace
    t1 = min(60*5,duration)               # s
    t2 = min(60*30,duration)              # s
    t3 = min(60*60*2,duration)            # s
    t4 = min(60*60*24*1,duration)         # s
    t5 = duration                         # s
    if (duration <= t1):
        steps = int(t1/dt1)
    elif (duration <= t2):
        steps = int((t1/dt1)+((t2-t1)/dt2))
    elif (duration <= t3):
        steps = int((t1/dt1)+((t2-t1)/dt2)+((t3-t2)/dt3))
    elif (duration <= t4):
        steps = int((t1/dt1)+((t2-t1)/dt2)+((t3-t2)/dt3)+((t4-t3)/dt4))
    else:
        steps = int((t1/dt1)+((t2-t1)/dt2)+((t3-t2)/dt3)+((t4-t3)/dt4)+((t5-t4)/dt5))
    xAxisTime = array([])
    for i in range(steps+1):
        if (len(xAxisTime) == 0):
            xAxisTime = append(xAxisTime,0)
        elif (xAxisTime[-1] < t1):
            xAxisTime = append(xAxisTime,xAxisTime[-1]+dt1)
        elif (xAxisTime[-1] < t2):
            xAxisTime = append(xAxisTime,xAxisTime[-1]+dt2)
        elif (xAxisTime[-1] < t3):
            xAxisTime = append(xAxisTime,xAxisTime[-1]+dt3)
        elif (xAxisTime[-1] < t4):
            xAxisTime = append(xAxisTime,xAxisTime[-1]+dt4)
        elif (xAxisTime[-1] < t5):
            xAxisTime = append(xAxisTime,xAxisTime[-1]+dt5)
    xAxisTimeMax = int(duration)



    T1old = Ti0
    T4old = Te
    while True:
        Ri = clcRprestup(Lair,emiss,sbc,T1old,Ti0)
        Re = clcRprestup(Lair,emiss,sbc,T4old,Te)
        Tstart = statHeatTrans()
        T1new = Tstart[1]
        T4new = Tstart[-1]
        if ((chyba(T1new,T1old)<0.001) and (chyba(T4new,T4old)<0.001)):
            break
        T1old = T1new
        T4old = T4new

    if (provozniNapeti):
        Tzaklad = ones(nNode)*Tstavba   # Základní teplota (bez napětí) je teplota před uvedením do provozu
    else:
        Tzaklad = Tstart                # Základní teplota (bez napětí) je teplota před LOCA (tj. před LOCA není žádné napětí od teplot)


    # Průřez
    beznyMetr = float(1)                                # m; rozměr kolmý na řešený průřez
    ratioE = modulusSteel/modulusConc
    bTot = ratioE*beznyMetr
    b1 = beznyMetr
    h1 = steelThick
    h2 = length-(steelThick+steelThickOut)
    h3 = steelThickOut
    a1 = h1*bTot
    a2 = h2*b1
    a3 = h3*bTot
    c1 = h1/2
    c2 = h1+h2/2
    c3 = h1+h2+h3/2
    aTot = a1+a2+a3
    cTot = (a1*c1+a2*c2+a3*c3)/aTot
    i1 = bTot*h1*h1*h1/12
    i2 = b1*h2*h2*h2/12
    i3 = bTot*h3*h3*h3/12
    d1 = a1*power((cTot-c1), 2)
    d2 = a2*power((cTot-c2), 2)
    d3 = a3*power((cTot-c3), 2)
    inertiaNahr = i1+i2+i3+d1+d2+d3             # m4; moment setrvačnosti
    areaNahr = aTot                             # m2; plocha náhradního průřezu
    center = cTot                               # m; těžiště



    # Popis výstupů
    lblTemp = str(Ti0)
    if (duration <= 60*60):
        lblTime = str(int(duration/60))+'m'
    else:
        lblTime = str(int(duration/(60*60)))+'h'


    # Funkce popisující teplotu uvnitř obálky
    nazevSouboru = 'ujvAirTepl'+str(Ti0)+'C.xlsx'       # Název souboru
    if not os.path.isfile(nazevSouboru):
        print('ERROR: File temperature.xlsx is missing in the folder.')
        return None
    excelData = read_excel(nazevSouboru,header=None)    # Načtení souboru
    vectCas = excelData.iloc[:,0].to_numpy()            # Vytvoření vektoru hodnot osy x (čas)
    vectHodn = excelData.iloc[:,1].to_numpy()           # Vytvoření vektoru hodnot osy y (čas)
    airTemp = interpolate.interp1d(vectCas, vectHodn)   # Vytvoření funkce airTemp(time) popisující teplotu v závislosti na čase

    # Funkce popisující tlak uvnitř obálky
    nazevSouboru = 'ujvAirTlak'+str(Ti0)+'C.xlsx'       # Název souboru
    if not os.path.isfile(nazevSouboru):
        print('ERROR: File pressure.xlsx is missing in the folder.')
        return None
    excelData = read_excel(nazevSouboru,header=None)    # Načtení souboru
    vectCas = excelData.iloc[:,0].to_numpy()            # Vytvoření vektoru hodnot osy x (čas)
    vectHodn = excelData.iloc[:,1].to_numpy()           # Vytvoření vektoru hodnot osy y (čas)
    fnctAirTlak = interpolate.interp1d(vectCas, vectHodn)   # Vytvoření funkce airTemp(time) popisující teplotu v závislosti na čase



    ############### VÝPOČET ###############
    print('Calculation version: '+lblVersion+'.')
    print('===== Temperature started. =====')
    airTempVect = zeros((steps+1))
    tempGradVect = zeros((steps+1))

    airTempVect[0] = Ti0
    tempGradVect[0] = Tstart[0]-Tstart[-1]

    airPressVect = zeros((steps+1))
    airPressVect[0] = fnctAirTlak(0)    # Tady není součinitel prssSftCoeff, protože před LOCA se uvažuje skutečný tlak - ne zvýšený.

    htcFvect = zeros((steps+1))
    htcCvect = zeros((steps+1))
    htcFvect[0] = 1/Ri
    htcCvect[0] = 1/Re

    matrixT = zeros((steps+1, nNode))           # Inicializace matice teplot v čase
    matrixT[0] = Tstart                         # Vložení teplot v t=0

    time = 0
    dt = dt1



    for i in range(steps):
        # 0) Začínáme na i=0.
        # 1) V daném časovém kroku (nyní) stanovíme teploty Tn a matice K a C.
        # 1) Stanovíme čas a teplotu požáru v budoucím časovém kroku
        # 2) Do Tnts (Temperature next time step) vkládáme jako první odhad teploty v t=0.
        # 4) Vypočítáme F a dF pro teploty v budoucím časovém kroku.
        # 5) Určíme residuum (chybu).
        # 6) Pokud je chyba velká, uděláme úpravu Tnts (Tnts = Tnts - R/dR) a přepočítáme 4) a 5).
        # 7) Pokračujeme na i=i+1.

        Tn = matrixT[i]          # tempertature now; teploty nyní (v čase i)
        K = matrixK(Tn)          # matice vodivosti nyní (v čase i)
        C = matrixC(Tn)          # matice kapacity nyní (v čase i)

        if (time < t1):
            dt = dt1
        elif (time < t2):
            dt = dt2
        elif (time < t3):
            dt = dt3
        elif (time < t4):
            dt = dt4
        elif (time < t5):
            dt = dt5

        time = time+dt

        if (i % (int(steps/100)*printPercent) == 0):
            print('[Temperature] Time: '+str(int(time/6)/10)+' minutes ('+str(int(100*i/steps))+' %)')
            eel.print_status('[Temperature] Time: ' + str(int(time / 6) / 10) + ' minutes (' + str(int(100 * i / steps)) + ' %)')()
        ftrAirTemp = airTemp(time)      # teplota požáru v budoucím časovém kroku i+1
        airTempVect[i+1] = ftrAirTemp
        airPressVect[i+1] = prssSftCoeff*fnctAirTlak(time)

        htcFvect[i+1] = 1/clcRprestup(Lair,emiss,sbc,Tn[0],ftrAirTemp)
        htcCvect[i+1] = 1/clcRprestup(Lair,emiss,sbc,Tn[-1],Te)

        Tnts = matrixT[i]               # první odhad budoucích teplot (teploty v čase i+1)
        F = vectorF(Tnts,ftrAirTemp)    # vektor toků v čase i+1
        dF = vectordF(Tnts,ftrAirTemp)             # derivace vektoru toků v čase i+1

        R = C.dot((Tnts-Tn)/dt) + K.dot(Tnts) + F # vektor residua
        dR = C/dt + K + dF                        # derivace vektoru residua
        maximum = amax(absolute(R))               # největší hodnodnota ve vektoru residua

        pocitadlo = 0
        while maximum > (1/1000):
            dX = linalg.solve(dR,R) # změna teplot Tnts (Newton method)
            Tnts = Tnts - dX        # oprava teplot Tnts

            F = vectorF(Tnts,ftrAirTemp)       # vektor toků v čase i+1
            dF = vectordF(Tnts,ftrAirTemp)     # derivace vektoru toků v čase i+1

            R = C.dot((Tnts-Tn)/dt) + K.dot(Tnts) + F
            dR = C/dt + K + dF
            maximum = amax(absolute(R))
        matrixT[i+1] = Tnts           # Uložení teplot v čase i+1 do matice teplot
        tempGradVect[i+1] = Tnts[0]-Tnts[-1]
    print('===== Temperature ended. =====')
    eel.print_status('===== Temperature ended. =====')()


    eel.print_status('===== Stress started. =====')()
    print('===== Stress started. =====')
    matrixStrnF = zeros((steps+1, nNode))
    matrixStrnR = zeros((steps+1, nNode))
    matrixStrnD = zeros((steps+1, nNode))   # Strain for free displacement (and fixed rotation)
    matrixStrssF = zeros((steps+1, nNode))
    matrixStrss = zeros((steps+1, nNode))
    matrixStrssD = zeros((steps+1, nNode))   # Stress for free displacement (and fixed rotation)

    for i in range(steps+1):
        temp = matrixT[i]           # Vektor teplot v case t
        tempDiff = zeros(nNode)     # Vektor zmen teplot od t=0
        for j in nodesRange:
            tempDiff[j] = (temp[j] - Tzaklad[j])

        freeStrains = zeros(nNode)      # Volné pretvoreni od teploty
        fixedStresses = zeros(nNode)     # Volné napeti od teploty
        for j in nodesRange:
            temperature = tempDiff[j]
            freeStrain = temperature*thermExpan
            freeStrains[j] = freeStrain
            fixedStresses[j] = -freeStrain*modulusConc
            if ((steelThick != 0) and (j < indSl)):
                fixedStresses[j] = -freeStrain*modulusSteel
            elif ((steelThickOut != 0) and (j >= indSl2)):
                fixedStresses[j] = -freeStrain*modulusSteel

        # Normalova sila a Moment
        normalForce = 0 # MN
        bendingMoment = 0 # MNm
        for j in elemRange:
            stress = (fixedStresses[j]+fixedStresses[j+1])/2 # Pro Fixed-end
            force = stress*Le*beznyMetr
            normalForce = normalForce + force
            bendingMoment = bendingMoment + force*(center-elemCenters[j])

        # Skutecne pretvoreni a krivost
        strain = -normalForce/(modulusConc*areaNahr) # - = MN / (MPa*m2)
        curvature = -bendingMoment/(modulusConc*inertiaNahr) # 1/m = MNm / (MPa*m4)

        # Real strain (free-end a free-displacement)
        realStrains = zeros(nNode)
        realStrainsDispl = zeros(nNode)
        for j in nodesRange:
            position = j*Le
            realStrains[j] = strain - curvature*(position-center)   # free-end
            realStrainsDispl[j] = strain                            # free-displacement

        # Real stress
        realStress = zeros(nNode)
        realStressDispl = zeros(nNode)
        for j in nodesRange:
            strainDiff = realStrains[j] - freeStrains[j]
            strainDispl = realStrainsDispl[j] - freeStrains[j]  # Strain difference for free displacement (and fixed rotation)
            realStress[j] = strainDiff * modulusConc
            realStressDispl[j] = strainDispl * modulusConc
            if ((steelThick != 0) and (j < indSl)):
                realStress[j] = strainDiff * modulusSteel
                realStressDispl[j] = strainDispl * modulusSteel
            elif ((steelThickOut != 0) and (j >= indSl2)):
                realStress[j] = strainDiff * modulusSteel
                realStressDispl[j] = strainDispl * modulusSteel





        # Lokalizace do matic
        matrixStrnF[i] = freeStrains
        matrixStrnR[i] = realStrains
        matrixStrnD[i] = realStrainsDispl
        matrixStrssF[i] = fixedStresses
        matrixStrss[i] = realStress
        matrixStrssD[i] = realStressDispl

        if (i % (int(steps/100)*printPercent) == 0):
            print('[Stress] Time: '+str(int(100*i/steps))+' %')

    matrixTempStrssF = matrixStrssF
    matrixTempStrss = matrixStrss
    matrixTempStrssD = matrixStrssD
    print('===== Stress ended. =====')




    print('===== Overpressure started. =====')
    ri = polomerVnitrni
    re = ri+length
    ri2 = ri*ri
    re2 = re*re
    matrixStrnOvrprss = zeros((steps+1, nNode))
    matrixStrssOvrprss = zeros((steps+1, nNode))
    vectStrnOvrprss = zeros((nNode))
    vectStrssOvrprss = zeros((nNode))
    vectMaxStrssOvrprss = zeros((steps+1))
    for i in range(steps+1):
        pi = airPressVect[i]
        maxStrssOvrprss = 0
        for j in range(nNode):
            rAct = ri+(j*Le)
            strssAct = ((pi*ri2-pe*re2)/(re2-ri2))+((pi-pe)*re2*ri2/((re2-ri2)*rAct*rAct))
            strnAct = strssAct/modulusTotal
            vectStrnOvrprss[j] = strnAct
            if ((steelThick != 0) and (j < indSl)):
                vectStrssOvrprss[j] = strnAct * modulusSteel
            elif ((steelThickOut != 0) and (j >= indSl2)):
                vectStrssOvrprss[j] = strnAct * modulusSteel
            else:
                vectStrssOvrprss[j] = strnAct * modulusConc
                maxStrssOvrprss = max(strnAct * modulusConc, maxStrssOvrprss)
        matrixStrnOvrprss[i] = vectStrnOvrprss
        matrixStrssOvrprss[i] = vectStrssOvrprss
        vectMaxStrssOvrprss[i] = maxStrssOvrprss

    matrixStrnOvrprss = pretlak*matrixStrnOvrprss
    matrixStrssOvrprss = pretlak*matrixStrssOvrprss
    vectMaxStrssOvrprss = pretlak*vectMaxStrssOvrprss
    print('===== Overpressure ended. =====')



    print('===== Prestress started. =====')
    # prstrssStrn = prestressStrss/modulusTotal
    pi = prestressStrss*length/((ri+re)/2)
    pe = 0
    matrixStrnPrstrss = zeros((steps+1, nNode))
    matrixStrssPrstrss = zeros((steps+1, nNode))
    vectMaxStrssPrstrss = zeros((steps+1))
    for i in range(steps+1):
        valMaxStrss = 0 # inicializace
        for j in range(nNode):
            rAct = ri+(j*Le)
            strssAct = ((pi*ri2-pe*re2)/(re2-ri2))+((pi-pe)*re2*ri2/((re2-ri2)*rAct*rAct))
            valMaxStrss = min(strssAct, valMaxStrss)
            strnAct = strssAct/modulusTotal
            matrixStrnPrstrss[i][j] = strnAct
            if ((steelThick != 0) and (j < indSl)):
                matrixStrssPrstrss[i][j] =  strnAct * modulusSteel
            elif ((steelThickOut != 0) and (j >= indSl2)):
                matrixStrssPrstrss[i][j] =  strnAct * modulusSteel
            else:
                matrixStrssPrstrss[i][j] =  strnAct * modulusConc
        vectMaxStrssPrstrss[i] = valMaxStrss

    matrixStrnPrstrss = predpeti*matrixStrnPrstrss
    matrixStrssPrstrss = predpeti*matrixStrssPrstrss
    vectMaxStrssPrstrss = predpeti*vectMaxStrssPrstrss
    print('===== Prestress ended. =====')



    print('===== Non-temperature started. =====')
    vectMaxStrssNontemp = zeros((steps+1))
    for i in range(steps+1):
        vectMaxStrssNontemp[i] = vectMaxStrssPrstrss[i]+vectMaxStrssOvrprss[i]
    print('===== Non-temperature ended. =====')



    print('===== All stresses started. =====')
    if (pretlak):
        print('Adding inner overpressure to overall stress.')
        matrixStrnR = matrixStrnR + matrixStrnOvrprss
        matrixStrnD = matrixStrnD + matrixStrnOvrprss
        matrixStrssF = matrixStrssF + matrixStrssOvrprss
        matrixStrss = matrixStrss + matrixStrssOvrprss
        matrixStrssD = matrixStrssD + matrixStrssOvrprss
    if (predpeti):
        print('Adding Prestress stress stress to overall stress.')
        matrixStrnR = matrixStrnR + matrixStrnPrstrss
        matrixStrnD = matrixStrnD + matrixStrnPrstrss
        matrixStrssF = matrixStrssF + matrixStrssPrstrss
        matrixStrss = matrixStrss + matrixStrssPrstrss
        matrixStrssD = matrixStrssD + matrixStrssPrstrss
    print('===== All stresses ended. =====')


    print('===== Vnitrni sily started. =====')
    vectNormalFxE = zeros(steps+1)
    vectBendingFxE = zeros(steps+1)
    vectNormalFrD = zeros(steps+1)
    vectBendingFrD = zeros(steps+1)
    vectNormalFrE = zeros(steps+1)
    vectBendingFrE = zeros(steps+1)

    vectNormalConcrFxE = zeros(steps+1)
    vectBendingConcrFxE = zeros(steps+1)
    vectNormalConcrFrD = zeros(steps+1)
    vectBendingConcrFrD = zeros(steps+1)
    vectNormalConcrFrE = zeros(steps+1)
    vectBendingConcrFrE = zeros(steps+1)

    vectNormalSteelFxE = zeros(steps+1)
    vectBendingSteelFxE = zeros(steps+1)
    vectNormalSteelFrD = zeros(steps+1)
    vectBendingSteelFrD = zeros(steps+1)
    vectNormalSteelFrE = zeros(steps+1)
    vectBendingSteelFrE = zeros(steps+1)

    for i in range(steps+1):
        normalFxE = 0 # MN
        bendingFxE = 0 # MNm
        normalFrD = 0 # MN
        bendingFrD = 0 # MNm
        normalFrE = 0 # MN
        bendingFrE = 0 # MNm

        normalSteelFxE = 0 # MN
        bendingSteelFxE = 0 # MNm
        normalSteelFrD = 0 # MN
        bendingSteelFrD = 0 # MNm
        normalSteelFrE = 0 # MN
        bendingSteelFrE = 0 # MNm

        normalConcrFxE = 0 # MN
        bendingConcrFxE = 0 # MNm
        normalConcrFrD = 0 # MN
        bendingConcrFrD = 0 # MNm
        normalConcrFrE = 0 # MN
        bendingConcrFrE = 0 # MNm

        fixedStresses = matrixStrssF[i] # Fixed-end napeti
        realStressDispl = matrixStrssD[i] # Free-displacement napeti
        realStress = matrixStrss[i] # Free-end napeti
        for j in elemRange:
            stressFxE = (fixedStresses[j]+fixedStresses[j+1])/2 # Pro Fixed-end
            force = stressFxE*Le*beznyMetr
            normalFxE = normalFxE + force
            bendingFxE = bendingFxE + force*(center-elemCenters[j])
            if ((steelThick != 0) and (j < indSl)):
                normalSteelFxE = normalSteelFxE + force
                bendingSteelFxE = bendingSteelFxE + force*(center-elemCenters[j])
            elif ((steelThickOut != 0) and (j >= indSl2)):
                normalSteelFxE = normalSteelFxE + force
                bendingSteelFxE = bendingSteelFxE + force*(center-elemCenters[j])
            else:
                normalConcrFxE = normalConcrFxE + force
                bendingConcrFxE = bendingConcrFxE + force*(center-elemCenters[j])

            stressFrD = (realStressDispl[j]+realStressDispl[j+1])/2 # Pro Fixed-displacement
            force = stressFrD*Le*beznyMetr
            normalFrD = normalFrD + force
            bendingFrD = bendingFrD + force*(center-elemCenters[j])
            if ((steelThick != 0) and (j < indSl)):
                normalSteelFrD = normalSteelFrD + force
                bendingSteelFrD = bendingSteelFrD + force*(center-elemCenters[j])
            elif ((steelThickOut != 0) and (j >= indSl2)):
                normalSteelFrD = normalSteelFrD + force
                bendingSteelFrD = bendingSteelFrD + force*(center-elemCenters[j])
            else:
                normalConcrFrD = normalConcrFrD + force
                bendingConcrFrD = bendingConcrFrD + force*(center-elemCenters[j])

            stressFrE = (realStress[j]+realStress[j+1])/2 # Pro Free-end
            force = stressFrE*Le*beznyMetr
            normalFrE = normalFrE + force
            bendingFrE = bendingFrE + force*(center-elemCenters[j])
            if ((steelThick != 0) and (j < indSl)):
                normalSteelFrE = normalSteelFrE + force
                bendingSteelFrE = bendingSteelFrE + force*(center-elemCenters[j])
            elif ((steelThickOut != 0) and (j >= indSl2)):
                normalSteelFrE = normalSteelFrE + force
                bendingSteelFrE = bendingSteelFrE + force*(center-elemCenters[j])
            else:
                normalConcrFrE = normalConcrFrE + force
                bendingConcrFrE = bendingConcrFrE + force*(center-elemCenters[j])


        vectNormalFxE[i] = normalFxE
        vectBendingFxE[i] = bendingFxE
        vectNormalFrD[i] = normalFrD
        vectBendingFrD[i] = bendingFrD
        vectNormalFrE[i] = normalFrE
        vectBendingFrE[i] = bendingFrE

        vectNormalSteelFxE[i] = normalSteelFxE
        vectBendingSteelFxE[i] = bendingSteelFxE
        vectNormalSteelFrD[i] = normalSteelFrD
        vectBendingSteelFrD[i] = bendingSteelFrD
        vectNormalSteelFrE[i] = normalSteelFrE
        vectBendingSteelFrE[i] = bendingSteelFrE

        vectNormalConcrFxE[i] = normalConcrFxE
        vectBendingConcrFxE[i] = bendingConcrFxE
        vectNormalConcrFrD[i] = normalConcrFrD
        vectBendingConcrFrD[i] = bendingConcrFrD
        vectNormalConcrFrE[i] = normalConcrFrE
        vectBendingConcrFrE[i] = bendingConcrFrE
    print('===== Vnitrni sily ended. =====')



    print('===== Results started. =====')
    eel.print_status('Progress: Compilation of results started.')()

    # Uložit výpočtová data
    inputs = array([])
    inputs = append(inputs,length)
    inputs = append(inputs,steelThick)
    inputs = append(inputs,steelThickOut)
    inputs = append(inputs,Le)

    dataSavePath = 'savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion
    Path(dataSavePath).mkdir(parents=True, exist_ok=True)
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/inputs.csv', inputs, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/xAxisTime.csv', xAxisTime, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/airTempVect.csv', airTempVect, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/airPressVect.csv', airPressVect, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/matrixT.csv', matrixT, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strnFxE.csv', matrixStrnF, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strnFrE.csv', matrixStrnR, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strnFrD.csv', matrixStrnD, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strssFxE.csv', matrixStrssF, delimiter=";") # Napeti celkova
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strssFrE.csv', matrixStrss, delimiter=";") # Napeti celkova
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strssFrD.csv', matrixStrssD, delimiter=";") # Napeti celkova
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strssTempFxE.csv', matrixTempStrssF, delimiter=";") # Napeti jen od teplot
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strssTempFrE.csv', matrixTempStrss, delimiter=";") # Napeti jen od teplot
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/strssTempFrD.csv', matrixTempStrssD, delimiter=";") # Napeti jen od teplot
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/vectPrestress.csv', matrixStrssPrstrss[0], delimiter=";") # Vektor rozlozeni napeti od predpeti
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/matOverpressure.csv', matrixStrssOvrprss, delimiter=";") # Napeti od pretlaku
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalFxE.csv', vectNormalFxE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalFrD.csv', vectNormalFrD, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalFrE.csv', vectNormalFrE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingFxE.csv', vectBendingFxE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingFrD.csv', vectBendingFrD, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingFrE.csv', vectBendingFrE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalSteelFxE.csv', vectNormalSteelFxE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalSteelFrD.csv', vectNormalSteelFrD, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalSteelFrE.csv', vectNormalSteelFrE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingSteelFxE.csv', vectBendingSteelFxE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingSteelFrD.csv', vectBendingSteelFrD, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingSteelFrE.csv', vectBendingSteelFrE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalConcrFxE.csv', vectNormalConcrFxE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalConcrFrD.csv', vectNormalConcrFrD, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/normalConcrFrE.csv', vectNormalConcrFrE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingConcrFxE.csv', vectBendingConcrFxE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingConcrFrD.csv', vectBendingConcrFrD, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/bendingConcrFrE.csv', vectBendingConcrFrE, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/soucinitelPrestupuFire.csv', htcFvect, delimiter=";")
    savetxt('savedData/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/soucinitelPrestupuCool.csv', htcCvect, delimiter=";")




    # Funkce pro oříznutí celkové matice a vytvoření matic pro beton a ocel
    def cropMatConcr(fullMat):
        if ((steelThick != 0) and (steelThickOut != 0)):
            cropedMat = fullMat[:,indSl:indSl2]
        elif (steelThick != 0):
            cropedMat = fullMat[:,indSl:]
        elif (steelThickOut != 0):
            cropedMat = fullMat[:,:indSl2]
        else:
            cropedMat = fullMat
        return cropedMat

    def cropMatSteel(fullMat):
        if ((steelThick != 0) and (steelThickOut != 0)):
            cropedMat = delete(fullMat, s_[indSl:indSl2], 1)
        elif (steelThick != 0):
            cropedMat = fullMat[:,:indSl]
        elif (steelThickOut != 0):
            cropedMat = fullMat[:,indSl2:]
        else:
            cropedMat = []
        return cropedMat

    # Matice teplot v betonu a v oceli
    matTempC = cropMatConcr(matrixT)
    matTempS = cropMatSteel(matrixT)

    # Matice napětí (real stress) v betonu a oceli
    matStrssFc = cropMatConcr(matrixStrssF) # for fixed end
    matStrssDc = cropMatConcr(matrixStrssD) # for free displacement
    matStrssRc = cropMatConcr(matrixStrss)  # for free end
    matStrssFs = cropMatSteel(matrixStrssF) # for fixed end
    matStrssDs = cropMatSteel(matrixStrssD) # for free displacement
    matStrssRs = cropMatSteel(matrixStrss)  # for free end

    # Kdy (v časovém jakém kroku) je dosaženo: maximální přetlak
    maxPressAirStep = argmax(airPressVect)

    # Kdy (v časovém jakém kroku) je dosaženo: maximální teploty
    maxTempAirStep = argmax(airTempVect)
    maxTempCstep = rowWithMaxValue(matTempC)
    maxTempSstep = rowWithMaxValue(matTempS)

    # Kdy (v časovém jakém kroku) je dosaženo: maximální tahy v betonu
    maxStrssFcStep = rowWithMaxValue(matStrssFc)
    maxStrssDcStep = rowWithMaxValue(matStrssDc)
    maxStrssRcStep = rowWithMaxValue(matStrssRc)

    # Kdy (v časovém jakém kroku) je dosaženo: maximální tlaky v betonu
    minStrssFcStep = rowWithMinValue(matStrssFc)
    minStrssDcStep = rowWithMinValue(matStrssDc)
    minStrssRcStep = rowWithMinValue(matStrssRc)

    # Kdy (v časovém jakém kroku) je dosaženo: maximální tahy ve výztuži
    maxStrssFsStep = rowWithMaxValue(matStrssFs)
    maxStrssDsStep = rowWithMaxValue(matStrssDs)
    maxStrssRsStep = rowWithMaxValue(matStrssRs)

    # Kdy (v časovém jakém kroku) je dosaženo: maximální tlaky ve výztuži
    minStrssFsStep = rowWithMinValue(matStrssFs)
    minStrssDsStep = rowWithMinValue(matStrssDs)
    minStrssRsStep = rowWithMinValue(matStrssRs)



    # Global min and max values
    strssMax = max(amax(matrixStrssF),amax(matrixStrss),amax(matrixStrssD))
    strssMin = min(amin(matrixStrssF),amin(matrixStrss),amin(matrixStrssD))
    strssMaxC = max(amax(matStrssFc),amax(matStrssDc),amax(matStrssRc))
    strssMinC = min(amin(matStrssFc),amin(matStrssDc),amin(matStrssRc))

    strnMax = max(amax(matrixStrnF),amax(matrixStrnR),amax(matrixStrnD),amax(matrixStrnPrstrss),amax(matrixStrnOvrprss))
    strnMin = min(amin(matrixStrnF),amin(matrixStrnR),amin(matrixStrnD),amin(matrixStrnPrstrss),amin(matrixStrnOvrprss))
    print('===== Results ended. =====')
    eel.print_status('Progress: Compilation of results ended.')()





    print('===== Plot started. =====')
    eel.print_status('Progress: Plotting of figures started.')()
    plotSavePath = 'figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion
    Path(plotSavePath).mkdir(parents=True, exist_ok=True)

    # Plot graf časový vývoj teploty ve vnitřní atmosféře
    plotLOCA(xAxisTimeMax)
    # plotLOCA(60)
    # plotLOCA(600)

    # Plot graf časový vývoj tlaku ve vnitřní atmosféře
    plotPressTime(xAxisTimeMax)
    # plotPressTime(60)
    # plotPressTime(600)





    # Plot graf (3x subplot) průběhů teplot (pro max. teploty vzduchu, betonu, oceli)
    vectTmaxTa = matrixT[maxTempAirStep]
    vectTmaxC = matrixT[maxTempCstep]
    vectTmaxS = matrixT[maxTempSstep]
    plotTemp3fig()



    # Plot graf (3x subplot) průběhů napětí (pro max. tah v betonu)
    figLabel = 'Stress distribution in the wall when maximal tension in concrete is reached'
    strssType = 'T'
    stepFxE = maxStrssFcStep
    strssFxE = matrixStrssF[stepFxE]
    stepFrL = maxStrssDcStep
    strssFrL = matrixStrssD[stepFrL]
    stepFrE = maxStrssRcStep
    strssFrE = matrixStrss[stepFrE]
    plotStrss3fig(figLabel,strssFxE,stepFxE,strssFrL,stepFrL,strssFrE,stepFrE,strssType)

    # Plot graf (3x subplot) průběhů napětí (pro max. tlak v betonu)
    figLabel = 'Stress distribution in the wall when maximal compression in concrete is reached'
    strssType = 'C'
    stepFxE = minStrssFcStep
    strssFxE = matrixStrssF[stepFxE]
    stepFrL = minStrssDcStep
    strssFrL = matrixStrssD[stepFrL]
    stepFrE = minStrssRcStep
    strssFrE = matrixStrss[stepFrE]
    plotStrss3fig(figLabel,strssFxE,stepFxE,strssFrL,stepFrL,strssFrE,stepFrE,strssType)

    # Plot graf (3x subplot) průběhů napětí (pro max. tah v oceli)
    figLabel = 'Stress distribution in the wall when maximal tension in steel is reached'
    strssType = 'ST'
    stepFxE = maxStrssFsStep
    strssFxE = matrixStrssF[stepFxE]
    stepFrL = maxStrssDsStep
    strssFrL = matrixStrssD[stepFrL]
    stepFrE = maxStrssRsStep
    strssFrE = matrixStrss[stepFrE]
    plotStrss3fig(figLabel,strssFxE,stepFxE,strssFrL,stepFrL,strssFrE,stepFrE,strssType)

    # Plot graf (3x subplot) průběhů napětí (pro max. tlak v oceli)
    figLabel = 'Stress distribution in the wall when maximal compression in steel is reached'
    strssType = 'SC'
    stepFxE = minStrssFsStep
    strssFxE = matrixStrssF[stepFxE]
    stepFrL = minStrssDsStep
    strssFrL = matrixStrssD[stepFrL]
    stepFrE = minStrssRsStep
    strssFrE = matrixStrss[stepFrE]
    plotStrss3fig(figLabel,strssFxE,stepFxE,strssFrL,stepFrL,strssFrE,stepFrE,strssType)

    # Plot graf (3x subplot) průběhů napětí (pro max. přetlak)
    figLabel = 'Stress distribution in the wall when maximal inner pressure is reached'
    strssType = 'PRESURE'
    plotStrss3fig(figLabel,matrixStrssF[maxPressAirStep],maxPressAirStep,matrixStrssD[maxPressAirStep],maxPressAirStep,matrixStrss[maxPressAirStep],maxPressAirStep,strssType)

    # Plot graf (3x subplot) průběhů napětí (pro t=0)
    figLabel = 'Stress distribution in the wall before LOCA'
    strssType = 'TIME0'
    plotStrss3fig(figLabel,matrixStrssF[0],0,matrixStrssD[0],0,matrixStrss[0],0,strssType)

    # Plot graf (3x subplot) průběhů napětí (pro t=0 a jen teplota)
    figLabel = 'Temperature-stress distribution in the wall before LOCA'
    strssType = 'TIME0temp'
    plotStrss3fig(figLabel,matrixTempStrssF[0],0,matrixTempStrssD[0],0,matrixTempStrss[0],0,strssType)



    # # Plot graf (3x subplot) průběhů napětí (pro pevně daný čas)
    # strssType = 'FIXED'
    # figLabel = 'Stress distribution in the wall at 5 minutes'
    # stepFixed = 300
    # plotStrss3fig(figLabel,matrixStrssF[stepFixed],stepFixed,matrixStrssD[stepFixed],stepFixed,matrixStrss[stepFixed],stepFixed,strssType)
    # figLabel = 'Stress distribution in the wall at 30 minutes'
    # stepFixed = 400
    # plotStrss3fig(figLabel,matrixStrssF[stepFixed],stepFixed,matrixStrssD[stepFixed],stepFixed,matrixStrss[stepFixed],stepFixed,strssType)
    # figLabel = 'Stress distribution in the wall at 2 hours'
    # stepFixed = 490
    # plotStrss3fig(figLabel,matrixStrssF[stepFixed],stepFixed,matrixStrssD[stepFixed],stepFixed,matrixStrss[stepFixed],stepFixed,strssType)

    # Plot Normal
    plotNormalTime(xAxisTimeMax,vectNormalFxE*1000,vectNormalSteelFxE*1000,vectNormalConcrFxE*1000,'FxE','fixed-end')
    plotNormalTime(xAxisTimeMax,vectNormalFrD*1000,vectNormalSteelFrD*1000,vectNormalConcrFrD*1000,'FrD','free-displacement')
    plotNormalTime(xAxisTimeMax,vectNormalFrE*1000,vectNormalSteelFrE*1000,vectNormalConcrFrE*1000,'FrE','free-end')
    plotNormalTime(int(60*60),vectNormalFxE*1000,vectNormalSteelFxE*1000,vectNormalConcrFxE*1000,'FxE','fixed-end')
    plotNormalTime(int(60*60),vectNormalFrD*1000,vectNormalSteelFrD*1000,vectNormalConcrFrD*1000,'FrD','free-displacement')
    plotNormalTime(int(60*60),vectNormalFrE*1000,vectNormalSteelFrE*1000,vectNormalConcrFrE*1000,'FrE','free-end')

    # Plot Bending
    plotBendingTime(xAxisTimeMax,vectBendingFxE*1000,vectBendingSteelFxE*1000,vectBendingConcrFxE*1000,'FxE','fixed-end')
    plotBendingTime(xAxisTimeMax,vectBendingFrD*1000,vectBendingSteelFrD*1000,vectBendingConcrFrD*1000,'FrD','free-displacement')
    plotBendingTime(xAxisTimeMax,vectBendingFrE*1000,vectBendingSteelFrE*1000,vectBendingConcrFrE*1000,'FrE','free-end')
    plotBendingTime(int(60*60),vectBendingFxE*1000,vectBendingSteelFxE*1000,vectBendingConcrFxE*1000,'FxE','fixed-end')
    plotBendingTime(int(60*60),vectBendingFrD*1000,vectBendingSteelFrD*1000,vectBendingConcrFrD*1000,'FrD','free-displacement')
    plotBendingTime(int(60*60),vectBendingFrE*1000,vectBendingSteelFrE*1000,vectBendingConcrFrE*1000,'FrE','free-end')







    # Plot graf časový vývoj maximálního tlaku v betonu, tlaku v oceli a tahu v betonu
        # (1) Časový vývoj maximálního tlaku v betonu
    strssCcTimeFxE = array([amin(matStrssFc[i]) for i in range(steps+1)])
    strssCcTimeFrL = array([amin(matStrssDc[i]) for i in range(steps+1)])
    strssCcTimeFrE = array([amin(matStrssRc[i]) for i in range(steps+1)])
        # (2) Časový vývoj maximálního tahu v betonu
    strssCtTimeFxE = array([amax(matStrssFc[i]) for i in range(steps+1)])
    strssCtTimeFrL = array([amax(matStrssDc[i]) for i in range(steps+1)])
    strssCtTimeFrE = array([amax(matStrssRc[i]) for i in range(steps+1)])
        # (3) Časový vývoj maximálního napětí v oceli
            # Časový vývoj maximálního tlaku v oceli
    strssScTimeFxE = array([amin(matStrssFs[i]) for i in range(steps+1)])
    strssScTimeFrL = array([amin(matStrssDs[i]) for i in range(steps+1)])
    strssScTimeFrE = array([amin(matStrssRs[i]) for i in range(steps+1)])
            # Časový vývoj maximálního tahu v oceli
    strssStTimeFxE = array([amax(matStrssFs[i]) for i in range(steps+1)])
    strssStTimeFrL = array([amax(matStrssDs[i]) for i in range(steps+1)])
    strssStTimeFrE = array([amax(matStrssRs[i]) for i in range(steps+1)])
            # Časový vývoj maximálního napětí v oceli
    strssSTimeFxE = strssScTimeFxE
    strssSTimeFrL = strssScTimeFrL
    strssSTimeFrE = strssScTimeFrE
    for i in range(steps+1):
        if (strssStTimeFxE[i] >= 0):
            strssSTimeFxE[i] = strssStTimeFxE[i]
        if (strssStTimeFrL[i] >= 0):
            strssSTimeFrL[i] = strssStTimeFrL[i]
        if (strssStTimeFrE[i] >= 0):
            strssSTimeFrE[i] = strssStTimeFrE[i]
        # (4) Plot
    plotStrssTime(xAxisTimeMax)
    # plotStrssTime(600)



    # Plot graf časový vývoj maximální teploty v air/concrete/steel
    tempConcrMaxTime = array([amax(matTempC[i]) for i in range(steps+1)])
    tempSteelMaxTime = array([amax(matTempS[i]) for i in range(steps+1)])
    plotTempTime(xAxisTimeMax)
    # plotTempTime(60)
    # plotTempTime(600)



    # Plot graf časový vývoj napětí od předpětí a přetlaku
    plotNonTempTime(xAxisTimeMax)
    # plotNonTempTime(60)
    # plotNonTempTime(600)







    # Plot graf časový vývoj teplotního gradientu
    plotGradient(xAxisTimeMax)



    # Plot rozložení napětí od předpětí
    vectStrssPrstrss = matrixStrssPrstrss[0]
    plotPresstressStress()



    # Plot temp distribution in t=0
    plotTempTstart()


    # Plot heat transfer coeff evolution
    plotSoucinitele()



    # Plot gif FRAME - průběh teplot, přetvoření a napětí ve stěně
    timeStep = 0
    fakeFrames = []
    fakeFrame = plotGifFrame(xAxisThick, matrixT[timeStep], airTempVect[timeStep], matrixStrnF[timeStep], matrixStrnR[timeStep], matrixStrnD[timeStep], matrixStrss[timeStep], matrixStrssF[timeStep], matrixStrssD[timeStep], matrixStrnOvrprss[timeStep], matrixStrnPrstrss[timeStep], xAxisTime[timeStep])
    fakeFrames.append(fakeFrame)
    gifName = str('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'beforeLOCA'+'.gif')
    gif.save(fakeFrames, gifName, duration=10000)

    timeStep = maxPressAirStep
    fakeFrames = []
    fakeFrame = plotGifFrame(xAxisThick, matrixT[timeStep], airTempVect[timeStep], matrixStrnF[timeStep], matrixStrnR[timeStep], matrixStrnD[timeStep], matrixStrss[timeStep], matrixStrssF[timeStep], matrixStrssD[timeStep], matrixStrnOvrprss[timeStep], matrixStrnPrstrss[timeStep], xAxisTime[timeStep])
    fakeFrames.append(fakeFrame)
    gifName = str('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'maxInnerPress'+'.gif')
    gif.save(fakeFrames, gifName, duration=10000)

    timeStep = maxTempAirStep
    fakeFrames = []
    fakeFrame = plotGifFrame(xAxisThick, matrixT[timeStep], airTempVect[timeStep], matrixStrnF[timeStep], matrixStrnR[timeStep], matrixStrnD[timeStep], matrixStrss[timeStep], matrixStrssF[timeStep], matrixStrssD[timeStep], matrixStrnOvrprss[timeStep], matrixStrnPrstrss[timeStep], xAxisTime[timeStep])
    fakeFrames.append(fakeFrame)
    gifName = str('figs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'_'+'maxInnerTemp'+'.gif')
    gif.save(fakeFrames, gifName, duration=10000)

    eel.print_status('Progress: Plotting of figures ended.')()


    # Plot gif - časová evoluce (gif) průběhů teplot ve stěně (frame)
    eel.print_status('Progress: Plotting of gif animation started.')()

    gifSavePath = 'gifs/'
    Path(gifSavePath).mkdir(parents=True, exist_ok=True)

    frames = []
    for i in range(steps+1):
        frameTime = xAxisTime[i]
        if (i % (int(steps/100)*printPercent) == 0):
            print('Plot status: '+str(int(frameTime/6)/10)+' min ('+str(int(100*i/steps))+' %).')
            eel.print_status('Animation plot status: '+str(int(frameTime/6)/10)+' min ('+str(int(100*i/steps))+' %).')()
        # if ((frameTime <= 180) or ((frameTime <= 1800) and ((frameTime % 30) == 0)) or ((frameTime <= 7200) and ((frameTime % 60) == 0)) or ((frameTime <= 172800) and ((frameTime % 300) == 0)) or ((frameTime % 1800) == 0)):
        if (True):
            vectT = matrixT[i]                      # Vektor rozložení teploty v daném čase
            airTempVal = airTempVect[i]             # Teplota ve vnitřní atmosféře v daném čase
            vectOvrprssStrn = matrixStrnOvrprss[i]  # Vektor přetvoření od přetlaku v daném čase
            vectPrestrssStrn = matrixStrnPrstrss[i]  # Vektor přetvoření od předpětí v daném čase
            strnT = matrixStrnF[i]    # Vektor Theoretical strain from temperature v daném čase
            strnR = matrixStrnR[i]    # Vektor Real strain free-end v daném čase
            strnD = matrixStrnD[i]    # Vektor Real strain free-displacement v daném čase
            strssFE = matrixStrssF[i] # Vektor Real stress fixed-end v daném čase
            strssSE = matrixStrss[i]  # Vektor Real stress free-end v daném čase
            strssD = matrixStrssD[i]  # Vektor Real stress free-displacement v daném čase
            frame = plotGifFrame(xAxisThick, vectT, airTempVal, strnT, strnR, strnD, strssSE, strssFE, strssD, vectOvrprssStrn, vectPrestrssStrn, frameTime)
            frames.append(frame)
    gifName = str('gifs/'+lblTemp+'C_'+lblType+'_'+lblTime+'_v'+lblVersion+'.gif')
    # gifName = str('gifs/'+lblTemp+'C_'+lblTime+'_v'+lblVersion+'.gif')
    gif.save(frames, gifName, duration=100)



    print('===== Plot ended. =====')
    eel.print_status('Progress: Plotting of gif animation ended.')()
    eel.print_status('Progress: Python calculation completed.')()







import eel

# initializing the application
eel.init("static_web_folder")

@eel.expose
def get_python_result(Tstavba, Ti0, Te, duration, concrThick, steelThick, polomerVnitrni, sigmaP0, plochaPredp, density0, waterCont, thermExpan, modulusConc, modulusSteel, emiss, Lair, pe, Le, dt1, dt2, dt3, dt4, dt5):
    eel.print_status('Vypocet spusten.')()
    hlavni_vypocet(Tstavba, Ti0, Te, duration, concrThick, steelThick, polomerVnitrni, sigmaP0, plochaPredp, density0, waterCont, thermExpan, modulusConc, modulusSteel, emiss, Lair, pe, Le, dt1, dt2, dt3, dt4, dt5)
    return ('COMPLETED: Calculation has ended. For results see figures and data in folders.')

# starting the application
eel.start("index.html", size=(1270, 800))