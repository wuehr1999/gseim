"""
Copyright (C) 2021 - Ruchita Korgaonkar <ruchita@iitb.ac.in>
This file is part of GSEIM.

GSEIM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import time
import numpy as np
import csv
import matplotlib
import matplotlib.pylab as plt
from os.path import expanduser
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QRadioButton,
    QFrame, QMenu, QAction, QScrollArea, QLineEdit, QFileDialog,
    QComboBox, QMainWindow, QSizePolicy, QVBoxLayout, QListWidget,
    QCheckBox, QHeaderView, QListWidgetItem, QAbstractItemView,
    QHBoxLayout, QColorDialog, QPlainTextEdit, QMessageBox,
    QTableWidget, QTableWidgetItem
    )
import numpy as np
import random
import os.path
from os import path
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    print('QT5')
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    print('Qt4')
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

class NavigationToolbar(NavigationToolbar):
    toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Home', 'Pan', 'Zoom', 'Save')]

class ScriptObject(QWidget):
    def __init__(self):
        super().__init__()
        self.title = "generatePlot.py"
        self.InitWindow()
 
    def InitWindow(self):
        self.setWindowTitle(self.title)
        vbox = QVBoxLayout()
        self.frame16 = QFrame(self)
        self.frame16.setGeometry(QRect(1,1,100,25))
        self.saveBtn = QPushButton("Save Script",self.frame16)
        vbox.addWidget(self.saveBtn)
        self.saveBtn.clicked.connect(self.openFileSaveDialog)
        self.saveBtn.setEnabled(1)

        self.frame1 = QFrame(self);self.frame1.setGeometry(QRect(1, 25, 600, 300))
        self.plainText = QPlainTextEdit(self.frame1)
        self.plainText.resize(600,300)
        self.plainText.setPlaceholderText("This is some text for our plaintextedit                                               ")
 
        self.plainText.setReadOnly(True)
 
        self.text = "import matplotlib.pylab as plt"##\nfrom matplotlib.figure import Figure"
        self.text = self.text + "\nfig, ax = plt.subplots()"
        
        if mainWin.ext == '.csv':
            self.text = self.text + "\nimport pandas as pd"
            self.text = self.text + "\nreader = pd.read_csv('"+mainWin.fileName+"')"
        if (mainWin.ext == '.dat')|(mainWin.ext == '.txt'):
            self.text = self.text + "\nimport numpy as np"
            self.text = self.text + "\nreader = np.loadtxt('"+mainWin.fileName+"')"
            self.text = self.text + "\nls=[]"
            self.text = self.text + "\nlbs=[]"
        if mainWin.YIndex != None:
            for i in range(0,len(mainWin.YIndex)):
                if mainWin.ext == '.csv':
                    self.text = self.text + "\nax.plot(reader.iloc[:," \
                        +str(mainWin.XIndex)+"],reader.iloc[:,"+str(mainWin.YIndex[i]) +"],"
                elif (mainWin.ext == '.dat') | (mainWin.ext == '.txt') :
                    self.text = self.text + "\nl=ax.plot(reader[:,"+str(mainWin.XIndex)+"],reader[:,"+str(mainWin.YIndex[i]) +"],"
                self.text = self.text + "\n\tcolor = '" +str(mainWin.colPlotObject[mainWin.YIndex[i]].lineColor)
                self.text = self.text +"', linestyle ='" + str(mainWin.colPlotObject[mainWin.YIndex[i]].lineStyle)
                self.text = self.text +"', linewidth = " + str(mainWin.colPlotObject[mainWin.YIndex[i]].width)
                self.text = self.text +",\n\tdrawstyle = '" + str(mainWin.colPlotObject[mainWin.YIndex[i]].drawStyle)
                self.text = self.text +"', label = '" + str(mainWin.colPlotObject[mainWin.YIndex[i]].label)
                self.text = self.text +"',\n\t marker = '" + str(mainWin.colPlotObject[mainWin.YIndex[i]]. marker)
                self.text = self.text +"', markersize = " + str(mainWin.colPlotObject[mainWin.YIndex[i]].size)
                self.text = self.text +",markeredgecolor = '" + str(mainWin.colPlotObject[mainWin.YIndex[i]].edgeColor)
                self.text = self.text +"',markerfacecolor = '" + str(mainWin.colPlotObject[mainWin.YIndex[i]].faceColor)+"')"
                self.text = self.text +"\nls.append(l[0])"
                self.text = self.text +"\nlbs.append('"+str(mainWin.colPlotObject[mainWin.YIndex[i]].label)+"')"
        if mainWin.YIndexR != None:
            self.text = self.text+ "\nax2 = ax.twinx()"
            for i in range(0,len(mainWin.YIndexR)):
                if mainWin.ext == '.csv':
                    self.text = self.text + "\nax2.plot(reader.iloc[:," \
                        +str(mainWin.XIndex)+"],reader.iloc[:,"+str(mainWin.YIndexR[i]) +"],"
                elif (mainWin.ext == '.dat') | (mainWin.ext == '.txt') :
                    self.text = self.text + "\nl=ax2.plot(reader[:,"+str(mainWin.XIndex)+"],reader[:,"+str(mainWin.YIndexR[i]) +"],"
                self.text = self.text + "\n\tcolor = '" +str(mainWin.colPlotObject[mainWin.YIndexR[i]].lineColor)
                self.text = self.text +"', linestyle ='" + str(mainWin.colPlotObject[mainWin.YIndexR[i]].lineStyle)
                self.text = self.text +"', linewidth = " + str(mainWin.colPlotObject[mainWin.YIndexR[i]].width)
                self.text = self.text +",\n\tdrawstyle = '" + str(mainWin.colPlotObject[mainWin.YIndexR[i]].drawStyle)
                self.text = self.text +"', label = '" + str(mainWin.colPlotObject[mainWin.YIndexR[i]].label)
                self.text = self.text +"',\n\t marker = '" + str(mainWin.colPlotObject[mainWin.YIndexR[i]]. marker)
                self.text = self.text +"', markersize = " + str(mainWin.colPlotObject[mainWin.YIndexR[i]].size)
                self.text = self.text +",markeredgecolor = '" + str(mainWin.colPlotObject[mainWin.YIndexR[i]].edgeColor)
                self.text = self.text +"',markerfacecolor = '" + str(mainWin.colPlotObject[mainWin.YIndexR[i]].faceColor)+"')"
                self.text = self.text +"\nls.append(l[0])"
                self.text = self.text +"\nlbs.append('"+str(mainWin.colPlotObject[mainWin.YIndexR[i]].label)+"')"
            if mainWin.plotAxesProp.yScale2 == 'linear':
                if len(mainWin.YIndexR)!=0:
                    self.text = self.text + "\nax2.ticklabel_format(axis='y', style='sci', scilimits=(" \
                    +str(mainWin.plotAxesProp.ySNL2)+","
                    self.text = self.text + str(mainWin.plotAxesProp.ySNU2)+"), useMathText="+str(mainWin.plotAxesProp.ySNM2)+")"
                else:
                    self.text = self.text + "\nax2.ticklabel_format(axis='y', style='sci', scilimits=(" \
                    +str(mainWin.plotAxesProp.ySNL)+","
                    self.text = self.text + str(mainWin.plotAxesProp.ySNU)+"), useMathText="+str(mainWin.plotAxesProp.ySNM)+")"
            if mainWin.plotAxesProp.yScale2 != None and len(mainWin.YIndexR)!=0:
                self.text = self.text + "\nax2.set_yscale('" + str(mainWin.plotAxesProp.yScale2) + "')"
            else:
                self.text = self.text + "\nax2.set_yscale('" + str(mainWin.plotAxesProp.yScale) + "')"
            if mainWin.plotAxesProp.yLabel2 != None:
                self.text = self.text + "\nax2.set_ylabel('" + str(mainWin.plotAxesProp.yLabel2) + "')"

            if (mainWin.plotAxesProp.yMin2 != None) & (mainWin.plotAxesProp.yMax2 != None) and len(mainWin.YIndexR)!=0:
                self.text = self.text + "\nax2.set_ylim(bottom = " + str(mainWin.plotAxesProp.yMin2)
                self.text = self.text + ", top = " + str(mainWin.plotAxesProp.yMax2)+")"
            else:
                self.text = self.text + "\nax2.set_ylim(bottom = " + str(mainWin.plotAxesProp.yMin)
                self.text = self.text + ", top = " + str(mainWin.plotAxesProp.yMax)+")"
        if mainWin.plotAxesProp.xScale != None:
            self.text = self.text + "\nax.set_xscale('" + str(mainWin.plotAxesProp.xScale) + "')"
        if mainWin.plotAxesProp.yScale != None and len(mainWin.YIndex)!=0:
            self.text = self.text + "\nax.set_yscale('" + str(mainWin.plotAxesProp.yScale) + "')"
        else:
            self.text = self.text + "\nax.set_yscale('" + str(mainWin.plotAxesProp.yScale2) + "')"
        if mainWin.plotAxesProp.xLabel != None:
            self.text = self.text + "\nax.set_xlabel('" + str(mainWin.plotAxesProp.xLabel) + "')"
        if mainWin.plotAxesProp.yLabel != None:
            self.text = self.text + "\nax.set_ylabel('" + str(mainWin.plotAxesProp.yLabel) + "')"
            
        if (mainWin.plotAxesProp.xMin != None) & (mainWin.plotAxesProp.xMax != None):
            self.text = self.text + "\nax.set_xlim(left = " + str(mainWin.plotAxesProp.xMin)
            self.text = self.text + ", right = " + str(mainWin.plotAxesProp.xMax) + ")"

        if (mainWin.plotAxesProp.yMin != None) & (mainWin.plotAxesProp.yMax != None) and len(mainWin.YIndex)!=0:
            self.text = self.text + "\nax.set_ylim(bottom = " + str(mainWin.plotAxesProp.yMin)
            self.text = self.text + ", top = " + str(mainWin.plotAxesProp.yMax)+")" 
        else:
            self.text = self.text + "\nax.set_ylim(bottom = " + str(mainWin.plotAxesProp.yMin2)
            self.text = self.text + ", top = " + str(mainWin.plotAxesProp.yMax2)+")" 
        if mainWin.gridObject.gridEnable:
            self.text =  self.text + "\nax.grid(color ='" + str(mainWin.gridObject.lineColor)
            self.text = self.text + "', linestyle = '" + str(mainWin.gridObject.lineStyle)
            self.text = self.text +"',axis = '" + str(mainWin.gridObject.axis)
            self.text = self.text +"',which = '" + str(mainWin.gridObject.which)
            self.text = self.text + "', linewidth = "+ str(mainWin.gridObject.width) + ")"   
                
        if mainWin.legendEnable:
            if len(mainWin.YIndex)>0 and len(mainWin.YIndexR)==0:
                self.text = self.text + "\nax.legend(loc = '" + mainWin.legProp.location
                self.text = self.text + "',frameon = " + str(mainWin.legProp.frameon)
                self.text = self.text + ", fontsize = " + str(mainWin.legProp.fontsize)
                if mainWin.legProp.title == None:
                    self.text = self.text + ", title = " + 'None'
                else:
                    self.text = self.text + ", title = '" + mainWin.legProp.title + "'"
                self.text = self.text + ",\n\tmarkerfirst = " + str(mainWin.legProp.markerfirst)
                self.text = self.text + ", markerscale = " + str(mainWin.legProp.markerscale)
                self.text = self.text + ", labelspacing = " + str(mainWin.legProp.labelspacing)
                self.text = self.text + ", columnspacing = "+ str(mainWin.legProp.columnspacing)+")"
            elif len(mainWin.YIndexR)>0 and len(mainWin.YIndex)==0:
                self.text = self.text + "\nax2.legend(loc = '" + mainWin.legProp.location
                self.text = self.text + "',frameon = " + str(mainWin.legProp.frameon)
                self.text = self.text + ", fontsize = " + str(mainWin.legProp.fontsize)
                if mainWin.legProp.title == None:
                    self.text = self.text + ", title = " + 'None'
                else:
                    self.text = self.text + ", title = '" + mainWin.legProp.title + "'"
                self.text = self.text + ",\n\tmarkerfirst = " + str(mainWin.legProp.markerfirst)
                self.text = self.text + ", markerscale = " + str(mainWin.legProp.markerscale)
                self.text = self.text + ", labelspacing = " + str(mainWin.legProp.labelspacing)
                self.text = self.text + ", columnspacing = "+ str(mainWin.legProp.columnspacing)+")"
            elif len(mainWin.YIndex) > 0 and len(mainWin.YIndexR) > 0:
                
                self.text = self.text+ "\nax.legend(ls,lbs"
                if mainWin.legProp.title == None:
                    self.text = self.text + ", title = " + 'None'
                else:
                    self.text = self.text + ", title = '" + mainWin.legProp.title + "'"
                self.text = self.text + ",\n\tmarkerfirst = " + str(mainWin.legProp.markerfirst)
                self.text = self.text + ", markerscale = " + str(mainWin.legProp.markerscale)
                self.text = self.text + ", labelspacing = " + str(mainWin.legProp.labelspacing)
                self.text = self.text + ", columnspacing = "+ str(mainWin.legProp.columnspacing)+")"
        if mainWin.titleObject.titleEnable:
            self.text = self.text + "\nax.set_title(label = '" + str(mainWin.titleObject.label)
            self.text = self.text + "',loc = '" +str(mainWin.titleObject.loc) + "')" 
        if mainWin.plotAxesProp.xScale == 'linear':
            self.text = self.text + "\nax.ticklabel_format(axis='x', style='sci', scilimits=(" + str(mainWin.plotAxesProp.xSNL)+","
            self.text = self.text + str(mainWin.plotAxesProp.xSNU)+"), useMathText="+str(mainWin.plotAxesProp.xSNM)+")"
        if mainWin.plotAxesProp.yScale == 'linear' and len(mainWin.YIndex)!=0: 
            self.text = self.text + "\nax.ticklabel_format(axis='y', style='sci', scilimits=("+str(mainWin.plotAxesProp.ySNL)+","
            self.text = self.text + str(mainWin.plotAxesProp.ySNU)+"), useMathText="+str(mainWin.plotAxesProp.ySNM)+")"
        else:
            self.text = self.text + "\nax.ticklabel_format(axis='y', style='sci', scilimits=("+str(mainWin.plotAxesProp.ySNL2)+","
            self.text = self.text + str(mainWin.plotAxesProp.ySNU2)+"), useMathText="+str(mainWin.plotAxesProp.ySNM2)+")"
        self.text = self.text + "\nplt.show()"
        self.text = self.text +"\n"
        self.plainText.appendPlainText(self.text)
        self.plainText.setUndoRedoEnabled(False)
        vbox.addWidget(self.plainText)
        
    def openFileSaveDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()", 
                            "","*.py ", options=options)
        print(self.fileName)
        f= open(os.path.expanduser(self.fileName+'.py'),"w+")
        f.write(self.text)
        f.close();

class PlotObject(object):
    def __init__(self, label = '',lineStyle = 'solid', drawStyle = 'default', width = 0.7, lineColor = None,
                        marker = '', size = 5,edgeColor = 'red',faceColor = 'red'):
        self.lineStyle =lineStyle; self.drawStyle = drawStyle; 
        self.width = width; self.lineColor = lineColor;
        self.marker = marker; self.size = size;
        self.edgeColor = edgeColor;self.faceColor = faceColor;
        self.label = label;
    def setLineStyle(self,lineStyle):
        self.lineStyle = lineStyle;
    def setDrawStyle(self,drawStyle):
        self.drawStyle = drawStyle;
    def setWidth(self,width):
        self.width =  width;
    def setLineColor(self,lineColor):
        self.lineColor = lineColor;
    def setMarker(self,marker):
        self.marker = marker;
    def setSize(self,size):
        self.size = size;
    def setEdgeColor(self,edgeColor):
        self.edgeColor = edgeColor
    def setFaceColor(self,faceColor):
        self.faceColor = faceColor
    def setLabel(self,label):
        self.label = label;
        
class TitleObject(object):
    def __init__(self, label = "", loc = 'center',titleEnable=True):
        self.label = label; self.loc = loc; 
        self.titleEnable = titleEnable
    def setLabel(self,label):
        self.label = label;
    def setTitleEnable(self,titleEnable):
        self.titleEnable = titleEnable
    def setLoc(self,loc):
        self.loc = loc;

class GridObject(object):
    def __init__(self, lineStyle = 'solid', width = 0.7, lineColor = 'silver',
                        which = 'both', axis = 'both', gridEnable = True):
        self.lineStyle =lineStyle;  
        self.width = width; self.lineColor = lineColor;
        self.which = which; self.axis = axis;
        self.gridEnable = True;
    def setLineStyle(self,lineStyle):
        self.lineStyle = lineStyle;
    def setWidth(self,width):
        self.width =  width;
    def setLineColor(self,lineColor):
        self.lineColor = lineColor;
    def setWhich(self,which):
        self.which = which;
    def setAxis(self,axis):
        self.axis = axis
    def setGridEnable(self,gridEnable):
        self.gridEnable = gridEnable;

class xTicksPropObject(object):
    def __init__(self,xTicksEnable = True, direction='inout',rotation=0,xticks = [],xtickslabels=[]):
        self.direction = direction;self.rotation = rotation; 
        self.xticks = xticks; self.xtickslabels = xtickslabels;
        self.xTicksEnable = xTicksEnable
    def setxTicksEnable(self,xTicksEnable):
        self.xTicksEnable = xTicksEnable;
    def setDirection(self,direction):
        self.direction = direction;
    def setRotation(self,rotation):
        self.rotation = rotation;
    def setXTicks(self,xticks):
        self.xticks = xticks;
    def setXTicksLabels(self,xtickslabels):
        self.xtickslabels = xtickslabels;

class AxesPropObject(object):
    def __init__(self, xScale = 'linear',xLabel = 'X-Axis', xMin = None, xMax = None, xSNL=-3,xSNU=3,xSNM=True,
                       yScale = 'linear',yLabel = 'Y-Axis', yMin = None, yMax = None,ySNL=-3,ySNU=3,ySNM=True,
                       yScale2 = 'linear',yLabel2 = 'Y-Axis', yMin2 = None, yMax2 = None,ySNL2=-3,ySNU2=3,ySNM2=True):
        self.xScale = xScale;self.xLabel = xLabel; 
        self.xMin = xMin; self.xMax = xMax;
        self.yScale = yScale; self.yLabel = yLabel;
        self.yMin = yMin; self.yMax = yMax
        self.yScale2 = yScale2; self.yLabel2 = yLabel2;
        self.yMin2 = yMin2; self.yMax2 = yMax2
        self.ySNL = ySNL; self.xSNL = xSNL
        self.ySNU = ySNU; self.xSNU = xSNU
        self.ySNM = ySNM; self.xSNM = xSNM
        self.ySNL2 = ySNL2; 
        self.ySNU2 = ySNU2;
        self.ySNM2 = ySNM2;
    def setxScale(self,xScale):
        self.xScale = xScale;
    def setxLabel(self,xLabel):
        self.xLabel = xLabel;
    def setxMin(self,xMin):
        self.xMin = xMin;
    def setxMax(self,xMax):
        self.xMax = xMax;
    def setyScale(self,yScale):
        self.yScale = yScale;
    def setyLabel(self,yLabel):
        self.yLabel = yLabel;
    def setyMin(self,yMin):
        self.yMin = yMin;
    def setyMax(self,yMax):
        self.yMax = yMax;
    def setyScale2(self,yScale2):
        self.yScale2 = yScale2;
    def setyLabel2(self,yLabel2):
        self.yLabel2 = yLabel2;
    def setyMin2(self,yMin2):
        self.yMin2 = yMin2;
    def setyMax2(self,yMax2):
        self.yMax2 = yMax2;
    def setySNM(self,ySNM):
        self.ySNM = ySNM;
    def setySNM2(self,ySNM2):
        self.ySNM2 = ySNM2;
    def setxSNM(self,xSNM):
        self.xSNM2 = xSNM;
    def setySNL(self,ySNL):
        self.ySNL = ySNL;
    def setySNL2(self,ySNL2):
        self.ySNL2 = ySNL2;
    def setxSNL(self,xSNL):
        self.xSNL = xSNL;
    def setySNU(self,ySNU):
        self.ySNU = ySNU;
    def setySNU2(self,ySNU2):
        self.ySNU2 = ySNU2;
    def setxSNU(self,xSNU):
        self.xSNU = xSNU;
        
class LegendObject(object):
    def __init__(self, location = 'best',frameon = True, fontsize = 10, title = None,
                        markerfirst = True, markerscale = 1.0, labelspacing = 0.5, columnspacing = 2.0):
        self.location = location; self.frameon = frameon;
        self.fontsize = fontsize; self.title = title;
        self.markerfirst = markerfirst; self.markerscale = markerscale; 
        self.labelspacing = labelspacing; self.columnspacing = columnspacing
    def setlocation(self,location):
        self.location = location;
    def setframeon(self,frameon):
        self.frameon = frameon;
    def setfontsize(self,fontsize):
        self.fontsize =  fontsize;
    def settitle(self,title):
        self.title = title;
    def setmarkerfirst(self,markerfirst):
        self.markerfirst = markerfirst;
    def setmarkerscale(self,markerscale):
        self.markerscale = markerscale;
    def setlabelspacing(self,labelspacing):
        self.labelspacing = labelspacing
    def setcolumnspacing(self,columnspacing):
        self.columnspacing = columnspacing
        
class LinePropPopup(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        canvas = QPixmap(185, 130)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 40, 185, 130))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)
        self.frame19 = QFrame(self);self.frame19.setGeometry(QRect(210, 40, 180, 130))
        self.labelC2 = QLabel(self.frame19)
        self.labelC2.setPixmap(canvas)

        self.frame1 = QFrame(self);self.frame1.setGeometry(QRect(10, 10, 250, 25))
        self.combo1 = QComboBox(self.frame1)
        self.combo1.currentIndexChanged.connect(self.yLineProp)

        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 175, 50, 25))
        self.def1 = QLabel("Default:",self.frameD);self.def1.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(65, 175, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.cb1.stateChanged.connect(self.dafaultLineProp)

        self.frame29 = QFrame(self);self.frame29.setGeometry(QRect(220,1,175,25))
        self.leglabelT = QLabel("Label:(used for legend)",self.frame29);
        self.frame28 = QFrame(self);self.frame28.setGeometry(QRect(220,15,175,25))
        self.leglabel = QLineEdit("0.5",self.frame28);
        self.leglabel.setFixedWidth(120)

        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(20, 45, 250, 25))
        self.Label1 = QLabel("Line Properties:",self.frame2);
        myFont = QFont(); myFont.setBold(True);self.Label1.setFont(myFont)

        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 65, 250, 25))
        self.Label2 = QLabel("Line Style:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(90, 65, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("solid");self.combo2.addItem("None");self.combo2.addItem("dotted")
        self.combo2.addItem("dashed");self.combo2.addItem("dashdot");

        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 90, 250, 25))
        self.Label3 = QLabel("Draw Style:",self.frame4)
        self.frame8 = QFrame(self);self.frame8.setGeometry(QRect(90, 90, 250, 25))
        self.combo3 = QComboBox(self.frame8)
        self.combo3.addItem("default");self.combo3.addItem("steps-post");
        self.combo3.addItem("steps-pre");self.combo3.addItem("steps-mid")

        self.frame5 = QFrame(self);self.frame5.setGeometry(QRect(20, 115, 250, 25))
        self.Label4 = QLabel("Width:",self.frame5)
        self.frame26 = QFrame(self);self.frame26.setGeometry(QRect(90,115,115,25))
        self.wdBtn = QLineEdit("0.5",self.frame26);
        self.wdBtn.setFixedWidth(95)
        validator = QDoubleValidator()
        self.wdBtn.setValidator(validator)

        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(20, 140, 250, 25))
        self.Label5 = QLabel("Color:",self.frame6)
        self.frame25 = QFrame(self);self.frame25.setGeometry(QRect(100,140,30,25))
        self.lnBtn = QPushButton(self.frame25);self.lnBtn.clicked.connect(self.openlnColorDlg)
        self.pixmap = QPixmap(10,10);self.pixmap.fill(QColor("red"));
        self.lnCIcon= QIcon(self.pixmap);self.lnBtn.setIcon(self.lnCIcon);

        self.frame11 = QFrame(self);self.frame11.setGeometry(QRect(220, 45, 250, 25))
        self.Label6 = QLabel("Marker Properties:",self.frame11);self.Label6.setFont(myFont)

        self.frame12 = QFrame(self);self.frame12.setGeometry(QRect(220, 65, 250, 25))
        self.Label7 = QLabel("Style:",self.frame12)
        self.frame13 = QFrame(self);self.frame13.setGeometry(QRect(300, 65, 250, 25))
        self.combo4 = QComboBox(self.frame13)
        mStyles = ["",".",",","o","v","^","<",">","1","2","3","4","8","s","p","P","*",
            "h","H","+","x","X","D","d","|","_","0","1","2","3","4","5","6","7","8","9","10","11"]
        for i in range(0,len(mStyles)):
             self.combo4.addItem(mStyles[i]);

             self.frame14 = QFrame(self);self.frame14.setGeometry(QRect(220, 90, 250, 25))
             self.Label8 = QLabel("Size:",self.frame14)
             self.frame27 = QFrame(self);self.frame27.setGeometry(QRect(300,90,115,25))
             self.sizeBtn = QLineEdit("0.5",self.frame27);
             self.sizeBtn.setFixedWidth(75)
             validator = QDoubleValidator()
             self.sizeBtn.setValidator(validator)

             self.frame16 = QFrame(self);self.frame16.setGeometry(QRect(220, 115, 250, 25))
             self.Label9 = QLabel("Edge Color:",self.frame16)
             self.frame24 = QFrame(self);self.frame24.setGeometry(QRect(325,115,30,25))
             self.edBtn = QPushButton(self.frame24);self.edBtn.clicked.connect(self.openEdColorDlg)
             self.pixmap = QPixmap(10,10);self.pixmap.fill(QColor("black"));
             self.edCIcon= QIcon(self.pixmap);self.edBtn.setIcon(self.edCIcon);

             self.frame17 = QFrame(self);self.frame17.setGeometry(QRect(220, 140, 250, 25))
             self.Label10 = QLabel("Face Color:",self.frame17)
             self.frame23 = QFrame(self);self.frame23.setGeometry(QRect(325,140,30,25))
             self.fcBtn = QPushButton(self.frame23);self.fcBtn.clicked.connect(self.openFcColorDlg)
             self.pixmap.fill(QColor("red"));
             self.redIcon= QIcon(self.pixmap);self.fcBtn.setIcon(self.redIcon);

             painter = QPainter(self.labelC.pixmap())
             painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
             painter.drawLine(QPoint(5, 10), QPoint(5, 125))
             painter.drawLine(QPoint(5, 125), QPoint(180, 125))
             painter.drawLine(QPoint(180, 10), QPoint(180, 125))
             painter.drawLine(QPoint(120, 10), QPoint(180, 10))
             painter.drawLine(QPoint(5, 10), QPoint(7, 10))             
             painter.end()

             painter = QPainter(self.labelC2.pixmap())
             painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
             painter.drawLine(QPoint(5, 10), QPoint(5, 125))
             painter.drawLine(QPoint(5, 125), QPoint(175, 125))
             painter.drawLine(QPoint(175, 10), QPoint(175, 125))
             painter.drawLine(QPoint(138, 10), QPoint(175, 10))
             painter.drawLine(QPoint(5, 10), QPoint(7, 10))             
             painter.end()

             self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(120, 170, 90, 25))
             self.applyBtn = QPushButton("Apply",self.frame20);
             self.applyBtn.clicked.connect(self.applyBtnAction)

             self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(210, 170, 90, 25))
             self.cancelBtn = QPushButton("Cancel",self.frame21)
             self.cancelBtn.clicked.connect(self.cancelBtnAction)

             self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(300, 170, 90, 25))
             self.okBtn = QPushButton("Ok",self.frame22)
             self.okBtn.clicked.connect(self.okBtnAction)

             self.linecolor = QColor('red');
             self.edgecolor = QColor('black');
             self.facecolor = QColor('red');
    def applyBtnAction(self):
        if self.combo1.currentText()!= '':
            linePropSelect = self.combo1.currentText()
            selIndexL = mainWin.YCols.indexFromItem(mainWin.YCols.findItems(linePropSelect,Qt.MatchContains)[0])
            selIndex = int(selIndexL.row())
            mainWin.colPlotObject[selIndex].setLineStyle(self.combo2.currentText())
            mainWin.colPlotObject[selIndex].setDrawStyle(self.combo3.currentText())
            mainWin.colPlotObject[selIndex].setWidth(float(self.wdBtn.text()))
            mainWin.colPlotObject[selIndex].setMarker(self.combo4.currentText())
            mainWin.colPlotObject[selIndex].setLineColor(self.linecolor.name())
            mainWin.colPlotObject[selIndex].setSize(float(self.sizeBtn.text()))
            mainWin.colPlotObject[selIndex].setFaceColor(self.facecolor.name())
            mainWin.colPlotObject[selIndex].setEdgeColor(self.edgecolor.name())
            mainWin.colPlotObject[selIndex].setLabel(self.leglabel.text())
            mainWin.plotDataWithChangedOptions()
    def yLineProp(self):
        linePropSelect = self.combo1.currentText()    
        selIndexL = mainWin.YCols.indexFromItem(mainWin.YCols.findItems(linePropSelect,Qt.MatchContains)[0])
        selIndex = int(selIndexL.row())
        self.combo2.setCurrentText(mainWin.colPlotObject[selIndex].lineStyle)
        self.combo3.setCurrentText(mainWin.colPlotObject[selIndex].drawStyle)
        self.wdBtn.setText(str(mainWin.colPlotObject[selIndex].width))
        self.sizeBtn.setText(str(mainWin.colPlotObject[selIndex].size))
        self.pixmap.fill(QColor(mainWin.colPlotObject[selIndex].lineColor));
        self.redIcon= QIcon(self.pixmap);
        self.lnBtn.setIcon(self.redIcon);  
        self.combo4.setCurrentText(mainWin.colPlotObject[selIndex].marker)
        self.pixmap.fill(QColor(mainWin.colPlotObject[selIndex].faceColor));
        self.redIcon= QIcon(self.pixmap);
        self.fcBtn.setIcon(self.redIcon);
        self.pixmap.fill(QColor(mainWin.colPlotObject[selIndex].edgeColor));
        self.redIcon= QIcon(self.pixmap);
        self.edBtn.setIcon(self.redIcon);  
        self.leglabel.setText(str(mainWin.colPlotObject[selIndex].label))
        self.linecolor = QColor(mainWin.colPlotObject[selIndex].lineColor);
        self.edgecolor = QColor(mainWin.colPlotObject[selIndex].edgeColor);
        self.facecolor = QColor(mainWin.colPlotObject[selIndex].faceColor);
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
    def dafaultLineProp(self):
        linePropSelect = self.combo1.currentText()    
        selIndexL = mainWin.YCols.indexFromItem(mainWin.YCols.findItems(linePropSelect,Qt.MatchContains)[0])
        selIndex = int(selIndexL.row())
        mainWin.colPlotObject[selIndex] = PlotObject(label = linePropSelect, lineColor = mainWin.colorSet[selIndex])
        self.yLineProp()

    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();
    def openFcColorDlg(self):
        self.facecolor = QColorDialog.getColor()
        self.pixmap.fill(self.facecolor);
        self.redIcon= QIcon(self.pixmap);
        self.fcBtn.setIcon(self.redIcon);
    def openEdColorDlg(self):
        self.edgecolor = QColorDialog.getColor()
        self.pixmap.fill(self.edgecolor);
        self.edCIcon= QIcon(self.pixmap);
        self.edBtn.setIcon(self.edCIcon); 
    def openlnColorDlg(self):
        self.linecolor = QColorDialog.getColor()
        self.pixmap.fill(self.linecolor);
        self.lnCIcon= QIcon(self.pixmap);
        self.lnBtn.setIcon(self.lnCIcon);

class Header(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 25, 150, 25))
        self.def1 = QLabel("No of lines in Header:",self.frameD);
        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(200, 25, 75, 25))
        self.hBtn = QLineEdit("0",self.frame2);
        self.hBtn.setFixedWidth(75)
        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(200, 70, 90, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)

    def okBtnAction(self):
        self.close();
        mainWin.header = int(self.hBtn.text())
        mainWin.readFile();
        
class GridPopup(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        canvas = QPixmap(380, 130)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 20, 380, 130))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)

        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 175, 50, 25))
        self.def1 = QLabel("Default:",self.frameD);self.def1.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(65, 175, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.cb1.stateChanged.connect(self.dafaultLineProp)

        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(20, 25, 250, 25))
        self.Label1 = QLabel("Grid:",self.frame2);
        myFont = QFont(); myFont.setBold(True);self.Label1.setFont(myFont)
        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 50, 250, 25))
        self.Label3 = QLabel("Grid on-off:",self.frame4)
        self.frameCB2 = QFrame(self);self.frameCB2.setGeometry(QRect(90, 50, 25, 25))
        self.cb2 = QCheckBox(self.frameCB2)
        self.cb2.setCheckState(QtCore.Qt.Checked)

        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 85, 250, 25))
        self.Label2 = QLabel("Line Style:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(90, 85, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("solid");self.combo2.addItem("None");self.combo2.addItem("dotted")
        self.combo2.addItem("dashed");self.combo2.addItem("dashdot");

        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(20, 120, 250, 25))
        self.Label5 = QLabel("Color:",self.frame6)
        self.frame25 = QFrame(self);self.frame25.setGeometry(QRect(100,120,30,25))
        self.lnBtn = QPushButton(self.frame25);self.lnBtn.clicked.connect(self.openlnColorDlg)
        self.pixmap = QPixmap(10,10);self.pixmap.fill(QColor("red"));
        self.lnCIcon= QIcon(self.pixmap);self.lnBtn.setIcon(self.lnCIcon);

        self.frame12 = QFrame(self);self.frame12.setGeometry(QRect(220, 50, 250, 25))
        self.Label7 = QLabel("Which:",self.frame12)
        self.frame13 = QFrame(self);self.frame13.setGeometry(QRect(300, 50, 250, 25))
        self.combo4 = QComboBox(self.frame13)
        self.combo4.addItem("major");self.combo4.addItem("minor");self.combo4.addItem("both");

        self.frame14 = QFrame(self);self.frame14.setGeometry(QRect(220, 85, 250, 25))
        self.Label8 = QLabel("Width:",self.frame14)
        self.frame27 = QFrame(self);self.frame27.setGeometry(QRect(300,85,115,25))
        self.sizeBtn = QLineEdit("0.5",self.frame27);
        self.sizeBtn.setFixedWidth(75)
        validator = QDoubleValidator()
        self.sizeBtn.setValidator(validator)

        self.frame16 = QFrame(self);self.frame16.setGeometry(QRect(220, 120, 250, 25))
        self.Label9 = QLabel("Axis:",self.frame16)
        self.frame24 = QFrame(self);self.frame24.setGeometry(QRect(300,120,150,25))
        self.combo3 = QComboBox(self.frame24)
        self.combo3.addItem("both");self.combo3.addItem("x");self.combo3.addItem("y");

        painter = QPainter(self.labelC.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 10), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(375, 125))
        painter.drawLine(QPoint(375, 10), QPoint(375, 125))
        painter.drawLine(QPoint(50, 10), QPoint(375, 10))
        painter.drawLine(QPoint(5, 10), QPoint(7, 10))
        painter.end()

        self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(120, 170, 90, 25))
        self.applyBtn = QPushButton("Apply",self.frame20);
        self.applyBtn.clicked.connect(self.applyBtnAction)

        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(210, 170, 90, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)

        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(300, 170, 90, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)
        self.linecolor = QColor('black');

    def applyBtnAction(self):
        mainWin.gridObject.setLineStyle(self.combo2.currentText())
        mainWin.gridObject.setAxis(self.combo3.currentText())
        mainWin.gridObject.setWidth(float(self.sizeBtn.text()))
        mainWin.gridObject.setWhich(self.combo4.currentText())
        mainWin.gridObject.setLineColor(self.linecolor.name())
        mainWin.gridObject.gridEnable = self.cb2.isChecked()
        mainWin.m.changeGridProps()
        mainWin.plotW.changeGridProps()
    def GridPropShow(self):
        self.combo2.setCurrentText(mainWin.gridObject.lineStyle)
        self.combo3.setCurrentText(mainWin.gridObject.axis)
        self.sizeBtn.setText(str(mainWin.gridObject.width))

        self.pixmap.fill(QColor(mainWin.gridObject.lineColor));
        self.redIcon= QIcon(self.pixmap);
        self.lnBtn.setIcon(self.redIcon);  

        self.combo4.setCurrentText(mainWin.gridObject.which)
        self.linecolor = QColor(mainWin.gridObject.lineColor);
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        bs = QtCore.Qt.Unchecked
        if mainWin.gridObject.gridEnable:
            bs = QtCore.Qt.Checked
        self.cb2.setCheckState(bs)
    def dafaultLineProp(self):
        mainWin.gridObject = GridObject()
        self.GridPropShow()

    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();
 
    def openlnColorDlg(self):
        self.linecolor = QColorDialog.getColor()
        self.pixmap.fill(self.linecolor);
        self.lnCIcon= QIcon(self.pixmap);
        self.lnBtn.setIcon(self.lnCIcon);

class AxesPopup(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        canvas = QPixmap(185, 130)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 10, 185, 130))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)
        self.frame19 = QFrame(self);self.frame19.setGeometry(QRect(210, 10, 185, 130))
        self.labelC2 = QLabel(self.frame19)
        self.labelC2.setPixmap(canvas)
        self.frameY2 = QFrame(self);self.frameY2.setGeometry(QRect(410, 10, 185, 130))
        self.labelY2 = QLabel(self.frameY2)
        self.labelY2.setPixmap(canvas)
        
        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 145, 50, 25))
        self.def1 = QLabel("Default:",self.frameD);self.def1.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(65, 145, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.cb1.stateChanged.connect(self.dafaultLineProp)

        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(20, 15, 250, 25))
        self.Label1 = QLabel("X-Axis:",self.frame2);self.Label1.setFont(myFont)
        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 35, 250, 25))
        self.Label2 = QLabel("Scale:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(60, 35, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("linear");self.combo2.addItem("log");self.combo2.addItem("logit")

        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 60, 250, 25))
        self.Label3 = QLabel("Label:",self.frame4)
        self.frame8 = QFrame(self);self.frame8.setGeometry(QRect(60, 60, 250, 25))
        self.labelEdit = QLineEdit("X-Axis",self.frame8);
        self.labelEdit.setFixedWidth(120)

        self.frame5 = QFrame(self);self.frame5.setGeometry(QRect(20, 85, 250, 25))
        self.Label4 = QLabel("Left:",self.frame5)
        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(60, 85, 250, 25))
        self.XLimL = QLineEdit("0",self.frame6);
        self.XLimL.setFixedWidth(120)
        validator = QDoubleValidator()
        self.XLimL.setValidator(validator)

        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(20, 110, 250, 25))
        self.Label5 = QLabel("Right:",self.frame6)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(60, 110, 250, 25))
        self.XLimR = QLineEdit("10",self.frame7);
        self.XLimR.setFixedWidth(120)
        validator = QDoubleValidator()
        self.XLimR.setValidator(validator)

        self.frame11 = QFrame(self);self.frame11.setGeometry(QRect(220, 15, 250, 25))
        self.Label6 = QLabel("Y-Axis:",self.frame11);self.Label6.setFont(myFont)

        self.frame12 = QFrame(self);self.frame12.setGeometry(QRect(220, 35, 250, 25))
        self.Label7 = QLabel("Scale:",self.frame12)
        self.frame13 = QFrame(self);self.frame13.setGeometry(QRect(270, 35, 250, 25))
        self.combo4 = QComboBox(self.frame13)
        self.combo4.addItem("linear");self.combo4.addItem("log");self.combo4.addItem("logit")

        self.frame14 = QFrame(self);self.frame14.setGeometry(QRect(220, 60, 250, 25))
        self.Label8 = QLabel("Label:",self.frame14)
        self.frame15 = QFrame(self);self.frame15.setGeometry(QRect(270, 60, 250, 25))
        self.YlabelEdit = QLineEdit("Y-Axis",self.frame15);
        self.YlabelEdit.setFixedWidth(120)

        self.frame16 = QFrame(self);self.frame16.setGeometry(QRect(220, 85, 250, 25))
        self.Label9 = QLabel("Bottom:",self.frame16)
        self.frame24 = QFrame(self);self.frame24.setGeometry(QRect(270,85,250,25))
        self.YlimB = QLineEdit("0",self.frame24);
        self.YlimB.setFixedWidth(120)
        validator = QDoubleValidator()
        self.YlimB.setValidator(validator)

        self.frame17 = QFrame(self);self.frame17.setGeometry(QRect(220, 110, 250, 25))
        self.Label10 = QLabel("Top:",self.frame17)
        self.frame23 = QFrame(self);self.frame23.setGeometry(QRect(270,110,250,25))  
        self.YlimT = QLineEdit("10",self.frame23);
        self.YlimT.setFixedWidth(120)
        validator = QDoubleValidator()
        self.YlimT.setValidator(validator)

        self.frameY11 = QFrame(self);self.frameY11.setGeometry(QRect(420, 15, 250, 25))
        self.LabelY6 = QLabel("Y-Axis (Right):",self.frameY11);self.LabelY6.setFont(myFont)

        self.frameY12 = QFrame(self);self.frameY12.setGeometry(QRect(420, 35, 250, 25))
        self.LabelY7 = QLabel("Scale:",self.frameY12)
        self.frameY13 = QFrame(self);self.frameY13.setGeometry(QRect(470, 35, 250, 25))
        self.comboY4 = QComboBox(self.frameY13)
        self.comboY4.addItem("linear");self.comboY4.addItem("log");self.comboY4.addItem("logit")

        self.frameY14 = QFrame(self);self.frameY14.setGeometry(QRect(420, 60, 250, 25))
        self.LabelY8 = QLabel("Label:",self.frameY14)
        self.frameY15 = QFrame(self);self.frameY15.setGeometry(QRect(470, 60, 250, 25))
        self.YlabelEdit2 = QLineEdit("Y-Axis",self.frameY15);
        self.YlabelEdit2.setFixedWidth(120)

        self.frameY16 = QFrame(self);self.frameY16.setGeometry(QRect(420, 85, 250, 25))
        self.LabelY9 = QLabel("Bottom:",self.frameY16)
        self.frameY24 = QFrame(self);self.frameY24.setGeometry(QRect(470,85,250,25))
        self.YlimB2 = QLineEdit("0",self.frameY24);
        self.YlimB2.setFixedWidth(120)
        validator = QDoubleValidator()
        self.YlimB2.setValidator(validator)

        self.frameY17 = QFrame(self);self.frameY17.setGeometry(QRect(420, 110, 250, 25))
        self.LabelY10 = QLabel("Top:",self.frameY17)
        self.frameY23 = QFrame(self);self.frameY23.setGeometry(QRect(470,110,250,25))  
        self.YlimT2 = QLineEdit("10",self.frameY23);
        self.YlimT2.setFixedWidth(120)
        validator = QDoubleValidator()
        self.YlimT2.setValidator(validator)
        painter = QPainter(self.labelC.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 10), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(180, 125))
        painter.drawLine(QPoint(180, 10), QPoint(180, 125))
        painter.drawLine(QPoint(60, 10), QPoint(180, 10))
        painter.drawLine(QPoint(5, 10), QPoint(7, 10))
        painter.end()

        painter = QPainter(self.labelC2.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 10), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(180, 125))
        painter.drawLine(QPoint(180, 10), QPoint(180, 125))
        painter.drawLine(QPoint(60, 10), QPoint(180, 10))
        painter.drawLine(QPoint(5, 10), QPoint(7, 10))
        painter.end()

        painter = QPainter(self.labelY2.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 10), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(180, 125))
        painter.drawLine(QPoint(180, 10), QPoint(180, 125))
        painter.drawLine(QPoint(110, 10), QPoint(180, 10))
        painter.drawLine(QPoint(5, 10), QPoint(7, 10))
        painter.end()

        self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(120, 140, 90, 25))
        self.applyBtn = QPushButton("Apply",self.frame20);
        self.applyBtn.clicked.connect(self.applyBtnAction)

        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(210, 140, 90, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)

        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(300, 140, 90, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)

    def applyBtnAction(self):
        if self.combo2.currentText()!='None':
            mainWin.plotAxesProp.setxScale(self.combo2.currentText())
        if self.labelEdit.text() !='None':
            mainWin.plotAxesProp.setxLabel(self.labelEdit.text())
        if self.XLimL.text()!='None':
            mainWin.plotAxesProp.setxMin(float(self.XLimL.text()))
        if self.XLimR.text() !='None':
            mainWin.plotAxesProp.setxMax(float(self.XLimR.text()))
        if self.combo4.currentText()!='None':
            mainWin.plotAxesProp.setyScale(self.combo4.currentText())
        if self.YlabelEdit.text()!='None':
            mainWin.plotAxesProp.setyLabel(self.YlabelEdit.text())
        if self.YlimB.text()!='None':
            mainWin.plotAxesProp.setyMin(float(self.YlimB.text()))
        if self.YlimT.text()!='None':
            mainWin.plotAxesProp.setyMax(float(self.YlimT.text()))
        if self.comboY4.currentText()!='None':
            mainWin.plotAxesProp.setyScale2(self.comboY4.currentText())
        if self.YlabelEdit2.text()!='None':
            mainWin.plotAxesProp.setyLabel2(self.YlabelEdit2.text())
        if self.YlimB2.text()!='None':
            mainWin.plotAxesProp.setyMin2(float(self.YlimB2.text()))
        if self.YlimT2.text()!='None':
            mainWin.plotAxesProp.setyMax2(float(self.YlimT2.text()))
        mainWin.m.changeAxesProps()
        mainWin.plotW.changeAxesProps()

    def AxesPropShow(self):
        self.combo2.setCurrentText(str(mainWin.plotAxesProp.xScale))
        self.combo4.setCurrentText(mainWin.plotAxesProp.yScale)
        self.comboY4.setCurrentText(mainWin.plotAxesProp.yScale2)
        self.labelEdit.setText(str(mainWin.plotAxesProp.xLabel))
        if mainWin.plotAxesProp.xMin!= None:
            self.XLimL.setText(str(format(mainWin.plotAxesProp.xMin)))
        else:
            self.XLimL.setText(str(mainWin.plotAxesProp.xMin))

        if mainWin.plotAxesProp.xMax!= None:           
            self.XLimR.setText(str(format(mainWin.plotAxesProp.xMax)))
        else:
            self.XLimR.setText(str(mainWin.plotAxesProp.xMax))      

        self.YlabelEdit.setText(str(mainWin.plotAxesProp.yLabel))

        if mainWin.plotAxesProp.yMin != None:
            self.YlimB.setText(str(format(mainWin.plotAxesProp.yMin)))
        else:
            self.YlimB.setText(str(mainWin.plotAxesProp.yMin))
        if mainWin.plotAxesProp.yMax != None:
            self.YlimT.setText(str(format(mainWin.plotAxesProp.yMax)))
        else:
            self.YlimT.setText(str(mainWin.plotAxesProp.yMax))
        
        if mainWin.plotAxesProp.yMin2 != None:
            self.YlimB2.setText(str(format(mainWin.plotAxesProp.yMin2)))
        else:
            self.YlimB2.setText(str(mainWin.plotAxesProp.yMin2))
        if mainWin.plotAxesProp.yMax2 != None:
            self.YlimT2.setText(str(format(mainWin.plotAxesProp.yMax2)))
        else:
            self.YlimT2.setText(str(mainWin.plotAxesProp.yMax2))

        self.cb1.setCheckState(QtCore.Qt.Unchecked)
    def dafaultLineProp(self):       
        mainWin.plotAxesProp = AxesPropObject() 
        self.AxesPropShow()

    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();  

class TitlePopup(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        canvas = QPixmap(250, 100)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 10, 250, 100))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)

        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 145, 50, 25))
        self.def1 = QLabel("Default:",self.frameD);self.def1.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(65, 145, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.cb1.stateChanged.connect(self.dafaultLineProp)

        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(20, 15, 250, 25))
        self.Label1 = QLabel("Title:",self.frame2);self.Label1.setFont(myFont)

        self.frameT = QFrame(self);self.frameT.setGeometry(QRect(20, 30, 50, 25))
        self.def1 = QLabel("Title:",self.frameT);
        self.frameTE = QFrame(self);self.frameTE.setGeometry(QRect(80, 30, 25, 25))
        self.cb3 = QCheckBox(self.frameTE)
        self.cb3.setCheckState(QtCore.Qt.Checked)

        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 55, 250, 25))
        self.Label2 = QLabel("Position:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(80, 55, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("left");self.combo2.addItem("center");self.combo2.addItem("right")

        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 80, 250, 25))
        self.Label3 = QLabel("Label:",self.frame4)
        self.frame8 = QFrame(self);self.frame8.setGeometry(QRect(80, 80, 250, 25))
        self.labelEdit = QLineEdit("",self.frame8);
        self.labelEdit.setFixedWidth(175)

        painter = QPainter(self.labelC.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 5), QPoint(5, 95))
        painter.drawLine(QPoint(5, 95), QPoint(245, 95))
        painter.drawLine(QPoint(245, 5), QPoint(245, 95))
        painter.drawLine(QPoint(5, 5), QPoint(245, 5))
        painter.end()

        self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(20, 170, 80, 25))
        self.applyBtn = QPushButton("Apply",self.frame20);
        self.applyBtn.clicked.connect(self.applyBtnAction)

        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(120, 170,80, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)

        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(220, 170, 80, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)

    def applyBtnAction(self):
        mainWin.titleObject.setLoc(self.combo2.currentText())
        mainWin.titleObject.setLabel(self.labelEdit.text())
        mainWin.titleObject.setTitleEnable(self.cb3.isChecked())
        mainWin.m.changeTitle()
        mainWin.plotW.changeTitle()
    def TitlePropShow(self):
        self.combo2.setCurrentText(str(mainWin.titleObject.loc))
        self.labelEdit.setText(str(mainWin.titleObject.label))
        bs = QtCore.Qt.Unchecked
        if mainWin.titleObject.titleEnable:
            bs = QtCore.Qt.Checked
        self.cb3.setCheckState(bs)

    def dafaultLineProp(self):
        mainWin.titleObject = TitleObject()
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.TitlePropShow()

    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();
class TickLabelNotation(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();
    def widget(self):
        canvas = QPixmap(460, 130)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 10, 460, 130))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)
        
        self.frameEg = QFrame(self);self.frameEg.setGeometry(QRect(20, 15, 450, 40))
        self.LabelEg = QLabel("E.g.: m=-3, n=3 ticklabels beyond 10\u207b\u00b3 - 10\u00b3 will be converted to\n       scientific notation",self.frameEg)
        self.LabelEg.setFont(myFont)
        
        self.frameSN = QFrame(self);self.frameSN.setGeometry(QRect(20, 50, 250, 25))
        self.LabelSN = QLabel("Limits For Scientific Notation:",self.frameSN);self.LabelSN.setFont(myFont)
        self.frameSNLb = QFrame(self);self.frameSNLb.setGeometry(QRect(20, 75, 170, 25))
        self.SNLb = QLabel("X-TickLabels: Lower(m):",self.frameSNLb);self.SNLb.setFont(myFont)
        self.frameSNL = QFrame(self);self.frameSNL.setGeometry(QRect(185, 75, 50, 25))
        self.SNLEdit = QLineEdit("-3",self.frameSNL)
        self.SNLEdit.setFixedWidth(40)
        self.frameSNLu = QFrame(self);self.frameSNLu.setGeometry(QRect(230, 75, 150, 25))
        self.SNLu = QLabel("Upper(n):",self.frameSNLu);self.SNLu.setFont(myFont)
        self.frameSNU = QFrame(self);self.frameSNU.setGeometry(QRect(295, 75, 50, 25))
        self.SNUEdit = QLineEdit("3",self.frameSNU)
        self.SNUEdit.setFixedWidth(40)
        self.frameSNXM = QFrame(self);self.frameSNXM.setGeometry(QRect(340, 75, 150, 25))
        self.SNXM = QLabel("UseMathsText:",self.frameSNXM);self.SNXM.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(445, 75, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Checked)
        
        self.frameYSNLb = QFrame(self);self.frameYSNLb.setGeometry(QRect(20, 100, 170, 25))
        self.YSNLb = QLabel("Y-TickLabels: Lower(m):",self.frameYSNLb);self.YSNLb.setFont(myFont)
        self.frameYSNL = QFrame(self);self.frameYSNL.setGeometry(QRect(185, 100, 50, 25))
        self.YSNLEdit = QLineEdit("-3",self.frameYSNL)
        self.YSNLEdit.setFixedWidth(40)
        self.frameYSNLu = QFrame(self);self.frameYSNLu.setGeometry(QRect(230, 100, 150, 25))
        self.YSNLu = QLabel("Upper(n):",self.frameYSNLu);self.YSNLu.setFont(myFont)
        self.frameYSNU = QFrame(self);self.frameYSNU.setGeometry(QRect(295, 100, 50, 25))
        self.YSNUEdit = QLineEdit("3",self.frameYSNU)
        self.YSNUEdit.setFixedWidth(40)
        self.frameSNYM = QFrame(self);self.frameSNYM.setGeometry(QRect(340, 100, 150, 25))
        self.SNYM = QLabel("UseMathsText:",self.frameSNYM);self.SNYM.setFont(myFont)
        self.frameCB2 = QFrame(self);self.frameCB2.setGeometry(QRect(445, 100, 25, 25))
        self.cb2 = QCheckBox(self.frameCB2)
        self.cb2.setCheckState(QtCore.Qt.Checked)
        
        painter = QPainter(self.labelC.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 5), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(455, 125))
        painter.drawLine(QPoint(455, 5), QPoint(455, 125))
        painter.drawLine(QPoint(5, 5), QPoint(455, 5))
        painter.end()
        
        self.frame23 =  QFrame(self);self.frame23.setGeometry(QRect(20, 170, 80, 25))
        self.applyBtn = QPushButton("Apply",self.frame23)
        self.applyBtn.clicked.connect(self.applyBtnAction)
        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(120, 170,80, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)
        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(220, 170, 80, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)

    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();
    def applyBtnAction(self):
        if self.YSNUEdit.text() != 'None':
            mainWin.plotAxesProp.setySNU(float(self.YSNUEdit.text()))
        if self.YSNLEdit.text() != 'None':
            mainWin.plotAxesProp.setySNL(float(self.YSNLEdit.text()))    
        if self.SNUEdit.text() != 'None':
            mainWin.plotAxesProp.setxSNU(float(self.SNUEdit.text()))
        if self.SNLEdit.text() != 'None':
            mainWin.plotAxesProp.setxSNL(float(self.SNLEdit.text()))
        bs = QtCore.Qt.Unchecked
        if self.cb2.isChecked():
            bs = QtCore.Qt.Checked
            mainWin.plotAxesProp.setySNM(True)
        else:
            mainWin.plotAxesProp.setySNM(False)
            
        if self.cb1.isChecked():
            bs = QtCore.Qt.Checked
            mainWin.plotAxesProp.setxSNM(True)
        else:
            mainWin.plotAxesProp.setxSNM(False)
        mainWin.m.changeAxesProps()
        mainWin.plotW.changeAxesProps()
class xTicksPopup(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        canvas = QPixmap(460, 130)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 10, 460, 130))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)

        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 145, 50, 25))
        self.def1 = QLabel("Default:",self.frameD);self.def1.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(65, 145, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.cb1.stateChanged.connect(self.dafaultLineProp)

        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(20, 15, 250, 25))
        self.Label1 = QLabel("x-Ticks:",self.frame2);self.Label1.setFont(myFont)

        self.frameT = QFrame(self);self.frameT.setGeometry(QRect(20, 30, 50, 25))
        self.def1 = QLabel("x-Ticks:",self.frameT);
        self.frameTE = QFrame(self);self.frameTE.setGeometry(QRect(105, 30, 25, 25))
        self.cb3 = QCheckBox(self.frameTE)
        self.cb3.setCheckState(QtCore.Qt.Checked)

        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 55, 250, 25))
        self.Label2 = QLabel("Direction:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(105, 55, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("in");self.combo2.addItem("out");self.combo2.addItem("inout")

        self.frame5 = QFrame(self);self.frame5.setGeometry(QRect(210, 55, 250, 25))
        self.Label4 = QLabel("Rotation:",self.frame5)
        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(270, 55, 50, 25))
        self.rotEdit = QLineEdit("0",self.frame6)
        self.rotEdit.setFixedWidth(50)
        self.frame25 = QFrame(self);self.frame25.setGeometry(QRect(320, 58, 250, 25))
        self.Label7 = QLabel("(in degrees:)",self.frame25)

        self.frame9 = QFrame(self);self.frame9.setGeometry(QRect(20, 80, 280, 25))
        self.Label5 = QLabel("xTicks:",self.frame9)
        self.frame10 = QFrame(self);self.frame10.setGeometry(QRect(105, 80, 280, 25))
        self.ticksEdit = QLineEdit("",self.frame10);
        self.ticksEdit.setFixedWidth(280)

        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 105, 280, 25))
        self.Label3 = QLabel("xTicksLabels:",self.frame4)
        self.frame8 = QFrame(self);self.frame8.setGeometry(QRect(105, 105, 280, 25))
        self.labelEdit = QLineEdit("",self.frame8);
        self.labelEdit.setFixedWidth(280)

        painter = QPainter(self.labelC.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 5), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(455, 125))
        painter.drawLine(QPoint(455, 5), QPoint(455, 125))
        painter.drawLine(QPoint(5, 5), QPoint(455, 5))
        painter.end()

        self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(385, 80, 80, 25))
        self.applyBtn = QPushButton("Apply",self.frame20);
        self.applyBtn.clicked.connect(self.applyBtnAction)

        self.frame23 =  QFrame(self);self.frame23.setGeometry(QRect(385, 105, 80, 25))
        self.applyBtn2 = QPushButton("Apply",self.frame23);
        self.applyBtn2.clicked.connect(self.applyBtnAction2)

        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(120, 170,80, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)

        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(220, 170, 80, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)

    def applyBtnAction(self):    
        mainWin.xTicksObject.setDirection(self.combo2.currentText())
        mainWin.xTicksObject.setRotation(self.rotEdit.text())
        mainWin.xTicksObject.setxTicksEnable(self.cb3.isChecked())
        if self.ticksEdit.text() != 'None':
            xticklist = []
            xtickl = self.ticksEdit.text();
            xtickstr = xtickl[1:len(xtickl)-1]
            if xtickstr.find(",")>=0:
                for item in xtickstr.split(","):
                    xticklist.append(float(item))
            else:
                for item in xtickstr.split():
                    xticklist.append(float(item))

            mainWin.xTicksObject.setXTicks(xticklist)

        if self.labelEdit.text() != 'None':
            xticklabel = self.labelEdit.text();
            xticklabelstr = xticklabel[1:len(xticklabel)-1]
            xticklabelstr = xticklabelstr.replace("'","")
            mainWin.xTicksObject.setXTicksLabels(list(xticklabelstr.split(", ")))
        mainWin.m.changeXTicks()
        mainWin.plotW.changeXTicks()
    def applyBtnAction2(self):
        mainWin.xTicksObject.setDirection(self.combo2.currentText())
        mainWin.xTicksObject.setRotation(self.rotEdit.text())
        mainWin.xTicksObject.setxTicksEnable(self.cb3.isChecked())
        if self.labelEdit.text() != 'None':
            xticklabel = self.labelEdit.text();
            xticklabelstr = xticklabel[1:len(xticklabel)-1]
            xticklabelstr = xticklabelstr.replace("'","")
            xticklabelstr = xticklabelstr.replace("\\\\","\\")
            mainWin.xTicksObject.setXTicksLabels(list(xticklabelstr.split(", ")))

        mainWin.m.changeXTicksLabels()
        mainWin.plotW.changeXTicksLabels()
    def xTicksPropShow(self):
        self.combo2.setCurrentText(str(mainWin.xTicksObject.direction))
        self.rotEdit.setText(str(mainWin.xTicksObject.rotation))
        self.ticksEdit.setText(str(mainWin.xTicksObject.xticks))
        self.labelEdit.setText(str(mainWin.xTicksObject.xtickslabels))
        bs = QtCore.Qt.Unchecked
        if mainWin.xTicksObject.xTicksEnable:
            bs = QtCore.Qt.Checked
        self.cb3.setCheckState(bs)

    def dafaultLineProp(self):
        mainWin.xTicksObject.direction = 'inout'
        mainWin.xTicksObject.rotation = 0
        self.xTicksPropShow()
    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();

class yTicksPopup(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        canvas = QPixmap(460, 130)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 10, 460, 130))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)

        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 145, 50, 25))
        self.def1 = QLabel("Default:",self.frameD);self.def1.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(65, 145, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.cb1.stateChanged.connect(self.dafaultLineProp)

        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(20, 15, 250, 25))
        self.Label1 = QLabel("y-Ticks:",self.frame2);self.Label1.setFont(myFont)

        self.frameT = QFrame(self);self.frameT.setGeometry(QRect(20, 30, 50, 25))
        self.def1 = QLabel("y-Ticks:",self.frameT);
        self.frameTE = QFrame(self);self.frameTE.setGeometry(QRect(105, 30, 25, 25))
        self.cb3 = QCheckBox(self.frameTE)
        self.cb3.setCheckState(QtCore.Qt.Checked)

        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 55, 250, 25))
        self.Label2 = QLabel("Direction:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(105, 55, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("in");self.combo2.addItem("out");self.combo2.addItem("inout")

        self.frame5 = QFrame(self);self.frame5.setGeometry(QRect(210, 55, 250, 25))
        self.Label4 = QLabel("Rotation:",self.frame5)
        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(270, 55, 50, 25))
        self.rotEdit = QLineEdit("0",self.frame6)
        self.rotEdit.setFixedWidth(50)
        self.frame25 = QFrame(self);self.frame25.setGeometry(QRect(320, 58, 250, 25))
        self.Label7 = QLabel("(in degrees:)",self.frame25)

        self.frame9 = QFrame(self);self.frame9.setGeometry(QRect(20, 80, 280, 25))
        self.Label5 = QLabel("yTicks:",self.frame9)
        self.frame10 = QFrame(self);self.frame10.setGeometry(QRect(105, 80, 280, 25))
        self.ticksEdit = QLineEdit("",self.frame10);
        self.ticksEdit.setFixedWidth(280)

        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 105, 280, 25))
        self.Label3 = QLabel("yTicksLabels:",self.frame4)
        self.frame8 = QFrame(self);self.frame8.setGeometry(QRect(105, 105, 280, 25))
        self.labelEdit = QLineEdit("",self.frame8);
        self.labelEdit.setFixedWidth(280)

        painter = QPainter(self.labelC.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 5), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(455, 125))
        painter.drawLine(QPoint(455, 5), QPoint(455, 125))
        painter.drawLine(QPoint(5, 5), QPoint(455, 5))
        painter.end()

        self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(385, 80, 80, 25))
        self.applyBtn = QPushButton("Apply",self.frame20);
        self.applyBtn.clicked.connect(self.applyBtnAction)

        self.frame23 =  QFrame(self);self.frame23.setGeometry(QRect(385, 105, 80, 25))
        self.applyBtn2 = QPushButton("Apply",self.frame23);
        self.applyBtn2.clicked.connect(self.applyBtnAction2)

        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(120, 170,80, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)

        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(220, 170, 80, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)

    def applyBtnAction(self):    
        mainWin.yTicksObject.setDirection(self.combo2.currentText())
        mainWin.yTicksObject.setRotation(self.rotEdit.text())
        mainWin.yTicksObject.setxTicksEnable(self.cb3.isChecked())
        if self.ticksEdit.text() != 'None':
            xticklist = []
            xtickl = self.ticksEdit.text();
            xtickstr = xtickl[1:len(xtickl)-1]
            if xtickstr.find(",")>=0:
                for item in xtickstr.split(","):
                    xticklist.append(float(item))
            else:
                for item in xtickstr.split():
                    xticklist.append(float(item))

            mainWin.yTicksObject.setXTicks(list(xticklist))

        if self.labelEdit.text() != 'None':
            xticklabel = self.labelEdit.text();
            xticklabelstr = xticklabel[1:len(xticklabel)-1]
            xticklabelstr = xticklabelstr.replace("'","")
            mainWin.yTicksObject.setXTicksLabels(list(xticklabelstr.split(", ")))
        mainWin.m.changeYTicks()
        mainWin.plotW.changeYTicks()
    def applyBtnAction2(self):
        mainWin.yTicksObject.setDirection(self.combo2.currentText())
        mainWin.yTicksObject.setRotation(self.rotEdit.text())
        mainWin.yTicksObject.setxTicksEnable(self.cb3.isChecked())
        if self.labelEdit.text() != 'None':
            xticklabel = self.labelEdit.text();
            xticklabelstr = xticklabel[1:len(xticklabel)-1]
            xticklabelstr = xticklabelstr.replace("'","")
            xticklabelstr = xticklabelstr.replace("\\\\","\\")
            mainWin.yTicksObject.setXTicksLabels(list(xticklabelstr.split(", ")))

        mainWin.m.changeYTicksLabels()
        mainWin.plotW.changeYTicksLabels()
    def yTicksPropShow(self):
        self.combo2.setCurrentText(str(mainWin.yTicksObject.direction))
        self.rotEdit.setText(str(mainWin.yTicksObject.rotation))
        self.ticksEdit.setText(str(mainWin.yTicksObject.xticks))
        self.labelEdit.setText(str(mainWin.yTicksObject.xtickslabels))
        bs = QtCore.Qt.Unchecked
        if mainWin.yTicksObject.xTicksEnable:
            bs = QtCore.Qt.Checked
        self.cb3.setCheckState(bs)

    def dafaultLineProp(self):
        mainWin.yTicksObject.direction = 'inout'
        mainWin.yTicksObject.rotation = 0
        self.yTicksPropShow()
    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();
class yTicksPopup2(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.widget();

    def widget(self):
        canvas = QPixmap(460, 130)
        canvas.fill(QColor("#FFFFFF"));
        myFont = QFont(); myFont.setBold(True);
        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 10, 460, 130))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)

        self.frameD = QFrame(self);self.frameD.setGeometry(QRect(10, 145, 50, 25))
        self.def1 = QLabel("Default:",self.frameD);self.def1.setFont(myFont)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(65, 145, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Unchecked)
        self.cb1.stateChanged.connect(self.dafaultLineProp)

        self.frame2 = QFrame(self);self.frame2.setGeometry(QRect(20, 15, 250, 25))
        self.Label1 = QLabel("y-Ticks:",self.frame2);self.Label1.setFont(myFont)

        self.frameT = QFrame(self);self.frameT.setGeometry(QRect(20, 30, 50, 25))
        self.def1 = QLabel("y-Ticks:",self.frameT);
        self.frameTE = QFrame(self);self.frameTE.setGeometry(QRect(105, 30, 25, 25))
        self.cb3 = QCheckBox(self.frameTE)
        self.cb3.setCheckState(QtCore.Qt.Checked)

        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 55, 250, 25))
        self.Label2 = QLabel("Direction:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(105, 55, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("in");self.combo2.addItem("out");self.combo2.addItem("inout")

        self.frame5 = QFrame(self);self.frame5.setGeometry(QRect(210, 55, 250, 25))
        self.Label4 = QLabel("Rotation:",self.frame5)
        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(270, 55, 50, 25))
        self.rotEdit = QLineEdit("0",self.frame6)
        self.rotEdit.setFixedWidth(50)
        self.frame25 = QFrame(self);self.frame25.setGeometry(QRect(320, 58, 250, 25))
        self.Label7 = QLabel("(in degrees:)",self.frame25)

        self.frame9 = QFrame(self);self.frame9.setGeometry(QRect(20, 80, 280, 25))
        self.Label5 = QLabel("yTicks:",self.frame9)
        self.frame10 = QFrame(self);self.frame10.setGeometry(QRect(105, 80, 280, 25))
        self.ticksEdit = QLineEdit("",self.frame10);
        self.ticksEdit.setFixedWidth(280)

        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 105, 280, 25))
        self.Label3 = QLabel("yTicksLabels:",self.frame4)
        self.frame8 = QFrame(self);self.frame8.setGeometry(QRect(105, 105, 280, 25))
        self.labelEdit = QLineEdit("",self.frame8);
        self.labelEdit.setFixedWidth(280)

        painter = QPainter(self.labelC.pixmap())
        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 5), QPoint(5, 125))
        painter.drawLine(QPoint(5, 125), QPoint(455, 125))
        painter.drawLine(QPoint(455, 5), QPoint(455, 125))
        painter.drawLine(QPoint(5, 5), QPoint(455, 5))
        painter.end()

        self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(385, 80, 80, 25))
        self.applyBtn = QPushButton("Apply",self.frame20);
        self.applyBtn.clicked.connect(self.applyBtnAction)

        self.frame23 =  QFrame(self);self.frame23.setGeometry(QRect(385, 105, 80, 25))
        self.applyBtn2 = QPushButton("Apply",self.frame23);
        self.applyBtn2.clicked.connect(self.applyBtnAction2)

        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(120, 170,80, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)

        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(220, 170, 80, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)

    def applyBtnAction(self):    
        mainWin.yTicksObject2.setDirection(self.combo2.currentText())
        mainWin.yTicksObject2.setRotation(self.rotEdit.text())
        mainWin.yTicksObject2.setxTicksEnable(self.cb3.isChecked())
        if self.ticksEdit.text() != 'None':
            xticklist = []
            xtickl = self.ticksEdit.text();
            xtickstr = xtickl[1:len(xtickl)-1]
            if xtickstr.find(",")>=0:
                for item in xtickstr.split(","):
                    xticklist.append(float(item))
            else:
                for item in xtickstr.split():
                    xticklist.append(float(item))

            mainWin.yTicksObject2.setXTicks(list(xticklist))

        if self.labelEdit.text() != 'None':
            xticklabel = self.labelEdit.text();
            xticklabelstr = xticklabel[1:len(xticklabel)-1]
            xticklabelstr = xticklabelstr.replace("'","")
            mainWin.yTicksObject2.setXTicksLabels(list(xticklabelstr.split(", ")))
        mainWin.m.changeYTicks2()
        mainWin.plotW.changeYTicks2()
    def applyBtnAction2(self):
        mainWin.yTicksObject2.setDirection(self.combo2.currentText())
        mainWin.yTicksObject2.setRotation(self.rotEdit.text())
        mainWin.yTicksObject2.setxTicksEnable(self.cb3.isChecked())
        if self.labelEdit.text() != 'None':
            xticklabel = self.labelEdit.text();
            xticklabelstr = xticklabel[1:len(xticklabel)-1]
            xticklabelstr = xticklabelstr.replace("'","")
            xticklabelstr = xticklabelstr.replace("\\\\","\\")
            mainWin.yTicksObject2.setXTicksLabels(list(xticklabelstr.split(", ")))
        mainWin.m.changeYTicksLabels2()
        mainWin.plotW.changeYTicksLabels2()
    def yTicksPropShow(self):
        self.combo2.setCurrentText(str(mainWin.yTicksObject2.direction))
        self.rotEdit.setText(str(mainWin.yTicksObject2.rotation))
        self.ticksEdit.setText(str(mainWin.yTicksObject2.xticks))
        self.labelEdit.setText(str(mainWin.yTicksObject2.xtickslabels))
        bs = QtCore.Qt.Unchecked
        if mainWin.yTicksObject2.xTicksEnable:
            bs = QtCore.Qt.Checked
        self.cb3.setCheckState(bs)

    def dafaultLineProp(self):
        mainWin.yTicksObject2.direction = 'inout'
        mainWin.yTicksObject2.rotation = 0
        self.yTicksPropShow2()
    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();

class LegendPopup(QMainWindow):

    def __init__(self):
        QMainWindow.__init__(self) 
        self.widget();

    def widget(self):
        canvas = QPixmap(380, 125)
        canvas.fill(QColor("#FFFFFF"));

        self.frame18 = QFrame(self);self.frame18.setGeometry(QRect(10, 30, 380, 125))
        self.labelC = QLabel(self.frame18)      
        self.labelC.setPixmap(canvas)

        self.frame1 = QFrame(self);self.frame1.setGeometry(QRect(20, 10, 250, 25))
        self.combo1 = QLabel("Show Legend",self.frame1)
        self.frameCB1 = QFrame(self);self.frameCB1.setGeometry(QRect(150, 10, 25, 25))
        self.cb1 = QCheckBox(self.frameCB1)
        self.cb1.setCheckState(QtCore.Qt.Checked)

        myFont = QFont(); myFont.setBold(True);
        self.frame3 = QFrame(self);self.frame3.setGeometry(QRect(20, 40, 250, 25))
        self.Label2 = QLabel("Location:",self.frame3)
        self.frame7 = QFrame(self);self.frame7.setGeometry(QRect(90, 40, 250, 25))
        self.combo2 = QComboBox(self.frame7)
        self.combo2.addItem("best");self.combo2.addItem("upper right");self.combo2.addItem("upper left")
        self.combo2.addItem("lower left");self.combo2.addItem("lower right");
        self.combo2.addItem("lower center");self.combo2.addItem("upper center");
        self.combo2.addItem("center left");self.combo2.addItem("center right");self.combo2.addItem("center");

        self.frame4 = QFrame(self);self.frame4.setGeometry(QRect(20, 70, 250, 25))
        self.Label3 = QLabel("Font Size:",self.frame4)
        self.frame8 = QFrame(self);self.frame8.setGeometry(QRect(90, 70, 250, 25))
        self.fSize = QLineEdit("10",self.frame8);
        self.fSize.setFixedWidth(75)
        validator = QDoubleValidator()
        self.fSize.setValidator(validator)
        self.frame5 = QFrame(self);self.frame5.setGeometry(QRect(20, 100, 250, 25))
        self.Label4 = QLabel("Title:",self.frame5)
        self.frame26 = QFrame(self);self.frame26.setGeometry(QRect(90,100,115,25))
        self.legTitle = QLineEdit(None,self.frame26);
        self.legTitle.setFixedWidth(95)
        self.frame6 = QFrame(self);self.frame6.setGeometry(QRect(20, 130, 250, 25))
        self.Label5 = QLabel("Frame:",self.frame6)
        self.frame25 = QFrame(self);self.frame25.setGeometry(QRect(100,130,30,25))
        self.fr1 = QCheckBox(self.frame25)
        self.fr1.setCheckState(QtCore.Qt.Checked)

        self.frame12 = QFrame(self);self.frame12.setGeometry(QRect(205, 40, 250, 25))
        self.Label7 = QLabel("Label Spacing:",self.frame12)
        self.frame13 = QFrame(self);self.frame13.setGeometry(QRect(308, 40, 250, 25))
        self.lblSpc = QLineEdit("0.5",self.frame13);
        self.lblSpc.setFixedWidth(75)
        validator = QDoubleValidator()
        self.lblSpc.setValidator(validator)
        self.frame14 = QFrame(self);self.frame14.setGeometry(QRect(205, 70, 250, 25))
        self.Label8 = QLabel("Marker Scale:",self.frame14)
        self.frame27 = QFrame(self);self.frame27.setGeometry(QRect(308,70,115,25))
        self.mscale = QLineEdit("1.0",self.frame27);
        self.mscale.setFixedWidth(75)
        validator = QDoubleValidator()
        self.mscale.setValidator(validator)
        self.frame16 = QFrame(self);self.frame16.setGeometry(QRect(205, 100, 250, 25))
        self.Label9 = QLabel("Marker First:",self.frame16)
        self.frame24 = QFrame(self);self.frame24.setGeometry(QRect(308,100,30,25))
        self.m1 = QCheckBox(self.frame24)
        self.m1.setCheckState(QtCore.Qt.Checked)
        self.frame17 = QFrame(self);self.frame17.setGeometry(QRect(205, 130, 250, 25))
        self.Label10 = QLabel("Column Spacing:",self.frame17)
        self.frame23 = QFrame(self);self.frame23.setGeometry(QRect(308,125,250,25))
        self.clmSpc = QLineEdit("2.0",self.frame23);
        self.clmSpc.setFixedWidth(75)
        validator = QDoubleValidator()
        self.clmSpc.setValidator(validator)

        painter = QPainter(self.labelC.pixmap())

        painter.setPen(QPen(Qt.black, 1, Qt.SolidLine))
        painter.drawLine(QPoint(5, 5), QPoint(5, 120))
        painter.drawLine(QPoint(5, 120), QPoint(375, 120))
        painter.drawLine(QPoint(375, 120), QPoint(375, 5))
        painter.drawLine(QPoint(375, 5), QPoint(5, 5))

        painter.end()
        self.frame20 =  QFrame(self);self.frame20.setGeometry(QRect(120, 170, 90, 25))
        self.applyBtn = QPushButton("Apply",self.frame20);
        self.applyBtn.clicked.connect(self.applyBtnAction)
        self.frame21 =  QFrame(self);self.frame21.setGeometry(QRect(210, 170, 90, 25))
        self.cancelBtn = QPushButton("Cancel",self.frame21)
        self.cancelBtn.clicked.connect(self.cancelBtnAction)
        self.frame22 =  QFrame(self);self.frame22.setGeometry(QRect(300, 170, 90, 25))
        self.okBtn = QPushButton("Ok",self.frame22)
        self.okBtn.clicked.connect(self.okBtnAction)
        self.linecolor = QColor('red');
        self.edgecolor = QColor('black');
        self.facecolor = QColor('red');

    def applyBtnAction(self):
        mainWin.legProp.setlocation(self.combo2.currentText())
        mainWin.legProp.setfontsize(float(self.fSize.text()))
        mainWin.legProp.settitle((self.legTitle.text()))
        mainWin.legProp.setlabelspacing(float(self.lblSpc.text()))
        mainWin.legProp.setmarkerscale(float(self.mscale.text()))
        mainWin.legProp.setcolumnspacing(float(self.clmSpc.text()))
        mainWin.legProp.setmarkerfirst(self.m1.isChecked())
        mainWin.legProp.setframeon(self.fr1.isChecked())
        mainWin.legendEnable = self.cb1.isChecked()
        if mainWin.legendEnable:
            mainWin.m.showLegend();
            mainWin.plotW.showLegend();
        else:
            mainWin.m.removeLegend()
            mainWin.plotW.removeLegend()

    def LegendPropShow(self):
        self.combo2.setCurrentText(mainWin.legProp.location)
        self.fSize.setText(str(mainWin.legProp.fontsize))
        self.legTitle.setText(mainWin.legProp.title)
        self.lblSpc.setText(str(mainWin.legProp.labelspacing))
        self.mscale.setText(str(mainWin.legProp.markerscale))
        self.clmSpc.setText(str(mainWin.legProp.columnspacing))
        bs = QtCore.Qt.Unchecked
        if mainWin.legProp.markerfirst :
            bs = QtCore.Qt.Checked            
        self.m1.setCheckState(bs)
        bs = QtCore.Qt.Unchecked
        if mainWin.legProp.frameon:
            bs = QtCore.Qt.Checked      
        self.fr1.setCheckState(bs)
        bs = QtCore.Qt.Unchecked
        if mainWin.legendEnable:
            bs = QtCore.Qt.Checked   
        self.cb1.setCheckState(bs)

    def cancelBtnAction(self):
        self.close();
    def okBtnAction(self):
        self.close();

class ApplicationWindow(QMainWindow):
    
    NoOfCol = 0;total_rows=0;

    def __init__(self, dir_output_files, file_last_opened):
        print('App')
        super().__init__()

        self.dir_output_files = dir_output_files
        self.file_last_opened = file_last_opened

        self.title = "Demo"
        self.left = 300
        self.top = 100
        self.width = 1000
        self.height = 620
        self.resize(self.width, self.height)
        self.move(self.left, self.top)
        
        self.widgetR();
        
    def widgetR(self):       
        self.colorSet = ["red","green","blue","magenta","black","orange","violet","brown"]
        self._main = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout()#._main)     
        self.DefaultPath = ""    
        self.msg = QMessageBox()
        self.frame1 = QFrame(self)
        self.frame1.setGeometry(QRect(10, 40, 100, 25))
        self.LegBtn1 = QPushButton("BrowseFile",self.frame1)
        self.LegBtn1.setChecked(True)
        self.LegBtn1.clicked.connect(self.openFileNameDialog)
        layout.addWidget(self.LegBtn1)
        
        self.frame6 = QFrame(self)
        self.frame6.setGeometry(QRect(300, 570, 1000, 30))
        self.labelFile = QLabel("Placeholder for file path and file name ",self.frame6)
        
        self.frame2 = QFrame(self)
        self.frame2.setGeometry(QRect(10, 237, 350, 30))
        self.label = QLabel("Number of columns:     ",self.frame2)
        layout.addWidget(self.label)
        self.label.setEnabled(0)

        self.frame5 = QFrame(self)
        self.frame5.setGeometry(QRect(10, 255, 350, 30))
        self.labelR = QLabel("Number of rows:                ",self.frame5)
        layout.addWidget(self.labelR)
        self.labelR.setEnabled(0)
        
        self.framegf = QFrame(self)
        self.framegf.setGeometry(QRect(300, 555, 100, 30))
        self.Labelf = QLabel("Filepath:",self.framegf)
        layout.addWidget(self.Labelf)
        
        layout.addWidget(self.labelFile)
        self.labelFile.setEnabled(0) 
        self.scrollFile = QScrollArea(self)
        self.scrollFile.setWidget(self.labelFile)
        self.scrollFile.setWidgetResizable(True)
        self.scrollFile.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        self.scrollFile.setFixedHeight(45)
        self.scrollFile.setFixedWidth(680)
        #self.scrollFile.setFrameShape(QFrame.NoFrame)
        layout.addWidget(self.scrollFile)
        self.scrollFile.move(300,570)
           
        self.framefn = QFrame(self)
        self.framefn.setGeometry(QRect(10, 65, 272, 30))
        self.labelFilen = QLabel("Placeholder for file name ",self.framefn)
        self.labelFilen.setEnabled(0)
        layout.addWidget(self.labelFile)
        self.scrollFilen = QScrollArea(self)
        self.scrollFilen.setWidget(self.labelFilen)
        self.scrollFilen.setWidgetResizable(True)
        self.scrollFilen.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff);
        self.scrollFilen.setFixedHeight(45)
        self.scrollFilen.setFixedWidth(172)
        self.scrollFilen.setFrameShape(QFrame.NoFrame)
        layout.addWidget(self.scrollFilen)
        self.scrollFilen.move(10,65)
        
        self.frameg8 = QFrame(self)
        self.frameg8.setGeometry(QRect(10, 117, 250, 30))
        self.Labelg = QLabel("Output Files:",self.frameg8)
        layout.addWidget(self.Labelg)
        
        self.frameg9 = QFrame(self)
        self.frameg9.setGeometry(QRect(10, 133, 200, 150))    
        self.OPFiles = QListWidget(self.frameg9)
        self.OPFiles.clicked.connect(self.FileSelectionChange)
        self.OPFiles.setSelectionMode(QAbstractItemView.NoSelection)
        layout.addWidget(self.OPFiles)
        self.scrollg2 = QScrollArea(self)
        self.scrollg2.setWidget(self.OPFiles)
        self.scrollg2.setWidgetResizable(True)
        self.scrollg2.setFixedHeight(100)
        self.scrollg2.setFixedWidth(150)
        layout.addWidget(self.scrollg2)
        self.scrollg2.move(10,133)
        self.FileIndex = None;
        
        self.frame3 = QFrame(self)
        self.frame3.setGeometry(QRect(10, 272, 250, 30))
        self.LabelX = QLabel("X-axis:",self.frame3)
        layout.addWidget(self.LabelX)
        self.frame7 = QFrame(self)
        self.frame7.setGeometry(QRect(10, 290, 200, 200))   
        self.senList = QListWidget(self.frame7)
        self.senList.clicked.connect(self.XSelectionChange)

        self.senList.setSelectionMode(QAbstractItemView.NoSelection)
        layout.addWidget(self.senList)
        self.scroll = QScrollArea(self)
        self.scroll.setWidget(self.senList)
        self.scroll.setWidgetResizable(True)
        self.scroll.setFixedHeight(150)
        self.scroll.setFixedWidth(150)
        layout.addWidget(self.scroll)
        self.scroll.move(10,290)
        self.XIndex = None;
        
        self.frame8 = QFrame(self)
        self.frame8.setGeometry(QRect(10, 447, 250, 30))
        self.LabelY = QLabel("Y-axis:",self.frame8)
        layout.addWidget(self.LabelY)
        
        self.frame9 = QFrame(self)
        self.frame9.setGeometry(QRect(10, 465, 300, 200))    
        self.YCols = QTableWidget(self.frame9)
        self.YCols.setColumnCount(3)
        self.YCols.setHorizontalHeaderLabels(["L","R",""])
        header = self.YCols.horizontalHeader()
        header.setSectionResizeMode(2, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.YCols.clicked.connect(self.YSelectionChange)
        self.YCols.setSelectionMode(QAbstractItemView.NoSelection)
        layout.addWidget(self.YCols)
        self.scroll2 = QScrollArea(self)
        self.scroll2.setWidget(self.YCols)
        self.scroll2.setWidgetResizable(True)
        self.scroll2.setFixedHeight(150)
        self.scroll2.setFixedWidth(270)
        layout.addWidget(self.scroll2)
        self.scroll2.move(10,465)
        self.YIndex = None;
        self.YIndexR = None;
        
        self.frame10 = QFrame(self)
        self.frame10.setGeometry(QRect(160, 125, 123, 25))
        self.linePropBtn = QPushButton("Plot Properties",self.frame10)
        self.linePropBtn.setChecked(True)
        self.linePropBtn.clicked.connect(self.plotProp)
        self.linePropBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.linePropBtn)
        
        self.frame11 = QFrame(self)
        self.frame11.setGeometry(QRect(160,150,123,25))
        self.axesBtn = QPushButton("Axes",self.frame11)
        self.axesBtn.clicked.connect(self.axesProp)
        self.axesBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.axesBtn)

        self.frame12 = QFrame(self)
        self.frame12.setGeometry(QRect(160,175,123,25))
        self.legendBtn = QPushButton("Legend",self.frame12)
        self.legendBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.legendBtn)
        self.legendBtn.clicked.connect(self.legendProp)
        
        self.frame15 = QFrame(self)
        self.frame15.setGeometry(QRect(160,200,123,25))
        self.gridBtn = QPushButton("Grid",self.frame15)
        self.gridBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.gridBtn)
        self.gridBtn.clicked.connect(self.gridProp)
        
        self.frame16 = QFrame(self)
        self.frame16.setGeometry(QRect(160,350,123,25))
        self.saveBtn = QPushButton("Save Figure",self.frame16)
        self.saveBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.saveBtn)
        self.saveBtn.clicked.connect(self.openFileSaveDialog)
        self.saveBtn.setEnabled(0)
        
        self.frame17 = QFrame(self)
        self.frame17.setGeometry(QRect(160,225,123,25))
        self.titleBtn = QPushButton("Title",self.frame17)
        self.titleBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.titleBtn)
        self.titleBtn.clicked.connect(self.titleProp)
        
        self.framexTicks = QFrame(self)
        self.framexTicks.setGeometry(QRect(160,250,123,25))
        self.xTicksBtn = QPushButton("xTicks",self.framexTicks)
        self.xTicksBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.xTicksBtn)
        self.xTicksBtn.clicked.connect(self.xTicksWProp)

        self.frameyTicks = QFrame(self)
        self.frameyTicks.setGeometry(QRect(160,275,123,25))
        self.yTicksBtn = QPushButton("yTicks",self.frameyTicks)
        self.yTicksBtn.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.yTicksBtn)
        self.yTicksBtn.clicked.connect(self.yTicksWProp)
        self.frameyTicks2 = QFrame(self)
        self.frameyTicks2.setGeometry(QRect(160,300,123,25))
        self.yTicksBtn2 = QPushButton("yTicks(Right)",self.frameyTicks2)
        self.yTicksBtn2.setFixedSize(QtCore.QSize(123, 25))
        layout.addWidget(self.yTicksBtn2)
        self.yTicksBtn2.clicked.connect(self.yTicksWProp2)
        self.frameTicksNt = QFrame(self)
        self.frameTicksNt.setGeometry(QRect(160,325,123,25))
        self.TicksNtBtn = QPushButton("TickLabelsNotation",self.frameTicksNt)
        layout.addWidget(self.TicksNtBtn)
        self.TicksNtBtn.clicked.connect(self.TicksNtWProp)
        
        self.m = PlotCanvas(self, width=5, height=3)
        self.m.move(285,40)
        self.navi = self.addToolBar(NavigationToolbar(self.m,self))
        
        self.plotW = PlotWindow()#self, width=5, height=3)
        self.plotW.move(1,1)
        self.navi2 = self.plotW.addToolBar(NavigationToolbar(self.plotW,self.plotW))
        
        self.frame13 = QFrame(self)
        self.frame13.setGeometry(QRect(600,520,120,25))
        self.scriptBtn = QPushButton("Generate Script",self.frame13)
        layout.addWidget(self.scriptBtn)
        self.scriptBtn.clicked.connect(self.genScript)
        self.scriptBtn.setEnabled(0)
        
        self.frame14 = QFrame(self)
        self.frame14.setGeometry(QRect(400,520,120,25))
        self.pltBtn = QPushButton("Plot Data",self.frame14)
        layout.addWidget(self.pltBtn)
        self.pltBtn.clicked.connect(self.plotWindow)
        self.pltBtn.setEnabled(0)
        self.header = 0
        self.legendEnable = True
        self.w=LinePropPopup();  
        self.w.setGeometry(QRect(700,530,400,200)) 
        self.plotAxesProp = AxesPropObject()           
        self.axesW = AxesPopup()
        self.axesW.setGeometry(QRect(700,550,600,180))  
        self.plotAxesProp.setxSNU(3); 
        self.plotAxesProp.setxSNL(-3);
        self.plotAxesProp.setySNU(3);
        self.plotAxesProp.setySNL(-3);
        self.plotAxesProp.setxSNM(True); 
        self.plotAxesProp.setySNM(True);
        
        self.legProp = LegendObject();             
        self.legendW = LegendPopup();
        self.legendW.setGeometry(QRect(700,530,400,200))
        self.gridObject = GridObject();             
        self.gridW = GridPopup();
        self.gridW.setGeometry(QRect(700,530,400,200))
        self.titleObject = TitleObject();             
        self.titleW = TitlePopup();
        self.titleW.setGeometry(QRect(700,530,300,200))
        
        self.headerW = Header();
        self.headerW.setGeometry(QRect(200,200,300,100))
        self.xTicksObject = xTicksPropObject();
        self.yTicksObject = xTicksPropObject();
        self.yTicksObject2 = xTicksPropObject();
        self.xTicksW =  xTicksPopup();
        self.xTicksW.setGeometry(QRect(700,530,480,200))
        self.yTicksW =  yTicksPopup();
        self.yTicksW.setGeometry(QRect(700,530,480,200))
        self.yTicksW2 =  yTicksPopup2();
        self.yTicksW2.setGeometry(QRect(700,530,480,200))
        self.TickLabelNotationW =  TickLabelNotation();
        self.TickLabelNotationW.setGeometry(QRect(700,530,480,200))
    def plotWindow(self):
        self.plotW.setGeometry(QRect(10,10,700,450))
        self.plotW.show()       
    def titleProp(self):
        self.titleW.TitlePropShow()
        self.titleW.show();       
    def yTicksWProp(self):
        self.yTicksW.yTicksPropShow()        
        self.yTicksW.show();
    def yTicksWProp2(self):
        self.yTicksW2.yTicksPropShow()        
        self.yTicksW2.show();
    def yTicksWProp2(self):
        self.yTicksW2.yTicksPropShow()        
        self.yTicksW2.show();
    def xTicksWProp(self):
        self.xTicksW.xTicksPropShow()        
        self.xTicksW.show();
    def genScript(self):
        self.scriptW = ScriptObject()
        self.scriptW.setGeometry(QRect(10,10,600,300))
        self.scriptW.show()
    def TicksNtWProp(self):
        self.TickLabelNotationW.show();
    def legendProp(self):
        self.legendW.LegendPropShow();
        self.legendW.show();
    def axesProp(self):
        self.axesW.AxesPropShow();
        self.axesW.show();
    def gridProp(self):
        self.gridW.GridPropShow();
        self.gridW.show();
    def plotProp(self):
        self.w.show() 
    def plotDataWithChangedOptions(self):
        self.m.clearPlot();self.plotW.clearPlot()
        yColor = []
        self.w.combo1.clear()
        if self.ext == '.csv':
            x = self.reader.iloc[:,self.XIndex]
        if (self.ext == '.dat') | (self.ext == '.txt'):
            x = self.reader[:,self.XIndex] 
        self.pls  = [];self.plbs = [];
        for i in range(0,len(self.YIndex)):
            pltColor = self.colorSet[(self.YIndex[i])%8]
            if self.colPlotObject[self.YIndex[i]].lineStyle == None: 
                self.colPlotObject[self.YIndex[i]].setLineStyle('solid')
            if self.colPlotObject[self.YIndex[i]].lineColor == None:
                self.colPlotObject[self.YIndex[i]].setLineColor(pltColor) 
            yColor.insert(i,(pltColor))
            if self.ext == '.csv':
                y = self.reader.iloc[:,self.YIndex[i]]
            if (self.ext == '.dat') | (self.ext == '.txt'):
                y = self.reader[:,self.YIndex[i]]
            pl = self.m.plotData(x,y,
                 self.colPlotObject[self.YIndex[i]].lineColor,
                 self.colPlotObject[self.YIndex[i]].lineStyle,
                 self.colPlotObject[self.YIndex[i]].width,
                 self.colPlotObject[self.YIndex[i]].drawStyle,
                 self.colPlotObject[self.YIndex[i]].label,
                 self.colPlotObject[self.YIndex[i]].marker,
                 self.colPlotObject[self.YIndex[i]].size, 
                 self.colPlotObject[self.YIndex[i]].edgeColor,
                 self.colPlotObject[self.YIndex[i]].faceColor);
            self.pls.append(pl);self.plbs.append(self.colPlotObject[self.YIndex[i]].label);
            self.plotW.plotData(x,y,
                 self.colPlotObject[self.YIndex[i]].lineColor,
                 self.colPlotObject[self.YIndex[i]].lineStyle,
                 self.colPlotObject[self.YIndex[i]].width,
                 self.colPlotObject[self.YIndex[i]].drawStyle,
                 self.colPlotObject[self.YIndex[i]].label,
                 self.colPlotObject[self.YIndex[i]].marker,
                 self.colPlotObject[self.YIndex[i]].size, 
                 self.colPlotObject[self.YIndex[i]].edgeColor,
                 self.colPlotObject[self.YIndex[i]].faceColor);
            self.w.combo1.addItem(self.YCols.item(self.YIndex[i],2).text())
            pixmap = QPixmap(10,10);
            pixmap.fill(QColor(self.colPlotObject[self.YIndex[i]].lineColor));
            self.redIcon= QIcon(pixmap);
            self.YCols.item(self.YIndex[i],0).setIcon(QIcon(self.redIcon))
        
        for i in range(0,len(self.YIndexR)):
            pltColor = self.colorSet[(self.YIndexR[i])%8]
            if self.colPlotObject[self.YIndexR[i]].lineStyle == None: 
                self.colPlotObject[self.YIndexR[i]].setLineStyle('solid')
            if self.colPlotObject[self.YIndexR[i]].lineColor == None:
                self.colPlotObject[self.YIndexR[i]].setLineColor(pltColor) 
            yColor.insert(i,(pltColor))
            if self.ext == '.csv':
                y = self.reader.iloc[:,self.YIndexR[i]]
            if (self.ext == '.dat') | (self.ext == '.txt'):
                y = self.reader[:,self.YIndexR[i]]
            pl = self.m.plotData(x,y,
                self.colPlotObject[self.YIndexR[i]].lineColor,
                self.colPlotObject[self.YIndexR[i]].lineStyle,
                self.colPlotObject[self.YIndexR[i]].width,
                self.colPlotObject[self.YIndexR[i]].drawStyle,
                self.colPlotObject[self.YIndexR[i]].label,
                self.colPlotObject[self.YIndexR[i]].marker,
                self.colPlotObject[self.YIndexR[i]].size, 
                self.colPlotObject[self.YIndexR[i]].edgeColor,
                self.colPlotObject[self.YIndexR[i]].faceColor,yax='right');
            self.pls.append(pl);self.plbs.append(self.colPlotObject[self.YIndexR[i]].label);
            self.plotW.plotData(x,y,
                self.colPlotObject[self.YIndexR[i]].lineColor,
                self.colPlotObject[self.YIndexR[i]].lineStyle,
                self.colPlotObject[self.YIndexR[i]].width,
                self.colPlotObject[self.YIndexR[i]].drawStyle,
                self.colPlotObject[self.YIndexR[i]].label,
                self.colPlotObject[self.YIndexR[i]].marker,
                self.colPlotObject[self.YIndexR[i]].size, 
                self.colPlotObject[self.YIndexR[i]].edgeColor,
                self.colPlotObject[self.YIndexR[i]].faceColor,yax='right');
            self.w.combo1.addItem(self.YCols.item(self.YIndexR[i],2).text())
            pixmap = QPixmap(10,10);
            pixmap.fill(QColor(self.colPlotObject[self.YIndexR[i]].lineColor));
            self.redIcon= QIcon(pixmap);
            self.YCols.item(self.YIndexR[i],1).setIcon(QIcon(self.redIcon))
        if self.legendEnable:
            self.m.showLegend();
            self.plotW.showLegend();
        self.plotAxesProp.setxLabel(self.senList.item(self.XIndex).text())
        
        self.m.getAxesProps();
        self.m.changeAxesProps();
        self.plotW.changeAxesProps();
        self.linePropBtn.setEnabled(1)               
        self.axesBtn.setEnabled(1)
        self.legendBtn.setEnabled(1)
        self.gridBtn.setEnabled(1)
        self.saveBtn.setEnabled(1)
        self.m.changeGridProps();
        self.plotW.changeGridProps();
        self.titleBtn.setEnabled(1)
        self.m.changeTitle()
        self.plotW.changeTitle(); 
        self.m.getxTicks();  
    def FileSelectionChange(self):
        NoOfItems = self.OPFiles.count()
        prev = self.FileIndex;
        listOfFileIndex = []
        for i in range(0,NoOfItems):
            item=self.OPFiles.item(i) 
            if item.checkState():
                listOfFileIndex.append(i)
        if len(listOfFileIndex)>1:
            for i in range(0,len(listOfFileIndex)):
                if self.FileIndex == listOfFileIndex[i]:
                    self.OPFiles.item(listOfFileIndex[i]).setCheckState(QtCore.Qt.Unchecked)
                else:
                    prev = listOfFileIndex[i]
        else:
            if len(listOfFileIndex) ==0:
                prev = None
            else:
                prev = listOfFileIndex[0]
        if self.FileIndex != prev:
            self.FileIndex = prev
            if self.FileIndex == None:
                NoOfItem =self.senList.count()
                self.label.setEnabled(0)
                self.labelR.setEnabled(0)
                self.scriptBtn.setEnabled(0)
                #self.labelFile.setEnabled(0)
                for ic in range(0,NoOfItem):
                    self.senList.takeItem(0)#self.senList.item(ic));
                    self.YCols.takeItem(ic,2);
                    self.YCols.takeItem(ic,1);
                    self.YCols.takeItem(ic,0);
                self.YCols.setRowCount(0)
                self.m.clearPlot()
                self.plotW.clearPlot()
                self.fileName = None
            else:
                self.fileName = self.dir_output_files + self.DataFiles[prev]
                if path.exists(self.fileName):
                    self.readFile();
                else:
                    self.labelR.setEnabled(0)
                    self.label.setEnabled(0)
                    self.msg.setText("Data File Does not Exist!")
                    self.showWarningDialog();
                    
    def XSelectionChange(self):
        NoOfItems = self.senList.count()
        prev = self.XIndex;
        listOfXIndex = []
        for i in range(0,NoOfItems):
            item=self.senList.item(i) 
            if item.checkState():
                listOfXIndex.append(i)
        if len(listOfXIndex)>1:
            for i in range(0,len(listOfXIndex)):
                if self.XIndex == listOfXIndex[i]:
                    self.senList.item(listOfXIndex[i]).setCheckState(QtCore.Qt.Unchecked)
                else:
                    prev = listOfXIndex[i]
              
        else:
            if len(listOfXIndex) ==0:
                prev = None
            else:
                prev = listOfXIndex[0]
        
        if self.XIndex != prev:
            self.XIndex = prev
            if self.XIndex == None:
                self.m.clearPlot()
                self.plotW.clearPlot()
            else:
                self.plotDataWithChangedOptions();
    def YSelectionChange(self):
        NoOfItems = self.NoOfCol#YCols.count()
        prev = self.YIndex;
        listOfYIndex = []
        listOfYIndexR = []
        for i in range(0,NoOfItems):
            item=self.YCols.item(i,0) 
            item2=self.YCols.item(i,1)
            pixmap = QPixmap(10,10);
            pixmap.fill(QColor("white"));
            self.redIcon= QIcon(pixmap);
            
            if (item.checkState()==2 and (item2.checkState()!=2)) or (item.checkState()==2 and \
                (self.FlagLR[i] == 'N' or self.FlagLR[i] == 'R')):
                
                listOfYIndex.append(i)
                self.FlagLR[i]='L'
                self.YCols.item(i,1).setCheckState(QtCore.Qt.Unchecked)
                self.YCols.item(i,1).setIcon(QIcon(self.redIcon))
               
            else:
                self.YCols.item(i,0).setIcon(QIcon(self.redIcon)) 
                
            if (item2.checkState()==2 and (item.checkState()!=2)) or \
               (item2.checkState()==2 and (self.FlagLR[i] == 'N' or self.FlagLR[i] == 'L')):
                listOfYIndexR.append(i)
                self.FlagLR[i]='R'
                self.YCols.item(i,0).setCheckState(QtCore.Qt.Unchecked)
                self.YCols.item(i,0).setIcon(QIcon(self.redIcon))
            else:
                self.YCols.item(i,1).setIcon(QIcon(self.redIcon)) 
                self.YCols.item(i,1).setCheckState(QtCore.Qt.Unchecked)
               
        if self.YIndex != listOfYIndex or self.YIndexR != listOfYIndexR:
            
            self.YIndex = listOfYIndex
            self.YIndexR = listOfYIndexR
            yColor = []
            if self.XIndex == None:
                self.m.clearPlot()
                self.plotW.clearPlot()
                yColor = []
                self.w.combo1.clear()
            else:
                self.plotDataWithChangedOptions();

    def openFileNameDialog(self):
        
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog#"All Files (*);;Python Files (*.py)"*.xlsx *.csv

        if path.exists(self.file_last_opened):
            with open(os.path.expanduser(self.file_last_opened),'r') as f:
                    netlistlines=f.readlines()
            self.DefaultPath=netlistlines[0]
        else:
            self.DefaultPath=""

        self.ProjectfileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", \
            self.DefaultPath," *.in", options=options)
        print(self.ProjectfileName)
        if path.exists(self.ProjectfileName):
            self.labelFile.setText(self.ProjectfileName)#'File: ...'+self.ProjectfileName[len(self.ProjectfileName)-25:])
            self.labelFile.setEnabled(1)
            self.labelFilen.setText(self.ProjectfileName[self.ProjectfileName.rindex('/')+1:])
            self.labelFilen.setEnabled(1)
            self.labelR.setEnabled(0)
            self.label.setEnabled(0)
            self.m.clearPlot()
            self.SolveBlkLine = []
            self.DefaultPath = self.ProjectfileName[:self.ProjectfileName.rindex('/')]
            with open(os.path.expanduser(self.file_last_opened), 'w+') as filetowrite:
                filetowrite.write(self.DefaultPath)
            self.netlist = self.dir_output_files + \
                self.ProjectfileName[self.ProjectfileName.rindex('/'):self.ProjectfileName.rindex('.')]+'.in'
            if path.exists(self.netlist):
                NoOfOPFiles =self.OPFiles.count()
                self.FileIndex = None
                for ic in range(0,NoOfOPFiles):
                    self.OPFiles.takeItem(0)

                with open(os.path.expanduser(self.netlist),'r') as f:
                    netlistlines=f.readlines()
                    self.SolveBlkLine = []
                    self.SolveBlkEndLine = []
                    self.OPBlkLine = []
                    self.OPBlkEndLine = []
                    self.DataFiles =[]
                     
                    LineNumber = 0
                    for line in netlistlines:
                        if 'begin_solve' in line:
                            self.SolveBlkLine.append(LineNumber)
                        if 'end_solve' in line:
                            self.SolveBlkEndLine.append(LineNumber)
                        if 'begin_output' in line:
                            self.OPBlkLine.append(LineNumber)
                            opFileLine = netlistlines[LineNumber+1]
                            dfile = opFileLine[opFileLine.index('=')+1:opFileLine.index('.dat')+4]
                            self.DataFiles.append(dfile)

                            itm = QListWidgetItem()
                            itm.setText(dfile)#+'.dat')
                            itm.setFlags(itm.flags() | QtCore.Qt.ItemIsUserCheckable)
                            itm.setCheckState(QtCore.Qt.Unchecked)

                            self.OPFiles.addItem(itm)
                        if 'end_output' in line:
                            self.OPBlkEndLine.append(LineNumber)
                        LineNumber+=1
                        
                    iblk =0;
                    self.variables=[]
                    self.varCnt = []
                    if len(self.SolveBlkLine)==0:
                        self.msg.setText("No Solve blocks for selected project!")
                        self.showWarningDialog();
                    else:
                        for blk in self.OPBlkLine:
                            for ln in range(blk, self.OPBlkEndLine[iblk]):
                                lin = netlistlines[ln]
                                if 'variables:' in lin:
                                    var = lin[lin.index(':')+2:-1].split(' ')
                                    self.variables.append(var)
                                    self.varCnt.append(len(var))
                            iblk+=1        

                NoOfItem =self.senList.count()
                for ic in range(0,NoOfItem):
                    self.senList.takeItem(0)#self.senList.item(ic));
                    self.YCols.takeItem(ic,2);
                    self.YCols.takeItem(ic,1);
                    self.YCols.takeItem(ic,0);
                self.YCols.setRowCount(0)
            else:
                self.msg.setText("Netlist file doesn't exist for selected project!")
                self.showWarningDialog();
                
    def openFileSaveDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.fileName, _ = QFileDialog.getSaveFileName(self,"QFileDialog.getSaveFileName()", 
            "","*.eps *.dat *.pgf *.pdf *.ps *.raw *.svg *.tiff", options=options)
        ax = self.m.fig.axes;
        self.m.fig.savefig(self.fileName)
    def showWarningDialog(self):
        
        self.msg.setIcon(QMessageBox.Information)

        self.msg.setWindowTitle("Warning!")
        self.msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        returnValue = self.msg.exec()
    def readFile(self):
        NoOfItem =self.senList.count()
        for ic in range(0,NoOfItem):
            self.senList.takeItem(0)#self.senList.item(ic));
            self.YCols.takeItem(ic,2);
            self.YCols.takeItem(ic,1);
            self.YCols.takeItem(ic,0);
        self.YCols.setRowCount(0)
        if self.fileName and path.exists(self.fileName):
            self.m.clearPlot()
            print(self.fileName)
            self.scriptBtn.setEnabled(1)
            self.ext = self.fileName[self.fileName.rindex('.'):]
            self.colPlotObject = [];
            if (self.ext == '.dat') | (self.ext == '.txt'):
                self.reader = np.loadtxt(self.fileName)    
            if str(self.reader)=='[]':
                self.msg.setText("Output file has no data!")
                self.showWarningDialog();
            else:           
                self.NoOfCol = self.reader.shape[1]
                self.total_rows = self.reader.shape[0]
                if self.NoOfCol == 0 or self.total_rows == 0:
                    self.msg.setText("Output file has no data!")
                    self.showWarningDialog();
                else:
                    fileI = self.DataFiles.index(self.fileName[self.fileName.rindex('/')+1:])                  
                    self.fileVariables=self.variables[fileI] 
                    if len(self.fileVariables)<self.NoOfCol:
                        self.fileVariables.insert(0,'time')
                    self.YCols.setRowCount(self.NoOfCol)
                    
                    self.YCols.setVerticalHeaderLabels(list([None]*self.NoOfCol))
                    self.FlagLR = []
                    for cln in range(0,self.NoOfCol):
                        if (self.ext == '.dat') | (self.ext == '.txt') | (self.header == 0):  
                            cl = self.fileVariables[cln]#"x"+str(cln)
                        self.FlagLR.append('N')
                        self.colPlotObject.append(PlotObject(label=cl))
                        #print(cl,self.colPlotObject)
                        itm =QListWidgetItem()
                        itm.setText(cl)
                        itm.setFlags(itm.flags() | QtCore.Qt.ItemIsUserCheckable)
                        itm.setCheckState(QtCore.Qt.Unchecked)
                        self.senList.addItem(itm)
                        pixmap = QPixmap(10,10);
                        pixmap.fill(QColor("white"));
                        self.redIcon= QIcon(pixmap);
                        itm =  QTableWidgetItem("")##QListWidgetItem()
                        itm.setFlags(itm.flags() | QtCore.Qt.ItemIsUserCheckable)
                        itm.setCheckState(QtCore.Qt.Unchecked)
                        itm.setIcon(self.redIcon)
                        self.YCols.setItem(cln,0,itm)
                        itm =  QTableWidgetItem("")##QListWidgetItem()
                        itm.setFlags(itm.flags() | QtCore.Qt.ItemIsUserCheckable)
                        itm.setCheckState(QtCore.Qt.Unchecked)                        
                        itm.setIcon(self.redIcon)
                        self.YCols.setItem(cln,1,itm)
                        self.YCols.verticalHeader().setSectionResizeMode(cln, QHeaderView.ResizeToContents)
                        itm =  QTableWidgetItem(cl)##QListWidgetItem()
                        self.YCols.setItem(cln,2,itm)
                        self.YCols.verticalHeader().setSectionResizeMode(cln, QHeaderView.ResizeToContents)
                    
                texC = "Number of columns:"+str(self.NoOfCol)                
                texR = "Number of Rows:"+str(self.total_rows) 
                self.label.setText(texC)
                self.label.setEnabled(1)
                self.labelR.setText(texR)
                self.labelR.setEnabled(1)
                #self.labelFile.setText(self.ProjectfileName)#('File: ...'+self.ProjectfileName[len(self.ProjectfileName)-25:])
                #self.labelFile.setEnabled(1)
                #self.labelFilen.setText(self.ProjectfileName[self.ProjectfileName.rindex('/')+1:])
                #self.labelFilen.setEnabled(1)
                if (self.NoOfCol) >=2:
                    self.XIndex = 0;
                    self.YIndex = [];#[1];
                    self.YIndexR = [];
                    self.senList.item(0).setCheckState(QtCore.Qt.Checked)
                    self.pltBtn.setEnabled(1)               

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=3, dpi=100):
        self.fig = Figure(figsize=(width, height))#,tight_layout=True)#, dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Fixed, QSizePolicy.Fixed)
        FigureCanvas.setGeometry(self,QRect(695,200,700,450))#self)
        FigureCanvas.updateGeometry(self)
        self.plot()

    def plot(self):
        data = [random.random() for i in range(250)]
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()
        self.ax.set_xlabel("X-axis")
        self.ax.set_ylabel("Y-axis")
        self.ax2 = self.ax.twinx()
        self.draw()
        
    def clearPlot(self):
        self.ax.cla()
        self.ax2.cla()
        self.ax.grid()
        self.draw()
        
    def plotData(self,x,y,pltColor='red',lineStyle='solid',lineWidth=0.7,drawStyle="default",
        pltLabel = '',markerStyle ='',markerSize = 7, markerEdgeColor = 'black',
        markerFaceColor='red',yax='left'):
        
        if yax =='right':
            pl=self.ax2.plot(x,y, color = pltColor,linestyle=lineStyle,
                linewidth=lineWidth,drawstyle=drawStyle,
                label=pltLabel, marker = markerStyle, markersize = markerSize,
                markeredgecolor = markerEdgeColor,markerfacecolor = markerFaceColor)
            self.ax.autoscale(enable=True,axis='x',tight=True)
            self.draw()
        else:
            pl=self.ax.plot(x,y, color = pltColor,linestyle=lineStyle,
                linewidth=lineWidth,drawstyle=drawStyle,
                label=pltLabel, marker = markerStyle, markersize = markerSize,
                markeredgecolor = markerEdgeColor,markerfacecolor = markerFaceColor)
            self.ax.autoscale(enable=True,axis='x',tight=True)
            self.draw()
        self.getAxesProps()
        return pl;
    def showLegend(self):
        ax = self.fig.axes
        ax = ax[0]
        if len(mainWin.YIndex) > 0 and len(mainWin.YIndexR) == 0:
            self.ax.legend(loc = mainWin.legProp.location,frameon = mainWin.legProp.frameon, 
                fontsize = mainWin.legProp.fontsize, title = mainWin.legProp.title,
                markerfirst = mainWin.legProp.markerfirst, markerscale = mainWin.legProp.markerscale, 
                labelspacing = mainWin.legProp.labelspacing, columnspacing = mainWin.legProp.columnspacing)
        elif len(mainWin.YIndexR)>0 and len(mainWin.YIndex) == 0:
            self.ax2.legend(loc = mainWin.legProp.location,frameon = mainWin.legProp.frameon, 
                fontsize = mainWin.legProp.fontsize, title = mainWin.legProp.title,
                markerfirst = mainWin.legProp.markerfirst, markerscale = mainWin.legProp.markerscale, 
                labelspacing = mainWin.legProp.labelspacing, columnspacing = mainWin.legProp.columnspacing)
        elif len(mainWin.YIndex) > 0 and len(mainWin.YIndexR) > 0:
            lns = mainWin.pls[0]
            for l in range(1,len(mainWin.pls)):
                lns = lns+mainWin.pls[l]
            self.ax.legend(lns,mainWin.plbs)
            
        self.draw()
        
    def removeLegend(self):
        ax = self.fig.axes
        ax = ax[0]
        ax.get_legend().remove()
        self.draw()
        
    def changeAxesProps(self):
        ax = self.fig.axes
        ax = ax[0]
        if mainWin.plotAxesProp.xScale != None:
            self.ax.set_xscale(mainWin.plotAxesProp.xScale)
        if mainWin.plotAxesProp.yScale != None and len(mainWin.YIndex)!=0:
            self.ax.set_yscale(mainWin.plotAxesProp.yScale)
        else:
            self.ax.set_yscale(mainWin.plotAxesProp.yScale2)
        if mainWin.plotAxesProp.yScale2 != None and len(mainWin.YIndexR)!=0:
            self.ax2.set_yscale(mainWin.plotAxesProp.yScale2) 
        else:
            self.ax2.set_yscale(mainWin.plotAxesProp.yScale) 
            
        if mainWin.plotAxesProp.xLabel != None:
            self.ax.set_xlabel(mainWin.plotAxesProp.xLabel)
        if mainWin.plotAxesProp.yLabel != None and len(mainWin.YIndex)!=0:
            self.ax.set_ylabel(mainWin.plotAxesProp.yLabel)
        else:
            self.ax.set_ylabel(mainWin.plotAxesProp.yLabel2)
        if mainWin.plotAxesProp.yLabel2 != None and len(mainWin.YIndexR)!=0:
            self.ax2.set_ylabel(mainWin.plotAxesProp.yLabel2)
        else:
            self.ax2.set_ylabel(mainWin.plotAxesProp.yLabel)
               
        if (mainWin.plotAxesProp.xMin != None) & (mainWin.plotAxesProp.xMax != None):
            self.ax.set_xlim(left = mainWin.plotAxesProp.xMin, right = mainWin.plotAxesProp.xMax)

        if (mainWin.plotAxesProp.yMin != None) & (mainWin.plotAxesProp.yMax != None) and len(mainWin.YIndex)!=0:
            self.ax.set_ylim(bottom = mainWin.plotAxesProp.yMin, top = mainWin.plotAxesProp.yMax)
        else:
            self.ax.set_ylim(bottom = mainWin.plotAxesProp.yMin2, top = mainWin.plotAxesProp.yMax2)
            
        if (mainWin.plotAxesProp.yMin2 != None) & (mainWin.plotAxesProp.yMax2 != None) and len(mainWin.YIndexR)!=0:
            self.ax2.set_ylim(bottom = mainWin.plotAxesProp.yMin2, top = mainWin.plotAxesProp.yMax2)
        else:
            self.ax2.set_ylim(bottom = mainWin.plotAxesProp.yMin, top = mainWin.plotAxesProp.yMax)
        if mainWin.plotAxesProp.xScale == 'linear':
            self.ax.ticklabel_format(axis="x", style="sci", scilimits=(mainWin.plotAxesProp.xSNL,\
                    mainWin.plotAxesProp.xSNU),useMathText=mainWin.plotAxesProp.xSNM)
        if mainWin.plotAxesProp.yScale == 'linear' and len(mainWin.YIndex)!=0:
            self.ax.ticklabel_format(axis="y", style="sci", scilimits=(mainWin.plotAxesProp.ySNL,\
            mainWin.plotAxesProp.ySNU),useMathText=mainWin.plotAxesProp.ySNM)
        else:
            self.ax.ticklabel_format(axis="y", style="sci", scilimits=(mainWin.plotAxesProp.ySNL2,\
            mainWin.plotAxesProp.ySNU2),useMathText=mainWin.plotAxesProp.ySNM2)
        if mainWin.plotAxesProp.yScale2 == 'linear' and len(mainWin.YIndexR)!=0:
            self.ax2.ticklabel_format(axis="y", style="sci", scilimits=(mainWin.plotAxesProp.ySNL2,\
            mainWin.plotAxesProp.ySNU2),useMathText=mainWin.plotAxesProp.ySNM2)
            
        else:
            self.ax2.ticklabel_format(axis="y", style="sci", scilimits=(mainWin.plotAxesProp.ySNL,\
            mainWin.plotAxesProp.ySNU),useMathText=mainWin.plotAxesProp.ySNM)
            
        self.draw()
            
    def getAxesProps(self):
        xMin, xMax = self.ax.get_xlim()  
        yMin, yMax = self.ax.get_ylim()
        yMin2, yMax2 = self.ax2.get_ylim()
        mainWin.plotAxesProp.setxMin(xMin); 
        mainWin.plotAxesProp.setxMax(xMax);
        mainWin.plotAxesProp.setyMin(yMin);
        mainWin.plotAxesProp.setyMax(yMax);
        mainWin.plotAxesProp.setyMin2(yMin2);
        mainWin.plotAxesProp.setyMax2(yMax2);

    def changeGridProps(self):
        if mainWin.gridObject.gridEnable:
            self.ax.grid(b=True,color = mainWin.gridObject.lineColor,
                    linestyle = mainWin.gridObject.lineStyle,
                    axis = mainWin.gridObject.axis,
                    which = mainWin.gridObject.which,
                    linewidth = mainWin.gridObject.width)
        else:
            self.ax.grid(b=None)

        self.draw()

    def changeTitle(self):
        if mainWin.titleObject.titleEnable:
            self.ax.set_title(label='')
            self.ax.set_title(label='',loc='left')
            self.ax.set_title(label='', loc='right')
            self.ax.set_title(label=mainWin.titleObject.label,loc=mainWin.titleObject.loc)
        else:
            self.ax.set_title(label='')
        
        self.draw()
        
    def getxTicks(self):
        mainWin.xTicksObject.setXTicks(list(self.ax.get_xticks()))
        mainWin.yTicksObject.setXTicks(list(self.ax.get_yticks()))
        mainWin.yTicksObject2.setXTicks(list(self.ax2.get_yticks()))
        bn = []
        for tick in self.ax.get_xticklabels():
            bn+= [str(tick.get_text())]
  
        mainWin.xTicksObject.setXTicksLabels(bn)
        bn = []
        for tick in self.ax.get_yticklabels():
            bn+= [str(tick.get_text())]
        mainWin.yTicksObject.setXTicksLabels(bn)
        
        bn = []
        for tick in self.ax2.get_yticklabels():
            bn+= [str(tick.get_text())]
        mainWin.yTicksObject2.setXTicksLabels(bn)
    def changeXTicks(self):
        if mainWin.xTicksObject.xTicksEnable:
            self.ax.set_xticks(mainWin.xTicksObject.xticks)#,mainWin.xTicksObject.direction)
            self.ax.tick_params(axis='x',direction=mainWin.xTicksObject.direction)
        else:
            self.ax.set_xticks([])
        self.draw()
    def changeXTicksLabels(self):
        if mainWin.xTicksObject.xTicksEnable:
            self.ax.set_xticks(mainWin.xTicksObject.xticks)
            self.ax.tick_params(axis='x',direction=mainWin.xTicksObject.direction)
            self.ax.set_xticklabels(mainWin.xTicksObject.xtickslabels,rotation = mainWin.xTicksObject.rotation)
        else:
            self.ax.set_xticks([])
        self.draw()
    def changeYTicks(self):
        if mainWin.yTicksObject.xTicksEnable:
            self.ax.set_yticks(mainWin.yTicksObject.xticks)#,mainWin.xTicksObject.direction)
            self.ax.tick_params(axis='y',direction=mainWin.yTicksObject.direction)
        else:
            self.ax.set_yticks([])
        self.draw()
    def changeYTicksLabels(self):
        if mainWin.yTicksObject.xTicksEnable:
            self.ax.set_yticks(mainWin.yTicksObject.xticks)
            self.ax.tick_params(axis='y',direction=mainWin.yTicksObject.direction)
            self.ax.set_yticklabels(mainWin.yTicksObject.xtickslabels,rotation = mainWin.yTicksObject.rotation)
        else:
            self.ax.set_yticks([])
        self.draw()    
    def changeYTicks2(self):
        if mainWin.yTicksObject2.xTicksEnable:
            self.ax2.set_yticks(mainWin.yTicksObject2.xticks)#,mainWin.xTicksObject.direction)
            self.ax2.tick_params(axis='y',direction=mainWin.yTicksObject2.direction)
        else:
            self.ax.set_yticks([])
        self.draw()
    def changeYTicksLabels2(self):
        if mainWin.yTicksObject2.xTicksEnable:
            self.ax2.set_yticks(mainWin.yTicksObject2.xticks)
            self.ax2.tick_params(axis='y',direction=mainWin.yTicksObject2.direction)
            self.ax2.set_yticklabels(mainWin.yTicksObject2.xtickslabels,
                rotation = mainWin.yTicksObject2.rotation)
        else:
            self.ax2.set_yticks([])
        self.draw()    
class PlotWindow(QMainWindow, PlotCanvas):
    def __init__(self):
        super().__init__()
if __name__ == "__main__":
    QGuiApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True);
    app = QtWidgets.QApplication([])

#   Note: include / at the end in passing dir_output_files
    if len(sys.argv) != 3:
        print('plot_main: need 3 arguments. Halting...')
        sys.exit()
    dir_output_files = sys.argv[1]
    file_last_opened = sys.argv[2]

    mainWin = ApplicationWindow(dir_output_files, file_last_opened)

    mainWin.show()
    sys.exit( app.exec_())
