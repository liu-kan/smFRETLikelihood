# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 12:42:04 2017

@author: liuk
"""

#from pyqtgraph.Qt import QtGui, QtCore,QtWidgets  # (the example applies equally well to PySide)
#import pyqtgraph as pg
import matplotlib
#from maxlikh import MaxLikehood


from PyQt5 import QtGui, QtCore,QtWidgets
import sys,os
import multiprocessing
import numpy as np
from ui.mprogress import *
from ui.pyqtMatplot import pyqtMatplot
from algo.calcThread import *
class Window(QtWidgets.QWidget):
    def __init__(self):
        super(Window, self).__init__()
        self.Sth=1.0
        self.burst=None
        self.fretBtn = QtWidgets.QPushButton('E-FRET Histogram')
        self.thresholdEBtn = QtWidgets.QPushButton('select threshold of E')
        self.maxLikBtn = QtWidgets.QPushButton('MaxLikiHood')
        self.filename = QtWidgets.QLineEdit('/home/liuk/proj/data/86c_224c.sqlite')
        logList = QtWidgets.QListWidget()
        self.n_states = QtWidgets.QLineEdit('2')

        self.fretBtn.clicked.connect(self.handleE_FRET)
        self.thresholdEBtn.clicked.connect(self.handleE_threshold)

        self.maxLikBtn.clicked.connect(self.handleE_maxlh)
        layout = QtWidgets.QGridLayout()
        sb=QtWidgets.QVBoxLayout()
        hb=QtWidgets.QHBoxLayout()
        hb.setSpacing(0)
        self.setLayout(layout)
        sb.addWidget(self.filename)
        sb.addWidget(self.fretBtn)
        sb.addWidget(self.thresholdEBtn)
        sb.addWidget(self.n_states)
        sb.addWidget(self.maxLikBtn)

        sb.addWidget(logList)

        wl = QtWidgets.QWidget()
        wl.setLayout(sb)
        layout.addWidget(wl, 0, 0,1,1)   # button goes in upper-left

        # v = pg.GraphicsView()
        # self.vb = pg.ViewBox()
        # self.vb.setAspectLocked()
        # v.setCentralItem(self.vb)

        # hb.addWidget(v)
        # wr = QtGui.QWidget()
        # wr.setLayout(hb)
        wr=pyqtMatplot(parent=self)
        layout.addWidget(wr, 0, 1,5,1)
        self.matplot=wr

    def handleE_FRET(self):
        if os.path.exists(self.filename.text()):
            t= e_sHistThread(self.filename.text(),self.matplot,self)
            progressWidget = progressDlg(t,10,parent=self,title="calc E-S Histogram ...")
            progressWidget.exec_()
        else:
            QtWidgets.QMessageBox.warning(self,'Warning','The File does not exist!')
    def handleE_threshold(self):
        self.matplot.canvas.updatingLine=not self.matplot.canvas.updatingLine
        print(self.matplot.canvas.updatingLine)
    def handleE_maxlh(self):
        print("SyncResolution",self.burst["SyncResolution"])
        t= maxlikThread(self.burst,int(self.n_states.text()))
        progressWidget = progressDlg(t,10,parent=self,title="calc maxLikelihood of photons ...")
        progressWidget.exec_()

if __name__ == '__main__':
    if os.name!='posix':
        print("no fork,freeze_support need.")
        multiprocessing.freeze_support()
    import sys
    app = QtWidgets.QApplication(sys.argv)
    #pg.setConfigOptions(antialias=True)
    window = Window()
    window.setWindowTitle("SMFRET")
    window.show()
    sys.exit(app.exec_())
