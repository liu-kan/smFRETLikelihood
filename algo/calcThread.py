from PyQt5 import QtCore
import algo.BurstSearch as BurstSearch
import algo.BGrate as BGrate
import algo.fretAndS as fretAndS
import algo.maxLik as maxLik
#import pyqtgraph as pg
import ui.pyqtMatplot
import copy
class e_sHistThread(QtCore.QThread):
    total = QtCore.pyqtSignal(int)
    update = QtCore.pyqtSignal()
    def __init__(self, dbname,viewbox, parent=None):
        super(e_sHistThread, self).__init__(parent)
        self.dbname = dbname
        self.count = 0
        self.timer=QtCore.QTimer(self)
        self.timer.timeout.connect(self.updateCount)
        self.timer.start(1000)
        self.viewbox=viewbox
        self.mainui=None
        if parent is not None:
            self.mainui=parent
    def updateCount(self):
        self.count = self.count + 1
        self.update.emit()
    def run(self):
        self.br=BGrate.calcBGrate(self.dbname,20,400)
        self.burst=BurstSearch.findBurst(self.br,self.dbname,["All"])
        self.burstSeff, self.burstFRET,self.wei,self.H,self.xedges, \
            self.yedges=fretAndS.FretAndS(self.dbname,self.burst,(27,27),self.br)
        #self.timer.stop()
        self.viewbox.canvas.drawHist(self.H,self.xedges,self.yedges)
        if self.mainui is not None:
            print ("copy burst")
            self.mainui.burst=copy.deepcopy(self.burst)
            print(self.mainui.burst["SyncResolution"])

class maxlikThread(QtCore.QThread):
    """docstring for maxlikThread """
    total = QtCore.pyqtSignal(int)
    update = QtCore.pyqtSignal()
    def __init__(self, burst,n_states,parent=None):
        super(maxlikThread,self).__init__()
        self.count = 0
        self.timer=QtCore.QTimer(self)
        self.timer.timeout.connect(self.updateCount)
        self.timer.start(1000)
        self.burst=burst
        self.n_states=n_states
        #self.viewbox=viewbox
    def updateCount(self):
        self.count = self.count + 1
        self.update.emit()
    def run (self):
        gsml=maxLik.GS_MLE(self.burst,0.891)
        gsml.n_states=self.n_states
        gsml.MaxLikehood([0.3,0.7,0.2, 3,3,3, 3,3,3])
        self.timer.stop()
