from PyQt5 import QtGui, QtCore,QtWidgets
from PyQt5.QtWidgets import QVBoxLayout

import matplotlib as mpl
mpl.use('Qt5Agg')
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import matplotlib.cm as cm
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

class pyqtMatplot(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(pyqtMatplot, self).__init__(parent)
        self.mainui=parent
        self.setupUi()
    def setupUi(self):
        vbox = QVBoxLayout()
        figure = mpl.figure.Figure(figsize=(10, 10))
        self.canvas = MatplotlibWidget(figure,self.mainui)
        #vbox.addWidget(self.mpl_toolbar)
        self.canvas.setParent(self)
        vbox.addWidget(self.canvas)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        vbox.addWidget(self.mpl_toolbar)
        self.setLayout(vbox)


class MatplotlibWidget(FigureCanvasQTAgg):
    def __init__(self,fig,mainui):
        super(MatplotlibWidget, self).__init__(fig)
        #-- set up an axe artist --
        self.fig=fig
        self.mainui=mainui
        self.ax = self.fig.add_axes([0.05, 0.05, 0.9, 0.9])
        # im=self.dax.imshow(H.transpose()[::-1], interpolation='bessel',
        #               cmap=cm.jet,extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        # fig.colorbar(im)
        # self.draw()
        self.updatingLine=False
        self.mpl_connect('motion_notify_event', self.onMouseMove)
        self.mpl_connect('button_press_event', self.onclick)
    def drawHist(self,H,xedges,yedges):
        im=self.ax.imshow(H.transpose()[::-1], interpolation='bessel',
                       cmap=cm.jet,extent=[0,1,0,1])
                       #[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        self.fig.colorbar(im)
        self.draw()
    def onMouseMove(self,event):
        if not event.inaxes:
            return
        if self.updatingLine :#and len(self.dax.lines)>0:
            if len(self.ax.lines)<1:
                self.ax.axhline(y=event.ydata, color="k")
            else:
                self.ax.lines = []#self.ax.lines[-1:] # keep the first two lines
                self.ax.axhline(y=event.ydata, color="k")
                self.draw()

    def onclick(self, event):
        if not event.inaxes:
            return
        if self.updatingLine:
            self.updatingLine=not self.updatingLine
        if len(self.ax.lines)>0:
            print(self.ax.lines[-1].get_ydata()[0])
            if self.mainui is not None:
                self.mainui.Sth=self.ax.lines[-1].get_ydata()[0]
