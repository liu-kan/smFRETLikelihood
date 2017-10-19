import matplotlib
#import matplotlib as mpl
matplotlib.use('Qt5Agg')

from scipy.optimize import leastsq

import matplotlib.cm as cm

from PyQt5.QtWidgets import QVBoxLayout
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
import sys

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import matplotlib.pyplot as plt
from PyQt5 import QtWidgets

class ControlMainWindow(QtWidgets.QMainWindow):
    def __init__(self,H,xedges,yedges, parent=None):
        super(ControlMainWindow, self).__init__(parent)
        self.H=H
        self.xedges=xedges
        self.yedges=yedges
        self.main_frame = QtWidgets.QWidget()
        self.setupUi()

    def setupUi(self):
        vbox = QVBoxLayout()
        figure = matplotlib.figure.Figure(figsize=(10, 10))
        leftcanvas = MatplotlibWidget(figure,self.H,self.xedges,self.yedges)
        #vbox.addWidget(self.mpl_toolbar)
        leftcanvas.setParent(self.main_frame)
        vbox.addWidget(leftcanvas)
        self.mpl_toolbar = NavigationToolbar(leftcanvas, self.main_frame)
        vbox.addWidget(self.mpl_toolbar)
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

class MatplotlibWidget(FigureCanvasQTAgg):
    def __init__(self, fig,H,xedges,yedges):
        super(MatplotlibWidget, self).__init__(fig)

        #-- set up an axe artist --

        self.dax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        im=self.dax.imshow(H.transpose()[::-1], interpolation='bessel',
                      cmap=cm.jet,extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        fig.colorbar(im)
        self.rs=RectSect(self.dax,xedges,yedges)
        self.rb=RectBuilder(self.dax,xedges,yedges,self.rs.toggle_selectorRS)
        self.draw()
        self.updatingLine=True

        #---- setup event ----
        #self.mpl_connect('key_press_event', self.onkeypress)
        #self.mpl_connect('motion_notify_event', self.onMouseMove)
        #self.mpl_connect('button_press_event', self.onclick)
    def onkeypress(self, event):
        print(' Key pressed.')
        if event.key in ['D', 'd'] and self.rs.active:
            print(' RectangleSelector deactivated.')
            self.rs.set_active(False)
        if event.key in ['A', 'a'] and not self.rs.active:
            print(' RectangleSelector activated.')
            self.rs.set_active(True)
    def onMouseMove(self,event):
        if not event.inaxes:
            return
        if self.updatingLine:
            self.dax.lines = [self.dax.lines[0]] # keep the first two lines
            self.dax.axhline(y=event.ydata, color="k")
            self.draw()

    def onclick(self, event):
        if not event.inaxes:
            return
        self.updatingLine=not self.updatingLine


def betweenXY(v,arr):
    si=-1
    ei=-1
    if v==None:
        return si,ei
    for i, av in enumerate(arr):
        if ei==-1 and av<=v:
            si=i
        if si!=-1 and av>=v:
            ei=i
            break
    return si,ei

class toggle_selector:
    def __init__(self, rs,pltc):

        self.rs=rs
        rs.set_active(False)
        pltc.connect('key_press_event', self)
    def __call__(self, event):
        print(' Key pressed.')
        if event.key in ['D', 'd'] and self.rs.active:
            print(' RectangleSelector deactivated.')
            self.rs.set_active(False)
        if event.key in ['A', 'a'] and not self.rs.active:
            print(' RectangleSelector activated.')
            self.rs.set_active(True)

class RectSect:
    def __init__(self, ax,xedges,yedges):

        self.xedges=xedges
        self.yedges=yedges
        self.rects = ax
        self.toggle_selectorRS = matplotlib.widgets.RectangleSelector(ax, self,
                                       drawtype='box', useblit=False,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=False)
        self.toggle_selectorRS.set_active(False)
        self.n_m1=lambda x: True if x!=-1 else False

    def __call__(self, eclick,erelease):
        'eclick and erelease are the press and release events'
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        #print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
        #print(" The button you used were: %s %s" % (eclick.button, erelease.button))

        ax=betweenXY(x1,self.xedges)
        ay=betweenXY(y1,self.yedges)

        ax2=betweenXY(x2,self.xedges)
        ay2=betweenXY(y2,self.yedges)

        if all(map(self.n_m1,ax)) and all(map(self.n_m1,ay)) and all(map(self.n_m1,ax2)) and all(map(self.n_m1,ay2)):
            rx=self.xedges[ax2[1]]-self.xedges[ax[0]]
            ry=self.yedges[ay2[1]]-self.yedges[ay[0]]
            print(rx,ry)
            b=matplotlib.patches.Rectangle((self.xedges[ax[0]],
                                        self.yedges[ay[0]]), rx,
                                       ry,color='black')
            self.rects.add_patch(b)
            self.rects.figure.canvas.draw()
        print(len(self.rects.patches))

class RectBuilder:
    def __init__(self, rect,xedges,yedges,rs):
        self.rects = rect
        self.rx=xedges[1]-xedges[0]
        self.ry=yedges[1]-yedges[0]
        self.rs=rs
        self.xedges=xedges
        self.yedges=yedges
        self.cid = rect.figure.canvas.mpl_connect('button_press_event', self)
        self.n_m1=lambda x: True if x!=-1 else False

    def __call__(self, event):
        #print('click', event)
        if self.rs.active: return
        if event.inaxes!=self.rects.axes: return

        ax=betweenXY(event.xdata,self.xedges)
        ay=betweenXY(event.ydata,self.yedges)

        if all(map(self.n_m1,ax)) and all(map(self.n_m1,ay)):
            b=matplotlib.patches.Rectangle((self.xedges[ax[0]],
                                        self.yedges[ay[0]]), self.rx,
                                       self.ry,color='black')
            #print((self.xedges[ax[0]],self.yedges[ay[0]]))
            self.rects.add_patch(b)
            self.rects.figure.canvas.draw()
        print(len(self.rects.patches))
