# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 11:22:22 2017

@author: liuk
"""

from PyQt5 import QtCore, QtGui,QtWidgets
import sys, time

class mythread(QtCore.QThread):

    total = QtCore.pyqtSignal(int)
    update = QtCore.pyqtSignal()

    def __init__(self, n, parent=None):
        super(mythread, self).__init__(parent)
        self.n = n

    def run(self):
        i = 0
        #self.total.emit(self.n)
        oldtime=0
        while (i<self.n):
            if (time.time()-oldtime>=1):
                i+=1
                print( str(i))
                self.update.emit()
                oldtime=time.time()
# create the dialog for zoom to point
class progress(QtWidgets.QProgressBar):

    def __init__(self,thread, len=10,textVisible=False ,parent=None):
        super(progress, self).__init__(parent)
        # Set up the user interface from Designer.
        self.setTextVisible(textVisible)
        self.setValue(0)

        self.thread =thread
        self.setMaximum(len)
        self.thread.total.connect(self.setMaximum)
        self.thread.update.connect(self.update)
        self.thread.finished.connect(self.close)

        self.n = 0
        self.thread.start()

    def update(self):
        self.n += 1
        #print( self.n)
        self.setValue(self.n%self.maximum())
class progressDlg(QtWidgets.QDialog):
    def __init__(self,thread, len=10,textVisible=False,parent=None,title=None):
        super().__init__(parent=parent)
        self.resize(300, 60)
        self.thread=thread
        vb=QtWidgets.QVBoxLayout()
        progressWidget = progress(thread, len,textVisible,parent=self)
        vb.addWidget(progressWidget)
        self.setLayout(vb)
        self.setWindowModality(QtCore.Qt.ApplicationModal)
        self.thread.finished.connect(self.close)
        if(title!=None):
            self.setWindowTitle(title)
        if(parent!=None):
            resultsY = parent.geometry().center().y()-self.geometry().height()/2 ;
            resultsX = parent.geometry().center().x()-self.geometry().width()/2;
            self.move(resultsX, resultsY);
if __name__=="__main__":
    app = QtWidgets.QApplication([])
    t= mythread( 20)
    progressWidget = progress(t,10)
    progressWidget.move(300, 300)
    progressWidget.show()
    sys.exit(app.exec_())
