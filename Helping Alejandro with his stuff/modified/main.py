# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 13:53:01 2017

@author: gerar
"""
from PyQt5 import QtGui, QtWidgets
import sys
import interfaz_prueba

class ExampleApp(QtWidgets.QMainWindow, interfaz_prueba.Ui_MainWindow):
    def __init__(self,parent=None):
        super(ExampleApp,self).__init__(parent)
        self.setupUi(self)

def main():
    app = QtWidgets.QApplication(sys.argv)
    form =ExampleApp()
    form.show()
    app.exec_()

if __name__ == "__main__":
    main()