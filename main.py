import os
import sys
import time

import PySide2

from pterasoftware import *


print("Builtin modules imported")
from PySide2.QtCore import Qt, QRect, Slot, QThreadPool

print("QTCore imported")
from PySide2.QtGui import QPixmap, QCloseEvent, QPalette, QColor

print("QtGUI imported")
from PySide2.QtWidgets import QDialog, QMainWindow, QApplication, QListWidgetItem, QSplashScreen, QMessageBox, QSlider, \
    QLabel, QAbstractItemView, QWidget, QColorDialog, QMenu

from pterasoftware.ui_resources.mainWindow import Ui_MainWindowDesign

class MainWindow(QMainWindow, Ui_MainWindowDesign):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)


if __name__ == '__main__':

    app = QApplication(sys.argv)
    pixmap = QPixmap('docs/logo.png')
    splash = QSplashScreen(pixmap)
    splash.setWindowFlags(Qt.WindowStaysOnTopHint)
    splash.setEnabled(False)
    splash.setMask(pixmap.mask())
    splash.show()
    time.sleep(0.1)  # This seems to fix the splash mask displaying but not the actual image
    app.processEvents()

    window = MainWindow()
    window.show()
    window.raise_()
    window.activateWindow()
    splash.finish(window)
    sys.exit(app.exec_())

    print("Yup, it works")





