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

from pterasoftware.ui_resources import mainWindow

class MainWindow(QMainWindow, mainWindow):
    if 0==0:
        print("Duh")


if __name__ == '__main__':


    print("Yup, it works")





