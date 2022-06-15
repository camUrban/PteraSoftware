import os
import sys
import time

print("Builtin modules imported")
from PySide2.QtCore import Qt

print("QTCore imported")
from PySide2.QtGui import QPixmap

print("QtGUI imported")
from PySide2.QtWidgets import QMainWindow, QApplication, QSplashScreen

from pterasoftware.ui_resources.mainWindow import Ui_MainWindowDesign


class MainWindow(QMainWindow, Ui_MainWindowDesign):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setupUi(self)

        self.actionExample_1.triggered.connect(lambda x: self.exampleMenu(0))
        self.actionExample_2.triggered.connect(lambda x: self.exampleMenu(1))
        self.actionExample_3.triggered.connect(lambda x: self.exampleMenu(2))
        self.actionExample_4.triggered.connect(lambda x: self.exampleMenu(3))
        self.actionExample_5.triggered.connect(lambda x: self.exampleMenu(4))
        self.actionExample_6.triggered.connect(lambda x: self.exampleMenu(5))
        self.actionExample_7.triggered.connect(lambda x: self.exampleMenu(6))
        self.actionExample_8.triggered.connect(lambda x: self.exampleMenu(7))
        self.actionExample_9.triggered.connect(lambda x: self.exampleMenu(8))
        self.actionExample_10.triggered.connect(lambda x: self.exampleMenu(9))

    def exampleMenu(self, ex_num):
        import importlib
        files = []
        for i, filename in enumerate(os.listdir("examples")):
            f = "examples." + filename
            files.append(f)
        file = files[ex_num]
        file = file.replace(".py", "")
        print(file)
        importlib.import_module(file)

    def printTerminalOutput(self):
        pass


if __name__ == '__main__':
    app = QApplication(sys.argv)
    pixmap = QPixmap('docs/logo.png')
    splash = QSplashScreen(pixmap)
    splash.setWindowFlags(Qt.WindowStaysOnTopHint)
    splash.setEnabled(False)
    splash.setMask(pixmap.mask())
    splash.show()
    time.sleep(1)
    app.processEvents()

    window = MainWindow()
    window.show()
    window.raise_()
    window.activateWindow()
    splash.finish(window)
    sys.exit(app.exec_())

    print("Yup, it works")
