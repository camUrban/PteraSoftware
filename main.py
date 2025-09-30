# NOTE: I haven't yet started refactoring this module.
"""This script opens the GUI to demonstrate Ptera Software analyzing example models.
It is in development and will be able to run custom models in the future"""

import os
import sys
import time
import importlib

print("Builtin modules imported")
from PySide6.QtCore import Qt

print("QTCore imported")
from PySide6.QtGui import QPixmap

print("QtGUI imported")
from PySide6.QtWidgets import QMainWindow, QApplication, QSplashScreen, QDialog

from pterasoftware.ui_resources.main_window import Ui_MainWindowDesign
from pterasoftware.ui_resources.textdialog import Ui_TextAboutDialog


class TextAboutDialog(QDialog):
    def __init__(self, title):
        super(TextAboutDialog, self).__init__()
        self.ui = Ui_TextAboutDialog()
        self.ui.setupUi(self)
        self.setWindowTitle(title)


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

        self.actionAbout.triggered.connect(self.menuREADME)

        self.displayText = ""

    def exampleMenu(self, ex_num):
        files = []
        for i, filename in enumerate(os.listdir("examples")):
            f = "examples." + filename
            files.append(f)
        file = files[ex_num]
        file = file.replace(".py", "")
        self.printTerminalOutput(f"Example {ex_num + 1} executed")
        print(file)
        importlib.import_module(file)

    def printTerminalOutput(self, text):
        print("Printing terminal output")
        self.terminalOutput.addItem(f"{text}")
        self.terminalOutput.scrollToBottom()

    def updateDisplayText(self, text):
        self.displayText = text
        self.printTerminalOutput(self, text)

    def menuREADME(self):
        from PySide6.QtGui import QTextDocument

        self.dialog = TextAboutDialog("About Ptera Software")
        doc = QTextDocument()
        doc.setMarkdown(self._read_file("README.md"))
        self.dialog.ui.textEdit.setDocument(doc)
        self.dialog.show()

    def _read_file(self, file_path: str) -> str:
        from PySide6.QtCore import QFile
        from PySide6.QtCore import QTextStream
        from PySide6.QtCore import QIODevice

        file = QFile(file_path)
        file.open(QIODevice.ReadOnly)
        ts = QTextStream(file)
        string = ts.readAll()
        return string


if __name__ == "__main__":
    app = QApplication(sys.argv)
    pixmap = QPixmap("docs/logo.png")
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
    sys.exit(app.exec())
    print("Exiting")
