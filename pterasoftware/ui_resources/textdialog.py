# NOTE: I haven't yet started refactoring this module.
from PySide6 import QtCore, QtWidgets


class Ui_TextAboutDialog(object):
    def setupUi(self, TextAboutDialog):
        TextAboutDialog.setObjectName("TextAboutDialog")
        TextAboutDialog.resize(1000, 900)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.MinimumExpanding,
            QtWidgets.QSizePolicy.MinimumExpanding,
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(TextAboutDialog.sizePolicy().hasHeightForWidth())
        TextAboutDialog.setSizePolicy(sizePolicy)
        TextAboutDialog.setMinimumSize(QtCore.QSize(1000, 900))
        TextAboutDialog.setSizeGripEnabled(False)
        TextAboutDialog.setModal(True)
        self.textEdit = QtWidgets.QTextEdit(TextAboutDialog)
        self.textEdit.setEnabled(True)
        self.textEdit.setGeometry(QtCore.QRect(0, 0, 1000, 900))
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.MinimumExpanding,
            QtWidgets.QSizePolicy.MinimumExpanding,
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.textEdit.sizePolicy().hasHeightForWidth())
        self.textEdit.setSizePolicy(sizePolicy)
        self.textEdit.setMinimumSize(QtCore.QSize(1000, 900))
        self.textEdit.setMaximumSize(QtCore.QSize(0, 0))
        self.textEdit.setAutoFormatting(QtWidgets.QTextEdit.AutoAll)
        self.textEdit.setTextInteractionFlags(QtCore.Qt.NoTextInteraction)
        self.textEdit.setObjectName("textEdit")

        self.retranslateUi(TextAboutDialog)
        QtCore.QMetaObject.connectSlotsByName(TextAboutDialog)

    def retranslateUi(self, TextAboutDialog):
        TextAboutDialog.setWindowTitle(
            QtWidgets.QApplication.translate("TextAboutDialog", "Dialog", None, -1)
        )
