# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'mainWindow.ui'
##
## Created by: Qt User Interface Compiler version 5.15.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *


class Ui_MainWindowDesign(object):
    def setupUi(self, MainWindowDesign):
        if not MainWindowDesign.objectName():
            MainWindowDesign.setObjectName(u"MainWindowDesign")
        MainWindowDesign.resize(1134, 887)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindowDesign.sizePolicy().hasHeightForWidth())
        MainWindowDesign.setSizePolicy(sizePolicy)
        MainWindowDesign.setMinimumSize(QSize(1052, 880))
        self.actionAbout = QAction(MainWindowDesign)
        self.actionAbout.setObjectName(u"actionAbout")
        self.actionManual = QAction(MainWindowDesign)
        self.actionManual.setObjectName(u"actionManual")
        self.actionQuit = QAction(MainWindowDesign)
        self.actionQuit.setObjectName(u"actionQuit")
        self.actionImport_Data = QAction(MainWindowDesign)
        self.actionImport_Data.setObjectName(u"actionImport_Data")
        self.actionExport_Data = QAction(MainWindowDesign)
        self.actionExport_Data.setObjectName(u"actionExport_Data")
        self.actionPreset_Aircraft_List = QAction(MainWindowDesign)
        self.actionPreset_Aircraft_List.setObjectName(u"actionPreset_Aircraft_List")
        self.centralWidget = QWidget(MainWindowDesign)
        self.centralWidget.setObjectName(u"centralWidget")
        self.gridLayout_4 = QGridLayout(self.centralWidget)
        self.gridLayout_4.setSpacing(6)
        self.gridLayout_4.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.frame = QFrame(self.centralWidget)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.horizontalLayout = QHBoxLayout(self.frame)
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout.setObjectName(u"horizontalLayout")

        self.gridLayout_4.addWidget(self.frame, 2, 0, 1, 1)

        self.tabWidget = QTabWidget(self.centralWidget)
        self.tabWidget.setObjectName(u"tabWidget")
        sizePolicy1 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy1)
        self.tabWidget.setToolTipDuration(-1)
        self.tab_2 = QWidget()
        self.tab_2.setObjectName(u"tab_2")
        self.horizontalLayout_2 = QHBoxLayout(self.tab_2)
        self.horizontalLayout_2.setSpacing(6)
        self.horizontalLayout_2.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setSpacing(6)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.label_13 = QLabel(self.tab_2)
        self.label_13.setObjectName(u"label_13")
        sizePolicy2 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.label_13.sizePolicy().hasHeightForWidth())
        self.label_13.setSizePolicy(sizePolicy2)
        font = QFont()
        font.setPointSize(12)
        font.setBold(True)
        font.setWeight(75)
        self.label_13.setFont(font)
        self.label_13.setAlignment(Qt.AlignCenter)

        self.verticalLayout_2.addWidget(self.label_13)

        self.label_14 = QLabel(self.tab_2)
        self.label_14.setObjectName(u"label_14")
        sizePolicy2.setHeightForWidth(self.label_14.sizePolicy().hasHeightForWidth())
        self.label_14.setSizePolicy(sizePolicy2)
        font1 = QFont()
        font1.setPointSize(10)
        font1.setBold(True)
        font1.setWeight(75)
        self.label_14.setFont(font1)

        self.verticalLayout_2.addWidget(self.label_14)

        self.lineEdit = QLineEdit(self.tab_2)
        self.lineEdit.setObjectName(u"lineEdit")
        sizePolicy3 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.lineEdit.sizePolicy().hasHeightForWidth())
        self.lineEdit.setSizePolicy(sizePolicy3)

        self.verticalLayout_2.addWidget(self.lineEdit)

        self.label_15 = QLabel(self.tab_2)
        self.label_15.setObjectName(u"label_15")
        sizePolicy2.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy2)
        self.label_15.setFont(font1)

        self.verticalLayout_2.addWidget(self.label_15)

        self.label_16 = QLabel(self.tab_2)
        self.label_16.setObjectName(u"label_16")
        sizePolicy4 = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy4)
        font2 = QFont()
        font2.setPointSize(10)
        self.label_16.setFont(font2)

        self.verticalLayout_2.addWidget(self.label_16)

        self.doubleSpinBox = QDoubleSpinBox(self.tab_2)
        self.doubleSpinBox.setObjectName(u"doubleSpinBox")
        sizePolicy5 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Maximum)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.doubleSpinBox.sizePolicy().hasHeightForWidth())
        self.doubleSpinBox.setSizePolicy(sizePolicy5)
        self.doubleSpinBox.setMinimumSize(QSize(149, 0))

        self.verticalLayout_2.addWidget(self.doubleSpinBox)

        self.label_17 = QLabel(self.tab_2)
        self.label_17.setObjectName(u"label_17")
        sizePolicy4.setHeightForWidth(self.label_17.sizePolicy().hasHeightForWidth())
        self.label_17.setSizePolicy(sizePolicy4)
        self.label_17.setFont(font2)

        self.verticalLayout_2.addWidget(self.label_17)

        self.doubleSpinBox_2 = QDoubleSpinBox(self.tab_2)
        self.doubleSpinBox_2.setObjectName(u"doubleSpinBox_2")
        sizePolicy6 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)
        sizePolicy6.setHeightForWidth(self.doubleSpinBox_2.sizePolicy().hasHeightForWidth())
        self.doubleSpinBox_2.setSizePolicy(sizePolicy6)
        self.doubleSpinBox_2.setMinimumSize(QSize(149, 0))

        self.verticalLayout_2.addWidget(self.doubleSpinBox_2)

        self.label_18 = QLabel(self.tab_2)
        self.label_18.setObjectName(u"label_18")
        sizePolicy4.setHeightForWidth(self.label_18.sizePolicy().hasHeightForWidth())
        self.label_18.setSizePolicy(sizePolicy4)

        self.verticalLayout_2.addWidget(self.label_18)

        self.doubleSpinBox_3 = QDoubleSpinBox(self.tab_2)
        self.doubleSpinBox_3.setObjectName(u"doubleSpinBox_3")
        sizePolicy6.setHeightForWidth(self.doubleSpinBox_3.sizePolicy().hasHeightForWidth())
        self.doubleSpinBox_3.setSizePolicy(sizePolicy6)
        self.doubleSpinBox_3.setMinimumSize(QSize(149, 0))

        self.verticalLayout_2.addWidget(self.doubleSpinBox_3)

        self.label_19 = QLabel(self.tab_2)
        self.label_19.setObjectName(u"label_19")
        self.label_19.setAlignment(Qt.AlignCenter)

        self.verticalLayout_2.addWidget(self.label_19)

        self.q_injection = QPushButton(self.tab_2)
        self.q_injection.setObjectName(u"q_injection")

        self.verticalLayout_2.addWidget(self.q_injection)


        self.horizontalLayout_2.addLayout(self.verticalLayout_2)

        self.groupBox = QGroupBox(self.tab_2)
        self.groupBox.setObjectName(u"groupBox")
        self.gridLayout = QGridLayout(self.groupBox)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setObjectName(u"gridLayout")
        self.label = QLabel(self.groupBox)
        self.label.setObjectName(u"label")

        self.gridLayout.addWidget(self.label, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_2 = QLabel(self.groupBox)
        self.label_2.setObjectName(u"label_2")

        self.gridLayout.addWidget(self.label_2, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_11 = QLabel(self.groupBox)
        self.label_11.setObjectName(u"label_11")
        sizePolicy2.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy2)
        self.label_11.setFont(font1)

        self.gridLayout.addWidget(self.label_11, 2, 0, 1, 2)

        self.q_alpha = QDoubleSpinBox(self.groupBox)
        self.q_alpha.setObjectName(u"q_alpha")
        self.q_alpha.setMinimum(-33.000000000000000)
        self.q_alpha.setMaximum(45.000000000000000)

        self.gridLayout.addWidget(self.q_alpha, 5, 0, 1, 1)

        self.q_x = QDoubleSpinBox(self.groupBox)
        self.q_x.setObjectName(u"q_x")
        self.q_x.setMaximum(49.000000000000000)

        self.gridLayout.addWidget(self.q_x, 3, 0, 1, 1)

        self.label_5 = QLabel(self.groupBox)
        self.label_5.setObjectName(u"label_5")
        sizePolicy7 = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.Preferred)
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy7)

        self.gridLayout.addWidget(self.label_5, 8, 1, 1, 1, Qt.AlignHCenter)

        self.q_y = QDoubleSpinBox(self.groupBox)
        self.q_y.setObjectName(u"q_y")
        self.q_y.setMaximum(49.000000000000000)

        self.gridLayout.addWidget(self.q_y, 4, 0, 1, 1)

        self.checkBox = QCheckBox(self.groupBox)
        self.checkBox.setObjectName(u"checkBox")

        self.gridLayout.addWidget(self.checkBox, 11, 0, 1, 1)

        self.label_4 = QLabel(self.groupBox)
        self.label_4.setObjectName(u"label_4")

        self.gridLayout.addWidget(self.label_4, 10, 1, 1, 1, Qt.AlignHCenter)

        self.label_12 = QLabel(self.groupBox)
        self.label_12.setObjectName(u"label_12")
        sizePolicy2.setHeightForWidth(self.label_12.sizePolicy().hasHeightForWidth())
        self.label_12.setSizePolicy(sizePolicy2)
        self.label_12.setFont(font1)

        self.gridLayout.addWidget(self.label_12, 7, 0, 1, 1)

        self.comboBox = QComboBox(self.groupBox)
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.setObjectName(u"comboBox")

        self.gridLayout.addWidget(self.comboBox, 10, 0, 1, 1)

        self.label_3 = QLabel(self.groupBox)
        self.label_3.setObjectName(u"label_3")

        self.gridLayout.addWidget(self.label_3, 5, 1, 1, 1, Qt.AlignHCenter)

        self.q_beta = QDoubleSpinBox(self.groupBox)
        self.q_beta.setObjectName(u"q_beta")
        self.q_beta.setMinimum(-50.000000000000000)
        self.q_beta.setMaximum(44.000000000000000)

        self.gridLayout.addWidget(self.q_beta, 8, 0, 1, 1)

        self.MainWingLabel = QLabel(self.groupBox)
        self.MainWingLabel.setObjectName(u"MainWingLabel")
        sizePolicy8 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Maximum)
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(0)
        sizePolicy8.setHeightForWidth(self.MainWingLabel.sizePolicy().hasHeightForWidth())
        self.MainWingLabel.setSizePolicy(sizePolicy8)
        self.MainWingLabel.setFont(font)
        self.MainWingLabel.setContextMenuPolicy(Qt.PreventContextMenu)
        self.MainWingLabel.setFrameShadow(QFrame.Sunken)
        self.MainWingLabel.setLineWidth(1)
        self.MainWingLabel.setAlignment(Qt.AlignCenter)

        self.gridLayout.addWidget(self.MainWingLabel, 0, 0, 1, 2)


        self.horizontalLayout_2.addWidget(self.groupBox)

        self.tabWidget_2 = QTabWidget(self.tab_2)
        self.tabWidget_2.setObjectName(u"tabWidget_2")
        sizePolicy9 = QSizePolicy(QSizePolicy.Ignored, QSizePolicy.Expanding)
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.tabWidget_2.sizePolicy().hasHeightForWidth())
        self.tabWidget_2.setSizePolicy(sizePolicy9)
        self.tab_7 = QWidget()
        self.tab_7.setObjectName(u"tab_7")
        self.tab_7.setMinimumSize(QSize(318, 0))
        self.groupBox_3 = QGroupBox(self.tab_7)
        self.groupBox_3.setObjectName(u"groupBox_3")
        self.groupBox_3.setGeometry(QRect(-7, -1, 331, 351))
        self.gridLayout_7 = QGridLayout(self.groupBox_3)
        self.gridLayout_7.setSpacing(6)
        self.gridLayout_7.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.label_32 = QLabel(self.groupBox_3)
        self.label_32.setObjectName(u"label_32")
        self.label_32.setFont(font1)

        self.gridLayout_7.addWidget(self.label_32, 12, 0, 1, 1)

        self.label_30 = QLabel(self.groupBox_3)
        self.label_30.setObjectName(u"label_30")
        sizePolicy2.setHeightForWidth(self.label_30.sizePolicy().hasHeightForWidth())
        self.label_30.setSizePolicy(sizePolicy2)
        self.label_30.setFont(font1)

        self.gridLayout_7.addWidget(self.label_30, 7, 0, 1, 1)

        self.checkBox_2 = QCheckBox(self.groupBox_3)
        self.checkBox_2.setObjectName(u"checkBox_2")

        self.gridLayout_7.addWidget(self.checkBox_2, 11, 0, 1, 1)

        self.q_alpha_3 = QDoubleSpinBox(self.groupBox_3)
        self.q_alpha_3.setObjectName(u"q_alpha_3")
        self.q_alpha_3.setMinimum(-33.000000000000000)
        self.q_alpha_3.setMaximum(45.000000000000000)

        self.gridLayout_7.addWidget(self.q_alpha_3, 5, 0, 1, 1)

        self.q_x_3 = QDoubleSpinBox(self.groupBox_3)
        self.q_x_3.setObjectName(u"q_x_3")
        self.q_x_3.setMaximum(49.000000000000000)

        self.gridLayout_7.addWidget(self.q_x_3, 3, 0, 1, 1)

        self.label_29 = QLabel(self.groupBox_3)
        self.label_29.setObjectName(u"label_29")

        self.gridLayout_7.addWidget(self.label_29, 10, 1, 1, 1, Qt.AlignHCenter)

        self.label_31 = QLabel(self.groupBox_3)
        self.label_31.setObjectName(u"label_31")

        self.gridLayout_7.addWidget(self.label_31, 5, 1, 1, 1, Qt.AlignHCenter)

        self.label_23 = QLabel(self.groupBox_3)
        self.label_23.setObjectName(u"label_23")

        self.gridLayout_7.addWidget(self.label_23, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_25 = QLabel(self.groupBox_3)
        self.label_25.setObjectName(u"label_25")
        sizePolicy7.setHeightForWidth(self.label_25.sizePolicy().hasHeightForWidth())
        self.label_25.setSizePolicy(sizePolicy7)
        self.label_25.setMinimumSize(QSize(0, 0))

        self.gridLayout_7.addWidget(self.label_25, 8, 1, 1, 1, Qt.AlignHCenter)

        self.q_beta_3 = QDoubleSpinBox(self.groupBox_3)
        self.q_beta_3.setObjectName(u"q_beta_3")
        self.q_beta_3.setMinimum(-50.000000000000000)
        self.q_beta_3.setMaximum(44.000000000000000)

        self.gridLayout_7.addWidget(self.q_beta_3, 8, 0, 1, 1)

        self.comboBox_2 = QComboBox(self.groupBox_3)
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.addItem("")
        self.comboBox_2.setObjectName(u"comboBox_2")

        self.gridLayout_7.addWidget(self.comboBox_2, 10, 0, 1, 1)

        self.doubleSpinBox_5 = QDoubleSpinBox(self.groupBox_3)
        self.doubleSpinBox_5.setObjectName(u"doubleSpinBox_5")

        self.gridLayout_7.addWidget(self.doubleSpinBox_5, 14, 0, 1, 1)

        self.label_33 = QLabel(self.groupBox_3)
        self.label_33.setObjectName(u"label_33")
        self.label_33.setAlignment(Qt.AlignCenter)

        self.gridLayout_7.addWidget(self.label_33, 13, 1, 1, 1)

        self.doubleSpinBox_4 = QDoubleSpinBox(self.groupBox_3)
        self.doubleSpinBox_4.setObjectName(u"doubleSpinBox_4")

        self.gridLayout_7.addWidget(self.doubleSpinBox_4, 13, 0, 1, 1)

        self.q_y_3 = QDoubleSpinBox(self.groupBox_3)
        self.q_y_3.setObjectName(u"q_y_3")
        self.q_y_3.setMaximum(49.000000000000000)

        self.gridLayout_7.addWidget(self.q_y_3, 4, 0, 1, 1)

        self.label_34 = QLabel(self.groupBox_3)
        self.label_34.setObjectName(u"label_34")

        self.gridLayout_7.addWidget(self.label_34, 14, 1, 1, 1)

        self.MainWingLabel_2 = QLabel(self.groupBox_3)
        self.MainWingLabel_2.setObjectName(u"MainWingLabel_2")
        sizePolicy8.setHeightForWidth(self.MainWingLabel_2.sizePolicy().hasHeightForWidth())
        self.MainWingLabel_2.setSizePolicy(sizePolicy8)
        self.MainWingLabel_2.setFont(font)
        self.MainWingLabel_2.setContextMenuPolicy(Qt.PreventContextMenu)
        self.MainWingLabel_2.setFrameShadow(QFrame.Sunken)
        self.MainWingLabel_2.setLineWidth(1)
        self.MainWingLabel_2.setAlignment(Qt.AlignCenter)

        self.gridLayout_7.addWidget(self.MainWingLabel_2, 0, 0, 1, 2)

        self.label_24 = QLabel(self.groupBox_3)
        self.label_24.setObjectName(u"label_24")
        sizePolicy2.setHeightForWidth(self.label_24.sizePolicy().hasHeightForWidth())
        self.label_24.setSizePolicy(sizePolicy2)
        self.label_24.setFont(font1)

        self.gridLayout_7.addWidget(self.label_24, 2, 0, 1, 2)

        self.label_20 = QLabel(self.groupBox_3)
        self.label_20.setObjectName(u"label_20")

        self.gridLayout_7.addWidget(self.label_20, 3, 1, 1, 1, Qt.AlignHCenter)

        self.doubleSpinBox_6 = QDoubleSpinBox(self.groupBox_3)
        self.doubleSpinBox_6.setObjectName(u"doubleSpinBox_6")

        self.gridLayout_7.addWidget(self.doubleSpinBox_6, 15, 0, 1, 1)

        self.label_35 = QLabel(self.groupBox_3)
        self.label_35.setObjectName(u"label_35")
        self.label_35.setAlignment(Qt.AlignCenter)

        self.gridLayout_7.addWidget(self.label_35, 15, 1, 1, 1)

        self.tabWidget_2.addTab(self.tab_7, "")
        self.tab_8 = QWidget()
        self.tab_8.setObjectName(u"tab_8")
        self.tabWidget_3 = QTabWidget(self.tab_8)
        self.tabWidget_3.setObjectName(u"tabWidget_3")
        self.tabWidget_3.setGeometry(QRect(0, -20, 324, 371))
        sizePolicy9.setHeightForWidth(self.tabWidget_3.sizePolicy().hasHeightForWidth())
        self.tabWidget_3.setSizePolicy(sizePolicy9)
        self.tab_13 = QWidget()
        self.tab_13.setObjectName(u"tab_13")
        sizePolicy7.setHeightForWidth(self.tab_13.sizePolicy().hasHeightForWidth())
        self.tab_13.setSizePolicy(sizePolicy7)
        self.tab_13.setMinimumSize(QSize(318, 0))
        self.groupBox_7 = QGroupBox(self.tab_13)
        self.groupBox_7.setObjectName(u"groupBox_7")
        self.groupBox_7.setGeometry(QRect(-7, -1, 331, 351))
        self.gridLayout_10 = QGridLayout(self.groupBox_7)
        self.gridLayout_10.setSpacing(6)
        self.gridLayout_10.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.label_58 = QLabel(self.groupBox_7)
        self.label_58.setObjectName(u"label_58")
        self.label_58.setFont(font1)

        self.gridLayout_10.addWidget(self.label_58, 12, 0, 1, 1)

        self.label_59 = QLabel(self.groupBox_7)
        self.label_59.setObjectName(u"label_59")
        sizePolicy2.setHeightForWidth(self.label_59.sizePolicy().hasHeightForWidth())
        self.label_59.setSizePolicy(sizePolicy2)
        self.label_59.setFont(font1)

        self.gridLayout_10.addWidget(self.label_59, 7, 0, 1, 1)

        self.checkBox_5 = QCheckBox(self.groupBox_7)
        self.checkBox_5.setObjectName(u"checkBox_5")

        self.gridLayout_10.addWidget(self.checkBox_5, 11, 0, 1, 1)

        self.q_alpha_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_alpha_6.setObjectName(u"q_alpha_6")
        self.q_alpha_6.setMinimum(-33.000000000000000)
        self.q_alpha_6.setMaximum(45.000000000000000)

        self.gridLayout_10.addWidget(self.q_alpha_6, 5, 0, 1, 1)

        self.q_x_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_x_6.setObjectName(u"q_x_6")
        self.q_x_6.setMaximum(49.000000000000000)

        self.gridLayout_10.addWidget(self.q_x_6, 3, 0, 1, 1)

        self.label_60 = QLabel(self.groupBox_7)
        self.label_60.setObjectName(u"label_60")

        self.gridLayout_10.addWidget(self.label_60, 10, 1, 1, 1, Qt.AlignHCenter)

        self.label_61 = QLabel(self.groupBox_7)
        self.label_61.setObjectName(u"label_61")

        self.gridLayout_10.addWidget(self.label_61, 5, 1, 1, 1, Qt.AlignHCenter)

        self.label_62 = QLabel(self.groupBox_7)
        self.label_62.setObjectName(u"label_62")

        self.gridLayout_10.addWidget(self.label_62, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_63 = QLabel(self.groupBox_7)
        self.label_63.setObjectName(u"label_63")
        sizePolicy7.setHeightForWidth(self.label_63.sizePolicy().hasHeightForWidth())
        self.label_63.setSizePolicy(sizePolicy7)
        self.label_63.setMinimumSize(QSize(0, 0))

        self.gridLayout_10.addWidget(self.label_63, 8, 1, 1, 1, Qt.AlignHCenter)

        self.q_beta_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_beta_6.setObjectName(u"q_beta_6")
        self.q_beta_6.setMinimum(-50.000000000000000)
        self.q_beta_6.setMaximum(44.000000000000000)

        self.gridLayout_10.addWidget(self.q_beta_6, 8, 0, 1, 1)

        self.comboBox_5 = QComboBox(self.groupBox_7)
        self.comboBox_5.addItem("")
        self.comboBox_5.addItem("")
        self.comboBox_5.addItem("")
        self.comboBox_5.setObjectName(u"comboBox_5")

        self.gridLayout_10.addWidget(self.comboBox_5, 10, 0, 1, 1)

        self.doubleSpinBox_13 = QDoubleSpinBox(self.groupBox_7)
        self.doubleSpinBox_13.setObjectName(u"doubleSpinBox_13")

        self.gridLayout_10.addWidget(self.doubleSpinBox_13, 14, 0, 1, 1)

        self.label_64 = QLabel(self.groupBox_7)
        self.label_64.setObjectName(u"label_64")
        self.label_64.setAlignment(Qt.AlignCenter)

        self.gridLayout_10.addWidget(self.label_64, 13, 1, 1, 1)

        self.doubleSpinBox_14 = QDoubleSpinBox(self.groupBox_7)
        self.doubleSpinBox_14.setObjectName(u"doubleSpinBox_14")

        self.gridLayout_10.addWidget(self.doubleSpinBox_14, 13, 0, 1, 1)

        self.q_y_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_y_6.setObjectName(u"q_y_6")
        self.q_y_6.setMaximum(49.000000000000000)

        self.gridLayout_10.addWidget(self.q_y_6, 4, 0, 1, 1)

        self.label_65 = QLabel(self.groupBox_7)
        self.label_65.setObjectName(u"label_65")

        self.gridLayout_10.addWidget(self.label_65, 14, 1, 1, 1)

        self.MainWingLabel_5 = QLabel(self.groupBox_7)
        self.MainWingLabel_5.setObjectName(u"MainWingLabel_5")
        sizePolicy8.setHeightForWidth(self.MainWingLabel_5.sizePolicy().hasHeightForWidth())
        self.MainWingLabel_5.setSizePolicy(sizePolicy8)
        self.MainWingLabel_5.setFont(font)
        self.MainWingLabel_5.setContextMenuPolicy(Qt.PreventContextMenu)
        self.MainWingLabel_5.setFrameShadow(QFrame.Sunken)
        self.MainWingLabel_5.setLineWidth(1)
        self.MainWingLabel_5.setAlignment(Qt.AlignCenter)

        self.gridLayout_10.addWidget(self.MainWingLabel_5, 0, 0, 1, 2)

        self.label_66 = QLabel(self.groupBox_7)
        self.label_66.setObjectName(u"label_66")
        sizePolicy2.setHeightForWidth(self.label_66.sizePolicy().hasHeightForWidth())
        self.label_66.setSizePolicy(sizePolicy2)
        self.label_66.setFont(font1)

        self.gridLayout_10.addWidget(self.label_66, 2, 0, 1, 2)

        self.label_67 = QLabel(self.groupBox_7)
        self.label_67.setObjectName(u"label_67")

        self.gridLayout_10.addWidget(self.label_67, 3, 1, 1, 1, Qt.AlignHCenter)

        self.doubleSpinBox_15 = QDoubleSpinBox(self.groupBox_7)
        self.doubleSpinBox_15.setObjectName(u"doubleSpinBox_15")

        self.gridLayout_10.addWidget(self.doubleSpinBox_15, 15, 0, 1, 1)

        self.label_68 = QLabel(self.groupBox_7)
        self.label_68.setObjectName(u"label_68")
        self.label_68.setAlignment(Qt.AlignCenter)

        self.gridLayout_10.addWidget(self.label_68, 15, 1, 1, 1)

        self.tabWidget_3.addTab(self.tab_13, "")
        self.tab_14 = QWidget()
        self.tab_14.setObjectName(u"tab_14")
        self.tabWidget_3.addTab(self.tab_14, "")
        self.tabWidget_2.addTab(self.tab_8, "")

        self.horizontalLayout_2.addWidget(self.tabWidget_2)

        self.tabWidget.addTab(self.tab_2, "")
        self.tab = QWidget()
        self.tab.setObjectName(u"tab")
        self.horizontalLayout_3 = QHBoxLayout(self.tab)
        self.horizontalLayout_3.setSpacing(6)
        self.horizontalLayout_3.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.groupBox_2 = QGroupBox(self.tab)
        self.groupBox_2.setObjectName(u"groupBox_2")
        self.gridLayout_2 = QGridLayout(self.groupBox_2)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.label_10 = QLabel(self.groupBox_2)
        self.label_10.setObjectName(u"label_10")

        self.gridLayout_2.addWidget(self.label_10, 4, 1, 1, 1, Qt.AlignHCenter)

        self.ee_phi = QDoubleSpinBox(self.groupBox_2)
        self.ee_phi.setObjectName(u"ee_phi")
        self.ee_phi.setMinimum(57.000000000000000)
        self.ee_phi.setMaximum(135.000000000000000)

        self.gridLayout_2.addWidget(self.ee_phi, 3, 0, 1, 1)

        self.ee_z = QDoubleSpinBox(self.groupBox_2)
        self.ee_z.setObjectName(u"ee_z")
        self.ee_z.setMinimum(-300.000000000000000)
        self.ee_z.setMaximum(300.000000000000000)

        self.gridLayout_2.addWidget(self.ee_z, 2, 0, 1, 1)

        self.ee_start = QPushButton(self.groupBox_2)
        self.ee_start.setObjectName(u"ee_start")

        self.gridLayout_2.addWidget(self.ee_start, 5, 0, 1, 1)

        self.ee_x = QDoubleSpinBox(self.groupBox_2)
        self.ee_x.setObjectName(u"ee_x")
        self.ee_x.setMinimum(-300.000000000000000)
        self.ee_x.setMaximum(300.000000000000000)

        self.gridLayout_2.addWidget(self.ee_x, 0, 0, 1, 1)

        self.ee_y = QDoubleSpinBox(self.groupBox_2)
        self.ee_y.setObjectName(u"ee_y")
        self.ee_y.setMinimum(-300.000000000000000)
        self.ee_y.setMaximum(300.000000000000000)

        self.gridLayout_2.addWidget(self.ee_y, 1, 0, 1, 1)

        self.ee_theta = QDoubleSpinBox(self.groupBox_2)
        self.ee_theta.setObjectName(u"ee_theta")
        self.ee_theta.setMinimum(-50.000000000000000)
        self.ee_theta.setMaximum(44.000000000000000)

        self.gridLayout_2.addWidget(self.ee_theta, 4, 0, 1, 1)

        self.ee_homing = QPushButton(self.groupBox_2)
        self.ee_homing.setObjectName(u"ee_homing")

        self.gridLayout_2.addWidget(self.ee_homing, 5, 2, 1, 1)

        self.label_8 = QLabel(self.groupBox_2)
        self.label_8.setObjectName(u"label_8")

        self.gridLayout_2.addWidget(self.label_8, 2, 1, 1, 1, Qt.AlignHCenter)

        self.label_6 = QLabel(self.groupBox_2)
        self.label_6.setObjectName(u"label_6")

        self.gridLayout_2.addWidget(self.label_6, 0, 1, 1, 1, Qt.AlignHCenter)

        self.label_9 = QLabel(self.groupBox_2)
        self.label_9.setObjectName(u"label_9")

        self.gridLayout_2.addWidget(self.label_9, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_7 = QLabel(self.groupBox_2)
        self.label_7.setObjectName(u"label_7")

        self.gridLayout_2.addWidget(self.label_7, 1, 1, 1, 1, Qt.AlignHCenter)

        self.ee_injection = QPushButton(self.groupBox_2)
        self.ee_injection.setObjectName(u"ee_injection")

        self.gridLayout_2.addWidget(self.ee_injection, 5, 1, 1, 1)


        self.horizontalLayout_3.addWidget(self.groupBox_2)

        self.tabWidget.addTab(self.tab, "")
        self.tab_4 = QWidget()
        self.tab_4.setObjectName(u"tab_4")
        self.tabWidget.addTab(self.tab_4, "")
        self.tab_5 = QWidget()
        self.tab_5.setObjectName(u"tab_5")
        self.tabWidget.addTab(self.tab_5, "")
        self.tab_6 = QWidget()
        self.tab_6.setObjectName(u"tab_6")
        self.tabWidget.addTab(self.tab_6, "")
        self.tab_3 = QWidget()
        self.tab_3.setObjectName(u"tab_3")
        self.groupBox_4 = QGroupBox(self.tab_3)
        self.groupBox_4.setObjectName(u"groupBox_4")
        self.groupBox_4.setGeometry(QRect(0, 0, 331, 235))
        self.gridLayout_6 = QGridLayout(self.groupBox_4)
        self.gridLayout_6.setSpacing(6)
        self.gridLayout_6.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.q_start_2 = QPushButton(self.groupBox_4)
        self.q_start_2.setObjectName(u"q_start_2")

        self.gridLayout_6.addWidget(self.q_start_2, 5, 0, 1, 1)

        self.q_y_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_y_2.setObjectName(u"q_y_2")
        self.q_y_2.setMaximum(49.000000000000000)

        self.gridLayout_6.addWidget(self.q_y_2, 1, 0, 1, 1)

        self.q_homing_2 = QPushButton(self.groupBox_4)
        self.q_homing_2.setObjectName(u"q_homing_2")

        self.gridLayout_6.addWidget(self.q_homing_2, 5, 2, 1, 1)

        self.q_beta_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_beta_2.setObjectName(u"q_beta_2")
        self.q_beta_2.setMinimum(-50.000000000000000)
        self.q_beta_2.setMaximum(44.000000000000000)

        self.gridLayout_6.addWidget(self.q_beta_2, 3, 0, 1, 1)

        self.q_alpha_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_alpha_2.setObjectName(u"q_alpha_2")
        self.q_alpha_2.setMinimum(-33.000000000000000)
        self.q_alpha_2.setMaximum(45.000000000000000)

        self.gridLayout_6.addWidget(self.q_alpha_2, 2, 0, 1, 1)

        self.q_x_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_x_2.setObjectName(u"q_x_2")
        self.q_x_2.setMaximum(49.000000000000000)

        self.gridLayout_6.addWidget(self.q_x_2, 0, 0, 1, 1)

        self.q_d_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_d_2.setObjectName(u"q_d_2")
        self.q_d_2.setMinimum(-45.000000000000000)
        self.q_d_2.setMaximum(0.000000000000000)

        self.gridLayout_6.addWidget(self.q_d_2, 4, 0, 1, 1)

        self.q_injection_2 = QPushButton(self.groupBox_4)
        self.q_injection_2.setObjectName(u"q_injection_2")

        self.gridLayout_6.addWidget(self.q_injection_2, 5, 1, 1, 1)

        self.label_21 = QLabel(self.groupBox_4)
        self.label_21.setObjectName(u"label_21")

        self.gridLayout_6.addWidget(self.label_21, 0, 1, 1, 1, Qt.AlignHCenter)

        self.label_22 = QLabel(self.groupBox_4)
        self.label_22.setObjectName(u"label_22")

        self.gridLayout_6.addWidget(self.label_22, 1, 1, 1, 1, Qt.AlignHCenter)

        self.label_26 = QLabel(self.groupBox_4)
        self.label_26.setObjectName(u"label_26")

        self.gridLayout_6.addWidget(self.label_26, 2, 1, 1, 1, Qt.AlignHCenter)

        self.label_27 = QLabel(self.groupBox_4)
        self.label_27.setObjectName(u"label_27")

        self.gridLayout_6.addWidget(self.label_27, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_28 = QLabel(self.groupBox_4)
        self.label_28.setObjectName(u"label_28")

        self.gridLayout_6.addWidget(self.label_28, 4, 1, 1, 1, Qt.AlignHCenter)

        self.tabWidget.addTab(self.tab_3, "")

        self.gridLayout_4.addWidget(self.tabWidget, 2, 1, 1, 1)

        self.frame_3 = QFrame(self.centralWidget)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setAutoFillBackground(False)
        self.frame_3.setFrameShape(QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QFrame.Raised)
        self.gridLayout_5 = QGridLayout(self.frame_3)
        self.gridLayout_5.setSpacing(6)
        self.gridLayout_5.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_5.setObjectName(u"gridLayout_5")
        self.Logo = QLabel(self.frame_3)
        self.Logo.setObjectName(u"Logo")
        self.Logo.setAlignment(Qt.AlignCenter)

        self.gridLayout_5.addWidget(self.Logo, 0, 0, 1, 1)


        self.gridLayout_4.addWidget(self.frame_3, 1, 0, 1, 2)

        self.progressBar = QProgressBar(self.centralWidget)
        self.progressBar.setObjectName(u"progressBar")
        self.progressBar.setMinimum(0)
        self.progressBar.setValue(24)

        self.gridLayout_4.addWidget(self.progressBar, 5, 0, 1, 2)

        self.frame_2 = QFrame(self.centralWidget)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Raised)
        self.verticalLayout = QVBoxLayout(self.frame_2)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.q_start = QPushButton(self.frame_2)
        self.q_start.setObjectName(u"q_start")
        self.q_start.setContextMenuPolicy(Qt.NoContextMenu)

        self.verticalLayout.addWidget(self.q_start)

        self.q_homing = QPushButton(self.frame_2)
        self.q_homing.setObjectName(u"q_homing")

        self.verticalLayout.addWidget(self.q_homing)

        self.terminalOutpotLogo = QLabel(self.frame_2)
        self.terminalOutpotLogo.setObjectName(u"terminalOutpotLogo")

        self.verticalLayout.addWidget(self.terminalOutpotLogo)

        self.view_logging = QListView(self.frame_2)
        self.view_logging.setObjectName(u"view_logging")
        sizePolicy10 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        sizePolicy10.setHorizontalStretch(0)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.view_logging.sizePolicy().hasHeightForWidth())
        self.view_logging.setSizePolicy(sizePolicy10)

        self.verticalLayout.addWidget(self.view_logging)


        self.gridLayout_4.addWidget(self.frame_2, 4, 0, 1, 2)

        MainWindowDesign.setCentralWidget(self.centralWidget)
        self.menuBar = QMenuBar(MainWindowDesign)
        self.menuBar.setObjectName(u"menuBar")
        self.menuBar.setGeometry(QRect(0, 0, 1134, 21))
        self.menuHelp = QMenu(self.menuBar)
        self.menuHelp.setObjectName(u"menuHelp")
        self.menuFile = QMenu(self.menuBar)
        self.menuFile.setObjectName(u"menuFile")
        self.menuAircraft = QMenu(self.menuBar)
        self.menuAircraft.setObjectName(u"menuAircraft")
        MainWindowDesign.setMenuBar(self.menuBar)
        self.statusBar = QStatusBar(MainWindowDesign)
        self.statusBar.setObjectName(u"statusBar")
        MainWindowDesign.setStatusBar(self.statusBar)

        self.menuBar.addAction(self.menuFile.menuAction())
        self.menuBar.addAction(self.menuAircraft.menuAction())
        self.menuBar.addAction(self.menuHelp.menuAction())
        self.menuHelp.addAction(self.actionAbout)
        self.menuHelp.addAction(self.actionManual)
        self.menuHelp.addSeparator()
        self.menuHelp.addAction(self.actionQuit)
        self.menuFile.addAction(self.actionImport_Data)
        self.menuFile.addAction(self.actionExport_Data)
        self.menuAircraft.addAction(self.actionPreset_Aircraft_List)

        self.retranslateUi(MainWindowDesign)

        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(1)
        self.tabWidget_3.setCurrentIndex(0)


        QMetaObject.connectSlotsByName(MainWindowDesign)
    # setupUi

    def retranslateUi(self, MainWindowDesign):
        MainWindowDesign.setWindowTitle(QCoreApplication.translate("MainWindowDesign", u"MainWindow", None))
        self.actionAbout.setText(QCoreApplication.translate("MainWindowDesign", u"About", None))
        self.actionManual.setText(QCoreApplication.translate("MainWindowDesign", u"Manual", None))
        self.actionQuit.setText(QCoreApplication.translate("MainWindowDesign", u"Quit", None))
        self.actionImport_Data.setText(QCoreApplication.translate("MainWindowDesign", u"Import Data", None))
        self.actionExport_Data.setText(QCoreApplication.translate("MainWindowDesign", u"Export Data", None))
        self.actionPreset_Aircraft_List.setText(QCoreApplication.translate("MainWindowDesign", u"Preset Aircraft List", None))
        self.label_13.setText(QCoreApplication.translate("MainWindowDesign", u"Aircraft", None))
        self.label_14.setText(QCoreApplication.translate("MainWindowDesign", u"Name", None))
        self.label_15.setText(QCoreApplication.translate("MainWindowDesign", u"CofG Position", None))
        self.label_16.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_17.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_18.setText(QCoreApplication.translate("MainWindowDesign", u"Z", None))
        self.label_19.setText(QCoreApplication.translate("MainWindowDesign", u"Other Options go here", None))
        self.q_injection.setText(QCoreApplication.translate("MainWindowDesign", u"Save Aircraft", None))
        self.groupBox.setTitle("")
        self.label.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_2.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_11.setText(QCoreApplication.translate("MainWindowDesign", u"Leading Edge Location", None))
        self.q_x.setSuffix("")
        self.label_5.setText(QCoreApplication.translate("MainWindowDesign", u"Number of Panels", None))
        self.checkBox.setText(QCoreApplication.translate("MainWindowDesign", u"Symmetric?", None))
        self.label_4.setText(QCoreApplication.translate("MainWindowDesign", u"Panel Spacing Option", None))
        self.label_12.setText(QCoreApplication.translate("MainWindowDesign", u"Panels", None))
        self.comboBox.setItemText(0, QCoreApplication.translate("MainWindowDesign", u"Uniform", None))
        self.comboBox.setItemText(1, QCoreApplication.translate("MainWindowDesign", u"Cosine", None))
        self.comboBox.setItemText(2, QCoreApplication.translate("MainWindowDesign", u"Reverse Cosine", None))

        self.label_3.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.MainWingLabel.setText(QCoreApplication.translate("MainWindowDesign", u"Main Wing", None))
        self.groupBox_3.setTitle("")
        self.label_32.setText(QCoreApplication.translate("MainWindowDesign", u"Chord and Control Surface", None))
        self.label_30.setText(QCoreApplication.translate("MainWindowDesign", u"Panels", None))
        self.checkBox_2.setText(QCoreApplication.translate("MainWindowDesign", u"Symmetric?", None))
        self.q_x_3.setSuffix("")
        self.label_29.setText(QCoreApplication.translate("MainWindowDesign", u"Panel Spacing Option", None))
        self.label_31.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.label_23.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_25.setText(QCoreApplication.translate("MainWindowDesign", u"Number of Panels", None))
        self.comboBox_2.setItemText(0, QCoreApplication.translate("MainWindowDesign", u"Uniform", None))
        self.comboBox_2.setItemText(1, QCoreApplication.translate("MainWindowDesign", u"Cosine", None))
        self.comboBox_2.setItemText(2, QCoreApplication.translate("MainWindowDesign", u"Reverse Cosine", None))

        self.label_33.setText(QCoreApplication.translate("MainWindowDesign", u"Chord", None))
        self.label_34.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Hinge Point", None))
        self.MainWingLabel_2.setText(QCoreApplication.translate("MainWindowDesign", u"Root Cross Section", None))
        self.label_24.setText(QCoreApplication.translate("MainWindowDesign", u"Leading Edge Location", None))
        self.label_20.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_35.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Deflection", None))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_7), QCoreApplication.translate("MainWindowDesign", u"Wing Root Cross Section", None))
        self.groupBox_7.setTitle("")
        self.label_58.setText(QCoreApplication.translate("MainWindowDesign", u"Chord and Control Surface", None))
        self.label_59.setText(QCoreApplication.translate("MainWindowDesign", u"Panels", None))
        self.checkBox_5.setText(QCoreApplication.translate("MainWindowDesign", u"Symmetric?", None))
        self.q_x_6.setSuffix("")
        self.label_60.setText(QCoreApplication.translate("MainWindowDesign", u"Panel Spacing Option", None))
        self.label_61.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.label_62.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_63.setText(QCoreApplication.translate("MainWindowDesign", u"Number of Panels", None))
        self.comboBox_5.setItemText(0, QCoreApplication.translate("MainWindowDesign", u"Uniform", None))
        self.comboBox_5.setItemText(1, QCoreApplication.translate("MainWindowDesign", u"Cosine", None))
        self.comboBox_5.setItemText(2, QCoreApplication.translate("MainWindowDesign", u"Reverse Cosine", None))

        self.label_64.setText(QCoreApplication.translate("MainWindowDesign", u"Chord", None))
        self.label_65.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Hinge Point", None))
        self.MainWingLabel_5.setText(QCoreApplication.translate("MainWindowDesign", u"Root Cross Section", None))
        self.label_66.setText(QCoreApplication.translate("MainWindowDesign", u"Leading Edge Location", None))
        self.label_67.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_68.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Deflection", None))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_13), QCoreApplication.translate("MainWindowDesign", u"Wing Root Cross Section", None))
        self.tabWidget_3.setTabText(self.tabWidget_3.indexOf(self.tab_14), QCoreApplication.translate("MainWindowDesign", u"Wing Tip Cross Section", None))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_8), QCoreApplication.translate("MainWindowDesign", u"Wing Tip Cross Section", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), QCoreApplication.translate("MainWindowDesign", u"Aircraft Parameters", None))
        self.groupBox_2.setTitle(QCoreApplication.translate("MainWindowDesign", u"Controller", None))
        self.label_10.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03b8</p></body></html>", None))
        self.ee_start.setText(QCoreApplication.translate("MainWindowDesign", u"Move", None))
        self.ee_homing.setText(QCoreApplication.translate("MainWindowDesign", u"Homing", None))
        self.label_8.setText(QCoreApplication.translate("MainWindowDesign", u"z", None))
        self.label_6.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_9.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03c6</p></body></html>", None))
        self.label_7.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.ee_injection.setText(QCoreApplication.translate("MainWindowDesign", u"Injection", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), QCoreApplication.translate("MainWindowDesign", u"Model ", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4), QCoreApplication.translate("MainWindowDesign", u"Materials", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5), QCoreApplication.translate("MainWindowDesign", u"Boundary Conditions", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6), QCoreApplication.translate("MainWindowDesign", u"Solution", None))
        self.groupBox_4.setTitle("")
        self.q_start_2.setText(QCoreApplication.translate("MainWindowDesign", u"Move", None))
        self.q_homing_2.setText(QCoreApplication.translate("MainWindowDesign", u"Homing", None))
        self.q_x_2.setSuffix("")
        self.q_injection_2.setText(QCoreApplication.translate("MainWindowDesign", u"Injection", None))
        self.label_21.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_22.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_26.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.label_27.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03b2</p></body></html>", None))
        self.label_28.setText(QCoreApplication.translate("MainWindowDesign", u"d", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), QCoreApplication.translate("MainWindowDesign", u"Visualiser Settings", None))
        self.Logo.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p><img src=\"docs/Black_Text_Logo.png\"/></p></body></html>", None))
        self.q_start.setText(QCoreApplication.translate("MainWindowDesign", u"Generate Results", None))
        self.q_homing.setText(QCoreApplication.translate("MainWindowDesign", u"Plot Visualisation", None))
        self.terminalOutpotLogo.setText(QCoreApplication.translate("MainWindowDesign", u"Terminal Output", None))
        self.menuHelp.setTitle(QCoreApplication.translate("MainWindowDesign", u"Help", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindowDesign", u"File", None))
        self.menuAircraft.setTitle(QCoreApplication.translate("MainWindowDesign", u"Aircraft", None))
    # retranslateUi

