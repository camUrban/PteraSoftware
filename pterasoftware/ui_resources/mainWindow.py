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
        MainWindowDesign.resize(1145, 880)
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
        self.actionExample_1 = QAction(MainWindowDesign)
        self.actionExample_1.setObjectName(u"actionExample_1")
        self.actionExample_2 = QAction(MainWindowDesign)
        self.actionExample_2.setObjectName(u"actionExample_2")
        self.actionExample_3 = QAction(MainWindowDesign)
        self.actionExample_3.setObjectName(u"actionExample_3")
        self.actionExample_4 = QAction(MainWindowDesign)
        self.actionExample_4.setObjectName(u"actionExample_4")
        self.actionExample_5 = QAction(MainWindowDesign)
        self.actionExample_5.setObjectName(u"actionExample_5")
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
        self.gridLayout_3 = QGridLayout(self.tab_7)
        self.gridLayout_3.setSpacing(6)
        self.gridLayout_3.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.groupBox_8 = QGroupBox(self.tab_7)
        self.groupBox_8.setObjectName(u"groupBox_8")
        self.gridLayout_11 = QGridLayout(self.groupBox_8)
        self.gridLayout_11.setSpacing(6)
        self.gridLayout_11.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_11.setObjectName(u"gridLayout_11")
        self.label_69 = QLabel(self.groupBox_8)
        self.label_69.setObjectName(u"label_69")
        self.label_69.setFont(font1)

        self.gridLayout_11.addWidget(self.label_69, 12, 0, 1, 1)

        self.label_70 = QLabel(self.groupBox_8)
        self.label_70.setObjectName(u"label_70")
        sizePolicy2.setHeightForWidth(self.label_70.sizePolicy().hasHeightForWidth())
        self.label_70.setSizePolicy(sizePolicy2)
        self.label_70.setFont(font1)

        self.gridLayout_11.addWidget(self.label_70, 7, 0, 1, 1)

        self.checkBox_6 = QCheckBox(self.groupBox_8)
        self.checkBox_6.setObjectName(u"checkBox_6")

        self.gridLayout_11.addWidget(self.checkBox_6, 11, 0, 1, 1)

        self.q_alpha_7 = QDoubleSpinBox(self.groupBox_8)
        self.q_alpha_7.setObjectName(u"q_alpha_7")
        self.q_alpha_7.setMinimum(-33.000000000000000)
        self.q_alpha_7.setMaximum(45.000000000000000)

        self.gridLayout_11.addWidget(self.q_alpha_7, 5, 0, 1, 1)

        self.q_x_7 = QDoubleSpinBox(self.groupBox_8)
        self.q_x_7.setObjectName(u"q_x_7")
        self.q_x_7.setMaximum(49.000000000000000)

        self.gridLayout_11.addWidget(self.q_x_7, 3, 0, 1, 1)

        self.label_71 = QLabel(self.groupBox_8)
        self.label_71.setObjectName(u"label_71")

        self.gridLayout_11.addWidget(self.label_71, 10, 1, 1, 1, Qt.AlignHCenter)

        self.label_72 = QLabel(self.groupBox_8)
        self.label_72.setObjectName(u"label_72")

        self.gridLayout_11.addWidget(self.label_72, 5, 1, 1, 1, Qt.AlignHCenter)

        self.label_73 = QLabel(self.groupBox_8)
        self.label_73.setObjectName(u"label_73")

        self.gridLayout_11.addWidget(self.label_73, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_74 = QLabel(self.groupBox_8)
        self.label_74.setObjectName(u"label_74")
        sizePolicy7.setHeightForWidth(self.label_74.sizePolicy().hasHeightForWidth())
        self.label_74.setSizePolicy(sizePolicy7)
        self.label_74.setMinimumSize(QSize(0, 0))

        self.gridLayout_11.addWidget(self.label_74, 8, 1, 1, 1, Qt.AlignHCenter)

        self.q_beta_7 = QDoubleSpinBox(self.groupBox_8)
        self.q_beta_7.setObjectName(u"q_beta_7")
        self.q_beta_7.setMinimum(-50.000000000000000)
        self.q_beta_7.setMaximum(44.000000000000000)

        self.gridLayout_11.addWidget(self.q_beta_7, 8, 0, 1, 1)

        self.comboBox_6 = QComboBox(self.groupBox_8)
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.addItem("")
        self.comboBox_6.setObjectName(u"comboBox_6")

        self.gridLayout_11.addWidget(self.comboBox_6, 10, 0, 1, 1)

        self.doubleSpinBox_16 = QDoubleSpinBox(self.groupBox_8)
        self.doubleSpinBox_16.setObjectName(u"doubleSpinBox_16")

        self.gridLayout_11.addWidget(self.doubleSpinBox_16, 14, 0, 1, 1)

        self.label_75 = QLabel(self.groupBox_8)
        self.label_75.setObjectName(u"label_75")
        self.label_75.setAlignment(Qt.AlignCenter)

        self.gridLayout_11.addWidget(self.label_75, 13, 1, 1, 1)

        self.doubleSpinBox_17 = QDoubleSpinBox(self.groupBox_8)
        self.doubleSpinBox_17.setObjectName(u"doubleSpinBox_17")

        self.gridLayout_11.addWidget(self.doubleSpinBox_17, 13, 0, 1, 1)

        self.q_y_7 = QDoubleSpinBox(self.groupBox_8)
        self.q_y_7.setObjectName(u"q_y_7")
        self.q_y_7.setMaximum(49.000000000000000)

        self.gridLayout_11.addWidget(self.q_y_7, 4, 0, 1, 1)

        self.label_76 = QLabel(self.groupBox_8)
        self.label_76.setObjectName(u"label_76")

        self.gridLayout_11.addWidget(self.label_76, 14, 1, 1, 1)

        self.MainWingLabel_6 = QLabel(self.groupBox_8)
        self.MainWingLabel_6.setObjectName(u"MainWingLabel_6")
        sizePolicy8.setHeightForWidth(self.MainWingLabel_6.sizePolicy().hasHeightForWidth())
        self.MainWingLabel_6.setSizePolicy(sizePolicy8)
        self.MainWingLabel_6.setFont(font)
        self.MainWingLabel_6.setContextMenuPolicy(Qt.PreventContextMenu)
        self.MainWingLabel_6.setFrameShadow(QFrame.Sunken)
        self.MainWingLabel_6.setLineWidth(1)
        self.MainWingLabel_6.setAlignment(Qt.AlignCenter)

        self.gridLayout_11.addWidget(self.MainWingLabel_6, 0, 0, 1, 2)

        self.label_77 = QLabel(self.groupBox_8)
        self.label_77.setObjectName(u"label_77")
        sizePolicy2.setHeightForWidth(self.label_77.sizePolicy().hasHeightForWidth())
        self.label_77.setSizePolicy(sizePolicy2)
        self.label_77.setFont(font1)

        self.gridLayout_11.addWidget(self.label_77, 2, 0, 1, 2)

        self.label_78 = QLabel(self.groupBox_8)
        self.label_78.setObjectName(u"label_78")

        self.gridLayout_11.addWidget(self.label_78, 3, 1, 1, 1, Qt.AlignHCenter)

        self.doubleSpinBox_18 = QDoubleSpinBox(self.groupBox_8)
        self.doubleSpinBox_18.setObjectName(u"doubleSpinBox_18")

        self.gridLayout_11.addWidget(self.doubleSpinBox_18, 15, 0, 1, 1)

        self.label_79 = QLabel(self.groupBox_8)
        self.label_79.setObjectName(u"label_79")
        self.label_79.setAlignment(Qt.AlignCenter)

        self.gridLayout_11.addWidget(self.label_79, 15, 1, 1, 1)


        self.gridLayout_3.addWidget(self.groupBox_8, 0, 0, 1, 1)

        self.tabWidget_2.addTab(self.tab_7, "")
        self.tab_8 = QWidget()
        self.tab_8.setObjectName(u"tab_8")
        self.gridLayout_7 = QGridLayout(self.tab_8)
        self.gridLayout_7.setSpacing(6)
        self.gridLayout_7.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.groupBox_7 = QGroupBox(self.tab_8)
        self.groupBox_7.setObjectName(u"groupBox_7")
        self.gridLayout_10 = QGridLayout(self.groupBox_7)
        self.gridLayout_10.setSpacing(6)
        self.gridLayout_10.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.doubleSpinBox_14 = QDoubleSpinBox(self.groupBox_7)
        self.doubleSpinBox_14.setObjectName(u"doubleSpinBox_14")

        self.gridLayout_10.addWidget(self.doubleSpinBox_14, 13, 0, 1, 1)

        self.q_y_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_y_6.setObjectName(u"q_y_6")
        self.q_y_6.setMaximum(49.000000000000000)

        self.gridLayout_10.addWidget(self.q_y_6, 4, 0, 1, 1)

        self.label_62 = QLabel(self.groupBox_7)
        self.label_62.setObjectName(u"label_62")

        self.gridLayout_10.addWidget(self.label_62, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_60 = QLabel(self.groupBox_7)
        self.label_60.setObjectName(u"label_60")

        self.gridLayout_10.addWidget(self.label_60, 10, 1, 1, 1, Qt.AlignHCenter)

        self.q_x_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_x_6.setObjectName(u"q_x_6")
        self.q_x_6.setMaximum(49.000000000000000)

        self.gridLayout_10.addWidget(self.q_x_6, 3, 0, 1, 1)

        self.label_63 = QLabel(self.groupBox_7)
        self.label_63.setObjectName(u"label_63")
        sizePolicy7.setHeightForWidth(self.label_63.sizePolicy().hasHeightForWidth())
        self.label_63.setSizePolicy(sizePolicy7)
        self.label_63.setMinimumSize(QSize(0, 0))

        self.gridLayout_10.addWidget(self.label_63, 8, 1, 1, 1, Qt.AlignHCenter)

        self.comboBox_5 = QComboBox(self.groupBox_7)
        self.comboBox_5.addItem("")
        self.comboBox_5.addItem("")
        self.comboBox_5.addItem("")
        self.comboBox_5.setObjectName(u"comboBox_5")

        self.gridLayout_10.addWidget(self.comboBox_5, 10, 0, 1, 1)

        self.label_58 = QLabel(self.groupBox_7)
        self.label_58.setObjectName(u"label_58")
        self.label_58.setFont(font1)

        self.gridLayout_10.addWidget(self.label_58, 12, 0, 1, 1)

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

        self.doubleSpinBox_13 = QDoubleSpinBox(self.groupBox_7)
        self.doubleSpinBox_13.setObjectName(u"doubleSpinBox_13")

        self.gridLayout_10.addWidget(self.doubleSpinBox_13, 14, 0, 1, 1)

        self.label_65 = QLabel(self.groupBox_7)
        self.label_65.setObjectName(u"label_65")

        self.gridLayout_10.addWidget(self.label_65, 14, 1, 1, 1)

        self.q_beta_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_beta_6.setObjectName(u"q_beta_6")
        self.q_beta_6.setMinimum(-50.000000000000000)
        self.q_beta_6.setMaximum(44.000000000000000)

        self.gridLayout_10.addWidget(self.q_beta_6, 8, 0, 1, 1)

        self.label_66 = QLabel(self.groupBox_7)
        self.label_66.setObjectName(u"label_66")
        sizePolicy2.setHeightForWidth(self.label_66.sizePolicy().hasHeightForWidth())
        self.label_66.setSizePolicy(sizePolicy2)
        self.label_66.setFont(font1)

        self.gridLayout_10.addWidget(self.label_66, 2, 0, 1, 2)

        self.checkBox_5 = QCheckBox(self.groupBox_7)
        self.checkBox_5.setObjectName(u"checkBox_5")

        self.gridLayout_10.addWidget(self.checkBox_5, 11, 0, 1, 1)

        self.label_67 = QLabel(self.groupBox_7)
        self.label_67.setObjectName(u"label_67")

        self.gridLayout_10.addWidget(self.label_67, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_64 = QLabel(self.groupBox_7)
        self.label_64.setObjectName(u"label_64")
        self.label_64.setAlignment(Qt.AlignCenter)

        self.gridLayout_10.addWidget(self.label_64, 13, 1, 1, 1)

        self.label_59 = QLabel(self.groupBox_7)
        self.label_59.setObjectName(u"label_59")
        sizePolicy2.setHeightForWidth(self.label_59.sizePolicy().hasHeightForWidth())
        self.label_59.setSizePolicy(sizePolicy2)
        self.label_59.setFont(font1)

        self.gridLayout_10.addWidget(self.label_59, 7, 0, 1, 1)

        self.label_61 = QLabel(self.groupBox_7)
        self.label_61.setObjectName(u"label_61")

        self.gridLayout_10.addWidget(self.label_61, 5, 1, 1, 1, Qt.AlignHCenter)

        self.doubleSpinBox_15 = QDoubleSpinBox(self.groupBox_7)
        self.doubleSpinBox_15.setObjectName(u"doubleSpinBox_15")

        self.gridLayout_10.addWidget(self.doubleSpinBox_15, 15, 0, 1, 1)

        self.q_alpha_6 = QDoubleSpinBox(self.groupBox_7)
        self.q_alpha_6.setObjectName(u"q_alpha_6")
        self.q_alpha_6.setMinimum(-33.000000000000000)
        self.q_alpha_6.setMaximum(45.000000000000000)

        self.gridLayout_10.addWidget(self.q_alpha_6, 5, 0, 1, 1)

        self.aerofoilSelection = QComboBox(self.groupBox_7)
        self.aerofoilSelection.setObjectName(u"aerofoilSelection")

        self.gridLayout_10.addWidget(self.aerofoilSelection, 16, 0, 1, 1)

        self.label_68 = QLabel(self.groupBox_7)
        self.label_68.setObjectName(u"label_68")
        self.label_68.setAlignment(Qt.AlignCenter)

        self.gridLayout_10.addWidget(self.label_68, 15, 1, 1, 1)

        self.label_20 = QLabel(self.groupBox_7)
        self.label_20.setObjectName(u"label_20")
        self.label_20.setAlignment(Qt.AlignCenter)

        self.gridLayout_10.addWidget(self.label_20, 16, 1, 1, 1)


        self.gridLayout_7.addWidget(self.groupBox_7, 0, 0, 1, 1)

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

        self.terminalOutput = QListView(self.frame_2)
        self.terminalOutput.setObjectName(u"terminalOutput")
        sizePolicy10 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        sizePolicy10.setHorizontalStretch(0)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.terminalOutput.sizePolicy().hasHeightForWidth())
        self.terminalOutput.setSizePolicy(sizePolicy10)

        self.verticalLayout.addWidget(self.terminalOutput)


        self.gridLayout_4.addWidget(self.frame_2, 4, 0, 1, 2)

        MainWindowDesign.setCentralWidget(self.centralWidget)
        self.menuBar = QMenuBar(MainWindowDesign)
        self.menuBar.setObjectName(u"menuBar")
        self.menuBar.setGeometry(QRect(0, 0, 1145, 21))
        self.menuHelp = QMenu(self.menuBar)
        self.menuHelp.setObjectName(u"menuHelp")
        self.menuExamples = QMenu(self.menuHelp)
        self.menuExamples.setObjectName(u"menuExamples")
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
        self.menuHelp.addAction(self.menuExamples.menuAction())
        self.menuHelp.addAction(self.actionManual)
        self.menuHelp.addSeparator()
        self.menuHelp.addAction(self.actionQuit)
        self.menuExamples.addAction(self.actionExample_1)
        self.menuExamples.addAction(self.actionExample_2)
        self.menuExamples.addAction(self.actionExample_3)
        self.menuExamples.addAction(self.actionExample_4)
        self.menuExamples.addAction(self.actionExample_5)
        self.menuFile.addAction(self.actionImport_Data)
        self.menuFile.addAction(self.actionExport_Data)
        self.menuAircraft.addAction(self.actionPreset_Aircraft_List)

        self.retranslateUi(MainWindowDesign)

        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(1)


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
        self.actionExample_1.setText(QCoreApplication.translate("MainWindowDesign", u"Example 1: Steady horseshoe vortex lattice method solver", None))
        self.actionExample_2.setText(QCoreApplication.translate("MainWindowDesign", u"Example 2: Steady ring vortex lattice method solver", None))
        self.actionExample_3.setText(QCoreApplication.translate("MainWindowDesign", u"Example 3: Unsteady ring vortex lattice method solver static", None))
        self.actionExample_4.setText(QCoreApplication.translate("MainWindowDesign", u"Example 4: Unsteady ring vortex lattice method solver variable", None))
        self.actionExample_5.setText(QCoreApplication.translate("MainWindowDesign", u"Example 5: Unsteady ring vortex lattice method solver variable formation", None))
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
        self.label_12.setText(QCoreApplication.translate("MainWindowDesign", u"Chord-wise Panels", None))
        self.comboBox.setItemText(0, QCoreApplication.translate("MainWindowDesign", u"Uniform", None))
        self.comboBox.setItemText(1, QCoreApplication.translate("MainWindowDesign", u"Cosine", None))
        self.comboBox.setItemText(2, QCoreApplication.translate("MainWindowDesign", u"Reverse Cosine", None))

        self.label_3.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.MainWingLabel.setText(QCoreApplication.translate("MainWindowDesign", u"Main Wing", None))
        self.groupBox_8.setTitle("")
        self.label_69.setText(QCoreApplication.translate("MainWindowDesign", u"Chord and Control Surface", None))
        self.label_70.setText(QCoreApplication.translate("MainWindowDesign", u"Span-wise panels", None))
        self.checkBox_6.setText(QCoreApplication.translate("MainWindowDesign", u"Symmetric?", None))
        self.q_x_7.setSuffix("")
        self.label_71.setText(QCoreApplication.translate("MainWindowDesign", u"Panel Spacing Option", None))
        self.label_72.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.label_73.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_74.setText(QCoreApplication.translate("MainWindowDesign", u"Number of Panels", None))
        self.comboBox_6.setItemText(0, QCoreApplication.translate("MainWindowDesign", u"Uniform", None))
        self.comboBox_6.setItemText(1, QCoreApplication.translate("MainWindowDesign", u"Cosine", None))
        self.comboBox_6.setItemText(2, QCoreApplication.translate("MainWindowDesign", u"Reverse Cosine", None))

        self.label_75.setText(QCoreApplication.translate("MainWindowDesign", u"Chord", None))
        self.label_76.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Hinge Point", None))
        self.MainWingLabel_6.setText(QCoreApplication.translate("MainWindowDesign", u"Root Cross Section", None))
        self.label_77.setText(QCoreApplication.translate("MainWindowDesign", u"Leading Edge Location", None))
        self.label_78.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_79.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Deflection", None))
        self.tabWidget_2.setTabText(self.tabWidget_2.indexOf(self.tab_7), QCoreApplication.translate("MainWindowDesign", u"Wing Root Cross Section", None))
        self.groupBox_7.setTitle("")
        self.label_62.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_60.setText(QCoreApplication.translate("MainWindowDesign", u"Panel Spacing Option", None))
        self.q_x_6.setSuffix("")
        self.label_63.setText(QCoreApplication.translate("MainWindowDesign", u"Number of Panels", None))
        self.comboBox_5.setItemText(0, QCoreApplication.translate("MainWindowDesign", u"Uniform", None))
        self.comboBox_5.setItemText(1, QCoreApplication.translate("MainWindowDesign", u"Cosine", None))
        self.comboBox_5.setItemText(2, QCoreApplication.translate("MainWindowDesign", u"Reverse Cosine", None))

        self.label_58.setText(QCoreApplication.translate("MainWindowDesign", u"Chord and Control Surface", None))
        self.MainWingLabel_5.setText(QCoreApplication.translate("MainWindowDesign", u"Root Cross Section", None))
        self.label_65.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Hinge Point", None))
        self.label_66.setText(QCoreApplication.translate("MainWindowDesign", u"Leading Edge Location", None))
        self.checkBox_5.setText(QCoreApplication.translate("MainWindowDesign", u"Symmetric?", None))
        self.label_67.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_64.setText(QCoreApplication.translate("MainWindowDesign", u"Chord", None))
        self.label_59.setText(QCoreApplication.translate("MainWindowDesign", u"Span-wise Panels", None))
        self.label_61.setText(QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.label_68.setText(QCoreApplication.translate("MainWindowDesign", u"Control Surface Deflection", None))
        self.label_20.setText(QCoreApplication.translate("MainWindowDesign", u"Aerofoil", None))
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
        self.menuExamples.setTitle(QCoreApplication.translate("MainWindowDesign", u"Examples", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindowDesign", u"File", None))
        self.menuAircraft.setTitle(QCoreApplication.translate("MainWindowDesign", u"Aircraft", None))
    # retranslateUi

