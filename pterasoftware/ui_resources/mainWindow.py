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
        MainWindowDesign.resize(837, 731)
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
        self.frame_2 = QFrame(self.centralWidget)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Raised)
        self.verticalLayout = QVBoxLayout(self.frame_2)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.label_25 = QLabel(self.frame_2)
        self.label_25.setObjectName(u"label_25")

        self.verticalLayout.addWidget(self.label_25)

        self.view_logging = QListView(self.frame_2)
        self.view_logging.setObjectName(u"view_logging")

        self.verticalLayout.addWidget(self.view_logging)

        self.gridLayout_4.addWidget(self.frame_2, 4, 0, 1, 2)

        self.frame = QFrame(self.centralWidget)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.horizontalLayout = QHBoxLayout(self.frame)
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout.setObjectName(u"horizontalLayout")

        self.gridLayout_4.addWidget(self.frame, 2, 0, 1, 1)

        self.frame_3 = QFrame(self.centralWidget)
        self.frame_3.setObjectName(u"frame_3")
        self.frame_3.setAutoFillBackground(False)
        self.frame_3.setFrameShape(QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QFrame.Raised)
        self.gridLayout_5 = QGridLayout(self.frame_3)
        self.gridLayout_5.setSpacing(6)
        self.gridLayout_5.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_5.setObjectName(u"gridLayout_5")
        self.label_17 = QLabel(self.frame_3)
        self.label_17.setObjectName(u"label_17")
        self.label_17.setEnabled(True)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_17.sizePolicy().hasHeightForWidth())
        self.label_17.setSizePolicy(sizePolicy)
        self.label_17.setMinimumSize(QSize(500, 0))
        self.label_17.setMaximumSize(QSize(799, 16777215))

        self.gridLayout_5.addWidget(self.label_17, 0, 0, 1, 1)

        self.label_16 = QLabel(self.frame_3)
        self.label_16.setObjectName(u"label_16")

        self.gridLayout_5.addWidget(self.label_16, 0, 1, 1, 1)

        self.gridLayout_4.addWidget(self.frame_3, 1, 0, 1, 2)

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
        self.groupBox = QGroupBox(self.tab_2)
        self.groupBox.setObjectName(u"groupBox")
        self.gridLayout = QGridLayout(self.groupBox)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setObjectName(u"gridLayout")
        self.label_3 = QLabel(self.groupBox)
        self.label_3.setObjectName(u"label_3")

        self.gridLayout.addWidget(self.label_3, 2, 1, 1, 1, Qt.AlignHCenter)

        self.q_alpha = QDoubleSpinBox(self.groupBox)
        self.q_alpha.setObjectName(u"q_alpha")
        self.q_alpha.setMinimum(-33.000000000000000)
        self.q_alpha.setMaximum(45.000000000000000)

        self.gridLayout.addWidget(self.q_alpha, 2, 0, 1, 1)

        self.label = QLabel(self.groupBox)
        self.label.setObjectName(u"label")

        self.gridLayout.addWidget(self.label, 0, 1, 1, 1, Qt.AlignHCenter)

        self.q_x = QDoubleSpinBox(self.groupBox)
        self.q_x.setObjectName(u"q_x")
        self.q_x.setMaximum(49.000000000000000)

        self.gridLayout.addWidget(self.q_x, 0, 0, 1, 1)

        self.q_y = QDoubleSpinBox(self.groupBox)
        self.q_y.setObjectName(u"q_y")
        self.q_y.setMaximum(49.000000000000000)

        self.gridLayout.addWidget(self.q_y, 1, 0, 1, 1)

        self.q_beta = QDoubleSpinBox(self.groupBox)
        self.q_beta.setObjectName(u"q_beta")
        self.q_beta.setMinimum(-50.000000000000000)
        self.q_beta.setMaximum(44.000000000000000)

        self.gridLayout.addWidget(self.q_beta, 3, 0, 1, 1)

        self.q_injection = QPushButton(self.groupBox)
        self.q_injection.setObjectName(u"q_injection")

        self.gridLayout.addWidget(self.q_injection, 5, 1, 1, 1)

        self.label_4 = QLabel(self.groupBox)
        self.label_4.setObjectName(u"label_4")

        self.gridLayout.addWidget(self.label_4, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_2 = QLabel(self.groupBox)
        self.label_2.setObjectName(u"label_2")

        self.gridLayout.addWidget(self.label_2, 1, 1, 1, 1, Qt.AlignHCenter)

        self.q_d = QDoubleSpinBox(self.groupBox)
        self.q_d.setObjectName(u"q_d")
        self.q_d.setMinimum(-45.000000000000000)
        self.q_d.setMaximum(0.000000000000000)

        self.gridLayout.addWidget(self.q_d, 4, 0, 1, 1)

        self.q_homing = QPushButton(self.groupBox)
        self.q_homing.setObjectName(u"q_homing")

        self.gridLayout.addWidget(self.q_homing, 5, 2, 1, 1)

        self.label_5 = QLabel(self.groupBox)
        self.label_5.setObjectName(u"label_5")

        self.gridLayout.addWidget(self.label_5, 3, 1, 1, 1, Qt.AlignHCenter)

        self.q_start = QPushButton(self.groupBox)
        self.q_start.setObjectName(u"q_start")
        self.q_start.setContextMenuPolicy(Qt.NoContextMenu)

        self.gridLayout.addWidget(self.q_start, 5, 0, 1, 1)

        self.horizontalLayout_2.addWidget(self.groupBox)

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

        MainWindowDesign.setCentralWidget(self.centralWidget)
        self.menuBar = QMenuBar(MainWindowDesign)
        self.menuBar.setObjectName(u"menuBar")
        self.menuBar.setGeometry(QRect(0, 0, 837, 21))
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

        self.tabWidget.setCurrentIndex(5)

        QMetaObject.connectSlotsByName(MainWindowDesign)

    # setupUi

    def retranslateUi(self, MainWindowDesign):
        MainWindowDesign.setWindowTitle(QCoreApplication.translate("MainWindowDesign", u"MainWindow", None))
        self.actionAbout.setText(QCoreApplication.translate("MainWindowDesign", u"About", None))
        self.actionManual.setText(QCoreApplication.translate("MainWindowDesign", u"Manual", None))
        self.actionQuit.setText(QCoreApplication.translate("MainWindowDesign", u"Quit", None))
        self.actionImport_Data.setText(QCoreApplication.translate("MainWindowDesign", u"Import Data", None))
        self.actionExport_Data.setText(QCoreApplication.translate("MainWindowDesign", u"Export Data", None))
        self.actionPreset_Aircraft_List.setText(
            QCoreApplication.translate("MainWindowDesign", u"Preset Aircraft List", None))
        self.label_25.setText(QCoreApplication.translate("MainWindowDesign", u"Terminal Output", None))
        self.label_17.setText(QCoreApplication.translate("MainWindowDesign",
                                                         u"<html><head/><body><p><span style=\" font-size:36pt; font-weight:600;\">PteraSoftware</span></p></body></html>",
                                                         None))
        self.label_16.setText(QCoreApplication.translate("MainWindowDesign",
                                                         u"<html><head/><body><p><img src=\":/images/icon.png\"/></p></body></html>",
                                                         None))
        self.groupBox.setTitle("")
        self.label_3.setText(
            QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03b1</p></body></html>", None))
        self.label.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.q_x.setSuffix("")
        self.q_injection.setText(QCoreApplication.translate("MainWindowDesign", u"Plot graphs", None))
        self.label_4.setText(QCoreApplication.translate("MainWindowDesign", u"d", None))
        self.label_2.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.q_homing.setText(QCoreApplication.translate("MainWindowDesign", u"Plot Visualisation", None))
        self.label_5.setText(
            QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03b2</p></body></html>", None))
        self.q_start.setText(QCoreApplication.translate("MainWindowDesign", u"Generate Results", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2),
                                  QCoreApplication.translate("MainWindowDesign", u"Aircraft Parameters", None))
        self.groupBox_2.setTitle(QCoreApplication.translate("MainWindowDesign", u"Controller", None))
        self.label_10.setText(
            QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03b8</p></body></html>", None))
        self.ee_start.setText(QCoreApplication.translate("MainWindowDesign", u"Move", None))
        self.ee_homing.setText(QCoreApplication.translate("MainWindowDesign", u"Homing", None))
        self.label_8.setText(QCoreApplication.translate("MainWindowDesign", u"z", None))
        self.label_6.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_9.setText(
            QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03c6</p></body></html>", None))
        self.label_7.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.ee_injection.setText(QCoreApplication.translate("MainWindowDesign", u"Injection", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab),
                                  QCoreApplication.translate("MainWindowDesign", u"Model ", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_4),
                                  QCoreApplication.translate("MainWindowDesign", u"Materials", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_5),
                                  QCoreApplication.translate("MainWindowDesign", u"Boundary Conditions", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_6),
                                  QCoreApplication.translate("MainWindowDesign", u"Solution", None))
        self.groupBox_4.setTitle("")
        self.q_start_2.setText(QCoreApplication.translate("MainWindowDesign", u"Move", None))
        self.q_homing_2.setText(QCoreApplication.translate("MainWindowDesign", u"Homing", None))
        self.q_x_2.setSuffix("")
        self.q_injection_2.setText(QCoreApplication.translate("MainWindowDesign", u"Injection", None))
        self.label_21.setText(QCoreApplication.translate("MainWindowDesign", u"x", None))
        self.label_22.setText(QCoreApplication.translate("MainWindowDesign", u"y", None))
        self.label_26.setText(
            QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>z</p></body></html>", None))
        self.label_27.setText(
            QCoreApplication.translate("MainWindowDesign", u"<html><head/><body><p>\u03b2</p></body></html>", None))
        self.label_28.setText(QCoreApplication.translate("MainWindowDesign", u"d", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3),
                                  QCoreApplication.translate("MainWindowDesign", u"Visualiser Settings", None))
        self.menuHelp.setTitle(QCoreApplication.translate("MainWindowDesign", u"Help", None))
        self.menuFile.setTitle(QCoreApplication.translate("MainWindowDesign", u"File", None))
        self.menuAircraft.setTitle(QCoreApplication.translate("MainWindowDesign", u"Aircraft", None))
    # retranslateUi
