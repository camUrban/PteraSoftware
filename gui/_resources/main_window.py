# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'main_window.ui'
##
## Created by: Qt User Interface Compiler version 6.9.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (
    QCoreApplication,
    QDate,
    QDateTime,
    QLocale,
    QMetaObject,
    QObject,
    QPoint,
    QRect,
    QSize,
    Qt,
    QTime,
    QUrl,
)
from PySide6.QtGui import (
    QAction,
    QBrush,
    QColor,
    QConicalGradient,
    QCursor,
    QFont,
    QFontDatabase,
    QGradient,
    QIcon,
    QImage,
    QKeySequence,
    QLinearGradient,
    QPainter,
    QPalette,
    QPixmap,
    QRadialGradient,
    QTransform,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMenu,
    QMenuBar,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QStatusBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)


class Ui_MainWindowDesign(object):
    def setupUi(self, MainWindowDesign):
        if not MainWindowDesign.objectName():
            MainWindowDesign.setObjectName("MainWindowDesign")
        MainWindowDesign.resize(1145, 880)
        sizePolicy = QSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindowDesign.sizePolicy().hasHeightForWidth())
        MainWindowDesign.setSizePolicy(sizePolicy)
        MainWindowDesign.setMinimumSize(QSize(1052, 880))
        self.actionAbout = QAction(MainWindowDesign)
        self.actionAbout.setObjectName("actionAbout")
        self.actionManual = QAction(MainWindowDesign)
        self.actionManual.setObjectName("actionManual")
        self.actionQuit = QAction(MainWindowDesign)
        self.actionQuit.setObjectName("actionQuit")
        self.actionImport_Data = QAction(MainWindowDesign)
        self.actionImport_Data.setObjectName("actionImport_Data")
        self.actionExport_Data = QAction(MainWindowDesign)
        self.actionExport_Data.setObjectName("actionExport_Data")
        self.actionPreset_Aircraft_List = QAction(MainWindowDesign)
        self.actionPreset_Aircraft_List.setObjectName("actionPreset_Aircraft_List")
        self.actionExample_1 = QAction(MainWindowDesign)
        self.actionExample_1.setObjectName("actionExample_1")
        self.actionExample_2 = QAction(MainWindowDesign)
        self.actionExample_2.setObjectName("actionExample_2")
        self.actionExample_3 = QAction(MainWindowDesign)
        self.actionExample_3.setObjectName("actionExample_3")
        self.actionExample_4 = QAction(MainWindowDesign)
        self.actionExample_4.setObjectName("actionExample_4")
        self.actionExample_5 = QAction(MainWindowDesign)
        self.actionExample_5.setObjectName("actionExample_5")
        self.actionExample_6 = QAction(MainWindowDesign)
        self.actionExample_6.setObjectName("actionExample_6")
        self.actionExample_7 = QAction(MainWindowDesign)
        self.actionExample_7.setObjectName("actionExample_7")
        self.actionExample_8 = QAction(MainWindowDesign)
        self.actionExample_8.setObjectName("actionExample_8")
        self.actionExample_9 = QAction(MainWindowDesign)
        self.actionExample_9.setObjectName("actionExample_9")
        self.actionExample_10 = QAction(MainWindowDesign)
        self.actionExample_10.setObjectName("actionExample_10")
        self.centralWidget = QWidget(MainWindowDesign)
        self.centralWidget.setObjectName("centralWidget")
        self.gridLayout_4 = QGridLayout(self.centralWidget)
        self.gridLayout_4.setSpacing(6)
        self.gridLayout_4.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.frame = QFrame(self.centralWidget)
        self.frame.setObjectName("frame")
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.horizontalLayout = QHBoxLayout(self.frame)
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout.setObjectName("horizontalLayout")

        self.gridLayout_4.addWidget(self.frame, 2, 0, 1, 1)

        self.tabWidget = QTabWidget(self.centralWidget)
        self.tabWidget.setObjectName("tabWidget")
        sizePolicy1 = QSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy1)
        self.tabWidget.setToolTipDuration(-1)
        self.aircraftParametersTab = QWidget()
        self.aircraftParametersTab.setObjectName("aircraftParametersTab")
        self.horizontalLayout_2 = QHBoxLayout(self.aircraftParametersTab)
        self.horizontalLayout_2.setSpacing(6)
        self.horizontalLayout_2.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.verticalLayout_2 = QVBoxLayout()
        self.verticalLayout_2.setSpacing(6)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_13 = QLabel(self.aircraftParametersTab)
        self.label_13.setObjectName("label_13")
        sizePolicy2 = QSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.label_13.sizePolicy().hasHeightForWidth())
        self.label_13.setSizePolicy(sizePolicy2)
        font = QFont()
        font.setPointSize(12)
        font.setBold(True)
        self.label_13.setFont(font)
        self.label_13.setAlignment(Qt.AlignCenter)

        self.verticalLayout_2.addWidget(self.label_13)

        self.label_14 = QLabel(self.aircraftParametersTab)
        self.label_14.setObjectName("label_14")
        sizePolicy2.setHeightForWidth(self.label_14.sizePolicy().hasHeightForWidth())
        self.label_14.setSizePolicy(sizePolicy2)
        font1 = QFont()
        font1.setPointSize(10)
        font1.setBold(True)
        self.label_14.setFont(font1)

        self.verticalLayout_2.addWidget(self.label_14)

        self.nameInput = QLineEdit(self.aircraftParametersTab)
        self.nameInput.setObjectName("nameInput")
        sizePolicy3 = QSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed
        )
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.nameInput.sizePolicy().hasHeightForWidth())
        self.nameInput.setSizePolicy(sizePolicy3)

        self.verticalLayout_2.addWidget(self.nameInput)

        self.label_15 = QLabel(self.aircraftParametersTab)
        self.label_15.setObjectName("label_15")
        sizePolicy2.setHeightForWidth(self.label_15.sizePolicy().hasHeightForWidth())
        self.label_15.setSizePolicy(sizePolicy2)
        self.label_15.setFont(font1)

        self.verticalLayout_2.addWidget(self.label_15)

        self.label_16 = QLabel(self.aircraftParametersTab)
        self.label_16.setObjectName("label_16")
        sizePolicy4 = QSizePolicy(
            QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum
        )
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy4)
        font2 = QFont()
        font2.setPointSize(10)
        self.label_16.setFont(font2)

        self.verticalLayout_2.addWidget(self.label_16)

        self.x_CofG = QDoubleSpinBox(self.aircraftParametersTab)
        self.x_CofG.setObjectName("x_CofG")
        sizePolicy5 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Maximum)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.x_CofG.sizePolicy().hasHeightForWidth())
        self.x_CofG.setSizePolicy(sizePolicy5)
        self.x_CofG.setMinimumSize(QSize(149, 0))

        self.verticalLayout_2.addWidget(self.x_CofG)

        self.label_17 = QLabel(self.aircraftParametersTab)
        self.label_17.setObjectName("label_17")
        sizePolicy4.setHeightForWidth(self.label_17.sizePolicy().hasHeightForWidth())
        self.label_17.setSizePolicy(sizePolicy4)
        self.label_17.setFont(font2)

        self.verticalLayout_2.addWidget(self.label_17)

        self.y_CofG = QDoubleSpinBox(self.aircraftParametersTab)
        self.y_CofG.setObjectName("y_CofG")
        sizePolicy6 = QSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)
        sizePolicy6.setHeightForWidth(self.y_CofG.sizePolicy().hasHeightForWidth())
        self.y_CofG.setSizePolicy(sizePolicy6)
        self.y_CofG.setMinimumSize(QSize(149, 0))

        self.verticalLayout_2.addWidget(self.y_CofG)

        self.label_18 = QLabel(self.aircraftParametersTab)
        self.label_18.setObjectName("label_18")
        sizePolicy4.setHeightForWidth(self.label_18.sizePolicy().hasHeightForWidth())
        self.label_18.setSizePolicy(sizePolicy4)

        self.verticalLayout_2.addWidget(self.label_18)

        self.z_CofG = QDoubleSpinBox(self.aircraftParametersTab)
        self.z_CofG.setObjectName("z_CofG")
        sizePolicy6.setHeightForWidth(self.z_CofG.sizePolicy().hasHeightForWidth())
        self.z_CofG.setSizePolicy(sizePolicy6)
        self.z_CofG.setMinimumSize(QSize(149, 0))

        self.verticalLayout_2.addWidget(self.z_CofG)

        self.label_19 = QLabel(self.aircraftParametersTab)
        self.label_19.setObjectName("label_19")
        self.label_19.setAlignment(Qt.AlignCenter)

        self.verticalLayout_2.addWidget(self.label_19)

        self.saveAircraftButtom = QPushButton(self.aircraftParametersTab)
        self.saveAircraftButtom.setObjectName("saveAircraftButtom")

        self.verticalLayout_2.addWidget(self.saveAircraftButtom)

        self.horizontalLayout_2.addLayout(self.verticalLayout_2)

        self.groupBox = QGroupBox(self.aircraftParametersTab)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QGridLayout(self.groupBox)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QLabel(self.groupBox)
        self.label.setObjectName("label")

        self.gridLayout.addWidget(self.label, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_2 = QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")

        self.gridLayout.addWidget(self.label_2, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_11 = QLabel(self.groupBox)
        self.label_11.setObjectName("label_11")
        sizePolicy2.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy2)
        self.label_11.setFont(font1)

        self.gridLayout.addWidget(self.label_11, 2, 0, 1, 2)

        self.LE_loc_z = QDoubleSpinBox(self.groupBox)
        self.LE_loc_z.setObjectName("LE_loc_z")
        self.LE_loc_z.setMinimum(-33.000000000000000)
        self.LE_loc_z.setMaximum(45.000000000000000)

        self.gridLayout.addWidget(self.LE_loc_z, 5, 0, 1, 1)

        self.LE_loc_x = QDoubleSpinBox(self.groupBox)
        self.LE_loc_x.setObjectName("LE_loc_x")
        self.LE_loc_x.setMaximum(49.000000000000000)

        self.gridLayout.addWidget(self.LE_loc_x, 3, 0, 1, 1)

        self.label_5 = QLabel(self.groupBox)
        self.label_5.setObjectName("label_5")
        sizePolicy7 = QSizePolicy(
            QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Preferred
        )
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy7)

        self.gridLayout.addWidget(self.label_5, 8, 1, 1, 1, Qt.AlignHCenter)

        self.LE_loc_y = QDoubleSpinBox(self.groupBox)
        self.LE_loc_y.setObjectName("LE_loc_y")
        self.LE_loc_y.setMaximum(49.000000000000000)

        self.gridLayout.addWidget(self.LE_loc_y, 4, 0, 1, 1)

        self.mainWingSymmetryCheckbox = QCheckBox(self.groupBox)
        self.mainWingSymmetryCheckbox.setObjectName("mainWingSymmetryCheckbox")

        self.gridLayout.addWidget(self.mainWingSymmetryCheckbox, 11, 0, 1, 1)

        self.label_4 = QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")

        self.gridLayout.addWidget(self.label_4, 10, 1, 1, 1, Qt.AlignHCenter)

        self.label_12 = QLabel(self.groupBox)
        self.label_12.setObjectName("label_12")
        sizePolicy2.setHeightForWidth(self.label_12.sizePolicy().hasHeightForWidth())
        self.label_12.setSizePolicy(sizePolicy2)
        self.label_12.setFont(font1)

        self.gridLayout.addWidget(self.label_12, 7, 0, 1, 1)

        self.CWPanelSpcaing = QComboBox(self.groupBox)
        self.CWPanelSpcaing.addItem("")
        self.CWPanelSpcaing.addItem("")
        self.CWPanelSpcaing.addItem("")
        self.CWPanelSpcaing.setObjectName("CWPanelSpcaing")

        self.gridLayout.addWidget(self.CWPanelSpcaing, 10, 0, 1, 1)

        self.label_3 = QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")

        self.gridLayout.addWidget(self.label_3, 5, 1, 1, 1, Qt.AlignHCenter)

        self.panelNum_CWPanels = QDoubleSpinBox(self.groupBox)
        self.panelNum_CWPanels.setObjectName("panelNum_CWPanels")
        self.panelNum_CWPanels.setMinimum(-50.000000000000000)
        self.panelNum_CWPanels.setMaximum(44.000000000000000)

        self.gridLayout.addWidget(self.panelNum_CWPanels, 8, 0, 1, 1)

        self.MainWingLabel = QLabel(self.groupBox)
        self.MainWingLabel.setObjectName("MainWingLabel")
        sizePolicy8 = QSizePolicy(
            QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Maximum
        )
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(0)
        sizePolicy8.setHeightForWidth(
            self.MainWingLabel.sizePolicy().hasHeightForWidth()
        )
        self.MainWingLabel.setSizePolicy(sizePolicy8)
        self.MainWingLabel.setFont(font)
        self.MainWingLabel.setFrameShadow(QFrame.Sunken)
        self.MainWingLabel.setLineWidth(1)
        self.MainWingLabel.setAlignment(Qt.AlignCenter)

        self.gridLayout.addWidget(self.MainWingLabel, 0, 0, 1, 2)

        self.horizontalLayout_2.addWidget(self.groupBox)

        self.tabWidget_2 = QTabWidget(self.aircraftParametersTab)
        self.tabWidget_2.setObjectName("tabWidget_2")
        sizePolicy9 = QSizePolicy(
            QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Expanding
        )
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.tabWidget_2.sizePolicy().hasHeightForWidth())
        self.tabWidget_2.setSizePolicy(sizePolicy9)
        self.tab_7 = QWidget()
        self.tab_7.setObjectName("tab_7")
        self.tab_7.setMinimumSize(QSize(318, 0))
        self.gridLayout_3 = QGridLayout(self.tab_7)
        self.gridLayout_3.setSpacing(6)
        self.gridLayout_3.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.groupBox_8 = QGroupBox(self.tab_7)
        self.groupBox_8.setObjectName("groupBox_8")
        self.gridLayout_11 = QGridLayout(self.groupBox_8)
        self.gridLayout_11.setSpacing(6)
        self.gridLayout_11.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_11.setObjectName("gridLayout_11")
        self.WRCW_LE_loc_y = QDoubleSpinBox(self.groupBox_8)
        self.WRCW_LE_loc_y.setObjectName("WRCW_LE_loc_y")
        self.WRCW_LE_loc_y.setMaximum(49.000000000000000)

        self.gridLayout_11.addWidget(self.WRCW_LE_loc_y, 4, 0, 1, 1)

        self.panelNum_SWPanels = QDoubleSpinBox(self.groupBox_8)
        self.panelNum_SWPanels.setObjectName("panelNum_SWPanels")
        self.panelNum_SWPanels.setMinimum(-50.000000000000000)
        self.panelNum_SWPanels.setMaximum(44.000000000000000)

        self.gridLayout_11.addWidget(self.panelNum_SWPanels, 8, 0, 1, 1)

        self.WRCS_CSHP = QDoubleSpinBox(self.groupBox_8)
        self.WRCS_CSHP.setObjectName("WRCS_CSHP")

        self.gridLayout_11.addWidget(self.WRCS_CSHP, 14, 0, 1, 1)

        self.label_74 = QLabel(self.groupBox_8)
        self.label_74.setObjectName("label_74")
        sizePolicy7.setHeightForWidth(self.label_74.sizePolicy().hasHeightForWidth())
        self.label_74.setSizePolicy(sizePolicy7)
        self.label_74.setMinimumSize(QSize(0, 0))

        self.gridLayout_11.addWidget(self.label_74, 8, 1, 1, 1, Qt.AlignHCenter)

        self.label_78 = QLabel(self.groupBox_8)
        self.label_78.setObjectName("label_78")

        self.gridLayout_11.addWidget(self.label_78, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_69 = QLabel(self.groupBox_8)
        self.label_69.setObjectName("label_69")
        font3 = QFont()
        font3.setPointSize(10)
        font3.setBold(False)
        self.label_69.setFont(font3)

        self.gridLayout_11.addWidget(self.label_69, 12, 0, 1, 1)

        self.label_75 = QLabel(self.groupBox_8)
        self.label_75.setObjectName("label_75")

        self.gridLayout_11.addWidget(self.label_75, 13, 1, 1, 1)

        self.WRCS_CSD = QDoubleSpinBox(self.groupBox_8)
        self.WRCS_CSD.setObjectName("WRCS_CSD")

        self.gridLayout_11.addWidget(self.WRCS_CSD, 15, 0, 1, 1)

        self.label_79 = QLabel(self.groupBox_8)
        self.label_79.setObjectName("label_79")

        self.gridLayout_11.addWidget(self.label_79, 15, 1, 1, 1)

        self.label_73 = QLabel(self.groupBox_8)
        self.label_73.setObjectName("label_73")

        self.gridLayout_11.addWidget(self.label_73, 4, 1, 1, 1, Qt.AlignHCenter)

        self.WRCS_Chord = QDoubleSpinBox(self.groupBox_8)
        self.WRCS_Chord.setObjectName("WRCS_Chord")

        self.gridLayout_11.addWidget(self.WRCS_Chord, 13, 0, 1, 1)

        self.label_70 = QLabel(self.groupBox_8)
        self.label_70.setObjectName("label_70")
        sizePolicy2.setHeightForWidth(self.label_70.sizePolicy().hasHeightForWidth())
        self.label_70.setSizePolicy(sizePolicy2)
        self.label_70.setFont(font3)

        self.gridLayout_11.addWidget(self.label_70, 7, 0, 1, 1)

        self.label_72 = QLabel(self.groupBox_8)
        self.label_72.setObjectName("label_72")

        self.gridLayout_11.addWidget(self.label_72, 5, 1, 1, 1, Qt.AlignHCenter)

        self.label_71 = QLabel(self.groupBox_8)
        self.label_71.setObjectName("label_71")

        self.gridLayout_11.addWidget(self.label_71, 10, 1, 1, 1, Qt.AlignHCenter)

        self.WRCW_LE_loc_z = QDoubleSpinBox(self.groupBox_8)
        self.WRCW_LE_loc_z.setObjectName("WRCW_LE_loc_z")
        self.WRCW_LE_loc_z.setMinimum(-33.000000000000000)
        self.WRCW_LE_loc_z.setMaximum(45.000000000000000)

        self.gridLayout_11.addWidget(self.WRCW_LE_loc_z, 5, 0, 1, 1)

        self.label_77 = QLabel(self.groupBox_8)
        self.label_77.setObjectName("label_77")
        sizePolicy2.setHeightForWidth(self.label_77.sizePolicy().hasHeightForWidth())
        self.label_77.setSizePolicy(sizePolicy2)
        self.label_77.setFont(font3)

        self.gridLayout_11.addWidget(self.label_77, 2, 0, 1, 2)

        self.label_76 = QLabel(self.groupBox_8)
        self.label_76.setObjectName("label_76")

        self.gridLayout_11.addWidget(self.label_76, 14, 1, 1, 1)

        self.WRCW_LE_loc_x = QDoubleSpinBox(self.groupBox_8)
        self.WRCW_LE_loc_x.setObjectName("WRCW_LE_loc_x")
        self.WRCW_LE_loc_x.setMaximum(49.000000000000000)

        self.gridLayout_11.addWidget(self.WRCW_LE_loc_x, 3, 0, 1, 1)

        self.SWPanelSpcaing = QComboBox(self.groupBox_8)
        self.SWPanelSpcaing.addItem("")
        self.SWPanelSpcaing.addItem("")
        self.SWPanelSpcaing.addItem("")
        self.SWPanelSpcaing.setObjectName("SWPanelSpcaing")

        self.gridLayout_11.addWidget(self.SWPanelSpcaing, 10, 0, 1, 1)

        self.WRCS_S_symCheckbox = QCheckBox(self.groupBox_8)
        self.WRCS_S_symCheckbox.setObjectName("WRCS_S_symCheckbox")

        self.gridLayout_11.addWidget(self.WRCS_S_symCheckbox, 11, 0, 1, 1)

        self.MainWingLabel_6 = QLabel(self.groupBox_8)
        self.MainWingLabel_6.setObjectName("MainWingLabel_6")
        sizePolicy8.setHeightForWidth(
            self.MainWingLabel_6.sizePolicy().hasHeightForWidth()
        )
        self.MainWingLabel_6.setSizePolicy(sizePolicy8)
        font4 = QFont()
        font4.setPointSize(12)
        font4.setBold(False)
        self.MainWingLabel_6.setFont(font4)
        self.MainWingLabel_6.setLineWidth(1)

        self.gridLayout_11.addWidget(self.MainWingLabel_6, 0, 0, 1, 2)

        self.WRCS_aerofoil = QComboBox(self.groupBox_8)
        self.WRCS_aerofoil.setObjectName("WRCS_aerofoil")

        self.gridLayout_11.addWidget(self.WRCS_aerofoil, 16, 0, 1, 1)

        self.label_23 = QLabel(self.groupBox_8)
        self.label_23.setObjectName("label_23")
        sizePolicy.setHeightForWidth(self.label_23.sizePolicy().hasHeightForWidth())
        self.label_23.setSizePolicy(sizePolicy)
        self.label_23.setToolTipDuration(0)

        self.gridLayout_11.addWidget(self.label_23, 16, 1, 1, 1)

        self.gridLayout_3.addWidget(self.groupBox_8, 0, 0, 1, 1)

        self.tabWidget_2.addTab(self.tab_7, "")
        self.tab_8 = QWidget()
        self.tab_8.setObjectName("tab_8")
        self.gridLayout_7 = QGridLayout(self.tab_8)
        self.gridLayout_7.setSpacing(6)
        self.gridLayout_7.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.groupBox_7 = QGroupBox(self.tab_8)
        self.groupBox_7.setObjectName("groupBox_7")
        self.gridLayout_10 = QGridLayout(self.groupBox_7)
        self.gridLayout_10.setSpacing(6)
        self.gridLayout_10.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_10.setObjectName("gridLayout_10")
        self.WTCS_Chord = QDoubleSpinBox(self.groupBox_7)
        self.WTCS_Chord.setObjectName("WTCS_Chord")

        self.gridLayout_10.addWidget(self.WTCS_Chord, 13, 0, 1, 1)

        self.WTCW_LE_loc_y = QDoubleSpinBox(self.groupBox_7)
        self.WTCW_LE_loc_y.setObjectName("WTCW_LE_loc_y")
        self.WTCW_LE_loc_y.setMaximum(49.000000000000000)

        self.gridLayout_10.addWidget(self.WTCW_LE_loc_y, 4, 0, 1, 1)

        self.label_62 = QLabel(self.groupBox_7)
        self.label_62.setObjectName("label_62")

        self.gridLayout_10.addWidget(self.label_62, 4, 1, 1, 1, Qt.AlignHCenter)

        self.label_60 = QLabel(self.groupBox_7)
        self.label_60.setObjectName("label_60")

        self.gridLayout_10.addWidget(self.label_60, 10, 1, 1, 1, Qt.AlignHCenter)

        self.WTCW_LE_loc_x = QDoubleSpinBox(self.groupBox_7)
        self.WTCW_LE_loc_x.setObjectName("WTCW_LE_loc_x")
        self.WTCW_LE_loc_x.setMaximum(49.000000000000000)

        self.gridLayout_10.addWidget(self.WTCW_LE_loc_x, 3, 0, 1, 1)

        self.label_63 = QLabel(self.groupBox_7)
        self.label_63.setObjectName("label_63")
        sizePolicy7.setHeightForWidth(self.label_63.sizePolicy().hasHeightForWidth())
        self.label_63.setSizePolicy(sizePolicy7)
        self.label_63.setMinimumSize(QSize(0, 0))

        self.gridLayout_10.addWidget(self.label_63, 8, 1, 1, 1, Qt.AlignHCenter)

        self.SWPanelSpcaing_WTCS = QComboBox(self.groupBox_7)
        self.SWPanelSpcaing_WTCS.addItem("")
        self.SWPanelSpcaing_WTCS.addItem("")
        self.SWPanelSpcaing_WTCS.addItem("")
        self.SWPanelSpcaing_WTCS.setObjectName("SWPanelSpcaing_WTCS")

        self.gridLayout_10.addWidget(self.SWPanelSpcaing_WTCS, 10, 0, 1, 1)

        self.label_58 = QLabel(self.groupBox_7)
        self.label_58.setObjectName("label_58")
        self.label_58.setFont(font3)

        self.gridLayout_10.addWidget(self.label_58, 12, 0, 1, 1)

        self.MainWingLabel_5 = QLabel(self.groupBox_7)
        self.MainWingLabel_5.setObjectName("MainWingLabel_5")
        sizePolicy8.setHeightForWidth(
            self.MainWingLabel_5.sizePolicy().hasHeightForWidth()
        )
        self.MainWingLabel_5.setSizePolicy(sizePolicy8)
        self.MainWingLabel_5.setFont(font4)
        self.MainWingLabel_5.setLineWidth(1)

        self.gridLayout_10.addWidget(self.MainWingLabel_5, 0, 0, 1, 2)

        self.WTCS_CSHP = QDoubleSpinBox(self.groupBox_7)
        self.WTCS_CSHP.setObjectName("WTCS_CSHP")

        self.gridLayout_10.addWidget(self.WTCS_CSHP, 14, 0, 1, 1)

        self.label_65 = QLabel(self.groupBox_7)
        self.label_65.setObjectName("label_65")

        self.gridLayout_10.addWidget(self.label_65, 14, 1, 1, 1)

        self.panelNum_SWPanels_WTCS = QDoubleSpinBox(self.groupBox_7)
        self.panelNum_SWPanels_WTCS.setObjectName("panelNum_SWPanels_WTCS")
        self.panelNum_SWPanels_WTCS.setMinimum(-50.000000000000000)
        self.panelNum_SWPanels_WTCS.setMaximum(44.000000000000000)

        self.gridLayout_10.addWidget(self.panelNum_SWPanels_WTCS, 8, 0, 1, 1)

        self.label_66 = QLabel(self.groupBox_7)
        self.label_66.setObjectName("label_66")
        sizePolicy2.setHeightForWidth(self.label_66.sizePolicy().hasHeightForWidth())
        self.label_66.setSizePolicy(sizePolicy2)
        self.label_66.setFont(font3)

        self.gridLayout_10.addWidget(self.label_66, 2, 0, 1, 2)

        self.WTCS_S_symCheckbox = QCheckBox(self.groupBox_7)
        self.WTCS_S_symCheckbox.setObjectName("WTCS_S_symCheckbox")

        self.gridLayout_10.addWidget(self.WTCS_S_symCheckbox, 11, 0, 1, 1)

        self.label_67 = QLabel(self.groupBox_7)
        self.label_67.setObjectName("label_67")

        self.gridLayout_10.addWidget(self.label_67, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_64 = QLabel(self.groupBox_7)
        self.label_64.setObjectName("label_64")

        self.gridLayout_10.addWidget(self.label_64, 13, 1, 1, 1)

        self.label_59 = QLabel(self.groupBox_7)
        self.label_59.setObjectName("label_59")
        sizePolicy2.setHeightForWidth(self.label_59.sizePolicy().hasHeightForWidth())
        self.label_59.setSizePolicy(sizePolicy2)
        self.label_59.setFont(font3)

        self.gridLayout_10.addWidget(self.label_59, 7, 0, 1, 1)

        self.label_61 = QLabel(self.groupBox_7)
        self.label_61.setObjectName("label_61")

        self.gridLayout_10.addWidget(self.label_61, 5, 1, 1, 1, Qt.AlignHCenter)

        self.WTCS_CSD = QDoubleSpinBox(self.groupBox_7)
        self.WTCS_CSD.setObjectName("WTCS_CSD")

        self.gridLayout_10.addWidget(self.WTCS_CSD, 15, 0, 1, 1)

        self.WTCW_LE_loc_z = QDoubleSpinBox(self.groupBox_7)
        self.WTCW_LE_loc_z.setObjectName("WTCW_LE_loc_z")
        self.WTCW_LE_loc_z.setMinimum(-33.000000000000000)
        self.WTCW_LE_loc_z.setMaximum(45.000000000000000)

        self.gridLayout_10.addWidget(self.WTCW_LE_loc_z, 5, 0, 1, 1)

        self.WTCS_aerofoil = QComboBox(self.groupBox_7)
        self.WTCS_aerofoil.setObjectName("WTCS_aerofoil")

        self.gridLayout_10.addWidget(self.WTCS_aerofoil, 16, 0, 1, 1)

        self.label_68 = QLabel(self.groupBox_7)
        self.label_68.setObjectName("label_68")

        self.gridLayout_10.addWidget(self.label_68, 15, 1, 1, 1)

        self.label_20 = QLabel(self.groupBox_7)
        self.label_20.setObjectName("label_20")

        self.gridLayout_10.addWidget(self.label_20, 16, 1, 1, 1)

        self.gridLayout_7.addWidget(self.groupBox_7, 0, 0, 1, 1)

        self.tabWidget_2.addTab(self.tab_8, "")

        self.horizontalLayout_2.addWidget(self.tabWidget_2)

        self.tabWidget.addTab(self.aircraftParametersTab, "")
        self.modellingTab = QWidget()
        self.modellingTab.setObjectName("modellingTab")
        self.horizontalLayout_3 = QHBoxLayout(self.modellingTab)
        self.horizontalLayout_3.setSpacing(6)
        self.horizontalLayout_3.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.groupBox_2 = QGroupBox(self.modellingTab)
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_2 = QGridLayout(self.groupBox_2)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_2.setObjectName("gridLayout_2")

        self.horizontalLayout_3.addWidget(self.groupBox_2)

        self.tabWidget.addTab(self.modellingTab, "")
        self.tab_4 = QWidget()
        self.tab_4.setObjectName("tab_4")
        self.tabWidget.addTab(self.tab_4, "")
        self.tab_5 = QWidget()
        self.tab_5.setObjectName("tab_5")
        self.tabWidget.addTab(self.tab_5, "")
        self.tab_6 = QWidget()
        self.tab_6.setObjectName("tab_6")
        self.tabWidget.addTab(self.tab_6, "")
        self.visualiserSettingsTab = QWidget()
        self.visualiserSettingsTab.setObjectName("visualiserSettingsTab")
        self.groupBox_4 = QGroupBox(self.visualiserSettingsTab)
        self.groupBox_4.setObjectName("groupBox_4")
        self.groupBox_4.setGeometry(QRect(0, 0, 331, 235))
        self.gridLayout_6 = QGridLayout(self.groupBox_4)
        self.gridLayout_6.setSpacing(6)
        self.gridLayout_6.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.q_start_2 = QPushButton(self.groupBox_4)
        self.q_start_2.setObjectName("q_start_2")

        self.gridLayout_6.addWidget(self.q_start_2, 5, 0, 1, 1)

        self.q_y_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_y_2.setObjectName("q_y_2")
        self.q_y_2.setMaximum(49.000000000000000)

        self.gridLayout_6.addWidget(self.q_y_2, 1, 0, 1, 1)

        self.q_homing_2 = QPushButton(self.groupBox_4)
        self.q_homing_2.setObjectName("q_homing_2")

        self.gridLayout_6.addWidget(self.q_homing_2, 5, 2, 1, 1)

        self.q_beta_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_beta_2.setObjectName("q_beta_2")
        self.q_beta_2.setMinimum(-50.000000000000000)
        self.q_beta_2.setMaximum(44.000000000000000)

        self.gridLayout_6.addWidget(self.q_beta_2, 3, 0, 1, 1)

        self.q_alpha_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_alpha_2.setObjectName("q_alpha_2")
        self.q_alpha_2.setMinimum(-33.000000000000000)
        self.q_alpha_2.setMaximum(45.000000000000000)

        self.gridLayout_6.addWidget(self.q_alpha_2, 2, 0, 1, 1)

        self.q_x_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_x_2.setObjectName("q_x_2")
        self.q_x_2.setMaximum(49.000000000000000)

        self.gridLayout_6.addWidget(self.q_x_2, 0, 0, 1, 1)

        self.q_d_2 = QDoubleSpinBox(self.groupBox_4)
        self.q_d_2.setObjectName("q_d_2")
        self.q_d_2.setMinimum(-45.000000000000000)
        self.q_d_2.setMaximum(0.000000000000000)

        self.gridLayout_6.addWidget(self.q_d_2, 4, 0, 1, 1)

        self.q_injection_2 = QPushButton(self.groupBox_4)
        self.q_injection_2.setObjectName("q_injection_2")

        self.gridLayout_6.addWidget(self.q_injection_2, 5, 1, 1, 1)

        self.label_21 = QLabel(self.groupBox_4)
        self.label_21.setObjectName("label_21")

        self.gridLayout_6.addWidget(self.label_21, 0, 1, 1, 1, Qt.AlignHCenter)

        self.label_22 = QLabel(self.groupBox_4)
        self.label_22.setObjectName("label_22")

        self.gridLayout_6.addWidget(self.label_22, 1, 1, 1, 1, Qt.AlignHCenter)

        self.label_26 = QLabel(self.groupBox_4)
        self.label_26.setObjectName("label_26")

        self.gridLayout_6.addWidget(self.label_26, 2, 1, 1, 1, Qt.AlignHCenter)

        self.label_27 = QLabel(self.groupBox_4)
        self.label_27.setObjectName("label_27")

        self.gridLayout_6.addWidget(self.label_27, 3, 1, 1, 1, Qt.AlignHCenter)

        self.label_28 = QLabel(self.groupBox_4)
        self.label_28.setObjectName("label_28")

        self.gridLayout_6.addWidget(self.label_28, 4, 1, 1, 1, Qt.AlignHCenter)

        self.tabWidget.addTab(self.visualiserSettingsTab, "")

        self.gridLayout_4.addWidget(self.tabWidget, 2, 1, 1, 1)

        self.frame_3 = QFrame(self.centralWidget)
        self.frame_3.setObjectName("frame_3")
        self.frame_3.setAutoFillBackground(False)
        self.frame_3.setFrameShape(QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QFrame.Raised)
        self.gridLayout_5 = QGridLayout(self.frame_3)
        self.gridLayout_5.setSpacing(6)
        self.gridLayout_5.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.Logo = QLabel(self.frame_3)
        self.Logo.setObjectName("Logo")
        self.Logo.setAlignment(Qt.AlignCenter)

        self.gridLayout_5.addWidget(self.Logo, 0, 0, 1, 1)

        self.gridLayout_4.addWidget(self.frame_3, 1, 0, 1, 2)

        self.progressBar = QProgressBar(self.centralWidget)
        self.progressBar.setObjectName("progressBar")
        self.progressBar.setMinimum(0)
        self.progressBar.setValue(24)

        self.gridLayout_4.addWidget(self.progressBar, 5, 0, 1, 2)

        self.frame_2 = QFrame(self.centralWidget)
        self.frame_2.setObjectName("frame_2")
        self.frame_2.setFrameShape(QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Raised)
        self.verticalLayout = QVBoxLayout(self.frame_2)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout.setObjectName("verticalLayout")
        self.genResultsButton = QPushButton(self.frame_2)
        self.genResultsButton.setObjectName("genResultsButton")
        self.genResultsButton.setContextMenuPolicy(Qt.NoContextMenu)

        self.verticalLayout.addWidget(self.genResultsButton)

        self.plotVisualsButton = QPushButton(self.frame_2)
        self.plotVisualsButton.setObjectName("plotVisualsButton")

        self.verticalLayout.addWidget(self.plotVisualsButton)

        self.terminalOutpotLogo = QLabel(self.frame_2)
        self.terminalOutpotLogo.setObjectName("terminalOutpotLogo")

        self.verticalLayout.addWidget(self.terminalOutpotLogo)

        self.terminalOutput = QListWidget(self.frame_2)
        self.terminalOutput.setObjectName("terminalOutput")
        sizePolicy10 = QSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        sizePolicy10.setHorizontalStretch(5)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(
            self.terminalOutput.sizePolicy().hasHeightForWidth()
        )
        self.terminalOutput.setSizePolicy(sizePolicy10)

        self.verticalLayout.addWidget(self.terminalOutput)

        self.gridLayout_4.addWidget(self.frame_2, 4, 0, 1, 2)

        MainWindowDesign.setCentralWidget(self.centralWidget)
        self.menuBar = QMenuBar(MainWindowDesign)
        self.menuBar.setObjectName("menuBar")
        self.menuBar.setGeometry(QRect(0, 0, 1145, 21))
        self.menuHelp = QMenu(self.menuBar)
        self.menuHelp.setObjectName("menuHelp")
        self.menuExamples = QMenu(self.menuHelp)
        self.menuExamples.setObjectName("menuExamples")
        self.menuFile = QMenu(self.menuBar)
        self.menuFile.setObjectName("menuFile")
        self.menuAircraft = QMenu(self.menuBar)
        self.menuAircraft.setObjectName("menuAircraft")
        MainWindowDesign.setMenuBar(self.menuBar)
        self.statusBar = QStatusBar(MainWindowDesign)
        self.statusBar.setObjectName("statusBar")
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
        self.menuExamples.addAction(self.actionExample_6)
        self.menuExamples.addAction(self.actionExample_7)
        self.menuExamples.addAction(self.actionExample_8)
        self.menuExamples.addAction(self.actionExample_9)
        self.menuExamples.addAction(self.actionExample_10)
        self.menuFile.addAction(self.actionImport_Data)
        self.menuFile.addAction(self.actionExport_Data)
        self.menuAircraft.addAction(self.actionPreset_Aircraft_List)

        self.retranslateUi(MainWindowDesign)

        self.tabWidget.setCurrentIndex(0)
        self.tabWidget_2.setCurrentIndex(1)

        QMetaObject.connectSlotsByName(MainWindowDesign)

    # setupUi

    def retranslateUi(self, MainWindowDesign):
        MainWindowDesign.setWindowTitle(
            QCoreApplication.translate("MainWindowDesign", "MainWindow", None)
        )
        self.actionAbout.setText(
            QCoreApplication.translate("MainWindowDesign", "About", None)
        )
        self.actionManual.setText(
            QCoreApplication.translate("MainWindowDesign", "Manual", None)
        )
        self.actionQuit.setText(
            QCoreApplication.translate("MainWindowDesign", "Quit", None)
        )
        self.actionImport_Data.setText(
            QCoreApplication.translate("MainWindowDesign", "Import Data", None)
        )
        self.actionExport_Data.setText(
            QCoreApplication.translate("MainWindowDesign", "Export Data", None)
        )
        self.actionPreset_Aircraft_List.setText(
            QCoreApplication.translate("MainWindowDesign", "Preset Aircraft List", None)
        )
        self.actionExample_1.setText(
            QCoreApplication.translate(
                "MainWindowDesign", "Example 1: Analyze steady trim example", None
            )
        )
        self.actionExample_2.setText(
            QCoreApplication.translate(
                "MainWindowDesign", "Example 2: Analyze unsteady trim example", None
            )
        )
        self.actionExample_3.setText(
            QCoreApplication.translate(
                "MainWindowDesign", "Example 3: Steady convergence example", None
            )
        )
        self.actionExample_4.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Example 4: Steady horseshoe vortex lattice method solver\n"
                "                ",
                None,
            )
        )
        self.actionExample_5.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Example 5: Steady ring vortex lattice method solver static\n"
                "                ",
                None,
            )
        )
        self.actionExample_6.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Example 6: Unsteady ring vortex lattice method solver static\n"
                "                ",
                None,
            )
        )
        self.actionExample_7.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Example 7: Unsteady ring vortex lattice method solver variable\n"
                "                ",
                None,
            )
        )
        self.actionExample_8.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Example 8: Unsteady ring vortex lattice method solver variable\n"
                "                    formation\n"
                "                ",
                None,
            )
        )
        self.actionExample_9.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Example 9: Unsteady static convergence example",
                None,
            )
        )
        self.actionExample_10.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Example 10: Unsteady variable convergence example",
                None,
            )
        )
        self.label_13.setText(
            QCoreApplication.translate("MainWindowDesign", "Aircraft", None)
        )
        self.label_14.setText(
            QCoreApplication.translate("MainWindowDesign", "Name", None)
        )
        self.label_15.setText(
            QCoreApplication.translate("MainWindowDesign", "CofG Position", None)
        )
        self.label_16.setText(QCoreApplication.translate("MainWindowDesign", "x", None))
        self.label_17.setText(QCoreApplication.translate("MainWindowDesign", "y", None))
        self.label_18.setText(QCoreApplication.translate("MainWindowDesign", "Z", None))
        self.label_19.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Other Options go here\n"
                "                                                    ",
                None,
            )
        )
        self.saveAircraftButtom.setText(
            QCoreApplication.translate("MainWindowDesign", "Save Aircraft", None)
        )
        self.groupBox.setTitle("")
        self.label.setText(QCoreApplication.translate("MainWindowDesign", "x", None))
        self.label_2.setText(QCoreApplication.translate("MainWindowDesign", "y", None))
        self.label_11.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Leading Edge Location\n"
                "                                                        ",
                None,
            )
        )
        self.LE_loc_x.setSuffix("")
        self.label_5.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Number of Panels\n"
                "                                                        ",
                None,
            )
        )
        self.mainWingSymmetryCheckbox.setText(
            QCoreApplication.translate("MainWindowDesign", "Symmetric?", None)
        )
        self.label_4.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Panel Spacing Option\n"
                "                                                        ",
                None,
            )
        )
        self.label_12.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Chord-wise Panels\n"
                "                                                        ",
                None,
            )
        )
        self.CWPanelSpcaing.setItemText(
            0, QCoreApplication.translate("MainWindowDesign", "Uniform", None)
        )
        self.CWPanelSpcaing.setItemText(
            1, QCoreApplication.translate("MainWindowDesign", "Cosine", None)
        )
        self.CWPanelSpcaing.setItemText(
            2,
            QCoreApplication.translate(
                "MainWindowDesign",
                "Reverse Cosine\n"
                "                                                            ",
                None,
            ),
        )

        self.label_3.setText(
            QCoreApplication.translate(
                "MainWindowDesign", "<html><head/><body><p>z</p></body></html>", None
            )
        )
        self.MainWingLabel.setText(
            QCoreApplication.translate("MainWindowDesign", "Main Wing", None)
        )
        self.groupBox_8.setTitle("")
        self.label_74.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Number\n"
                "                                                                            of Panels\n"
                "                                                                        ",
                None,
            )
        )
        self.label_78.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "x\n"
                "                                                                        ",
                None,
            )
        )
        self.label_69.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Chord\n"
                "                                                                            and Control\n"
                "                                                                            Surface\n"
                "                                                                        ",
                None,
            )
        )
        self.label_75.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Chord\n"
                "                                                                        ",
                None,
            )
        )
        self.label_79.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Control\n"
                "                                                                            Surface\n"
                "                                                                            Deflection\n"
                "                                                                        ",
                None,
            )
        )
        self.label_73.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "y\n"
                "                                                                        ",
                None,
            )
        )
        self.label_70.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                            Span-wise\n"
                "                                                                            panels\n"
                "                                                                        ",
                None,
            )
        )
        self.label_72.setText(
            QCoreApplication.translate(
                "MainWindowDesign", "<html><head/><body><p>z</p></body></html>", None
            )
        )
        self.label_71.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Panel\n"
                "                                                                            Spacing\n"
                "                                                                            Option\n"
                "                                                                        ",
                None,
            )
        )
        self.label_77.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Leading\n"
                "                                                                            Edge\n"
                "                                                                            Location\n"
                "                                                                        ",
                None,
            )
        )
        self.label_76.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Control\n"
                "                                                                            Surface\n"
                "                                                                            Hinge Point\n"
                "                                                                        ",
                None,
            )
        )
        self.WRCW_LE_loc_x.setSuffix("")
        self.SWPanelSpcaing.setItemText(
            0,
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                                Uniform\n"
                "                                                                            ",
                None,
            ),
        )
        self.SWPanelSpcaing.setItemText(
            1,
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                                Cosine\n"
                "                                                                            ",
                None,
            ),
        )
        self.SWPanelSpcaing.setItemText(
            2,
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                                Reverse\n"
                "                                                                                Cosine\n"
                "                                                                            ",
                None,
            ),
        )

        self.WRCS_S_symCheckbox.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                            Symmetric?\n"
                "                                                                        ",
                None,
            )
        )
        self.MainWingLabel_6.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Root\n"
                "                                                                            Cross\n"
                "                                                                            Section\n"
                "                                                                        ",
                None,
            )
        )
        self.label_23.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                            Aerofoil\n"
                "                                                                        ",
                None,
            )
        )
        self.tabWidget_2.setTabText(
            self.tabWidget_2.indexOf(self.tab_7),
            QCoreApplication.translate(
                "MainWindowDesign", "Wing Root Cross Section", None
            ),
        )
        self.groupBox_7.setTitle("")
        self.label_62.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "y\n"
                "                                                                        ",
                None,
            )
        )
        self.label_60.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Panel\n"
                "                                                                            Spacing\n"
                "                                                                            Option\n"
                "                                                                        ",
                None,
            )
        )
        self.WTCW_LE_loc_x.setSuffix("")
        self.label_63.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Number\n"
                "                                                                            of Panels\n"
                "                                                                        ",
                None,
            )
        )
        self.SWPanelSpcaing_WTCS.setItemText(
            0,
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                                Uniform\n"
                "                                                                            ",
                None,
            ),
        )
        self.SWPanelSpcaing_WTCS.setItemText(
            1,
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                                Cosine\n"
                "                                                                            ",
                None,
            ),
        )
        self.SWPanelSpcaing_WTCS.setItemText(
            2,
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                                Reverse\n"
                "                                                                                Cosine\n"
                "                                                                            ",
                None,
            ),
        )

        self.label_58.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Chord\n"
                "                                                                            and Control\n"
                "                                                                            Surface\n"
                "                                                                        ",
                None,
            )
        )
        self.MainWingLabel_5.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Root\n"
                "                                                                            Cross\n"
                "                                                                            Section\n"
                "                                                                        ",
                None,
            )
        )
        self.label_65.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Control\n"
                "                                                                            Surface\n"
                "                                                                            Hinge Point\n"
                "                                                                        ",
                None,
            )
        )
        self.label_66.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Leading\n"
                "                                                                            Edge\n"
                "                                                                            Location\n"
                "                                                                        ",
                None,
            )
        )
        self.WTCS_S_symCheckbox.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                            Symmetric?\n"
                "                                                                        ",
                None,
            )
        )
        self.label_67.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "x\n"
                "                                                                        ",
                None,
            )
        )
        self.label_64.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Chord\n"
                "                                                                        ",
                None,
            )
        )
        self.label_59.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                            Span-wise\n"
                "                                                                            Panels\n"
                "                                                                        ",
                None,
            )
        )
        self.label_61.setText(
            QCoreApplication.translate(
                "MainWindowDesign", "<html><head/><body><p>z</p></body></html>", None
            )
        )
        self.label_68.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "Control\n"
                "                                                                            Surface\n"
                "                                                                            Deflection\n"
                "                                                                        ",
                None,
            )
        )
        self.label_20.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "\n"
                "                                                                            Aerofoil\n"
                "                                                                        ",
                None,
            )
        )
        self.tabWidget_2.setTabText(
            self.tabWidget_2.indexOf(self.tab_8),
            QCoreApplication.translate(
                "MainWindowDesign", "Wing Tip Cross Section", None
            ),
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.aircraftParametersTab),
            QCoreApplication.translate("MainWindowDesign", "Aircraft Parameters", None),
        )
        self.groupBox_2.setTitle(
            QCoreApplication.translate("MainWindowDesign", "Controller", None)
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.modellingTab),
            QCoreApplication.translate("MainWindowDesign", "Model", None),
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_4),
            QCoreApplication.translate("MainWindowDesign", "Materials", None),
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_5),
            QCoreApplication.translate("MainWindowDesign", "Boundary Conditions", None),
        )
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.tab_6),
            QCoreApplication.translate("MainWindowDesign", "Solution", None),
        )
        self.groupBox_4.setTitle("")
        self.q_start_2.setText(
            QCoreApplication.translate("MainWindowDesign", "Move", None)
        )
        self.q_homing_2.setText(
            QCoreApplication.translate("MainWindowDesign", "Homing", None)
        )
        self.q_x_2.setSuffix("")
        self.q_injection_2.setText(
            QCoreApplication.translate("MainWindowDesign", "Injection", None)
        )
        self.label_21.setText(QCoreApplication.translate("MainWindowDesign", "x", None))
        self.label_22.setText(QCoreApplication.translate("MainWindowDesign", "y", None))
        self.label_26.setText(
            QCoreApplication.translate(
                "MainWindowDesign", "<html><head/><body><p>z</p></body></html>", None
            )
        )
        self.label_27.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "<html><head/><body><p>\u03b2</p></body></html>",
                None,
            )
        )
        self.label_28.setText(QCoreApplication.translate("MainWindowDesign", "d", None))
        self.tabWidget.setTabText(
            self.tabWidget.indexOf(self.visualiserSettingsTab),
            QCoreApplication.translate("MainWindowDesign", "Visualiser Settings", None),
        )
        self.Logo.setText(
            QCoreApplication.translate(
                "MainWindowDesign",
                "<html><head/><body><p><img\n"
                '                                            src="docs/Black_Text_Logo.png"/></p></body></html>\n'
                "                                        ",
                None,
            )
        )
        self.genResultsButton.setText(
            QCoreApplication.translate("MainWindowDesign", "Generate Results", None)
        )
        self.plotVisualsButton.setText(
            QCoreApplication.translate("MainWindowDesign", "Plot Visualisation", None)
        )
        self.terminalOutpotLogo.setText(
            QCoreApplication.translate("MainWindowDesign", "Terminal Output", None)
        )
        self.menuHelp.setTitle(
            QCoreApplication.translate("MainWindowDesign", "Help", None)
        )
        self.menuExamples.setTitle(
            QCoreApplication.translate("MainWindowDesign", "Examples", None)
        )
        self.menuFile.setTitle(
            QCoreApplication.translate("MainWindowDesign", "File", None)
        )
        self.menuAircraft.setTitle(
            QCoreApplication.translate("MainWindowDesign", "Aircraft", None)
        )

    # retranslateUi
