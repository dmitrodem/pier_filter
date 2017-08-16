#!/usr/bin/env python2

import sys
from PyQt5.QtWidgets import (QApplication, QWidget,
                             QHBoxLayout, QVBoxLayout,
                             QPushButton, QLabel, QSpinBox,
                             QTableWidget, QTableWidgetItem)
from PyQt5.Qt import QSizePolicy, QSize, QColor
from PyQt5.QtGui import (QPen, QColor, QFont, QPainter)

from pyqtgraph import PlotWidget

class RowWidget(QWidget):
    def __init__(self):
        super(RowWidget, self).__init__()
        self.initUI()

    def initUI(self):
        length_label = QLabel("Length [mm]")
        self.length_field = QSpinBox()
        width_label = QLabel("Width [mm]")
        self.width_field = QSpinBox()
        hbox = QHBoxLayout()
        hbox.addWidget(length_label)
        hbox.addWidget(self.length_field)
        hbox.addWidget(width_label)
        hbox.addWidget(self.width_field)
        self.setLayout(hbox)

class DrawingWidget(QWidget):
    def __init__(self):
        super(DrawingWidget, self).__init__()
        self.setMinimumSize(100, 100)
        self.show()

    def paintEvent(self, e):
        print "paint event"
        qp = QPainter()
        qp.begin(self)
        self.drawWidget(qp)
        qp.end()

    def drawWidget(self, qp):
        size = self.size()
        w = size.width()
        h = size.height()

        off_w = w * 0.05
        off_h = h * 0.05

        width = w * 0.9
        height = h * 0.9 
        
        qp.setBrush(QColor(100, 200, 0, 255))
        qp.drawRect(off_w, off_h, width, height)
        
class FilterGUI(QWidget):
    def __init__(self):
        super(FilterGUI, self).__init__()
        self.initUI()

    def initUI(self):
        self.drawing_area = DrawingWidget()
        self.plotting_area = PlotWidget()
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(['Width, mm', 'Length, mm'])
        add_section_button = QPushButton("Add &Section")
        add_section_button.clicked.connect(self.addRow)
        vbox_left = QVBoxLayout()
        vbox_left.addWidget(self.drawing_area)
        vbox_left.addWidget(self.plotting_area)

        vbox_right = QVBoxLayout()
        vbox_right.addWidget(self.table)
        vbox_right.addWidget(add_section_button)
        vbox_right.addStretch()
        
        hbox = QHBoxLayout()
        hbox.addLayout(vbox_left)
        hbox.addLayout(vbox_right)

        self.setLayout(hbox)
        self.setWindowTitle("FilterGUI")
        self.show()

    def addRow(self):
        self.table.insertRow(self.table.rowCount())
        self.dumpConfig()

    def dumpConfig(self):
        for row in xrange(self.table.rowCount()):
            items = []
            for column in xrange(self.table.columnCount()):
                item = self.table.item(row, column)
                if hasattr(item, "text"):
                    item = item.text()
                items.append(item)
            print items

def run():
    app = QApplication(sys.argv)
    gui = FilterGUI()
    gui.show()
    sys.exit(app.exec_())

run()

