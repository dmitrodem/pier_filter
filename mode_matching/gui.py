#!/usr/bin/env python2

import sys
from PyQt5.QtWidgets import (QApplication, QWidget,
                             QHBoxLayout, QVBoxLayout,
                             QPushButton, QLabel, QSpinBox,
                             QTableWidget, QTableWidgetItem,
                             QMenu, QAction)
from PyQt5.Qt import QSizePolicy, QSize, QColor, QItemSelection
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from PyQt5.QtGui import (QPen, QColor, QFont, QPainter)

from pyqtgraph import PlotWidget

class DrawingWidget(QWidget):
    def __init__(self):
        super(DrawingWidget, self).__init__()
        self.setMinimumSize(200, 200)
        self.show()

    def paintEvent(self, e):
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

class TableWidget(QTableWidget):
    def __init__(self, rows, columns):
        super(TableWidget, self).__init__(rows, columns)

        self.insertAboveAction = QAction("Insert Row &Above")
        self.insertAboveAction.triggered.connect(self.insertAbove)
        
        self.insertBelowAction = QAction("Insert Row &Below")
        self.insertBelowAction.triggered.connect(self.insertBelow)

        self.deleteAction = QAction("&Delete Selection")
        self.deleteAction.triggered.connect(self.deleteRows)

    def contextMenuEvent(self, event):
        menu = QMenu()
        menu.addAction(self.insertAboveAction)
        menu.addAction(self.insertBelowAction)
        menu.addAction(self.deleteAction)
        action = menu.exec_(self.mapToGlobal(event.pos()))

    def getSelectedRows(self):
        indices = self.selectionModel().selection().indexes()
        rows = []
        for i in indices:
            rows.append(i.row())
        rows = list(set(rows))

        return rows
    
    def insertAbove(self):
        rows = self.getSelectedRows()
        if rows:
            self.insertRow(rows[0])
        else:
            self.insertRow(0)
            
    def insertBelow(self):
        rows = self.getSelectedRows()
        if rows:
            self.insertRow(rows[-1])
        else:
            self.insertRow(self.rowCount())
            
    def deleteRows(self):
        rows = self.getSelectedRows()
        for r in rows:
            self.removeRow(r)

class FilterGUI(QWidget):
    valueChanged = pyqtSignal([int], ['QString'])

    @pyqtSlot()
    def foo(self):
        print "simple signal"

    def __init__(self):
        super(FilterGUI, self).__init__()
        self.valueChanged.connect(self.foo)
        self.valueChanged.emit('')
        self.initUI()

    def initUI(self):
        self.drawing_area = DrawingWidget()
        self.plotting_area = PlotWidget()
        #self.table = QTableWidget(0, 2)
        self.table = TableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(['Width, mm', 'Length, mm'])
        add_section_button = QPushButton("Add &Section")
        add_section_button.clicked.connect(self.addRow)
        calculate_button = QPushButton("&Calculate");
        calculate_button.clicked.connect(self.dumpConfig)
        vbox_left = QVBoxLayout()
        vbox_left.addWidget(self.drawing_area)
        vbox_left.addWidget(self.plotting_area)

        vbox_right = QVBoxLayout()
        vbox_right.addWidget(self.table)
        vbox_right.addWidget(add_section_button)
        vbox_right.addWidget(calculate_button)
        vbox_right.addStretch()
        
        hbox = QHBoxLayout()
        hbox.addLayout(vbox_left)
        hbox.addLayout(vbox_right)

        self.setLayout(hbox)
        self.setWindowTitle("FilterGUI")
        self.show()

    def addRow(self):
        self.table.insertRow(self.table.rowCount())
        

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

