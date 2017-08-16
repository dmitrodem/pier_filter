#!/usr/bin/env python2

import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QMainWindow,
                             QPushButton, QHBoxLayout, QVBoxLayout)

class MyWidget(QWidget):
    def __init__(self):
        super(MyWidget, self).__init__()
        self.initUI()
        
    def initUI(self):
        self.resize(250, 150)
        self.move(300, 300)
        self.setWindowTitle('Simple example')
        okButton = QPushButton('OK')
        cancelButton = QPushButton('Cancel')
        hbox = QHBoxLayout()
        hbox.addStretch(1)
        hbox.addWidget(okButton)
        hbox.addWidget(cancelButton)

        vbox = QVBoxLayout()
        vbox.addStretch(1)
        vbox.addLayout(hbox)

        self.setLayout(vbox)
        self.show()

def run():
    app = QApplication(sys.argv)
    w = MyWidget()
    sys.exit(app.exec_())


run()
