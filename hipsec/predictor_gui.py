#!/usr/bin/env python

'''
Hipsec predicts if gene overexpression of an extracellular protein will lead 
to successful high-level production in Aspergillus niger. The prediction is 
based on a protein's amino acid sequence. More information about the method
can be found in:

Exploring sequence characteristics related to high-level production of secreted 
proteins in Aspergillus niger. B.A. van den Berg, M.J.T. Reinders, M. Hulsman, 
L. Wu, H.J. Pel, J.A. Roubos, D. de Ridder (2012) PLoS ONE (accepted).

Copyright (C) 2013  B.A. van den Berg

Hipsec is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>. 
'''

import os
import sys
from PyQt4.QtCore import QSize
from PyQt4.QtGui import QApplication 
from PyQt4.QtGui import QMainWindow
from PyQt4.QtGui import QIcon
from PyQt4.QtGui import QTextEdit
from PyQt4.QtGui import QAction
from PyQt4.QtGui import QFileDialog
from PyQt4.QtGui import QTextCursor
from PyQt4.QtGui import qApp
import predictor

__author__ = "Bastiaan van den Berg"
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Bastiaan van den Berg"
__email__ = "b.a.vandenberg@tudelft.nl"

class PredictorGUI(QMainWindow):
    
    def __init__(self):
        super(PredictorGUI, self).__init__()
        self.default_location = os.getenv("HOME")
        self.current_load_location = self.default_location
        self.current_save_location = self.default_location
        self.initUI()
        
    def initUI(self):      

        self.textEdit = QTextEdit()
        self.textEdit.setReadOnly(True)
        
        self.setCentralWidget(self.textEdit)
        self.statusBar()

        icon_open = QIcon.fromTheme("document-open", QIcon(":/open.png"))
        icon_save = QIcon.fromTheme("document-save-as")
        icon_copy = QIcon.fromTheme("edit-copy")
        icon_clear = QIcon.fromTheme("edit-clear")

        openFile = QAction(QIcon(icon_open), 'Open', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new File')
        openFile.triggered.connect(self.showLoadDialog)

        self.saveFile = QAction(QIcon(icon_save), 'Save', self)
        self.saveFile.setShortcut('Ctrl+S')
        self.saveFile.setStatusTip('Save data to File')
        self.saveFile.triggered.connect(self.showSaveDialog)
        self.saveFile.setDisabled(True)

        self.clearData = QAction(QIcon(icon_clear), 'Clear', self)
        self.clearData.setShortcut('Ctrl+D')
        self.clearData.setStatusTip('Clear data')
        self.clearData.triggered.connect(self.clear_event)
        self.clearData.setDisabled(True)

        self.copyData = QAction(QIcon(icon_copy), 'Copy', self)
        self.copyData.setShortcut('Ctrl+C')
        self.copyData.setStatusTip('Copy data to clipboard')
        self.copyData.triggered.connect(self.copy_event)
        self.copyData.setDisabled(True)

        exitAction = QAction(QIcon(), 'Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.triggered.connect(qApp.quit)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)
        fileMenu.addAction(self.saveFile)   
        fileMenu.addAction(self.copyData)
        fileMenu.addAction(self.clearData)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)

        fileToolBar = self.addToolBar("Load and Save")
        fileToolBar.setIconSize(QSize(22, 22))
        fileToolBar.addAction(openFile)
        fileToolBar.addAction(self.saveFile)
        fileToolBar.addAction(self.copyData)
        fileToolBar.addAction(self.clearData)

        self.setGeometry(300, 300, 400, 600)
        self.setWindowTitle('File dialog')
        self.show()
    
    def showLoadDialog(self):

        fname = QFileDialog.getOpenFileName(self, 'Open file',
                    self.current_load_location)
        self.current_load_location = os.path.dirname(str(fname))
        data = predictor.get_output_string(fname)                        
        self.textEdit.setText(data)
        self.saveFile.setDisabled(False)
        self.clearData.setDisabled(False)
        self.copyData.setDisabled(False)
        
    def showSaveDialog(self):

        if(self.current_save_location == self.default_location):
            self.current_save_location = self.current_load_location

        fname = QFileDialog.getSaveFileName(self, 'Save file', 
                    self.current_save_location)
        self.current_save_location = os.path.dirname(str(fname))
        
        with open(str(fname), 'w') as fout:
            fout.write(self.textEdit.toPlainText())
    
    def clear_event(self):
        self.textEdit.setText('')
        self.saveFile.setDisabled(True)
        self.clearData.setDisabled(True)
        self.copyData.setDisabled(True)
        
    def copy_event(self):
        self.textEdit.selectAll()
        self.textEdit.copy()
        self.textEdit.moveCursor(QTextCursor.Start)

def main():
    
    app = QApplication(sys.argv)
    ex = PredictorGUI()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
