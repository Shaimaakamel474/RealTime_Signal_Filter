import sys
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget,QVBoxLayout,QCheckBox ,QPushButton, QHBoxLayout
from PyQt5 import QtWidgets, uic, QtGui
from matplotlib.pyplot import figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
import matplotlib.pyplot as plt
from cmath import *
from numpy import *
import numpy as np
import matplotlib
import pyqtgraph as pg
import pandas as pd
import re
# from Main import  plot_frequency_response_phase
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy import signal
matplotlib.use('Qt5Agg')
#____________________________________________________________________________________________________________________#

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=7, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(1, 1, 1)
        super(MplCanvas, self).__init__(fig)
#____________________________________________________________________________________________________________________#


class Create_checkbox_with_button:
    def __init__(self , checkbox_val , layout_widget):
        self.container_widget = QWidget()
        self.checkbox = QCheckBox(checkbox_val)
        self.checkbox.setStyleSheet("""
            color: Black;
            height: 35px; 
            font-family: "MS Shell Dlg 2"; /* Set the font */
            font-size: 9pt; /* Set the font size */
        """)
        self.button = QPushButton('Delete')
        self.button.setStyleSheet("""
            background-color: #0062cc;
            border-radius: 10px;
            border-style: solid;
            border-width: 1px;
            border-color: #0062cc;
            color: white;
            height: 35px; 
            font-family: "MS Shell Dlg 2"; /* Set the font */
            font-size: 9pt; /* Set the font size */
        """)

        self.container_layout = QHBoxLayout()
        self.container_layout.addWidget(self.checkbox)
        self.container_layout.addWidget(self.button)
        self.container_widget.setLayout(self.container_layout)
        layout_widget.addWidget(self.container_widget)
        self.button.clicked.connect(lambda: self.delete_checkbox_with_button(self.container_widget,layout_widget))
        self.checkbox.clicked.connect(lambda:  self.take_checkbox_txt(self.checkbox))
        
    
    def delete_checkbox_with_button(self, container_widget,layout_widget):
        layout_widget.removeWidget(container_widget)
        container_widget.deleteLater()
        checkbox = container_widget.findChild(QCheckBox)
        MyDialog.checked_checkbox_names.remove(checkbox.text())

    def take_checkbox_txt(self,checkbox):
        if checkbox.isChecked():
            MyDialog.checked_checkbox_names.append(checkbox.text())
        else:
            MyDialog.checked_checkbox_names.remove(checkbox.text())
#____________________________________________________________________________________________________________________#

class MyDialog(QtWidgets.QDialog):
    checked_checkbox_names = []
    def __init__(self, main):
        super(MyDialog, self).__init__()
        # Load the UI Page
        uic.loadUi(r'untitled_allPassFilter.ui', self)
        #-----------------creats canvas--------------------------------------


        self.comboBox.addItem("Select an option")
        self.comboBox.model().item(0).setEnabled(False)
        self.Filter_graph.setBackground('w')
        self.phase_filter_graph.setBackground('w')

        self.canvas_phase = MplCanvas(self, width=10, height=10, dpi=100)
        self.canvas_phase.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.layout_phase = QtWidgets.QVBoxLayout()
        self.layout_phase.addWidget(self.canvas_phase)
        self.layout_phase.setContentsMargins(0, 0, 0, 0)
        self.graph = pg.PlotItem() 
        self.canvas1 = MplCanvas(self, width=10, height=10, dpi=100)
        self.canvas1.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.layout1 = QtWidgets.QVBoxLayout()
        self.layout1.addWidget(self.canvas1)
        self.layout1.setContentsMargins(0, 0, 0, 0)
        self.main =main

        self.pushButton_4.clicked.connect(lambda: self.main.plot_frequency_response_phase(flag = True))

        self.canvas2 = MplCanvas(self, width=5, height=4, dpi=100)
        self.layout2 = QtWidgets.QVBoxLayout()
        self.layout2.addWidget(self.canvas2)

        # self.a = []
        # self.b = []
        #---------------initializations-------------------------------------------
        self.file_path= r"Filters.csv"
        self.read_csv_filters(self.file_path)
        self.all_pass_filter_csv_real= []
        self.all_pass_filter_csv_imag= []
        self.add_complex_to_lib_from_combo = None
        self.selected_mode_str=None
        self.lineEdit_val = None
        self.angle_filtered = None
        self.a = []
        self.b = []
        #------------------connections-------------------------------------------
        self.comboBox.activated.connect(lambda: self.get_combobox_value())
        self.Done_btn.clicked.connect(self.get_lineEdit_values)
        self.Add_filter_combo_btn.clicked.connect(  lambda : self.add_to_lib(self.selected_mode_str,self.added_filters))
        self.Add_filter_from_line_btn.clicked.connect(lambda: self.add_to_lib(self.lineEdit_val,self.added_filters))
        self.Apply_btn.clicked.connect(lambda: self.plot_final_phase_gain())

    #------------------------------------------------------------------------------------------------------------------------
    def plot_canvas(self, canvas, widget, layout, x,y ):
        canvas.axes.cla()
        canvas.axes.plot(x, y)
        canvas.draw()
        widget.setCentralItem(self.graph)
        widget.setLayout(layout)

# function to get all pass filters coef from a csv file
    def read_csv_filters (self, file_path):
        data = pd.read_csv(file_path)
        self.all_pass_filter_csv_real = data['real']
        self.all_pass_filter_csv_imag = data['imag']
        self.write_on_combobox()

# fill the combobox  with csv values 
    def write_on_combobox(self):
        for first, second in zip(self.all_pass_filter_csv_real, self.all_pass_filter_csv_imag):
            complex_number = complex(first, second)
            string_complex_number = str(complex_number)
            self.add_complex_to_lib_from_combo = string_complex_number
            sanitized_complex_number = string_complex_number.translate(str.maketrans('', '', '()'))
            self.comboBox.addItem(sanitized_complex_number)

# displaying a filter when selecting a combobox value 
    def get_combobox_value(self):
        self.selected_mode_str = self.comboBox.currentText()  
        complex_number_list = self.from_str_to_float(self.selected_mode_str)
        self.plot_all_pass_filter(complex_number_list)

# function to create a libirary
    def add_to_lib(self, checkbox , widget):
        Create_checkbox_with_button(checkbox,widget)
        print("Checked Checkboxes:", self.checked_checkbox_names)

# displaying a filter when write on a lineEdite
    def get_lineEdit_values(self):
        real = float(self.Real_line.text())
        imag = float(self.Img_line.text())
        complex_value = complex(real , imag)
        complex_value = str(complex_value)
        self.lineEdit_val = complex_value
        complex_number_list = self.from_str_to_float(complex_value)
        self.plot_all_pass_filter(complex_number_list)

# function to create all_pass_filter
    def phase_response(self,a):
        a= self.list_of_lists_to_complex(a)
        self.b = [-np.conj(a), 1.0]
        self.a = [1.0, -a]
        self.x_b = [-np.conj(a)]
        self.y_a = [-a]
        self.x_b_complex = [complex(x) for x in self.x_b]
        self.y_a_complex = [complex(y) for y in self.y_a]
        w, h = signal.freqz(self.b, self.a)
        angels = np.zeros(512) if a==1 else np.unwrap(np.angle(h))
        w = w/max(w)
        h = np.around(angels, decimals=3)
        return w , h 
    
# plot final_phase_gain 
    def plot_final_phase_gain(self):
        filter_angles = np.zeros(512)
        w = np.zeros(512)
        if self.checked_checkbox_names :
            for a in self.checked_checkbox_names :
                a = self.from_str_to_float(a)
                w , angeles = self.phase_response(a)
                filter_angles = np.add(filter_angles, angeles)
                self.angle_filtered = filter_angles
                pass
            self.plot_canvas(self.canvas2 , self.phase_filter_graph,
                        self.layout2, w , filter_angles  )
            


    # plot a choosen all-pass filter from either a combobox or lineEdite
    def plot_all_pass_filter(self,a):
        w , h = self.phase_response(a)
        self.plot_canvas(self.canvas1 , self.Filter_graph,
                    self.layout1, w , h  )
        
#--------------------utility functions ----------------------
    def list_of_lists_to_complex(self,complex_data):
        if isinstance(complex_data[0], list):
            complex_numbers = [complex(r, imag) for r, imag in complex_data]
        else:
            complex_numbers = complex(complex_data[0], complex_data[1])
        return complex_numbers
    
    def from_str_to_float(self,selected_mode_str):
        minus_sign_count = selected_mode_str.count('-')
        plus_sign_count = selected_mode_str.count('+')
        parts = re.split(pattern=r"([-()+j])", string=selected_mode_str)
        numeric_parts = [part for part in parts if part and part not in {'-', '(', '+', ')', 'j'}]
        if minus_sign_count == 2:
            numeric_floats = [-float(part) for part in numeric_parts]
        else:
            if minus_sign_count == 1 and plus_sign_count:
                numeric_floats = [float(part) for part in numeric_parts]
                numeric_floats[0]=  -abs(numeric_floats[0])
            else:
                numeric_floats = [float(part) for part in numeric_parts]
                numeric_floats[-1]=  -abs(numeric_floats[-1]) 
        
        if minus_sign_count == 0:
            numeric_floats = [float(part) for part in numeric_parts]

        if len(numeric_floats) == 1:
            numeric_floats.insert(0, 0.0)
        return numeric_floats

def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MyDialog(QtWidgets.QDialog)
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()