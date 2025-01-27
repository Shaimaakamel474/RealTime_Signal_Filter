from matplotlib.figure import Figure
from PyQt5.QtCore import Qt
from PyQt5 import QtCore
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QGraphicsScene ,QLabel , QHBoxLayout
from PyQt5 import QtWidgets, uic, QtGui
from matplotlib.pyplot import figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from PyQt5.QtWidgets import QWidget
import matplotlib.pyplot as plt
from cmath import *
from numpy import *
from scipy.signal import tf2zpk
import sys
from PyQt5.QtWidgets import QDialog
from scipy import signal
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyqtgraph as pg
import pandas as pd
matplotlib.use('Qt5Agg')
import numpy as np
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import simpleSplit
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from AllPass import MyDialog
import csv
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=7, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(1, 1, 1)
        super(MplCanvas, self).__init__(fig)
    def clear_all_plots(self):
        self.axes.cla()  # Clear the axes
        self.draw()  # Redraw the canvas to reflect the changes
       



class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Load the UI Page
        uic.loadUi(r'MainWindow - untitled.ui', self)
        self.graph = pg.PlotItem()
        self.setup_layouts()
        self.Widget_Touchpad.setMouseTracking(True)
        self.Widget_Touchpad.mouseMoveEvent = self.calculateMousePosition

        self.load_button.clicked.connect(self.load_signal)  

        self.Delete_combox.currentIndexChanged.connect(self.Delete_Options)

        self.comboBox_swap.currentIndexChanged.connect(lambda indx : self.swap_markers_and_positions(indx))
        self.button_Save.clicked.connect(self.save_poles_zeros_to_csv)
        self.Button_Export.clicked.connect(self.open_file_dialog)

        # self.Save_Method_Combox.currentIndexChanged.connect(self.save_and_reload_filter)

        self.Combox_Filters.currentIndexChanged.connect(lambda : self.Display_Filters())
        self.Combox_Filters_type.currentIndexChanged.connect(lambda : self.Display_Filters())
        
        self.Second_Window.clicked.connect(self.Open_Second_Window)


        self.Combox_Generat_C_code.currentIndexChanged.connect(lambda : self.realize_and_export())
                # Plot unit circle
        self.plot_unit_circle()
        
        self.zeros = []     
        self.poles = []

        self.zeros_postions = []
        self.poles_postions = []

        self.new_zeros = []
        self.new_poles = []
        self.filtered_signal = [0]
        self.lowpass = False

        self.last_cursor_pos = None
        self.mag_list = [0]
        self.time_list = [0]

        # Track the selected marker for dragging
        self.selected_marker = None
        self.dragging_marker = None
        self.selected_primary_marker = None
        self.selected_conjugate_marker = None
        self.offset = (0, 0)

        self.canvas1.mpl_connect('button_press_event', self.on_press)
        self.canvas1.mpl_connect('motion_notify_event', self.on_motion)
        self.canvas1.mpl_connect('button_release_event', self.on_release)
       
        self.plot_frequency_response()




    def setup_layouts(self):

        self.Real_Signal.setBackground('w')
        self.Filtered_Signal.setBackground('w')
        self.canvas1 = MplCanvas(self, width=10, height=10, dpi=100)
        self.canvas1.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.layout1 = QtWidgets.QVBoxLayout()
        self.layout1.addWidget(self.canvas1)



        self.canvas_mag = MplCanvas(self, width=5, height=4, dpi=100)
        self.canvas1.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.layout_mag = QtWidgets.QVBoxLayout()
        self.layout_mag.addWidget(self.canvas_mag)

        self.canvas_phase = MplCanvas(self, width=10, height=10, dpi=100)
        self.canvas1.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.layout_phase = QtWidgets.QVBoxLayout()
        self.layout_phase.addWidget(self.canvas_phase)
        
        
                        # Set the layout for the phase graph
        self.Phase_Graph.setCentralItem(self.graph)
        self.layout_phase.setContentsMargins(0, 0, 0, 0)
        self.Phase_Graph.setLayout(self.layout_phase)


              # Set the layout for the magnitude graph
        self.Mag_Graph.setCentralItem(self.graph)
        self.layout_mag.setContentsMargins(0, 0, 0, 0)
        self.Mag_Graph.setLayout(self.layout_mag)


        self.circle_graph.setCentralItem(self.graph)
        self.layout1.setContentsMargins(0, 0, 0, 0)
        self.circle_graph.setLayout(self.layout1)
        




    def plot_unit_circle(self):
        # Define the angles for the unit circle
        angles = np.linspace(0, 2 * np.pi, 1000)
        x = np.cos(angles)
        y = np.sin(angles)

        # Plot the unit circle
        self.canvas1.axes.plot(x, y, 'k')  # Black line for the unit circle

        # Plot the x-axis (horizontal line)
        self.canvas1.axes.axhline(0, color='k')

        # Plot the y-axis (vertical line)
        self.canvas1.axes.axvline(0, color='k')

        # Number of angular divisions
        angular_divisions = 24  # 24 sectors, 15° apart

        # Number of radial divisions
        radial_divisions = 8  # Concentric circles starting further from the center

        # Plot angular lines (radial lines)
        for theta in np.linspace(0, 2 * np.pi, angular_divisions, endpoint=False):
            x_line = [0, np.cos(theta)]  # Scale radial lines to match the max radius (1)
            y_line = [0, np.sin(theta)]
            self.canvas1.axes.plot(x_line, y_line, 'k', linewidth=0.6)  # Black solid lines

        # Plot radial gridlines (concentric circles)
        # Start the circles from 0.2 instead of 0 to avoid clutter in the center
        radii = np.linspace(0.2, 1, radial_divisions)
        for radius in radii:
            x_circle = radius * np.cos(angles)
            y_circle = radius * np.sin(angles)
            self.canvas1.axes.plot(x_circle, y_circle, 'k', linewidth=0.6)  # Black solid lines

        
        #         # Add numbers to the unit circle
        # for theta in np.linspace(0, 2 * np.pi, angular_divisions, endpoint=False):
        #     x_label = 1.1 * np.cos(theta)  # Position labels slightly outside the circle
        #     y_label = 1.1 * np.sin(theta)
        #     angle_deg = int(np.degrees(theta))  # Convert angle to degrees
        #     self.canvas1.axes.text(
        #         x_label, y_label, f"{angle_deg/360:.2f}", ha='center', va='center', fontsize=
        #     )
            
        
        # Hide the plot widget's axes
        self.circle_graph.setAxisItems({})
        self.circle_graph.setRange(QtCore.QRectF(-1.1, -1.1, 2.2, 2.2))

        # Set the canvas background color to white
        self.canvas1.setStyleSheet("background-color: white;")

        # Resize the canvas to fill the entire plot widget
        self.canvas1.setGeometry(self.circle_graph.geometry())

        # Set the limits to ensure the origin is at the center of the circle
        self.canvas1.axes.set_xlim(-1.1, 1.1)
        self.canvas1.axes.set_ylim(-1.1, 1.1)

        # Refresh the canvas to display the updated plot
        self.canvas1.draw()



    def plot_frequency_response(self):
        # Convert zeros and poles to complex numbers
        zeros = [complex(x, y) for x, y in self.zeros_postions]
        poles = [complex(x, y) for x, y in self.poles_postions]
        
        # Create frequency response using freqz_zpk from scipy
        w, h = signal.freqz_zpk(zeros, poles, 1.0)  # Frequency and response
        
        # Magnitude Response
        self.canvas_mag.axes.set_facecolor('white')
        self.canvas_mag.axes.clear()  # Clear the previous plot
        self.canvas_mag.axes.semilogx(w, abs(h))  # Magnitude in dB
        self.canvas_mag.axes.set_xlabel('Frequency (Hz)')
        self.canvas_mag.axes.set_ylabel('Magnitude (dB)')
        self.canvas_mag.axes.set_title('Magnitude Response')
        self.canvas_mag.axes.grid(True)
        self.canvas_mag.axes.set_xlim(1e-2, 10)  # Log scale for frequency
        self.canvas_mag.draw()  # Redraw the canvas
        
  
        
        # Phase Response
        self.canvas_phase.axes.clear()  # Clear previous phase plot
        self.canvas_phase.axes.semilogx(w, np.angle(h, deg=True))  # Phase in degrees
        self.canvas_phase.axes.set_xlabel('Frequency (Hz)')
        self.canvas_phase.axes.set_ylabel('Phase (degrees)')
        self.canvas_phase.axes.set_title('Phase Response')
        self.canvas_phase.axes.grid(True)
        self.canvas_phase.axes.set_xlim(1e-2, 10)  # Log scale for frequency
        self.canvas_phase.draw()  # Redraw the phase canvas
        



    # to handle all actions on the unit circle >> move , delete specific object , draw new object
    def on_press(self , event ):
        if event.button == 1:  # Left mouse button
            if event.inaxes == self.canvas1.axes:
                x, y = event.xdata, event.ydata
                if x is not None and y is not None:
                    for primary_marker, conjugate_marker in self.zeros + self.poles:
                        if primary_marker.contains(event)[0]:
                            self.selected_primary_marker = primary_marker
                            self.selected_conjugate_marker = conjugate_marker
                            self.offset = (primary_marker.get_offsets()[0][0] - x, primary_marker.get_offsets()[0][1] - y)
                            self.dragging_marker = True  # Set dragging mode
                            break
                    
                    
                    else:self.draw_New_object(x,y)

        elif event.button == 3:  # Right mouse button for deletion
            self.delete_selected(event)        
                    

#draw new zero , pole in the pressed position
    def draw_New_object(self , x,y):
        if self.Button_Zero.isChecked():
            if self.add_conjugates.isChecked():
                if y != 0:
                    # Add conjugate zero to the list
                    conjugate_marker = self.canvas1.axes.scatter(x, -y, color='b', marker='o', s=100,
                                                                facecolors='none', edgecolors='b',
                                                                linewidths=2)
                    primary_marker = self.canvas1.axes.scatter(x, y, color='b', marker='o', s=100,
                                                            facecolors='none', edgecolors='b',
                                                            linewidths=2)
                    self.zeros.append((primary_marker, conjugate_marker))
                    self.zeros_postions.append((x, y))
            else:
                # Add zero to the list
                marker = self.canvas1.axes.scatter(x, y, color='b', marker='o', s=100,
                                                facecolors='none', edgecolors='b', linewidths=2)
                self.zeros.append((marker, None))
                self.zeros_postions.append((x, y))
        
        elif self.Button_Pole.isChecked():
            if self.add_conjugates.isChecked():
                if y != 0:
                    # Add conjugate pole to the list
                    conjugate_marker = self.canvas1.axes.scatter(x, -y, color='r', marker='x', s=100,
                                                                linewidths=2)
                    primary_marker = self.canvas1.axes.scatter(x, y, color='r', marker='x', s=100,
                                                            linewidths=2)
                    self.poles.append((primary_marker, conjugate_marker))
                    self.poles_postions.append((x, y))
            else: 
                # Add pole to the list
                marker = self.canvas1.axes.scatter(x, y, color='r', marker='x', s=100, linewidths=2)
                self.poles.append((marker, None))
                self.poles_postions.append((x, y))
        self.canvas1.draw()
        self.plot_frequency_response()


    def delete_selected(self, event):
        x, y = event.xdata, event.ydata

        if x is not None and y is not None:
            clicked_position = np.array([x, y])
            min_distance = float('inf')
            closest_marker = None

            all_markers = self.zeros + self.poles  # Combine both zero and pole markers
            # print(all_markers)
            # loop in all poles , zeros to get the cloest obj to the pressed location

            for primary_marker, conjugate_marker in all_markers:
                
                if primary_marker is not None:
                    marker_position = primary_marker.get_offsets()[0]

                    distance = np.linalg.norm(marker_position - clicked_position)
                    if distance < min_distance:
                        min_distance = distance
                        
                        closest_marker = (primary_marker, conjugate_marker)

            # delete the cloest obj from our data 
            if closest_marker is not None:
                if closest_marker in self.zeros:
                    self.zeros.remove(closest_marker)
                    self.zeros_postions = [(primary_marker.get_offsets()[0][0], primary_marker.get_offsets()[0][1]) for primary_marker, _ in self.zeros]
                elif closest_marker in self.poles:
                    self.poles.remove(closest_marker)
                    self.poles_postions = [(primary_marker.get_offsets()[0][0], primary_marker.get_offsets()[0][1]) for primary_marker, _ in self.poles]

                closest_marker[0].remove()  # Remove primary marker

                if closest_marker[1] is not None:
                    closest_marker[1].remove()  # Remove conjugate marker if it exists

        self.canvas1.draw()
        self.plot_frequency_response()
        # self.plot_frequency_response_phase(self.canvas_phase, self.phase_responce_graph, self.layout_phase, flag=False)
    
 
# update the graph after delete the all poles ,  all zeros , both
    def redraw_plot(self, delete_zeros=False, delete_poles=False):
        # Redraws the plot after removing zeros or poles
        self.canvas1.axes.clear()
        self.plot_unit_circle()
     
        
        for zero in self.zeros:
            if not delete_zeros:
                # Unpack tuple if necessary and add the artist object
                if isinstance(zero, tuple):
                    self.canvas1.axes.add_artist(zero[0])  # Assuming the artist object is at index 0
                    if zero[1] is not None :self.canvas1.axes.add_artist(zero[1])  #check if there conj
        for pole in self.poles:
            if not delete_poles:
                # Unpack tuple if necessary and add the artist object
                if isinstance(pole, tuple):
                    self.canvas1.axes.add_artist(pole[0])  # Assuming the artist object is at index 0
                    if pole[1] is not None :self.canvas1.axes.add_artist(pole[1])  #check if there conj
        
        self.canvas1.draw()
        self.plot_frequency_response()

    def delete_zeros(self):
        # Delete all zeros on the plot      
        self.zeros = []
        self.zeros_postions = []
        self.redraw_plot(delete_zeros=True)
        self.plot_frequency_response()
       

    def delete_poles(self):
        # Delete all poles on the plot  
        self.poles = []
        self.poles_postions = []
        self.redraw_plot(delete_poles=True)
        self.plot_frequency_response()
       

    def delete_all(self):
        # Delete all poles and zeros on the plot
        self.zeros = []
        self.zeros_postions = []
        self.poles = []
        self.poles_postions = []
        self.redraw_plot()
        self.plot_frequency_response()
    

    def Delete_Options(self):
        state=self.Delete_combox.currentText()
        if state=="ALL_Poles":
            self.delete_poles()
        elif state=="ALL_Zeros":
            self.delete_zeros()
        elif state=="Delete Both":
            self.delete_all()
        self.Delete_combox.setCurrentText("Delete")



    def swap_markers_and_positions(self , index):
        if index == 1 and len(self.poles) > 0:
            zeros, zeros_postions = self.convert_poles_to_zeros()
            self.zeros +=zeros
            self.zeros_postions+=zeros_postions
            self.poles, self.poles_postions=[],[]

        elif index == 2 and len(self.zeros) > 0:
            poles, poles_postions = self.convert_zeros_to_poles()
            self.poles+=poles
            self.poles_postions+=poles_postions
            
            self.zeros, self.zeros_postions =[],[]
        elif index == 3 and (len(self.zeros) > 0 or len(self.poles) > 0):
            zeros, zeros_postions = self.convert_poles_to_zeros()
            poles, poles_postions = self.convert_zeros_to_poles()
            self.zeros, self.zeros_postions =zeros, zeros_postions
            self.poles, self.poles_postions=poles, poles_postions
        self.comboBox_swap.setCurrentIndex(0)
        self.canvas1.draw()  # Redraw the canvas with updated markers
        self.plot_frequency_response()

    def convert_poles_to_zeros(self):
                # Clear existing markers
        for obj in self.poles :
            obj[0].remove()
            if obj[1] is not None:
                obj[1].remove()
        New_Zeros, New_Zeros_Pos = [], []
        for indx, (_, conj) in enumerate(self.poles):            
            x, y = self.poles_postions[indx] if indx < len(self.poles_postions) else self.poles_postions[-1]
            
            # Plot conjugate marker if it exists
            conjugate_marker = None
            if conj is not None:
                conjugate_marker = self.canvas1.axes.scatter(x, -y, color='b', marker='o', s=100, facecolors='none', edgecolors='b', linewidths=2)
            
            primary_marker = self.canvas1.axes.scatter(x, y, color='b', marker='o', s=100, facecolors='none', edgecolors='b', linewidths=2)
            
            New_Zeros.append((primary_marker, conjugate_marker))
            New_Zeros_Pos.append((x, y))
        return New_Zeros, New_Zeros_Pos

    def convert_zeros_to_poles(self):

        for obj in self.zeros :
            obj[0].remove()
            if obj[1] is not None:
                obj[1].remove()
        New_Poles, New_poles_pos = [], []
        for indx, (_, conj) in enumerate(self.zeros):
            x, y = self.zeros_postions[indx] if indx < len(self.zeros_postions) else self.zeros_postions[-1]
            
            # Plot conjugate marker if it exists
            conjugate_marker = None
            if conj is not None:
                conjugate_marker = self.canvas1.axes.scatter(x, -y, color='r', marker='x', s=100, linewidths=2)
            
            primary_marker = self.canvas1.axes.scatter(x, y, color='r', marker='x', s=100, linewidths=2)
            
            New_Poles.append((primary_marker, conjugate_marker))
            New_poles_pos.append((x, y))

        return New_Poles, New_poles_pos
    



            
    def save_poles_zeros_to_csv(self):
        """
        Save poles and zeros data to a CSV file with a dynamic filename.
        
        The filename format is: Filter_{num_poles}Poles_{num_zeros}Zeros.csv
        Stores x, y coordinates for poles and zeros
        """
        filename = f"Filter_{len(self.poles)}Poles_{len(self.zeros)}Zeros.csv"
        # filename="tst.csv"
        with open(filename, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            
            # Write header
            csvwriter.writerow(['Type', 'X', 'Y', 'Is_Conjugate'])
            
            # Track added conjugate points to avoid duplicates
            added_conjugates = set()
            
            # Write poles data
            for i, (pole, conjugate) in enumerate(self.poles):
                x, y = self.poles_postions[i]
                
                # Add primary point
                csvwriter.writerow(['Pole', x, y, conjugate is not None])
                
                # Add conjugate point if it exists and hasn't been added before
                if conjugate is not None and (x, -y) not in added_conjugates:
                    csvwriter.writerow(['Pole', x, -y, True])
                    added_conjugates.add((x, -y))
            
            # Reset for zeros
            added_conjugates.clear()
            
            # Write zeros data
            for i, (zero, conjugate) in enumerate(self.zeros):
                x, y = self.zeros_postions[i]
                
                # Add primary point
                csvwriter.writerow(['Zero', x, y, conjugate is not None])
                
                # Add conjugate point if it exists and hasn't been added before
                if conjugate is not None and (x, -y) not in added_conjugates:
                    csvwriter.writerow(['Zero', x, -y, True])
                    added_conjugates.add((x, -y))
        
        print(f"Saved data to {filename}")
      



    def open_file_dialog(self):
        # Open the file dialog and get the file path
        file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "CSV Files (*.csv);;All Files (*)")
        
        if file_path:
            self.load_poles_zeros_from_csv(file_path)
    

    def load_poles_zeros_from_csv(self, filename):
        """
        Load poles and zeros data from a CSV file and replot on the canvas.
        
        Clears existing plots before creating new ones.
        """
        # Clear existing plots
        self.canvas1.axes.clear()
        self.plot_unit_circle()
        
        # Reset lists
        self.poles = []
        self.poles_postions = []
        self.zeros = []
        self.zeros_postions = []
        
        # Track added conjugate points to avoid duplicates
        added_conjugates = set()
        
        # Read CSV file
        with open(filename, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            next(csvreader)  # Skip header
            
            for row in csvreader:
                point_type, x, y, is_conjugate = row
                x = float(x)
                y = float(y)
                is_conjugate = is_conjugate == 'True'
                
                # Skip if this conjugate point has already been added
                if is_conjugate and (x, y) in added_conjugates:
                    continue
                
                if point_type == 'Pole':
                    if is_conjugate:
                        primary_marker = self.canvas1.axes.scatter(x, y, color='r', marker='x', s=100, linewidths=2)
                        conjugate_marker = self.canvas1.axes.scatter(x, -y, color='r', marker='x', s=100, linewidths=2)
                        self.poles.append((primary_marker, conjugate_marker))
                        self.poles_postions.append((x, y))
                        
                        # Mark both points as added to prevent duplicates
                        added_conjugates.add((x, y))
                        added_conjugates.add((x, -y))
                    else:
                        marker = self.canvas1.axes.scatter(x, y, color='r', marker='x', s=100, linewidths=2)
                        self.poles.append((marker, None))
                        self.poles_postions.append((x, y))
                
                elif point_type == 'Zero':
                    if is_conjugate:
                        primary_marker = self.canvas1.axes.scatter(x, y, color='b', marker='o', s=100, 
                                                                facecolors='none', edgecolors='b', linewidths=2)
                        conjugate_marker = self.canvas1.axes.scatter(x, -y, color='b', marker='o', s=100, 
                                                                    facecolors='none', edgecolors='b', linewidths=2)
                        self.zeros.append((primary_marker, conjugate_marker))
                        self.zeros_postions.append((x, y))
                        
                        # Mark both points as added to prevent duplicates
                        added_conjugates.add((x, y))
                        added_conjugates.add((x, -y))
                    else:
                        marker = self.canvas1.axes.scatter(x, y, color='b', marker='o', s=100, 
                                                           facecolors='none', edgecolors='b', linewidths=2)
                        self.zeros.append((marker, None))
                        self.zeros_postions.append((x, y))
        
        # Redraw the canvas
        self.canvas1.draw()
        print(f"Loaded data from {filename}")

     


    
    def on_motion(self, event):
        if self.dragging_marker and self.selected_primary_marker is not None:
            if event.inaxes == self.canvas1.axes:
                x, y = event.xdata, event.ydata
                if x is not None and y is not None:
                    # Update primary marker position while dragging
                    new_x_primary = x + self.offset[0]
                    new_y_primary = y + self.offset[1]
                    self.selected_primary_marker.set_offsets([(new_x_primary, new_y_primary)])

                    # Update conjugate marker position if it exists
                    if self.selected_conjugate_marker:
                        new_x_conjugate = new_x_primary
                        new_y_conjugate = -new_y_primary
                        self.selected_conjugate_marker.set_offsets([(new_x_conjugate, new_y_conjugate)])

                    # Update the zeros_list and poles_list outside the if condition
                    self.zeros_postions = [tuple(primary.get_offsets()[0]) if primary is not None else None for primary, _
                                    in self.zeros]
                    self.poles_postions = [tuple(primary.get_offsets()[0]) if primary is not None else None for primary, _
                                    in self.poles]

                    # Clear the axes before plotting
                    self.canvas_mag.axes.clear()
                    self.canvas_phase.axes.clear()

                    # Plot magnitude and phase responses
                    self.plot_frequency_response()
                    

                    # Redraw canvas
                    self.canvas1.draw()




     
    def on_release(self , event):
        if self.selected_primary_marker is not None:
            # Reset selected marker after release
            self.selected_primary_marker = None
            self.offset = (0, 0)
            self.dragging_marker = False
        self.plot_frequency_response()
        
    
    def Display_Filters(self):
        # Clear previous filter data
                # Clear existing plots
        self.canvas1.axes.clear()
        self.plot_unit_circle()
        
        self.poles, self.zeros = [], []
        self.poles_postions, self.zeros_postions = [], []
        exist=set()
        b,a=self.get_filter_data()
       
        # Convert filter coefficients to zeros and poles
        zeros, poles, _ = signal.tf2zpk(b, a)
        # for zero in zeros:
        #     print(zero)
        # print("//////////////////////////////////////////////////////////////")
        # for pole in poles :
        #     print(pole)
        # Plot the zeros (in blue)
        for zero in zeros:
            
            if zero.imag == 0:
                zero = zero.real

          
            conj_z = zero.conjugate() if isinstance(zero, complex) else None
            # If conjugate exists in the list, pair them; otherwise, mark "None"
            if conj_z in zeros and (zero, conj_z) not in exist and (conj_z, zero) not in exist:
                primary_marker = self.canvas1.axes.scatter(np.real(zero), np.imag(zero), color='b', marker='o', s=100,
                                                    facecolors='none', edgecolors='b', linewidths=2) 
                conj=self.canvas1.axes.scatter(np.real(zero), -np.imag(zero), color='b', marker='o', s=100,
                                                    facecolors='none', edgecolors='b', linewidths=2)
                
                self.zeros.append((primary_marker, conj))  # Append to zeros
                self.zeros_postions.append((np.real(zero), np.imag(zero)))
                exist.add((zero, conj_z))
            
            
            elif conj_z not in zeros:
                primary_marker = self.canvas1.axes.scatter(np.real(zero), np.imag(zero), color='b', marker='o', s=100,
                                                    facecolors='none', edgecolors='b', linewidths=2) 

                self.zeros.append((primary_marker, None))  # Append to zeros
                self.zeros_postions.append((np.real(zero), np.imag(zero)))
                exist.add((zero, conj_z))  # Mark this zero as processed
       
        for pole in poles:
            
            if pole.imag == 0:
                pole = pole.real

            # العثور على الكونجوجيت (conjugate)
            conj_p = pole.conjugate() if isinstance(pole, complex) else None
            # If conjugate exists in the list, pair them; otherwise, mark "None"
            if conj_p in poles and (pole, conj_p) not in exist and (conj_p, pole) not in exist:
                primary_marker = self.canvas1.axes.scatter(np.real(pole), np.imag(pole), color='r', marker='x', s=100, linewidths=2) 
                conj=self.canvas1.axes.scatter(np.real(pole), -np.imag(pole), color='r', marker='x', s=100, linewidths=2)
                
                self.poles.append((primary_marker, conj))  # Append to zeros
                self.poles_postions.append((np.real(pole), np.imag(pole)))
                exist.add((pole, conj_p))
            
            
            elif conj_p not in poles:
                primary_marker = self.canvas1.axes.scatter(np.real(pole), np.imag(pole), color='r', marker='x', s=100, linewidths=2)  

                self.poles.append((primary_marker, None))  # Append to zeros
                self.poles_postions.append((np.real(pole), np.imag(pole)))
                exist.add((pole, conj_p))  # Mark this zero as processed
            



        self.canvas1.draw()

        # Plot the frequency response (make sure this function exists and works as expected)
        self.plot_frequency_response()

    
    def get_filter_data(self):
        # Get the selected filter type from the combo box
        curr_filter = self.Combox_Filters.currentText()
        curr_type=self.Combox_Filters_type.currentText()
            # Define the filter coefficients based on the selected filter type
        if curr_filter == 'Butterworth':
            if curr_type == "LPF":
                # Butterworth Low Pass Filter
                b, a = signal.butter(4, 0.2, btype='low', analog=False)
            elif curr_type == "HPF":
                # Butterworth High Pass Filter
                b, a = signal.butter(4, 0.2, btype='high', analog=False)
 
        elif curr_filter == 'Chebyshev_I':
            if curr_type == "LPF":
                # Chebyshev Type I Low Pass Filter
                b, a = signal.cheby1(4, 1, 0.2, btype='low', analog=False)
            elif curr_type == "HPF":
                # Chebyshev Type I High Pass Filter
                b, a = signal.cheby1(4, 1, 0.2, btype='high', analog=False)

        elif curr_filter == 'Chebyshev_II':
            if curr_type == "LPF":
                # Chebyshev Type II Low Pass Filter (Inverse Chebyshev)
                b, a = signal.cheby2(4, 40, 0.2, btype='low', analog=False)
            elif curr_type == "HPF":
                # Chebyshev Type II High Pass Filter (Inverse Chebyshev)
                b, a = signal.cheby2(4, 40, 0.2, btype='high', analog=False)

        elif curr_filter == 'Bessel':
            if curr_type == "LPF":
                # Bessel Low Pass Filter
                b, a = signal.bessel(4, 0.2, btype='low', analog=False)
            elif curr_type == "HPF":
                # Bessel High Pass Filter
                b, a = signal.bessel(4, 0.2, btype='high', analog=False)

        elif curr_filter == 'Elliptic':
            if curr_type == "LPF":
                b, a = signal.ellip(4, 0.8, 40, 0.2, btype='low', analog=False)
            elif curr_type == "HPF":
                # Elliptic High Pass Filter
                b, a = signal.ellip(4, 1, 40, 0.2, btype='high', analog=False)
        return b,a




    def ensure_conjugates(self):
        """
        Ensure that every pole and zero has a conjugate. 
        If the conjugate marker is not None, add the conjugate position.
        If the conjugate marker is None, do nothing.
        """
        # Ensure conjugates for poles
        updated_poles = []
        poles_positions = self.poles_postions
        zeros_positions = self.zeros_postions

        for (x, y), marker in zip(poles_positions, self.poles):
            if marker[1] is not None:  # If conjugate marker is not None
                conjugate = (x, -y)  # The conjugate position
                if conjugate not in poles_positions:  # Add conjugate position if it doesn't exist
                    updated_poles.append(conjugate)
            # If marker[1] is None, do nothing

        # Add new conjugates to poles
        poles_positions.extend(updated_poles)

        # Ensure conjugates for zeros
        updated_zeros = []
        for (x, y), marker in zip(zeros_positions, self.zeros):
            if marker[1] is not None:  # If conjugate marker is not None
                conjugate = (x, -y)  # The conjugate position
                if conjugate not in zeros_positions:  # Add conjugate position if it doesn't exist
                    updated_zeros.append(conjugate)
            # If marker[1] is None, do nothing

        # Add new conjugates to zeros
        zeros_positions.extend(updated_zeros)

        return poles_positions, zeros_positions


    def direct_form_II(self, poles, zeros):
        """
        Perform Direct Form II realization of the filter.
        """
        numerator = [1]
        denominator = [1]
        for x, y in poles:
            denominator = np.convolve(denominator, [1, -complex(x, y)])
        for x, y in zeros:
            numerator = np.convolve(numerator, [1, -complex(x, y)])
        return numerator, denominator

    def cascade_form(self, poles, zeros):
        """
        Perform Cascade realization (splitting into second-order sections).
        """
        sections = []
        if len(poles) % 2 != 0 or len(zeros) % 2 != 0:
            raise ValueError("Number of poles and zeros must be even for cascade form.")

        for i in range(0, len(poles), 2):
            section_poles = poles[i:i+2]
            section_zeros = zeros[i:i+2]

            numerator = [1]
            denominator = [1]

            for x, y in section_poles:
                denominator = np.convolve(denominator, [1, -complex(x, y)])
            for x, y in section_zeros:
                numerator = np.convolve(numerator, [1, -complex(x, y)])

            sections.append((numerator, denominator))
        return sections

    def export_to_C_code(self, numerator, denominator, method, section_index=None):
        """
        Export the filter realization as C code.
        """
        if method == 'direct_form_II':
            code = f"""
#include <complex.h>  // For complex numbers

double complex b[] = {{{', '.join(f'{coef.real} + {coef.imag}*I' for coef in numerator)}}};
double complex a[] = {{{', '.join(f'{coef.real} + {coef.imag}*I' for coef in denominator[1:])}}};

double complex x[{len(numerator)}] = {{0}};  // Input buffer
double complex y[{len(denominator)}] = {{0}};  // Output buffer

double complex filter(double complex input) {{
    for (int i = {len(numerator)-1}; i > 0; i--) {{
        x[i] = x[i-1];  // Shift input
        y[i] = y[i-1];  // Shift output
    }}

    x[0] = input;
    y[0] = 0 + 0*I;

    for (int i = 0; i < {len(numerator)}; i++) {{
        y[0] += b[i] * x[i];
    }}
    for (int i = 1; i < {len(denominator)}; i++) {{
        y[0] -= a[i] * y[i];
    }}

    return y[0];
}}

int main() {{
    double complex input_signal[] = {{1 + 0*I, 2 + 1*I, 3 - 1*I, 0 + 0*I}};
    int n = sizeof(input_signal) / sizeof(input_signal[0]);

    printf("Filtered output:\\n");
    for (int i = 0; i < n; i++) {{
        double complex output = filter(input_signal[i]);
        printf("Input: %.2f%+.2fi, Output: %.2f%+.2fi\\n",
               creal(input_signal[i]), cimag(input_signal[i]),
               creal(output), cimag(output));
    }}

    return 0;
}}
"""
        elif method == 'cascade_form':
            code = f"""
#include <complex.h>  // For complex numbers

// Section {section_index + 1}
double complex b{section_index + 1}[] = {{{', '.join(f'{coef.real} + {coef.imag}*I' for coef in numerator)}}};
double complex a{section_index + 1}[] = {{{', '.join(f'{coef.real} + {coef.imag}*I' for coef in denominator[1:])}}};

double complex x{section_index + 1}[{len(numerator)}] = {{0}};  // Input buffer
double complex y{section_index + 1}[{len(denominator)}] = {{0}};  // Output buffer

double complex filter{section_index + 1}(double complex input) {{
    for (int i = {len(numerator)-1}; i > 0; i--) {{
        x{section_index + 1}[i] = x{section_index + 1}[i-1];  // Shift input
        y{section_index + 1}[i] = y{section_index + 1}[i-1];  // Shift output
    }}

    x{section_index + 1}[0] = input;
    y{section_index + 1}[0] = 0 + 0*I;

    for (int i = 0; i < {len(numerator)}; i++) {{
        y{section_index + 1}[0] += b{section_index + 1}[i] * x{section_index + 1}[i];
    }}
    for (int i = 1; i < {len(denominator)}; i++) {{
        y{section_index + 1}[0] -= a{section_index + 1}[i] * y{section_index + 1}[i];
    }}

    return y{section_index + 1}[0];
}}

int main() {{
    double complex input_signal[] = {{1 + 0*I, 2 + 1*I, 3 - 1*I, 0 + 0*I}};
    int n = sizeof(input_signal) / sizeof(input_signal[0]);

    printf("Filtered output:\\n");
    for (int i = 0; i < n; i++) {{
        double complex output = filter{section_index + 1}(input_signal[i]);
        printf("Input: %.2f%+.2fi, Output: %.2f%+.2fi\\n",
               creal(input_signal[i]), cimag(input_signal[i]),
               creal(output), cimag(output));
    }}

    return 0;
}}
"""
        return code

    def realize_and_export(self):
        """
        Realize the filter and export the C code.
        """
        poles_postions , zeros_postions =self.ensure_conjugates()
        self.filter_type=self.Combox_Generat_C_code.currentText()
        # Collect the C code for the entire filter
        full_c_code = f"Generated C Code for Filter Realization\nFilter Type: {self.filter_type}\n\n"
       
        if self.filter_type == 'direct_form_II':
            numerator, denominator = self.direct_form_II(poles_postions, zeros_postions)
            full_c_code += self.export_to_C_code(numerator, denominator, 'direct_form_II')
        
                    # Now save the full C code into a single PDF
            self.save_pdf(full_c_code, "direct_form_II_filter_C_code.pdf", self.filter_type)
        elif self.filter_type == 'cascade_form':
            sections = self.cascade_form(poles_postions, zeros_postions)
            # Generate the code for each section and append to the full C code string
            for section_index, (numerator, denominator) in enumerate(sections):
                section_code = self.export_to_C_code(numerator, denominator, 'cascade_form', section_index)
                full_c_code += f"\n\n{section_code}"

            # Now save the full C code into a single PDF
            self.save_pdf(full_c_code, "cascade_filter_C_code.pdf", self.filter_type)
            self.Combox_Generat_C_code.setCurrentText("Generate_C_Code")

    def save_pdf(self, c_code, filename, filter_type):
        """
        Save the generated C code to a PDF file.
        """
        # Create a PDF canvas
        c = canvas.Canvas(filename, pagesize=letter)
        width, height = letter  # Get the dimensions of the page

        # Set up some styling
        c.setFont("Times-Roman", 12)  # Use Times-Roman font, size 12 (larger and clearer)
        margin = 40  # Margin around the text
        line_height = 16  # Line height for the text
        text_object = c.beginText(margin, height - margin - 40)  # Start text at the top-left corner
        text_object.setFont("Times-Roman", 12)

        # Add a title at the top of the page
        c.setFont("Times-Bold", 16)  # Larger bold font for the title
        c.drawString(margin, height - margin - 20, f"Filter Design C Code - {filter_type} Form")

        # Set text object font back to regular
        c.setFont("Times-Roman", 12)

        # Split the C code into lines to fit the page width
        lines = simpleSplit(c_code, "Times-Roman", 12, width - 2 * margin)

        # Add each line of code to the PDF
        for line in lines:
            # If the line doesn't fit in the current page, create a new page
            if text_object.getY() < margin + line_height:
                c.drawText(text_object)
                c.showPage()  # Create a new page
                text_object = c.beginText(margin, height - margin - 40)
                text_object.setFont("Times-Roman", 12)
            text_object.textLine(line)

        # Draw the remaining text on the canvas
        c.drawText(text_object)

        # Save the PDF file
        c.showPage()  # Ensure all content is drawn on the last page
        c.save()

        print(f"PDF saved as {filter_type}")


    def Open_Second_Window(self):
        main = MyDialog(QtWidgets.QDialog)
        main.exec_()










# ////////////////////////////////////////////////////////////////////////////////////////////////





            
    def calculateMousePosition(self, event):
        cursor_pos = event.pos()
        # print(f"Mouse Position: ({cursor_pos.x()}, {cursor_pos.y()})")
        self.magnitude = (cursor_pos.x() + cursor_pos.y() ) /2
        self.drawGraph()


    def drawGraph(self):
        self.mag_list.append (self.magnitude)
        self.time_list.append(self.time_list[-1] + 0.0005)
        # print("Contents of self.mag_list:", self.mag_list)

        if self.time_list[-1] > 1.5:
            self.Real_Signal.plotItem.setXRange(
                self.time_list[-1] - 1.5,
                self.time_list[-1] + 1.5,
                padding=0
            )
            self.hoveredOutput.plotItem.setXRange(
                self.time_list[-1] - 1.5,
                self.time_list[-1] + 1.5,
                padding=0
            )

        self.Real_Signal.plotItem.clear()
        self.Real_Signal.plotItem.plot(
            x=self.time_list[-500:],
            y=self.mag_list[-500:],
           pen={'color': 'b', 'width': 2}# Use 'b' for blue color
            

        )

        zeros_array = np.array(self.zeros_postions).flatten()
        poles_array = np.array(self.poles_postions).flatten()
        
        numerator, denominator = signal.zpk2tf(zeros_array, poles_array, 1)
        self.filtered_signal = signal.lfilter(numerator, denominator, self.mag_list.copy())

        # print(self.filtered_signal)

        self.Filtered_Signal.plotItem.clear()
        self.Filtered_Signal.plotItem.plot(
            x=self.time_list[-500:],
            y= np.real(self.filtered_signal[-500:]),
            pen={'color': 'b', 'width': 2}# Use 'b' for blue color  
        )

    def load_signal(self):
        self.clear_graphs()
        
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Signal", "", "CSV Files (*.csv);;All Files (*)")
        if not file_path:
            return  # User canceled the file dialog

        try:
            # Load the signal data from the CSV file
            data = pd.read_csv(file_path)
            if 'Time' not in data.columns or 'Signal' not in data.columns:
                raise ValueError("CSV file must contain 'Time' and 'Signal' columns.")

            time = data['Time'].to_numpy()
            signal_data = data['Signal'].to_numpy()

            # Plot the original signal
            self.Real_Signal.plotItem.clear()
            self.Real_Signal.plotItem.plot(x=time, y=signal_data, pen='b')

            # Apply the filter to the loaded signal
            zeros_array = np.array(self.zeros_postions).flatten()
            poles_array = np.array(self.poles_postions).flatten()
            numerator, denominator = signal.zpk2tf(zeros_array, poles_array, 1)
            filtered_signal = signal.lfilter(numerator, denominator, signal_data)

            # Plot the filtered signal
            self.Filtered_Signal.plotItem.clear()
            self.Filtered_Signal.plotItem.plot(x=time, y=np.real(filtered_signal), pen='g')

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load signal: {str(e)}")



    
    def clear_graphs(self):
        self.mag_list = [0]
        self.time_list = [0]
        self.Real_Signal.plotItem.clear()
        self.Filtered_Signal.plotItem.clear()          



def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()