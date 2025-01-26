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
import sys
from PyQt5.QtWidgets import QDialog
from scipy import signal
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyqtgraph as pg
matplotlib.use('Qt5Agg')



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
        uic.loadUi(r'RealTime_Signal_Filter\MainWindow - untitled.ui', self)
        self.graph = pg.PlotItem()

        

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

        

        self.button_ClearZeros.clicked.connect(lambda : self.delete_zeros())
        self.button_ClearPoles.clicked.connect(lambda : self.delete_poles())
        self.Button_Clear.clicked.connect(lambda : self.delete_all())

        self.Button_Swap.clicked.connect(self.swap_markers_and_positions)

                # Plot unit circle
        self.plot_unit_circle()
        
        self.zeros = []     
        self.poles = []

        self.zeros_postions = []
        self.poles_postions = []

        self.canvas1.mpl_connect('button_press_event', self.on_press)
        self.plot_frequency_response_mag()


# plot the unit circle widget
    def plot_unit_circle(self):
    # Plot unit circle using Matplotlib on the canvas
        self.draw_circle()

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
        self.canvas1.draw()
        self.circle_graph.setCentralItem(self.graph)
        self.layout1.setContentsMargins(0, 0, 0, 0)
        self.circle_graph.setLayout(self.layout1)
    
    def draw_circle(self):
            # Plot unit circle using Matplotlib on the canvas
        angles = np.linspace(0, 2 * np.pi, 1000)
        x = np.cos(angles)
        y = np.sin(angles)
        self.canvas1.axes.plot(x, y, 'k')
        self.canvas1.axes.axhline(0, color='k')  
        
        # Plot the y-axis (vertical line)
        self.canvas1.axes.axvline(0, color='k')  
        
        
    def plot_frequency_response_mag(self):
            # Create a system from zeros and poles
            zeros = [complex(x, y) for x, y in self.zeros_postions]
            poles = [complex(x, y) for x, y in self.poles_postions]
            self.canvas_mag.axes.set_facecolor('white')
            w, h = signal.freqz_zpk(zeros, poles, 1.0)
           
            # Plot the magnitude response
            self.canvas_mag.axes.clear()
            self.canvas_mag.axes.semilogx(w, abs(h))
            self.canvas_mag.axes.set_xlabel('Frequency')
            self.canvas_mag.axes.set_ylabel('Magnitude (dB)')
            self.canvas_mag.axes.set_title('Magnitude Response')
            self.canvas_mag.axes.grid()
            self.canvas_mag.axes.set_xlim(1e-2, 10)
            self.canvas_mag.draw()
            self.Mag_Graph.setCentralItem(self.graph)
            self.layout_mag.setContentsMargins(0, 0, 0, 0)
            self.Mag_Graph.setLayout(self.layout_mag)


            # Plot the phase response
            self.canvas_phase.axes.clear()
            self.canvas_phase.axes.semilogx(w, np.angle(h, deg=True))
            self.canvas_phase.axes.set_xlabel('Frequency')
            self.canvas_phase.axes.set_ylabel('Phase (degrees)')
            self.canvas_phase.axes.set_title('Phase Response')
            self.canvas_phase.axes.grid()
            self.canvas_phase.axes.set_xlim(1e-2, 10)
            self.canvas_phase.draw()
            self.Phase_Graph.setCentralItem(self.graph)
            self.layout_phase.setContentsMargins(0, 0, 0, 0)
            self.Phase_Graph.setLayout(self.layout_phase)
    


    def plot_frequency_response_phase(self, canvas_phase, widget, layout, flag):
        # Create a system from zeros and poles
        zeros = [complex(x, y) for x, y in self.zeros_list]
        poles = [complex(x, y) for x, y in self.poles_list]
        print("zeros complex : ", zeros)
        print("poles complex : ", poles)

                # Compute the frequency response
        w, h = signal.freqz_zpk(zeros, poles, 1.0)






    # to handle all actions on the unit circle >> move , delete specific object , draw new object
    def on_press(self , event ):
        if event.button == 1:  # Left mouse button
            if event.inaxes == self.canvas1.axes:
                x, y = event.xdata, event.ydata
                if x is not None and y is not None:
                    self.draw_New_object(x,y)

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
        self.plot_frequency_response_mag()


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
        self.plot_frequency_response_mag()
        # self.plot_frequency_response_phase(self.canvas_phase, self.phase_responce_graph, self.layout_phase, flag=False)
    
    



# update the graph after delete the all poles ,  all zeros , both
    def redraw_plot(self, delete_zeros=False, delete_poles=False):
        # Redraws the plot after removing zeros or poles
        self.canvas1.axes.clear()
        self.draw_circle()
     
        
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
        self.plot_frequency_response_mag()

    def delete_zeros(self):
        # Delete all zeros on the plot      
        self.zeros = []
        self.zeros_postions = []
        self.redraw_plot(delete_zeros=True)
        self.plot_frequency_response_mag()
        # self.plot_frequency_response_phase(self.canvas_phase, self.phase_responce_graph, self.layout_phase, flag=False)

    def delete_poles(self):
        # Delete all poles on the plot  
        self.poles = []
        self.poles_postions = []
        self.redraw_plot(delete_poles=True)
        self.plot_frequency_response_mag()
        # self.plot_frequency_response_phase(self.canvas_phase, self.phase_responce_graph, self.layout_phase, flag=False)

    def delete_all(self):
        # Delete all poles and zeros on the plot
        self.zeros = []
        self.zeros_postions = []
        self.poles = []
        self.poles_postions = []
        self.redraw_plot()
        self.plot_frequency_response_mag()
        # self.plot_frequency_response_phase(self.canvas_phase, self.phase_responce_graph, self.layout_phase, flag=False)

    def swap_markers_and_positions(self):
        # Clear existing 
        for obj in self.poles + self.zeros:
            obj[0].remove()
            if obj[1] is not None :obj[1].remove()

        New_Poles, New_Zeros = [], []  # Lists to store new poles and zeros
        New_poles_pos, New_Zeros_Pos = [], []  # Lists to store new pole and zero positions

        # Swap poles to zeros
        for indx, (_, conj) in enumerate(self.poles):            
            # Use the corresponding position, or last position if index is out of range
            x, y = self.poles_postions[indx] if indx < len(self.poles_postions) else self.poles_postions[-1]
            
            if conj is not None:
                conjugate_marker = self.canvas1.axes.scatter(x, -y, color='b', marker='o', s=100, facecolors='none', edgecolors='b', linewidths=2)
            else:
                conjugate_marker = None
            primary_marker = self.canvas1.axes.scatter(x, y, color='b', marker='o', s=100, facecolors='none', edgecolors='b', linewidths=2)
            
            New_Zeros.append((primary_marker, conjugate_marker))
            New_Zeros_Pos.append((x, y))
        
        # Swap zeros to poles
        for indx, (_, conj) in enumerate(self.zeros):
            # Use the corresponding position, or last position if index is out of range
            x, y = self.zeros_postions[indx] if indx < len(self.zeros_postions) else self.zeros_postions[-1]
            
            if conj is not None:
                conjugate_marker = self.canvas1.axes.scatter(x, -y, color='r', marker='x', s=100, linewidths=2)
            else:
                conjugate_marker = None
            primary_marker = self.canvas1.axes.scatter(x, y, color='r', marker='x', s=100, linewidths=2)
            
            New_Poles.append((primary_marker, conjugate_marker))
            New_poles_pos.append((x, y))

        # Update the poles and zeros with the new swapped lists
        self.poles, self.zeros = New_Poles, New_Zeros
        self.poles_postions, self.zeros_postions = New_poles_pos, New_Zeros_Pos
        
        self.canvas1.draw()  # Redraw the canvas with updated markers
        self.canvas1.flush_events()
        self.plot_frequency_response_mag()
            
            



def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()