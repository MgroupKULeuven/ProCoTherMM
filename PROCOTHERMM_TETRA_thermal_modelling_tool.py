# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 08:32:48 2019

@author: frederik


Contents of GUI
1D transient problem
2D transient problem

Ability to change material parameters
Ability to change boundary conditions
Visualisation of data 
Extraction of data

"""

from PyQt5 import QtGui, uic 
from PyQt5.QtWidgets import QMainWindow, QApplication,QDesktopWidget, QWidget, QSizePolicy,QPushButton,QLabel,QVBoxLayout, QSlider, QLineEdit, QCheckBox 
from PyQt5.QtCore import QThread, pyqtSignal, Qt
import sys
import numpy as np
import matplotlib.cm as pltcm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import time
import csv

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5.5, height=4.5, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,QSizePolicy.Expanding,QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.xlabel = ''
        self.ylabel = ''
        self.title = ''
        self.plot(0,0)        
    def set_xlabel(self,tex):
        self.xlabel = tex
    def set_ylabel(self,tex):
        self.ylabel = tex
    def set_title(self,tex):
        self.title = tex
    def plot(self,x,y):
        self.axes.clear()
        self.axes.plot(x,y)
        self.axes.set_xlabel(self.xlabel)
        self.axes.set_ylabel(self.ylabel)
        self.axes.set_title(self.title)
        self.draw()
    def plot_extra(self,x,y):
        self.axes.plot(x,y)
        self.axes.set_xlabel(self.xlabel)
        self.axes.set_ylabel(self.ylabel)
        self.axes.set_title(self.title)
        self.draw() 
    def contourf(self,val):
        self.axes.clear()
        self.axes.contourf(val,20)
        self.axes.set_xlabel(self.xlabel)
        self.axes.set_ylabel(self.ylabel)
        self.axes.set_title(self.title)
        self.draw()
    def plot_point_extra(self,x,y):
        self.axes.scatter(x,y)
        self.draw()

##############################################################################
def deg2K(T):
    return T + 273

##############################################################################
def K2deg(T):
    return T - 273

############################################################################## 
def lookuptableinterpol(default,x,y,X):
    result = default
    for i in range(len(x)-1):
        if (X>=x[i])and(X<x[i+1]):
            result = linpol(x[i],x[i+1],y[i],y[i+1],X)
    return result

##############################################################################    
def linpol(x0,x1,k0,k1,x):
    return ((k1-k0)/(x1-x0))*x + (k0 - ((k1-k0)/(x1-x0))*x0)
 

##############################################################################
class Materiaal_Steel():
    def __init__(self):
        None
    def rho(self,T):
        return 7800
    def k(self,T):
        T = K2deg(T)
        temps = [400,500,600,700,800,900,2000] #°C
        ks = [50,40,36,31.8,27.7,27.2,27.2]
        default_k = 50
        return lookuptableinterpol(default_k,temps,ks,T)  
    def cp(self,T):
        T = K2deg(T)
        cp = 500
        temps = [0,300,600,700,800,1000,1600] #°C
        cps = [450,550,750,950,800,650,630]
        return lookuptableinterpol(cp,temps,cps,T)
    def epsilon(self,T):
        T = K2deg(T)
        return 0.775 + 0.0001*T
    def sigma(self,T): 
        return 5.669*(10**(-8)) # W/(m^2 K^4)
    def contactconductance(self,T):
        return 7

    
class simulation_part():
    def __init__(self,dimensions,ID,T_init_degreesCelsius, simulation_time):
        self.sim_data = dict()
        self.boundary_conditions = dict()
        self.sim_data['T_init_degreesCelsius'] = T_init_degreesCelsius
        self.sim_data['LX'] = 0.22
        self.sim_data['A_x'] = 10
        self.sim_data['Delta_X'] = 0.02
        self.sim_data['NX'] = int(self.sim_data['LX']/self.sim_data['Delta_X'])
        self.sim_data['NY'] = 1 # set second dimension off
        self.sim_data['A_y'] = self.sim_data['Delta_X']
        self.dimensions = dimensions
       
        if self.dimensions==2:
            self.sim_data['LY'] = 0.22
            self.sim_data['Delta_Y'] = 0.02
            self.sim_data['NY'] = int(self.sim_data['LY']/self.sim_data['Delta_Y'])
            self.sim_data['A_x'] = 10/self.sim_data['NX'] # particle surface size on bounds in x
            self.sim_data['A_y'] = 10/self.sim_data['NY'] # particle surface size on bounds in y
            self.y = np.linspace(0,self.sim_data['LY'],self.sim_data['NY'])
            
        self.x = np.linspace(0,self.sim_data['LX'],self.sim_data['NX'])
        
        self.TimeIndex = 0
        self.material = Materiaal_Steel()
        self.init_simulation_time(simulation_time)
        self.init_temperature()
        
    def init_temperature(self,temp='None'):
        if not(temp=='None'):
            self.sim_data['T_init_degreesCelsius'] = temp
        print('new initial temeprature = ', self.sim_data['T_init_degreesCelsius'])
        self.Temperature_init = deg2K(self.sim_data['T_init_degreesCelsius'])*np.ones((self.sim_data['NX'],self.sim_data['NY']))
        self.Temperature = np.zeros((self.sim_data['Nt'],self.sim_data['NX'],self.sim_data['NY'])) # TODO: create dynamic size!
        self.Temperature[0,:,:] = self.Temperature_init
        self.average_Temperature = np.zeros((self.sim_data['Nt']))
        self.average_Temperature[0] = self.get_average_temperature(0)

    def init_simulation_time(self,time):
        self.sim_data['simulation_time'] = time 
        self.sim_data['Nt'] =  1000
        self.sim_data['Delta_T'] = float(self.sim_data['simulation_time']/self.sim_data['Nt']) # delta_time
        self.time = np.linspace(0,self.sim_data['simulation_time'],self.sim_data['Nt'])
        self.init_temperature()
        
    def set_boundary_condition(self,entry,value):
        self.boundary_conditions[entry] = value

    def simulate_time_step(self):
        j = self.TimeIndex
        # GET BOUNDARY CONDITIONS FROM NEIGBOURS            
        #print('ID_  ', self.ID)
        integration = Euler_(self.Temperature[j,:,:],self.sim_data,self.material,self.boundary_conditions)
       # print(integration)
        if self.dimensions ==1:
            self.Temperature[j+1,:] = np.vstack(integration)
        if self.dimensions ==2:
            self.Temperature[j+1,:,:] = integration
        # compute average_temperature
        self.average_Temperature[j+1] = self.get_average_temperature(j+1)
        #print(self.Temperature)
        self.TimeIndex = self.TimeIndex+1
        

    def get_average_temperature(self,time_index):
        return np.mean(self.Temperature[time_index,:,:])
        

##############################################################################
def ODE_(T,k_list,sim_data,material,boundary_conditions):
    dim = len(k_list)
    k_x = k_list[0]
    k_y = 0
    if dim==2:
        k_y = k_list[1]
        
    nx = sim_data['NX']
    ny = sim_data['NY']

    m = sim_data['A_x']*sim_data['A_y']*material.rho(T[k_x,k_y])
    ## x-direction  = VERTICAL
    if k_x < nx-1: 
        conduction_top = -material.k(T[k_x,k_y])*sim_data['A_x']*( T[k_x+1,k_y]- T[k_x,k_y] )/sim_data['Delta_X']
    if k_x==nx-1: # top
        conduction_top = 0
        convection_top = boundary_conditions['h_top']*sim_data['A_x']*(T[k_x,k_y] - boundary_conditions['T_env_top'])
        topradiation = 0
        topevaporation = 0
        # TODO: radiation and evaporation
    if k_x>0:
        conduction_bottom = -material.k(T[k_x,k_y])*sim_data['A_x']*( T[k_x,k_y]- T[k_x-1,k_y] )/sim_data['Delta_X']
    if k_x==0: #Bottom
        convection_bottom = boundary_conditions['h_bottom']*sim_data['A_x']*(T[k_x,k_y] - boundary_conditions['T_env_bottom'])
        conduction_bottom = 0  

     
    if dim ==2:    ## y-direction = HORIZONTAL
        if k_y < ny-1: 
            conduction_right = -material.k(T[k_x,k_y])*sim_data['A_y']*( T[k_x,k_y+1]- T[k_x,k_y] )/sim_data['Delta_Y']
        if k_y==ny-1: # right
            conduction_right = 0
            convection_right = boundary_conditions['h_right']*sim_data['A_y']*(T[k_x,k_y] - boundary_conditions['T_env_right'])
        if k_y>0:
            conduction_left = -material.k(T[k_x,k_y])*sim_data['A_y']*( T[k_x,k_y]- T[k_x,k_y-1] )/sim_data['Delta_Y']
        if k_y==0: #left
            convection_left = boundary_conditions['h_left']*sim_data['A_y']*(T[k_x,k_y] - boundary_conditions['T_env_left'])
            conduction_left = 0  
    

    # compute right part of ODE
    dT_dt = 1/( m*material.cp(T[k_x,k_y]) ) * (conduction_bottom - conduction_top  ) 
    if k_x == nx-1: # top
        dT_dt = dT_dt + 1/( m*material.cp(T[k_x,k_y]))*(- convection_top - topradiation - topevaporation)
    if k_x==0: # bottom
        dT_dt = dT_dt + 1/( m*material.cp(T[k_x,k_y]))*(- convection_bottom)
    if dim==2:
        dT_dt = dT_dt + 1/( m*material.cp(T[k_x,k_y]) ) * (conduction_left - conduction_right)
        if k_y == ny-1: # right
            dT_dt = dT_dt + 1/( m*material.cp(T[k_x,k_y]))*(- convection_right )
        if k_y==0: # left
            dT_dt = dT_dt + 1/( m*material.cp(T[k_x,k_y]))*(- convection_left )
    
    return dT_dt

##############################################################################
def Euler_(T,sim_data,material,boundary_conditions):
    '''
    simulate thermal ODE one time step using Euler integration
    '''
    nx = sim_data['NX']
    ny = sim_data['NY']
    if not(sim_data['NY'] == 1):
        T_ = np.zeros((nx,ny))
        for k_x in range(nx):
            for k_y in range(ny):
                T_[k_x,k_y] = T[k_x,k_y] + sim_data['Delta_T']*ODE_(T,[k_x,k_y],sim_data,material,boundary_conditions)
    else:
        T_ = np.zeros((nx))
        for k in range(nx):
            T_[k] = T[k] + sim_data['Delta_T']*ODE_(T,[k],sim_data,material,boundary_conditions)     
    return T_

##############################################################################
def RK4(T,sim_data,material,boundary_conditions):
    '''
    simulate thermal ODE one time step using Runge Kutta 4 integration
    '''
    nx = sim_data['NX']
    T_hat_j_plus_half = np.zeros((nx))
    T_hat_hat_j_plus_half = np.zeros((nx))
    T_hat_j_plus_1 = np.zeros((nx))
    Kutta1 = np.zeros((nx))
    Kutta2 = np.zeros((nx))
    Kutta3 = np.zeros((nx))
    Kutta4 = np.zeros((nx))
    T_ = np.zeros((nx))
    for k in range(nx):
        Kutta1[k] = ODE(T,k,sim_data,material,boundary_conditions)
        T_hat_j_plus_half[k] = T[k] + (sim_data['Delta_T']/2.)*Kutta1[k]
    for k in range(nx):
        Kutta2[k] = ODE(T_hat_j_plus_half,k,sim_data,material,boundary_conditions)
        T_hat_hat_j_plus_half[k] = T[k] + (sim_data['Delta_T']/2.)*Kutta2[k]
    for k in range(nx):
        Kutta3[k] = ODE(T_hat_hat_j_plus_half,k,sim_data,material,boundary_conditions)
        T_hat_j_plus_1[k] = T[k] + (sim_data['Delta_T'])*Kutta3[k]
    for k in range(nx):
        Kutta4[k] = ODE(T_hat_j_plus_1,k,sim_data,material,boundary_conditions)
        T_[k] = T[k] + (sim_data['Delta_T']/6.)*( Kutta1[k] + 2.*Kutta2[k] + 2.*Kutta3[k] + Kutta4[k] )        
    return T_


###############################################################################
class Gui(QMainWindow):
    def __init__(self, parent = None):
        super(Gui, self).__init__(parent)
        self.setWindowTitle('Gui for Thermal Simulation - [Frederik Debrouwere]')
        self.centralWidget = QWidget(self)
        sizeObject = QDesktopWidget().screenGeometry(-1)
        self.GEOMETRICS = {}
        self.GEOMETRICS['app_width'] = 1500 # sizeObject.width() #RectScreen.width()#1920 #1500
        self.GEOMETRICS['app_height'] = 900#sizeObject.height()-100
        self.resize(self.GEOMETRICS['app_width'],self.GEOMETRICS['app_height'])
        self.centralWidget.setGeometry(0, 0, self.GEOMETRICS['app_width'], self.GEOMETRICS['app_height'])
        self.solver = dict()

        self.dimension = 2

        self.part = simulation_part(2,'ID1',100,1000) # create dummy part
        
        DeltaSliderX = 45
        DeltaSliderY = 55
        self.widgetContents = {}
        self.slab_widgets = {}
        self.widgetContents['QPushButton'] = {}
        self.widgetContents['QPushButton']['startbtn'] = QPushButton('simulate',self.centralWidget)
        self.widgetContents['QPushButton']['startbtn'].clicked.connect(self.simulate_and_plot)
        self.widgetContents['QPushButton']['startbtn'].setGeometry(700, 20, 200, 30)
        self.widgetContents['QPushButton']['export'] = QPushButton('export T to csv',self.centralWidget)
        self.widgetContents['QPushButton']['export'].clicked.connect(self.write_to_csv)
        self.widgetContents['QPushButton']['export'].setGeometry(700, 20+35, 200, 30)
        self.widgetContents['QLabel'] = {}
        self.widgetContents['QLabel']['log'] = QLabel('log',self.centralWidget)
        self.widgetContents['QLabel']['log'].setGeometry(0, 0, 400, 20)

        self.widgetContents['graphs'] = {}
        self.widgetContents['graphs']['controuf'] = PlotCanvas(self.centralWidget)
        self.widgetContents['graphs']['controuf'].move(100,300)
        self.widgetContents['graphs']['controuf'].set_xlabel('y position [.]')
        self.widgetContents['graphs']['controuf'].set_ylabel('x position [.]')
        self.widgetContents['graphs']['controuf'].set_title('Contour plot: z = temperature [°C]')
        self.widgetContents['graphs']['xy_t'] = PlotCanvas(self.centralWidget)
        self.widgetContents['graphs']['xy_t'].move(900,300)
        self.widgetContents['graphs']['xy_t'].set_xlabel('time [s]')
        self.widgetContents['graphs']['xy_t'].set_ylabel('temperature [°C]')
        self.widgetContents['QSlider'] = {}
        self.widgetContents['QSlider']['time'] = QSlider(Qt.Horizontal,self.centralWidget)
        self.widgetContents['QSlider']['time'].setGeometry(900+DeltaSliderY, 800, 550-2*DeltaSliderY, 30)
       # self.widgetContents['QSlider']['time'].valueChanged[int].connect(self.moved_slider_time)
        self.widgetContents['QSlider']['time'].sliderReleased.connect(self.moved_slider_time)
        self.widgetContents['QSlider']['time'].setMinimum(0)
        self.widgetContents['QSlider']['time'].setMaximum(self.part.sim_data['Nt']-1)
        self.widgetContents['QSlider']['time'].setValue(0)
        self.widgetContents['QSlider']['time'].setTickPosition(QSlider.TicksBelow)
        self.widgetContents['QSlider']['time'].setTickInterval(int(self.part.sim_data['Nt']/10))
        self.widgetContents['QLabel']['timeSlider'] = QLabel('Time [s]',self.centralWidget)
        self.widgetContents['QLabel']['timeSlider'].move(850,800)
        self.widgetContents['QLabel']['timeValue'] = QLabel('0',self.centralWidget)
        self.widgetContents['QLabel']['timeValue'].setGeometry(850,800+20, 100, 20)
        self.widgetContents['QSlider']['x'] = QSlider(Qt.Vertical,self.centralWidget)
        self.widgetContents['QSlider']['x'].setGeometry(20, 300+DeltaSliderX, 30, 450-2*DeltaSliderX)
        self.widgetContents['QSlider']['x'].sliderReleased.connect(self.moved_slider_x)
        self.widgetContents['QSlider']['x'].setMinimum(0)
        self.widgetContents['QSlider']['x'].setMaximum(self.part.sim_data['NX']-1)
        self.widgetContents['QSlider']['x'].setValue(0)
        self.widgetContents['QSlider']['x'].setTickPosition(QSlider.TicksBelow)
        self.widgetContents['QSlider']['x'].setTickInterval(int(self.part.sim_data['NX']/10))
        self.widgetContents['QSlider']['y'] = QSlider(Qt.Horizontal,self.centralWidget)
        self.widgetContents['QSlider']['y'].setGeometry(100+DeltaSliderY, 800, 550-2*DeltaSliderY, 30)
        self.widgetContents['QSlider']['y'].sliderReleased.connect(self.moved_slider_y)
        self.widgetContents['QSlider']['y'].setMinimum(0)
        self.widgetContents['QSlider']['y'].setMaximum(self.part.sim_data['NY']-1)
        self.widgetContents['QSlider']['y'].setValue(0)
        self.widgetContents['QSlider']['y'].setTickPosition(QSlider.TicksBelow)
        self.widgetContents['QSlider']['y'].setTickInterval(int(self.part.sim_data['NY']/10))
        self.widgetContents['QLabel']['xSlider'] = QLabel('X [m]',self.centralWidget)
        self.widgetContents['QLabel']['xSlider'].move(75,800-60)
        self.widgetContents['QLabel']['xValue'] = QLabel('0',self.centralWidget)
        self.widgetContents['QLabel']['xValue'].setGeometry(75,800+20-60, 100, 20)
        self.widgetContents['QLabel']['ySlider'] = QLabel('Y [m]',self.centralWidget)
        self.widgetContents['QLabel']['ySlider'].move(75,800)
        self.widgetContents['QLabel']['yValue'] = QLabel('0',self.centralWidget)
        self.widgetContents['QLabel']['yValue'].setGeometry(75,800+20, 100, 20)
        self.widgetContents['QLineEdit'] = {}
        self.widgetContents['QLineEdit']['SimulationTime'] = QLineEdit('',self.centralWidget)
        self.widgetContents['QLineEdit']['SimulationTime'].setGeometry(1080, 20, 100, 30)
        self.widgetContents['QLineEdit']['SimulationTime'].setText('10000')
        self.widgetContents['QLabel']['SimulationTime'] = QLabel('Simulation time [s]',self.centralWidget)
        self.widgetContents['QLabel']['SimulationTime'].setGeometry(930, 20, 150, 30)
        self.widgetContents['QLineEdit']['InitialTemperature'] = QLineEdit('',self.centralWidget)
        self.widgetContents['QLineEdit']['InitialTemperature'].setGeometry(1080, 20+30, 100, 30)
        self.widgetContents['QLineEdit']['InitialTemperature'].setText('100')
        self.widgetContents['QLabel']['InitialTemperature'] = QLabel('Initial Temp. [°C]',self.centralWidget)
        self.widgetContents['QLabel']['InitialTemperature'].setGeometry(930, 20+30, 150, 30)
        self.widgetContents['QCheckBox'] = {}
        self.widgetContents['QCheckBox']['1D'] = QCheckBox('1D',self.centralWidget)
        self.widgetContents['QCheckBox']['2D'] = QCheckBox('2D',self.centralWidget)
        self.widgetContents['QCheckBox']['1D'].stateChanged.connect(self.checkbox1Dtoggled)
        self.widgetContents['QCheckBox']['2D'].stateChanged.connect(self.checkbox2Dtoggled)
        self.widgetContents['QCheckBox']['1D'].move(700, 20+2*35)
        self.widgetContents['QCheckBox']['2D'].move(760, 20+2*35)
        self.widgetContents['QCheckBox']['2D'].toggle()
        self.widgetContents['QLabel']['Temp_xyt'] = QLabel(' ',self.centralWidget)
        self.widgetContents['QLabel']['Temp_xyt'].setGeometry(500-150,300, 100, 20)
        self.widgetContents['QLabel']['Temp_xyt2'] = QLabel(' ',self.centralWidget)
        self.widgetContents['QLabel']['Temp_xyt2'].setGeometry(1300-150,300, 100, 20)
        self.widgetContents['QLabel']['Temp_average'] = QLabel(' ',self.centralWidget)
        self.widgetContents['QLabel']['Temp_average'].setGeometry(1300-150,300+30, 100, 20)
        self.widgetContents['QLabel']['Temp1'] = QLabel('Temperature [°C] = ',self.centralWidget)
        self.widgetContents['QLabel']['Temp1'].setGeometry(500-300,300, 150, 20)
        self.widgetContents['QLabel']['Temp2'] = QLabel('Temperature [°C] = ',self.centralWidget)
        self.widgetContents['QLabel']['Temp2'].setGeometry(1300-300,300, 150, 20)
        self.widgetContents['QLabel']['TempAVG'] = QLabel('Average Temp [°C] = ',self.centralWidget)
        self.widgetContents['QLabel']['TempAVG'].setGeometry(1300-300,300+30, 150, 20)
       # self.widgetContents['QTextEdit'][]
       
        self.widgetContents['QCheckBox']['BC_ConvTop'] = QCheckBox('Top : convection = h [W/(m²K)]',self.centralWidget)
        self.widgetContents['QCheckBox']['BC_ConvTop'].setGeometry(50, 50, 270, 25)
        self.widgetContents['QCheckBox']['BC_ConvTop'].toggle()
        self.widgetContents['QLineEdit']['BC_ConvTop_h'] = QLineEdit('7',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvTop_h'].setGeometry(330, 50, 150, 25)
        self.widgetContents['QLabel']['BC_ConvTop_Temp'] = QLabel('Temperature [°C] = ',self.centralWidget)
        self.widgetContents['QLabel']['BC_ConvTop_Temp'].setGeometry(80, 50+25, 230, 25)
        self.widgetContents['QLineEdit']['BC_ConvTop_Temp'] = QLineEdit('20',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvTop_Temp'].setGeometry(330, 50+25, 150, 25)

        self.widgetContents['QCheckBox']['BC_ConvBottom'] = QCheckBox('Bottom : convection = h [W/(m²K)]',self.centralWidget)
        self.widgetContents['QCheckBox']['BC_ConvBottom'].setGeometry(50, 50  + 60, 270, 25)
        self.widgetContents['QCheckBox']['BC_ConvBottom'].toggle()
        self.widgetContents['QLineEdit']['BC_ConvBottom_h'] = QLineEdit('7',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvBottom_h'].setGeometry(330, 50 + 60, 150, 25)
        self.widgetContents['QLabel']['BC_ConvBottom_Temp'] = QLabel('Temperature [°C] = ',self.centralWidget)
        self.widgetContents['QLabel']['BC_ConvBottom_Temp'].setGeometry(80, 50+25  + 60, 230, 25)
        self.widgetContents['QLineEdit']['BC_ConvBottom_Temp'] = QLineEdit('20',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvBottom_Temp'].setGeometry(330, 50+25  + 60, 150, 25)

        self.widgetContents['QCheckBox']['BC_ConvLeft'] = QCheckBox('Left : convection = h [W/(m²K)]',self.centralWidget)
        self.widgetContents['QCheckBox']['BC_ConvLeft'].setGeometry(50, 50  + 2*60, 270, 25)
        self.widgetContents['QCheckBox']['BC_ConvLeft'].toggle()
        self.widgetContents['QLineEdit']['BC_ConvLeft_h'] = QLineEdit('7',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvLeft_h'].setGeometry(330, 50 + 2*60, 150, 25)
        self.widgetContents['QLabel']['BC_ConvLeft_Temp'] = QLabel('Temperature [°C] = ',self.centralWidget)
        self.widgetContents['QLabel']['BC_ConvLeft_Temp'].setGeometry(80, 50+25  + 2*60, 230, 25)
        self.widgetContents['QLineEdit']['BC_ConvLeft_Temp'] = QLineEdit('20',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvLeft_Temp'].setGeometry(330, 50+25  + 2*60, 150, 25)

        self.widgetContents['QCheckBox']['BC_ConvRight'] = QCheckBox('Right : convection = h [W/(m²K)]',self.centralWidget)
        self.widgetContents['QCheckBox']['BC_ConvRight'].setGeometry(50, 50  + 3*60, 270, 25)
        self.widgetContents['QCheckBox']['BC_ConvRight'].toggle()
        self.widgetContents['QLineEdit']['BC_ConvRight_h'] = QLineEdit('7',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvRight_h'].setGeometry(330, 50 + 3*60, 150, 25)
        self.widgetContents['QLabel']['BC_ConvRight_Temp'] = QLabel('Temperature [°C] = ',self.centralWidget)
        self.widgetContents['QLabel']['BC_ConvRight_Temp'].setGeometry(80, 50+25  + 3*60, 230, 25)
        self.widgetContents['QLineEdit']['BC_ConvRight_Temp'] = QLineEdit('20',self.centralWidget)
        self.widgetContents['QLineEdit']['BC_ConvRight_Temp'].setGeometry(330, 50+25  + 3*60, 150, 25)


    def is_number(self,inp):
        print('got ',inp)
        try:
            number = float(inp)
        except Exception:
            QtGui.QMessageBox.about(self, 'Error','Input can only be a number')
            pass
        return number

    def update_parameters(self):
        self.part = simulation_part(self.dimension,'ID1',self.is_number(self.widgetContents['QLineEdit']['InitialTemperature'].text()),self.is_number(self.widgetContents['QLineEdit']['SimulationTime'].text()))
        
        if self.widgetContents['QCheckBox']['BC_ConvTop'].isChecked():
            self.part.set_boundary_condition('T_env_top', deg2K( self.is_number(  self.widgetContents['QLineEdit']['BC_ConvTop_Temp'].text()  ) ))
            self.part.set_boundary_condition('h_top',self.is_number(  self.widgetContents['QLineEdit']['BC_ConvTop_h'].text()  ) )
        else:
            self.part.set_boundary_condition('T_env_top', 0.0)
            self.part.set_boundary_condition('h_top',0 )

        if self.widgetContents['QCheckBox']['BC_ConvBottom'].isChecked():
            self.part.set_boundary_condition('T_env_bottom', deg2K( self.is_number(  self.widgetContents['QLineEdit']['BC_ConvBottom_Temp'].text()  ) ))
            self.part.set_boundary_condition('h_bottom',self.is_number(  self.widgetContents['QLineEdit']['BC_ConvBottom_h'].text()  ) )
        else:
            self.part.set_boundary_condition('T_env_bottom', 0.0)
            self.part.set_boundary_condition('h_bottom',0 )
        
        if self.widgetContents['QCheckBox']['BC_ConvLeft'].isChecked():
            self.part.set_boundary_condition('T_env_left', deg2K( self.is_number(  self.widgetContents['QLineEdit']['BC_ConvLeft_Temp'].text()  ) ))
            self.part.set_boundary_condition('h_left',self.is_number(  self.widgetContents['QLineEdit']['BC_ConvLeft_h'].text()  ) )
        else:
            self.part.set_boundary_condition('T_env_left', 0.0)
            self.part.set_boundary_condition('h_left',0 )
            
        if self.widgetContents['QCheckBox']['BC_ConvRight'].isChecked():
            self.part.set_boundary_condition('T_env_right', deg2K( self.is_number(  self.widgetContents['QLineEdit']['BC_ConvRight_Temp'].text()  ) ))
            self.part.set_boundary_condition('h_right',self.is_number(  self.widgetContents['QLineEdit']['BC_ConvRight_h'].text()  ) )
        else:
            self.part.set_boundary_condition('T_env_right', 0.0)
            self.part.set_boundary_condition('h_right',0 )
        
    #    self.part.init_temperature(temp = self.is_number(self.widgetContents['QLineEdit']['InitialTemperature'].text()))
    #    self.part.init_simulation_time(self.is_number(self.widgetContents['QLineEdit']['SimulationTime'].text()))

    def log(self,tex):
        self.widgetContents['QLabel']['log'].setText(tex)

    def setting_up_visuals_for_1D(self):
        self.widgetContents['QSlider']['y'].setMinimum(0)
        self.widgetContents['QSlider']['y'].setMaximum(0)
        self.widgetContents['QSlider']['y'].setValue(0)
        self.widgetContents['QSlider']['y'].setTickPosition(QSlider.TicksBelow)
        self.widgetContents['QSlider']['y'].setTickInterval(1)
        self.widgetContents['graphs']['controuf'].set_xlabel('time [s]')
        self.widgetContents['QSlider']['y'].hide()
        self.widgetContents['QLabel']['ySlider'].hide()
        self.widgetContents['QLabel']['yValue'].hide()
        '''
        self.widgetContents['QCheckBox']['BC_ConvLeft'].hide()
        self.widgetContents['QLineEdit']['BC_ConvLeft_h'].hide()
        self.widgetContents['QLabel']['BC_ConvLeft_Temp'].hide()
        self.widgetContents['QLineEdit']['BC_ConvLeft_Temp'].hide()
        self.widgetContents['QCheckBox']['BC_ConvRight'].hide()
        self.widgetContents['QLineEdit']['BC_ConvRight_h'].hide()
        self.widgetContents['QLabel']['BC_ConvRight_Temp'].hide()
        self.widgetContents['QLineEdit']['BC_ConvRight_Temp'].hide()
        '''

    def setting_up_visuals_for_2D(self):
        self.widgetContents['QSlider']['y'].setMinimum(0)
        self.widgetContents['QSlider']['y'].setMaximum(self.part.sim_data['NY']-1)
        self.widgetContents['QSlider']['y'].setValue(0)
        self.widgetContents['QSlider']['y'].setTickPosition(QSlider.TicksBelow)
        self.widgetContents['QSlider']['y'].setTickInterval(int(self.part.sim_data['NY']/10))
        self.widgetContents['QSlider']['y'].show()
        self.widgetContents['QLabel']['ySlider'].show()
        self.widgetContents['QLabel']['yValue'].show()
        '''
        self.widgetContents['QCheckBox']['BC_ConvLeft'].show()
        self.widgetContents['QLineEdit']['BC_ConvLeft_h'].show()
        self.widgetContents['QLabel']['BC_ConvLeft_Temp'].show()
        self.widgetContents['QLineEdit']['BC_ConvLeft_Temp'].show()
        self.widgetContents['QCheckBox']['BC_ConvRight'].show()
        self.widgetContents['QLineEdit']['BC_ConvRight_h'].show()
        self.widgetContents['QLabel']['BC_ConvRight_Temp'].show()
        self.widgetContents['QLineEdit']['BC_ConvRight_Temp'].show()
        '''

    def checkbox1Dtoggled(self,state):
        if state == Qt.Checked:
            self.widgetContents['QCheckBox']['2D'].setChecked(False)
            self.dimension = 1
            self.setting_up_visuals_for_1D()
        else:
            self.widgetContents['QCheckBox']['2D'].setChecked(True)
            self.dimension = 2
            self.setting_up_visuals_for_2D()
    def checkbox2Dtoggled(self,state):
        if state == Qt.Checked:
            self.widgetContents['QCheckBox']['1D'].setChecked(False)
            self.dimension = 2
            self.setting_up_visuals_for_2D()
        else:
            self.widgetContents['QCheckBox']['1D'].setChecked(True)
            self.dimension = 1
            self.setting_up_visuals_for_1D()

    def moved_slider_x(self):
        
        x_ = int(self.widgetContents['QSlider']['x'].value())
        self.widgetContents['QLabel']['xValue'].setText(str(round(self.part.x[x_],2)))
        self.showTemperature_at_xyz()
        self.plot()

        
    def moved_slider_y(self):
       # self.plot_xy_t(int(self.widgetContents['QSlider']['x'].value()),int(self.widgetContents['QSlider']['y'].value()))
        y_ = int(self.widgetContents['QSlider']['y'].value())
        self.widgetContents['QLabel']['yValue'].setText(str(round(self.part.y[y_],2)))
        self.showTemperature_at_xyz()
        self.plot()

    def moved_slider_time(self):
        print('plotting now')
       # self.plot_contour_2D(int(self.widgetContents['QSlider']['time'].value()))
        self.showTemperature_at_xyz()
        t_ = int(self.widgetContents['QSlider']['time'].value())
        self.widgetContents['QLabel']['timeValue'].setText(str(round(self.part.time[t_],2)))
        self.plot()

    def showTemperature_at_xyz(self):
        t_ = int(self.widgetContents['QSlider']['time'].value())
        x_ = int(self.widgetContents['QSlider']['x'].value())
        y_ = int(self.widgetContents['QSlider']['y'].value())
        print('Trying to get temp at xyt : ', x_, y_, t_)
        print(self.part.Temperature[t_,x_,y_])
        text_ = str(self.part.Temperature[t_,x_,y_])
        print(text_)
        self.widgetContents['QLabel']['Temp_xyt'].setText(text_)
        self.widgetContents['QLabel']['Temp_xyt2'].setText(text_)
        self.widgetContents['QLabel']['Temp_average'].setText(str(self.part.average_Temperature[t_]))

    def simulate_and_plot(self):
        self.update_parameters()
        self.simulate()
        self.plot()
        self.moved_slider_time()
        self.moved_slider_x()
        if self.part.dimensions==2:
            self.moved_slider_y()
       # self.showTemperature_at_xyz()
    def plot(self):
        if self.part.dimensions==1:
            print('got here')
            self.plot_contour_1D()
        if self.part.dimensions==2:
            self.plot_contour_2D(int(self.widgetContents['QSlider']['time'].value()))
        self.plot_xy_t(int(self.widgetContents['QSlider']['x'].value()),int(self.widgetContents['QSlider']['y'].value()))
        #print(self.part.get_average_temperature(10))
    
        
    def simulate(self):
        self.log('Starting simulation ...')
        print('starting simulation')
        for i in range(self.part.sim_data['Nt']-1):
            self.part.simulate_time_step()
        self.log('Simulation finished')

    def plot_contour_1D(self):
        self.widgetContents['graphs']['controuf'].contourf(np.transpose(self.part.Temperature[:,:,0]))

    def plot_contour_2D(self,time):
        self.widgetContents['graphs']['controuf'].contourf(self.part.Temperature[time,:,:])
       # self.plot_point_extra( int(self.widgetContents['QSlider']['x'].value()) , int(self.widgetContents['QSlider']['y'].value()) )

    def plot_xy_t(self,x,y):
        self.widgetContents['graphs']['xy_t'].plot(self.part.time,self.part.Temperature[:,x,y])
        self.widgetContents['graphs']['xy_t'].plot_extra(self.part.time,self.part.average_Temperature)

    def write_to_csv(self):
        self.log('Starting export ...')
        # todo / ALSO MAKE 1D VERSION   !!!!
        with open('test.csv', mode='w') as test:
            employee_writer = csv.writer(test, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            if self.part.dimensions==1:
                first_line = ['x']
            if self.part.dimensions==2:
                first_line = ['x', 'y']
            for i in range(self.part.sim_data['Nt']):           
                first_line.append('Temp at t='+str(self.part.sim_data['Delta_T']*i)+'s')
                
            employee_writer.writerow(first_line)
            for i in range(self.part.sim_data['NX']):
                for j in range(self.part.sim_data['NY']):
                    if self.part.dimensions==1:
                        line_temp = [str(self.part.sim_data['Delta_X']*i)]
                    if self.part.dimensions==2:
                        line_temp = [str(self.part.sim_data['Delta_X']*i),str(self.part.sim_data['Delta_Y']*j)]
                    for k in range(self.part.sim_data['Nt']-1):
                        line_temp.append(str(self.part.Temperature[k,i,j]))
                    employee_writer.writerow(line_temp)
        self.log('export finished')
        

##############################################################################
##############################################################################
app = QApplication(sys.argv)
app.setStyleSheet('QMainWindow{background-color: darkgray;}')
app.setStyleSheet('QMainWindow{background-color: white;}')
myWidget = Gui()
myWidget.show()
app.exec_()
