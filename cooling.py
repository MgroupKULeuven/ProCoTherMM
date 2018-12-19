# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 09:29:20 2018

@author: robot
"""
from PyQt4 import QtGui, uic
import sys
import numpy as np
import matplotlib.pyplot as plt
import time
# material properties

class Example(QtGui.QWidget):
    def __init__(self):
        super(Example, self).__init__()
        self.rho_ = 2700; # kg/m3
        self.lambda_= 237; # W/mK
        self.cp_ = 880;
        self.solver = dict()
        
        self.ui = uic.loadUi('Gui.ui',self)        
        self.ui.pushButton.clicked.connect(self.solver_FD_matrix)
        
        self.init_stuff()
        self.init_geometry()
        self.init_solver_conditions()
        self.init_matrix_construction()
        self.init_boundary_conditions()
        self.ui.label.setText('finished initializing')

        
        
        #self.setGeometry(100,100,1000,600)
        self.ui.show()
        
        
    def init_stuff(self):
        self.a = self.lambda_/(self.rho_*self.cp_)
        self.T0 = 600
        self.HTCleft =0
        self.Tleft=15
        self.HTCright = 3000
        self.Tright = 15
        self.HTCtop=800
        self.Ttop = 15
        self.HTCbottom = 800
        self.Tbottom = 15
        self.t_cool = 50
        
    def init_geometry(self):
        self.LengthX = 0.1
        self.LengthY = 0.1;
        self.A = 1;
        
    def init_solver_conditions(self):
        self.solver['tau'] = 0.005 # s
        self.solver['e'] = 0.005 #m
        self.solver['f'] = 0.005; 
        self.solver['CFL'] = self.a*self.solver['tau']/(self.solver['e']**2)
        
    def init_matrix_construction(self):
        self.NOFtimesteps=int(self.t_cool/self.solver['tau'])
        self.NOFXlayers = int(self.LengthX/self.solver['e'])
        self.NOFYlayers = int(self.LengthY/self.solver['f'])
        self.T = np.ones((self.NOFYlayers+1,self.NOFXlayers+1,self.NOFtimesteps))*self.T0
        self.Temp = np.ones((self.NOFYlayers,self.NOFXlayers,self.NOFtimesteps))*self.T0
        self.BC_type = np.zeros((2,self.NOFXlayers)) # eerste rij convectie of niet, tweede rij geleiding of niet
        self.BC_type[0,0] = 1
        self.BC_type[0,-1]=1
        self.BC_type[1,:] = 1
        print self.BC_type
        
    def init_boundary_conditions(self):
        self.BC = np.zeros((self.NOFYlayers,self.NOFXlayers,4))
        # eerste tab --> convectie links
        # tweede tab --> convectie rechts
        # derde tab --> convectie onder
        # vierde tab --> convectie boven
        self.BC[:,0,0] = 1
        self.BC[:,-1,1] = 1
        self.BC[-1,:,2] = 1
        self.BC[0,:,3] = 1
        
    def solver_FD_matrix(self):
        start=time.time()
        print 'starting solver', 
        self.T_array = dict()
        self.T_array[str(0)] = np.ones((self.NOFYlayers*self.NOFXlayers,1))*self.T0
        A__ = np.zeros((self.NOFYlayers*self.NOFXlayers,self.NOFYlayers*self.NOFXlayers))
        B__ = np.zeros((self.NOFYlayers*self.NOFXlayers,1))
        N_ = self.NOFYlayers
        mul = self.solver['tau']/(self.A*self.solver['e']*self.rho_*self.cp_)
        mul2 = -self.lambda_/self.solver['e']*self.A*mul
        for j in range(self.NOFXlayers):
            for k in range(self.NOFYlayers):
                A__[j+N_*k, j+N_*k] +=  self.BC[k,j,0] * self.HTCleft*self.A*mul
                B__[j+N_*k] = -self.Tleft * self.BC[k,j,0] * self.HTCleft*self.A*mul
                A__[j+N_*k, j+N_*k] +=  self.BC[k,j,1] * self.HTCright*self.A*mul
                B__[j+N_*k] = -self.Tright * self.BC[k,j,1] * self.HTCright*self.A*mul
                A__[j+N_*k, j+N_*k] +=  self.BC[k,j,2] * self.HTCbottom*self.A*mul
                B__[j+N_*k] = -self.Tbottom * self.BC[k,j,2] * self.HTCbottom*self.A*mul
                A__[j+N_*k, j+N_*k] +=  self.BC[k,j,3] * self.HTCtop*self.A*mul
                B__[j+N_*k] = -self.Ttop * self.BC[k,j,3] * self.HTCtop*self.A*mul
                if (j<self.NOFXlayers-1):
                    A__[j+N_*k, j+N_*k] += mul2
                    A__[j+N_*k, j+1+N_*k] += -mul2
                if (j>0):
                    A__[j+N_*k, j+N_*k] += mul2
                    A__[j+N_*k, j-1+N_*k] += -mul2
                if (k<self.NOFYlayers-1):
                    A__[j+N_*k, j+N_*k] += mul2
                    A__[j+N_*k, j+N_*(k+1)] += -mul2
                if (k>0):
                    A__[j+N_*k, j+N_*k] += mul2
                    A__[j+N_*k, j+N_*(k-1)] += -mul2
        print A__
        for i in range(self.NOFtimesteps-1):
           # print np.dot( A__  , self.T_array[str(i)] )
            self.T_array[str(i+1)] = self.T_array[str(i)] + np.dot( A__  , self.T_array[str(i)] ) + B__
                    
        print 'solver done', time.time()-start
        print self.T_array[str(500)]
        for k in range(self.NOFYlayers):       
            self.Temp[k,:,500] =self.T_array[str(500)][N_*k:N_*k+(N_),0]
        print self.Temp[:,:,500]
        print self.T_array[str(500)][N_*k:N_*k+(N_),0]
        plt.figure(1)
        plt.contour(self.Temp[:,:,500])
        plt.show()
    def solver_FD(self):
        start=time.time()
        print 'starting solver'
       # self.ui.label.setText('starting solver')
        mul = self.solver['tau']/(self.A*self.solver['e']*self.rho_*self.cp_)
        mul2 = -self.lambda_/self.solver['e']*self.A*mul
        for i in range(self.NOFtimesteps-1):
            print i
            for j in range(self.NOFXlayers):
                for k in range(self.NOFYlayers):                           
                    QconvL = self.BC[k,j,0] * self.HTCleft*self.A*(self.T[k,j,i]-self.Tleft)*mul
                    QconvR = self.BC[k,j,1] * self.HTCright*self.A*(self.T[k,j,i]-self.Tright)*mul
                    QconvB = self.BC[k,j,2] * self.HTCbottom*self.A*(self.T[k,j,i]-self.Tbottom)*mul
                    QconvT = self.BC[k,j,3] * self.HTCtop*self.A*(self.T[k,j,i]-self.Ttop)*mul
                    # conduction
                    QinX = (j<self.NOFXlayers-1) * mul2*(self.T[k,j,i]-self.T[k,j+1,i])
                    QoutX = (j>0) * mul2*(self.T[k,j-1,i]-self.T[k,j,i])
                    QinY = (k<self.NOFYlayers-1) * mul2*(self.T[k,j,i]-self.T[k+1,j,i])
                    QoutY = (k>0) * mul2*(self.T[k-1,j,i]-self.T[k,j,i])
                    Qnet = QinX + QinY -QoutX - QoutY - QconvL - QconvR - QconvT - QconvB
                     # time integration of heat equation Euler integration
                    self.T[k,j,i+1] = self.T[k,j,i] + Qnet 
        print 'solver done', time.time()-start
        plt.figure(1)
        plt.contour(self.T[1:-2,1:-2,500])
        plt.show()
       # self.ui.label.setText('solver done')
def main():
    app = QtGui.QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
    
if __name__ =='__main__':
    main()


