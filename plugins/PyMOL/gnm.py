from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
from pymol.cgo import *
from pymol import cmd
from numpy import *


class net_model():
   def __init__(self):
      self.num_nodes=0
      self.num_modes=0
      self.modes=zeros(shape=(self.num_modes,self.num_nodes))
      self.nodes=zeros(shape=(self.num_nodes))


def __init__(self):

   net=net_model()

   self.menuBar.addmenu('GNM', 'GNM Menu')

   self.menuBar.addmenuitem('GNM', 'command',
                      'Contact Map',
                      label='Contact Map',
                     command = lambda : load_cmap())

   self.menuBar.addmenuitem('GNM', 'command',
                      'Load Modes',
                      label='Load Modes',
                     command = lambda : load_modes(net))

   self.menuBar.addmenuitem('GNM', 'command',
                      'Choose Mode',
                      label='Choose Mode',
                     command = lambda : plot_mode(net))

   self.menuBar.addmenuitem('GNM', 'command',
                      'check',
                      label='Help',
                     command = lambda : check(net))

def load_cmap():
   file_cmap= open_it()

   m=cmd.get_model('all')
   trans={}
   for ii in range(len(m.atom)):
      trans[m.atom[ii].id]=ii


   ElasticNet = [ BEGIN, LINES, COLOR, 1.0, 0.4, 0.0 ]

   for line in open(file_cmap,'r'):
      contact=map(float,line.split())

      aa=trans[contact[0]]
      bb=trans[contact[1]]

      ElasticNet.append(VERTEX)
      ElasticNet.append(m.atom[aa].coord[0])
      ElasticNet.append(m.atom[aa].coord[1])
      ElasticNet.append(m.atom[aa].coord[2])
      ElasticNet.append(VERTEX)
      ElasticNet.append(m.atom[bb].coord[0])
      ElasticNet.append(m.atom[bb].coord[1])
      ElasticNet.append(m.atom[bb].coord[2])



   ElasticNet.append(END)

   cmd.load_cgo(ElasticNet,'ElasticNet')

def load_modes(net):
   file_modes= open_it()
   f = open(file_modes, 'r')

   m=cmd.get_model('all')
   trans={}
   for ii in range(len(m.atom)):
      trans[m.atom[ii].id]=ii


   data=f.readline().split()

   net.num_modes=int(data[0])
   net.num_nodes=int(data[2])

   net.modes=zeros(shape=(net.num_modes,net.num_nodes))
   net.nodes=zeros(shape=(net.num_nodes))

   for ii in range(net.num_modes):
      data=f.readline().split()
      for jj in range(net.num_nodes):
         data=f.readline().split()
         net.modes[ii][jj]=data[1]
         net.nodes[jj]=data[0]

def plot_mode(net):

   mode= tkSimpleDialog.askstring('Modes of the GNM',
                                                      'Please enter the index of the mode:')
   mode=int(mode)

   for ii in range(net.num_nodes):
      selec='id '+str(net.nodes[ii])
      cmd.alter(selec,"q="+str(net.modes[mode][ii]))
   maxx=ma.maximum(net.modes[mode][:])
   cmd.spectrum("q", 'red_white_blue',minimum=str(-maxx),maximum=str(maxx))

def check(net):

#   m=cmd.get_model('all')
#   trans={}
#   for ii in range(len(m.atom)):
#      trans[m.atom[ii].pdb_index]=ii
#   print m.atom[0].id,m.atom[0].index

   print net.nodes

def open_it():
   f_name = tkFileDialog.askopenfilename()
   return f_name










###import tkFileDialog
### 
#### examples:
### 
###def open_it():
### 
###filename = tkFileDialog.askopenfilename()
### 
###print filename # test
### 
### 
### 
###def save_it():
### 
###filename = tkFileDialog.askopenfilename()
### 
###print filename # test
### 
### 
### 
###def save_as():
### 
###filename = tkFileDialog.asksaveasfilename()
### 
###print filename # test
