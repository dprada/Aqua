from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
from pymol import cmd
from numpy import *

def __init__(self):
   self.menuBar.addmenu('Network', 'Network Menu')

   self.menuBar.addmenuitem('Network', 'command',
                      'Add Network',
                      label='Add Network',
                     command = lambda : load_net())

   self.menuBar.addmenuitem('Network', 'command',
                     'Add Links',
                      label='Add Links',
                     command = lambda : load_links())

   self.menuBar.addmenuitem('Network', 'command',
                     'Add Weights',
                      label='Add Weights',
                     command = lambda : load_weights())

   self.menuBar.addmenuitem('Network', 'command',
                     'Add sets',
                      label='Add sets',
                     command = lambda : load_sets())

   self.menuBar.addmenuitem('Network', 'command',
                     'Add labels',
                      label='Add labels',
                     command = lambda : load_labels())

def load_net():
   file_coors= open_it()
   i=0
   for line in open(file_coors,'r'):
      coors_net=map(float,line.split())
      i=i+1
      cmd.pseudoatom('net',pos=coors_net,name='node',resi=i,chain='A')

def load_links():
   file_links= open_it()
   for line in open(file_links,'r'):
      edge=line.split()
      for i in range (0,len(edge)): edge[i]='resi '+edge[i]
      cmd.bond(edge[0],edge[1])

def load_weights():
   file_weights= open_it()
   i=0
   for line in open(file_weights,'r'):
      i=i+1
      selec='resi '+str(i)
      ocup=line.split()
      cmd.alter(selec,"b="+ocup[0])

   cmd.spectrum("b", 'blue_red')

def load_sets():
   dict1={'1':'00-0',
          '2':'00-1',
          '3':'00-2',
          '4':'10-0',
          '5':'02-0',
          '6':'10-1',
          '7':'02-1',
          '8':'10-2',
          '9':'02-2',
          '10':'12-0',
          '11':'12-1',
          '12':'12-2',
          '13':'00-K',
          '14':'10-K',
          '15':'02-K',
          '16':'12-K',
          '17':'12-K1',
          '18':'00-1K',
          '19':'10-1K',
          '20':'01-1K',
          '21':'12-1K',
          '22':'0K-0',
          '23':'0K-1',
          '24':'0K-2',
          '25':'1K-0',
          '26':'1K-1',
          '27':'1K-2',
          '28':'K0-0',
          '29':'K0-1',
          '30':'K0-2',
          '31':'K2-0',
          '32':'K2-1',
          '33':'K2-2'}
   file_sets= open_it()
   i=0
   for line in open(file_sets,'r'):
      i=i+1
      selec='resi '+str(i)
      ocup=line.split()
      cmd.alter(selec,"resn="+'"'+dict1[ocup[0]]+'"')
      cmd.alter(selec,"q="+ocup[0])
      
   for ii in sorted(dict1.iterkeys()):
      cmd.select('"'+dict1[ii]+'"','resn "'+dict1[ii]+'"')
   cmd.spectrum("q", 'rainbow')
   
def load_labels():
   file_labels= open_it()
   i=0
   for line in open(file_labels,'r'):
      i=i+1
      selec='resi '+str(i)
      txt=line.split()
      cmd.label(selec,txt[0])

#def load_mss():
#   file_mss= open_it()
#   i=0
#   for line in open(file_mss,'r')
   

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
