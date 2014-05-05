from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
from pymol import cmd
from numpy import *
import aqua as aqua

def __init__(self):
   self.menuBar.addmenu('Network', 'Network Menu')

   self.menuBar.addmenuitem('Network', 'command',
                      'Add Network',
                      label='Add Network',
                     command = lambda : load_net())

   self.menuBar.addmenuitem('Network', 'command',
                      'Help',
                      label='Help',
                     command = lambda : help_net())

   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add Links',
   #                   label='Add Links',
   #                  command = lambda : load_links())
   # 
   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add Weights',
   #                   label='Add Weights',
   #                  command = lambda : load_weights())
   # 
   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add sets',
   #                   label='Add sets',
   #                  command = lambda : load_sets())
   # 
   #self.menuBar.addmenuitem('Network', 'command',
   #                  'Add labels',
   #                   label='Add labels',
   #                  command = lambda : load_labels())

def open_it():
   f_name = tkFileDialog.askopenfilename()
   return f_name

def help_net():
   tkMessageBox.showinfo("Node \t --->\t Atom", "id \t --->\t id \n"+
                         "label \t --->\t name \n"+"cluster \t --->\t resi \n"+
                         "weight \t --->\t b \n"+ "att1 \t --->\t q")

def load_net():
   name_net= tkSimpleDialog.askstring('Loading Network',
                                                      'Please enter the name of the new object:')
   file_coors= open_it()
   net=aqua.network(file_coors)
   for ii in range(net.num_nodes):
      cmd.pseudoatom(name_net,pos=net.node[ii].coors)
      selec='id '+str(ii)
      cmd.alter(selec,'name='+'"'+net.node[ii].label+'"')
      #cmd.label(selec,'"'+net.node[ii].label+'"')
      cmd.alter(selec,'resi='+str(net.node[ii].cluster))
      cmd.alter(selec,'b='+str(net.node[ii].weight))
      cmd.alter(selec,'q='+str(net.node[ii].att1))
   for ii in range(net.num_nodes):
      for jj in net.node[ii].link.keys():
         cmd.bond('id '+str(ii),'id '+str(jj))
   #cmd.spectrum('b','blue_red')
   #cmd.spectrum('q','rainbow')
   #cmd.select('"'+dict1[ii]+'"','resn "'+dict1[ii]+'"')


