from Tkinter import *
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
from pymol import cmd
from numpy import *
import aqua as aqua
from os import remove as remove_file
import string
import random

def __init__(self):
   self.menuBar.addmenu('Network', 'Network Menu')

   self.menuBar.addmenuitem('Network', 'command',
                      'Add Network',
                      label='Add Network',
                     command = lambda : load_net())

   self.menuBar.addcascademenu('Network', 'Edges')

   self.menuBar.addmenuitem('Edges', 'command',
                      'Show Edges',
                      label='Show Edges',
                     command = lambda : show_edges_net())

   self.menuBar.addmenuitem('Edges', 'command',
                      'Hide Edges',
                      label='Hide Edges',
                     command = lambda : hide_edges_net())


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
   cmd.set("connect_mode", 1)
   nfile='/tmp/fake'+''.join(random.choice(string.ascii_uppercase) for i in range(6))+'.xyz'
   f = open(nfile, 'wb')
   f.write(str(net.num_nodes)+"\n")
   #f.write("FAKE\n")
   for ii in net.node:
      f.write("X\t%f\t%f\t%f\n" % (ii.coors[0],ii.coors[1],ii.coors[2]))
   f.close()
   cmd.load(nfile,name_net)
   remove_file(nfile)
   for ii in xrange(net.num_nodes):
      selec=name_net+' & id '+str(ii+1)
      #cmd.alter(selec,'name='+net.node[ii].label)
      cmd.alter(selec,'resi='+str(net.node[ii].cluster))
      cmd.alter(selec,'b='+str(net.node[ii].weight))
      cmd.alter(selec,'q='+str(net.node[ii].att1))
      cmd.alter(selec,'ID='+str(ii))
   #for ii in range(net.num_nodes):
   #   for jj in net.node[ii].link.keys():
   #      cmd.bond(name_net+' & id '+str(ii),name_net+' & id '+str(jj))

#def show_edges_net():
   

#def load_net():
#   name_net= tkSimpleDialog.askstring('Loading Network',
#                                                      'Please enter the name of the new object:')
#   file_coors= open_it()
#   net=aqua.network(file_coors)
#   contador=10
#   incremento=net.num_nodes/10
#   voy=incremento
#   #print '# Loading nodes:'
#   for ii in range(net.num_nodes):
#      #if ii==voy:
#      #   print '#  ',contador,'%'
#      #   contador=contador+10
#      #   voy+=incremento
#      #cmd.pseudoatom(name_net,pos=net.node[ii].coors,name=net.node[ii].label,resi=net.node[ii].cluster,b=net.node[ii].weight,q=net.node[ii].att1)
#      cmd.pseudoatom(name_net,pos=net.node[ii].coors)
#   #contador=10
#   #voy=incremento
#   #for ii in range(net.num_nodes):
#   #   if ii==voy:
#   #      print '#  ',contador,'%'
#   #      contador=contador+10
#   #      voy+=incremento
#   #   for jj in net.node[ii].link.keys():
#   #      cmd.bond('id '+str(ii),'id '+str(jj))
# 
# 
# 
# 
# 
# 
##cmd.pseudoatom(name_net,pos=net.node[ii].coors,'name='+'"'+net.node[ii].label+'"','resi='+str(net.node[ii].cluster),'b='+str(net.node[ii].weight),'q='+str(net.node[ii].att1))
##cmd.spectrum('b','blue_red')
##cmd.spectrum('q','rainbow')
##cmd.select('"'+dict1[ii]+'"','resn "'+dict1[ii]+'"')
##selec='id '+str(ii)
##cmd.alter(selec,'name='+'"'+net.node[ii].label+'"')
###cmd.label(selec,'"'+net.node[ii].label+'"')
##cmd.alter(selec,'resi='+str(net.node[ii].cluster))
##cmd.alter(selec,'b='+str(net.node[ii].weight))
##cmd.alter(selec,'q='+str(net.node[ii].att1))

