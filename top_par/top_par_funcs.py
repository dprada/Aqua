## Aux. function to load new topologies:

user_topol=[]

def add_topol(self,file_topol,verbose=False):

   newtop= __import__(file_topol)
   
   self.residue[newtop.residue_name]=newtop.residue_name
   self.residue_type[newtop.residue_name]=newtop.residue_type
   self.residue_atoms[newtop.residue_name]=[]
   
   new_list_ats=[]
   for ii,jj in newtop.atoms.iteritems():
      try:
         new_list_ats.append(self.atom[ii])
      except:
         xx='nat'+ii
         self.atom[ii]=xx
         self.atom_type[xx]=jj
         new_list_ats.append(xx)
      
   new_cov_bonds=[]
      
   for ii in newtop.covalent_bonds:
      new_cov_bonds.append([self.atom[ii[0]],self.atom[ii[1]]])

   for ii,jj in newtop.charge.iteritems():
      self.charge[self.atom[ii]]=jj

   for ii in newtop.donors:
      try:
         self.donors_exception[self.atom[ii]][newtop.residue_name]=['Always',True]
      except:
         self.donors_exception[self.atom[ii]]={ newtop.residue_name : [ 'Always'   , True  ]}

   for ii in newtop.acceptors:
      try:
         self.acceptors_exception[self.atom[ii]][newtop.residue_name]=['Always',True]
      except:
         self.acceptors_exception[self.atom[ii]]={ newtop.residue_name : [ 'Always'   , True  ]}


   self.residue_atoms[newtop.residue_name]=new_list_ats
   self.covalent_bonds[newtop.residue_name]=new_cov_bonds

   if verbose:
      print '# New topology:',newtop.residue_name
      print '#',newtop.residue_type,'with', len(newtop.atoms),'atoms and', len(newtop.covalent_bonds),'covalent bonds.'

