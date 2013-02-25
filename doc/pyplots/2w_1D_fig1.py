import pylab as pylab

fff=open('2w_1D_fig1.dat','r')
xx=[]
yy=[]
for ii in fff.readlines():
    aa=ii.split()
    xx.append(float(aa[0])); yy.append(float(aa[1]))

fff.close()
pylab.plot(xx,yy,'bo')
pylab.show()

