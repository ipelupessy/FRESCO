import os
import numpy

from amuse.units import units

data_directory = os.path.dirname(os.path.abspath(__file__))

f=open(os.path.join(data_directory,"ciexyz31.csv"),"r")

lines=f.readlines()

cielam=[] | units.nano(units.m)
x= []
y= []
z= []

for l in lines:
  line=l.split(",")
  cielam.append( float(line[0]) | units.nano(units.m) )
  x.append( float(line[1]) )
  y.append( float(line[2]) )
  z.append( float(line[3]) )
x=numpy.array(x)
y=numpy.array(y)
z=numpy.array(z)

xyz_data=dict()
xyz_data['x']=dict(wavelength=cielam, throughput=x)
xyz_data['y']=dict(wavelength=cielam, throughput=y)
xyz_data['z']=dict(wavelength=cielam, throughput=z)
