import os
import numpy 

from amuse.units import units,quantities

from matplotlib import pyplot

from blackbody import energy_flux2,total_bolometric_flux,B_lambda

filter_data=dict()

filters=[ "bess-u.pass",
          "bess-b.pass",
          "bess-v.pass",
          "bess-r.pass",
          "bess-i.pass" ] # Bessell 1990 filters

data_directory = os.path.dirname(os.path.abspath(__file__))


for fil in filters:

  f=open( os.path.join(data_directory,fil) ,"r")
  
  lines=f.readlines()
  
  lam=[] | units.angstrom
  throughput= []
  
  for l in lines:
    line=l.split()
    lam.append( float(line[0]) | units.angstrom )
    throughput.append( float(line[1]) )
  throughput=numpy.array(throughput)
  
  filter_data[ fil ]=dict(wavelength=lam, throughput=throughput)

def filter_band_flux(fdata, source, lambdas=None):
  if fdata.__class__.__name__=="str":
    fdata=filter_data[fdata]
  if lambdas is None:
    lambdas=fdata['wavelength']
    f=fdata['throughput']
  else:
    xp=fdata['wavelength']
    fp=fdata['throughput']
    f=numpy.interp(lambdas.value_in(units.angstrom),xp=xp.value_in(units.angstrom),fp=fp,left=0.,right=0.)
  src=f*source(lambdas)
  src=quantities.to_quantity(src)
  return numpy.trapz(src.number,x=lambdas.number) | (src.unit*lambdas.unit)

def filter_band_lambda(fdata):
  return (filter_band_flux(fdata, lambda x: x)/ filter_band_flux(fdata, lambda x: 1.)).in_(units.angstrom)

def plot_filters():

  n=1000
  for fil in filters:
    lam=(3000.+(9200.-3000.)*numpy.array(range(n+1))/n) | units.angstrom
    xp=filter_data[fil]['wavelength']
    fp=filter_data[fil]['throughput']
    f=numpy.interp(lam.value_in(units.nano(units.m)),xp=xp.value_in(units.nano(units.m)),fp=fp,left=0.,right=0.)
  
    pyplot.plot(lam.value_in(units.nano(units.m)), f)

  pyplot.show()

if __name__=="__main__":

  plot_filters()
  
  T=5000. | units.K
  
  fb=filter_band_flux(filter_data['bess-u.pass'], lambda x: B_lambda(x,T)).in_(units.W/units.m**2)
  
  print fb
  print (fb*(1.| units.RSun)**2/(1.| units.AU)**2).in_(units.W/units.m**2)
  print (energy_flux2(5778.|units.K)*(1.| units.RSun)**2/(1.| units.AU)**2).in_(units.W/units.m**2)  
