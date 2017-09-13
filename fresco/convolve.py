import numpy
from numpy.fft import fft2, ifft2

def convolve_fft(image1, image2, MinPad=True, pad=True):
  """ Not so simple convolution """
  #Just for comfort:
  FFt = fft2
  iFFt = ifft2
    
	#The size of the images:
  r1,c1 = image1.shape
  r2,c2 = image2.shape
  
	#MinPad results simpler padding,smaller images:
  if MinPad:
    r = r1+r2
    c = c1+c2
  else:
    #if the Numerical Recipies says so:
    r = 2*max(r1,r2)
    c = 2*max(c1,c2)
    
	#For nice FFT, we need the power of 2:
  if pad:
    pr2 = int(numpy.log(r)/numpy.log(2.0) + 1.0 )
    pc2 = int(numpy.log(c)/numpy.log(2.0) + 1.0 )
    rOrig = r
    cOrig = c
    r = 2**pr2
    c = 2**pc2
	
#numpy fft has the padding built in, which can save us some steps
#here. The thing is the s(hape) parameter:
#fftimage = FFt(image1,s=(r,c)) * FFt(image2,s=(r,c))
  fftimage = FFt(image1, s=(r,c))*FFt(image2[::-1,::-1],s=(r,c))
  
  if pad:
    return (iFFt(fftimage))[:rOrig,:cOrig].real
  else:
    return (iFFt(fftimage)).real

def Convolve(
        image,
        kernel,
        ):
    xpad,ypad=kernel.shape
    result = convolve_fft(image, kernel)
    return result[xpad/2:-xpad/2,ypad/2:-ypad/2]


def gaussian_filter(image, sigma, order=0, truncate=4.0):
    size=sigma*truncate
    N=numpy.ceil(size)*2
    x,y=numpy.mgrid[-size:size:-1j*N,-size:size:-1j*N]
    kernel=numpy.exp(-(x**2+y**2)/sigma**2)
    kernel=kernel/kernel.sum()
    #~ from matplotlib import pyplot
    #~ f=pyplot.figure()
    #~ pyplot.imshow(Convolve(image,kernel))
    #~ pyplot.show()
    print kernel.shape
    return Convolve(image,kernel)
