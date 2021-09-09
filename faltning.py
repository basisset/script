import numpy as np
import sys

x_in = np.zeros(500000)
y1_in = np.zeros(500000)
ytot = np.zeros(500000)
x = np.zeros(500000)
y = np.zeros(500000)

npeaks,start,end,fwhm = sys.stdin.readline().strip().split()
npeaks = int(npeaks)
start = float(start)
end = float(end)
fwhm = float(fwhm)

c = 2.772588722/(fwhm**2)

for i in range(npeaks):
	x_in[i],y1_in[i] = sys.stdin.readline().strip().split() 
	x_in[i] = float(x_in[i])
	y1_in[i] = float(y1_in[i])
	ytot[i] = y1_in[i]

step = 0.005
i = 0
x[0] = start
A = 1.0
while x[i] < end:
	i = i+1
	x[i] = x[i-1] + step
	y[i] = 0
	for j in range(npeaks):
		y[i]=y[i]+ytot[j]*A*np.exp(-c*(x_in[j]-x[i])**2)
	print("{} {}".format(x[i],y[i]))





