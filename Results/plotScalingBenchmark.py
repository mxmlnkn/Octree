from numpy import *
from matplotlib.pyplot import *
import argparse
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser()
parser.add_argument("file", help="Filename") # e.g. benchmark-yeesolver.dat
parser.add_argument("-n", "--measurements", help="Number of Measurements in a 1-column file. Measurements are preceded by one line containing number of cores used")
# e.g. 24
args = parser.parse_args()

data = genfromtxt( args.file, comments='#' )
ncols = int(args.measurements)
print "Len: ",size(data),", ncols: ",ncols
if len(data) == size(data):
    assert(len(data) % (ncols+1) == 0)
    data = data.reshape( len(data) / (ncols+1), ncols+1 )
    data = data[:,0:-4].transpose()

fig = figure( figsize=(16,6.5) );
subplot(131)
xlabel('Cores N')
ylabel('Execution time for one timestep / s')
#xscale('log')
#yscale('log')
x    = data[0,:]
y    = empty( len( x ) )
yerr = empty( len( x ) )
for i in range( len( x ) ):
    y[i]   = data[1:,i].mean()
    yerr[i] = data[1:,i].std()
errorbar( x,y, yerr=yerr, fmt='o' )

subplot(132)
xlabel('Cores N')
ylabel('Speedup')
#xscale('log')
#yscale('log')
xlim(0.5,196)

def fitfunc( N, a, e, b ):
    return a*N**e + b
popt, pcov = curve_fit(fitfunc, x, y[0]/y, sigma=yerr/y * y[0]/y, p0=[1.0,1.0,0.0] ) # 
perr  = np.sqrt(np.diag(pcov))
print "Sp(N) = (",popt[0],"+-",perr[0],") * N ** (",popt[1],"+-",perr[1],")"
fitx = exp( linspace( log(min(x)), 1.1*log(max(x)), 100 ) )
plot( fitx, fitfunc( fitx, *popt ), 'r-', label=r'Exponential Fit: $%.3f'%popt[2]+'+%.3f'%popt[0]+"\cdot N^{"+'%.3f'%popt[1]+"}$" )

def gustafson( N, a ):
    return a + N*(1.0-a)
poptg, pcovg = curve_fit(gustafson, x, y[0]/y,  p0=[0.5] ) # sigma=yerr/y * y[0]/y,
plot( fitx, gustafson( fitx, poptg[0] ), 'g-', label=r"Gustafson's Law: $"+'%.3f'%poptg[0]+' + (1-'+'%.3f'%poptg[0]+')\cdot N$' )

def amdahl( N, b ):
    return 1.0 / ( b + (1.0-b)/N )
popta, pcova = curve_fit(amdahl, x, y[0]/y,  p0=[0.5] ) # sigma=yerr/y * y[0]/y,
plot( fitx, amdahl( fitx, popta[0] ), 'k-', label=r"Amdahl's Law: $\frac{1}{%.3f"%popta[0]+' + \frac{1-%.3f'%popta[0]+'}{N}$' )

errorbar( x, y[0]/y, yerr = yerr/y * y[0]/y, fmt='bo' )
legend(loc='best')

subplot(133)
xlabel('Cores N')
ylabel('Parallel Efficiency')
#xscale('log')
#yscale('log')
xlim(0.5,196)

plot( fitx, fitfunc( fitx, *popt ) / fitx, 'r-', label='Exponential Fit' )

def fitfunc2( N, t_comm, t_work_0 ):
    return 1.0 * ( t_comm/t_work_0 + 1.0 ) / ( N*t_comm/t_work_0 + 1.0 )
popt2, pcov2 = curve_fit( fitfunc2, x, y[0]/y / x, sigma=yerr/y * y[0]/y / x ) # 
perr2  = np.sqrt(np.diag(pcov2))
#plot( fitx, fitfunc2( fitx, popt2[0], popt2[1] ), label=r"$\frac{t_\mathrm{comm}+t_\mathrm{work}^0 }{ t_\mathrm{comm} + \frac{t_\mathrm{work}^0}{N} }$" )

plot( fitx, gustafson( fitx, poptg[0] ) / fitx, 'g-', label=r"Gustafson's Law" )

errorbar( x, y[0]/y / x, yerr = yerr/y * y[0]/y / x, fmt='bo' )
legend(loc='best')

tight_layout()
show()
