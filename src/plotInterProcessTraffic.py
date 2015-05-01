from numpy import *
from matplotlib.pyplot import *
import os.path
from scipy.optimize import curve_fit

orderings = [ "GrayCode", "Morton", "Hilbert", "Rows" ]
"""basestrings = [ "./output/2015-01-27_04-28/Octree-Setup-6_Initial-4_Max-Refinement-6_",
                "./output/2015-01-27_04-21/Octree-Setup-6_Initial-3_Max-Refinement-4_",
                "./output/2015-01-27_04-26/Octree-Setup-4_Initial-3_Max-Refinement-4_",
                "./output/2015-01-27_04-46/Octree-Setup-5_Initial-4_Max-Refinement-6_"]
setupnames = [ "Circle Shape Refinement",
               "Circle Shape Refinement, less cells",
               "Refinement at only one point",
               "Uniform randomly refined mesh" ]"""

"""basestrings = [ "./output/2015-01-27_07-11/Octree-Setup-6_Initial-3_Max-Refinement-4_",
                "./output/2015-01-27_07-05/Octree-Setup-4_Initial-4_Max-Refinement-6_",
                "./output/2015-01-27_07-19/Octree-Setup-5_Initial-3_Max-Refinement-6_",]
setupnames = [ "3D Sphere Refinement, less cells",
               "3D Refinement at only one point",
               "3D Uniform randomly refined mesh" ]"""
"""basestrings = [ "./output/2015-02-01_05-59_Octree_Benchmark_3D_and_2D/Octree-Setup-6_Initial-6_Max-Refinement-10_",
                "./output/2015-02-01_05-59_Octree_Benchmark_3D_and_2D/Quadtree-Setup-6_Initial-8_Max-Refinement-12_" ]
setupnames = [ u"Sphäre, 278272 Zellen",
               "Kreis, 93184 Zellen"  ]"""
"""basestrings = [ "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Quadtree-Setup-5_Initial-8_Max-Refinement-6_" ,
                "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Quadtree-Setup-6_Initial-7_Max-Refinement-12_",
                "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Octree-Setup-5_Initial-5_Max-Refinement-6_"   ,
                "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Octree-Setup-6_Initial-3_Max-Refinement-7_"   ]"""
basestrings = [ "./output/2015-02-14_21-46-31_Octree_Benchmark_3D_and_2D_preinterpolated/Quadtree-Setup-5_Initial-8_Max-Refinement-6_" ,
                "./output/2015-02-14_21-46-31_Octree_Benchmark_3D_and_2D_preinterpolated/Quadtree-Setup-6_Initial-7_Max-Refinement-12_",
                "./output/2015-02-14_21-46-31_Octree_Benchmark_3D_and_2D_preinterpolated/Octree-Setup-5_Initial-5_Max-Refinement-6_"   ,
                "./output/2015-02-14_21-46-31_Octree_Benchmark_3D_and_2D_preinterpolated/Octree-Setup-6_Initial-3_Max-Refinement-7_"   ]
setupnames = [ u"Zufällig verfeinert 2D, 65536 Zellen",
               u"Kreis 2D, 44956 Zellen"  ,
               u"Zufällig verfeinert 3D, 32768 Zellen",
               u"Sphäre 3D, 48952 Zellen"]

def powerfit(x,b,c):
     return b*x**(1./c)

for isetup in range(len(basestrings)):
    fig = figure( figsize=(6,4) );
    #title( setupnames[isetup] )
    print setupnames[isetup]

    for i in range(len(orderings)):
        fname = basestrings[isetup] + orderings[i] + "_Ordering.dat"
        if os.path.isfile( fname ):
            data = genfromtxt( fname, comments='#' )
            if orderings[i] == "Rows":
                labelname = "Spaltenweise"
            elif orderings[i] == "Morton":
                labelname = "Z-Kurve"
            elif orderings[i] == "GrayCode":
                labelname = "mod. Gray-Code"
            else:
                labelname = orderings[i]
            plot( data[20:,0], data[20:,2], 'o', ms=5, label=labelname )
            # Plot scaling
            if orderings[i] == "Rows":
                x = data[-15:-1,0]
                poptg, pcovg = curve_fit(powerfit, x, data[-15:-1,2], p0=[1e5,1.0] )
                # sigma=yerr/y * y[0]/y,
                print "Fit for Rows: ",poptg, sqrt(diag(pcovg))

                poptg[0] *= 4.
                poptg[1] = 1.0
                plot( x, powerfit(x , *poptg ) , 'k', linestyle='dashed', linewidth=2 )
                tx = 5.0
                ty = 1.8*powerfit(tx,*poptg)
                text( 0.7*tx, ty, r"$\propto n$", fontsize=18 )

            if orderings[i] == "Hilbert":
                x = data[20:-20,0]
                poptg, pcovg = curve_fit(powerfit, x, data[20:-20,2], p0=[1.0,1.0] )
                print "Fit for Hilbert: ",poptg, sqrt(diag(pcovg))

                if poptg[1] < 2.0:
                    poptg[1] = 2.0
                else:
                    poptg[1] = 3.0
                poptg[0] /= 2.
                plot( x, powerfit(x , *poptg ) , 'k-', linestyle='dashed', linewidth=2 )

                tx = 1.5*x[ len(x)/2 ]
                ty = 1.5*powerfit(tx,*poptg)
                text( tx, ty, r"$\propto n^\frac{1}{%i}$" % int(poptg[1]), fontsize=18 )

        else:
            print "Can't find ", fname
    plot( data[20:,0], data[20:,1], 'b-', lw=2, label = "Maximum" )

    xscale('log')
    yscale('log')
    ylim( ( 0, max(data[:,1])*1.1 ) )
    if ( isetup == 1 ):
        legend( loc='lower right' )
    xlabel( 'Prozesse' )
    ylabel( 'Zu sendende Daten / Byte' );
    legend(loc='lower right')
    tight_layout()
    fig.savefig( basestrings[isetup] + 'Plot.pdf', format='PDF' )

show()
