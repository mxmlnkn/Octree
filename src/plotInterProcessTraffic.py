from numpy import *
from matplotlib.pyplot import *
import os.path

orderings = [ "Morton", "GrayCode", "Hilbert", "Rows" ]
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
basestrings = [ "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Quadtree-Setup-5_Initial-8_Max-Refinement-6_" ,
                "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Quadtree-Setup-6_Initial-7_Max-Refinement-12_",
                "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Octree-Setup-5_Initial-5_Max-Refinement-6_"   ,
                "./output/2015-02-14_01-44-53_Octree_Benchmark_3D_and_2D_multiple_Setups/Octree-Setup-6_Initial-5_Max-Refinement-10_"  ]
setupnames = [ u"Zufällig verfeinert 2D, 65536 Zellen",
               u"Kreis 2D, 44956 Zellen"  ,
               u"Zufällig verfeinert 3D, 32768 Zellen",
               u"Sphäre 3D, 49428 Zellen"  ]

for isetup in range(len(basestrings)):
    fig = figure( figsize=(8,5.5) );
    title( setupnames[isetup] )
    for i in range(len(orderings)):
        fname = basestrings[isetup] + orderings[i] + "_Ordering.dat"
        if os.path.isfile( fname ):
            data = genfromtxt( fname, comments='#' )
            plot( data[20:,0], data[20:,2]/1000., 'o', label = orderings[i] )
        else:
            print "Can't find ", fname
    plot( data[20:,0], data[20:,1]/1000., 'b-', lw=2, label = "Maximum" )
    xscale('log')
    yscale('log')
    ylim( ( 0, max(data[:,1]/1000.)*1.1 ) )
    if ( isetup == 1 ):
        legend( loc='lower right' )
    xlabel( 'Prozesse' )
    ylabel( 'Zu sendende Daten / kB' );
    legend(loc='best')
    tight_layout()
    fig.savefig( basestrings[isetup] + 'Plot.pdf', format='PDF' )

show()
