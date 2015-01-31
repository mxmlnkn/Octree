from numpy import *
from matplotlib.pyplot import *

orderings = [ "Morton", "GrayCode", "Hilbert" ] #, "Rows"
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
basestrings = [ "./output/2015-01-31_05-12/Octree-Setup-6_Initial-8_Max-Refinement-12_"]
setupnames = [ "Circle Refinement 93184 Cells" ]

fig = figure( figsize=(12,8) );

for isetup in range(len(basestrings)):
    subplot( 2,2,isetup )
    title( setupnames[isetup] )
    for i in range(len(orderings)):
        data = genfromtxt( basestrings[isetup] + orderings[i] + ".dat", comments='#' )
        plot( data[:,0], data[:,2]/1024., '-o', label = orderings[i] )
    plot( data[:,0], data[:,1]/1024., '-o', label = "Maximum" )
    xlim( ( min(data[:,0])-0.5, max(data[:,0])+0.5 ) )
    ylim( ( 0, max(data[:,1]/1024.)*1.1 ) )
    if ( isetup == 1 ):
        legend( loc='lower right' )
    xlabel( 'Processors / #' )
    ylabel( 'Communication Data / kiB' );
tight_layout()
fig.savefig( 'Octree-SFC-Benchmarks.pdf', format='PDF' )

for isetup in range(len(basestrings)):
    fig = figure( figsize=(6,5) );
    title( setupnames[isetup] )
    for i in range(len(orderings)):
        data = genfromtxt( basestrings[isetup] + orderings[i] + ".dat", comments='#' )
        plot( data[:,0], data[:,2]/1024., '-o', label = orderings[i] )
    plot( data[:,0], data[:,1]/1024., '-o', label = "Maximum" )
    xlim( ( min(data[:,0])-0.5, max(data[:,0])+0.5 ) )
    ylim( ( 0, max(data[:,1]/1024.)*1.1 ) )
    legend( loc='best' )
    xlabel( 'Processors / #' )
    ylabel( 'Communication Data / kiB' );
    tight_layout()
    fig.savefig( 'Octree-SFC-Benchmark_'+str(isetup)+'.pdf', format='PDF' )

show()
