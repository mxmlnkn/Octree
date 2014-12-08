from numpy import *
from matplotlib.pyplot import *

fig = figure( figsize=(16,6) );

subplot( 1,2,1 )
title( "Circular gradient" )
data = genfromtxt( "Octree_Benchmark_circle_minRecursion4_Morton_Ordering_.dat", comments='#' )
plot( data[:,0], data[:,1], label = "Total Data accessed" )
plot( data[:,0], data[:,2], label = "Data communicated Morton" )
data = genfromtxt( "Octree_Benchmark_circle_minRecursion4_GrayCode_Ordering_.dat", comments='#' )
plot( data[:,0], data[:,2], label = "Data communicated GrayCode" )
data = genfromtxt( "Octree_Benchmark_circle_minRecursion4_Hilbert_Ordering_.dat", comments='#' )
plot( data[:,0], data[:,2], label = "Data communicated Hilbert" )

subplot( 1,2,2 )
title( "Point Refinement Setup" )
data = genfromtxt( "Octree_Benchmark_point_minRecursion4_Morton_Ordering_.dat", comments='#' )
plot( data[:,0], data[:,1], label = "Total Data accessed" )
plot( data[:,0], data[:,2], label = "Data communicated Morton" )
data = genfromtxt( "Octree_Benchmark_point_minRecursion4_GrayCode_Ordering_.dat", comments='#' )
plot( data[:,0], data[:,2], label = "Data communicated GrayCode" )
data = genfromtxt( "Octree_Benchmark_point_minRecursion4_Hilbert_Ordering_.dat", comments='#' )
plot( data[:,0], data[:,2], label = "Data communicated Hilbert" )

legend( loc='lower right' )
tight_layout()

show()
