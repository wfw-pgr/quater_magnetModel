import numpy as np
import os, sys
import gmsh

# ------------------------------------------------- #
# --- [1] initialization of the gmsh            --- #
# ------------------------------------------------- #
gmsh.initialize()
gmsh.option.setNumber( "General.Terminal", 1 )
gmsh.model.add( "model" )

# ------------------------------------------------- #
# --- [2] initialize settings                   --- #
# ------------------------------------------------- #
ptsDim , lineDim , surfDim , voluDim  =  0,  1,  2,  3
pts    , line    , surf    , volu     = {}, {}, {}, {}
ptsPhys, linePhys, surfPhys, voluPhys = {}, {}, {}, {}
lc                                    = 0.1
x_, y_, z_, lc_, tag_                 = 0, 1, 2, 3, 4


# ------------------------------------------------- #
# --- [3] Modeling                              --- #
# ------------------------------------------------- #

lc   = 0.1
side = "+"
import generate__magnet as mag
mag.generate__magnet( lc=lc, side=side )

# ------------------------------------------------- #
# --- [4] Physical Grouping                     --- #
# ------------------------------------------------- #
gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

# ------------------------------------------------- #
# --- [5] Mesh settings                         --- #
# ------------------------------------------------- #
# dimtags1      = [ (3,1) ]
# dimtags2      = [ (3,2) ]
# dimtags3      = [ (3,3) ]
# element_size1 = 10
# element_size2 = 20
# element_size3 = 5
# points1  = gmsh.model.getBoundary( dimtags1, recursive=True )
# points2  = gmsh.model.getBoundary( dimtags2, recursive=True )
# points3  = gmsh.model.getBoundary( dimtags3, recursive=True )
# gmsh.model.occ.synchronize()
# gmsh.model.occ.mesh.setSize( points1, element_size1 )
# gmsh.model.mesh.setSize( points2, element_size2 )
# gmsh.model.mesh.setSize( points3, element_size3 )


# ------------------------------------------------- #
# --- [2] post process                          --- #
# ------------------------------------------------- #
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write( "msh/model.geo_unrolled" )
gmsh.write( "msh/model.msh" )
gmsh.finalize()

