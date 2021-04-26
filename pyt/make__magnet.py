import numpy as np
import os, sys
import gmsh

# ------------------------------------------------- #
# --- [1] initialization of the gmsh            --- #
# ------------------------------------------------- #
gmsh.initialize()
gmsh.option.setNumber( "General.Terminal", 1 )
gmsh.option.setNumber( "Mesh.Algorithm"  , 1 )
gmsh.option.setNumber( "Mesh.Algorithm3D", 1 )
gmsh.option.setNumber( "Mesh.SubdivisionAlgorithm", 1 )
gmsh.model.add( "model" )


# ------------------------------------------------- #
# --- [2] Modeling                              --- #
# ------------------------------------------------- #
make_model = True
if ( make_model ):
    side = "+"
    import generate__magnetParts as mag
    mag.generate__magnetParts( side=side )
else:
    stpFile = "msh/model.step"
    gmsh.model.occ.importShapes( stpFile )
gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()


# ------------------------------------------------- #
# --- [3] Mesh settings                         --- #
# ------------------------------------------------- #
meshFile = "dat/mesh.conf"
import nkGmshRoutines.assign__meshsize as ams
meshes = ams.assign__meshsize( meshFile=meshFile )
gmsh.option.setNumber( "Mesh.CharacteristicLengthMin", np.min( meshes["meshsize_list"] ) )
gmsh.option.setNumber( "Mesh.CharacteristicLengthMax", np.max( meshes["meshsize_list"] ) )


# ------------------------------------------------- #
# --- [4] post process                          --- #
# ------------------------------------------------- #
if ( make_model ):
    gmsh.write( "msh/model.step" )
    gmsh.write( "msh/model.geo_unrolled" )
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write( "msh/model.msh" )
gmsh.finalize()

