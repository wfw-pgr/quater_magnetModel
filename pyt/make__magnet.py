import os, sys
import gmsh
import numpy   as np

# ========================================================= #
# ===  make__magnet routine                             === #
# ========================================================= #

def make__magnet():

    # ------------------------------------------------- #
    # --- [1] load config                           --- #
    # ------------------------------------------------- #
    cnsFile = "dat/parameter.conf"
    import nkUtilities.load__constants as lcn
    const   = lcn.load__constants( inpFile=cnsFile )
    side    = const["geometry.side"]

    
    # ------------------------------------------------- #
    # --- [2] initialization of the gmsh            --- #
    # ------------------------------------------------- #
    gmsh.initialize()
    gmsh.option.setNumber( "General.Terminal", 1 )
    gmsh.option.setNumber( "Mesh.Algorithm"  , const["mesh.algorithm2D"] )
    gmsh.option.setNumber( "Mesh.Algorithm3D", const["mesh.algorithm3D"] )
    gmsh.option.setNumber( "Mesh.SubdivisionAlgorithm", const["mesh.subdivision"] )
    gmsh.model.add( "model" )

    
    # ------------------------------------------------- #
    # --- [3] Modeling                              --- #
    # ------------------------------------------------- #
    if ( const["geometry.import_model"] ):
        stpFile = "msh/model.step"
        gmsh.model.occ.importShapes( stpFile )
        const["geometry.save_step"] = False
    else:
        import generate__magnetParts as mag
        mag.generate__magnetParts( side=side )
        
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    
    # ------------------------------------------------- #
    # --- [4] save model                            --- #
    # ------------------------------------------------- #
    if ( const["geometry.save_step"] ):
        gmsh.write( "msh/model.step" )

        
    # ------------------------------------------------- #
    # --- [5] Mesh settings                         --- #
    # ------------------------------------------------- #
    meshFile = "dat/mesh.conf"
    if   ( side == "+" ):
        physFile = "dat/phys_right.conf"
    elif ( side == "-" ):
        physFile = "dat/phys_left.conf"
    elif ( side in ["+-","-+"] ):
        physFile = "dat/phys_both.conf"
    else:
        sys.exit( "[make__magnet.py] side == {0} ??? ".format( side ) )

    if ( const["mesh.uniform"] ):
        gmsh.option.setNumber( "Mesh.CharacteristicLengthMin", 0.2 )
        gmsh.option.setNumber( "Mesh.CharacteristicLengthMax", 0.2 )
    else:
        import nkGmshRoutines.assign__meshsize as ams
        meshes = ams.assign__meshsize( meshFile=meshFile, physFile=physFile )

    
    # ------------------------------------------------- #
    # --- [6] Meshing / save mesh                   --- #
    # ------------------------------------------------- #
    #  -- [6-1] meshing                             --  #
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    
    #  -- [6-2] save mesh                           --  #
    gmsh.write( "msh/model.msh" )
    if ( const["mesh.save_bdf"] ):
        gmsh.write( "msh/model.bdf"  )


    # ------------------------------------------------- #
    # --- [7] post-process                          --- #
    # ------------------------------------------------- #
    gmsh.finalize()
    return()



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #
if ( __name__=="__main__" ):
    make__magnet()
