import os, sys
import numpy as np
import gmsh
# import gmsh_api.gmsh as gmsh

# ========================================================= #
# ===  define ports for drill magnet                    === #
# ========================================================= #

def define__ports( inpFile="dat/ports.conf" ):

    # ------------------------------------------------- #
    # --- [1] load parameters                       --- #
    # ------------------------------------------------- #
    import nkUtilities.load__table2dictarr as ltd
    params = ltd.load__table2dictarr( inpFile=inpFile )
    
    # ------------------------------------------------- #
    # --- [2] make ports                            --- #
    # ------------------------------------------------- #
    nums   = []
    for param in params:
        #  -- [2-1] object type    --  #
        if ( param["type"].lower() == "pipe" ):
            ret = gmsh.model.occ.addCylinder( 0.0, 0.0, -0.5*param["wz"], \
                                              0.0, 0.0, +1.0*param["wz"], param["r1"] )
        if ( param["type"].lower() == "cone" ):
            ret = gmsh.model.occ.addCOne    ( -0.5*param["wx"], 0.0, 0.0, \
                                              +1.0*param["wx"], 0.0, 0.0, param["r1"], param["r2"] )
        if ( param["type"].lower() == "cube" ):
            ret = gmsh.model.occ.addBox     ( -0.5*param["wx"], -0.5*param["wy"], -0.5*param["wz"], \
                                              +1.0*param["wx"], +1.0*param["wy"], +1.0*param["wz"], )
        #  -- [2-2] rotate object  --  #
        ptheta     = param["theta"] / 180.0 * np.pi
        gmsh.model.occ.rotate( [(3,ret)], 0,0,0, 0,0,1, ptheta )
        #  -- [2-3] translate      --  #
        dx, dy, dz = param["r_pos"]*np.cos( ptheta ), param["r_pos"]*np.sin( ptheta ), 0.0
        dx, dy, dz = dx+param["dx"], dy+param["dy"], dz+param["dz"]
        gmsh.model.occ.translate( [(3,ret)], dx, dy, dz )
        #  -- [2-4] append number  --  #
        nums.append( ret )
    return( nums )



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    gmsh.initialize()
    gmsh.option.setNumber( "General.Terminal", 1 )
    gmsh.model.add( "model" )

    portFile = "dat/ports.conf"
    define__ports( inpFile=portFile )
    
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber( "Mesh.CharacteristicLengthMin", 0.1 )
    gmsh.option.setNumber( "Mesh.CharacteristicLengthMax", 0.1 )
    gmsh.model.mesh.generate(3)
    gmsh.write( "msh/ports.msh" )
    gmsh.finalize()
    
