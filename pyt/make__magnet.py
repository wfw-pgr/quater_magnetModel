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
    # --- [4] define port                           --- #
    # ------------------------------------------------- #
    if ( const["geometry.add_port"] ):
        # -- [4-1] save wo port model               --  #
        gmsh.write( "msh/model_woport.step" )
        # -- [4-2] define ports                     --  #
        import nkGmshRoutines.define__ports as dfp
        inpFile      = "dat/ports.conf"
        portNums     = dfp.define__ports( inpFile=inpFile )
        gmsh.model.occ.synchronize()
        # -- [4-3] boolean cut from yoke            --  #
        tools        = [ (3,tool  ) for tool   in portNums                 ]
        targets      = [ (3,target) for target in const["geometry.yoke_tobecut"] ]
        copy         = gmsh.model.occ.copy( targets )
        yoke_p       = gmsh.model.occ.cut ( targets, tools, removeObject=True, removeTool=False )
        holes        = gmsh.model.occ.intersect( tools, copy, removeObject=True, removeTool=False )
        gmsh.model.occ.synchronize()
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        gmsh.write( "msh/model.geo_unrolled" )

    
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

    if ( const["mesh.compound"] ):
        surfDim,voluDim = 2, 3
        physNum_gap     = 301
        physNum_pole    = 302
        volu_gap        = gmsh.model.getEntitiesForPhysicalGroup( 3, physNum_gap  )
        volu_pole       = gmsh.model.getEntitiesForPhysicalGroup( 3, physNum_pole )
        dimtag_gap      = [ (voluDim,vnum) for vnum in volu_gap  ]
        dimtag_pole     = [ (voluDim,vnum) for vnum in volu_pole ]
        surf_gap        = gmsh.model.getBoundary( dimtag_gap  )
        surf_pole       = gmsh.model.getBoundary( dimtag_pole )
        surf_gap        = [ dimtag[1] for dimtag in surf_gap  ]
        surf_pole       = [ dimtag[1] for dimtag in surf_pole ]
        surf_common     = list( set( surf_gap ) & set( surf_pole ) )
        gmsh.model.mesh.setCompound( surfDim, surf_common )

        
    # ------------------------------------------------- #
    # --- [6] Meshing / save mesh                   --- #
    # ------------------------------------------------- #
    #  -- [6-1] meshing                             --  #
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)

    #  -- [6-2] optimization                        --  #
    if ( const["mesh.optimize"] ):
        gmsh.option.setNumber( "Mesh.OptimizeThreshold", const["mesh.opt_threshold"] )
        gmsh.model.mesh.optimize( "Netgen" )
        gmsh.model.mesh.optimize( "Relocate3D" )
    
    #  -- [6-3] save mesh                           --  #
    gmsh.option.setNumber( "Mesh.SaveElementTagType", 2 )
    gmsh.option.setNumber( "Mesh.BdfFieldFormat"    , 0 )
    if ( const["mesh.save_bdf"] ):
        gmsh.write( "msh/model.bdf" )
    if ( const["mesh.save_msh"] ):
        gmsh.write( "msh/model.msh" )


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
