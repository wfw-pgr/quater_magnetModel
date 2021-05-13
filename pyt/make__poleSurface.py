import numpy   as np
import os, sys, subprocess
import gmsh


# ========================================================= #
# ===  interpolation onto mesh                          === #
# ========================================================= #
def make__poleSurface():

    # ------------------------------------------------- #
    # --- [1] load constants                        --- #
    # ------------------------------------------------- #
    cnsFile = "dat/parameter.conf"
    import nkUtilities.load__constants as lcn
    const   = lcn.load__constants( inpFile=cnsFile )
    lc1     = const["geometry.pole_lc_top"]
    lc2     = const["geometry.pole_lc_bot"]
    radius  = const["geometry.r_pole"]
    side    = const["geometry.side"]

    # ------------------------------------------------- #
    # --- [2] interpolation / gmsh <-> elmer        --- #
    # ------------------------------------------------- #
    nop     = True
    
    if ( side in ["+","+-","-+"] ):
        generate__mesh_to_interpolate( lc1=lc1, lc2=lc2, radius=radius, side="+" )
        convert__meshFormat( direction="gmsh->elmer" )
        ret     = interpolate__grid_to_mesh( side="+" )
        nop     = False

    if ( side in ["-","+-","-+"] ):
        generate__mesh_to_interpolate( lc1=lc1, lc2=lc2, radius=radius, side="-" )
        convert__meshFormat( direction="gmsh->elmer" )
        ret     = interpolate__grid_to_mesh( side="-" )
        nop     = False

    if ( side in ["+-","-+"]     ):
        generate__mesh_to_interpolate( lc1=lc1, lc2=lc2, radius=radius, side="+-" )
        convert__meshFormat( direction="gmsh->elmer" )
        ret     = interpolate__grid_to_mesh( side="+-" )
        nop     = False

    if ( nop ):
        print( "[make__poleSurface.py] no Operation !!! ERROR !!! " )
    else:
        convert__meshFormat( direction="elmer->gmsh", side=side )

    return()


# ========================================================= #
# ===   generate__mesh_to_interpolate                   === #
# ========================================================= #

def generate__mesh_to_interpolate( lc1=0.0, lc2=0.0, radius=1.0, side="+" ):

    # ------------------------------------------------- #
    # --- [1] initialization of the gmsh            --- #
    # ------------------------------------------------- #
    gmsh.initialize()
    gmsh.option.setNumber( "General.Terminal", 1 )
    gmsh.model.add( "model" )

    # ------------------------------------------------- #
    # --- [2] Modeling                              --- #
    # ------------------------------------------------- #
    r1, r2  = 0.0, radius
    origin  = [ 0.0, 0.0 ]
    import nkGmshRoutines.generate__sector90 as s90
    if   ( side == "+" ):
        ret1    = s90.generate__sector90( lc=lc1, r1=r1, r2=r2, quadrant=1, \
                                          origin=origin, defineSurf=True )
        ret2    = s90.generate__sector90( lc=lc2, r1=r1, r2=r2, quadrant=4, \
                                          origin=origin, defineSurf=True )
    elif ( side == "-" ):
        ret1    = s90.generate__sector90( lc=lc1, r1=r1, r2=r2, quadrant=2, \
                                          origin=origin, defineSurf=True )
        ret2    = s90.generate__sector90( lc=lc2, r1=r1, r2=r2, quadrant=3, \
                                          origin=origin, defineSurf=True )

    elif ( side in ["+-","-+"] ):
        ret1    = s90.generate__sector90( lc=lc1, r1=r1, r2=r2, quadrant=1, \
                                          origin=origin, defineSurf=True )
        ret2    = s90.generate__sector90( lc=lc2, r1=r1, r2=r2, quadrant=2, \
                                          origin=origin, defineSurf=True )
        ret3    = s90.generate__sector90( lc=lc2, r1=r1, r2=r2, quadrant=3, \
                                          origin=origin, defineSurf=True )
        ret4    = s90.generate__sector90( lc=lc2, r1=r1, r2=r2, quadrant=4, \
                                          origin=origin, defineSurf=True )
        
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    
    # ------------------------------------------------- #
    # --- [3] Mesh generation & save                --- #
    # ------------------------------------------------- #
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write( "msh/mesh2d.msh" )
    gmsh.finalize()


    
# ========================================================= #
# === convert mesh into easy-to-read format using Elmer === #
# ========================================================= #
def convert__meshFormat( direction=None, side="+" ):

    if   ( side == "+" ):
        inpFile = "dat/onmesh_right.dat"
    elif ( side == "-" ):
        inpFile = "dat/onmesh_left.dat"
    elif ( side in ["+-","-+"] ):
        inpFile = "dat/onmesh_both.dat"
    else:
        print( "[convert__meshFormat] side == {0} ??? ".format( side ) )
        sys.exit()
        
    if ( direction is None ):
        sys.exit( "[convert__meshFormat] direction ??? ( gmsh->elmer or elmer->gmsh ) " )
    if ( not( direction in [ "elmer->gmsh", "gmsh->elmer" ] ) ):
        sys.exit( "[convert__meshFormat] direction ??? ( gmsh->elmer or elmer->gmsh ) " )
    
    # ------------------------------------------------- #
    # --- [1] gmsh -> elmer mode                    --- #
    # ------------------------------------------------- #
    if   ( direction.lower() == "gmsh->elmer"  ):

        # -- use ElmerGrid command to convert -- #
        print()
        print( "[convert__meshFormat] convert .msh (Gmsh) File into Elmer-Format.... " )
        cmd = "ElmerGrid 14 2 msh/mesh2d.msh"
        print( cmd )
        subprocess.call( cmd.split() )

    # ------------------------------------------------- #
    # --- [2] elmer -> gmsh mode                    --- #
    # ------------------------------------------------- #
    if   ( direction.lower() == "elmer->gmsh" ):
        
        # -- copy original Elmer-Format       -- #
        print()
        print( "[reconvert__meshFormat] re-convert Elmer-Format into Gmsh.... " )
        cmd = "mkdir -p msh/mesh3d"
        print( cmd )
        subprocess.call( cmd.split() )
        cmd = "cp -r msh/mesh2d/* msh/mesh3d/"
        print( cmd )
        subprocess.call( cmd, shell=True )
        # -- modify nodes position            -- #
        import nkUtilities.load__pointFile as lpf
        nodFile     = "msh/mesh3d/mesh.nodes"
        nodes       = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
        mesh        = lpf.load__pointFile( inpFile=nodFile, returnType="point" )
        mesh[:,2:5] = nodes[:,:]
        # -- rewrite back into nodFile        -- #
        fmt         = [ "%d", "%d", "%15.8e", "%15.8e", "%15.8e" ]
        with open( nodFile, "w" ) as f:
            np.savetxt( f, mesh, fmt=fmt )
        # -- reconvert into gmsh Format       -- #
        os.chdir( "msh/" )
        cmd = "ElmerGrid 2 4 mesh3d"
        print( cmd )
        subprocess.call( cmd.split() )
        os.chdir( "../" )
    return()


# ========================================================= #
# ===  interpolate__grid_to_mesh                        === #
# ========================================================= #
def interpolate__grid_to_mesh( gridFile="dat/mshape_svd.dat", meshFile="msh/mesh2d/mesh.nodes", \
                               side="+", interpolation="cubic" ):
    
    if ( side == "+" ):
        omsFile = "dat/onmesh_right.dat"
        elmFile = "dat/mesh_right.elements"
    if ( side == "-" ):
        omsFile = "dat/onmesh_left.dat"
        elmFile = "dat/mesh_left.elements"
    if ( side in ["+-","-+"] ):
        omsFile = "dat/onmesh_both.dat"
        elmFile = "dat/mesh_both.elements"
    
    # ------------------------------------------------- #
    # --- [1] load grid & mesh Data                 --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    grid      = lpf.load__pointFile( inpFile=gridFile, returnType="structured" )
    mesh      = lpf.load__pointFile( inpFile=meshFile, returnType="point"      )
    gridData  = grid[0,:,:,0:3]
    meshData  = mesh[:,2:5]

    # ------------------------------------------------- #
    # --- [2] interpolation 2D                      --- #
    # ------------------------------------------------- #
    if   ( interpolation == "linear" ):
        import nkInterpolator.interpolate__linear2D as li2
        ret       = li2.interpolate__linear2D( gridData=gridData, pointData=meshData )
    elif ( interpolation == "cubic"  ):
        import nkInterpolator.interpolate__bicubic  as bic
        ret       = bic.interpolate__bicubic( gridData=gridData, pointData=meshData )
        
    # ------------------------------------------------- #
    # --- [3] save in a File                        --- #
    # ------------------------------------------------- #
    #  -- [3-1] saving data                         --  #
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=omsFile, Data=ret )
    #  -- [3-2] copy mesh.elements                  --  #
    cmd       = "cp msh/mesh2d/mesh.elements {0}".format( elmFile )
    print( "\n" + "[make__poleSurface] copy mesh.elements... " )
    print( cmd + "\n" )
    subprocess.call( cmd, shell=True )
    
    #  -- [3-3] save figure                         --  #
    import nkUtilities.cMapTri as cmt
    pngFile = "png/onmesh.png"
    cmt.cMapTri( xAxis=ret[:,0], yAxis=ret[:,1], cMap=ret[:,2], pngFile=pngFile  )
    
    return( ret )


# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #

if ( __name__=="__main__" ):
    make__poleSurface()
