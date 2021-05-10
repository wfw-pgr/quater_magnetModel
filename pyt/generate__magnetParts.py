import os, sys, subprocess
import gmsh
import numpy                              as np
import nkGmshRoutines.generate__coneShape as con
import nkGmshRoutines.generate__sector180 as sec


# ========================================================= #
# === generate magnet shape                             === #
# ========================================================= #

def generate__magnetParts( side="+" ):

    cnsFile    = "dat/parameter.conf"
    import nkUtilities.load__constants as lcn
    const      = lcn.load__constants( inpFile=cnsFile )

    hexahedral = const["mesh.recombine"]
    
    r_pole     = const["geometry.r_pole"]
    w_iair1    = const["geometry.w_iair1"]
    w_iair2    = const["geometry.w_iair2"]
    w_coil     = const["geometry.w_coil"]
    w_yoke     = const["geometry.w_yoke"]
    w_cut      = const["geometry.w_cut"]
    w_oair     = const["geometry.w_oair"]

    z_gap      = const["geometry.z_gap"]
    z_pole     = const["geometry.z_pole"]
    z_root     = const["geometry.z_root"]
    
    h_iair1    = const["geometry.h_iair1"]
    h_iair2    = const["geometry.h_iair2"]
    h_coil     = const["geometry.h_coil"]
    h_yoke     = const["geometry.h_yoke"]
    h_cut      = const["geometry.h_cut"]
    h_oair     = const["geometry.h_oair"]

    h_slot     = h_iair1 + h_coil + h_iair2
        
    
    # ------------------------------------------------- #
    # --- [1] pole making                           --- #
    # ------------------------------------------------- #
    if ( const["geometry.flat_pole"] ):
        if   ( side == "+" ):
            generate__pole( r1=0.0, r2=r_pole, z1=0.0, z2=z_gap, z3=z_pole, z4=z_root, \
                            side="+", hexahedral=hexahedral )
            
        elif ( side == "-" ):
            generate__pole( r1=0.0, r2=r_pole, z1=0.0, z2=z_gap, z3=z_pole, z4=z_root, \
                            side="-", hexahedral=hexahedral )

        elif ( side in ["+-","-+"] ):
            sys.exit( "[generate__magnetParts] side == both for flat_pole is not coded, now." )
            
    else:
        if ( z_root != h_slot ):
            sys.exit( "[generate__magnetParts] incompatible slot-depth and pole-root-length")

        import generate__poleLayer as gpl
        if   ( side == "+" ):
            gpl.generate__poleLayer( side="+", z1=z_pole, z2=z_root, radius=r_pole )

        elif ( side == "-" ):
            gpl.generate__poleLayer( side="-", z1=z_pole, z2=z_root, radius=r_pole )

        elif ( side in ["+-","-+"] ):
            # -- generate (+) side and save it -- #
            gmsh.model.add( "pole_right" )
            gpl.generate__poleLayer( side="+" , z1=z_pole, z2=z_root, radius=r_pole )
            gmsh.write( "msh/pole_right.step" )
            gmsh.model.remove()
            # -- generate (+) side and save it -- #
            gmsh.model.add( "pole_left" )
            gpl.generate__poleLayer( side="-" , z1=z_pole, z2=z_root, radius=r_pole )
            gmsh.write( "msh/pole_left.step" )
            gmsh.model.remove()
            # -- load each model again         -- #
            gmsh.model.setCurrent( "model" )
            gmsh.model.occ.importShapes( "msh/pole_right.step" )
            gmsh.model.occ.importShapes( "msh/pole_left.step" )
            gmsh.model.occ.synchronize()
            gmsh.model.occ.removeAllDuplicates()
            gmsh.model.occ.synchronize()
            
        
    # ------------------------------------------------- #
    # --- [2] coil making                           --- #
    # ------------------------------------------------- #
    r1 = r_pole
    r2 = r_pole + w_iair1
    r3 = r_pole + w_iair1 + w_coil
    r4 = r_pole + w_iair1 + w_coil + w_iair2
    z1 = 0.0
    z2 = h_iair1
    z3 = h_iair1 + h_coil
    z4 = h_iair1 + h_coil + h_iair2
    if ( side in ["+","+-","-+"] ):
        generate__coilslot( r1=r1, r2=r2, r3=r3, r4=r4, \
                            z1=z1, z2=z2, z3=z3, z4=z4, side="+", hexahedral=hexahedral )
    if ( side in ["-","+-","-+"] ):
        generate__coilslot( r1=r1, r2=r2, r3=r3, r4=r4, \
                            z1=z1, z2=z2, z3=z3, z4=z4, side="-", hexahedral=hexahedral )

    # ------------------------------------------------- #
    # --- [3]  yoke making                          --- #
    # ------------------------------------------------- #
    r1 = 0.0
    r2 = r_pole + w_iair1 + w_coil + w_iair2
    r3 = r2 + w_yoke-w_cut
    r4 = r2 + w_yoke
    z1 = 0.0
    z2 = h_slot
    z3 = h_slot + h_yoke - h_cut
    z4 = h_slot + h_yoke
    if ( side in ["+","+-","-+"] ):
        generate__yoke    ( r1=r1, r2=r2, r3=r3, r4=r4, \
                            z1=z1, z2=z2, z3=z3, z4=z4, side="+", hexahedral=hexahedral )
    if ( side in ["-","+-","-+"] ):
        generate__yoke    ( r1=r1, r2=r2, r3=r3, r4=r4, \
                            z1=z1, z2=z2, z3=z3, z4=z4, side="-", hexahedral=hexahedral )
    
    # ------------------------------------------------- #
    # --- [4]  outside Air making                   --- #
    # ------------------------------------------------- #
    r1 = 0.0
    r3 = r_pole + w_iair1 + w_coil + w_iair2 + w_yoke
    r2 = r3 - w_cut
    r4 = r3 + w_oair
    z1 = 0.0
    z3 = h_slot + h_yoke
    z2 = h_slot + h_yoke - h_cut
    z4 = h_slot + h_yoke + h_oair
    if ( side in ["+","+-","-+"] ):
        generate__outAir  ( r1=r1, r2=r2, r3=r3, r4=r4, \
                            z1=z1, z2=z2, z3=z3, z4=z4, side="+", hexahedral=hexahedral )
    if ( side in ["-","+-","-+"] ):
        generate__outAir  ( r1=r1, r2=r2, r3=r3, r4=r4, \
                            z1=z1, z2=z2, z3=z3, z4=z4, side="-", hexahedral=hexahedral )
    

    
# ========================================================= #
# ===  generate pole parts                              === #
# ========================================================= #
def generate__pole( lc=0.0, r1=0.0, r2=0.7, \
                    z1=0.0, z2=0.2, z3=0.7, z4=1.0, side="+", hexahedral=False ):
    # ------------------------------------------------- #
    # --- [1] generate pole parts                   --- #
    # ------------------------------------------------- #
    origin  = [ 0.0, 0.0 ]
    gap     = sec.generate__sector180( lc=lc, r1=r1, r2=r2, hexahedral=hexahedral, tag=-1, \
                                       fuse=True, \
                                       zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    pole    = sec.generate__sector180( lc=lc, r1=r1, r2=r2, hexahedral=hexahedral, tag=-1, \
                                       fuse=True, 
                                       zoffset=z2, height=z3-z2, defineVolu=True, side=side )
    body    = sec.generate__sector180( lc=lc, r1=r1, r2=r2, hexahedral=hexahedral, tag=-1, \
                                       fuse=True, 
                                       zoffset=z3, height=z4-z3, defineVolu=True, side=side )
    return()



# ========================================================= #
# ===  generate coil & inner air parts                  === #
# ========================================================= #
def generate__coilslot( lc=0.0, r1=0.8, r2=0.9, r3=1.0, r4=1.1, \
                        z1=0.0, z2=0.2, z3=0.4, z4=0.7, side="+", hexahedral=False ):
    # ------------------------------------------------- #
    # --- [1] generate coil parts                   --- #
    # ------------------------------------------------- #
    origin  = [ 0.0, 0.0 ]
    air_inn = sec.generate__sector180( lc=lc, r1=r1, r2=r2, zoffset=z1, height=z4-z1, \
                                       defineVolu=True, side=side, hexahedral=hexahedral )
    air_bot = sec.generate__sector180( lc=lc, r1=r2, r2=r3, zoffset=z1, height=z2-z1, \
                                       defineVolu=True, side=side, hexahedral=hexahedral )
    coil    = sec.generate__sector180( lc=lc, r1=r2, r2=r3, zoffset=z2, height=z3-z2, \
                                       defineVolu=True, side=side, hexahedral=hexahedral )
    air_top = sec.generate__sector180( lc=lc, r1=r2, r2=r3, zoffset=z3, height=z4-z3, \
                                       defineVolu=True, side=side, hexahedral=hexahedral )
    air_out = sec.generate__sector180( lc=lc, r1=r3, r2=r4, zoffset=z1, height=z4-z1, \
                                       defineVolu=True, side=side, hexahedral=hexahedral )
    return()


# ========================================================= #
# ===  generate york parts                              === #
# ========================================================= #
def generate__yoke( lc=0.0, r1=0.0, r2=1.1, r3=1.4, r4=1.5, \
                    z1=0.0, z2=0.7, z3=1.1, z4=1.2, side="+", hexahedral=False ):
    # ------------------------------------------------- #
    # --- [1] generate coil parts                   --- #
    # ------------------------------------------------- #
    origin      = [ 0.0, 0.0 ]
    origin_cone = [0.0,0.0,z3]
    th1,th2     = -90.0, 90.0
    yoke_h      = sec.generate__sector180( lc=lc, r1=r1, r2=r4, zoffset=z2, height=z3-z2, \
                                           defineVolu=True, side=side, hexahedral=hexahedral )
    yoke_v      = sec.generate__sector180( lc=lc, r1=r2, r2=r4, zoffset=z1, height=z2-z1, \
                                           defineVolu=True, side=side, hexahedral=hexahedral )
    york_c      = con.generate__coneShape( lc =lc , origin=origin_cone, r1=r4, r2=r3, \
                                           th1=th1, th2=th2, height=z4-z3, side=side )
    return()


# ========================================================= #
# ===  generate outside Air region                      === #
# ========================================================= #
def generate__outAir( lc=0.0, r1=0.0, r2=1.1, r3=1.4, r4=1.5, \
                      z1=0.0, z2=0.7, z3=1.1, z4=1.2, side="+", hexahedral=False ):
    # ------------------------------------------------- #
    # --- [1] main outside air cylinder             --- #
    # ------------------------------------------------- #
    th1,th2     = -90.0, 90.0
    origin      = [ 0.0, 0.0 ]
    origin_cone = [0.0,0.0,z2]
    oAir_h      = sec.generate__sector180( lc=lc, r1=r1, r2=r4, zoffset=z2, height=z4-z2, \
                                           defineVolu=True, side=side, hexahedral=hexahedral )
    oAir_v      = sec.generate__sector180( lc=lc, r1=r3, r2=r4, zoffset=z1, height=z2-z1, \
                                           defineVolu=True, side=side, hexahedral=hexahedral )
    oAir_c      = con.generate__coneShape( lc =lc , origin=origin_cone, r1=r3, r2=r2, \
                                           th1=th1, th2=th2, height=z3-z2, side=side )
    # ------------------------------------------------- #
    # --- [2] cut cone shape from cylinder          --- #
    # ------------------------------------------------- #
    target      = [ (3,( (oAir_h[0] )["volu"])["sector"]),(3,( (oAir_h[1] )["volu"])["sector"]) ]
    tool        = [(3,(oAir_c["volu"])["cone"])]
    ret         = gmsh.model.occ.cut( target, tool )
    return()
