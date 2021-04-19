import os, sys, subprocess
import gmsh
import numpy                              as np
import nkGmshRoutines.generate__coneShape as con
import nkGmshRoutines.generate__sector180 as sec


# ========================================================= #
# === generate magnet shape                             === #
# ========================================================= #

def generate__magnetParts( side="+" ):

    cnsFile = "dat/parameter.conf"
    import nkUtilities.load__constants as lcn
    const = lcn.load__constants( inpFile=cnsFile )

    # ------------------------------------------------- #
    # --- [1] pole making                           --- #
    # ------------------------------------------------- #
    generate__pole( lc=const["lc_pole"], r1=0.0  , r2=const["r_pole"], \
                    z1=0.0, z2=const["z_gap"], z3=const["z_pole"], side=side )

    # ------------------------------------------------- #
    # --- [2] coil making                           --- #
    # ------------------------------------------------- #
    r1 = const["r_pole"]
    r2 = const["r_pole"]+const["w_iair1"]
    r3 = const["r_pole"]+const["w_iair1"]+const["w_coil"]
    r4 = const["r_pole"]+const["w_iair1"]+const["w_coil"]+const["w_iair2"]
    z1 = 0.0
    z2 = const["h_iair1"]
    z3 = const["h_iair1"]+const["h_coil"]
    z4 = const["h_iair1"]+const["h_coil"]+const["h_iair2"]
    generate__coilslot( lc=const["lc_coil"], r1=r1, r2=r2, r3=r3, r4=r4, \
                        z1=z1, z2=z2, z3=z3, z4=z4, side=side )

    # ------------------------------------------------- #
    # --- [3]  yoke making                          --- #
    # ------------------------------------------------- #
    r1 = 0.0
    r2 = const["r_pole"]+const["w_iair1"]+const["w_coil"]+const["w_iair2"]
    r3 = r2+const["w_yoke"]-const["w_cut"]
    r4 = r2+const["w_yoke"]
    z1 = 0.0
    z2 = const["h_iair1"]+const["h_coil"]+const["h_iair2"]
    z3 = z2+const["h_yoke"]-const["h_cut"]
    z4 = z2+const["h_yoke"]
    generate__yoke    ( lc=const["lc_yoke"], r1=r1, r2=r2, r3=r3, r4=r4, \
                        z1=z1, z2=z2, z3=z3, z4=z4, side=side )
    
    # ------------------------------------------------- #
    # --- [4]  outside Air making                   --- #
    # ------------------------------------------------- #
    r1 = 0.0
    r3 = const["r_pole"]+const["w_iair1"]+const["w_coil"]+const["w_iair2"]+const["w_yoke"]
    r2 = r3-const["w_cut"]
    r4 = r3+const["w_oair"]
    z1 = 0.0
    z3 = const["h_iair1"]+const["h_coil"]+const["h_iair2"]+const["h_yoke"]
    z2 = z3-const["h_cut"]
    z4 = z3+const["h_oair"]
    generate__outAir  ( lc=const["lc_oair"], r1=r1, r2=r2, r3=r3, r4=r4, \
                        z1=z1, z2=z2, z3=z3, z4=z4, side=side )
    

    
# ========================================================= #
# ===  generate pole parts                              === #
# ========================================================= #
def generate__pole( lc=0.2, r1=0.0, r2=0.7, \
                    z1=0.0, z2=0.2, z3=0.7, side="+" ):
    # ------------------------------------------------- #
    # --- [1] generate pole parts                   --- #
    # ------------------------------------------------- #
    origin  = [ 0.0, 0.0 ]
    gap     = sec.generate__sector180( lc=lc, r1=r1, r2=r2, \
                                       zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    pole    = sec.generate__sector180( lc=lc, r1=r1, r2=r2, \
                                       zoffset=z2, height=z3-z2, defineVolu=True, side=side )
    return()



# ========================================================= #
# ===  generate coil & inner air parts                  === #
# ========================================================= #
def generate__coilslot( lc=0.2, r1=0.8, r2=0.9, r3=1.0, r4=1.1, \
                        z1=0.0, z2=0.2, z3=0.4, z4=0.7, side="+" ):
    # ------------------------------------------------- #
    # --- [1] generate coil parts                   --- #
    # ------------------------------------------------- #
    print( r1, r2, r3, r4 )
    print( z1, z2, z3, z4 )
    origin  = [ 0.0, 0.0 ]
    air_inn = sec.generate__sector180( lc=lc, r1=r1, r2=r2, \
                                       zoffset=z1, height=z4-z1, defineVolu=True, side=side )
    air_bot = sec.generate__sector180( lc=lc, r1=r2, r2=r3, \
                                       zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    coil    = sec.generate__sector180( lc=lc, r1=r2, r2=r3, \
                                       zoffset=z2, height=z3-z2, defineVolu=True, side=side )
    air_top = sec.generate__sector180( lc=lc, r1=r2, r2=r3, \
                                       zoffset=z3, height=z4-z3, defineVolu=True, side=side )
    air_out = sec.generate__sector180( lc=lc, r1=r3, r2=r4, \
                                       zoffset=z1, height=z4-z1, defineVolu=True, side=side )
    return()


# ========================================================= #
# ===  generate york parts                              === #
# ========================================================= #
def generate__yoke( lc=0.2, r1=0.0, r2=1.1, r3=1.4, r4=1.5, \
                    z1=0.0, z2=0.7, z3=1.1, z4=1.2, side="+" ):
    # ------------------------------------------------- #
    # --- [1] generate coil parts                   --- #
    # ------------------------------------------------- #
    origin      = [ 0.0, 0.0 ]
    origin_cone = [0.0,0.0,z3]
    th1,th2     = -90.0, 90.0
    yoke_h      = sec.generate__sector180( lc=lc, r1=r1, r2=r4, \
                                           zoffset=z2, height=z3-z2, defineVolu=True, side=side )
    yoke_v      = sec.generate__sector180( lc=lc, r1=r2, r2=r4, \
                                           zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    york_c      = con.generate__coneShape( lc =lc , origin=origin_cone, r1=r4, r2=r3, \
                                           th1=th1, th2=th2, height=z4-z3, side=side )
    return()


# ========================================================= #
# ===  generate outside Air region                      === #
# ========================================================= #
def generate__outAir( lc=0.2, r1=0.0, r2=1.1, r3=1.4, r4=1.5, \
                      z1=0.0, z2=0.7, z3=1.1, z4=1.2, side="+" ):
    # ------------------------------------------------- #
    # --- [1] main outside air cylinder             --- #
    # ------------------------------------------------- #
    th1,th2     = -90.0, 90.0
    origin      = [ 0.0, 0.0 ]
    origin_cone = [0.0,0.0,z2]
    oAir_h      = sec.generate__sector180( lc=lc, r1=r1, r2=r4, \
                                           zoffset=z2, height=z4-z2, defineVolu=True, side=side )
    oAir_v      = sec.generate__sector180( lc=lc, r1=r3, r2=r4, \
                                           zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    oAir_c      = con.generate__coneShape( lc =lc , origin=origin_cone, r1=r3, r2=r2, \
                                           th1=th1, th2=th2, height=z3-z2, side=side )
    # ------------------------------------------------- #
    # --- [2] cut cone shape from cylinder          --- #
    # ------------------------------------------------- #
    target      = [ (3,( (oAir_h[0] )["volu"])["sector"]),(3,( (oAir_h[1] )["volu"])["sector"]) ]
    tool        = [(3,(oAir_c["volu"])["cone"])]
    ret         = gmsh.model.occ.cut( target, tool )
    return()
