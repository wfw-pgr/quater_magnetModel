import os, sys, subprocess
import gmsh
import numpy                              as np
import nkGmshRoutines.generate__fanShape  as fan
import nkGmshRoutines.generate__coneShape as con


# ========================================================= #
# === generate magnet shape                             === #
# ========================================================= #

def generate__magnet( lc=0.4, side="+" ):

    cnsFile = "dat/parameter.conf"
    import nkUtilities.load__constants as lcn
    const = lcn.load__constants( inpFile=cnsFile )

    # ------------------------------------------------- #
    # --- [1] pole making                           --- #
    # ------------------------------------------------- #
    generate__pole( lc=lc , r1=0.0  , r2=const["r_pole"], \
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
    generate__coilslot( lc=lc, r1=r1, r2=r2, r3=r3, r4=r4, \
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
    generate__yoke    ( lc=lc, r1=r1, r2=r2, r3=r3, r4=r4, \
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
    generate__outAir  ( lc=lc, r1=r1, r2=r2, r3=r3, r4=r4, \
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
    th1,th2 = -90.0, 90.0
    gap     = fan.generate__fanShape( lc=lc, th1=th1, th2=th2, origin=origin, r1=r1, r2=r2, \
                                      zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    pole    = fan.generate__fanShape( lc=lc, th1=th1, th2=th2, origin=origin, r1=r1, r2=r2, \
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
    origin  = [ 0.0, 0.0 ]
    th1,th2 = -90.0, 90.0
    air_inn = fan.generate__fanShape( lc=lc, th1=th1, th2=th2, origin=origin, r1=r1, r2=r2, \
                                      zoffset=z1, height=z4-z1, defineVolu=True, side=side )
    air_bot = fan.generate__fanShape( lc=lc, th1=th1, th2=th2, origin=origin, r1=r2, r2=r3, \
                                      zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    coil    = fan.generate__fanShape( lc=lc, th1=th1, th2=th2, origin=origin, r1=r2, r2=r3, \
                                      zoffset=z2, height=z3-z2, defineVolu=True, side=side )
    air_top = fan.generate__fanShape( lc=lc, th1=th1, th2=th2, origin=origin, r1=r2, r2=r3, \
                                      zoffset=z3, height=z4-z3, defineVolu=True, side=side )
    air_out = fan.generate__fanShape( lc=lc, th1=th1, th2=th2, origin=origin, r1=r3, r2=r4, \
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
    th1,th2     = -90.0, 90.0
    origin_cone = [0.0,0.0,z3]
    york_h      = fan.generate__fanShape(  lc=lc, th1=th1, th2=th2, origin=origin, r1=r1, r2=r4, \
                                           zoffset=z2, height=z3-z2, defineVolu=True, side=side )
    york_v      = fan.generate__fanShape(  lc=lc, th1=th1, th2=th2, origin=origin, r1=r2, r2=r4, \
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
    oAir_h      = fan.generate__fanShape(  lc=lc, th1=th1, th2=th2, origin=origin, r1=r1, r2=r4, \
                                           zoffset=z2, height=z4-z2, defineVolu=True, side=side )
    oAir_v      = fan.generate__fanShape(  lc=lc, th1=th1, th2=th2, origin=origin, r1=r3, r2=r4, \
                                           zoffset=z1, height=z2-z1, defineVolu=True, side=side )
    oAir_c      = con.generate__coneShape( lc =lc , origin=origin_cone, r1=r3, r2=r2, \
                                           th1=th1, th2=th2, height=z3-z2, side=side )
    # ------------------------------------------------- #
    # --- [2] cut cone shape from cylinder          --- #
    # ------------------------------------------------- #
    target      = [(3,(oAir_h["volu"])["fan"])]
    tool        = [(3,(oAir_c["volu"])["cone"])]
    ret         = gmsh.model.occ.cut( target, tool )
    return()
