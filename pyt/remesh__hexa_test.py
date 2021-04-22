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
inpFile = "msh/model.step"
gmsh.model.occ.importShapes( inpFile )
gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()

# ------------------------------------------------- #
# --- [4] Physical Grouping                     --- #
# ------------------------------------------------- #
pntFile = "dat/physNumTable.conf"
import nkGmshRoutines.load__physNumTable as lpt
lpt.load__physNumTable( inpFile=pntFile )

# ------------------------------------------------- #
# --- [5] Mesh settings                         --- #
# ------------------------------------------------- #

#  -- [5-1] mesh config loading                 --  #
mshFile  = "dat/mesh.conf"
meshdict = {}
with open( mshFile, "r" ) as f:
    lines = f.readlines()
for line in lines:
    if   ( len( line.strip() ) == 0 ):
        pass
    elif ( ( line.strip() )[0] == "#" ):
        pass
    elif ( len( ( line.strip() ).split() ) == 3 ):
        aline              = ( line.strip() ).split()
        meshdict[aline[0]] = [ int( aline[1] ), float( aline[2] ) ]

# -- [5-2] mesh alignment                      --  #

physNumKeys = list( meshdict.keys() )
print( physNumKeys )
# print( [ meshdict[key] for key in physNumKeys ] )
# sys.exit()

physNumKeys = [ "gap", "pole", "coil", "outAir" ]
physNumKeys = [ "yoke" ]

for key in physNumKeys:
    # -- setting -- #
    physNum    = ( meshdict[key] )[0]
    meshsize   = ( meshdict[key] )[1]
    # -- volume  -- #
    v_entities = gmsh.model.getEntitiesForPhysicalGroup( voluDim, physNum )
    v_dimtags  = [ ( voluDim, v_enti ) for v_enti in v_entities ]
    # -- surface -- #
    s_dimtags  = gmsh.model.getBoundary( v_dimtags, oriented=False )
    s_entities = [ int( s_dimtag[1] ) for s_dimtag in s_dimtags ]
    # -- line    -- #
    l_dimtags  = gmsh.model.getBoundary( s_dimtags, combined=False, oriented=False )
    l_dimtags  = list( set( l_dimtags ) )
    l_entities = [ int( l_dimtag[1] ) for l_dimtag in l_dimtags ]
    
    # -- inifinite ==> line   -- #
    NN = 6
    is_rectangle  = True
    nSmooth       = 100
    for il,l_enti in enumerate(l_entities):
        gmsh.model.mesh.setTransfiniteCurve( l_enti, NN )
    # for iS,s_enti in enumerate(s_entities):
    #     hlines    = gmsh.model.getBoundary( (surfDim,s_enti), combined=False, oriented=False )
    #     if ( len( hlines ) <= 4 ):
    #         gmsh.model.mesh.setTransfiniteSurface( s_enti )
    #     if is_rectangle:
    #         gmsh.model.mesh.setRecombine( surfDim, s_enti )
    #         gmsh.model.mesh.setSmoothing( surfDim, s_enti, nSmooth)
    for iv,v_enti in enumerate(v_entities):
        hsurfs    = gmsh.model.getBoundary( (voluDim,v_enti), oriented=True )
        checker   = True
        print( "(key,v_enti,len(hsurfs)) = {0}, {1}, {2}".format( key, v_enti, len(hsurfs) )  )
        if ( len( hsurfs ) <= 6 ):
            for iS,s_dimtag in enumerate(hsurfs):
                hlines = gmsh.model.getBoundary( s_dimtag, combined=False, oriented=False )
                if ( len( hlines ) > 4 ):
                    checker = False
                    break
            if ( checker is True ):
                print( key, meshdict[key] )
                for iS,s_dimtag in enumerate(hsurfs):
                    gmsh.model.mesh.setTransfiniteSurface( s_dimtag[1] )
                    if is_rectangle:
                        gmsh.model.mesh.setRecombine( surfDim, s_dimtag[1] )
                        gmsh.model.mesh.setSmoothing( surfDim, s_dimtag[1], nSmooth)
                gmsh.model.mesh.setTransfiniteVolume( v_enti )
        
    # for iv,v_enti in enumerate(v_entities):
    #     v_dimtag   = [(voluDim,v_enti)]
    #     s_dimtags  = gmsh.model.getBoundary( v_dimtag, oriented=False )
    #     s_entities = [ int( s_dimtag[1] ) for s_dimtag in s_dimtags ]
        
    #     for iS,s_enti in enumerate( s_entities ):
    #         v_dimtag   = [(voluDim,v_enti)]
    #         s_dimtags  = gmsh.model.getBoundary( v_dimtag, oriented=False )
            
    #         # print( surf )
        
# print( physNumKeys )
# print( meshdict )

        
# ------------------------------------------------- #
# --- [6] post process                          --- #
# ------------------------------------------------- #

gmsh.option.setNumber( "Mesh.CharacteristicLengthMin", 0.1 )
gmsh.option.setNumber( "Mesh.CharacteristicLengthMax", 0.2 )

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write( "msh/imported.msh" )
gmsh.finalize()

