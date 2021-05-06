import numpy as np

# ========================================================= #
# ===  generatesampleSurface                            === #
# ========================================================= #

def generate__sampleSurface():

    x_,y_,z_     = 0, 1, 2
    xMin,xMax,LI = -1.0, +1.0, +31
    yMin,yMax,LJ = -1.0, +1.0, +31
    
    # ------------------------------------------------- #
    # --- [1] coordinate                            --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ xMin, xMax, LI ]
    x2MinMaxNum = [ yMin, yMax, LJ ]
    x3MinMaxNum = [  0.0,  0.0,  1 ]
    grid        = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )

    # ------------------------------------------------- #
    # --- [2] function                              --- #
    # ------------------------------------------------- #
    def surf_func( xpos, ypos, radius=1.0 ):
        radii_h = ( xpos**2 + ypos**2 ) / radius
        phi     = np.arctan2( ypos, xpos )
        zpos    = 0.1 * ( 1.0 - radii_h ) + 0.3
        return( zpos )

    # ------------------------------------------------- #
    # --- [3] save sample surface                   --- #
    # ------------------------------------------------- #
    grid[:,:,:,z_]  = surf_func( grid[:,:,:,x_], grid[:,:,:,y_] )
    grid            = np.reshape( grid, (1,LJ,LI,3) )
    mshape          = np.zeros( (1,LJ,LI,6) )
    mshape[...,0:3] = grid
    mshape[...,  3] = grid[...,2]
    mshape[...,  4] = 0
    mshape[...,  5] = 1.0
    outFile         = "dat/mshape_svd.dat"
    import nkUtilities.save__pointFile as spf
    spf.save__pointFile( outFile=outFile, Data=mshape )



# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #
if ( __name__=="__main__" ):
    generate__sampleSurface()
