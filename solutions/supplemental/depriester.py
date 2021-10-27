import numpy as np

def depriester( coefs, T, p ):
    # K = depreiester( coefs, T, p )
    #
    # Calculates K-values using the depriester charts
    #
    # INPUTS:
    #   coefs - the Depriester correlation coefficients, with each coefficient
    #           in a column and each species in a row.
    #   T     - temperature (Rankine)
    #   p     - pressure (psia)
    #
    # OUTPUT:
    #   K - the K-values
    #
    # Author: James C. Sutherland

    ns,nc = coefs.shape
    assert( nc==6 )  # 6 coefficients are required for depriester

    K = np.exp( coefs[:,0]/T**2 + coefs[:,1]/T + coefs[:,2] + coefs[:,3]*np.log(p) + coefs[:,4]/p**2 + coefs[:,5]/p )

    return K