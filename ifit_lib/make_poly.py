import numpy as np

#========================================================================================
#===================================== make_poly ========================================
#======================================================================================== 

def make_poly(grid, coefs):

    '''
    Function to construct a polynolial given a list of coefficients
    
    INPUTS
    ------
    grid, array
        X grid over which to calculate the polynomial
        
    coefs, list
        Polynomial coefficents (of acesnding rank)
    
    OUTPUTS
    -------
    poly, array
        Resulting polynomial, calculated as poly = p0 + p1*grid + p2*grid^2 + ...
    '''

    poly = np.zeros(len(grid))

    for i, p in enumerate(coefs):

        poly = np.add(poly, np.multiply(np.power(grid, i), p))
        
    return poly
