import numpy as np
import matplotlib.pyplot as plt

def batch_constant_r( xeq, yeq, Teq, Nstage, x0, xw, xdguess, LV, doPlot ):
    '''
    This function calculates the distillate composition that corresponds to
    the given residue composition for a batch distiller with a rectifying
    section that has a constant reflux ratio.
    
     INPUTS:
       xeq  - equilibrium liquid mole fractions
       yeq  - equilibrium vapor mole fractions
       Teq  - equilibrium temperature
       Nstage - number of equilibrium stages
       x0   - initial charge composition in the still
       xw   - the residue composition
       xdguess - a guess for the distillate composition
       Lv      - Liquid/Vapor ratio in the rectifying section
       doPlot  - a flag to trigger plotting the McCabe-Thiele diagram. 
                 (False - don't plot, True - plot)
    
     OUTPUT:
       The distillate composition corresponding to the supplied residue
       composition.  Note that this is only consistent when it is equal to the
       supplied xdguess.
    
    Author: James C. Sutherland
    '''
    # Build functions to interpolate the Txy data. 
    # These functions are essentially shortcuts that we will use below
    eval_y = lambda x: np.interp(x,xeq,yeq)
    eval_T = lambda x: np.interp(x,xeq,Teq) 

    xd = xdguess   # current guess for the distillate composition

    '''
    Create a function to evaluate the liquid composition from the operating line.
    This comes from solving the linear equation:
       y = slope * (x-x1) + y1
    for x, with (x1,y1) a point on the line (the distillate point in this case)
    '''
    eval_xop = lambda y: (y-xd)/LV+xd

    if doPlot:
        xx = np.linspace(0,xdguess)
        yy = eval_y(xx)
        plt.figure(figsize=(7,7))
        plt.plot(xx,yy,'b-')
        plt.plot([0,1],[0,1],'k-')
        plt.grid()
        plt.axis('equal')
        plt.axis([0,1,0,1])
        plt.xlabel('Liquid Benzene Mole Fraction')
        plt.ylabel('Vapor Benzene Mole Fraction')
        yy = np.linspace(0,xd)
        plt.plot(eval_xop(yy),yy,'r--')  # the operating line

    # we start the plot at the residue composition given.
    xop = xw;
    yop = 0;

    for i in range(0,Nstage):
        # set the equilibrium composition in the liquid. This is a vertical line
        # from the operating line (stepping up), so it is equal to xop.
        xstage = xop

        # set the vapor composition in equilibrium with xstage.  This comes from
        # interpolating the equilibrium data.
        ystage = eval_y( xstage )

        if doPlot:
            plt.plot( [xop,xstage],[yop,ystage],'g-')

        # step over to the operating line.  Here we know the vapor composition,
        # and need to find the liquid from the operating line.
        yop = ystage
        xop = eval_xop( yop )

        if doPlot:
            plt.plot( [xstage,xop],[ystage,yop],'g-')
            plt.text(xstage-.03,ystage+.03,str(i));

    if doPlot:
#         print( 'Distillate composition: {:.4f}'.format(xop))
        plt.show()

    return xop