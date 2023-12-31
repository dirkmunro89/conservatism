#
import numpy as np
#
def obj(x):
#
    a = 1.; b = 32.; c = -1.
    f = a*np.sin(b*x)*np.exp(c*x)
    df = a*c*np.sin(b*x)*np.exp(c*x) + a*b*np.exp(c*x)*np.cos(b*x)
#   ddf = 2*a*c*b*np.cos(b*x)*np.exp(c*x)+a*c*c*np.sin(b*x)*np.exp(c*x)-a*b*b*np.exp(c*x)*np.sin(b*x)
#
#   quad. approx. to function linearised in terms of reciprocal intervening variables
#   ddf = -2./x*df 
    ddf = abs(-2./x*df ) # nonconvex variant
#
    return [f, df, ddf]
#
def con(x):
#
    a = 1./4.; b = 32.
    g = a*np.cos(b*x)+0.1
    dg = -a*b*np.sin(b*x)
#   ddg = -a*b*b*np.cos(b*x)
#
#   quad. approx. to function linearised in terms of reciprocal intervening variables
#   ddg = -2./x*dg 
    ddg = abs(-2./x*dg)  # nonconvex variant
#
    return [g, dg, ddg]
