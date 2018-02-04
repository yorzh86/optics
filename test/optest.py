~ #check scipy evolution optimization

#~ example1
#~ ---------#
#~ minimize x1*x4*(x1+x2+x3)+x3
#~ when x[0]*x[1]*x[2]*x[3] <= 25.0
#~ when x1**2+x2**2+x3**2+x4**2 = 40
#~ when 1<= x[i] <=5
#~ starting guess 1,5,5,1

import numpy as np
from scipy.optimize import minimize

def objective(x):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    x4 = x[3]
    return x1*x4*(x1+x2+x3)+x3

def constraint1(x):
    return x[0]*x[1]*x[2]*x[3]-25.0

def constraint2(x):
    sum_sq = 40
    for i in range(4):
        sum_sq = sum_sq - x[i]**2
    return sum_sq  

x0 = [1,5,5,1]
print objective(x0)

b = (1.0, 5.0)
bnds = (b,b,b,b)
con1 = {'type':"ineq", 'fun': constraint1}
con2 = {'type':"eq", 'fun': constraint2}
cons = [con1, con2]

sol = minimize(objective, x0, method='SLSQP', bounds=bnds, constraints = cons)
print sol
print sol.x


#~ example2
#~ ---------#
#~ 1.3x**2+4x+0.6=0. Minimize this.
import scipy
objective2 = np.poly1d([1.3, 4.0, 0.6])
print objective2
x_ = scipy.optimize.fmin(objective2, [3])
print x_

#~ plot it:
import pylab as pl
x = np.linspace(-4, 1, 101)
pl.plot(x, objective2(x))
pl.plot(x_, objective2(x_), 'ro')



#~ example3
#~ ---------#
#~ minimize: f = 2x0 +x1
#~ when: 
#~ -x0+x1<=1
#~ x0+x1>=2
#~ x1>=0
#~ x0-2x1<=4
import cvxopt as cvx
from cvxopt import solvers as cvx_solvers

#1. Write matrix A - coefficients of x0, x1
A = cvx.matrix([  [-1.0, -1.0, 0, -1.0], [1.0, -1.0, -1.0, -2.0]  ])
#2. Write vector b - less than 1, more than 2, more than 0, less than4
b = cvx.matrix([1.0, -2.0, 0.0, 4.0])
#3. Write coefficients of fn equality
fn  = cvx.matrix([2.0, 1.0])
#4. Define and print solution
sol = cvx_solvers.lp(fn, A, b)
print sol['x']


#~ example 4
#~ solve system of equations(does not work for inequalities)
import numpy as np
import scipy.optimize as opt
# 1. Write equations
def system(x,a,b,c):
    x0,x1,x2 = x
    eqs = [
        3*x0-np.cos(x1*x2)+a, # == 0
        x0**2 -81*(x1+0.1)**2 + np.sin(x2)+b, # == 0
        np.exp(-x0*x1)+ 20*x2 + c, # == 0
    ]
    return eqs
#2. Define coefficients
a = -0.5
b = 1.06
c = (10*np.pi-3.0)/3

#3. Set initial guess
x0=[0.1, 0.1, -0.1]

#4. Solve system of nonlinear eqs and print results
result = opt.root(system, x0, args=(a,b,c))
print "root:", result.x
print "solution:", result.fun

