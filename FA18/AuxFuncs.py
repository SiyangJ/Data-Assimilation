import numpy as np

'''
l96rhs Summary:
    Input in the format for odeint to work
    pars: for now a single scalar that represents constant F
Note:
    Implemented by Amit Apte for MATH892 Spring 18
'''
def l96rhs(xin, tin, pars):
  F = pars
  ndim, = xin.shape
  # index of the current grid point
  k00 = np.arange(ndim)
  # index for the next grid point
  kp1 = np.mod(np.arange(1,ndim+1), ndim)
  # index for the previous grid point
  km1 = np.mod(np.arange(-1,ndim-1), ndim)
  # index for previous to previous grid point
  km2 = np.mod(np.arange(-2,ndim-2), ndim)
  frhs = (xin[kp1] - xin[km2]) * xin[km1] - xin[k00] + F
  return frhs

'''
ObsOp_40_20 Summary:
    X: each column is a state vector.
    H: R^40 -> R^10 -> R^20
    Projection followed by nonlinear transformation
Note:
    Modified from matlab file below
'''
def ObsOp_40_20(x):
    
    '''
    # NumPy is stupid
    # Its default is row vector, not column vector...
    # I mean WTF???
    # Alright I'll just treat a single vector somewhat differently.
    if len(x.shape)==1:
        proj_ind = np.arange(10)
        x_proj = x[proj_ind]
        p = x_proj[2*i];
        q = x_proj[2*i+1];
    
        p2 = p*p;
        q2 = q*q;
    
        y[4*i:4*i+4] = [p*q,p2+q2,p2,q2];
        return y
    '''
    
    n_col = x.shape[1]

    # projection step R^40 -> R^10
    proj_ind = np.arange(10)
    x_proj = x[proj_ind,:]

    # nonlinear transformation R^10 -> R^20
    y = np.zeros([20,n_col]);

    # TODO: Optimization in terms of efficiency.
    for j in range(n_col):
        for i in range(5):
            p = x_proj[2*i,j];
            q = x_proj[2*i+1,j];
    
            p2 = p*p;
            q2 = q*q;
    
            y[4*i:4*i+4,j] = np.array([p*q,p2+q2,p2,q2]).T;
    return y
'''
%OBS Summary of this function goes here
%   x: each column is a state vector.
%   H: R^40 -> R^10 -> R^20
%   Projection followed by nonlinear transformation

n_col = size(x,2);

% projection step R^40 -> R^10

proj_ind = 1:10;
x_proj = x(proj_ind,:);

% nonlinear transformation R^10 -> R^20

y = zeros(20,n_col);

% TODO: Optimization in terms of efficiency.
for j=1:n_col
for i=1:5
        
    p = x_proj(2*i-1,j);
    q = x_proj(2*i,j);
    
    p2 = p^2;
    q2 = q^2;
    
    y(4*i-3:4*i,j) = [p*q;p2+q2;p2;q2];   
end
end

end
'''

'''
Inv_20_10 Summary:
    Y: each column is an observation vector.
    X: 10 dimensional state vector
    H: R^40 -> R^10 -> R^20
    Always take the positive root
Note:
    Modified from matlab file below
'''
def Inv_20_10(y):
    n_col = y.shape[1];

    x = np.zeros([10,n_col])

    for j in range(n_col):
        for i in range(5):
            p2 = y[4*i+2,j];
            q2 = y[4*i+3,j];

            x[2*i:2*i+2,j]=np.array([np.sqrt(p2),np.sqrt(q2)]).T

    return x

'''    
%STUPID_INVERSE_20_10 Summary of this function goes here
%   Detailed explanation goes here

n_col = size(y,2);

% nonlinear transformation R^10 -> R^20

x = zeros(10,n_col);

% TODO: Optimization in terms of efficiency.
for j=1:n_col
for i=1:5
        
    p = y(4*i-1,j);
    q = y(4*i,j);
    
    x(2*i-1:2*i,j) = [sqrt(p);sqrt(q)];
'''
