syms A A0 a w t U Ut Ud B V X Y
A0 = [-1,-a;0,-1];
U = [cos(w*t),sin(w*t);-sin(w*t),cos(w*t)];
Ut = U.';
Ud = diff(U,t);
V = simplify(Ud*Ut);
%%
A = simplify(Ut*A0*U);
%%
B = A0+V;
Y = simplify(expm(B*t));
X = simplify(Ut*Y);
%%
syms Xd Xr
Xd = simplify(diff(X,t));
Xr = simplify(A*X);