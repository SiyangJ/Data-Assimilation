syms DF(x,y,z) a b c x y z
DF(x,y,z) = [-a,a,0;b-z,-1,-x;y,x,-c];
E = eig(DF);

E0 = eig(DF(0,0,0));

syms x1 x2 y1 y2 z1 z2

x1 = sqrt(c*(b-1));
y1 = x1;
x2 = -x1;
y2 = x2;

z1 = b-1;
z2 = z1;

syms DF1 d

DF1 = [-a,a,0;1,-1,-d;d,d,-c];

E1 = eig(DF(x1,y1,z1));
E2 = eig(DF(x2,y2,z2));

E1d = eig(DF1);