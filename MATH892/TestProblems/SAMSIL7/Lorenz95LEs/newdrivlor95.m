%This illustrates the use of the discrete QR method to approximate
%Lyapunov exponents and a orthonormal basis for the Lyapunov vectors
%associated with the first p Lyapunov exponents.
%The use of a general solution operator (formod) and finite differences
%allow for a Jacobian free method.
%
%Final time
TF=100;
%DIMENSION OF PROBLEM
DIM=40;
%FIXED TIME STEP
H = 1.E-2;
%NUMBER OF LYAPUNOV EXPONENTS
p=20;

%IC=e1
X = zeros(DIM,1);
X(1)=1;

%INITIALIZE
LE = zeros(p,1);
q = eye(DIM,p);
FEVAL = zeros(DIM,1);
NEWDIFF = zeros(DIM,p);
TIME = 0;
TAU = H/10;
 
%NUMBER OF TIMESTEPS
NUMSTEPS = floor(TF/H);

%INTEGRATE LORENZ '95 

for i = 1:NUMSTEPS

%UPDATE TRAJECTORY, I.E., FIND F(X).
   XNEW = formod(TIME,X,H);

   epsilon=sqrt(eps(1))*max(1,norm(XNEW));

%% OVER ALL THE LYAPUNOV EXPONENTS WE WANT, SOME p<=DIM
   for j=1:p

%EVALUATE F(X+H*Qj)
      NEWIC = X+epsilon*q(:,j);
      FEVAL = formod(TIME,NEWIC,H);

%APPROXIMATE F'(X_n)Q_n
      NEWDIFF(:,j) = (FEVAL - XNEW)/epsilon;

   end

%CALL mgs
   [q,r] = mgs(NEWDIFF);

%FORM LES
   for j=1:p
       LE(j) = LE(j) + log(r(j,j));
   end 

%UPDATE TIME AND DEPENDENT VARIABLES
   X = XNEW;
   TIME = TIME + H;

end

LE = LE/TIME

