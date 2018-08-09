function [xk1,Bk1]=predict(xk,Bk,A,Q)

xk1=A*xk+mvnrnd([0,0],Q)';
Bk1=A*Bk*A'+Q;

%  for i=1:1000
%  xk1=xk1+forward(A,Q,xk);
%  Bk1=Bk1+A*Bk*A'+Q;
%  end
%  
%  xk1=xk1*(1/1001);
%  Bk1=Bk1*(1/1001);
