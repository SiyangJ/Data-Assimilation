function [xak,Bak]=Kupdate(xk,yk,Bk,H,R)


K=Bk*H'*inv(H*Bk*H'+R);

xak=xk+K*(yk-H*xk);
Bak=Bk-K*H*Bk;
