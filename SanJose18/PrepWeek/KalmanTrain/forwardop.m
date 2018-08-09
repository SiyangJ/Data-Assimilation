function xf=forwardop(A,Q,x)

xf=A*x+mvnrnd([0,0],Q);

