function y=obs(H,x,R)

y=H*x+mvnrnd([0,0],R)';


