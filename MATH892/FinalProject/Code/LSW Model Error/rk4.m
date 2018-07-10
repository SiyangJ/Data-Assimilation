% 4th Order Runge-Kutta

function y_f = rk4(fname,t_i,t_f,h,ic)

t = t_i;
y_f = ic;

while (t < t_f)
    
    k_1 = feval(fname,t,y_f);
    k_2 = feval(fname,t+h/2,y_f+h*k_1/2);
    k_3 = feval(fname,t+h/2,y_f+h*k_2/2);
    k_4 = feval(fname,t+h,y_f+h*k_3);
    
    y_f = y_f+h*(k_1+2*k_2+2*k_3+k_4)/6;
    t = t+h;
end

