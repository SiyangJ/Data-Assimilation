% 4th order Runge-Kutta

function y_f = rk4_full(fname,t_i,t_f,h,ic)

t = t_i;
time = t_i:h:t_f;
y_f = zeros(length(ic),length(time));
y_f(:,1) = ic;
i = 1;

while (t < t_f)
    
    k_1 = feval(fname,t,y_f(:,i));
    k_2 = feval(fname,t+h/2,y_f(:,i)+h*k_1/2);
    k_3 = feval(fname,t+h/2,y_f(:,i)+h*k_2/2);
    k_4 = feval(fname,t+h,y_f(:,i)+h*k_3);
    
    y_f(:,i+1) = y_f(:,i)+h*(k_1+2*k_2+2*k_3+k_4)/6;
    t = t+h;
    i = i+1;
end