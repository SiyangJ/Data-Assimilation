%       Created by Dwight Nwaigwe  April 2011 Modified by Chirstian Sampson
%       2018
%	This program uses the ensemble kalman filter to estimate a system's state.
%	The state is x_new=f(x,u)+w, where u some input, w the
%	Gaussian distributed process noise, and f is a nonlinear function. The measurement 
%	is y_new=h(x)+v where h is a nonlinear function and v Gaussian distributed measurement noise.                 


%       The algorithm used in this code is referenced from the following:
%       S Gillijns et. al., "What Is the Ensemble Kalman Filter and How Well Does it Work?"
%       Proceedings of the 2006 American Control Conference,
%       Minneapolis, Minnesota, USA, June 14-16, 2006, pp 4448-4453.


function [Xtr,Ybar,Xest]= ensemblekfilter(f,h,x_tr,x_ini,w,z,num_iterations) 



[~,num_members]=size(x_ini);
p1=length(f);
m1=length(h);
%var_vector=[]
var_vector = cell([p1,1]); 
x_est=x_ini; % set first estimates to inital values guesses

for j=1:p1  %create vector containing variables x1 to xn
  eval(sprintf(' syms x%d', j));
  %var_vector=[var_vector sprintf('x%d ',j)];
  var_vector{j}=sprintf('x%d ',j);
  
  end
var_vector=strcat('[',var_vector);
var_vector=strcat(var_vector,']');

Zcov=eye(m1); %create measurement noise covariance matrix this will act a a preconditioner later, not always needed 
for j=1:m1
  Zcov(j,j)=z(j)^2;
  
end
%Each iteration is a time step

for i=1:num_iterations  
i
   x_tr=double(subs(f,str2sym(var_vector),x_tr))+w.*randn(p1,1); %compute true value of state at next time step
   Xtr(:,i)=x_tr; % save the true value
   for j=1:num_members
       
     W(:,j)=w.*randn(p1,1);                          %create process noise
     Z(:,j)=z.*randn(m1,1);                          %create measurement noise
     x_est(:,j)=double(subs(f,str2sym(var_vector),x_est(:,j))); %+W(:,j);      %forecast state
     y(:,j)=double(subs(h,str2sym(var_vector),x_tr))+Z(:,j);                 %make measurement
     y_for(:,j)=double(subs(h,str2sym(var_vector),x_est(:,j)));              %forecast measurement
   end

   x_estbar=mean(x_est,2); % calculate mean of ensemble members                    
   ybar=mean(y,2); 
   Ybar(:,i)=ybar;%store the observation
   y_forbar=mean(y_for,2); %calculate mean of elsemble observation maps
   
   for j=1:p1
     Ex(j,:)=[x_est(j,:)-x_estbar(j)]; % Create the ensemble forcast error matrix
     
   end

   for j=1:m1
     Ey(j,:)=[y_for(j,:)-y_forbar(j)]; % Create the ensemble forecast observation error matrix 
     
   end

   Bxy=Ex*Ey'/(num_members-1); %Estimate forecast/obs error covariance matrix
   Byy=Ey*Ey'/(num_members-1)+ Zcov;  %Estimate observation error covariance matrix with our observation operator Byy is pretty singular not always so.                    
   K=Bxy*inv(Byy); %calculate the Kalman Gain
   x_est=x_est+K*(y-y_for); % Preforem the analysis step on each forecast ensemble member.
   x_estbar=mean(x_est,2); % Calculate the averave to be our "analysis"
   Xest(:,i)=x_estbar; %Store mean the value for iteration i

end