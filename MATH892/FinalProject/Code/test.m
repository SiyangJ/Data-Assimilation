for i=2:3
    
    for j = 1:nens
        
        [~,xf] = ode45(TM,tobs(i-1:i),xens(:,j));
        xf = xf(end,:)';
        xens(:,j) = xf;
        
        
        
        
    end
    
    
    
end