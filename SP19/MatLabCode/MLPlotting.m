figure
for i = 1:100
    
    plot(1:20,Y(:,i),1:20,HML40_L_10_30(X(:,i)))
    pause()
    
end

%%

figure
for i = 1:100
    
    plot(1:40,X(:,i),1:40,MLinverse(Y(:,i)))
    pause()
    
end