% for i = 1:5
% 
%     DataGeneration(i);
%     
% end

%%

X = zeros(40,1000);
Y = zeros(20,1000);

for i=1:5
    
    [obs,states] = DataGeneration2(i);
    
    X(:,(i-1)*200+1:i*200) = states(:,2:end);
    Y(:,(i-1)*200+1:i*200) = obs(:,2:end);
    
end

save('Data/X_10_30_100.mat','X')
save('Data/Y_10_30_100.mat','Y')