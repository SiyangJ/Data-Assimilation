%% Plotting Identical Twin Experiment for All 40 Variables
%  With Error Bars

figure
for i=1:2
    subplot(1,2,i)
    hold on
    plot(tobs,truestate(i,:),'b-','LineWidth',3)
    ami = xam(i,:);
    fmi = xfm(i,:);
    ame = sqrt(reshape(xacov(i,i,:),[1,nobs+1]));
    fme = sqrt(reshape(xfcov(i,i,:),[1,nobs+1]));
    
    plot(tobs,ami,'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    plot(tobs,fmi,'r-s','LineWidth',2,'MarkerSize',8)
    
    fill([tobs,fliplr(tobs)],...
        [xfm(i,:) - 3* fme,...
         fliplr(xfm(i,:) + 3* fme)],...
        'r','EdgeColor','none','FaceColor','g','FaceAlpha','0.3')
    
    fill([tobs,fliplr(tobs)],...
        [xam(i,:) - 3* ame,...
         fliplr(xam(i,:) + 3* ame)],...
        'r','EdgeColor','none','FaceColor','r','FaceAlpha','0.5')
    
    plot(tobs,xmf(i,:),'k--*')
    hold off
end

%% Plotting Identical Twin Experiment for All 40 Variables

figure
for i=1:40
    subplot(8,5,i)
    hold on
    plot(ttraj,truetraj(i,:),'b-','LineWidth',4)
    ami = xam(i,:);
    fmi = xfm(i,:);
    plot(tobs,ami,'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    plot(tobs,fmi,'r-s','LineWidth',2,'MarkerSize',8)
    plot(tobs,xmf(i,:),'g--*')
    hold off
end
%% Plotting Identical Twin Experiment for All 40 Variables
figure
for i=1:2
    subplot(1,2,i)
    hold on
    plot(tobs,truestate(i,:),'b-','LineWidth',2)
    ami = xam(i,:);
    fmi = xfm(i,:);
    plot(tobs,ami,'r-','LineWidth',2)
    plot(tobs,fmi,'g-','LineWidth',0.5)
    %plot(tobs,xmf(i,:),'g--*')
    hold off
end

%% Plotting Identical Twin Experiment sequentially

erflag = true;

mag = [1;100];

figure
for i=1:nobs+1

    plot(1:ndim,truestate(1:ndim,i) .* mag,'b-','LineWidth',6)
    hold on
    plot(1:ndim,xam(1:ndim,i) .* mag,'r-o','LineWidth',3)
    plot(1:ndim,xfm(1:ndim,i) .* mag,'g-*','LineWidth',1)
    
    if erflag
        
        fill([1:ndim,fliplr(1:ndim)],...
        [(xfm(1:ndim,i).*mag)' - 3* (diag(xfcov(1:ndim,1:ndim,i)).*mag.^2)',...
         fliplr((xfm(1:ndim,i).*mag)' + 3* (diag(xfcov(1:ndim,1:ndim,i)).*mag.^2)')],...
        'r','EdgeColor','none','FaceColor','g','FaceAlpha','0.3')
    
        fill([1:ndim,fliplr(1:ndim)],...
        [(xam(1:ndim,i).*mag)' - 3* (diag(xacov(1:ndim,1:ndim,i)).*mag.^2)',...
         fliplr((xam(1:ndim,i).*mag)' + 3* (diag(xacov(1:ndim,1:ndim,i)).*mag.^2)')],...
        'r','EdgeColor','none','FaceColor','r','FaceAlpha','0.5')
    end
    
    hold off
    
    pause();
    
end


%% Plotting the forecast/analysis RMSE for each variable

frmss = zeros(40,1);
armss = zeros(40,1);

for i=1:40
    
    frmss(i) = rms(truestate(i,:)-xfm(i,:));
    armss(i) = rms(truestate(i,:)-xam(i,:));
    
end

plot(1:40,frmss,1:40,armss)
legend({'forecast','analysis'})
