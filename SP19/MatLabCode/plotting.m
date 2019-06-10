%% Plotting Identical Twin Experiment for All 40 Variables
%  With Error Bars

figure
for i=1:40
    subplot(8,5,i)
    hold on
    plot(ttraj,truetraj(i,:),'b-','LineWidth',4)
    ami = xam(i,:);
    fmi = xfm(i,:);
    ame = reshape(xacov(i,i,:),[1,nobs+1]);
    fme = reshape(xfcov(i,i,:),[1,nobs+1]);
    errorbar(tobs,ami,ame,'o','MarkerSize',10,'MarkerFaceColor',[.49 1 .63])
    errorbar(tobs,fmi,fme,'r-s','LineWidth',2,'MarkerSize',8)
    plot(tobs,xmf(i,:),'g--*')
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
for i=1:40
    subplot(8,5,i)
    hold on
    plot(ttraj,truetraj(i,:),'b-','LineWidth',4)
    ami = xam(i,:);
    fmi = xfm(i,:);
    plot(tobs,ami,'r-','LineWidth',2)
    plot(tobs,fmi,'g-','LineWidth',0.5)
    %plot(tobs,xmf(i,:),'g--*')
    hold off
end

%% Plotting Identical Twin Experiment sequentially

erflag = true;

figure
for i=1:nobs+1

    subplot(2,1,1)
    plot(1:10,truestate(1:10,i),'b-','LineWidth',6)
    hold on
    plot(1:10,xam(1:10,i),'r-o','LineWidth',3)
    plot(1:10,xfm(1:10,i),'g-*','LineWidth',1)
    
    if erflag
        
        fill([1:10,fliplr(1:10)],...
        [xfm(1:10,i)' - 3* diag(xfcov(1:10,1:10,i))',...
         fliplr(xfm(1:10,i)' + 3* diag(xfcov(1:10,1:10,i))')],...
        'r','EdgeColor','none','FaceColor','g','FaceAlpha','0.3')
    
        fill([1:10,fliplr(1:10)],...
        [xam(1:10,i)' - 3* diag(xacov(1:10,1:10,i))',...
         fliplr(xam(1:10,i)' + 3* diag(xacov(1:10,1:10,i))')],...
        'r','EdgeColor','none','FaceColor','r','FaceAlpha','0.5')
    end
    
    hold off
    
    subplot(2,1,2)
    plot(11:40,truestate(11:40,i),'b-','LineWidth',6)
    hold on
    plot(11:40,xam(11:40,i),'r-','LineWidth',3)
    plot(11:40,xfm(11:40,i),'g-','LineWidth',1)
    
    if erflag
        
        fill([11:40,fliplr(11:40)],...
        [xfm(11:40,i)' - 3* diag(xfcov(11:40,11:40,i))',...
         fliplr(xfm(11:40,i)' + 3* diag(xfcov(11:40,11:40,i))')],...
        'r','EdgeColor','none','FaceColor','g','FaceAlpha','0.3')
    
        fill([11:40,fliplr(11:40)],...
        [xam(11:40,i)' - 3* diag(xacov(11:40,11:40,i))',...
         fliplr(xam(11:40,i)' + 3* diag(xacov(11:40,11:40,i))')],...
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
