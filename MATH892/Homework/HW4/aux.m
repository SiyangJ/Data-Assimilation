figure
hold on
%plot(xt);
plot(yt');
hold off
legend('y1','y2')
%%
figure
hold on
errorbar(xp(1,:),squeeze(Pp(1,1,:)))
hold off
%%
figure
plot(xp(1,:))
%%
figure
hold on
plot(0:20,xt(1,:))
plot(0:20,xp(1,:))
errorbar(0:20,xf(1,:),squeeze(Pf(1,1,:)))
hold off
legend('true','predict','analysis')