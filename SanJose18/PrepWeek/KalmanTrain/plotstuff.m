plot(xtrue)
hold on
display('This is the true trajectory, Hit enter to continue plotting')
pause;
plot(xGPS)
hold on
display('These are the location observations, Hit enter to continue plotting')
pause;

plot(avgspeed100)
display('This is the position if the train went an average of 100mph as the train company claims, Hit enter to continue plotting.')
hold on
pause;

plot(xk)
display('These are each prediction from a precious analysis value, Hit enter to continue plotting')
hold on
pause;

plot(Xa)
display('These are the analysis values themselves, Hit enter to continue plotting')
hold on
pause;

figure
plot(vtrue)
hold on
display('This is the true velocities, must have ben dangerous track conditions, Hit enter to continue plotting')
pause;

plot(vGPS)
hold on
display('This is the velocity GPS observations, Hit enter to continue plotting')
pause;

plot(100*ones(1,60))
display('This is the if we had an average velocity of 100mph, Hit enter to continue plotting')
pause;

plot(vk)
display('These are the predictited velocities, from the last analysis step, Hit enter to continue plottting!')
hold on
pause;

plot(Va)
display('These are the analysis velocities, thats all!')
hold on
pause;




