function XDOT = FTHS(T,X)
%Couple Stommel-Lorenz (simplified ocean-atmosphere model), from
%Roebber PJ (1995) Climate variability in a low-order coupled
%atmosphere–ocean model. Tellus 47A:473–494,
%
%R. Tardif, G. J. Hakim, and C. Snyder. Coupled atmosphere-ocean data assimilation experiments with a
%low-order climate model. Climate Dynamics, pages 1–-13, 2013.

			aa=0.25E0;	
			bb=4.0E0;
			V1=0.832E0;
			V2=2.592E0;
			V3=10.3E0;
			KT=10.5E0;
			KZ=1.0E0;
			T0=128.15E0;
			alpha=1.0E0;
			beta=8.0E0;
			mu=4.0E0;
			F0=6.65E0;
			F1=2.0E0;
			F2=10.9E0;
			G0=-3.6E0;
			G1=1.24E0;
			G2=3.81E0;
			c1=1.25E0;
			c2=0.156E0;
			TA2=298.15E0; % 25 degrees celsius in kelvins
			gamm=0.06812E0;

                        omega=1.E0; % Time is in seconds, so this should be
                        omega=0.E0;  % omega = 2*pi/(365*24*60*60)
                        omega = 2*pi/(365*24*60*60);

% Definition of the coupling terms
			q=mu*(alpha*(X(5)-X(4))-beta*(X(8)-X(7)));
			FF=F0+F1*cos(omega*T)+F2*(X(5)-X(4))/T0;
			G=G0+G1*cos(omega*T)+G2*X(4)/T0;
			QS=c1+c2*(X(2)*X(2)+X(3)*X(3));
			TA1=TA2-gamm*X(1);
% Atmospheric Component
			XDOT(1)=-(X(2)*X(2)+X(3)*X(3))-aa*X(1)+aa*FF;
			XDOT(2)=X(1)*X(2)-bb*X(1)*X(3)-X(2)+G;
			XDOT(3)=X(1)*X(3)+bb*X(1)*X(2)-X(3);
% Ocean Component
			XDOT(4)=(0.5E0*q*(X(5)-X(6))+KT*(TA1-X(4))-KZ*(X(4)-X(6)))/V1;
			XDOT(5)=(0.5E0*q*(X(6)-X(4))+KT*(TA2-X(5))-KZ*(X(5)-X(6)))/V2;
			XDOT(6)=(0.5E0*q*(X(4)-X(5))+KZ*(X(4)-X(6))+KZ*(X(5)-X(6)))/V3;
			XDOT(7)=(0.5E0*q*(X(8)-X(9))-KZ*(X(7)-X(9))-QS)/V1;
			XDOT(8)=(0.5E0*q*(X(9)-X(7))-KZ*(X(8)-X(9))+QS)/V2;
			XDOT(9)=(0.5E0*q*(X(7)-X(8))+KZ*(X(7)-X(9))+KZ*(X(8)-X(9)))/V3;
