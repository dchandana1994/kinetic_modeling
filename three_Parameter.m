function dx = f(t,x)
global b;
dx(1,1) = -b(1)*( ( x(1)/b(3) ) - x(1) );
dx(2,1) = -b(1)*(50/49.78)*( ( x(1)/b(3) ) - x(1) ) - b(2)*x(2)*x(3);
dx(3,1) = -b(2)*x(2)*x(3);
end
function S = Sfun1D(b)
% computation of an error function for an ODE model
% INPUT: b - vector of parameters
global tdata xdata x0;
%% ODE model
% (nested function, uses parameters b(1) and b(2) of the main function)
%% numerical integration set up
tspan = [0:0.1:max(tdata)];
[tsol,xsol] = ode23(@f,tspan,x0);
%% plot result of the integration
plot(tdata,xdata,'x','MarkerSize',10);
hold on
plot(tsol,(xsol(:,1)));
hold off
drawnow
%% find predicted values x(tdata)
xpred = interp1(tsol,xsol(:,1),tdata);
%% compute total error
S = 0;	
for i = 1:length(tdata)
S = S + ( xpred(i)-xdata(i) )^2;
end
end
function paramfit1D
% main program for fitting parameters of an ODE model to data
% the model and the error function are defined in the file Sfun1D.m
clearvars -global
global tdata xdata x0 b;
%% declare data for the model
% time - x value
tdata(1) = 0; xdata(1) = 3.1;
tdata(2) = 5; xdata(2) = 2.35;
tdata(3) = 10; xdata(3) = 2;
tdata(4) = 20; xdata(4) = 1.6;
tdata(5) = 30; xdata(5) = 0.8;
tdata(6) = 40; xdata(6) = 0.45;
tdata(7) = 50; xdata(7) = 0.25;
%% initial condition
x0(1) = 3.1;
x0(2) = 0;
x0(3) = 7.57;
%% initial guess of parameter values
b(1) = 0.043;
b(2) = 1.5;
b(3) = 0.5;
%% minimization step
[bmin, Smin] = fminsearch(@Sfun1D,b);
disp('Estimated parameters b(i):');
disp(bmin)
disp('Smallest value of the error S:');
disp(Smin)
end
