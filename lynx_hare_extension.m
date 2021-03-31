function exitflag = lynx_hare_extension

%Enter the data and the initial guess.
td = [0:20];
hare = [30 47.2 70.2 77.4 36.3 20.6 18.1 21.4 22 25.4 27.1 40.3 57 76.6 52.3 19.5 11.2 7.6 14.6 16.2 24.7];
lynx = [4 6.1 9.8 35.2 59.4 41.7 19 13 8.3 9.1 7.4 8 12.3 19.5 45.7 51.1 29.7 15.8 9.7 10.1 8.6];
p = [34.913 3.856 0.397 0.018 0.786 0.023];

%Finally we use the fminsearch routine as follows:

[p,fval,exitflag] = fminsearch(@leastcomp,p,[],td,hare,lynx);
p
fval

function J = leastcomp(p,tdata,xdata,ydata)
%Create the least squares error function to be minimized.
n1 = length(tdata);
[t,y] = ode23(@lotvol,tdata,[p(1),p(2)],[],p(3),p(4),p(5),p(6));
errx = y(:,1)-xdata(1:n1)';
erry = y(:,2)-ydata(1:n1)';
J = errx'*errx + erry'*erry;

function dydt = lotvol(t,y,a1,a2,b1,b2)
% Predator and Prey Model
tmp1 = a1*y(1) - a2*y(1)*y(2);
tmp2 = -b1*y(2) + b2*y(1)*y(2);
dydt = [tmp1; tmp2];