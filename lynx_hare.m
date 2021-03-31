%LYNX HARE MODEL
%this script investigates the population dynamics of the predator-prey
%relationship within the populations of snoeshow hare adn cancada lynx as
%described by the Lotka-Volterra model

%the modelling function exists as a local function within this script
%% 1) EXPLORE THE LV MODEL
%set up code to simulate system usinf the ODE45 solver

% set up the parameters and intial conditions,
options = odeset('InitialStep',0.1,'MaxStep',0.1);
t_range= [0 20];
x_ini= [15000 1000];

k1= 1.1; %Rate of prey bith
k2= 0.001;  %Rate of predation
k3 = 10;  %Rate of predator death
p = [k1 k2 k3];

% Simulation: use the model_PP function to run the simulation
[t,x]=ode45(@model_PP,t_range,x_ini,options,p);
% Plot time series of all variables
figure(1);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level ');
title('Initial model simulation@Lotka-Volterra Predator-Prey Dynamics');

%% 1.2) CHECK THE MODEL WORKS AS EXPECTED
%try the model with only prey (no predators)
x_ini = [15000 0];

% Simulation: use the model_PP function to run the simulation
[t,x]=ode45(@model_PP,t_range,x_ini,options,p);
% Plot time series of all variables
figure(2);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level ');
title('Model Behaviour: no predator population present');
%the hare population rises exponentially - as expected without predator threat 

%try the model with only predators (no prey)
x_ini = [0 1000];

% Simulation: use the model_PP function to run the simulation
[t,x]=ode45(@model_PP,t_range,x_ini,options,p);
% Plot time series of all variables
figure(3);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level ');
title('Model Behaviour: no prey population present');
%the lynx population dies off exponentially - as expected without available prey 

%the model appears to be functioning as expected. 

%% 1.3) Identify the steady states of the model
%that is the values at which the populations remain stable. 

%Population X will be steady state at the following values:
%k1*X - k2*X*Y = 0
%factorises to X( k1 - k2Y ) = 0
%solving for X = 0 or when Y = k1/k2

%Population Y will be steady state at the following values:
%k2*X*Y - k3*Y = 0
%factorises to Y( k2X - k3 ) = 0
%solving for Y = 0 or when X = k3/k2

%check that the system stablises as expected  using those values
x_ini = [k3/k2 k1/k2];

% Simulation: use the model_PP function to run the simulation
[t,x]=ode45(@model_PP,t_range,x_ini,options,p);
% Plot time series of all variables
figure(4);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level ');
title('Model Behaviour: investigating steady state values')

%the system reaches a steady state with X at 10000 and Y at 1100 - steady
%state parameter functioning as expected. 


%% 1.4) Systematically varying the initial conditions from the steady state
variations = (1:100:1000);
for i = 1: length(variations)
    x_ini = [k3/k2+(variations(i)) k1/k2+(variations(i))];
    [t,x]=ode45(@model_PP,t_range,x_ini,options,p);
    
    figure(4+i);
    plot(t,x(:,:));
    xlabel('Time ');
    ylabel('Population Level ');
    title('Varying the Initial Conditions')
end
  
%with increaseed perturbation of the initial conditionss from the steady
%state values (leaving the rate constants unchanged) the following is
%observed: 
%the periods of the population cycles decrease(slowly)
% the amplitude of the population cycles increase (quickly)

%% 1.5) Systematically varying each of the rate constants
x_ini =[k3/k2 k1/k2]; %reset the initial values to run as steady state

%investigating k1 - rate of prey birth  - initial rate = 1.1
k1_vector = [0.8, 1.1, 1.4];
p = [k1 k2 k3];
for i = 1:length(k1_vector)
    p(1) = k1_vector(i);
    [t,x]=ode45(@model_PP,t_range,x_ini,options,p);
    
    figure(14+i);
    plot(t,x(:,:));
    xlabel('Time ');
    ylabel('Population Level ');
    title('Varying K1 - the rate of prey birth')
  
end

%changing k1 to drop it below 1 results in cycles of greater amplitude in
%the prey species, 9.5-11k and a reduced population in predators, oscillating
%between 0-1k
% changing k1 to bring it above 1.1 results in similar fluctuates in the
% prey species with a slightly lower amplitude - 9.5-10.5k but the predator population is not reduced instead
% oscilliatng between 1-2k


%investing k2 - predation rate - inital rate at 0.001
k2_vector = [0.0008, 0.001, 0.0012];
p = [k1 k2 k3];

for i = 1:length(k2_vector)
    p(2) = k2_vector(i);
    [t,x]=ode45(@model_PP,t_range,x_ini,options,p);
    
    figure(17+i);
    plot(t,x(:,:));
    xlabel('Time ');
    ylabel('Population Level ');
    title('Varying K2 - the predation rate')
end
%changing k2 to make it lower results in upward oscillations of prey,
%10-16k and of the predator population 800-2200
%changing k to make it larger results in downward oscillations of the prey
%population 7-10k and a decreased amplitude of cyle of the predator population (compared with the down swing)


%investing k3 - predator death - inital rate at 10
k3_vector = [8, 10, 12];
p = [k1 k2 k3];

for i = 1:length(k3_vector)
    p(3) = k3_vector(i);
    [t,x]=ode45(@model_PP,t_range,x_ini,options,p);
    
    figure(20+i);
    plot(t,x(:,:));
    xlabel('Time ');
    ylabel('Population Level ');
    title('Varying K3 - the predator death rate')
end
%changing k3 to make it lower results in upward oscillations of predator
%1-2k and downward oscillation of the prey population (10-6.5 k)
%changing k to make it larger results in upward oscillations of the prey
%population 10-14k and a reduction in predator population (500-1500)


%% FIT THE MODEL TO THE HARE-LYNX DATA
%2.1) estimating ecological k1 point range
%when predators are absent (0) prey growth should be exponential
%based on the exponential growth model X(t) = X_0e^(k1*t)
%use this and begin by finding an ecologically valid k1 range estimate

%N0 = N(hares at start of year) = 2 %one breeding pair to begin
%N0 = 2;
%n_newborn = 3; %average litter size
%litter_number = 2.6; %average number of litters per year
%new_offspring = n_newborn*litter_number; %number of offspring per year
%N1 = N0 + new_offspring; %number of off spring after 1 year
%N1 = N0 * exp(k1*t); where t=1 for 1 year.
% new_offspring +N0 = N0* exp(k1*t); %here t = 1 as 1 year has passed
% new_offspring + N0 = N0* exp(k1)
% logN0*new_offspring+N0 = k1 %max growth (no death, max litters, max number per litter)
% new offspring range from 0 -7.8 (3 per litter 2.6 per year)
new_offspring = (0:7); %range of ecologically valid possible offspring per year
N0 = 2; %an intial breeding pair, log N0 = log2
k1 = log2(new_offspring + N0); % point range of ecologically valid 

k1 = k1(5); %choose one of the ecological valid point estimates for sanity check
x_ini =[15000 1000]; %reset the initial values
p = [k1 k2 k3]; %reset params to hold the ecologically valid k1 value

%all values are above 1 (NB for prey birth) and all are  reasonaable)
% Simulation: use the model_PP function to run the simulation
[t,x]=ode45(@model_PP,t_range,x_ini,options,p);
% Plot time series of all variables
figure(24);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level');
title('Model behaviour when k1 is ecologically valid ');

% visual assessment suggests that:
%X is oscilating around 11800 and Y is around 3800


%2.2) finding the steady state of the X and Y populations from the data
%load in the behavioural file 
data=readtable('Fur_Pelts_1900_to_1920.csv');

year = table2array(data(:,1));
hare = table2array(data(:,2));
lynx = table2array(data(:,3));

figure(25)
plot(year,hare)
hold on
plot(year, lynx)
title('Fur Pelts 1900 to 1920') 
% visual assessment suggests that:
%X is oscilating around 50000 and Y is around 30000
k1 = median(k1); %take the median birth rate
Xss = mean(hare);
Yss = mean(lynx);
k2 = k1/Yss;
k3 = k1*(Xss/Yss);


x_ini = [hare(1), lynx(1)]; %reset the initial values
p = [k1 k2 k3]; %reset params to hold the ecologically valid k1 value
% Simulation: use the model_PP function to run the simulation
[t,x]=ode45(@model_PP,t_range,x_ini,options,p);
% Plot time series of all variables
figure(26);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level');
title('Model behaviour: ecologically valid parameters');

%% Try to run model again with better parameter estimates
%model seems to run with too mnay oscillations, period too short. 
%consider k1 again - this is not birth rate but population growth including
%death rates. 
%lower k1 to be a better estimation of popultion growth. 
k1 = 0.5; 
Xss = mean(hare);
Yss = mean(lynx);
k2 = k1/Yss;
k3 = k1*(Xss/Yss);

x_ini = [hare(1), lynx(1)]; %reset the initial values
p = [k1 k2 k3]; %reset params to hold the ecologically valid k1 value
% Simulation: use the model_PP function to run the simulation
t_range = [1900:1920];
[t,x]=ode45(@model_PP,t_range,x_ini,options,p);
% Plot time series of all variables
figure(27);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level');
title('Model behaviour: refitted ecologically valid parameters');

% Plot time series of all variables
figure(28);
plot(t,x(:,:));
xlabel('Time ');
ylabel('Population Level');
title('Model behaviourand histroical data');
hold on 
plot(year, lynx, 'Color', 'black')
hold on
plot(year, hare, 'Color', 'black')

%% Model extensions
function exitflag = lynx_hare

%Enter the data and the initial guess.
td = t;
p = [30 4 0.4 0.018 0.8 0.023];

%Finally we use the fminsearch routine as follows:

[p,fval,exitflag] = fminsearch(@leastcomp,p,[],td,hare,lynx);
p
fval
end
function J = leastcomp(p,tdata,xdata,ydata)
%Create the least squares error function to be minimized.
n1 = length(tdata);
[t,y] = ode23(@lotvol,tdata,[p(1),p(2)],[],p(3),p(4),p(5),p(6));
errx = y(:,1)-xdata(1:n1)';
erry = y(:,2)-ydata(1:n1)';
J = errx'*errx + erry'*erry;
end

function dydt = lotvol(t,y,a1,a2,b1,b2)
tmp1 = a1*y(1) - a2*y(1)*y(2);
tmp2 = -b1*y(2) + b2*y(1)*y(2);
dydt = [tmp1; tmp2];
end

%% Local function - the model fitting function for the ode45 solver
function dxdt= model_PP(t,x,p)
% Initialize model vector with zeroes,
dxdt=zeros(2,1);
% Load parameter values from vector p
k1 = p(1);  %Rate of prey bith
k2 = p(2);  %Rate of predation
k3 = p(3); %Rate of predator death

X = x(1); %initial prey population
Y = x(2); %inital predator population 

% Differential Equations:
% change in X
dxdt(1)=  k1*X - k2*X*Y;
% change in Y
dxdt(2)=  k2*X*Y - k3*Y;

end