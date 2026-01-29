%%
% Nume si prenume: Daju Bianca Teodora
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 5; 
n = 2; 

%% Process data (fixed, do not modify)
%se calculeaza constantele fizice ale rezervoarelor
c1 = (1000+n*300)/10000;
c2 = (1.15+2*(m+n/10)/20);
a1 = 2*c2*c1;
a2 = c1;
b0 = (1.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [2*(m/2+rand(1)*m/5); m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 10*c2/c1/1e4*1.5; % fundamental step size
Tfin = 30*c2/c1;

gain = 10;
umin = 0; umax = gain; % input saturation
ymin = 0; ymax = b0*gain/1.5; % output saturation

whtn_pow_in = 1e-6*5*(((m-1)*8+n/2)/5)/2*6/8; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(10); % input quantizer (DAC)

whtn_pow_out = 1e-5*5*(((m-1)*25+n/2)/5)*6/80*(0.5+0.3*(m-2)); % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = (m/2+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
u0 = 0;        % fixed
ust = 5;     % must be modified (saturation)
t1 = 10*c2/c1; % recommended 

%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_hidraulic_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

plot(t,u,t,y)
shg

%% System identification

%Indecsii la declansarea treptei(t=0)
i1=4976;
i2=6611;
%Indecsi la stabilizare
i3=12136;
i4=13255;

%Calcularea valorilor medii
y_st=mean(y(i3:i4));
y0=mean(y(i1:i2));

u_st=mean(u(i3:i4));
u0=mean(u(i1:i2));

k=(y_st-y0)/(u_st-u0)

figure
plot(t,u,t,y)
hold on

%%
%deducem doua constante de timp 
i5=6920;
i6=8566;
t_reg=t(i5:i6);
y_reg=log(y_st-y(i5:i6));
figure
plot(t_reg,y_reg)

Areg=[sum(t_reg.^2),sum(t_reg);
        sum(t_reg),length(t_reg)];
breg=[sum(y_reg.*t_reg);sum(y_reg)];
       
theta=inv(Areg)*breg  %solutia sistemului de regresie
T1=-1/theta(1)
%%
i7=10670;
i8=11089;
Ti=t(i8)-t(i7)

T2vec=0.1:0.01:T1;
Y_ecuatie=T1*T2vec.*log(T2vec)-T2vec*(Ti+T1*log(T1))+T1*Ti;

plot(T2vec,Y_ecuatie)
grid on

T2=3.03

%%
%validare 
%modelul in spatiul starilor a fost construit pe baza paramtrilor
%identificati K,T1,T2
H=tf(k,[T1*T2,T1+T2,1])

ysim=lsim(H,u,t)

figure
plot(t,u,t,y,t,ysim)

A=[0,1;-1/T1/T2,-(1/T1+1/T2)];
B=[0;k/T1/T2];
C=[1,0];
D=0; 
sys=ss(A,B,C,D); 
ysim2=lsim(sys,u,t,[y(1),1.13]);% aici imi modifica 10 le 

figure
plot(t,u,t,y,t,ysim2);


hold on


%J=1/sqrt(lenght(t)*norm(y-ysim2));
eMPN = norm(y-ysim2)/norm(y-mean(y))