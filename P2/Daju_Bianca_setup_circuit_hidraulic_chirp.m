%%
% Nume si prenume: Daju Bianca Teodora
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 5; 
n = 2; 

%% Process data (fixed, do not modify)
c1 = (1000+n*300)/10000;
c2 = (1.15+2*(m+n/10)/20);
a1 = 2*c2*c1;
a2 = c1;
b0 = (1.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [2*(m/2+rand(1)*m/5); m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 10*c2/c1/1e4*1.5; % fundamental step size
Tfin = 30*c2/c1*10; % simulation duration

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
wf=1/18.54;
fmin=wf/2/pi/10;
fmax=wf/2/pi*10;
Ain=1.4


%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_hidraulic_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

plot(t,u,t,y)
shg

%% System identification

yst=(2.34+1.94)/2
ust=1.45 %l am determinat cu u_op_region
K=yst/ust+0.02

 %%
 %Gasirea fazei de -90, determinarea parametrillor zeta(factor de amortizare)
 % si wn(pulsatia
 %naturala)
 w1=pi/(1011.57-993.42);
 deltaT1=(1002.92-993.42); %se ia un maxim de la rosu si unul de la albastru pe x
 phi1=-rad2deg(w1*deltaT1); %trb sa fie aprox -90

 w2=pi/(1045.68-1028.99)
 deltaT2=1037.54-1028.99
 phi2=-rad2deg(w2*deltaT2) %valoarea cea mai buna

 w3=pi/(1079.47-1063.47);
 deltaT3=(1071.87-1063.47);
 phi3=-rad2deg(w3*deltaT3); 

 w4=pi/(1045.27-1028.68);
 deltaT4=(1037.68-1028.68); 
 phi4=-rad2deg(w4*deltaT4);


 Ay=(2.69-1.65)/2
 Au=1.4    %din proiectare Ain
 Im=-Ay/Au
 wn=w2      
 zeta=-K/2/Im  %trebuie sa aiba valoare peste 1, IMPORTANT , sa nu fie nici apropae de 1

 
 %% validare 

H=tf([K*wn^2],[1,2*zeta*wn,wn^2])
zpk(H)
% 
% T1 si T2 ,se pune in command window zpk(H) si voi lua factorii de
% %jos supra 1 
% 

T1=1/0.05
T2=1/0.7
% 
A=[0,1;-1/T1/T2,-(1/T1+1/T2)];
B=[0;K/T1/T2];
C=[1,0];
D=0; 
sys=ss(A,B,C,D); 
ysim2=lsim(sys,u,t,[y(1),2.2]); %iesirea simulata a sistemului , conditiile initiale


figure
plot(t,u,t,y,t,ysim2);

%J=1/sqrt(lenght(t)*norm(y-ysim2));
eMPN = norm(y-ysim2)/norm(y-mean(y))*100