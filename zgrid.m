clear all; clc;
N = 100;
z_star = 0.6+0.37i;
sampling_freq = 20;
T = 1/sampling_freq;
G=tf([1.4432],[1 -0.2299],T);
%D_ang =pi+atan(imag(z_star)/(real(z_star)-0.229));
G_ang =angle(evalfr(G,z_star));
d_ang = pi - G_ang;
d_ang2 = d_ang - pi;
alpha_ang=linspace(d_ang2+10*pi/180, pi-10*pi/180,N);
beta1 = -1;
d_ang1 = angle(z_star+beta1);
alpha =zeros(N,1);
beta =zeros(N,1);
beta_ang =zeros(N,1);
myPoles =[];
for k = 1:N
    alpha(k) = imag(z_star)/tan(alpha_ang(k))-real(z_star);
    beta_ang(k) = alpha_ang(k)-d_ang2-d_ang1;
    beta1(k) = imag(z_star)/tan(beta_ang(k))-real(z_star);
    %beta2(k)=imag(z_star)/tan(alpha_ang(k))-real(z_star);
    D_z(k) = tf([1 alpha(k)], [1 beta1(k)-1 -beta1(k)], T);
    K(k) = 1/abs(evalfr(G*D_z(k), z_star));
    D(k) = K(k)*D_z(k);
    %controller_G = minreal(G*D_z/(1+G*D_z));
    controller_G2 = feedback(G*D(k),1);
    myPoles(k,:) = pole(controller_G2);
    L = tf([1 -1],[1 0],T);
    Kv = evalfr(minreal(L*G*D),1)/T
end
figure
plot(alpha, real(myPoles), 'kx', 'markersize',3, 'linewi',2)
grid on

Tfinal = 2;
t = (0:T:Tfinal);
stepReference = ones(length(t),1)*2*pi;
rampReference = t;
M = 5;
alpha_restrict = linspace(-0.8,-0.1,M);
mystepResponse = zeros(length(t),M);
myrampResponse = zeros(length(t),M);
mystepControl = zeros(length(t),M);
myrampcontrol = zeros(length(t),M);
clear alpha_ang beta_ang beta K D myPoles Kv D_z Kv
for i = 1:M
    alpha_ang(i) = angle(alpha_restrict(i)+z_star);
    beta_ang = alpha_ang(i) - d_ang2 - d_ang1;
    beta1(i) = imag(z_star)./tan(beta_ang)-real(z_star);
    D_z = tf([1 alpha_restrict(i)], [1 beta1(i)-1 -beta1(i)], T);
    K(i) = 1/abs(evalfr(G*D_z, z_star));
    D = K(i)*D_z;
    controller_G2 = feedback(G*D,1);
    myPoles(k,:) = pole(controller_G2);
    L = tf([1 -1],[1 0],T);
    Kv(i) = evalfr(minreal(L*G*D),1)/T;
    mystepResponse(:,i) = lsim(controller_G2,stepReference,t);
    myrampResponse(:,i) = lsim(controller_G2,rampReference,t);
    G_u = minreal(D/(1+G*D));
    myStepControl(:,i) = lsim(G_u, stepReference, t);
    myrampControl(:,i) = lsim(G_u, rampReference, t);
end

%Kv = rampReference-myrampResponse;
figure
plot(t,mystepResponse)
grid on
hold on
plot(t,stepReference)
%stepinfo(t,mystepReference)
grid on
hold off
title('Step Response')
figure
plot(t,myrampResponse)
grid on
hold on
plot(t,rampReference)
grid on
hold off
title('Ramp Response')
figure
stairs(t,myrampControl)
grid on
title('Ramp Control')
figure
stairs(t,myStepControl)
grid on
title('Step Control')