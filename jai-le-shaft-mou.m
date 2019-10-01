%%
clc
clear all
close all
 
B_moteur = 0.01;
J_moteur = 0.02;
K_b = 0.5;
K_i = 0.5;
T = 0.01;
K = 100;
J_load = 1;
B_load = 1;
La = 0.008;
Ra = 8;
N = 0.1;
Kp = 0.318;
 
 
% x1 = tetaL, x2 = wL, x3 = ia, x4 = ea
% u = tetaD
A = [0 1 0 0;
     0 -((B_moteur/N) + N*B_load)/((J_moteur/N) + N*J_load) K_i/((J_moteur/N) + N*J_load) 0;
     0 -K_b/(La*N) -Ra/La 1/La;
     0 0 0 -1/T];
B = [0;
     0;
     0;
     (K/T)];
% y = tetaL
C = [1 0 0 0];
D = 0;
 
sys = ss(A, B, C, D);
 
H = tf(sys)
 
figure
impulse(H)
%% mou
Kl = 10000

A = [0                  1               0               0       0               0;
     -Kl/J_load         -B_load/J_load  0               0       N*Kl/J_load     0;
     0                  0               -Ra/La          1/La    0               -K_b/La;
     0                  0               0               -1/T    0               0;
     0                  0               0               0       0               1;
     N*Kl/(J_moteur)    0               K_i/J_moteur    0       -Kl*N^2/J_moteur    -B_moteur/J_moteur];
B = [0;
     0;
     0;
     (K/T);
     0;
     0];
% y = tetaL
C = [1 0 0 0 0 0];
D = 0;
 
sys = ss(A, B, C, D);
 
H = tf(sys)
 
figure
impulse(H)


%% separation en 4 x 1 ordre
[B,A] = tfdata(H,'v')
[r,p,k] = residue(B,A)
integrator = tf(1,[1 0])
 
 
h1 = tf(r(1),[1 -p(1)])
h2 = tf(r(2),[1 -p(2)])
h3 = tf(r(3),[1 -p(3)])
h4 = tf(r(4),[1 -p(4)])
G = h1*h2*h3*h4
dcG = dcgain(G/integrator)
dcH = dcgain(H/integrator)
GAIN = dcH/dcG
figure()
hold on
% impulse(h1)
% impulse(h2)
% impulse(h3)
% impulse(h4)
[b1,a1] = residue([r(1) r(3) r(4)],[p(1) p(3) p(4)],[]);
impulse(tf(b1,a1));
[b2,a2] = residue([r(3) r(4)], [p(3) p(4)],[]);
impulse(tf(b2,a2));
[b3,a3] = residue([r(4)], [p(4)],[]);
impulse(tf(b3,a3));
impulse(H);
impulse(G*GAIN);
legend('1','2','3','4','H','G')
 
%% reponse a la consigne
 
rep = Kp*H/(Kp*H + 1)
 
figure()
% impulse(rep)
step(rep)
 
 
%% identification du moteur
% x = [-diff(vitesse)./diff(t),tension(1:end-1)]
% rep = pinv(x)*vitesse(1:end-1)
% T_calc = rep(1)
% K_calc = rep(2)
% 
% h = tf(K_calc,[T_calc 1])
% step(h)
E = 8;
tm = 0.52;
Ia = 1.09;
 
ra = E/Ia;
ki = tm/Ia;
kb = ki;
 
syms jm bm;
eq1 = (ki/ra)/(bm + ki*kb/ra) == K_calc;
bm = double(solve(eq1))
eq2 = (jm)/(bm + ki*kb/ra) == T_calc;
jm = double(solve(eq2))
 
%% toute la fonction de transfert
 
G0 = tf(K/T,[1 1/T]);
G1 = tf(1,[La Ra]);
G2 = tf((K_i/((J_moteur/N) + N*J_load)),[1 ((B_moteur/N) + N*B_load)/((J_moteur/N) + N*J_load)]);
G3 = -1/N * K_b;
G4 = tf(1,[1 0]);
 
M  = G0*G1*G2*G4/(1-G1*G2*G3)
% Reduction d'ordre
M_ord3  = G0*G2*G4/(1-G2*G3)
M_ord2  = G2*G4/(1-G2*G3)
m  = minreal(M)
figure()
impulse(M)
hold on
impulse(M_ord3)
impulse(M_ord2*80)
% impulse(G1)
% impulse(G2)
% impulse(G4)
legend('M','ord3','ord2')
%% reponse en boucle ferme

