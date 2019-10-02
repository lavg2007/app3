%%
clc
clear all
close all
% constante
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

%% mise en matrice des quatres fonction de transfert
A = [0 1 0 0;
     0 -((B_moteur/N) + N*B_load)/((J_moteur/N) + N*J_load) K_i/((J_moteur/N) + N*J_load) 0;
     0 -K_b/(La*N) -Ra/La 1/La;
     0 0 0 -1/T];
B = [0;
     0;
     0;
     (K/T)];
% la sortie est directement theta de la load
C = [1 0 0 0];
D = 0;

sys = ss(A, B, C, D);
% fonction de transfert du systeme en boucle ouverte (methode numerique)
H = tf(sys)

% reponse impulsionnel de la fonciton de transfert
figure()
impulse(H)
%% dans le cas du shaft mou
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
 
H_mou = tf(sys)
time = [0:0.001:5] % plus longue periode sur le graph

figure
impulse(H_mou,time)
figure
impulse(H-H_mou,time)


%% separation en 4 x 1 ordre
[B,A] = tfdata(H,'v')
[r,p,k] = residue(B,A) % separation en systeme d'ordre 1
integrator = tf(1,[1 0])
 
% creation des 4 fonctions de transfert
h1 = tf(r(1),[1 -p(1)])
h2 = tf(r(2),[1 -p(2)])
h3 = tf(r(3),[1 -p(3)])
h4 = tf(r(4),[1 -p(4)])
% assemblage de la fonction de transfert total
G = h1*h2*h3*h4
% calcul de la difference en gain
dcG = dcgain(G/integrator)
dcH = dcgain(H/integrator)
GAIN = dcH/dcG

figure()
hold on
% impulse(h1)
% impulse(h2)
% impulse(h3)
% impulse(h4)
[b1,a1] = residue([r(1) r(3) r(4)],[p(1) p(3) p(4)],[]); % ordre 3
impulse(tf(b1,a1));
[b2,a2] = residue([r(3) r(4)], [p(3) p(4)],[]);% ordre 2
impulse(tf(b2,a2));
% [b3,a3] = residue([r(4)], [p(4)],[]); % ordre 1
% impulse(tf(b3,a3));
impulse(H); % original
impulse(G*GAIN); % recompose avec gain ajuste
legend('3','2','H','G')
 
%% reponse a la consigne
% l'equation de la reponse est
rep_tet_d = Kp*H/(Kp*H + 1)
rep_tet_d_norm = minreal(rep_tet_d)
% sous forme complete

figure()
step(rep_tet_d) 
 
%% identification du moteur
load('donnees_moteur_2016.mat')
x = [-diff(vitesse)./diff(t),tension(1:end-1)]
rep = pinv(x)*vitesse(1:end-1)
T_calc = rep(1)
K_calc = rep(2)

h = tf(K_calc,[T_calc 1])
step(h)

E = 8;
tm = 0.52;
Ia = 1.09;
% trouver les caractristique du moteur a rotor bloque
ra = E/Ia;
ki = tm/Ia;
kb = ki;

% trouver les caractristique du moteur en mouvement
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
% reconstruction de  la boucle ouverte complete
M  = G0*G1*G2*G4/(1-G1*G2*G3)
m  = minreal(M) % normalisation
% Reduction d'ordre
M_ord3  = G0*G2*G4/(1-G2*G3)
M_ord2  = G2*G4/(1-G2*G3)

figure()
impulse(M)
hold on
impulse(M_ord3)
impulse(M_ord2*80)
legend('M','ord3','ord2')
