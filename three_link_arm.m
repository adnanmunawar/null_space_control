%% Homework 4 Solution
% Modified by Adnan Munawar


%% Problem 1: Position Kinematics
% Calculate the forward position kinematics of the tip position (x3, y3) with respect to
% the base F0. That is, define the 4x4 homogeneous transformation matrix for the base
% to the tip, T03. You are not required to use D-H parameters, but you must be clear in
% showing all the intermediate steps

%%
%
clc
clear all
global l1 l2 l3
global m1 m2 m3 ml
syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
syms l1 l2 l3 m1 m2 m3 ml I1 I2 I3 g real;

% l1 = 1; l2=1; l3=1; m1=0.5; m2=0.5; m3=0.5; g=9.81;
%%
%Define transformation matrices
T1 = [cos(q1)    -sin(q1)    0    l1*cos(q1);
      sin(q1)    cos(q1)     0    l1*sin(q1);
      0          0           1    0;
      0          0           0    1];
  
T2 = [cos(q2)    -sin(q2)    0    l2*cos(q2);
      sin(q2)    cos(q2)     0    l2*sin(q2);
      0          0           1    0;
      0          0           0    1];
  
T3 = [cos(q3)    -sin(q3)    0    l3*cos(q3);
      sin(q3)    cos(q3)     0    l3*sin(q3);
      0          0           1    0;
      0          0           0    1];

%Generate transformations from base to tip:
T01 = T1;
T02 = T1*T2;
T03 = T1*T2*T3;


%% Problem 2: Velocity Kinematics
% Calculate the forward velocity kinematics and write out in matrix form, clearly
% showing the full 6-DOF (6x3) Jacobian equation describing the 3 translational and 3
% angular velocities of the tip at the end of the last link as a function the 3 joint
% velocities.

%%
% Position Vector to Tip (Frame 3) is Base (Frame 0)
P03 = [T03(1,4);
       T03(2,4);
       T03(3,4)];

%%
% Jacobian
Jv = simplify([diff(P03(1),q1), diff(P03(1),q2), diff(P03(1),q3);
                diff(P03(2),q1), diff(P03(2),q2), diff(P03(2),q3);
                diff(P03(3),q1), diff(P03(3),q2), diff(P03(3),q3)]);

%%
% Alternately, Jv = simplify(jacobian(P03,[dq1;dq2;dq3]))
%
% Append Jw to get the full jacobian
Jw = [0 0 0 ; 0 0 0 ; 1 1 1];
J = [Jv ; Jw];

%% Problem 5: Lagrangian Dynamics
% For the configuration described above, symbolically derive the kinetic and potential
% energy of the arm by modeling the links as point masses at the center of each of the
% links, with an additional load mass at the end of the last link (4 masses total).
% Then solve for the Lagrangian for the 3-link arm shown. This should have the number
% provided above for the robot parameters, and leave the joint parameters as symbolic
% values.
%
% Derive the dynamics of the arm using the Euler-Lagrange approach. Arrange the
% dynamics into the matrix-vector standard form for representing multi-link arm
% dynamics, labeling key components of the equation. Be sure to clearly identify the
% equation and do your best to factor into the standard form discussed in class.

%%
% To solve for the lagrangian dynamics, we need the kinematic and potential
% energies of all these link masses:

%%
% Jacobian of link 1, 2 and 3
J_l1 = jacobian(T01(1:3,4),[q1]);
J_l2 = jacobian(T02(1:3,4),[q1 q2]);
J_l3 = jacobian(T03(1:3,4),[q1 q2 q3]);

J_l1 = simplify(J_l1);
J_l2 = simplify(J_l2);
J_l3 = simplify(J_l3);


%%
% Next, we calculate the velocity at each mass w.r.t base frame

v_m1 = subs(J_l1,l1,l1/2) * dq1;
v_m2 = subs(J_l2,[l1 l2],[l1 l2/2]) * [dq1 ; dq2];
v_m3 = subs(J_l3,[l1 l2 l3],[l1 l2 l3/2]) * [dq1 ; dq2 ; dq3];
v_ml = subs(J_l3,[l1 l2 l3],[l1 l2 l3]) * [dq1 ; dq2 ; dq3];

v_m1 = simplify(v_m1);
v_m2 = simplify(v_m2);
v_m3 = simplify(v_m3);
v_ml = simplify(v_ml);

%%
% Then we compute the Kinematic Energy of each mass
K1 = 0.5 * m1 * (v_m1.' * v_m1) + 0.5 * I1 * dq1^2;
K2 = 0.5 * m2 * (v_m2.' * v_m2) + 0.5 * I2 * (dq1 + dq2)^2;
K3 = 0.5 * m3 * (v_m3.' * v_m3) + 0.5 * I3 * (dq1 + dq2 + dq3)^2;
Kl = 0.5 * ml * (v_ml.' * v_ml);

K1 = simplify(K1);
K2 = simplify(K2);
K3 = simplify(K3);
Kl = simplify(Kl);

g = 9.8;

%%
% Also, we need the potential energy of each mass
P1 = m1 * g * subs(T01(2,4),l1,l1/2);
P2 = m2 * g * subs(T02(2,4),l2,l2/2);
P3 = m3 * g * subs(T03(2,4),l3,l3/2);
Pl = ml * g * T03(2,4);

P1 = simplify(P1);
P2 = simplify(P2);
P3 = simplify(P3);
Pl = simplify(Pl);

%%
% Summing up all the Kinematic energies and Potential energies

K = K1 + K2 + K3;
P = P1 + P2 + P3;

%%
% Computing the lagrangian
m1 = 3; m2 = 2; m3 = 1; ml = 0;
I1 = 100; I2 = 100; I3 = 100;
l1 = 1; l2 = 1; l3 = 1;
g = 9.8;
L = simplify( K - P );
L = subs(L);

%%
% This is the equation for Euler Lagrange Dynamics
% $$\tau_i = \frac{d}{dt}\frac{\partial L}{\partial \dot{q}_i} - \frac{\partial L}{\partial q_i}$ 
% where
% $$i \in (1,2,3)$

%%
% Treating $$\tau_i =
% \frac{d}{dt}\frac{\partial L}{\partial \dot{q}_i} - \frac{\partial
% L}{\partial q_i}$ as $$\tau_i =
% \frac{d}{dt}A - B$

syms th1(t) th2(t) th3(t);

%%
% First differentiate A1 wrt to $$\theta_i$
A1 = diff(L,dq1);
%%
% The replace variables of q1, q2, q3 and their derivatives with the
% time-dependant counter parts. This is done so that we can differentiate
% wrt time in the next step
A1t = subs(A1,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
%%
% Differentiate w.r.t time
dA1t = diff(A1t, t);
%%
% Replace the time-dependant varialbes with non-time dependant varialbes to
% make the output look cleaner
A1 = subs(dA1t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t) ... 
    diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
%%
% Compute Potential energy for first link
B1 = diff(L,q1);

%%
% Repeat for all the links

A2 = diff(L,dq2);
A2t = subs(A2,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
dA2t = diff(A2t, t);
A2 = subs(dA2t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t) ...
    diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
B2 = diff(L,q2);

A3 = diff(L,dq3);
A3t = subs(A3,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
dA3t = diff(A3t, t);
A3 = subs(dA3t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t) ...
    diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
B3 = diff(L,q3);

Tau_l1 = A1 - B1;
Tau_l2 = A2 - B2;
Tau_l3 = A3 - B3;

%%
% Now we want to split the dynamics into the M, C and G matrices.
M11 = simplify(Tau_l1 - subs(Tau_l1,ddq1,0)) /ddq1;
M12 = simplify(Tau_l1 - subs(Tau_l1,ddq2,0)) /ddq2;
M13 = simplify(Tau_l1 - subs(Tau_l1,ddq3,0)) /ddq3;

M21 = simplify(Tau_l2 - subs(Tau_l2,ddq1,0)) /ddq1;
M22 = simplify(Tau_l2 - subs(Tau_l2,ddq2,0)) /ddq2;
M23 = simplify(Tau_l2 - subs(Tau_l2,ddq3,0)) /ddq3;

M31 = simplify(Tau_l3 - subs(Tau_l3,ddq1,0)) /ddq1;
M32 = simplify(Tau_l3 - subs(Tau_l3,ddq2,0)) /ddq2;
M33 = simplify(Tau_l3 - subs(Tau_l3,ddq3,0)) /ddq3;

G11 = subs(Tau_l1, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
G21 = subs(Tau_l2, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
G31 = subs(Tau_l3, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);

%%
% The intertia matrix is as follows

M = [M11 M12 M13;
     M21 M22 M23;
     M31 M32 M33];
 
 
%%
% This is the gravity vector

G = [G11;
     G21;
     G31];

C11 = Tau_l1 - (M(1,:) * [ddq1 ddq2 ddq3].' + G11);
C21 = Tau_l2 - (M(2,:) * [ddq1 ddq2 ddq3].' + G21);
C31 = Tau_l3 - (M(3,:) * [ddq1 ddq2 ddq3].' + G31);

%%
% This is the Coriolis Vector
C = [C11;
     C21;
     C31];


