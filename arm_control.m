function arm_control(M,C,G,J, t0, tf)
syms q1 q2 q3
syms dq1 dq2 dq3
syms ddq1 ddq2 ddq3
global Me Ce Ge Je
global hAxes peAxes
global n_joints
global l1 l2 l3

n_joints = size(M,1);
Je = simplify(J(1:2,:));
Me = simplify(M);
Ce = simplify(C);
Ge = simplify(G);
Tau = zeros(n_joints, 1);
q0 = zeros(n_joints, 1);
dq0= zeros(n_joints, 1);
q0 = [-0.4,0.8]';
[t,y] = ode45(@(t,y) arm_model(t,y,Tau), [t0,tf], [q0,dq0]);
figure;
hold on;
grid on;
xlim([-4,4]);
ylim([-4,4]);
x_e = get_ee_pos(q0);
hAxes = plot([0], [0], '-ro','LineWidth',3,'MarkerSize',10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'blue');
peAxes = plot(x_e(1),x_e(2), '-m.');
sleep_time = 2*(tf - t0)/length(t);
for i=1:size(y,1)
   x_e = get_ee_pos(y(i,:));
   arm_plot(y(i,:));
   peAxes.XData = [peAxes.XData, x_e(1)];
   peAxes.YData = [peAxes.YData, x_e(2)];
   pause(sleep_time);
end
end

function [x_ee] = get_ee_pos(q)
global l1 l2 l3
global n_joints
if n_joints == 2
    q1 = q(1); q2 = q(2);
    x_ee = [l2*cos(q1 + q2) + l1*cos(q1), l2*sin(q1 + q2) + l1*sin(q1)].';
elseif n_joints == 3
    q1 = q(1); q2 = q(2); q3 = q(3);
    x_ee = [l3*cos(q1 + q2 + q3) + l2*cos(q1 + q2) + l1*cos(q1), l3*sin(q1 + q2 + q3) + l2*sin(q1 + q2) + l1*sin(q1)].';
end
end

function arm_plot(q)
global l1 l2 l3
global hAxes
global n_joints

if n_joints == 2
    q1 = q(1); q2 = q(2);
    P1 = [l1*cos(q1), l1*sin(q1), 0];
    P2 = [l2*cos(q1 + q2) + l1*cos(q1) ,l2*sin(q1 + q2) + l1*sin(q1), 0];
    x = [0 P1(1) P2(1)];
    y = [0 P1(2) P2(2)];
elseif n_joints == 3
    q1 = q(1); q2 = q(2); q3 = q(3);
    P1 = [l1*cos(q1), l1*sin(q1), 0];
    P2 = [l2*cos(q1 + q2) + l1*cos(q1) ,l2*sin(q1 + q2) + l1*sin(q1), 0];
    P3 = [l3*cos(q1 + q2 + q3) + l2*cos(q1 + q2) + l1*cos(q1), l3*sin(q1 + q2 + q3) + l2*sin(q1 + q2) + l1*sin(q1), 0];
    x = [0 P1(1) P2(1) P3(1)];
    y = [0 P1(2) P2(2) P3(2)];
end
hAxes.XData = x;
hAxes.YData = y;
end
