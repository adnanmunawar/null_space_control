function dydt = arm_model(t,y,Tau)
global Me Ce Ge Je
global ddq_pre
global n_joints
global q1 q2 q3 dq1 dq2 dq3
global l1 l2 l3

dy_idx_low = n_joints + 1;
dy_idx_high = 2*n_joints;
if n_joints == 2
    q1 = y(1); q2 = y(2);
    dq1 = y(3); dq2 = y(4);
elseif n_joints == 3
    q1 = y(1); q2 = y(2); q3 = y(3);
    dq1 = y(4); dq2 = y(5); dq3 = y(6);
end
M = eval(Me);
C = eval(Ce);
G = eval(Ge);
Fe = [-0.01,0.]';
J = eval(Je);
Tau = J.' * Fe;
dydt = zeros(n_joints,1);
dydt(1:n_joints) = y(dy_idx_low:dy_idx_high);
% dydt(dy_idx_low:dy_idx_high) = M\(Tau - C - G);
dydt(dy_idx_low:dy_idx_high) = Tau
end
