function arm_sim(M, C, G, J, v_filter)
syms q1 q2 q3
syms dq1 dq2 dq3
syms ddq1 ddq2 ddq3
syms t1 t2 t3
global l1 l2 l3
global x_d
global e de
global hAxes dAxes peAxes hText uText tText nText feAxes eeAxes
global fig
global Te T_null ctime Fe
global n_joints

 n_joints = size(M,1);
 J = eval(J(1:2,:));
 M = simplify(M);
 C = simplify(C);
 G = simplify(G);
 B = ones(diag(n_joints));
 
 q = zeros(n_joints, 1);
 dq = zeros(n_joints, 1);
 ddq = zeros(n_joints, 1);
 x_d = [-1.9, 1.5]';
 
 dt = 0.01;
 q_list = zeros(n_joints,2000);
 Kp = 10*diag([1,1]);
 Kd = sqrt(Kp);
 if n_joints == 2
     null_gains = ones(2,1);
     x_d_null = [0,0]';
     q = [1.0,0.3]';
 elseif n_joints == 3
     null_gains = ones(3,1);
     x_d_null = [0,0,0]';
     q = [0.6,0.4,-0.6]';
 end
 x_e = get_ee_pos(q);
 Kp_null = diag(null_gains);
 Kd_null = sqrt(Kp_null);
 fig = figure;
 grid on
 hold on
 hAxes = plot([0,0,1,2], [0,0,0,0], '-ro','LineWidth',3,'MarkerSize',10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'blue');
 dAxes = plot(0,0, '-bX','MarkerSize',10);
 peAxes = plot(x_e(1),x_e(2), '-m.');
 feAxes = quiver(x_e(1), x_e(2),0.1,0.1, '-g', 'LineWidth',3, 'MaxHeadSize', 0.9);
 eeAxes = quiver(x_e(1), x_e(2),0.1,0.1, '-r', 'LineWidth',3, 'MaxHeadSize', 0.9);
 hText = text(-3,3,'Error');
 uText = text(-3,2.5,'Input');
 nText = text(-3,2.0, 'Tnull');
 tText = text(-3,1.5,'t');
 qText = text(-3,-1.0,'t');
 xNullText = text(-3,-1.5,'t');
 eNullText = text(-3,-2.0,'t');
 xlim([-4,4]);
 ylim([-4,4]);
 ctime = 0;
 e = (x_d - x_e);
 e_null = Kp_null * (x_d_null - q) - Kd_null * dq;
 set(gca,'ButtonDownFcn', @click_callback);
 i=1;
 while norm(e_null) > 0.01 || norm(e) > 0.1
    q_list(:,i) = q;
    arm_plot(q)
    dq = ddq * dt;
    q = q + (dq * dt) + 0.5 * ddq * dt * dt;
    x_e_pre = x_e;
    x_e = get_ee_pos(q);
    e = (x_d - x_e);
    de = (x_e - x_e_pre) / dt;
    V_max = .3;
    if v_filter == false
        x_ee = Kp * e - Kd * de;
    else
        x_ee = velocity_filter(Kp, Kd, x_e, de, x_d, V_max);
    end
    e_null = Kp_null * wrapToPi(x_d_null - q) - Kd_null * dq;
    if n_joints == 2
        q1 = q(1); q2 = q(2);
        dq1 = dq(1); dq2 = dq(2);
    elseif n_joints == 3
        q1 = q(1); q2 = q(2); q3 = q(3);
        dq1 = dq(1); dq2 = dq(2); dq3 = dq(3);
    end
    Jee_numeric = eval(J);
    M_numeric = eval(M);
    C_numeric = eval(C);
    G_numeric = eval(G);
    M_ee_inv = Jee_numeric / M_numeric * Jee_numeric.';
    if (abs(det(Jee_numeric * Jee_numeric.'))) > 0.001
        M_ee = double(inv(M_ee_inv));
    else
        [U,S,V] = svd(M_ee_inv);
        for j=1:size(S,1)
            if S(j,j) < 0.00025
                S(j,j) = 0;
            else
                S(j,j) = 1/S(j,j);
            end
        end
        M_ee = V.' * S * U.';
    end
    xJee_numeric_inv = double(M_ee * Jee_numeric / M_numeric);
    Fe = M_ee * x_ee;
    Te = double(Jee_numeric.' * Fe);
    T_null = (eye(n_joints) - Jee_numeric.' * xJee_numeric_inv) * e_null;
    T_ee = (Te + G_numeric + T_null);
    ddq = double(M_numeric \ (T_ee - C_numeric - G_numeric - 0*B*dq));
    ctime = ctime + dt;
    i = i+1;
    peAxes.XData = [peAxes.XData, x_e(1)];
    peAxes.YData = [peAxes.YData, x_e(2)];
    fe = Fe/norm(Fe);
    temp_x_ee = e/norm(e);
    feAxes.XData = x_e(1);
    feAxes.YData = x_e(2);
    feAxes.UData = fe(1);
    feAxes.VData = fe(2);
    eeAxes.XData = x_e(1);
    eeAxes.YData = x_e(2);
    eeAxes.UData = temp_x_ee(1);
    eeAxes.VData = temp_x_ee(2);
    qText.String = strcat('q= [', num2str(q'), ']');
    xNullText.String = strcat('x_{null}= [', num2str(x_d_null'), ']');
    eNullText.String = strcat('e_{null}= [', num2str(wrapToPi(x_d_null - q)'), ']');
 end
 lim = i;
 for j=1:10
    for i=1:lim
        arm_plot(q_list(:,i));
        pause(0.01)
    end
 end
 
 
end

function arm_plot(q)
global l1 l2 l3
global x_d
global e de
global hAxes dAxes
global hText uText tText nText
global Te T_null ctime Fe
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
dAxes.XData = x_d(1);
dAxes.YData = x_d(2);
hText.String = strcat('e= [', num2str(e'), ']');
uText.String = strcat('Tx= [', num2str(Te'), ']');
nText.String = strcat('Tn= [', num2str(T_null'), ']');
tText.String = strcat('t= ', num2str(ctime));
pause(0.0001)
end

function [Ux] = velocity_filter(Kp, Kd, x, dx, tx, Vmax)
Kp = diag(Kp);
Kd = diag(Kd);
lamb = Kp ./ Kd;
x_tilde = x - tx;
sat = Vmax ./ (lamb .* abs(x_tilde));
scale = ones(2,1);
if any(sat < 1)
    [val, index] = min(sat);
    unclipped = Kp .* x_tilde(index);
    clipped = Kd .* Vmax .* sign(x_tilde(index));
    scale = ones(2,1) .* sat(index);
    scale(index) = 1;
end
clipped_sat = sat ./ scale;
clipped_sat(clipped_sat > 1) = 1;
Ux = -Kd .* (dx + clipped_sat .* scale .* lamb .* x_tilde);
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

function click_callback(gcbo, axes)
global x_d
% global fig
% cP = fig.CurrentPoint;
cP = get(gcbo, 'CurrentPoint');
x_d(1) = cP(2,1);
x_d(2) = cP(2,2);
end