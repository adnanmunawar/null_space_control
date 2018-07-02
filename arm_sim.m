function arm_sim(M, C, G, J_ee)
syms q1 q2 q3
syms dq1 dq2 dq3
syms ddq1 ddq2 ddq3
syms t1 t2 t3
global l1 l2 l3
global x_d
global e de

 J_ee = subs(J_ee);
 Me = simplify(M);
 Ce = simplify(C);
 Ge = simplify(G);
 
 q = [0,0,0]';
 dq = [0,0,0]';
 ddq = [0,0,0]';
 T_ext = [t1, t2, t3]';
 dt = 0.005;
 pos = zeros(3,2000);
 e = [0,0,0]';
 de = [0,0,0]';
 RHS = simplify(Me\(T_ext - Ce - Ge));
 x_d = [-3.0,0.0,3.18]';
 x_c = [0,0,0]';
 Kp = diag([3,3,1]);
 Kd = diag([3,3,1]);
 for i=1:2000
    arm_plot(q(1), q(2), q(3))
    plot(x_d(1), x_d(2), '-bX', 'MarkerSize', 10);
    dq = dq + (ddq * dt);
    q = q + (dq * dt);
    x_c_pre = x_c;
    x_c = get_ee_pos(q(1), q(2), q(3));
    e = (x_d - x_c);
    de = (x_c - x_c_pre)*(1/dt);
    x_ee = Kp*e - Kd*de;
    pos(:,i) = q;
    q1   = q(1); q2   = q(2); q3   = q(3);
    dq1  = dq(1); dq2  = dq(2); dq3  = dq(3);
    ddq1 = ddq(1); ddq2 = ddq(2); ddq3 = ddq(3);
     Jee_numeric = eval(J_ee);
     M_ee_inv = Jee_numeric / eval(Me) * Jee_numeric.';
     if (abs(det(Jee_numeric * Jee_numeric.'))) > 0.005
         M_ee = double(inv(M_ee_inv));
     else
         [U,S,V] = svd(M_ee_inv);
         for j=1:size(S,1)
             if S(j,j) < 0.005
                 S(j,j) = 0;
             else
                 S(j,j) = 1/S(j,j);
             end
         end
         M_ee = V * S * U.';
     end
     T_ee = double(Jee_numeric.' * M_ee * x_ee) + eval(Ge);
     t1 = T_ee(1); t2 = T_ee(2); t3 = T_ee(3);
     ddq= double(eval(RHS));
 end
 for i=1:size(pos,2)
     arm_plot(pos(1,i), pos(2,i), pos(3,i));
 end
end

function arm_plot(q1, q2, q3)
global l1 l2 l3
global x_d
global e de
global hText
P1 = [l1*cos(q1), l1*sin(q1), 0];
P2 = [l2*cos(q1 + q2) + l1*cos(q1) ,l2*sin(q1 + q2) + l1*sin(q1), 0];
P3 = [l3*cos(q1 + q2 + q3) + l2*cos(q1 + q2) + l1*cos(q1), l3*sin(q1 + q2 + q3) + l2*sin(q1 + q2) + l1*sin(q1), 0];
x = [0 P1(1) P2(1) P3(1)];
y = [0 P1(2) P2(2) P3(2)];
delete(hText);
hold on
plot(x, y, '-ro','LineWidth',3,'MarkerSize',10, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'blue')
plot(x_d(1), x_d(2), '-rX')
hText = text(3,3, num2str(e));
hold off
xlim([-4,4]);
ylim([-4,4]);
pause(0.001)
end

function [P_ee] = get_ee_pos(q1,q2,q3)
global l1 l2 l3
P_ee = [l3*cos(q1 + q2 + q3) + l2*cos(q1 + q2) + l1*cos(q1), l3*sin(q1 + q2 + q3) + l2*sin(q1 + q2) + l1*sin(q1), q1+q2+q3].';
end