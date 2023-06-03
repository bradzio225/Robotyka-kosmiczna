%--------------------------------------------------------------------------
% Function fot calculating symbolic M and C matrices
%--------------------------------------------------------------------------
clc
clear all
close all
tic

for satelliteMassScalerIndex = 1:1:3

p1 = 0; p2 = 0;
m0 = [35;35*10;35*200];
m1 = 2.5; m2 = 1.5;
L1 = 0.6; L2 = 0.5;
a1 = L1/2; a2 = L2/2;
I0 = [1.5;1.5*10;1.5*200];
tk = 3.5;
q1=0; q2=0; dq1=0; dq2=0;
q_p = [0 0]; dq_p = [0 0];
q_k = [pi/8 pi/4]; dq_k = [0 0];
V = [pi/24, pi/12];
t_a = [0.37 0.4]; t_b = [0.55 0.65]; t_k = [3.5 3.5];
dt = tk/2004; time = 0:dt:(tk(1)-dt);

joint_1 = Trajectory_Generation(q_p(1),dq_p(1),q_k(1),dq_k(1),t_k(1),t_a(1),t_b(1),dt,V(1));
joint_2 = Trajectory_Generation(q_p(2),dq_p(2),q_k(2),dq_k(2),t_k(2),t_a(2),t_b(2),dt,V(2));

q_s = [joint_1.q; joint_2.q];
dq_s = [joint_1.dq; joint_2.dq];

%% Preallocation for optimization purposes

U_1 = zeros(1,length(joint_1.q)-1);
U_2 = zeros(1,length(joint_1.q)-1);
earth_position_x = zeros(1,length(U_1)-1);
earth_position_y = zeros(1,length(U_1)-1);
earth_orientation = zeros(1,length(U_1)-1);
earth_q_1 = zeros(1,length(U_1)-1);
earth_q_2 = zeros(1,length(U_1)-1);
earth_Vel_x = zeros(1,length(U_1)-1);
earth_Vel_y = zeros(1,length(U_1)-1);
earth_Omee = zeros(1,length(U_1)-1);
space_position_x = zeros(1,length(U_1)-1);
space_position_y = zeros(1,length(U_1)-1);
space_orientation = zeros(1,length(U_1)-1);
space_q_1 = zeros(1,length(U_1)-1);
space_q_2 = zeros(1,length(U_1)-1);
space_Vel_x = zeros(1,length(U_1)-1);
space_Vel_y = zeros(1,length(U_1)-1);
space_Omee = zeros(1,length(U_1)-1);

for i = 1:1:(length(joint_1.q) - 1)
A = [(1/3 * m1 + m2) * L1^2 + 1/3 * m2 * L2^2 + m2 * L1 * L2 * cos(joint_2.q(i)),... 
    1/3 * m2 * L2^2 + 1/2 * m2 * L1 * L2 * cos(joint_2.q(i));...
    1/3 * m2 * L2^2 + 1/2  * m2 * L1 * L2 * cos(joint_2.q(i)), 1/3 * m2 * L2^2];

B = [joint_1.ddq(i); joint_2.ddq(i)];

C = [0, -m2 * L1 * L2 * (joint_1.dq(i) + 1/2 * joint_2.dq(i)) * sin(joint_2.q(i));...
    1/2 * m2 * L1 * L2 * joint_1.dq(i) * sin(joint_2.q(i)), 0];

D = [joint_1.dq(i); joint_2.dq(i)];

    Out = A * B + C * D;

    U_1(i) = Out(1);
    U_2(i) = Out(2);
end

for k = 1:1:length(q_s)-1
    q_send(1) = q_s(1,k);
    q_send(2) = q_s(2,k);
    dq_send(1) = dq_s(1,k);
    dq_send(2) = dq_s(2,k);
    [earth_position_x(k), earth_position_y(k), earth_orientation(k), earth_q_1(k), earth_q_2(k), earth_Vel_x(k), earth_Vel_y(k),earth_Omee(k)]  = direct_2DoF(q_send, dq_send, L1, L2);
end

earthDirect = [earth_position_x;earth_position_y;earth_orientation;earth_q_1;earth_q_2;earth_Vel_x;earth_Vel_y;earth_Omee];

%% Graphs for satellite:
[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10] = satelite_inverse_dynamics(p1, p2, m0(satelliteMassScalerIndex), m1, m2, I0(satelliteMassScalerIndex), L1, L2, U_1, U_2, q1, q2, dq1, dq2, tk);

q = [x4;x5];
dq =[x9;x10];

for k = 1:1:length(q)
    q_send(1) = q(1,k);
    q_send(2) = q(2,k);
    dq_send(1) = dq(1,k);
    dq_send(2) = dq(2,k);
    [space_position_x(k), space_position_y(k), space_orientation(k), space_q_1(k), space_q_2(k), space_Vel_x(k), space_Vel_y(k),space_Omee(k)]  = direct_2DoF(q_send, dq_send, L1, L2);
end

figure('Name', 'Satellite')
subplot(2,1,1)
plot(time,x1, time, x2)
title("Position of the satellite")
xlabel('Time [s]')
ylabel('Position [m]')
legend('Position x_S', 'Position y_S')
grid on;
subplot(2,1,2)
plot(time, x3)
title("Orientation of the satellite")
xlabel('Time [s]')
ylabel('Orientation [rad]')
legend('Orientation  \theta_S')
grid on;

figure('Name', 'Manipulator on Satellite')
plot(time,x4, time, x5)
title("Angular position of the joints")
legend('Angular position of the first joint \theta_1', 'Angular position of the second joint \theta_2')
xlabel('Time [s]')
ylabel('Angular position [rad]')
grid on;

%% Graphs for the effector 
figure()
subplot(2,2,1)
plot(time, space_position_x+x1)
hold on
plot(time, space_position_y+x2)
grid on;
title("Position of the endeffector - SPACE - for m_0 = " + m0(satelliteMassScalerIndex)  + "kg and I_0 = " + I0(satelliteMassScalerIndex) + "kgm^2.")
legend('x_e_e _S_P_A_C_E', 'y_e_e _S_P_A_C_E')
xlabel('Time [s]');
ylabel('Position [m]');
hold off
subplot(2,2,2)
plot(time, space_orientation+x3)
grid on;
title("Orientation of the endeffector - SPACE - for m_0 = " + m0(satelliteMassScalerIndex)  + "kg and I_0 = " + I0(satelliteMassScalerIndex) + "kgm^2.")
legend('Ψ_e_e _S_P_A_C_E')
xlabel('Time [s]');
ylabel('Orientation[rad]');

subplot(2,2,3)
plot(time, earth_position_x)
hold on
plot(time, earth_position_y)
grid on;
title('Position of the endeffector - EARTH')
legend('x_e_e', 'y_e_e')
xlabel('Time [s]');
ylabel('Position [m]');
hold off

subplot(2,2,4)
plot(time, earth_orientation)
grid on;
title('Orientation of the endeffector - EARTH')
legend('Ψ_e_e')
xlabel('Time [s]');
ylabel('Orientation[rad]');

figure('Name', 'Manipulator on Satellite/Earth Comparision')
subplot(1,2,1)
plot(time,x4, time, x5)
title("Angular position of the joints - Space - for m_0 = " + m0(satelliteMassScalerIndex)  + "kg and I_0 = " + I0(satelliteMassScalerIndex) + "kgm^2.")
legend('Angular position of the first joint \theta_1', 'Angular position of the second joint \theta_2')
xlabel('Time [s]')
ylabel('Angular position [rad]')
grid on;
subplot(1,2,2)
plot(time,earth_q_1, time, earth_q_2)
title("Angular position of the joints - Earth")
legend('Angular position of the first joint \theta_1', 'Angular position of the second joint \theta_2')
xlabel('Time [s]')
ylabel('Angular position [rad]')
grid on;

end

toc
%% Function Definitions:

function [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10] = satelite_inverse_dynamics...
(p1, p2, m0, m1, m2, I0, L1, L2, U_1, U_2, q1, q2, dq1, dq2, tk)

syms x1_s  % x-positon of the satellite's center of mass [m]
syms x2_s  % y-positon of the satellite's center of mass [m]
syms x3_s  % satellite orientation [rad] 
syms x4_s  % angular position of the first joint [rad]
syms x5_s  % angular position of the second joint [rad]
syms x6_s  % the first derivative of x1 [m/sec]
syms x7_s  % the first derivative of x2 [m/sec]
syms x8_s  % the first derivative of x3 [rad/sec]
syms x9_s % the first derivative of x4 [rad/sec]
syms x10_s % the first derivative of x5 [rad/sec]

% Initial Conditions:
a1 = L1/2;
a2 = L2/2;
I1 = 1/12 * m1 * L1^2;
I2 = 1/12 * m2 * L2^2;
h = tk/length(U_1);

% Base (Satellite)
J0v = [1 0 0 0 0 ;
       0 1 0 0 0 ];
J0vt = J0v';
J0w = [0 0 1 0 0 ];
J0wt = J0w';

% First link
var1 = - p1*sin(x3_s) -p2*cos(x3_s) - a1*sin(x3_s+x4_s);
var2 = -a1*sin(x3_s+x4_s);
var3 = -p2*sin(x3_s) + p1*cos(x3_s)+a1*cos(x4_s+x3_s);
var4 = a1*cos(x4_s+x3_s);
J1v = [1 0  var1  var2 0 ;
       0 1  var3  var4 0 ];
J1vt = J1v';  
J1w = [0 0 1 1 0];
J1wt = J1w';

% Second link
var5 = -p2*cos(x3_s) - L1*sin(x4_s+x3_s) - p1*sin(x3_s) - a2*sin(x4_s+x5_s+x3_s);
var6 = -L1*sin(x4_s+x3_s) - a2*sin(x4_s+x5_s+x3_s);
var7 = -a2*sin(x4_s+x5_s+x3_s);
var8 = -p2*sin(x3_s) + a2*cos(x4_s+x5_s+x3_s) + L1*cos(x4_s+x3_s)+p1*cos(x3_s);
var9 =  a2*cos(x4_s+x5_s+x3_s) + L1*cos(x4_s+x3_s);
var10 = a2*cos(x4_s+x5_s+x3_s);
J2v = [1 0   var5 var6  var7 ;
       0 1   var8 var9 var10 ];   
J2vt = J2v';     
J2w = [0 0 1 1 1 ];
J2wt = J2w';

% M matrix

Tv = m0*J0vt*J0v +  m1*J1vt*J1v +  m2*J2vt*J2v;
Tw = J0wt*I0*J0w +  J1wt*I1*J1w +  J2wt*I2*J2w;
M  = Tv+Tw;

% C matrix

Qdot = [x6_s; x7_s; x8_s; x9_s; x10_s];
Qdott =[x6_s x7_s x8_s x9_s x10_s ];
Q = [x1_s; x2_s; x3_s; x4_s; x5_s];

for i=1:1:size(Qdot,1)
    for j=1:1:size(Qdot,1)
        sum = 0;
        for k=1:1:size(Qdot,1)
           sum = sum + diff(M(i,j),Q(k))*Qdot(k);
        end
        N1(i,j)=sum;
    end
end
for i=1:1:size(Qdot,1)
    for j=1:1:size(Qdot,1)
        for k=1:1:size(Qdot,1)     
          P(j,k) = diff(M(j,k),Q(i));
        end
    end
    N2(i,1) = Qdott*P*Qdot;
    clear P;
end

C = N1*Qdot - 0.5*N2;

gFH_M = matlabFunction(M);
gFH_C = matlabFunction(C);

x1=zeros(1,length(U_1)-1);
x2=zeros(1,length(U_1)-1);
x3=zeros(1,length(U_1)-1);
x4=zeros(1,length(U_1)-1);
x5=zeros(1,length(U_1)-1);
x6=zeros(1,length(U_1)-1);
x7=zeros(1,length(U_1)-1);
x8=zeros(1,length(U_1)-1);
x9=zeros(1,length(U_1)-1);
x10=zeros(1,length(U_1)-1);

x4(1)=q1;
x5(1)=q2;
x9(1)=dq1;
x10(1)=dq2;

%% Ode 4

for iteration = 1:1:(length(U_1)-1)

U = [0; 0; 0; U_1(iteration); U_2(iteration)];

% n1

M_v = gFH_M(x3(iteration),x4(iteration),x5(iteration));

C_v1 = gFH_C(x10(iteration),x3(iteration),x4(iteration),x5(iteration),x6(iteration),x7(iteration),x8(iteration),x9(iteration));

Qdot_v = [x6(iteration),x7(iteration),x8(iteration),x9(iteration),x10(iteration)];

k11 = h*Qdot_v;
k12 = h*inv(M_v)*(U-C_v1);
% end n1

% n2
k21 = h*(Qdot_v + 1/2*k11);

x6_2 = x6(iteration) + 1/2 * k12(1);
x7_2 = x7(iteration) + 1/2 * k12(2);
x8_2 = x8(iteration) + 1/2 * k12(3);
x9_2 = x9(iteration) + 1/2 * k12(4);
x10_2 = x10(iteration) + 1/2 * k12(5);

C_v2 = gFH_C(x10_2,x3(iteration),x4(iteration),x5(iteration),x6_2,x7_2,x8_2,x9_2);

k22 = h*inv(M_v)*(U-C_v2);
% end n2

% n3
k31 = h*(Qdot_v + 1/2*k21);

x6_2 = x6(iteration) + 1/2 * k22(1);
x7_2 = x7(iteration) + 1/2 * k22(2);
x8_2 = x8(iteration) + 1/2 * k22(3);
x9_2 = x9(iteration) + 1/2 * k22(4);
x10_2 = x10(iteration) + 1/2 * k22(5);

C_v3 = gFH_C(x10_2,x3(iteration),x4(iteration),x5(iteration),x6_2,x7_2,x8_2,x9_2);
k32 = h*inv(M_v)*(U-C_v3);
% end n3

% n4
k41 = h*(Qdot_v + k31);

x6_2 = x6(iteration) + k32(1);
x7_2 = x7(iteration) + k32(2);
x8_2 = x8(iteration) + k32(3);
x9_2 = x9(iteration) + k32(4);
x10_2 = x10(iteration) + k32(5);

C_v4 = gFH_C(x10_2,x3(iteration),x4(iteration),x5(iteration),x6_2,x7_2,x8_2,x9_2);

k42 = h*inv(M_v)*(U-C_v4);

% end n4

k_1 = (1/6)*(k11 + 2*k21 + 2*k31 + k41);
k_2 = (1/6)*(k12 + 2*k22 + 2*k32 + k42);

x1(iteration+1) = x1(iteration) + k_1(1);
x2(iteration+1) = x2(iteration) + k_1(2);
x3(iteration+1) = x3(iteration) + k_1(3);
x4(iteration+1) = x4(iteration) + k_1(4);
x5(iteration+1) = x5(iteration) + k_1(5);

x6(iteration+1) = x6(iteration) + k_2(1);
x7(iteration+1) = x7(iteration) + k_2(2);
x8(iteration+1) = x8(iteration) + k_2(3);
x9(iteration+1) = x9(iteration) + k_2(4);
x10(iteration+1) = x10(iteration) + k_2(5);

end
end

function[Pee_x,Pee_y,Psiee,Psiee_q_x,Psiee_q_y, Vee_x, Vee_y,Omee] = direct_2DoF(q, dq, L1, L2)

Pee_x = L1 * cos(q(1)) + L2 * cos(q(1) + q(2));
Pee_y = L1 * sin(q(1)) + L2 * sin(q(1) + q(2));

Psiee = q(1) + q(2);

Psiee_q_x = q(1);
Psiee_q_y = q(2);

Vee_x = -L1 * dq(1) * sin(q(1)) - L2 * (dq(1) + dq(2)) * sin(q(1) + q(2));
Vee_y = L1 * dq(1) * cos(q(1)) + L2 * (dq(1) + dq(2)) * cos(q(1) + q(2));

Omee = dq(1) + dq(2);

end
 
function [r1] = Trajectory_Generation(q_p,dq_p,q_k,dq_k,t_k,t_a,t_b,dt,V)
syms a_0 a_1 a_2 a_3 b_0 b_1 c_0 c_1 c_2 c_3 t_c

[r1] = solve(...
    a_0 - q_p == 0, ...
        ...
    a_1 - dq_p == 0, ...
        ...
    a_0 + a_1 * t_a + a_2 * t_a^2 + a_3 * t_a^3 - b_0 - b_1 * t_a == 0, ...
        ...
    a_1 + 2 * a_2 * t_a + 3 * a_3 * t_a^2 - b_1 == 0,...
        ...
    a_1 + 2 * a_2 * t_a + 3 * a_3 * t_a^2 - c_1 - 2 * c_2 * t_c + 3 * c_3 * t_c^2 == 0, ...
        ...
    V - b_1 == 0, ...
        ...
    b_0 + b_1 * t_c - c_0 - c_1 * t_c - c_2 * t_c^2 - c_3 * t_c^3 == 0, ...
        ...
    b_1 - c_1 - 2 * c_2 * t_c - 3 * c_3 * t_c^2 == 0, ...
        ...
    c_0 + c_1 * t_k + c_2 * t_k^2 + c_3 * t_k^3 - q_k == 0, ...
        ...
    c_1 + 2 * c_2 * t_k + 3 * c_3 * t_k^2 - dq_k == 0, ...
        ...
    t_c - t_k + t_b == 0, ...
        ...
    a_0, a_1, a_2, a_3, b_0, b_1, c_0, c_1, c_2, c_3, t_c);

index = 1;

q = zeros(1,length(0:dt:t_k));
dq = zeros(1,length(0:dt:t_k));
ddq = zeros(1,length(0:dt:t_k));

for t = 0:dt:t_k
    if(t <= t_a)
        q(index) = r1.a_0 + r1.a_1 * t + r1.a_2 * t^2 + r1.a_3 * t^3;
        dq(index) = r1.a_1 + 2 * r1.a_2 * t + 3 * r1.a_3* t^2;
        ddq(index) = 2 * r1.a_2 + 6 * r1.a_3 * t;
        index = index + 1;
    else if (t < t_k - t_b && t > t_a)
        q(index) = r1.b_0 + r1.b_1 * t;
        dq(index) = r1.b_1;
        ddq(index) = 0;
        index = index + 1;
    else
        q(index) = r1.c_0 + r1.c_1 * t + r1.c_2 * t^2 + r1.c_3 * t^3;
        dq(index) = r1.c_1 + 2 * r1.c_2 * t + 3 * r1.c_3* t^2;
        ddq(index)= 2 * r1.c_2 + 6 * r1.c_3 * t;
        index = index + 1;
    end
    end

end
r1.q = q;
r1.dq = dq;
r1.ddq = ddq;
end