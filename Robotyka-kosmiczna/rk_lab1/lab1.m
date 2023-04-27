clear all
syms a_0 a_1 a_2 a_3 b_0 b_1 c_0 c_1 c_2 c_3
syms q_p q_k v
syms t_a t_c t_k t_b
syms d_qk d_qp

t_c = t_k-t_b;

A = a_0 - q_p

B = a_1 - d_qp;

C = a_0 + a_1*t_a + a_2*(t_a)^2 + a_3*(t_a)^3 - b_0 - b_1*t_a;

D = a_1 + 2*a_2*t_a + 3*a_3*(t_a)^2 - b_1;

E = a_1 + 2*a_2*t_a + 3*a_3*(t_a)^2 - c_1 - 2*c_2*t_c - 3*c_3*(t_c)^2;

F = b_1 - v;

G = b_0 + b_1*t_c - c_0 - c_1*t_c - c_2*(t_c)^2 - c_3*(t_c)^3;

H = b_1 - c_1 - 2*c_2*t_c - 3*c_3*(t_c)^2;

I = c_0 + c_1*t_k + c_2*(t_k)^2 + c_3*(t_k)^3 -q_k;

J = c_1 + 2*c_2*t_k + 3*c_3*(t_k)^2 - d_qk;

[Z] = solve(A==0, B==0, C==0, D==0, E==0, F==0, G==0, H==0, I==0, J==0 ,a_0, a_1, a_2, a_3, b_0, b_1, c_0, c_1, c_2, c_3);



solved = Z;
q_p_in = 0;
dq_p_in = 20;
q_k_in = 100;
dq_k_in = 2;
Tk_in = 12;
Ta_in = 4;
Tb_in = 6;
dt_in = 1;
V_in = 10;
out = Trajectory_Generation(solved, q_p_in, dq_p_in, q_k_in, dq_k_in, Tk_in, Ta_in, Tb_in, dt_in, V_in);


function [q,dq,ddq] = Trajectory_Generation(solved, q_p_in,dq_p_in,q_k_in,dq_k_in,Tk_in,Ta_in,Tb_in,dt_in,V_in)
    syms q_p q_k v
    syms t_a t_c t_k t_b
    syms d_qk d_qp
    val = subs(solved,[q_p, d_qp, q_k, d_qk, t_k, t_a, t_b, v],[q_p_in dq_p_in q_k_in dq_k_in Tk_in Ta_in Tb_in V_in]);

    t1 = 0:dt_in:Ta_in;
    q_1 = val.a_0 + val.a_1 * t1 + val.a_2 * (t1).^2 + val.a_3 * (t1).^3;
    d_q1 = val.a_1 + 2*val.a_2*t1 + 3*val.a_3 * (t1).^2;
    dd_q1 = 2*val.a_2 + 6*val.a_3*t1;
    
    t2 = Ta_in:dt_in:Tb_in;
    q_2 = val.b_0 + val.b_1 * t2;
    d_q2 = val.b_1 + zeros(1, length(t2));
    dd_q2 = zeros(1, length(t2));
    
    t3 = Tb_in:dt_in:Tk_in;
    q_3 = val.c_0 + val.c_1 * t3 + val.c_2 * (t3).^2 + val.c_3 * (t3).^3;
    d_q3 = val.c_1 + 2*val.c_2*t3 + 3*val.c_3 * (t3).^2;
    dd_q3 = 2*val.c_2 + 6*val.c_3*t3;

    t = [0:dt_in:Tk_in];
    q = [q_1, q_2, q_3];
    dq = [d_q1, d_q2, d_q3];
    ddq = [dd_q1, dd_q2, dd_q3];

    figure;
    subplot(3, 1, 1);
    plot(t1, q_1, t2, q_2, t3, q_3, t1(1), q_p_in, 'o', Tk_in, q_k_in, 'o')
    legend('Faza I', 'Faza II', 'Faza III', 'Pozycja początkowa', 'Pozycja końcowa')
    title('Pozycja')
    grid on
    axis auto

    subplot(3, 1, 2);
    plot(t1, d_q1, t2, d_q2, t3, d_q3, t1(1), dq_p_in, 'o', Tk_in, dq_k_in, 'o')
    legend('Faza I', 'Faza II', 'Faza III', 'Prędkość początkowa', 'Prędkość końcowa')
    title('Prędkość')
    grid on
    axis auto

    subplot(3, 1, 3);
    plot(t1, dd_q1, t2, dd_q2, t3, dd_q3)
    legend('Faza I', 'Faza II', 'Faza III')
    title('Przyspieszenie')
    grid on
    axis auto
end