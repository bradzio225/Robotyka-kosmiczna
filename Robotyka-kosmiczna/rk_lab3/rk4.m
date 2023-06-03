clear all

h = 0.1; % h = dt
N = 100;

% time = zeros(N,1);
% y = zeros(N,1);

time(1) = 0;
y(1) = 1;

for i=1:1:N
   k1 = h * funkcja(time(i), y(i));
   k2 = h * funkcja(time(i) + 0.5 * h, y(i) + 0.5 * k1);
   k3 = h * funkcja(time(i) + 0.5 * h, y(i) + 0.5 * k2);
   k4 = h * funkcja(time(i) + h, y(i) + k3);
   k = (1/6) * (k1 + 2 * k2 + 2 * k3 + k4);
   
   y(i+1) = y(i) + k;
   time(i+1) = time(i) + h;
   
end

figure
hold on
plot(time, y);

%sprawdzenie

% syms y(t) t
% rownanie = diff(y,t) == -3*y + t + 6;
% solution = dsolve(rownanie);

time2(1) = 0;
y2(1) = 1;
C3 = -8;

for i=1:1:N
   time2(i+1) = time(i) + h;
   y2(i + 1) = time2(i+1)/3 + (C3*exp(-3*time2(i+1)))/9 + 17/9;
end

plot(time2, y2)
