% Define the parameters for the first and second equations
m1 = 85; % Mass
c1 = 250; % Damping coefficient
k1 = 25000; % Spring constant
k2 = 13600; % Another spring constant
omega1 = 3.14; % Frequency
Y1 = 1; % Amplitude
l1 = 0.23; % Constant
l2 = 0.17; % Constant
l3 = 0.14; % Constant

m2 = 85; % Mass
c2 = 250; % Damping coefficient
k3 = 25000; % Spring constant
omega2 = 3.14; % Frequency
Y2 = 1; % Amplitude

% Time span and initial conditions
tspan = [0, 10];
initial_conditions = [0, 0]; 
h = 0.01; 

% Solve ODEs using RK4
[t1, z1] = RK4(@myODE3, tspan, initial_conditions, h);
[t2, z2] = RK4(@myODE4, tspan, initial_conditions, h);

% Plot the results
figure;
plot(t1, z1(:, 1), 'b-', 'LineWidth', 2, 'DisplayName', 'ODE 1');
hold on;
plot(t2, z2(:, 1), 'r-', 'LineWidth', 2, 'DisplayName', 'ODE 2');
hold off;

xlabel('Time');
ylabel('Displacement');
title('Comparison of Two ODEs using Runge-Kutta 4th Order Method');
legend;
grid on;

% Calculate root mean square displacement (RMSD) for the first ODE
T1 = t1(end);  
rmsd1 = sqrt(trapz(t1, z1(:, 1).^2) / T1);

% Calculate root mean square displacement (RMSD) for the second ODE
T2 = t2(end);  
rmsd2 = sqrt(trapz(t2, z2(:, 1).^2) / T2);

% Display the calculated RMSDs
disp(['Root Mean Square Displacement (RMSD) for ODE 1: ', num2str(rmsd1)]);
disp(['Root Mean Square Displacement (RMSD) for ODE 2: ', num2str(rmsd2)]);

% Define functions at the end of the script
% Solve the first differential equation: mz'' + cz' + k1z + 2k2pz = m(ω^2)Ysin(ωt)
function dzdt1 = myODE3(t, z)
    % Parameters for the first equation
    m1 = 85; % Mass
    c1 = 250; % Damping coefficient
    k1 = 25000; % Spring constant
    k2 = 13600; % Another spring constant
    omega1 = 3.14; % Frequency
    Y1 = 1; % Amplitude
    l1 = 0.23; % Constant
    l2 = 0.17; % Constant
    l3 = 0.14; % Constant

    % Calculate p
    p = (l3 - l1) / sqrt(l2^2 - z(1)^2) + 1;

    % First differential equation
    dzdt1 = [z(2); (m1 * (omega1^2) * Y1 * sin(omega1 * t) - c1 * z(2) - k1 * z(1) - 2 * k2 * p * z(1)) / m1];
end

% Solve the second differential equation: mz'' + cz' + k1z = m(ω^2)Ysin(ωt)
function dzdt2 = myODE4(t, z)
    % Parameters for the second equation
    m2 = 85; % Mass
    c2 = 250; % Damping coefficient
    k3 = 25000; % Spring constant
    omega2 = 3.14; % Frequency
    Y2 = 1; % Amplitude

    % Second differential equation
    dzdt2 = [z(2); (m2 * (omega2^2) * Y2 * sin(omega2 * t) - c2 * z(2) - k3 * z(1)) / m2];
end

% Runge-Kutta 4th Order Method
function [t, z] = RK4(myODE, tspan, initial_conditions, h)
    t_initial = tspan(1);
    t_final = tspan(2);
    
    num_steps = round((t_final - t_initial) / h);
    
    t = zeros(num_steps + 1, 1);
    z = zeros(num_steps + 1, length(initial_conditions));
    
    t(1) = t_initial;
    z(1, :) = initial_conditions;
    
    for i = 1:num_steps
        ti = t(i);
        zi = z(i, :);
        
        k1 = h * myODE(ti, zi)';
        k2 = h * myODE(ti + h/2, zi + k1/2)';
        k3 = h * myODE(ti + h/2, zi + k2/2)';
        k4 = h * myODE(ti + h, zi + k3)';
        
        t(i + 1) = ti + h;
        z(i + 1, :) = zi + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end
