% Multibody Joolian System Simulation
% Clearing
clear all
close all

% Constants
d2r = pi/180; % Degrees to Radians
r2d = 180/pi; % Radians to Degrees
G = 6.674*10^-11; % Newtonian Constant of Gravitation

% Joolian System Characteristics Array
% Radius (m), Mass (kg), Gravitational Parameter (m^3/s^2)
planet_characteristics = cell(1,6); % Preallocation
planet_characteristics{1} = [6000000, 4.2332127e24, 2.8252800e14]; % Jool
planet_characteristics{2} = [500000, 2.9397311e22, 1.9620000e12]; % Laythe
planet_characteristics{3} = [300000, 3.1087655e21, 2.0748150e11]; % Vall
planet_characteristics{4} = [6000000, 4.2332127e22, 2.8252800e12]; % Tylo
planet_characteristics{5} = [65000, 3.7261090e19, 2.4868349e9]; % Bop
planet_characteristics{6} = [44000, 1.0813507e19, 7.2170208e8]; % Pol

% Orbital Characteristics - Initial Conditions
% SMA, Ecc, Inc, AoP, RAAN, MNA, TA (Calculated Later)
orbits = cell(1, 5); % Preallocation
orbits{1} = [27184000, 0, 0, 0, 0, 3.140,0]; % Laythe
orbits{2} = [43152000, 0, 0, 0, 0, 0.900,0]; % Vall
orbits{3} = [68500000, 0, 0.025*d2r, 0, 0, 3.140,0]; % Tylo
orbits{4} = [128500000, 0.235, 15*d2r, 25*d2r, 10*d2r, 0.900,0]; % Bop
orbits{5} = [179890000, 0.171, 4.25*d2r, 15*d2r, 2*d2r, 0.900,0]; % Pol
for i = 1:5 % True Anomaly Calculation
    orbits{i}(1,7) = 2*atan(sqrt((1+orbits{i}(2))/(1-orbits{i}(2)))*tan(KTE(orbits{i}(2), orbits{i}(6))/2));
end

% Initial Position in Perifocal
r_P = cell(1, 6); % Preallocation
r_P{1} = [0; 0; 0]; % Jool's Initial Position in Perifocal
for i = 2:6
    r_P{i} = ((orbits{i-1}(1,1)*(1-orbits{i-1}(2)^2))/(1+orbits{i-1}(1,2)*cos(orbits{i-1}(1,7))))*[cos(orbits{i-1}(1,7)); sin(orbits{i-1}(1,7)); 0];
end

% Initial Velocity in Perifocal
v_P = cell(1, 6); % Preallocation
v_P{1} = [0; 0; 0]; % Jool's Initial Velocity in Perifocal
for i = 2:6
    v_P{i} = sqrt(planet_characteristics{1}(3)/(orbits{i-1}(1)*(1-orbits{i-1}(2)^2)))*[-sin(orbits{i-1}(1,7)); (orbits{i-1}(2) + cos(orbits{i-1}(1,7))); 0];
end

simulation_timestep = linspace(0, 1000000, 1000001);

% Position Matrix Preallocation
planet_positions = cell(1, 6); % Preallocation
for i = 1:6
    planet_positions{i} = zeros(length(simulation_timestep), 3);
end

% Velocity Matrix Preallocation
planet_velocities = cell(1, 6); % Preallocation
for i = 1:6
    planet_velocities{i} = zeros(length(simulation_timestep), 3);
end

% Initial Position and Velocity in JCI
planet_positions{1}(1, :) = [0,0,0]; % Jool's Initial Position in JCI
planet_velocities{1}(1, :) = [0,0,0]; % Jool's Initial Velocity in JCI
for i = 2:6
    C_PJ = create_rot_mtx(3,orbits{i-1}(1,4))*create_rot_mtx(1,orbits{i-1}(1,3))*create_rot_mtx(3,orbits{i-1}(1,5));
    C_JP = transpose(C_PJ);
    planet_positions{i}(1, :) = transpose(C_JP*r_P{i}); % Initial Position in JCI
    planet_velocities{i}(1, :) = transpose(C_JP*v_P{i}); % Initial Velocity in JCI
end

% Figure Setup
planet_names = {'Jool', 'Laythe', 'Vall', 'Tylo', 'Bop', 'Pol'};
planet_colors = [0 1 0; 0 0 1; 0.5 1 0.5; 0.7 0.7 0.7; 0.9 0.8 0.6; 1 1 0]; % Green, Dark Blue, Light Blue, Grey, Beige, Yellow
figure;

% Setup grid
ax = gca;
ax.Color = 'k'; % Set background color to black
ax.GridColor = 'w'; % Set grid color to white
hold on;

% Not really sure
planetMarkers = gobjects(1, 6);  % Store handles for planet markers
planetLabels = gobjects(1, 6);   % Store handles for planet labels
planetTrajectories = gobjects(1, 6);   % Store handles for planet Trajectories

% Customize the plot
axis equal;
axis([-2.15e8 2.15e8 -2.15e8 2.15e8 -2.15e8 2.15e8]);
title('Bodies Positions Over Time');
xlabel('X Position');
ylabel('Y Position');
zlabel('Z Position');
grid on;

% Simulation Loop
for t = 1:length(simulation_timestep)
    if mod(t, 1000) == 0
        disp(t)
        delete(planetMarkers); % Delete Previous Marker
        delete(planetLabels); % Delete Previous Label
        delete(planetTrajectories); % Delete Previous Planet Trajectory Points
        traj_points = 10000; % # of points plotted
        for i = 1:6
            if t > traj_points
                planetTrajectories(i) = plot3(planet_positions{i}(t-traj_points:1000:t, 1), planet_positions{i}(t-traj_points:1000:t, 2), planet_positions{i}(t-traj_points:1000:t, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
            else
                planetTrajectories(i) = plot3(planet_positions{i}(1:1000:t, 1), planet_positions{i}(1:1000:t, 2), planet_positions{i}(1:1000:t, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
            end
            planetMarkers(i) = plot3(planet_positions{i}(t, 1), planet_positions{i}(t, 2), planet_positions{i}(t, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', planet_colors(i, :));
            planetLabels(i) = text(planet_positions{i}(t, 1), planet_positions{i}(t, 2), planet_positions{i}(t, 3), planet_names{i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        end
        drawnow
    end
    for i = 1:6 % Current body
        tot_force_vect = [0, 0, 0]; % Initialize total force vector for each body

        for j = 1:6 % Interactions between current body and other bodies
            if i ~= j % Ignores interactions between current body and itself
                % Force vector calculation
                force_ij_mag = (G*planet_characteristics{i}(1, 2)*planet_characteristics{j}(1, 2)) / (norm(planet_positions{j}(t, :) - planet_positions{i}(t, :))^2);
                force_ij_vect = force_ij_mag * ((planet_positions{j}(t, :) - planet_positions{i}(t, :)) / norm(planet_positions{j}(t, :) - planet_positions{i}(t, :)));
                tot_force_vect = tot_force_vect + force_ij_vect; % Sum force vectors every loop
            end
        end

        % Calculate acceleration, update position and velocity
        a_i = tot_force_vect / planet_characteristics{i}(1, 2);
        planet_velocities{i}(t+1, :) = planet_velocities{i}(t, :) + a_i;
        planet_positions{i}(t+1, :) = planet_positions{i}(t, :) + planet_velocities{i}(t, :);
    end
end
hold off;

% Functions
function eccentric_anomaly = KTE(ecc, mna) % Iterative Method for Solving KTE
    guess_mtx = zeros(1:5); % Preallocation
    guess_mtx(1,1) = mna; % Mean Anomaly
    for i = 1:4
        guess_mtx(1, i+1) = guess_mtx(1, i) - (guess_mtx(1, i) - ecc * sin(guess_mtx(1, i)) - guess_mtx(1, 1)) / (1 - ecc * guess_mtx(1, i));
    end
    eccentric_anomaly = guess_mtx(1,5); % Output Eccentric Anomaly (radians)
end

function rot_mtx = create_rot_mtx(axis, angle) % Rotation Matrix (x-rot input is 1, y-rot input is 2, z-rot input is 3)
if axis==1
    rot_mtx=[1,0,0;0,cos(angle),sin(angle);0,-sin(angle),cos(angle)];  
elseif axis==2
    rot_mtx=[cos(angle),0,-sin(angle);0,1,0;sin(angle),0,cos(angle)];  
elseif axis==3
    rot_mtx=[cos(angle),sin(angle),0;-sin(angle),cos(angle),0;0,0,1];  
else
    disp('Unacceptable Function Argument')
end
end
