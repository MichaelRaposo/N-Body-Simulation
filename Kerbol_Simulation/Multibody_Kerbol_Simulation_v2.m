%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multibody Kerbol System Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clearing
clear all
close all

% Constants and Preallocation
d2r = pi/180; % Degrees to Radians
r2d = 180/pi; % Radians to Degrees
G = 6.674*10^-11; % Newtonian Constant of Gravitation

% Preallocation and Initialization
num_bodies = 17;
body_prop = cell(1, num_bodies);                                 % Preallocation
orbits = cell(1, num_bodies);                                   % Preallocation
r_P = cell(1, num_bodies);                                      % Preallocation
v_P = cell(1, num_bodies);                                      % Preallocation
planet_position = cell(1, num_bodies);                          % Preallocation
planet_velocity = cell(1, num_bodies);                          % Preallocation
planetMarkers = gobjects(1, num_bodies);         % Store handles for planet markers
planetLabels = gobjects(1, num_bodies);          % Store handles for planet labels
planetTrajectories = gobjects(1, num_bodies);    % Store handles for planet trajectories
planetMarkers2 = gobjects(1, num_bodies);        % Store handles for planet markers
planetLabels2 = gobjects(1, num_bodies);         % Store handles for planet labels
planetTrajectories2 = gobjects(1, num_bodies);   % Store handles for planet trajectories
timeCounter = gobjects(1); % Store time for display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Characteristics and Orbital Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Body names and colours - closest planet to star, followed by closest moons
body_names = {'Kerbol','Moho','Eve','Kerbin','Duna','Dres','Jool','Eeloo','Gilly','Mun','Minmus','Ike','Laythe', 'Vall', 'Tylo', 'Bop', 'Pol'};
planet_colors = [1 0.6471 0; 0.8235 0.7059 0.5490; 0.8627 0.6275 0.8627; 0 0 1; 1 0 0; 0.5882 0.5294 0.4706; 0 1 0; 1 1 1; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.7294 0.9255 0.6941; 0.5 0.5 0.5; 0 0 1; 0.5 1 0.5; 0.7 0.7 0.7; 0.9 0.8 0.6; 1 1 0];

% Body Properties: Radius (m), Mass (kg), Gravitational Parameter (m^3/s^2), SOI (m)
body_prop{1} = [261600000, 1.7565459e28, 1.1723328e18, 1e12];      % Kerbol (REF=1)
body_prop{2} = [250000, 2.5263314e21, 1.6860938e11, 9646663];      % Moho (REF=2)
body_prop{3} = [700000, 1.2243980e23, 8.1717302e12, 85109365];     % Eve (REF=3)
body_prop{4} = [600000, 5.2897088e22, 3.5303940e12, 84147790];     % Kerbin (REF=4)
body_prop{5} = [320000, 4.5154270e21, 3.0136321e11, 47921949];     % Duna (REF=5)
body_prop{6} = [138000, 3.2190937e20, 2.1484489e10, 32832840];     % Dres (REF=6)
body_prop{7} = [6000000, 4.2332127e24, 2.8252800e14, 2.4559852e9]; % Jool (REF=7)
body_prop{8} = [210000, 1.1149224e21, 7.4410815e10, 1.1908294e8];  % Eeloo (REF=8)
body_prop{9} = [13000, 1.2420363e17, 8289449.8, 126123.27];        % Gilly (REF=9)
body_prop{10} = [200000, 9.7599066e20, 6.5138398e10, 2429891.1];   % Mun (REF=10)
body_prop{11} = [60000, 2.6457580e19, 1.7658000e9, 2247735.4];     % Minmus (REF=11)
body_prop{12} = [130000, 2.7821615e20, 1.8568369e10, 1049598.9];   % Ike (REF=12)
body_prop{13} = [500000, 2.9397311e22, 1.9620000e12, 3723645.8];   % Laythe (REF=13)
body_prop{14} = [300000, 3.1087655e21, 2.0748150e11, 2406401.4];   % Vall (REF=14)
body_prop{15} = [600000, 4.2332127e22, 2.8252800e12, 10856518];    % Tylo (REF=15)
body_prop{16} = [65000, 3.7261090e19, 2.4868349e9, 1221060.9];     % Bop (REF=16)
body_prop{17} = [44000, 1.0813507e19, 7.2170208e8, 1042138.9];     % Pol (REF=17)

% Orbital Characteristics - Initial Conditions
% Reference Body (Parent Body), SMA, Ecc, Inc, AoP, RAAN, MNA, TA (Calculated Later)
orbits{1} = [1, 0, 0, 0, 0, 0, 0, 0];                                   % Kerbol (REF=1)
orbits{2} = [1, 5263138304, 0.2, d2r*7, d2r*15, d2r*70, 3.14, 0];       % Moho (REF=2)
orbits{3} = [1, 9832684544, 0.01, d2r*2.1, 0, d2r*15, 3.14, 0];         % Eve (REF=3)
orbits{4} = [1, 13599840256, 0, 0, 0, 0, 3.14, 0];                      % Kerbin (REF=4)
orbits{5} = [1, 20726155264, 0.051, d2r*0.06, 0, d2r*135.5, 3.14, 0];   % Duna (REF=5)
orbits{6} = [1, 40839348203, 0.145, d2r*5, d2r*90, d2r*280, 3.14, 0];   % Dres (REF=6)
orbits{7} = [1, 68773560320, 0.05, d2r*1.304, 0, d2r*52, 0.1, 0];       % Jool (REF=7)
orbits{8} = [1, 90118820000, 0.26, d2r*6.15, d2r*260, d2r*50, 3.14, 0]; % Eeloo (REF=8)
orbits{9} = [3, 31500000, 0.55, d2r*12, d2r*10, d2r*80, 0.9, 0];        % Gilly (REF=9)
orbits{10} = [4, 12000000, 0, 0, 0, 0, 1.7, 0];                         % Mun (REF=10)
orbits{11} = [4, 47000000, 0, d2r*6, d2r*38, d2r*78, 0.9, 0];           % Minmus (REF=11)
orbits{12} = [5, 3200000, 0.03, d2r*0.2, 0, 0, 1.7, 0];                 % Ike (REF=12)
orbits{13} = [7, 27184000, 0, 0, 0, 0, 3.14, 0];                        % Laythe (REF=13)
orbits{14} = [7, 43152000, 0, 0, 0, 0, 0.9, 0];                         % Vall (REF=14)
orbits{15} = [7, 68500000, 0, 0.025*d2r, 0, 0, 3.14, 0];                % Tylo (REF=15)
orbits{16} = [7, 128500000, 0.235, 15*d2r, 25*d2r, 10*d2r, 0.9, 0];     % Bop (REF=16)
orbits{17} = [7, 179890000, 0.171, 4.25*d2r, 15*d2r, 2*d2r, 0.9, 0];    % Pol (REF=17)

% True Anomaly Calculation
for i = 1:num_bodies
    orbits{i}(1,8) = 2*atan(sqrt((1+orbits{i}(3))/(1-orbits{i}(3)))*tan(KTE(orbits{i}(3), orbits{i}(7))/2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Positions and Velocities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial Position in Perifocal
r_P{1} = [0;0;0]; % Kerbol Initial Position in Perifocal
for i = 2:num_bodies
    r_P{i} = ((orbits{i}(1,2)*(1-orbits{i}(3)^2))/(1+orbits{i}(1,3)*cos(orbits{i}(1,8))))*[cos(orbits{i}(1,8)); sin(orbits{i}(1,8)); 0];
end

v_P{1} = [0;0;0]; % Kerbol Initial Velocity in Perifocal
for i = 2:num_bodies
    v_P{i} = sqrt(body_prop{orbits{i}(1,1)}(1,3)/(orbits{i}(2)*(1-orbits{i}(3)^2)))*[-sin(orbits{i}(1,8)); (orbits{i}(1,3) + cos(orbits{i}(1,8))); 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_time = linspace(0, 1, 1); % Simulation length in seconds
sim_percent = (length(sim_time)-1)/1000; % Simulation percentage constant
traj_points = 2500000; % Position/velocity stored at once

% Position Matrix Preallocation
for i = 1:num_bodies
    planet_position{i} = zeros(length(sim_time), 3);
end

% Velocity Matrix Preallocation
for i = 1:num_bodies
    planet_velocity{i} = zeros(length(sim_time), 3);
end

% Initial Position and Velocity in KoCI (Kerbol Centered Inertial)
planet_position{1}(1, :) = [0,0,0]; % Kerbol's Initial Position in KoCI
planet_velocity{1}(1, :) = [0,0,0]; % Kerbol's Initial Velocity in KoCI
for i = 2:num_bodies
    C_PKo = create_rot_mtx(3,orbits{i}(1,5))*create_rot_mtx(1,orbits{i}(1,4))*create_rot_mtx(3,orbits{i}(1,6));
    C_KoP = transpose(C_PKo);
    planet_position{i}(1, :) = transpose(C_KoP*r_P{i}) + planet_position{orbits{i}(1,1)}(1,:); % Initial Position in KoCI
    planet_velocity{i}(1, :) = transpose(C_KoP*v_P{i}) + planet_velocity{orbits{i}(1,1)}(1,:); % Initial Velocity in KoCI
end

for t = 1:length(sim_time) % Simulation Loop
    force_array = cell(num_bodies,num_bodies); 
    if mod(t,sim_percent)==0
        disp('Progress: ' + string(round((100*t)/length(sim_time),2))+'%');   
    end
    for i = 1:num_bodies % Current Body
        for j = i:num_bodies % Interacting Body
            if i ~= j % Ignores interaction between current body and itself
                % Force vector calculation
                force_ij_mag = (G*body_prop{i}(1, 2)*body_prop{j}(1, 2)) / (norm(planet_position{j}(t, :) - planet_position{i}(t, :))^2);
                force_ij_vect = force_ij_mag * ((planet_position{j}(t, :) - planet_position{i}(t, :)) / norm(planet_position{j}(t, :) - planet_position{i}(t, :)));
                force_array{i,j} = force_ij_vect;
            end
        end
    end
    %force_array = cellfun(@(x, y) [x, y], force_array, transpose(cellfun(@(x) x*-1,force_array,'UniformOutput',false)), 'UniformOutput', false);
    for i = 1:length(force_array)
        force_array{i,i} = [0,0,0];
    end
    for i = 1:num_bodies % Current Body
        total_force_vect = [0,0,0];
        for j = 1:i
            total_force_vect = total_force_vect - force_array{j,i};
        end
        for j = i:num_bodies
            total_force_vect = total_force_vect + force_array{i,j};
        end
        a_i = total_force_vect / body_prop{i}(1, 2);
        planet_velocity{i}(t+1, :) = planet_velocity{i}(t, :) + a_i;
        planet_position{i}(t+1, :) = planet_position{i}(t, :) + planet_velocity{i}(t, :);
    end
end

%%

% Figure Setup
figureSize = [100, 100, 1200, 800]; % Pixels
tiledlayout(2, 3, 'TileSpacing','compact'); % 2 rows, 3 columns
tl = gcf; % Get the current figure handle
tl.OuterPosition = figureSize;

% Create a VideoWriter object
movieFileName = 'Kerbol_System_Simulation.mp4';
movieWriter = VideoWriter(movieFileName, 'MPEG-4');
movieWriter.FrameRate = 144;  % Set the desired frame rate
spf = 1000; % Seconds per frame
open(movieWriter);

% In-Plane Plot (Tile 1)
nexttile([2, 2])
ax1 = gca;
ax1.Color = 'k'; % Set background color to black
ax1.GridColor = 'w'; % Set grid color to white
hold on;
axis(ax1, 'equal');
axis(ax1, [-1.2e11 1.2e11 -1.2e11 1.2e11]);
title(ax1, 'In-Plane Motion');
xlabel(ax1, 'X Position (meters)');
ylabel(ax1, 'Y Position (meters)');
grid(ax1, 'on');

% Out-of-Plane Plot (Tile 2)
nexttile([2, 1])
ax2 = gca;
ax2.Color = 'k'; % Set background color to black
ax2.GridColor = 'w'; % Set grid color to white
hold on;
axis(ax2, 'equal');
axis(ax2, [-0.42e11 0.42e11 -1.2e11 1.2e11]);
title(ax2, 'Out-of-Plane Motion');
xlabel(ax2, 'Z Position (meters)');
ylabel(ax2, 'Y Position (meters)');
grid(ax2, 'on');

for t = 1:length(sim_time)
    if mod(t, spf) == 0
        % Delete Previous Markers on ax1
        delete(planetMarkers); % Delete Previous Marker
        delete(planetLabels); % Delete Previous Label
        delete(planetTrajectories); % Delete Previous Planet Trajectory Points

        % Delete Previous Markers on ax2
        delete(planetMarkers2); % Delete Previous Marker on ax2getframe
        delete(planetLabels2); % Delete Previous Label on ax2
        delete(planetTrajectories2); % Delete Previous Planet Trajectory Points on ax2
        delete(timeCounter);

        traj_points = 2500000; % # of points plotted
        for i = 1:num_bodies % Planets
            if orbits{i}(1,1) == 1
                xlim(ax1, [planet_position{1}(t, 1) - 1.2e11, planet_position{1}(t, 1) + 1.2e11]);
                ylim(ax1, [planet_position{1}(t, 2) - 1.2e11, planet_position{1}(t, 2) + 1.2e11]);

                xlim(ax2, [planet_position{1}(t, 3) - 0.42e11, planet_position{1}(t, 3) + 0.42e11]);
                ylim(ax2, [planet_position{1}(t, 2) - 1.2e11, planet_position{1}(t, 2) + 1.2e11]);

                if t > traj_points
                    planetTrajectories(i) = plot(ax1, planet_position{i}((t-traj_points):1000:t, 1), planet_position{i}((t-traj_points):1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                    planetTrajectories2(i) = plot(ax2, planet_position{i}((t-traj_points):1000:t, 3), planet_position{i}((t-traj_points):1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                else
                    planetTrajectories(i) = plot(ax1, planet_position{i}(1:1000:t, 1), planet_position{i}(1:1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                    planetTrajectories2(i) = plot(ax2, planet_position{i}(1:1000:t, 3), planet_position{i}(1:1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                end

                planetMarkers(i) = plot(ax1, planet_position{i}(t, 1), planet_position{i}(t, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                planetLabels(i) = text(ax1, planet_position{i}(t, 1), planet_position{i}(t, 2), body_names{i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                planetMarkers2(i) = plot(ax2, planet_position{i}(t, 3), planet_position{i}(t, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                planetLabels2(i) = text(ax2, planet_position{i}(t, 3), planet_position{i}(t, 2), body_names{i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            
            elseif norm(planet_position{i}(t,:) - planet_position{orbits{i}(1,1)}(t,:)) > body_prop{orbits{i}(1,1)}(1,4)
                xlim(ax1, [planet_position{1}(t, 1) - 1.2e11, planet_position{1}(t, 1) + 1.2e11]);
                ylim(ax1, [planet_position{1}(t, 2) - 1.2e11, planet_position{1}(t, 2) + 1.2e11]);

                xlim(ax2, [planet_position{1}(t, 3) - 0.42e11, planet_position{1}(t, 3) + 0.42e11]);
                ylim(ax2, [planet_position{1}(t, 2) - 1.2e11, planet_position{1}(t, 2) + 1.2e11]);

                if t > traj_points
                    planetTrajectories(i) = plot(ax1, planet_position{i}((t-traj_points):1000:t, 1), planet_position{i}((t-traj_points):1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                    planetTrajectories2(i) = plot(ax2, planet_position{i}((t-traj_points):1000:t, 3), planet_position{i}((t-traj_points):1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                else
                    planetTrajectories(i) = plot(ax1, planet_position{i}(1:1000:t, 1), planet_position{i}(1:1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                    planetTrajectories2(i) = plot(ax2, planet_position{i}(1:1000:t, 3), planet_position{i}(1:1000:t, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                end

                planetMarkers(i) = plot(ax1, planet_position{i}(t, 1), planet_position{i}(t, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                planetLabels(i) = text(ax1, planet_position{i}(t, 1), planet_position{i}(t, 2), body_names{i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                planetMarkers2(i) = plot(ax2, planet_position{i}(t, 3), planet_position{i}(t, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                planetLabels2(i) = text(ax2, planet_position{i}(t, 3), planet_position{i}(t, 2), body_names{i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            
            end
        end
        timeCounter = annotation('textbox', [0, 0, 1, 0.05], 'String', ['Time: ' num2str(t) ' s elapsed'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'BackgroundColor', 'w');

        % Capture the current frame and write to the video file
        frame = getframe(gcf);
        writeVideo(movieWriter, frame);

        drawnow
    end
end
close(movieWriter); % Close the VideoWriter object
hold off;

% Functions
function eccentric_anomaly = KTE(ecc, mna) % Iterative Method for Solving KTE
    guess_mtx = zeros(1:10); % Preallocation
    guess_mtx(1,1) = mna; % Mean Anomaly
    for i = 1:9
        guess_mtx(1, i+1) = guess_mtx(1, i) - (guess_mtx(1, i) - ecc * sin(guess_mtx(1, i)) - guess_mtx(1, 1)) / (1 - ecc * guess_mtx(1, i));
    end
    eccentric_anomaly = guess_mtx(1,10); % Output Eccentric Anomaly (radians)
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
