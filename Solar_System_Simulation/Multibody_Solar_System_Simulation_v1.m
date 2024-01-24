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
num_bodies = 30;
body_prop = cell(1, num_bodies);                    % Preallocation
body_IDs = cell(2, num_bodies);                     % Preallocation
orbits = cell(1, num_bodies);                       % Preallocation
r_P = cell(1, num_bodies);                          % Preallocation
v_P = cell(1, num_bodies);                          % Preallocation
planet_position = cell(1, num_bodies);              % Preallocation
planet_velocity = cell(1, num_bodies);              % Preallocation
planetMarkers_1 = gobjects(1, num_bodies);          % Store handles for planet position
planetMarkers_2 = gobjects(1, num_bodies);          % Store handles for planet position
planetMarkers_3 = gobjects(1, num_bodies);          % Store handles for planet position
planetMarkers_4 = gobjects(1, num_bodies);          % Store handles for planet position
planetTrajectories_1 = gobjects(1, num_bodies);     % Store handles for planet trajectories
planetTrajectories_2 = gobjects(1, num_bodies);     % Store handles for planet trajectories
planetTrajectories_3 = gobjects(1, num_bodies);     % Store handles for planet trajectories
planetTrajectories_4 = gobjects(1, num_bodies);     % Store handles for planet trajectories
planetTrajectories_1_2 = gobjects(1, num_bodies);   % Store handles for planet trajectories for t > traj_points
planetTrajectories_2_2 = gobjects(1, num_bodies);   % Store handles for planet trajectories for t > traj_points
planetTrajectories_3_2 = gobjects(1, num_bodies);   % Store handles for planet trajectories for t > traj_points
planetTrajectories_4_2 = gobjects(1, num_bodies);   % Store handles for planet trajectories for t > traj_points
planetLabels_1 = gobjects(1, num_bodies);           % Store handles for planet labels
planetLabels_2 = gobjects(1, num_bodies);           % Store handles for planet labels
planetLabels_3 = gobjects(1, num_bodies);           % Store handles for planet labels
planetLabels_4 = gobjects(1, num_bodies);           % Store handles for planet labels
timeCounter = gobjects(1);                          % Store time for display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Characteristics and Orbital Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Body names and colours - closest planet to star, followed by closest moons
body_names = cell(2,num_bodies);
body_names(1,:) = {'Sun','Mercury','Venus','Earth','Mars','Ceres','Vesta','Jupiter','Saturn','Uranus','Neptune','Pluto','Moon','Phobos','Deimos', 'Io', 'Europa', 'Ganymede', 'Callisto','Enceladus','Tethys','Dione','Rhea','Titan','Miranda','Ariel','Umbriel','Titania','Oberon','Triton'};
url = 'https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND=%27MB%27';
API_data = webread(url);
API_lines = strsplit(API_data,'\n');
for i = 17:319 % Finds Planet's ID's, where not found copies name
    API_cells = strsplit(API_lines{i},' ');
    for j = 1:num_bodies
        if contains(lower(body_names{1,j}),lower(API_cells{3}))
            body_names{2,j} = API_cells{2};
        elseif isempty(body_names{2,j})
            body_names{2,j} = body_names{1,j};
        end
    end
end

planet_colors = zeros(length(body_names), 3);  % Initialize matrix with zeros
planet_colors(1, :) = [1 0.6471 0];  % Sun - Yellow
planet_colors(2:12, :) = [0.8235 0.7059 0.5490; 0.8627 0.6275 0.8627; 0 0.5 1; 1 0.2 0; 0.5882 0.5294 0.4706; 0.5882 0.5294 0.4706; 0.6471 0.5686 0.5255; 0.8471 0.7922 0.6157; 0.8 0.8 1; 0.4 0.4 1; 0.7882 0.5686 0.2235];    % Planet Colours
planet_colors(13:end, :) = repmat([0.5 0.5 0.5], numel(13:length(body_names)), 1);                                                      % Satelite Colours

% Body Properties: Radius (m), Mass (kg), Gravitational Parameter (m^3/s^2), SOI (m)
body_prop{1} = [696000000, 1.9885e30, 0, 0];    % Sol (REF=1)
body_prop{2} = [2439400, 0.330103e24, 0, 0];    % Mercury (REF=2)
body_prop{3} = [6051800, 4.86731e24, 0, 0];     % Venus (REF=3)
body_prop{4} = [6371008.4, 5.97217e24, 0, 0];   % Earth (REF=4)
body_prop{5} = [3389500, 0.641691e24, 0, 0];    % Mars (REF=5)
body_prop{6} = [3389500, 9.38416e20, 0, 0];     % Ceres (REF=6)
body_prop{7} = [285000, 2.59027e20, 0, 0];      % Vesta (REF=7)
body_prop{8} = [69911000, 1898.125e24, 0, 0];   % Jupiter (REF=8)
body_prop{9} = [58232000, 568.317e24, 0, 0];    % Saturn (REF=9)
body_prop{10} = [25362000, 86.8099e24, 0, 0];   % Uranus (REF=10)
body_prop{11} = [24622000, 102.4092e24, 0, 0];  % Neptune (REF=11)
body_prop{12} = [1188300, 13029.0e18, 0, 0];    % Pluto (REF=12)
body_prop{13} = [1737400, 0, 4.9028e12, 0];     % Luna (REF=13)
body_prop{14} = [11080, 0, 7.087e5, 0];         % Phobos (REF=14)
body_prop{15} = [6200, 0, 9.62e4, 0];           % Deimos (REF=15)
body_prop{16} = [1821490, 0, 5.9599155e12, 0];  % IO (REF=16)
body_prop{17} = [1560800, 0, 3.2027121e12, 0];  % Europa (REF=17)
body_prop{18} = [2631200, 0, 9.8878328e12, 0];  % Ganymede (REF=18)
body_prop{19} = [2410300, 0, 7.1792834e12, 0];  % Callisto (REF=19)
body_prop{20} = [252100, 0, 7.21037e9, 0];      % Enceladus (REF=20)
body_prop{21} = [531100, 0, 4.1213453e10, 0];   % Tethys (REF=21)
body_prop{22} = [561400, 0, 7.311607e10, 0];    % Dione (REF=22)
body_prop{23} = [763500, 0, 1.5394175e11, 0];   % Rhea (REF=23)
body_prop{24} = [2574760, 0, 8.9781371e12, 0];  % Titan (REF=24)
body_prop{25} = [235800, 0, 4.3e9, 0];          % Miranda (REF=25)
body_prop{26} = [578900, 0, 8.35e10, 0];        % Ariel (REF=26)
body_prop{27} = [584700, 0, 8.51e10, 0];        % Umbriel (REF=27)
body_prop{28} = [788900, 0, 2.269e11, 0];       % Titania (REF=28)
body_prop{29} = [761400, 0, 2.053e11, 0];       % Oberon (REF=29)
body_prop{30} = [1352600, 0, 1.4284955e12, 0];  % Triton (REF=30)

% GP from Mass Calculation
for i = 1:12
    body_prop{i}(1,3) = body_prop{i}(1,2) * G;
end

% Mass from GP Calculation
for i = 13:num_bodies
    body_prop{i}(1,2) = body_prop{i}(1,3) / G;
end

% Orbital Characteristics - Initial Conditions
% Reference Body (Parent Body), SMA, Ecc, Inc, AoP, RAAN, MNA, TA (Calculated Later)
orbits{1} = [1, 0, 0, 0, 0, 0, 0, 0];
for i = 2:12
    orbits{i}(1,1) = 1;
end
% Specifying Parent Bodies
orbits{13}(1,1) = 4;
orbits{14}(1,1) = 5;
orbits{15}(1,1) = 5;
orbits{16}(1,1) = 8;
orbits{17}(1,1) = 8;
orbits{18}(1,1) = 8;
orbits{19}(1,1) = 8;
orbits{20}(1,1) = 9;
orbits{21}(1,1) = 9;
orbits{22}(1,1) = 9;
orbits{23}(1,1) = 9;
orbits{24}(1,1) = 9;
orbits{25}(1,1) = 10;
orbits{26}(1,1) = 10;
orbits{27}(1,1) = 10;
orbits{28}(1,1) = 10;
orbits{29}(1,1) = 10;
orbits{30}(1,1) = 11;

start_time = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
end_time = start_time + (1/(365*24));
orbits_cells = cell(4,3);
for i = 2:num_bodies
    orbits_url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text";
    orbits_url = orbits_url + "&COMMAND='" + string(body_names{2,i}) + "'";
    orbits_url = orbits_url + "&OBJ_DATA='NO'";
    orbits_url = orbits_url + "&Make_EPHEM='YES'";
    orbits_url = orbits_url + "&EPHEM_TYPE='ELEMENTS'";
    orbits_url = orbits_url + "&CENTER=" + "'" + string(body_names{2,orbits{i}(1,1)}) + "'";
    orbits_url = orbits_url + "&START_TIME=" + "'" + string(start_time) + "'";
    orbits_url = orbits_url + "&STOP_TIME=" + "'" + string(end_time) + "'";
    orbits_data = strtrim(webread(orbits_url));
    orbits_lines = strsplit(orbits_data,'\n');

    for j = 1:length(orbits_lines)
        if ismember(orbits_lines{j},"$$SOE")
            for k = j+2:j+5
                temp_cell = strsplit(orbits_lines{k},'= ');
                orbits_cells{k-j-1,1} = temp_cell{1,2}(1:end-3);
                orbits_cells{k-j-1,2} = temp_cell{1,3}(1:end-3);
                orbits_cells{k-j-1,3} = temp_cell{1,4};
            end
        end
    end
    orbits{i}(1,2) = 1000*str2double(orbits_cells{4,1});
    orbits{i}(1,3) = str2double(orbits_cells{1,1});
    orbits{i}(1,4) = d2r*str2double(orbits_cells{1,3});
    orbits{i}(1,5) = d2r*str2double(orbits_cells{2,2});
    orbits{i}(1,6) = d2r*str2double(orbits_cells{2,1});
    orbits{i}(1,7) = d2r*str2double(orbits_cells{3,2});
    orbits{i}(1,8) = d2r*str2double(orbits_cells{3,3});
end

% Patched Conics SOI calculation

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
sim_time = linspace(0, 500000, 500001); % Simulation length in seconds
sim_percent = (length(sim_time)-1)/1000; % Simulation percentage constant
traj_points = 250000; % Position/velocity stored at once

% Position Matrix Preallocation
for i = 1:num_bodies
    planet_position{i} = zeros(traj_points, 3);
end

% Velocity Matrix Preallocation
for i = 1:num_bodies
    planet_velocity{i} = zeros(traj_points, 3);
end

% Initial Position and Velocity in KoCI (Kerbol Centered Inertial)
planet_position{1}(1, :) = [0,0,0]; % Sol's Initial Position in SCI
planet_velocity{1}(1, :) = [0,0,0]; % Sol's Initial Velocity in SCI
for i = 2:num_bodies
    C_PS = create_rot_mtx(3,orbits{i}(1,5))*create_rot_mtx(1,orbits{i}(1,4))*create_rot_mtx(3,orbits{i}(1,6));
    C_SP = transpose(C_PS);
    planet_position{i}(1, :) = transpose(C_SP*r_P{i}) + planet_position{orbits{i}(1,1)}(1,:); % Initial Position in SCI
    planet_velocity{i}(1, :) = transpose(C_SP*v_P{i}) + planet_velocity{orbits{i}(1,1)}(1,:); % Initial Velocity in SCI
end


% Figure Setup
figureSize = [100, 100, 1225, 975];        % Pixels
tiledlayout(3, 4, 'TileSpacing','compact'); % 3 rows, 4 columns
tl = gcf;                                   % Get the current figure handle
tl.OuterPosition = figureSize;              % Sets figure to specified size

% Create a VideoWriter object
movieFileName = 'Solar_System_Simulation.mp4';    % Output file name
movieWriter = VideoWriter(movieFileName, 'MPEG-4'); % Defines movieWriter
movieWriter.FrameRate = 60;                         % Set the desired frame rate
spf = 1000;                                         % Seconds simulated per frame
open(movieWriter);                                  % Opens movieWriter

% In-Plane Plot (Tile 1) - Inner Solar System
nexttile([2, 2])
ax1 = gca;
ax1.Color = 'k'; % Set background color to black
ax1.GridColor = 'k'; % Set grid color to white
hold on;
axis(ax1, 'equal');
axis(ax1, [-5.3e11 5.3e11 -5.3e11 5.3e11]);
title(ax1, 'In-Plane Motion');
xlabel(ax1, 'X Position (meters)');
ylabel(ax1, 'Y Position (meters)');
grid(ax1, 'on');

% In-Plane Plot (Tile 2) - Outer Solar System
nexttile([2,2])
ax2 = gca;
ax2.Color = 'k'; % Set background color to black
ax2.GridColor = 'k'; % Set grid color to white
hold on;
axis(ax2, 'equal');
axis(ax2, [-7.5e12 7.5e12 -7.5e12 7.5e12]);
title(ax2, 'In-Plane Motion');
xlabel(ax2, 'X Position (meters)');
grid(ax2, 'on');

% Out-of-Plane Plot (Tile 3) - Inner Solar System
nexttile([1, 2])
ax3 = gca;
ax3.Color = 'k'; % Set background color to black
ax3.GridColor = 'k'; % Set grid color to white
hold on;
axis(ax3, 'equal');
axis(ax3, [-5.3e11 5.3e11 -1.855e11 1.855e11]);
title(ax3, 'Out-of-Plane Motion');
xlabel(ax3, 'X Position (meters)');
ylabel(ax3, 'Z Position (meters)');
grid(ax3, 'on');

% Out-of-Plane Plot (Tile 4) - Outer Solar System
nexttile([1, 2])
ax4 = gca;
ax4.Color = 'k'; % Set background color to black
ax4.GridColor = 'k'; % Set grid color to white
hold on;
axis(ax4, 'equal');
axis(ax4, [-7.5e12 7.5e12 -2.625e12 2.625e12]);
title(ax4, 'Out-of-Plane Motion');
xlabel(ax4, 'X Position (meters)');
grid(ax4, 'on');

%%
% Simulation Loop
multiple = 0;
for t = 1:length(sim_time)
    if mod(t,sim_percent)==0
        disp('Progress: ' + string(round((100*t)/length(sim_time),2))+'%'); % Displays current progress
    end
    force_array = cell(num_bodies,num_bodies); 

    % Indexing Logic
    if t < traj_points
        k_cur = t;
        k_next = k_cur + 1; 
    elseif mod(t,traj_points) == 0
        multiple = multiple + 1;
        k_cur = traj_points;
        k_next = 1;
    else 
        k_cur = t - multiple*traj_points;
        k_next = k_cur + 1;
    end

    % Force Loop
    for i = 1:num_bodies % Current Body
        for j = i:num_bodies % Interacting Body
            if i ~= j % Ignores interaction between current body and itself
                % Force vector calculation
                force_ij_mag = (G*body_prop{i}(1, 2)*body_prop{j}(1, 2)) / (norm(planet_position{j}(k_cur, :) - planet_position{i}(k_cur, :))^2);
                force_ij_vect = force_ij_mag * ((planet_position{j}(k_cur, :) - planet_position{i}(k_cur, :)) / norm(planet_position{j}(k_cur, :) - planet_position{i}(k_cur, :)));
                force_array{i,j} = force_ij_vect;
            end
        end
    end
    for i = 1:length(force_array)
        force_array{i,i} = [0,0,0]; % Sets major diagonal to [0,0,0]
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
        planet_velocity{i}(k_next, :) = planet_velocity{i}(k_cur, :) + a_i;
        planet_position{i}(k_next, :) = planet_position{i}(k_cur, :) + planet_velocity{i}(k_cur, :);
    end
    if mod(t, spf) == 0
        delete(planetMarkers_1);
        delete(planetMarkers_2);
        delete(planetMarkers_3);
        delete(planetMarkers_4);

        delete(planetTrajectories_1);
        delete(planetTrajectories_2);
        delete(planetTrajectories_3);
        delete(planetTrajectories_4);
       
        delete(planetLabels_1);
        delete(planetLabels_2);
        delete(planetLabels_3);
        delete(planetLabels_4);

        delete(timeCounter);

        if multiple ~= 0
            delete(planetTrajectories_1_2);
            delete(planetTrajectories_2_2);
            delete(planetTrajectories_3_2);
            delete(planetTrajectories_4_2);
        end

        for i = 1:num_bodies % Planets
            if orbits{i}(1,1) == 1 || norm(planet_position{i}(k_next,:) - planet_position{orbits{i}(1,1)}(k_next,:)) > body_prop{orbits{i}(1,1)}(1,4)
                % Inner Solar System In-Plane Motion
                xlim(ax1, [planet_position{1}(k_next, 1) - 5.3e11, planet_position{1}(k_next, 1) + 5.3e11]);
                ylim(ax1, [planet_position{1}(k_next, 2) - 5.3e11, planet_position{1}(k_next, 2) + 5.3e11]);

                % Outer Solar System In-Plane Motion
                xlim(ax2, [planet_position{1}(k_next, 3) - 7.5e12, planet_position{1}(k_next, 3) + 7.5e12]);
                ylim(ax2, [planet_position{1}(k_next, 2) - 7.5e12, planet_position{1}(k_next, 2) + 7.5e12]);

                % Inner Solar System Out-of-Plane Motion
                xlim(ax3, [planet_position{1}(k_next, 1) - 5.3e11, planet_position{1}(k_next, 1) + 5.3e11]);
                ylim(ax3, [planet_position{1}(k_next, 2) - 1.855e11, planet_position{1}(k_next, 2) + 1.855e11]);

                % Outer Solar System Out-of-Plane Motion
                xlim(ax4, [planet_position{1}(k_next, 3) - 7.5e12, planet_position{1}(k_next, 3) + 7.5e12]);
                ylim(ax4, [planet_position{1}(k_next, 2) - 2.625e12, planet_position{1}(k_next, 2) + 2.625e12]);

                if t <= traj_points
                    if i == 1
                        planetMarkers_1(i) = plot(ax1, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_2(i) = plot(ax2, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_3(i) = plot(ax3, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_4(i) = plot(ax4, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                    elseif i <= 7 && i ~= 1
                        planetTrajectories_1(i) = plot(ax1, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_3(i) = plot(ax3, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetMarkers_1(i) = plot(ax1, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_3(i) = plot(ax3, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetLabels_1(i) = text(ax1, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                        planetLabels_3(i) = text(ax3, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                    elseif i <= 12 && i ~= 1
                        planetTrajectories_2(i) = plot(ax2, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_4(i) = plot(ax4, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetMarkers_2(i) = plot(ax2, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_4(i) = plot(ax4, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetLabels_2(i) = text(ax2, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                        planetLabels_4(i) = text(ax4, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                    end
                else
                    if i == 1
                        planetMarkers_1(i) = plot(ax1, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_2(i) = plot(ax2, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_3(i) = plot(ax3, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_4(i) = plot(ax4, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                    elseif i <= 7 && i ~= 1
                        planetTrajectories_1(i) = plot(ax1, planet_position{i}(k_next+1:1000:traj_points, 1), planet_position{i}(k_next+1:1000:traj_points, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_3(i) = plot(ax3, planet_position{i}(k_next+1:1000:traj_points, 1), planet_position{i}(k_next+1:1000:traj_points, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_1_2(i) = plot(ax1, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_3_2(i) = plot(ax3, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetMarkers_1(i) = plot(ax1, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_3(i) = plot(ax3, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetLabels_1(i) = text(ax1, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                        planetLabels_3(i) = text(ax3, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                    elseif i <= 12 && i ~= 1
                        planetTrajectories_2(i) = plot(ax2, planet_position{i}(k_next+1:1000:traj_points, 1), planet_position{i}(k_next+1:1000:traj_points, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_4(i) = plot(ax4, planet_position{i}(k_next+1:1000:traj_points, 1), planet_position{i}(k_next+1:1000:traj_points, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_2_2(i) = plot(ax2, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 2), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetTrajectories_4_2(i) = plot(ax4, planet_position{i}(1:1000:k_next, 1), planet_position{i}(1:1000:k_next, 3), 'LineWidth', 1, 'Color', planet_colors(i, :));
                        planetMarkers_2(i) = plot(ax2, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetMarkers_4(i) = plot(ax4, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), 'o', 'MarkerSize', 5, 'Color', planet_colors(i, :), 'MarkerFaceColor', 'none');
                        planetLabels_2(i) = text(ax2, planet_position{i}(k_next, 1), planet_position{i}(k_next, 2), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                        planetLabels_4(i) = text(ax4, planet_position{i}(k_next, 1), planet_position{i}(k_next, 3), body_names{1,i},'Color', planet_colors(i, :), 'FontSize', 8, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
                    end
                end      
            end
        end
        timeCounter = annotation('textbox', [0, 0, 1, 0.05], 'String', ['Time: ' num2str(t) ' s elapsed'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'BackgroundColor', 'w');

        % Capture the current frame and write to the video file

        drawnow
        frame = getframe(gcf);
        
        writeVideo(movieWriter, frame);
    end
end
close(movieWriter);

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