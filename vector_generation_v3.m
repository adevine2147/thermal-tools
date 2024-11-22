clc
clear
%------------------------------BACKGROUND----------------------------------
% the purpose of this code is do generate vectors to input into TD for
% accurate environmental heating. The user shall input all necessary 
% information about the orbit and the code will produce vectors to paste 
% into TD for analysis. This was originally created for SN2, orientations
% may have to be edited for other vehicles
% time in seconds, length in km
% Time zero is a quarter orbit before subsolar!
%
% ---------------------ASSUMPTIONS/ATTITUDE CONSTRAINTS-------------------
% slewing takes roughly 3 minutes ~GNC. this code won't capture slewing
% since positions are 90min/15positions ~6 minutes apart 
% SC appears to orbit earth in CCW according to TD
%
% SN2 has Z+ aligned with aperture
%         Y+ aligned with avionics tray
%         X+ aligned with torque rods side
% in hot case , SC X+ is pointing velocity (avionics tray facing nadir)
% in cold case, SC y+ is pointing velocity (avionics tray facing space)
%
% SUN POINTING, Spacecraft X- is aligned with
% geocentric Z+, and spacecraft Z- is aligned with geocentric Y+
%
% NADIR POINTING, Spacecraft X- is aligned with velocity, spacecraft Z+
% points towards GEOCENTRIC [0 0 0]

%----------------------------PROBLEM SETUP---------------------------------

h = 500; % altitude of spacecraft in kms
beta = 20*-1; % beta angle in DEGREES (double check definition, maybe change Z definition)
n_positions = 15; % number of positions minus 1 in orbit desired
nadir_time_fraction = 0.25; % fraction of orbit when spacecraft is pointed nadir
nadir_phase_shift = 0.7; % how long after beginning of orbit the spacecraft starts to point nadir
% 0.25
% 0.7

mu = 398600; %gravitational parameter of earth
R = 6378; %radius of earth km
a = R+h; %altitude from center of earth

% orbital period, and time at each position
T_orbit = floor(2*pi*mu^(-0.5)*a^1.5);
times = linspace(0,T_orbit,n_positions);
time_step = T_orbit/n_positions;
time_step_minutes = time_step/60;

% eclipse fraction/beta star, spacecraft thermal control handbook pg 41
beta_star = asind(R/(R+h)); 
if (beta < beta_star)
    eclipse_fraction = acosd(sqrt(h^2 + 2*R*h)/((R+h)*cosd(beta)))/180;
    t_eclipse_enter = T_orbit/2-eclipse_fraction*T_orbit/2;
    t_eclipse_exit = t_eclipse_enter + eclipse_fraction*T_orbit; 
else
    eclipse_fraction = 0; 
    t_eclipse_enter = 0;
    t_eclipse_exit = 0; 
end
% choosing what positions the spacecraft will be pointing towards nadir
if(nadir_time_fraction>0)
    n_nadir_positions = floor(floor(nadir_time_fraction*n_positions)); 
    nadir_positions = floor(nadir_phase_shift*n_positions):(floor(nadir_phase_shift*n_positions)+n_nadir_positions-1);
end

% forgot what i had in my other code, double check magnitudes
% phase shift of +pi/2 to start time at subsolar point
x_pos = @(t) a*cos(t*2*pi/T_orbit + pi/2); 
y_pos = @(t) a*sin(t*2*pi/T_orbit + pi/2)*cosd(beta); 
z_pos = @(t) a*sin(t*2*pi/T_orbit + pi/2)*sind(beta);

position_array = zeros(n_positions, 3); 

% creating cartesian positions, GCS origin (0,0,0) is center of earth
% [0,1,0] points to sun
for position=1:n_positions
        position_array(position,1) = x_pos(times(position));% x
        position_array(position,2) = y_pos(times(position));% y
        position_array(position,3) = z_pos(times(position));% z
end

% 1st column is times, 2-4 is sun vector, 5-7 is nadir vector, 8 is alt
% ratio
TD_output_array = zeros(n_positions,8);
TD_output_array(:,1) = times; 
TD_output_array(:,end) = a/R;

%Pointing arrays from n_positions to sun and nadir
geo_sun_point_array = repmat([0 1 0], n_positions, 1);

geo_nadir_point_array = zeros(n_positions,3);
for position=1:n_positions
     vec = [0 0 0] - position_array(position,:);
     geo_nadir_point_array(position,1) = vec(1)/norm(vec);
     geo_nadir_point_array(position,2) = vec(2)/norm(vec);
     geo_nadir_point_array(position,3) = vec(3)/norm(vec);
end

%-------------------------------MATH---------------------------------------
% math in radians, beta input in degrees

% basis vectors for geocentric coords, no touch
x0 = [1; 0; 0];
y0 = [0; 1; 0]; 
z0 = [0; 0; 1];
GCS_Basis_Matrix = [x0,y0,z0];

% vector normal to orbital plane, used for finding velocity vector
vno = [0, -sind(beta), cosd(beta)];

% creating all velocity vectors for each position and normalizing
% roundoff errors? 
velocities = zeros(n_positions, 3);
for position=1:n_positions
     velocity = cross(vno,geo_nadir_point_array(position,:))*-1;
     velocities(position,1) = velocity(1)/norm(velocity);
     velocities(position,2) = velocity(2)/norm(velocity);
     velocities(position,3) = velocity(3)/norm(velocity);
end

% basis matrices used for rotation matrices later
% n_positions pages of 3x3 matrices
% row 1 is x basis vector
% row 2 is y basis vector 
% row 3 is z basis vector
SCCS_basis_matrices = zeros(3,3,n_positions);

% looping through sun pointing positions and creating their basis
% matrices
%   Z- of SC points Y+ of GCS (solar array)
%   Y- of SC is pointed Z+ in GCS (heavens)
for position=1:n_positions
    SCCS_basis_matrices(2,:,position) = -1*[0, 0, 1]; % heavens
    SCCS_basis_matrices(3,:,position) = -1*[0, 1, 0]; % sun
    SCCS_basis_matrices(1,:,position) = cross( ...
        SCCS_basis_matrices(2,:,position), ...
        SCCS_basis_matrices(3,:,position));
end

% looping through nadir pointing positions and creating their basis
% matrices
%   Z+ of SC points at GCS origin (0,0,0)
%   X+ of SC points velocity (same plane as GCS nadir vector, 90 degrees
if(nadir_time_fraction>0)
    for position=nadir_positions(1):nadir_positions(end)
        SCCS_basis_matrices(1,:,position) = velocities(position,:);
        SCCS_basis_matrices(3,:,position) = (geo_nadir_point_array(position,:));
        SCCS_basis_matrices(2,:,position) = cross( ...
            SCCS_basis_matrices(1,:,position), ...
            SCCS_basis_matrices(3,:,position))*-1;
    end
end

% creating sun and nadir pointing vectors for TD
% 1 is time
% 234 is xyz for sun
% 567 is xyz for earth
% 8 is altitude ratio
for position=1:n_positions
    TD_output_array(position,2:4) = SCCS_basis_matrices(:,:,position)*transpose(geo_sun_point_array(position,:));
    TD_output_array(position,5:7) = SCCS_basis_matrices(:,:,position)*transpose(geo_nadir_point_array(position,:));
end


%-------------------------------VISUALIZATION------------------------------

axis equal;
hold on;
alpha(1);

shading flat

% orbital normal
quiver3(0,0,0, vno(1), vno(2), vno(3), 1000);
% earth to sun
quiver3(0,0,0, 0, 5000, 0, 'y', 'LineWidth',2, 'MaxHeadSize',1);
%plotting all nadir and sun pointing vectors from each position
model = stlread('aries.stl'); 
modelScaleFactor = 20;

for position=1:n_positions
    translation = position_array(position,:);
    rotation = SCCS_basis_matrices(:,:,position); 
    pos1 = model.Points*modelScaleFactor*rotation + translation; 
trisurf(model.ConnectivityList, pos1(:,1), pos1(:,2), pos1(:,3), ...
    'FaceColor', 'cyan', 'EdgeColor', 'none');

    %sun pointing
    % quiver3(position_array(position,1),position_array(position,2), ...
    %     position_array(position,3),geo_sun_point_array(position, 1), ...
    %     geo_sun_point_array(position, 2),geo_sun_point_array(position, 3),2000)
    % velocities
    % quiver3(position_array(position,1),position_array(position,2), ...
    %     position_array(position,3), velocities(position, 1),velocities ...
    %     (position, 2),velocities(position, 3),2000)
    %nadir pointing
    % quiver3(position_array(position,1),position_array(position,2), ...
    %     position_array(position,3),geo_nadir_point_array(position,1), ...
    %     geo_nadir_point_array(position,2), geo_nadir_point_array(position,3) )

    scale = 2000; 
    % spacecraft X, red 
    quiver3(position_array(position,1),position_array(position,2), ...
        position_array(position,3), ...
        SCCS_basis_matrices(1,1,position)*scale, ...
        SCCS_basis_matrices(1,2,position)*scale,  ...
        SCCS_basis_matrices(1,3,position)*scale, ...
        'r','LineWidth', 2, 'MaxHeadSize', 1);
    % spacecraft Y, green
    quiver3(position_array(position,1),position_array(position,2), ...
        position_array(position,3), ...
        SCCS_basis_matrices(2,1,position)*scale, ...
        SCCS_basis_matrices(2,2,position)*scale,  ...
        SCCS_basis_matrices(2,3,position)*scale, ...
        'g', 'LineWidth', 2, 'MaxHeadSize', 1);
    % spacecraft Z, blue
    quiver3(position_array(position,1),position_array(position,2), ...
        position_array(position,3), ...
        SCCS_basis_matrices(3,1,position)*scale, ...
        SCCS_basis_matrices(3,2,position)*scale,  ...
        SCCS_basis_matrices(3,3,position)*scale, ...
        'b','LineWidth',2, 'MaxHeadSize', 1);
    % quiver3(position_array(position,1),position_array(position,2),position_array(position,3), [u v w])
end

% 1st position dot
plot3(position_array(1,1),position_array(1,2),position_array(1,3),'-o','Color','b','MarkerSize',10,...
   'MarkerFaceColor','#D9FFFF')
% 1st position SCCS, Z- faces sun

%plotting orbital path
plot3(position_array(:,1),position_array(:,2),position_array(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
plot3(0,0,0,'-o','Color','b','MarkerSize',20,...
    'MarkerFaceColor','#008000')

% sun messed up the scales, maybe put in later
% plot3(0,40000,0,'-o','Color','b','MarkerSize',10,...
%     'MarkerFaceColor','#FFFF00')

