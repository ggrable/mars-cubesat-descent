% Mars Descent Simulaltion

% This code simulates, animates, analizes impact force and saves a video file
% of a custom descent of a lander on Mars. In the animation, blue represents
% parachute deployment, green represents impact survival and red represents
% impact accident.

% Rodrigo Santos Valente da Costa
% RodrigoS.V.daCosta@gmail.com
% June 2016

clear; close all;

% Lander and Planet Settings
L1= 1000; % Size of the lander in the animation (does not affect calculations), m
g= 3.711; % Acceleratiojn of gravity, m/s^2
A= 1; % Projected area of the cubesat, m^2
Cd= 1.05; % Coefficient of drag (1.05 for a cube), unitless
m= 5; % Mass of the lander, kg
rho= 0.0044; % Initial density of the atmosphere, kg/m^3
drag= 0; % Initial drag force, N

mars_mass= 6.39*10^23; % Mass of mars, kg
G= 6.67*10^(-11); % Gravitational constant
mars_radius= 3.39 * 10^6; % Mars radius, m

para_altitude= 12000; % Altitude of parachute deployment, m


% Initial Conditions
y0= 15000; % Initial altitude, m
v0= -975; % Initial vertical velocity, m/s

yAcc= 0; % Initial acceleration, m/s^2
yVel= v0; % Vertical velocity, m/s
y= y0; % Altitude, m

ydata= [];
yVeldata= [];
yAccdata= [];
time= [];
density= [];
gravity= [];
dragdata= [];

frame_size= y0 + L1; % External frame of the animation

% Time Settings
duration= 690; % Duration of simluation, s
delta_t= 0.01; % Delta time of the simulation, s
time_steps= round(duration/delta_t); % Number of time steps
frameRate= 2; % Animation frame rate, frames/s
framePeriod= 1./frameRate; % Period of one frame, s

% Draw initial figure
fig= figure(1);
axs=axes('Parent',fig);

% Open movie file
vidObj= VideoWriter('mars_descent_simulation.avi');
vidObj.FrameRate= frameRate;
open(vidObj);

% Draw Lander
cubesat= rectangle('Position', [0, -(frame_size - 1)*L1, L1, L1], 'Curvature', [0,0]);
border= line([-frame_size/2, -frame_size/2, frame_size/2, frame_size/2, -frame_size/2],...
    [0, frame_size, frame_size, 0, 0]);
axis(axs, [-frame_size/2, frame_size/2, 0, frame_size]);
axis(axs, 'square');
axis(axs, 'on'); % Turns axes values on
set(axs,'XTick',[]); % Turns x axis off
ylabel('Altitude (m)');
pressure_change= line([-frame_size, frame_size], [7000,7000], 'Color', 'red'); ...
    % Draws a line where the atmospheric density changes

% Video generation
for i= 1:time_steps
    t= (i-1)*delta_t;
    
    % Calculation of the gravity acceleration depending on the altitude
    g= G*mars_mass/(y + mars_radius)^2;
    
    if y<=0 % Stops when it hits the ground
        break;
    end
    
    if y>para_altitude
        set(cubesat, 'Position', [-L1/2, y, L1, L1], 'FaceColor', 'black');
        Cd= 1.05; % Coefficient of drag of Cubesat, unitless
        A= 0.01; % Surface area of the Cubesat only
    end 
    if y<=para_altitude
        set(cubesat, 'Position', [-L1/2, y, L1, L1], 'FaceColor', 'blue');
        Cd= 1.3; % Coefficient of drag of parachute, unitless
        A= 2; % Surface area of the parachute, m^2
    end     
    
    % Atmospheric density calculation depending on the altitude
    if y>7000
        temp=-23.4 - 0.00222*y;
        pressure= 0.699*exp(-0.00009*y);
    end
    if y<7000
        temp= -31 - 0.000998*y;
        pressure= 0.699*exp(-0.00009*y);
    end
    rho= pressure/(.1921*(temp+273.1));
    
    % Resettig video data and getting position/velocity
    if (mod(t,framePeriod) == 0)
        currFrame= getframe;
        writeVideo(vidObj, currFrame);
        time= [time t];
        ydata= [ydata y];
        yVeldata= [yVeldata yVel];
        yAccdata= [yAccdata yAcc];
        density= [density rho];
        gravity= [gravity g];
    end
        
    drag= A/2*Cd*rho*yVel^2; % Drag force, N
    
    yAcc= (drag - m*g)/m; % Acceleration considering drag, m/s^2
    
    yVel= yVel + yAcc*delta_t; % Velocity considering drag, m/s
    
    y= y + yVel*delta_t; % Altitude considering drag, m
end

% Impact calculation
s= 0.001; % Slow down distance after impact, m
impactForce= 1/2*m*yVel^2/s; % Impact force, N
impactStress= impactForce/0.01; % Impact stress, Pa

if impactStress<276*10^6
    set(cubesat, 'Position', [-L1/2, y, L1, L1], 'FaceColor', 'green');
    disp('Cubesat survived the land.');
end
if impactStress>=276*10^6
    set(cubesat, 'Position', [-L1/2, y, L1, L1], 'FaceColor', 'red');
    disp('Cubesat did not survive the land.');
end

elapsedTime= t;
disp('The descent and landing process took (in seconds):');
disp(elapsedTime);

%Close the movie file.
close(vidObj);

% Plot Time vs. Altitude
figure(2);
plot(time, ydata);
xlabel('Time (s)');
ylabel('Altitude (m)');
line([0, t], [7000,7000], 'Color', 'red');

% Plot Time vs. Velocity
figure(3);
plot(time, yVeldata);
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Plot Time vs. Acceleration
figure(4)
plot(time, yAccdata);
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

% Plot Time vs. Atmospheric Density
figure(5)
plot(time, density);
xlabel('Time (s)');
ylabel('Air density (kg/m^3)');

% Plot Time vs. Gravity
figure(6)
plot(time, gravity);
xlabel('Time (s)');
ylabel('Gravity Acceleration (m/s^2)');