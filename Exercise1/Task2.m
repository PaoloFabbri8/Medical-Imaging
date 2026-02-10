% TASK 2
% Implement the above formula in MATLAB or Python. Consider a spherical head with a radius of 10 cm. 
% If there is a single current dipole at ⃗r0 = (−7.5, 0, 0) cm with dipole moment Q = (0, 1, 0) · 50 nAm,
%  what is the measured magnetic field ⃗B at 
% i) the position ⃗r1 = (−10.1, 0, 2.7) cm?
% ii) the position ⃗r2 = (−12.0, 0, 3.2) cm?

% r1 roughly corresponds to the measurement distance of optically-pumped magnetometers (OPMs), 
% while ⃗r2 corresponds to that of SQUID magnetometers, which have more thermal insulation.
%Plot the magnitudes of the magnetic field component normals to the surface of the sphere on a spherical surface at radius 11 cm, 
% i.e. just outside the sphere. Return your code as well as the plot.


% Creating a function to compute thenmagnetic field (B)
%Dipole parameters
Q = [0; 1; 0] * 50e-9;
r0 = [-7.5; 0; 0] * 1e-2;

% Positions
r1 = [-10.1; 0; 2.7] * 1e-2;
r2 = [-12.0; 0; 3.2] * 1e-2;

% Constant
mu0 = 4*pi*1e-7;

B_r1 = sarvas_function(r1, r0, Q, mu0);
B_r2 = sarvas_function(r2, r0, Q, mu0);


disp('Magnetic field at r1 (OPM):');
disp(B_r1)

disp('Magnetic field at r2 (SQUID):');
disp(B_r2)


% Radius of the sphere
R = 0.11;

% Sphere creation
theta = linspace(0, pi, 50);        
phi   = linspace(0, 2*pi, 100);    
[Theta,Phi] = meshgrid(theta, phi);

X = R*sin(Theta).*cos(Phi);
Y = R*sin(Theta).*sin(Phi);
Z = R*cos(Theta);

% computation of B normal values for each point in the sphere
B_normal = zeros(size(X));

for i = 1:numel(X)
    r_point = [X(i); Y(i); Z(i)];  
    B = sarvas_function(r_point, r0, Q, mu0);  
    B_normal(i) = dot(B, r_point)/norm(r_point); 
end

% Plot 3D
figure;
surf(X, Y, Z, B_normal, 'EdgeColor', 'none');  
colormap jet;
colorbar;
axis equal;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Normal component of B on a sphere at 11 cm radius');