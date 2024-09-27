% Parameters
clear all;
rho = 8.238e3;      % Material density [kg/m^3]
Cp = 4.68e2;        % Heat capacity [J/kg*K]
Kappa = 1.34e1;     % Thermal conductivity [W/m*K]
alpha = 5.45e6;     % Absorption coefficient [m^-1]
Pav = 30;           % Average Power [W]
f = 1.02e5;         % Laser pulse frequency [Hz]
tau = 1.4e-7;       % Laser pulse duration [s]
w0 = 4e-5;          % Laser beam spot ∅ [m]
Tm = 1.670e3;          % Melting point [K]
Tb = 3.173e3;          % Boiling point [K]
Lf = 3.0e5;        % Latent heat of melting (kJ/kg)
melting_mask = false(200,200);  % Mask to apply temperature adjustment
melting_mask1 = false(200,200); % Mask to store last applied adjustment

R0 = w0 / 2;        % Radius of the laser beam spot at the focus [m]
I0 = Pav / (f * tau * pi * R0^2);   % Peak on-line laser intensity [W/m^2]

% Spatial and temporal discretization
dx = 2e-6;          % Discretization along x axis [m]
dz = 4e-7;          % Discretization along z axis [m]
dt1 = 0.5*(rho*Cp*dx^2*dz^2)/(Kappa*(dx^2+dz^2));    % Time step (s)
dt2 = tau/10;
dt = min([dt1 dt2]);

Lx = 4e-4;          % Length of the material in x [m]
Lz = 8e-5;          % Length of the material in z [m]

Nx = round(Lx / dx);        % Number of spatial grid points in X
Nz = round(Lz / dz);        % Number of spatial grid points in Z
Nt = round(tau / dt);       % Number of time steps

% Boundry Conditions
T0 = 293.15;        % Initial Temperature [K]
Tx = 293.15;        % Temperature Left and Right [K]
Tz = 293.15;        % Temperature Bottom [K]

% Initialize temperature field
T = T0 * ones(Nz, Nx);

% Transformation Variables
    x_pulse = Lx / 2;
    t_pulse = Nt*dt*0.55;

% Simulation loop
for n = 0:Nt
% Calculate heat source term (pulse at the center)
    A = alpha * I0 * exp(-((((1:Nx) * dx) - x_pulse) / R0) .^2)...
        * exp(-(((n * dt) - t_pulse) / tau)^2);
    for i = 1:Nz-1
        depth_decay = exp(-alpha*(i-1)*dz);
        for k = 2:Nx-1
            if i==1
                T(i, k) = T(i, k)...
                    + (((dt * Kappa) / (rho * Cp)) * ((T(i+1, k) - 2 * T(i, k) + T0) / (dz ^ 2) + (T(i, k+1) - 2 * T(i, k) + T(i, k-1)) / (dx ^ 2)))...
                    + (dt / (rho * Cp)) * A(1, k) * depth_decay;
            else
                T(i, k) = T(i, k)...
                    + (((dt * Kappa) / (rho * Cp)) * ((T(i+1, k) - 2 * T(i, k) + T(i-1, k)) / (dz ^ 2) + (T(i, k+1) - 2 * T(i, k) + T(i, k-1)) / (dx ^ 2)))...
                    + (dt / (rho * Cp)) * A(1, k) * depth_decay;
            end
        end
    end

% Boundary conditions (fixed temperature at boundaries)
    T(end, :) = Tz;
    T(:, 1) = Tx;
    T(:, end) = Tx;

% Adjust temperature for melting
    melting_mask = T > Tm;
    melting_mask =logical(melting_mask-melting_mask1);
    T = T - melting_mask * (Lf / Cp); 
    melting_mask1 = logical(melting_mask1 + melting_mask);

% Plot temperature distribution every time step
        figure;
        imagesc((1:Nx)*dx/1e-6, (1:Nx)*dz/1e-6, T)
        title(['Temperature Distribution at t = ' num2str(n*dt/1e-9) 'ns']);
        xlabel('X (micro meter)');
        ylabel('Z (micro meter)');
        colormap(gca,"default");
        colorbar;
        clim([290 82000]);
        cb = colorbar(); 
        ylabel(cb,'Temperature (K)','Rotation',270)
        drawnow;
end

% Crater
evaporatedMask = T>Tb;
colorbar;
figure;
imagesc((1:Nx)*2, (1:Nx)*0.4, evaporatedMask)
xlabel('X (micro meter)');
ylabel('Z (micro meter)');
colormap(gca,"gray");
colorbar;
clim([0 1]);
cb = colorbar(); 
ylabel(cb,'Ablated material','Rotation',270)
drawnow;