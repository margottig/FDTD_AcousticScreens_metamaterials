% Simulation parameters
Nx = 400;    % Number of grid points in x direction
Ny = 200;    % Number of grid points in y direction
dx = 1e-4;   % Grid spacing in x direction
dy = 1e-4;   % Grid spacing in y direction
c_air = 343; % Speed of sound in air (m/s)
freq = 1000; % Frequency (Hz)
lamda = c_air / freq; % Wavelength (m)
dt = 0.5 * min(dx, dy) / c_air; % Time step
Nt = 600;    % Number of time steps

% Lens parameters
lens_x = Nx / 2;  % Center of the lens in x direction
lens_y = Ny / 4;  % Position of the lens in y direction
h = 0.1625; % Height of the medium in meters
d = 0.25 * lamda; % Width of the medium in meters
n_o = 3;   % Index of refraction profile high
n_h = 1;   % Index of refraction profile low
a = (1 / h) * acosh(n_o / n_h); % Parameter for the hyperbolic secant profile

% Source parameters
source_x = Nx / 2; % Source position in x direction
source_y = 10;   % Source position in y direction
source_freq = 100e3; % Source frequency (Hz)

% Geometry of white blocks
blocks = [ 80 90 20 10;
          120 90 30 10;
          160 90 15 10;
          200 90 10 10;
          240 90 20 10;
          280 90 15 10;
          320 90 20 10];
      
% Initialize pressure fields
p = zeros(Nx, Ny);   % Pressure at current time step
p_old = zeros(Nx, Ny); % Pressure at previous time step
p_new = zeros(Nx, Ny); % Pressure at next time step

% Create lens (index of refraction profile)
lens = ones(Nx, Ny) * c_air;
for ix = 1:Nx
    for iy = 1:Ny
        x = (ix - lens_x) * dx;
        y = (iy - lens_y) * dy;
        if abs(x) <= d
            lens(ix, iy) = c_air / (n_h * cosh(a * x));
        end
    end
end

% FDTD loop
for n = 1:Nt
    % Update pressure field
    for ix = 2:Nx-1
        for iy = 2:Ny-1
            c = lens(ix, iy);
            p_new(ix, iy) = 2 * p(ix, iy) - p_old(ix, iy) + (c^2 * dt^2) * ...
                ((p(ix+1, iy) - 2 * p(ix, iy) + p(ix-1, iy)) / dx^2 + ...
                 (p(ix, iy+1) - 2 * p(ix, iy) + p(ix, iy-1)) / dy^2);
        end
    end
    
    % Apply source
    p_new(source_x, source_y) = p_new(source_x, source_y) + ...
        sin(2 * pi * source_freq * n * dt);
    
    % Apply white blocks (set pressure inside blocks to zero)
    for b = 1:size(blocks, 1)
        bx = blocks(b, 1);
        by = blocks(b, 2);
        bw = blocks(b, 3);
        bh = blocks(b, 4);
        p_new(bx:bx+bw, by:by+bh) = 0;
    end
    
    % Update fields
    p_old = p;
    p = p_new;
    
    % Visualization
    if mod(n, 10) == 0
        imagesc(p', [-1 1]);
        colorbar;
        hold on;
        % Plot the white blocks
        for b = 1:size(blocks, 1)
            bx = blocks(b, 1);
            by = blocks(b, 2);
            bw = blocks(b, 3);
            bh = blocks(b, 4);
            rectangle('Position', [bx, by, bw, bh], 'FaceColor', 'w');
        end
        title(['Time step: ', num2str(n)]);
        hold off;
        pause(0.01);
    end
end
