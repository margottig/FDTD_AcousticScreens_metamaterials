function [h,h2,h3]=fdtd2d_metalens(metalens, centralfrequency)

%%  --------------- GENERAL -----------------------------
% Initial parameters
dh = .01;       % Spatial definition [m]
dt = dh/341/2;  % Temporal definition [s]
ts = 15;         % Simulation time [ms]
% Medium characteristics
rho = 1.21;     % Medium density [kg/m3]
c = 341;        % Propagation speed [m/s]
K = (c^2)*rho;  % Compressibility modulus

%% -------------------- MAPPING MATRICES ----------------------------
% Map domain dimensions in meters and PML
lx = 5;         % Width [m]
ly = 4;         % Height [m]
lengthPML = 0.25;    % PML size/width/layer in m
% Map sizes matrix 
nPML = round(lengthPML/dh); 
nx   = round(lx/dh)+nPML*2; % Times 2 because PML is above and below
ny   = round(ly/dh)+nPML*2; % Times 2 because PML is to the left and right
% Create map matrix
p  = zeros(nx,ny);   % Pressure matrices
px = zeros(nx,ny);
py = zeros(nx,ny);
ux = zeros(nx+1,ny); % Particle velocity matrix in x
uy = zeros(nx,ny+1); % Particle velocity matrix in y

%% ------------- Receiver and source positions ------------
sound_source_position = [0.9,0.25]; %[m,m]
posNy = round((sound_source_position(2))/dh) + nPML; % source position in y
posNx = round(sound_source_position(1)/dh) + nPML; % source position in x
mic1 = [1.25,0.5]; % Receiver 1 position [m,m]
mic2 = [2,0.9];% Receiver 2 position [m,m]
mic3 = [1.5,3]; % Receiver 3 position [m,m]

%% ----------------------- PML ----------------------------
gammamax = 0.5; % Maximum impedance reduction / percentage
% Left and right PML
gammaux = zeros(nx+1,ny);
gammaux(1:nPML,:) = repmat(gammamax*((nPML:-1:1)'/nPML).^2,1,ny); % Left
gammaux(1+end-nPML:end,:) = repmat(gammamax*((1:1:nPML)'/nPML).^2,1,ny); % Right
gammax = (gammaux(1:end-1,:)+gammaux(2:end,:))/2; % Set
% Upper and lower PML
gammauy = zeros(nx,ny+1);
gammauy(:,1:nPML) = repmat(gammamax*((nPML:-1:1)/nPML).^2,nx,1); % Lower
gammauy(:,1+end-nPML:end) = repmat(gammamax*((1:1:nPML)/nPML).^2,nx,1); % Upper
gammay = (gammauy(:,1:end-1)+gammauy(:,2:end))/2; % Set

%% EXCITATION
maxtt = ts*10^-3/dt; % Length of the time vector
a    = centralfrequency/(sqrt(pi)/2)*4; 
t    = ((1:maxtt)/(1/dt)-4/a); % Time vector
w    = -(exp(-a^2*(t.^2)/2).*(a^2*(t.^2)-1)); % Ricker

%% Initialization of the listening point vectors
h = zeros(1,maxtt);
h2 = zeros(1,maxtt);
h3 = zeros(1,maxtt);

%% Boundary conditions
ux(1,:)   = -p(1,:)/rho/c;
ux(end,:) = p(end,:)/rho/c;
uy(:,1)   = -p(:,1)/rho/c;
uy(:,end) = p(:,end)/rho/c;

%% Read metamaterial from BMP File
metamaterialFile = metalens; % Path to the BMP file
metamaterialImage = imread(metamaterialFile); % Read the BMP file
metamaterialImage = imbinarize(metamaterialImage); % Convert to binary image (0 and 1)
metamaterialImage = imresize(metamaterialImage, 0.4,'nearest'); % Resize to fit the simulation grid
metamaterialIndices = find(metamaterialImage == 1);
[metamaterialX, metamaterialY] = ind2sub(size(metamaterialImage), metamaterialIndices);
metamaterialY = metamaterialY - 100; % Move Y indices

%% Calculation
fig1 = figure('Color',[1,1,1]);
colormap(jet(256))

for tt=1:maxtt
    % Pressure
    px = px.*(1-gammax)-K*dt/dh*diff(ux,1,1);
    py = py.*(1-gammay)-K*dt/dh*diff(uy,1,2);
    
    % Excitation
    px(posNx,posNy) = w(tt);
    py(posNx,posNy) = w(tt);
    p = px+py;
    
    % Velocity
    ux(2:nx,:) = ux(2:nx,:).*(1-gammaux(2:nx,:))-dt/rho/dh*diff(p,1,1);
    uy(:,2:ny) = uy(:,2:ny).*(1-gammauy(:,2:ny))-dt/rho/dh*diff(p,1,2);
    
    % setting the particle velocities to zero at the positions corresponding to the metamaterial
    for i = 1:length(metamaterialX)
        if metamaterialX(i) <= nx && metamaterialY(i) <= ny %     % Check if Indices are within Grid Bounds
            ux(metamaterialX(i), metamaterialY(i)) = 0;
            uy(metamaterialX(i), metamaterialY(i)) = 0;
        end
    end
    
    % Impulse response
    h(tt) = p(round(mic1(1)/dh)+nPML,round(mic1(2)/dh)+nPML);
    h2(tt) = p(round(mic2(1)/dh)+nPML,round(mic2(2)/dh)+nPML);
    h3(tt) = p(round(mic3(1)/dh)+nPML,round(mic3(2)/dh)+nPML);
    
    % Graphical representation
    if tt/100==round(tt/100)*100
        splmap = 10*log10(abs(p).^2/(2e-5)^2);
        %splmap = splmap - max(splmap(:));
        pcolor(((1:nx)-nPML)*dh,((1:ny)-nPML)*dh,splmap');
        hold on;
        % Show position of the microphones
        plot([mic1(1),mic2(1),mic3(1)],[mic1(2),mic2(2),mic3(2)], 'ro', 'MarkerSize', 2,'MarkerFaceColor','r');
        text([mic1(1),mic2(1),mic3(1)]+0.05,[mic1(2),mic2(2),mic3(2)],{'mic1','mic2','mic3'},'Color','white')

        % Plot metamaterial indices in white
        plot((metamaterialX-nPML)*dh, (metamaterialY-nPML)*dh, 'g.', 'MarkerSize', 1);

        hold off
        axis equal;axis([-nPML,nx-nPML,-nPML,ny-nPML]*dh)
        set(gca,'Clim',[50 110]);shading flat,cc = colorbar;
        title(['Time = ' num2str(round((tt)*1000*dt)) ' ms']);
        xlabel('X [m]'),ylabel('Y [m]'),cc.Label.String = 'SPL (dB)';
        drawnow
    end
end

%% ---- PLOTS --------------
time = seconds(t);
plot_received_signals(time, h, h2, h3);

return


% bmp files taken from Alagoz, Baris Baykant, and Serkan Alagoz. “Towards Earthquake Shields: A Numerical Investigation of Earthquake Shielding with Seismic Crystals.” Open Journal of Acoustics, vol. 01, no. 03, Scientific Research Publishing, Inc., 2011, pp. 63–69, doi:10.4236/oja.2011.13008.
