clc;
clear;
close all;

%% Plot for concentration cA

% constants
R = 8.314;    
T = 5000;   
% D0 = 1e-10;    
% Q = 62000;
% D = D0 * exp(-Q / (R * T)); % Diffusion coefficient
D = 10;

% grid and time parameters
Nx = 150;      % grid size in x-direction
Ny = 150;      % grid size in y-direction
dx = 1e-9;    % grid spacing (1 nm)
dy = 1e-9;    % grid spacing (1 nm)
dt = 0.01;    % time step

time = 10000; % time

% material matrix
cA = zeros(Ny, Nx);
cA(1, :) = 0.8;
cB = zeros(Ny, Nx);


% Blue = cB, Yellow = cA
num_colors = 256;
%custom_colormap = [linspace(0,1,num_colors)', linspace(0,1,num_colors)', linspace(1,0,num_colors)']; 

for iter = 1:time
    cA_new = cA;
    
    for i = 2:Ny-1
        for j = 2:Nx-1
            cA_new(i, j) = cA(i, j) + D * dt * ( ...
                (cA(i+1, j) - 2*cA(i, j) + cA(i-1, j)) + ...
                (cA(i, j+1) - 2*cA(i, j) + cA(i, j-1)));
        end
    end
    
        for j = 2:Nx-1
            cA_new(Ny, j) = cA(Ny, j) + D * dt * ( ...
                            (cA(Ny-1, j) - cA(Ny, j)) + ...
                            (cA(Ny, j+1) - 2 * cA(Ny, j) + cA(Ny, j-1)));
            
            % periodic boundary conditions
            cA_new(:, 1) = cA_new(:, Nx-1);
            cA_new(:, Nx) = cA_new(:, 2);
            cA = cA_new;
        end
end

cA(10, 30) = 1;
cA(10, 31) = 1;
cA(11, 30) = 1;
cA(11, 31) = 1;

cA(10, 60) = 1;
cA(10, 61) = 1;
cA(11, 60) = 1;
cA(11, 61) = 1;

cA(10, 90) = 1;
cA(10, 91) = 1;
cA(11, 90) = 1;
cA(11, 91) = 1;

cA(10, 120) = 1;
cA(10, 121) = 1;
cA(11, 120) = 1;
cA(11, 121) = 1;

cA(20, 30) = 0.8;
cA(20, 31) = 0.8;
cA(21, 30) = 0.8;
cA(21, 31) = 0.8;

cA(20, 60) = 0.8;
cA(20, 61) = 0.8;
cA(21, 60) = 0.8;
cA(21, 61) = 0.8;

cA(20, 90) = 0.8;
cA(20, 91) = 0.8;
cA(21, 90) = 0.8;
cA(21, 91) = 0.8;

cA(20, 120) = 0.8;
cA(20, 121) = 0.8;
cA(21, 120) = 0.8;
cA(21, 121) = 0.8;

cA(30, 30) = 0.6;
cA(30, 31) = 0.6;
cA(31, 30) = 0.6;
cA(31, 31) = 0.6;

cA(30, 60) = 0.6;
cA(30, 61) = 0.6;
cA(31, 60) = 0.6;
cA(31, 61) = 0.6;

cA(30, 90) = 0.6;
cA(30, 91) = 0.6;
cA(31, 90) = 0.6;
cA(31, 91) = 0.6;

cA(30, 120) = 0.6;
cA(30, 121) = 0.6;
cA(31, 120) = 0.6;
cA(31, 121) = 0.6;

% for i = 1:Ny
%     if mod(i, 2) == 0
%         for j = 1:Nx
%             if mod(j, 2) ~= 0
%                 cA(i, j) = 1;
%             end
%         end
%     end
% end

    % visualisation
    imagesc(cA);
    colormap(parula);
    colorbar;
    caxis([0, 1]);
    title(sprintf('cA at Time %d', iter));
    pause(0.01);

    %% Concentration dependence on depth for different y, if h is depth

% select specific y positions for plotting
y_positions = [30, 60, 90, 120]; % y indices
Ny = size(cA, 1);
figure;
colors = lines(length(y_positions)); % distinct colors for curves

for u = 1:length(y_positions)
    y = y_positions(u);
    subplot(2, 2, u);
    plot(1:Ny, cA(:, y), 'Color', colors(u, :), 'LineWidth', 1.5);

    xlabel('Depth (y)');
    ylabel('Concentration');
    title(sprintf('Concentration dependence on depth at y = %d', y));
    grid on;
end

sgtitle('Concentration dependence on depth at different y positions');

%% Cahn-Hilliard

% i is vertical
dt = 0.01;
h = 1;
a = 13;
b = 2.2;
gamma = 2.9; % gradient energy coefficient

time = 600;
D0 = 0.001;
D = 20 * D0;

miu = zeros(Ny, Nx);

phi = zeros(Ny, Nx);
phi_new = zeros(Ny, Nx);

for i = 1:Ny
    for j = 1:Nx
        phi_new(i,j) = ((2 * cA(i,j)) - 1);
        phi0(i, j) = phi_new(i, j);
    end
end


for t = 1:time
    i = 1;
    j = 1;
    
    miu(i, j) = b.^4 * phi_new(i, j).^3 - a * b.^2 * phi_new(i, j) - gamma * ...
                (phi_new(i+1, j) - phi_new(i, j) - (phi_new(i, j+1) - phi_new(i, j)));
    
    phi(i, j) = D * (miu(i+1, j) - miu(i, j)) + ...
                D * (miu(i, j+1) - miu(i, j));
    
    phi_new(i, j) = phi0(i, j) + dt * phi(i, j);

    for j = 2:Nx-1
        
        miu(i, j) = b.^4 * phi_new(i, j).^3 - a * b.^2 * phi_new(i, j) - gamma * ...
                    (phi_new(i+1, j) - phi_new(i, j) + phi_new(i, j+1) - 2 * phi_new(i, j) + phi_new(i, j-1));
        
        phi(i, j) = D * (miu(i+1, j) - miu(i, j)) + ...
                    D * (miu(i, j+1) - miu(i, j)) - ...
                    D * (miu(i, j) - miu(i, j-1));
        
        phi_new(i, j) = phi0(i, j) + dt * phi(i, j);

        %phi(i,j) = phi_new(i,j);

    end
    
    % phi_new(i, Nx) = phi_new(i, 1);

    for i = 2:Ny-1
        j = 1;
        
        miu(i, j) = b.^4 * phi_new(i, j).^3 - a * b.^2 * phi_new(i, j) - gamma * ...
                    (phi_new(i+1, j) - 2 * phi_new(i, j) + phi_new(i-1, j) + (phi_new(i, j+1) -  phi_new(i, j)));
        
        phi(i, j) = D * (miu(i+1, j) - miu(i, j)) - ...
                    D * (miu(i, j) - miu(i-1, j)) - ...
                    D * (miu(i, j+1) - miu(i, j));
        
        phi_new(i, j) = phi0(i, j) + dt * phi(i, j);
        
        % phi_new(i, Nx) = phi_new(i, 1);

    for j = 2:Nx-1
        
        miu(i, j) = (b.^4 * phi_new(i, j).^3) - (a * b.^2 * phi_new(i, j)) - (gamma) * ...
                    (phi_new(i+1, j) - 2 * phi_new(i, j) + phi_new(i-1, j) + phi_new(i, j+1) - 2 * phi_new(i, j) + phi_new(i, j-1));
        
        phi(i, j) = D * (miu(i+1, j) - 2 * miu(i, j) + miu(i-1, j)) + D * (miu(i, j+1) - 2 * miu(i, j) + miu(i, j-1));
        
        phi_new(i, j) = phi0(i, j) + dt * phi(i, j);
        
        %phi(i,j) = phi_new(i,j);

    end
    
    % phi_new(i, Nx) = phi_new(i, 1);

    end
    
    %phi_new(:, 1) = phi_new(:, Nx-1);
    %phi_new(:, Nx) = phi_new(:, 2);

    for i = 1:Ny
        phi_new(i, 1) = phi_new(i, Nx-1);
        phi_new(i, Nx) = phi_new(i, 2);
    end

    for j = 1:Nx
       
        phi_new(Ny, j) = phi_new(Ny-1, j);

    end
        
    phi = phi_new;
    phi0 = phi_new;

end

% visualisation
num_colors = 256;
hImg = imagesc(phi_new);  
colormap(parula); 
colorbar;
caxis([-1, 1]);
title('Diffusion');    
set(hImg, 'CData', phi_new);
title(sprintf('Diffusion at Time %d', t));
pause(0.01);

%% Phi dependence on depth for different y, if h is depth

% select specific y positions for plotting
y_positions = [30, 115]; % y indices
figure;
hold on;
colors = lines(length(y_positions)); % distinct colors for curves

for u = 1:length(y_positions)
    y = y_positions(u);
    plot(1:Ny, phi_new(:, y), 'Color', colors(u, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('y = %d', y));
end

xlabel('Depth (y)');
ylabel('Phi');
title('Phi Dependence on Depth');
legend('show');
grid on;
hold off;