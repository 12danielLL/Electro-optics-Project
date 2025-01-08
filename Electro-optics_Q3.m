%Q3
%%
%3-1
disp('3-1:')
disp('3-1B+C:');
R = 0.05; %radius of the circle [m]
L = 0.2; %viewing range [m]
N = 200; %number of samples in each dimension
delta_x = L / N;

%x,y axis
x = linspace(-L/2, L/2 - delta_x, N);
y = linspace(-L/2, L/2 - delta_x, N);

%circular key sample
circ_output = circ(L, N, R);

%plot
figure;
imagesc(x, y , circ_output);
axis xy;
axis image;
xlabel('x [m]');
ylabel('y [m]');
title('circular key try');

%%
%3-1D
disp('3-1D:');

%frequency axes
f_x = linspace(-1/(2*delta_x), 1/(2*delta_x) - 1/L, N); %[cycles/m]
f_y = linspace(-1/(2*delta_x), 1/(2*delta_x) - 1/L, N); %[cycles/m]

%Fourier transform
F_result = fftshift(fft2(circ_output));

axes=-L/2:L/N:L/2-L/N;
%plot
figure;
imagesc(axes,axes,abs(F_result));
axis xy;
axis image;
xlabel('x [m]');
ylabel('y [m]');
title('|fourier fransform|');
axis xy;
axis image;

%plot

[a, b] = meshgrid(-L/2:L/N:L/2-L/N, -L/2:L/N:L/2-L/N);
figure;
surf(a,b, abs(F_result));
axis square;
colormap(jet);
xlabel('x [m]');
ylabel('y[m]');
zlabel('|F|');
title('|fourier fransform| ');
camlight left; 
lighting phong; 
shading interp;

%%
%3-2
disp('3-2:');
%3-2C
disp('3-2C:');

lambda = 935.*10.^(-9); %[m]
z0 = 41999.90179; %[m]
k = 2*pi / lambda; %[rad/m]
i = sqrt(-1);
x = f_x * lambda * z0;
y = f_y * lambda * z0;

[X, Y] = meshgrid(x, y);
jinc_output = jinc(R * sqrt((X.^2 + Y.^2) / (lambda^2 * z0^2)));

%I~|E|^2
e1=exp(i * k * z0) ./ (lambda * z0 * i);
e2=exp(i *k * (X.^2 + Y.^2) / (2 * z0));
E =  e1*e2* R^2 .* jinc_output;

%plot intensity distribution
figure;
imagesc(x, y, abs(E).^2);%I~|E|^2
axis xy;
axis image;
colormap('jet');
colorbar;
xlabel('x[m]');
ylabel('y[m]');
title('Fraunhofer Pattern power distribution');

figure;
surf(x, y, abs(E).^2);
camlight left;
lighting phong;
shading interp
colorbar;
title('|E|^2');
xlabel('x[m]');
ylabel('y[m]');
zlabel('absolute value');

%%
%3-2D
disp('3-2D:');
mid_index = round(N / 2); %middle index
E_side_look =abs(E(:, mid_index));

%FWHM
half_peak = max(E_side_look) / 2; %value of half of the peak
above_half_peak = E_side_look > half_peak;%bigger valu from the half 
FWHM_range = find(diff(above_half_peak) ~= 0); 
FWHM = x(FWHM_range([1 end])); %FWHM size
fprintf('FWHM: %.4f to %.4f meters\n', FWHM(1), FWHM(2)); 

%plot
figure;
plot(x, E_side_look);
title('vertical cut');
xlabel('x[m]');
ylabel('abslute value');
hold on

plot(FWHM, [half_peak, half_peak], 'r');
text(FWHM(1), half_peak, sprintf('FWHM: %.4f m', FWHM(2) - FWHM(1)), 'VerticalAlignment', 'bottom');
hold off;

%%
%3-3
disp('3-3:');
%3-3B
disp('3-3B:');

% function [u2, x_prop] = propFresnel(u1, L, lambda, z)
%     % propagation – according to Fresnel
%     % assumes same x and y side lengths and
%     % uniform sampling
%     % u1 - source plane field
%     % L - source plane side length (FOV)
%     % lambda - wavelength
%     % z - propagation distance
%     % u2 - observation plane field
%     % x_prop – axis in observation plane field
% 
%     % Define the grid and sampling parameters
%     N = size(u1, 1); % Assume square grid
%     delta_x = L / N; % Sampling interval
%     x = (-N/2 : N/2-1) * delta_x; % Source plane coordinates
%     f_x = x / (lambda * z); % Frequency axis in the source plane
% 
%     % Fourier transform of the source plane field
%     U1 = F(u1);
% 
%     % Fresnel kernel
%     kernel = exp(1i * pi * lambda * z * (f_x.^2));
% 
%     % Propagation according to Fresnel
%     U2 = U1 .* kernel;
% 
%     % Inverse Fourier transform to get the observation plane field
%     u2 = iF(U2);
% 
%     % Define axis in the observation plane field
%     x_prop = x * (lambda * z);
% 
% end


%%
%3-3C
disp('3-3C:');

z1 = z0 ;%[m]
%circular key sample
u1 = circ(L, N, R);

%Fresnel pattern at z1
[u2, x_prop] = propFresnel(u1, L, lambda, z1);
%intensity distribution
intensity =abs(u2).^2;

figure;
imagesc(x_prop, x_prop, intensity);
axis xy;
colormap('jet');
colorbar;
xlabel('x[m]');
ylabel('y[m]');
title('intensity distribution');

figure;
surf(x_prop, x_prop, intensity);
camlight left;
lighting phong;
shading interp
colorbar;
title('|intensity|');
xlabel('x[m]');
ylabel('y[m]');
zlabel('absolute value');
camlight left;
lighting phong;
shading interp;

%FWHM of the intensity pattern
%x-direction
half_peak_X = max(intensity(:)) / 2;
X_range = find(intensity(size(intensity, 1) / 2, :) >= half_peak_X); %find indices
FWHM_X = x_prop(X_range(end)) - x_prop(X_range(1)); 

%y-direction
half_peak_Y = max(intensity(:)) / 2;
X_range = find(intensity(:, size(intensity, 2) / 2) >= half_peak_Y); %find indices 
FWHM_Y = x_prop(X_range(end)) - x_prop(X_range(1));

disp(['FWHM x-axis: ', num2str(FWHM_X), ' [m]']);
disp(['FWHM y-axis: ', num2str(FWHM_Y), ' [m]']);

%%
%3-3D
disp('3-3D:');
%same as 3-3C just change z1 from 'z0/50' to 'z0'
%%
%3-3E
disp('3-3E:');

distances = [10*z0, 2*z0, z0, 0.5*z0, 0.1*z0, 0.01*z0];
FWHM_X_val = zeros(size(distances));
FWHM_Y_val = zeros(size(distances));

figure;
for i = 1:length(distances)
   
    [u2, x_prop] = propFresnel(u1, L, lambda, distances(i));

    intensity = abs(u2).^2;
    
    %FWHM
    max_intensity = max(intensity(:));
    half_peak_intensity = max_intensity / 2;
    
    %indices where intensity is closest to half max
    [~, idx] = min(abs(intensity(:) - half_peak_intensity));
    [y_indx, x_indx] = ind2sub(size(intensity), idx);

    %FWHM in y-axis
    top_indx = max(1, y_indx - 1);
    down_indx = min(N, y_indx + 1);
    FWHM_Y = abs(x_prop(down_indx) - x_prop(top_indx));
    %FWHM in x-axis
    left_indx = max(1, x_indx - 1);
    right_indx = min(N, x_indx + 1);
    FWHM_X = abs(x_prop(right_indx) - x_prop(left_indx));
    
    %save FWHM values
    FWHM_X_val(i) = FWHM_X;
    FWHM_Y_val(i) = FWHM_Y;
    
    subplot(2, 6, i);
    imagesc(x_prop, x_prop, intensity);
    axis square;
    colormap('jet');
    colorbar;
    xlabel('x[m]');
    ylabel('y[m]');
    title(['Intensity at z = ', num2str(distances(i))]);
    
    subplot(2, 6, 6+i);
    plot(x_prop, intensity(round(N/2), :));
    xlabel('x[m]');
    ylabel('intensity');
    title(['intensity x-axis at z = ', num2str(distances(i))]);
end

%plot FWHM to distance
figure;
normalized_distances = distances / z0;
plot(normalized_distances, FWHM_X_val, '-o', 'DisplayName', 'FWHM x');
hold on;
plot(normalized_distances, FWHM_Y_val, '-o', 'DisplayName', 'FWHM y');
xlabel('Distance / z_0[m]');
ylabel('FWHM [m]');
title('FWHM to distance');
legend;





