function [u2, x_prop] = propFresnel(u1, L, lambda, z)
    % propagation – according to Fresnel
    % assumes same x and y side lengths and
    % uniform sampling
    % u1 - source plane field
    % L - source plane side length (FOV)
    % lambda - wavelength
    % z - propagation distance
    % u2 - observation plane field
    % x_prop – axis in observation plane field

    
    i=sqrt(-1);
    N = size(u1, 1); %square grid
    delta_x = L / N; 
    x = (-N/2 : N/2-1) * delta_x; % Source plane coordinates
    f_x = x / (lambda * z);

    U1 = F(u1);
    kernel = exp(i * pi * lambda * z * (f_x.^2));
    U2 = U1 .* kernel;
    u2 = iF(U2);

    x_prop = x * lambda * z;

end
