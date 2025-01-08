%1
disp('1:');
%1D
disp('1D:');

n0 = 1.694; 
a = [1e26, 1e27, 1e28]; %Hz^2
d = 410216.4882; %[nm]
theta_i = 36; %dagrees

lambda_margin = linspace(0, 800, 1000); %[nm]

T_value = zeros(length(a), length(lambda_margin));

%loop over each coefficient 'a'
for i = 1:length(a)
    %refractive index for each wavelength
    v = physconst('LightSpeed') ./ (lambda_margin * 1e-9); %Hz
    n = n0 - a(i) ./ v.^2;
    
    theta_t = asin(sin(theta_i) ./ n);
    
    r_s = sin(theta_t - theta_i) ./ sin(theta_t + theta_i);
    R = abs(r_s).^2;
    
    delta = 4 * pi * (n0 - a(i) ./ v.^2) .* d .* cos(theta_t) ./ lambda;
    
    %transfer coefficient
    T = ((1 - R).^2) ./ ((1 - R).^2 + 4 * R .* sin(delta/2).^2);
    
    T_value(i, :) = T;
end

%plot transfer coefficients for each 'a'
figure;
hold on;
for i = 1:length(a)
    plot(lambda_margin, T_value(i, :), 'LineWidth', 1.5);
end
hold off;

%plot
title('transfer coefficient to wavelength');
xlabel('wavelength[nm]');
ylabel('transfer coefficient');
legend('a = 1e26', 'a = 1e27', 'a = 1e28');
grid on;





