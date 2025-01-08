%2
disp('2:');
%2C1
disp('2C1:');
lambda_margin = (0.4:0.01:1);%𝜇m

f_len_new = zeros(size(lambda_margin));
n_val_of_n7 = zeros(size(lambda_margin));

%focal lengths for every wavelength
for j = 1:length(lambda_margin)
    wavelength = lambda_margin(j);

    %n_BK7 formula
    n = sqrt(1 + 1.03961212./(1-0.00600069867./wavelength^2) + 0.231792344./(1-0.0200179144./wavelength^2) + 1.01046945./(1-103.560653./wavelength^2));
    n_val_of_n7(j) = n;

    %lens maker's formula
    f = (n - 1) * ((1/(-44)) - (1/(-5))); % R1 = -44 and R2 = -5
    f_len_new(j) = 1 / f;
end


%plot n to λ
figure;
subplot(2, 1, 1);
plot(lambda_margin, n_val_of_n7, 'g', 'LineWidth', 2);
xlabel('λ(\mum)');
ylabel('n');
title('refractive index of lens to wavelength');
grid on;

%plot f to λ
subplot(2, 1, 2);
plot(lambda_margin, f_len_new, 'r', 'LineWidth', 2);
xlabel('λ[\mum]');
ylabel('f [mm]');
title('focal length to wavelength');
grid on;
%%
%2C2
disp('2C2:');
R = 34e-3; %radii of curvature
lambda1 = 744e-9;%[m]
lambda2 = 534e-9;%[m]
f = 33.203125e-3; %[m]
L1 = 70e-3;%d to lens [m]
L2 = 63.16348195e-3;%d from lens to figure[m]
theta_range = linspace(pi/5, -pi/5, 11);%angles 

%progression in space matrix
prog_mat = [1, L1; 0, 1];
%thin lens matrix
thin_lens_mat = [1, 0; -1/f, 1];

heights_arry = zeros(length(theta_range), 3);

%run over angles
for j = 1:length(theta_range)
    %incoming beam hight- 0.001[m]
    beam_initial = [0.001; tan(theta_range(j))];
    
    to_hit = prog_mat * beam_initial;
    hit = thin_lens_mat * to_hit;
    beam_finial = prog_mat * hit;
    
    heights_arry(j, :) = [beam_initial(1), hit(1), beam_finial(1)];
end

%plot
figure;
hold on;
for j = 1:length(theta_range)
    plot([0, L1, L1 + L2], heights_arry(j, :), '-');
end
xlabel('distance[m]');
ylabel('height[m]');
title('progression of the rays in the system');
grid on;

%%
%2C3
disp('2C3:');

f = 33.203125e-3; %[m]

%ray heights
heights_lambda1 = zeros(length(-10:1:10), 3);
heights_lambda2 = zeros(length(-10:1:10), 3);

%iterate over heights
h_range = -0.001:0.0001:0.001;
for j = 1:length(h_range)
    %incoming beam
    beam_initial = [h_range(j); 0];
    
    to_hit = prog_mat * beam_initial;
    hit = thin_lens_mat * to_hit;
    beam_finial = prog_mat * hit;
    
    heights_lambda1(j, :) = [beam_initial(1), hit(1), beam_finial(1)];
    
    %lambda2
    thin_lens_mat_lambda2 = [1, 0; -1/(f*(lambda2/lambda1)), 1];
    hit_lambda2 = thin_lens_mat_lambda2 * to_hit;
    beam_finial_lambda2 = prog_mat * hit_lambda2;
    
    heights_lambda2(j, :) = [beam_initial(1), hit_lambda2(1), beam_finial_lambda2(1)];
end

%plot
figure;
hold on;
%lambda1
plot([0, L1, L1 + L2], heights_lambda1(:, :), 'red');
%lambda2
plot([0, L1, L1 + L2], heights_lambda2(:, :), 'green');
xlabel('distance[m]');
ylabel('height[m]');
title('progression of the rays in the system');
hold off;

%%
%2D1
%%
%2D2
disp('2D2:');
R2=-34.*10.^(-3);%[m]
R1=34.*10.^(-3);%[m]
R3=136.503.*10.^(-3);%[m]
d1=0.22.*10.^(-3);%[m]
d2=2.59.*10.^(-3);%[m]
lambda_margin = (400*(10.^(-9)):10*(10.^-9):1000*(10.^(-9)));%[m]

f_len = zeros(size(lambda_margin));
n_val_of_n7 = zeros(size(lambda_margin));
n_val_of_n2 = zeros(size(lambda_margin));

for j = 1:length(lambda_margin)
    wavelength = lambda_margin(j);
%n_BK7 formula
n7=sqrt(1 + 1.03961212./(1-0.00600069867./(wavelength.*10^6)^2) + 0.231792344./(1-0.0200179144./(wavelength.*10^6)^2) + 1.01046945./(1-103.560653./(wavelength.*10^6)^2));
n_val_of_n7(j) = n7;
%n_F2 formula
n2=sqrt(1 + 1.34533359./(1-0.00997743871./(wavelength.*10^6)^2) + 0.209073176./(1-0.0470450767./(wavelength.*10^6)^2) + 0.937357162./(1-111.886764./(wavelength.*10^6)^2));
n_val_of_n2(j) = n2;

numerator = n2 .* R3 .* (n2 .* R2 .* (-n7 + 1) + (n7 - n2) .* (n7 .* R1 + d1 .* (-n7 + 1))) + ...
            (d2 .* (n2 * R2 .* (-n7 + 1) + (d1 .* (1 - n7) + n7 .* R1) .* (n7 - n2)) + n2 .* R2 .* (n7 .* R1 + d1 .* (-n7 + 1))) .* (n7 - n2);
            
denominator = n7 .* n2^2 .* R1 .* R2 * R3;

expression = numerator ./ denominator;
f_len(j)=-1./expression;
end

%refractive index to wavelength
figure;
subplot(2, 1, 1);
plot(lambda_margin, n_val_of_n7, 'r', 'LineWidth', 2);
xlabel('λ[m]');
ylabel('n_{BK7}');
title('refractive index - Wavelength (n_{BK7})');
grid on;

subplot(2, 1, 2);
plot(lambda_margin, n_val_of_n2, 'r', 'LineWidth', 2);
xlabel('λ[m]');
ylabel('n_{F2}');
title('refractive index - Wavelength (n_{F2})');
grid on;

%avarages
f_avg=mean(f_len);
disp(f_avg);
n7_avg=mean(n_val_of_n7);
disp(n7_avg);
R1_avg=2.*f_avg.*(n7_avg-1);
R2_avg=-2.*f_avg.*(n7_avg-1);

lambda_margin = (400*(10.^(-9)):10*(10.^-9):1000*(10.^(-9)));%[m]
f_len_new = zeros(size(lambda_margin));

for i = 1:length(lambda_margin)
    wavelength = lambda_margin(i);

    n7=sqrt(1 + 1.03961212./(1-0.00600069867./(wavelength.*10^6)^2) + 0.231792344./(1-0.0200179144./(wavelength.*10^6)^2) + 1.01046945./(1-103.560653./(wavelength.*10^6)^2));
    n_val_of_n7(i) = n7;
    C_ABCD = (n7 - 1) .* ((1./(R1_avg)) - (1./(R2_avg))); 
    f_len_new(i) = 1 ./ C_ABCD;
end
%plot focal length to wavelength
figure;
plot(lambda_margin, f_len, 'g', 'LineWidth', 2);
hold on;
plot(lambda_margin, f_len_new, 'r', 'LineWidth', 2);
hold off;
xlabel('λ[m]');
ylabel('f[m]');
title('focal length - wavelngth');
legend('with layer', 'without layer');

%%
%2D3
disp('2D3:');
%Principal planes
h1 = zeros(size(lambda_margin));
h2 = zeros(size(lambda_margin));

n_val_of_n7_D3 = zeros(size(lambda_margin));
n_val_of_n2_D3 = zeros(size(lambda_margin));

%ABCD matrix
A_bcd = zeros(size(lambda_margin));
a_B_cd = zeros(size(lambda_margin));
ab_C_d = zeros(size(lambda_margin));
abc_D = zeros(size(lambda_margin));

lambda_margin = (400*(10.^(-9)):10*(10.^-9):1000*(10.^(-9)));%[m]
for i = 1:length(lambda_margin)
    wavelength = lambda_margin(i);

   %n_BK7 formula
   n7=sqrt(1 + 1.03961212./(1-0.00600069867./(wavelength.*10^6)^2) + 0.231792344./(1-0.0200179144./(wavelength.*10^6)^2) + 1.01046945./(1-103.560653./(wavelength.*10^6)^2));
   n_val_of_n7_D3(i)=n7;
   %n_F2 formula
   n2=sqrt(1 + 1.34533359./(1-0.00997743871./(wavelength.*10^6)^2) + 0.209073176./(1-0.0470450767./(wavelength.*10^6)^2) + 0.937357162./(1-111.886764./(wavelength.*10^6)^2));
   n_val_of_n2_D3(i)=n2;

    %spherical surface hit matrix
    R1_mat = [1, 0; (1-n7)/(n7*R1_avg), 1/n7];
    R2_mat = [1, 0; (n7-n2)/(n2*R2_avg), n7/n2];
    R3_mat = [1, 0; (n2-1)/(1*R3), n2/1];
    %space progression matrix
    d1_mat = [1, d1; 0, 1];
    d2_mat = [1, d2; 0, 1];
    %total system matrix
    tot_systen_mat = R3_mat * d2_mat * R2_mat * d1_mat * R1_mat;
    
    A_bcd(i) = tot_systen_mat(1, 1);
    a_B_cd(i) = tot_systen_mat(1, 2);
    ab_C_d(i) = tot_systen_mat(2, 1);
    abc_D(i) = tot_systen_mat(2, 2);

    %h1 and h2 for each wavelength
    h1(i) = (abc_D(i) - 1 )/ab_C_d(i);
    h2(i) = (1 - A_bcd(i))/ab_C_d(i);
end

h1_avg = mean(h1);
h2_avg = mean(h2);
disp(['h1_average = ', num2str(h1_avg)]);
disp(['h2_average = ', num2str(h2_avg)]);


%plot h1 to wavelength 
figure;
subplot(2,1,1)
plot(lambda_margin, h1, 'b', 'DisplayName', 'h1');
xlabel('λ[m]');
ylabel('h1[m]');
title('h1 - λ');
legend('h1', 'Location', 'Best');
grid on;

%plot h2 to wavelength 
subplot(2,1,2)
plot(lambda_margin, h2, 'g', 'DisplayName', 'h2');
xlabel('λ[m]');
ylabel('h2[m]');
title('h2 - λ');
legend('h2', 'Location', 'Best');
grid on;














