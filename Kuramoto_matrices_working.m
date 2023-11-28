% Kuramoto Model for Hypergraphs (edges and triangles) 

tic
N = 5000; %% Number of Oscillators

A = importdata('edgelist_5000.txt'); % list of edges;
B = importdata('trianglelist_5000.txt'); % list of triangles;
k1 = importdata('k1_degree_5000.txt'); % degree of each node;
k2 = importdata('k2_degree_5000.txt'); % triangle degree of each node;
w = importdata('initial_frequency_5000.txt'); % frequency from Cauchy distribution;

K1 = 0:0.001:0.07; % Pairwise coupling strength
K2 = 0.05; % Triangular coupling strength

k_c = 2*mean(k1)/mean(k1.^2); % critical K_1 value that marks onset of synchronization
k_bicritical = (mean(k1.^4)*mean(k1)^2)/(mean(k1.^2)^2*mean(k1.^3)); % critical K_2 value for bistability in correlated links and triangles

% Initial distribution of oscillators
% a = 0; 
% b = 2*pi;
% d1_theta = (b-a)/N;
% theta_given = transpose(a + (0:1:N-1)*d1_theta);        
theta_given = pi*ones(N,1);

T = 10000; %% number of time steps for Euler Method
delta_t = 0.001;

%% Initialize Order Parameters

R1 = zeros(N,1);     %% order parameter of oscillators
R2 = zeros(N,1);     %% order parameter of oscillators

R1_out = zeros(length(K1),1);       %% global order parameter
R2_out = zeros(length(K1),1);       %% global order parameter 
z_out = zeros(length(K1),1);        %% e^{i\theta}/N

norm_R1 = zeros(T,1);       %% norm of R1 for each t
norm_R2 = zeros(T,1);       %% norm of R1 for each t
z = zeros(T,1);

norm_z = zeros(T,1);

theta_sim = zeros(N,1);     %% Simulated theta values

%% Simulation in increasing K1

for i = 1:length(K1)

    for t = 1 : T
        
        R1(:,1) = 0;
        R2(:,1) = 0;
        
        for l = 1 : length(A) % local order parameter R1
            R1(A(l,1)) = R1(A(l,1)) + exp(1j*theta_given(A(l,2)));
            R1(A(l,2)) = R1(A(l,2)) + exp(1j*theta_given(A(l,1)));
        end
        
        for m = 1 : length(B) % local order parameter R2
            R2(B(m,1)) = R2(B(m,1)) + exp(2j*theta_given(B(m,2)) - 1j*theta_given(B(m,3))) + exp(2j*theta_given(B(m,3)) - 1j*theta_given(B(m,2)));
            R2(B(m,2)) = R2(B(m,2)) + exp(2j*theta_given(B(m,1)) - 1j*theta_given(B(m,3))) + exp(2j*theta_given(B(m,3)) - 1j*theta_given(B(m,1)));
            R2(B(m,3)) = R2(B(m,3)) + exp(2j*theta_given(B(m,2)) - 1j*theta_given(B(m,1))) + exp(2j*theta_given(B(m,1)) - 1j*theta_given(B(m,2)));
        end
        
        d_theta = w + K1(i) * imag(R1.*exp(-1j*theta_given)) + K2 * imag(R2.*exp(-1j*theta_given));
        theta_sim = theta_given + d_theta.*delta_t;
        
        theta_given = theta_sim;
        norm_R1(t) = norm((sum(R1(:,1))/(2*length(A))));
        norm_R2(t) = norm((sum(R2(:,1))/(6*length(B))));
        z(t) = (sum(exp(1j*theta_given)))/N;
        norm_z(t) = norm(z(t));        
    end
    
    R1_out(i) = mean(norm_R1(T-1000:T));    % Global order parameter
    R2_out(i) = mean(norm_R2(T-1000:T));    % Global order parameter
    z_out(i) = mean(norm_z(T-1000:T));
    i
end

%% Initialize Order Parameters for decreasing K1

R1back = zeros(N,1);     %% order parameter for oscillators
R2back = zeros(N,1);     %% order parameter for oscillators

R1_out_back = zeros(length(K1),1);       %% global order parameter for K1
R2_out_back = zeros(length(K1),1);       %% global order parameter for K1
z_out_back = zeros(length(K1),1);        %% e^{i\theta}/N

norm_R1_back = zeros(T,1);       %% normalized R1 for each t
norm_R2_back = zeros(T,1);       %% normalized R1 for each t
z_back = zeros(T,1);

norm_z_back = zeros(T,1);

theta_sim_back = zeros(N,1);     %% Simulated theta values

% Simulation for decreasing K1

for j = length(K1):-1:1

    for t = 1 : T

        R1back(:,1) = 0;
        R2back(:,1) = 0;

        for l = 1 : length(A)
            R1back(A(l,1)) = R1back(A(l,1)) + exp(1j*theta_given(A(l,2)));
            R1back(A(l,2)) = R1back(A(l,2)) + exp(1j*theta_given(A(l,1)));
        end

        for m = 1 : length(B)
            R2back(B(m,1)) = R2back(B(m,1)) + exp(2j*theta_given(B(m,2)) - 1j*theta_given(B(m,3))) + exp(2j*theta_given(B(m,3)) - 1j*theta_given(B(m,2)));
            R2back(B(m,2)) = R2back(B(m,2)) + exp(2j*theta_given(B(m,1)) - 1j*theta_given(B(m,3))) + exp(2j*theta_given(B(m,3)) - 1j*theta_given(B(m,1)));
            R2back(B(m,3)) = R2back(B(m,3)) + exp(2j*theta_given(B(m,2)) - 1j*theta_given(B(m,1))) + exp(2j*theta_given(B(m,1)) - 1j*theta_given(B(m,2)));
        end

        d_theta = w + K1(j) * imag(R1back.*exp(-1j*theta_given)) + K2 * imag(R2back.*exp(-1j*theta_given));
        theta_sim_back = theta_given + d_theta.*delta_t;

        theta_given = theta_sim_back;
        norm_R1_back(t) = norm((sum(R1back(:,1))/(2*length(A))));
        norm_R2_back(t) = norm((sum(R2back(:,1))/(6*length(B))));
        z_back(t) = (sum(exp(1j*theta_given)))/N;
        norm_z_back(t) = norm(z_back(t));

    end

    R1_out_back(j) = mean(norm_R1_back(T-1000:T));
    R2_out_back(j) = mean(norm_R2_back(T-1000:T));
    z_out_back(j) = mean(norm_z_back(T-1000:T));
    j
end

plot(K1,R1_out, 'o-r');
hold on;
plot(K1,R2_out, 'o-b');
plot(K1, z_out,  'o-g');
plot(K1,R1_out_back,  'x-r');
hold on;
plot(K1,R2_out_back, 'x-b');
plot(K1, z_out_back,  'x-g');
legend('R1 for', 'R2 for', 'R for', 'R1 back', 'R2 back', 'R back');
xlabel('K1');
%ylim([0 0.5])
ylabel('order parameters');
title('order parameter R1 and R2 for varying K1 and K2 = 0.05');
toc