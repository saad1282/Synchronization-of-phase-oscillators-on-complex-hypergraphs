%% Create list of edges and triangles

N = 5000;
[A, k1] = generate_edge_list(2,N);  % list of links and its degree
[B, k2] = generate_edge_list(3,N);  % list of triangles and its degree

% Cauchy distribution for natural frequency
% 
% u = rand(N,1);
% x_0 = 0; % mean
% beta = 1; % scale parameter
% w = x_0 + beta.*(tan(pi*(u-1/2)));
% This has infinite range, so we will truncate it

w = zeros(N,1);
for i = 1 : N
    %w(i) = tan((i*pi)/N - (N+1)*pi/(2*N)); % This gives natural frequency between about -3e+3 and 3e+3
    w(i) = tan(pi*(2*i - N -1)*(N+1)); % This gives natural frequency close to 0
end

% adjacency matrix of the network, although we do not use it.
% adjacency_matrix = zeros(N,N);
 
% for i = 1:length(A)
%     adjacency_matrix(A(i,1),A(i,2)) = 1;
% end
 
%% save edgelists and natural frequencies as txt file

% save('edgelist.txt', 'A', '-ASCII','-append');
% type('edgelist.txt');
% save('k1_degree.txt', 'k1', '-ASCII','-append');
% type('k1_degree.txt');
% save('trianglelist.txt', 'B', '-ASCII','-append');
% type('trianglelist.txt');
% save('k2_degree.txt', 'k2', '-ASCII','-append');
% type('k2_degree.txt');
% save('initial_frequency.txt', 'w', '-ASCII','-append');
% type('initial_frequency.txt');
% save('adj_matrix.txt', 'adjacency_matrix', '-ASCII','-append');
% type('adj_matrix.txt');