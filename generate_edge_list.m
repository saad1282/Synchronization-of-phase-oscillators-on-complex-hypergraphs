%% HYPERGRAPH EDGELIST. This file generates the list of edges and triangles

function [edgelist, k] = generate_edge_list(hyperedge_size, N)
% hyperedge_size = size of a hyperedge; 2 for links and 3 for triangles
% N = total number of oscillators

%% degree distribution

% power law network

% k_min = 30; %minimum degree of a node
% k_max = 70; % maximum degree of a node
% u = rand(N,1); % random number on a cdf
% gamma = 3;    % gamma parameter for power law distribution
% a = k_min^(1-gamma);
% b = k_max^(1-gamma);
% k = ceil((a + u*(b-a)).^(1/(1-gamma)));

% k-regular network

% k_mean = 30;
% k = k_mean*(ones(N,1));

% uniform degree distribution network

min = 10;   % minimum degree
max = 30;   % maximum degree
prob = 1/(max - min);    % probability distribution
k = ceil((1/prob)*rand(N,1)+min);   % degree sequence

% Bimodal degree distribution
% 
% k_first = 20;
% k_second = 50;
% 
% k = rand(N,1);
% 
% for j = 1: length(k)
%     
%     if k(j) <= 0.5
%         k(j) = k_first;
%     else k(j) = k_second;
%     end
%     
% end

%% Cleaning up the degree sequence

k_sum = sum(k); % degree sequence sum
remainder = mod(k_sum, hyperedge_size); % as the sum of degree should be divisible by hyperedge_size

if (remainder ~= 0)
    update = hyperedge_size - remainder;
    ind = length(k);
    i = randsample(ind,1);
    k(i) = k(i) + update;
end

%% Create stublist

index = 1;
stub_num = sum(k);
stublist = zeros(stub_num, 1);

for j = 1 : length(k)
    for l = 1 : k(j)
        stublist(index) = j;
        index = index + 1;
    end
end

%% Use stublist to get hyperedges/edges

edgelist = zeros(stub_num/hyperedge_size, hyperedge_size);
index1 = 1;

while (stub_num > 0)
    
    indices = randsample(stub_num,hyperedge_size);
    edgelist(index1,:) = stublist(indices);
    stublist(indices) = [];
    stub_num = length(stublist);
    index1 = index1+1;
    
end
