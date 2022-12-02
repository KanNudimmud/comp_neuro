%% PCA/SVD EXAMPLE
%% MONKEY REACHING DATA
%% loading data, making matrix
clear all, close all, clc
load('Lab5_CenterOutTrain.mat');

% get binned spikes for neuron 61
ni = 61; %35, 50*, 61*, 63, 76, 77*,78, 139
dt = .1; 

tvec = 0:dt:4;
n = size(tvec,2)-1; % number of timepoints
m = size(go,1); % number of trials

X = zeros(n,m);
for tri=1:size(go,1)
    for dti = 2:size(tvec,2)
        X(dti-1, tri) = sum(unit(ni).times > instruction(tri) + tvec(dti-1) ...
            & unit(ni).times < instruction(tri) + tvec(dti));
    end
end

% sorting according to reach direction
[~,ind] = sort(direction); 
X = X(:,ind);

%% plotting data matrix
figure(1); 
imagesc(X);shg
ylabel('timestep'); xlabel('trial')
set(gcf, 'Color', [1 1 1], 'papersize', [10 4], 'paperposition', [0 0 10 4])

%% subtract mean (in both dimensions)
mu=sum(X,2)/m;	
MU=repmat(mu,1,m);	
Z1=X-MU;	
mu=sum(Z1',2)/n;	
MU2=repmat(mu,1,n);	
Z=(Z1'-MU2)';	
figure; imagesc(Z);
ylabel('timestep'); xlabel('trial')

%% FIRST OPTION: TIMESTEP BASIS
%% compute covariance matrix
cov=Z*Z';	
imagesc(cov); 
xlabel('timestep'); 
ylabel('timestep'); 
set(gcf, 'Color', [1 1 1], 'papersize', [10 4], 'paperposition', [0 0 10 4])

%% compute eigenvectors and eigenvalues of cov matrix, plot variance
[F,V]=eig(cov);	
F = fliplr(F);
var=fliplr(sum(V));
figure(4); h = subplot(2,1,1);
plot(var, 'o'); %set(gca, 'yscale', 'log')
set(gcf, 'Color', [1 1 1], 'papersize', [5 4], 'paperposition', [0 0 5 4])
Fn = F; 

%% plot first 2 eigenvectors
figure(5); hold on
plot(F(:,1),'r')	
plot(F(:,2),'b')
set(gcf, 'Color', [1 1 1], 'papersize', [10 4], 'paperposition', [0 0 10 4])
xlabel('timestep'); 

%% Filter the data using the first PCs
Zf = F'*Z;
Zffilt=Zf;	
Zffilt(3:end,:)=0;	
Zflt=F*Zffilt;	
figure(5); 
subplot(2,2,1);
imagesc(Z); title('original')
subplot(2,2,3);
imagesc(Zflt); title('filtered by timestep PCs')
set(gcf, 'Color', [1 1 1], 'papersize', [5 4], 'paperposition', [0 0 5 4])

%% SECOND OPTION: TRIAL BASIS
%% plotting data matrix
figure(1); imagesc(X');shg
xlabel('timestep'); ylabel('trial')

%% compute covariance matrix
cov=Z'*Z;	
imagesc(cov); % now rows for each trial

%% compute eigenvectors and eigenvalues of cov matrix, plot variance
[F,V]=eig(cov);	
F = fliplr(F);
var=fliplr(sum(V));
figure(4); g = subplot(2,1,2);
plot(var, 'o'); %set(gca, 'yscale', 'log')
linkaxes([h g], 'x'); xlim([0 158])
linkaxes([h g], 'y')
Fm = F; 

%% plot first 2 eigenvectors
figure(3);
hold on
plot(F(:,1),'r')	
plot(F(:,2),'b')
axis tight
xlabel('trial')

%% Filter the data using the first PCs
Zf = F'*Z';
Zffilt=Zf;	
Zffilt(3:end,:)=0;	
Zflt=F*Zffilt;	
figure(5); 
subplot(2,2,2);
imagesc(Z'); title('original')
subplot(2,2,4);
imagesc(Zflt); title('filtered by trial PCs')

%% SVD
%% perform SVD
[U,S,V] = svd(Z);

%% Compare U and V to eigenvector matrices from PCA
% (look at abs of because some eigenvectors can have a sign flip)
figure; 
subplot(2,2,1); imagesc(abs(U)); title('U')
subplot(2,2,2); imagesc(abs(V(1:n,1:n))); title('V(1:n,1:n)')
subplot(2,2,3); imagesc(abs(Fn));title('Fn')
subplot(2,2,4); imagesc(abs(Fm(1:n,1:n)));title('Fm(1:n,1:n)')

%% Plot the first four matrix components
% each matrix is an outer product of basis vectors from U and V
figure(6);
for PCi = 1:4
    Comp = U(:,PCi)*S(PCi,PCi)*V(:,PCi)';
    subplot(2,2,PCi)
    title(num2str(PCi))
    imagesc(Comp); 
end
set(gcf, 'Color', [1 1 1], 'papersize', [5 4], 'paperposition', [0 0 5 4])

%% end.