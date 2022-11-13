%% Exercise 6: Singular Value Decomposition (SVD) on images of faces
clear all, close all, clc
% We’ll perform SVD on images of faces. 
% SVD decomposes a matrix X into the product USV’,
% where U and V are square matrices and S is diagonal. 
% If X is n x m, columns of U are a basis in n-dimensional space, 
% and columns of V are a basis in m-dimensional space. 
% If you ve meansubtracted your data, 
% these are the same bases you get in PCA 
% (the eigenvectors of the covariance matrices, XX  and X’X respectively).
% S contains the squared eigenvalues of these covariance matrices 
% on the diagonal (XX  and X’X have the same eigenvalues). 
% Using SVD, we can look at the principal modes in both the pixel space

% We’ll use faces from the JAFFE database, which contains images 
% of several subjects each displaying several facial expressions 
% (Michael J. Lyons. Japanese Female Facial Expressions (JAFFE),
% Database of digital images (1997). http://www.kasrl.org/jaffe_info.html)
load( 'jaffe.mat' )

% i=1;
% figure()
% imagesc(reshape(IMS(i,:),137,86));
% colormap gray;

% Mean-subtract your data so that 
% the mean of each pixel is 0 and the mean of each image is zero.
meanImages = mean(IMS, 2);
IMS = IMS - repmat(meanImages, 1, size(IMS, 2));
for i = 1:size(IMS, 1)
    IMS(i, :) = IMS(i, :) ./ var(IMS(i, :));
end

% compute U, S and V
[U, S, V] = svd(IMS);

% confirm IMS = U*S*V'

figure
plot(diag(S))
xlabel('ranking singular value label')
ylabel('singular value')

figure
for i = 1:9
    subplot(3,3,i)
    imagesc(reshape(V(:, i), 137, 86));
    colormap gray
    title(['singular value ', num2str( i )])
    axis equal
end

n = 9;
VRedDim = V;
for i = (n + 1): size(V, 2)
    VRedDim(:, i) = 0;
end
I = U * S * VRedDim';

numImages = 5;
figure
count = 1;
for j = 1 : numImages
    subplot(numImages, 2, count)
    count = count + 1;
    imagesc(reshape(IMS( j, : ), 137, 86));
    axis equal
    colormap gray
    subplot(numImages, 2, count)
    count = count + 1;
    imagesc(reshape(I(j, :), 137, 86));
    colormap gray
    axis equal
    
end

% U is a basis for the space of images, of which there are 213. 
% The variable EMind indicates the emotion of each image, 
% and IDind indicates the subject id of each image.