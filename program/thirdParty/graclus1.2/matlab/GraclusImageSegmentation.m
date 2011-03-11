function GraclusImageSegmentation(imgFile,k)
%% function GraclusImageSegmentation(imgFile,k)
%
%  Inputs: imgFile (image to segment), k (number of segs)
%
%  Note: Make sure you add a path to Jianbo Shi's image segmentation
%  code within this function so that the graph can be created from
%  the image (via the function ICgraph).

I = imread(imgFile);
[Inr,Inc,nb] = size(I);

if (nb>1),
    I =double(rgb2gray(I));
else
    I = double(I);
end

%If you want to resize to nr x nc to make the image smaller:
%I = imresize(I,[nr, nc],'bicubic');
%Inr = nr; Inc = nc;

disp('Creating the graph...');
[W,imageEdges] = ICgraph(I);
maxW = max(max(W));
spectral = 1;
[n,n] = size(W);
W = sparsifyc(W,1e-6);
W2 = ceil(1e2*W);
[partition,ncutval]=graclus_mex(W2,nnz(W2),k,0,20,spectral);
GraclusDiscrete = zeros(n,k);
for i = 1:k
    idx = find(partition == i);
    for j = 1:length(idx)
        GraclusDiscrete(idx(j),i) = 1;
    end
end
GraclusLabel = zeros(Inr,Inc);
for j=1:size(GraclusDiscrete,2),
    GraclusLabel = GraclusLabel + j*reshape(GraclusDiscrete(:,j),Inr,Inc);
end

%Display the segmentation
figure(1);clf; imagesc(I);colormap(gray);axis off;
disp('Figure 1: The input image before segmentation....');

figure(2);clf
bw = edge(GraclusLabel,0.01);
J1=showmask(I,imdilate(bw,ones(2,2))); imagesc(J1);axis off
disp('Figure 2: The segmentation using Graclus.');
