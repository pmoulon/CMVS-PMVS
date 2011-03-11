function [A, m, n] = readMetis(fname)
% read from metix input graph

fid = fopen(fname,'r');
if fid == -1
   fprintf('readDense: fopen for file %s failed\n',fname); 
   error('Fatal Error');
end 

m = fscanf(fid,'%d', 1); % number of vertices
n = fscanf(fid,'%d', 1); % number of edges
%A = zeros(m,m);
A(1,1)=0;
A=sparse(A);
t = fgetl(fid);
format = sscanf(t,'%d');
if ((length(format) ==0) | (format ==0))
    for i=1:m
        t = fgetl(fid);
        l = sscanf(t, '%d');
        A(i,l) = 1;
    end
else
    for i=1:m
        t = fgetl(fid);
        l = sscanf(t, '%d');
        
        A(i,l(1:2:length(l))) = l(2:2:length(l));
    end
end
%A = sparse(A);