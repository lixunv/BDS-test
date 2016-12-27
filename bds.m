function [w, sig, c, c1, k] = bds (series, maxdim, distance, flag, maxram)

if nargin < 5
maxram = 150;
elseif maxram > 500
disp('Are you sure you have so much memory available?')
error('If so, you need to edit the code, otherwise try again with a lower value.')
end
if nargin < 4
flag = 0;
elseif ~any(flag == [0 1])
error('Unknown method for determining dimensional distance; try again with 0 or 1.')
end
if nargin < 3
distance = 1.5;
elseif distance < 0
error('The dimensional distance parameter must be positive.')
elseif flag == 1 & distance > 1
error('The correlation integral cannot exceed 1.')
end
if nargin < 2
maxdim = 2;
elseif maxdim < 1
error('The dimension needs to be at least 1.');
end
if nargin < 1
error('Cannot compute the BDS statistic on nothing.')
end
[rows,cols] = size(series);
if rows > 1 & cols == 1
n = rows;
series = series';
elseif cols > 1 & rows == 1
n = cols;
elseif cols > 1 & rows > 1
n = cols*rows;
series = series(:)'; % transformation into a row vector
disp(sprintf('\aTransformed matrix input into a single column.'))
else
error('Cannot compute the BDS statistic on a scalar!')
end
%%%%%%%%%%%% Determination of and preparations for fastest method given MAXRAM %%%%%%%%%%%
fastbuild = 0.000016 * (1:52) .* pow2(1:52); % memory requirements
slowbuild = 0.000045 * pow2(1:52); % for the various
holdinfo = 0.000005 * pow2(1:52); % algorithms in 
wordtable = 0.000008 * n^2 ./ (1:52); % megabytes for 
bitandop = 0.000024 * n^2 ./ (1:52); % given N
[ram1, bits1] = min(fastbuild + holdinfo + wordtable + bitandop); % number of bits for
[ram2, bits2] = min(fastbuild + holdinfo + wordtable); % which each of six
[ram3, bits3] = min(slowbuild + holdinfo + wordtable + bitandop); % methods uses minimum
[ram4, bits4] = min(slowbuild + holdinfo + wordtable); % memory; this memory
[ram5, bits5] = min( wordtable + bitandop); % is given by
[ram6, bits6] = min( wordtable); % ram1, ram2,..., ram6
if ram1 < maxram | ram2 < maxram
if ram1 < maxram
method = 1;
bits = bits1; ram = ram1;
else
method = 2; % maximum number of rows to put
bits = bits2; ram = ram2; % through BITAND and bit-counting
stepping = floor((maxram-ram)*bits/n/0.000024); % algorithm without exceeding MAXRAM
end
% Vector BITINFO lists the number of bits set for each integer between 0 and 2^bits
% (corresponding to the indices of the vector shifted by 1). See Kanzler (1998) for an
% explanation.
bitinfo = uint8(sum(rem(floor((0:pow2(bits)-1)'*pow2(1-bits:0)),2),2));
elseif ram3 < maxram | ram4 < maxram
if ram3 < maxram
method = 3;
bits = bits3; ram = ram3;
else
method = 4;
bits = bits4; ram = ram4;
stepping = floor((maxram - ram) * bits / n / 0.000024);
end
bitinfo(1:pow2(bits), :) = uint8(0); % the same as above, but created through
for bit = 1 : bits % a loop, which consumes less memory
bitinfo(1:pow2(bits)) = sum([bitinfo, ...
kron(ones(pow2(bits-bit),1), [zeros(pow2(bit-1),1); ones(pow2(bit-1),1)])],2);
end
elseif ram5 < maxram | ram6 < maxram
if ram5 < maxram
method = 5;
bits = bits5; ram = ram5;
else
method = 6;
bits = bits6; ram = ram6;
stepping = floor((maxram - ram) * bits / n / 0.000024);
end
else
disp('Insufficient amount of memory. Allocate more memory to the system')
disp('or reduce the number of observations, then try again.')
error(' ')
end
%
if ~flag
demeaned = series-sum(series)/n; % fastest algorithm for
epsilon = distance * sqrt(demeaned*demeaned'/(n-1)); % computing the standard
clear demeaned % to save memory % deviation of SERIES

elseif 0.000008 * 3 * sum(1:n-1) < maxram % check memory requirements for DIST and sorting
dist(1:sum(1:n-1)) = 0;
for i = 1 : n-1
dist(1+(i-1)*(n-1)-sum(0:i-2):i*(n-1)-sum(1:i-1)) = abs(series(i+1:n)-series(i));
end
sorted = sort(dist);
epsilon = sorted(round(distance*sum(1:n-1))); % DISTANCEth percentile of SORTED series
clear dist sorted
else
error('Insufficient RAM to compute EPSILON; allocate more memory or use METHOD = 1.')
end
colsum(1:n) = 1;
rowsum(1:n) = 0;
nwords = ceil((n-1)/bits);
wrdmtrx(1:n-1,1:nwords) = 0; % initialisation of bit-word table
for row = 1 : n-1
bitvec = abs(series(1+row:n) - series(row)) <= epsilon;
rowsum(row) = sum(bitvec);
colsum(1+row:n) = colsum(1+row:n) + bitvec;
nwords = ceil((n-row)/bits);
wrdmtrx(row,1:nwords) = (reshape([bitvec,zeros(1,nwords*bits-n+row)],...%transformation
bits, nwords)' *pow2(0:bits-1)')'; %into bit-words
end
clear series bitvec
bitsum(maxdim:-1:1) = cumsum([sum(rowsum(maxdim:n-1)), rowsum(maxdim-1:-1:1)]);
c1 (maxdim:-1:1) = bitsum(maxdim:-1:1) ./ cumsum([sum(1:n-maxdim), n-maxdim+1 : n-1]);
fullsum = rowsum + colsum;
k = (fullsum*fullsum' + 2*n - 3*(2*bitsum(1)+n)) / n/(n-1)/(n-2);
clear rowsum colsum fullsum bitsum
for m = 2 : maxdim
bitcount = 0;
if sum(method == [1 3])
wrdmtrx(m:n-1,:) = bitand(wrdmtrx(m:n-1,:),wrdmtrx(m-1:n-2,:)); % BITAND and bit
bitcount = sum(sum(bitinfo(wrdmtrx(m:n-1,:)+1))); % count all at once
elseif sum(method == [2 4])
for row = n-stepping : -stepping : m+1 % BITAND
wrdmtrx(row:row+stepping-1,:) = bitand(wrdmtrx(row:... % and bit
row+stepping-1,:), wrdmtrx(row-1:row+stepping-2,:)); % count in
bitcount=bitcount+sum(sum(bitinfo(wrdmtrx(row:row+stepping-1,:)+1))); % backward
end % loops
wrdmtrx(m:row-1,:) = bitand(wrdmtrx(m:row-1,:), wrdmtrx(m-1:row-2,:)); % through
bitcount = bitcount + sum(sum(bitinfo(wrdmtrx(m:row-1,:)+1))); % the table
elseif method == 5
wrdmtrx(m:n-1,:) = bitand(wrdmtrx(m:n-1,:), wrdmtrx(m-1:n-2,:)); % BITAND at once...
for col = 1 : ceil((n-1)/bits) % bit count
bitcount = bitcount + sum(sum(rem(floor(wrdmtrx(m:... % by brute force
n-1-(col-1)*bits, col) * pow2(1-bits:0)), 2))); % in loops
end
else
for row = n-stepping : -stepping : m+1
wrdmtrx(row:row+stepping-1,:) = bitand(wrdmtrx(row:... % BITAND
row+stepping-1,:), wrdmtrx(row-1:row+stepping-2,:)); % opera-
end % tions
wrdmtrx(m:row-1,:) = bitand(wrdmtrx(m:row-1,:), wrdmtrx(m-1:row-2,:)); % and brute-
for col = 1 : ceil((n-1)/bits) % force bit
bitcount = bitcount + sum(sum(rem(floor(wrdmtrx(m:... % counting
n-1-(col-1)*bits, col) * pow2(1-bits:0)),2))); % in loops
end
end
c(m-1) = bitcount / sum(1:n-m); % indexing of
sigma(m-1) = 2*sqrt(prod(ones(1,m)*k) + 2*ivp(k,m-(1:m-1),m-1)... % C and SIGMA
*(ivp(c1(1),2*(1:m-1),m-1))' + (m-1)*(m-1)... % runs from 1
*prod(ones(1,2*m)*c1(1)) - m*m*k*prod(ones(1,2*m-2)*c1(1))); % to MAXDIM-1
end
clear wrdmtrx
if maxdim > 1
w = sqrt(n-(2:maxdim)+1) .* (c - idvp(c1(2:maxdim), 2:maxdim, maxdim-1)) ./ sigma;
if exist('normcdf.m','file') & nargout > 1
sig = min(normcdf(w,0,1), 1-normcdf(w,0,1)) * 2;
elseif nargout > 1
sig(1:maxdim-1) = NaN;
end
else
w = [];
sig = [];
c = [];
end
%%%%%%%%%%%%%%%%%%%%%%% Sub-functions for computing integer powers %%%%%%%%%%%%%%%%%%%%%%%
function ipow = ivp (base, intpowvec, veclen)
ipow(1 : veclen) = 0;
for j = 1 : veclen
ipow(j) = prod(ones(1, intpowvec(j)) * base);
end
function ipow = idvp (basevec, intpowvec, veclen)
ipow(1 : veclen) = 0;
for j = 1 : veclen
ipow(j) = prod(ones(1, intpowvec(j)) * basevec(j));
end