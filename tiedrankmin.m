% tiedrankmin.m
%
% This function returns the sample ranks of the values in a vector. Tied
% values will be replaced by their minimum.
%
% Input parameter:
%   - x  : a vector
%
% This function is adapted from tiedrank2 function available at
% https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fww2.mathworks.cn%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fsubmissions%2F28528%2Fversions%2F5%2Fcontents%2Ftiedrank2.m&embed=web

function [r, tieadj,tiecount] = tiedrankmin(x)
%TIEDRANK Compute the ranks of a sample, adjusting for ties.
%   [R, TIEADJ] = TIEDRANK(X) computes the ranks of the values in the
%   vector X.  If any X values are tied, TIEDRANK computes their average
%   rank.  The return value TIEADJ is an adjustment for ties required by
%   the nonparametric tests SIGNRANK and RANKSUM, and for the computation
%   of Spearman's rank correlation.
%
%   [R, TIEADJ] = TIEDRANK(X,1) computes the ranks of the values in the
%   vector X.  TIEADJ is a vector of three adjustments for ties required
%   in the computation of Kendall's tau.  TIEDRANK(X,0) is the same as
%   TIEDRANK(X).
%
%   [R, TIEADJ] = TIEDRANK(X,0,1) computes the ranks from each end, so
%   that the smallest and largest values get rank 1, the next smallest
%   and largest get rank 2, etc.  These ranks are used for the
%   Ansari-Bradley test.
%
%   See also ANSARIBRADLEY, CORR, PARTIALCORR, RANKSUM, SIGNRANK.
%   Copyright 1993-2005 The MathWorks, Inc.
%   $Revision: 1.5.2.6 $  $Date: 2005/07/29 11:42:04 $

if isvector(x)
   [r,tieadj,tiecount] = tr(x);
else
   if isa(x,'single')
      outclass = 'single';
   else
      outclass = 'double';
   end
   % Operate on each column vector of the input (possibly > 2 dimensional)
   sz = size(x);
   ncols = sz(2:end);  % for 2x3x4, ncols will be [3 4]
   r = zeros(sz,outclass);
   
   tieadj = zeros([1,ncols],outclass);
   
   for j=1:prod(ncols)
      [r(:,j),tieadj(:,j)] = tr(x(:,j));
   end
end
% --------------------------------
function [r,tieadj,tiecount] = tr(x)
%TR Local tiedrank function to compute results for one column
% Sort, then leave the NaNs (which are sorted to the end) alone
[sx, rowidx] = sort(x(:));
numNaNs = sum(isnan(x));
xLen = numel(x) - numNaNs;

ranks = [1:xLen NaN(1,numNaNs)]';

tieadj = 0;

if isa(x,'single')
   ranks = single(ranks);
   tieadj = single(tieadj);
end
% Adjust for ties.  Avoid using diff(sx) here in case there are infs.
ties = (sx(1:xLen-1) == sx(2:xLen));
tieloc = [find(ties); xLen+2];
maxTies = numel(tieloc);
tiecount = 1;
while (tiecount < maxTies)
    tiestart = tieloc(tiecount);
    ntied = 2;
    while(tieloc(tiecount+1) == tieloc(tiecount)+1)
        tiecount = tiecount+1;
        ntied = ntied+1;
    end
    
    tieadj = tieadj + ntied*(ntied-1)*(ntied+1)/2;
    
    % Compute min of tied ranks
    ranks(tiestart:tiestart+ntied-1) = ranks(tiestart);
    tiecount = tiecount + 1;
end
% Broadcast the ranks back out, including NaN where required.
r(rowidx) = ranks;
r = reshape(r,size(x));
