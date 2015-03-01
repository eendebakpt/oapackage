function [  ] = Doptim( varargin )
%DOPTIM Stub for mex file
% Usage: Doptim(N, levels, [nrestarts], [alpha], [verbose], [maxtime])
%    calculates D-optimal designs),    written by Pieter.Eendebak <pieter.eendebak@gmail.com>
% 
% Input arguments:
%    N (integer): number of runs
%    s (array): levels of the factors
%    nrestarts (integer, default is 40): number of restarts
%    alpha (3x1 array, default is [1,1,0]): weights for the various efficiency measures
%    verbose (integer, default is 1): output level
%    maxtime (double, default: 150): maximum running time in seconds
% Output arguments:
%    A (Nxk or Nxkxm array): generated design(s)
%    d (3x1 or 3xm array): computed values
% 
% The function optimizes design according to the optimization function:
%
% 	F = alpha(1)*D + alpha(2) * Ds + alpha(2) * D1
%
% For more details see the webpage http://www.pietereendebak.nl/oapackage.
% 
% Example: [A, d] = Doptim(40, [2,2,2,2,2,2], 10, [1,2,0])

error('please compile the mex file Doptim'); 
return

%% Testing

tic;
[A, d] = Doptim(40, [2,2,2,2,2,2], 10, [1,2,0]);
toc;
d

%%
[A, d] = Doptim(48, [4,3,2,2,2], 10, [1,2,0], 2, 100,2); d
disp(A)
disp(d)

%%
[A, d] = Doptim(40, [2,2,2,2,2,2], 10, [1,2,0], 2, 100,2); d



