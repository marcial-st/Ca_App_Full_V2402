%% Thesis
%%% Function that evaluates the system (A), observation matrix (H), system randomness (Q), and
%%% observation error covariance (R) matrix.
%%% where: 
%%%     1.- T   : sample period
%%%     2.- sx  : x position standard deviation 
%%%     3.- sy  : y position standard deviation 
%%%     4.- ra  : x random jerk 
%%%     5.- rb  : y random jerk
function [A,H,Q,R]=k_eval(T,sx,sy,ra,rb)
%% Constant Acceleration and Random Acceleration Changes 
% System 
 A=[1  T T^2/2  0  0   0   ;...% system
    0  1   T    0  0   0   ;...
    0  0   1    0  0   0   ;...
    0  0   0    1  T T^2/2 ;...
    0  0   0    0  1   T   ;...
    0  0   0    0  0   1   ];   
 
H=[1 0 0 0 0 0;... % observation matrix
   0 0 0 1 0 0];

w=[ra*T^3/3 ra*T^2/2 ra*T rb*T^3/3 rb*T^2/2 rb*T];
Q=w'*w;

R=[ sx^2    0   ;...% observation error covariance matrix
     0    sy^2 ];  % assuming independence
 