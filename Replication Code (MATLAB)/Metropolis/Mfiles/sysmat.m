function [A,B,H,R,Se,Phi] = sysmat(T1,T0,para)

% This function computes the matrices of the state space representation.
% Input = para : Vector of Structural Parameters
%         T1 T0: 
% Output= Matrices of state space model. See kalman.m


nu       = para(1);
alp      = para(2);
del      = 0.025;
bet      = 0.99;
gamA     = para(3);
rhoA     = para(4);
rhoB     = para(5);
lnY0     = para(6);
lnhst    = para(7);
sigA     = para(8);
sigB     = para(9);

eq_y = 1;
eq_h = 2;

% /** number of observation variables **/

ny = 2;

% /** model variable indices **/

r_t   = 1;
w_t   = 2;
c_t   = 3;
k_t1  = 4;
h_t   = 5;
y_t   = 6;
i_t   = 7;
Er_t  = 8;
Ec_t  = 9;
Eda_t = 10;
da_t  = 11;
b_t   = 12;
a_t   = 13;

% /** shock indices **/

e_a  = 1;
e_b  = 2;

%=========================================================================
%                          TRANSITION EQUATION
%  
%           s(t) = Phi*s(t-1) + R*e(t)
%           e(t) ~ iid N(0,Se)
% 
%=========================================================================

nep = size(T0,2);

Phi = T1;

R   = T0;

Se  = zeros(nep,nep);

Se(e_a,e_a) = sigA^2;

Se(e_b,e_b) = sigB^2;

%=========================================================================
%                          MEASUREMENT EQUATION
%  
%           y(t) = a + b*s(t) + u(t) 
%           u(t) ~ N(0,HH)
% 
%=========================================================================

A         = zeros(ny,1);
A(eq_y,1) = lnY0;
A(eq_h,1) = lnhst;

nstate = size(Phi,2); 

B = zeros(ny,nstate);

B(eq_y,y_t) =  1;
B(eq_y,a_t) =  1; 

B(eq_h,h_t) =  1;

H = zeros(ny,ny);  

end

