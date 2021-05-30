%{   
 Suppose we have G(Z)
                          z^-1 (0.1459 + 0.344 z^-1 + 0.0445 z^-2)
                Gz =  -----------------------------------------------
                        1 - 1.806 z^-1 + 0.8964 z^-2 - 0.09072 z^-3
Solve using 2nd config. & with Dynamic Cancellation 
%}
clc; clear vars; close all;
%% Plant ( G(z) )
B = [0.1459 0.344 0.04451];           % numerator 
A = [1 -1.806 0.8964 -0.09072];       % denominator
d = 1;                                % order of delay of the system
%% Inputs
syms t z 
U(t) = exp(-0.04*t)*sin(0.2*t);
U(z) = vpa(ztrans(U(t)),2);
[N,M]= numden(U(z));
M    = sym2poly(M);
M    = M/M(1,1)
%% Poles 
Am = [1 -0.4177 0.0183];
A0 = [1 -0.4 0.04];
alpha = conv(Am,A0);
%% Call Diophantine solver function
[T,v] = Diophantine(M,B,d,alpha)
%% Diophantine solver function
function [S,R] = Diophantine(A,B,d,alpha)
%{
This function is intended to solve the Diphantine equation in the form of 
                            AR + z^(-d) BS = A0 Am = alpha;
where
- A = 1 + a_1 z^-1 + a_2 z^-1 + ... + a_na z^(-na) 
- B = b_0 + b_1 z^-1 + b_2 z^-1 + ... + b_nb z^(-nb) 
- R = 1 + r_1 z^-1 + r_2 z^-1 + ... + r_nr z^(-nr) 
- S = s_0 + s_1 z^-1 + s_2 z^-1 + ... + s_ns z^(-ns) 
- d : delay in the system. Notice that this form of the Diaphontaing solution
        is available for systems with d>=1
- alpha = 1 + alpha1 z^-1 + alpha2 z^-1 + ... + alpha_(nalpha z)^(-nalpha) = Am*A0,  required characteristic polynomial
- Am = required polynomial of the model;
- A0 = observer polynomail for compensation of the order
The function input outputs are given in the following
                        function [ S, R ] = Diophantine( A, B, d, alpha )
Inputs
A = [1,  a_1,  a_2,  a_3,  ..., a_na] 
B = [b_0,  b_1, b _2,  b_3,  ..., a_nb] 
d = delay time, a number.
alpha =  [1,  alpha_1, alpha _2,  alpha_3,  ..., alpha_nalpha], nalpha is
the final order of the closed loop transfer function
Outputs
S = [s_0,  s_1, s _2,  s_3,  ..., s_ns]
R =  [1,  r_1,  r_2,  r_3,  ..., r_nr]
to find the oreders of the polynomials we use these equations 
nr = nb + d - 1
ns = na - 1 
nalpha = na + nb + d - 1
the functions is used to estimate the polynomials S and R which are the
numerator and the denomenator of the controller transfer function,
respectively.
The Solution is given in matrix form by solving a linear system of
equations such as 
          ####### M*theta = (V-Y) --> theta = M^(-1)*(V-Y) #######
-- M : Sylvester matrix
-- V: vector contains the "alpha" polynomail coefficients without "1" at the
        first of it. 
        
                            V = transpose([alpha_1, alpha _2,  alpha_3,  ..., alpha_nalpha])
                            size(nalpha, 1)
-- Y: vector contains the "A" polynomail coefficients without "1" at the
        first of it.
                            Y = transpose([a_1,  a_2,  a_3,  ..., a_na, 0, 0, ..., 0])
                            size(nalpha, 1)
-- theta: vector contains the unknowns. That is, the coefficients of the R
        polynomial and the coefficients of the S polynomial
    
                            theta = tranpose([r_1, r_2, ..., r_nr, s_0, s1, s2, ..., s_ns])
 %}
%% Calculating orders
na = length(A)-1;
nb = length(B)-1;
nr = nb + d - 1;
ns = na - 1;
nalpha = na + nb +d - 1;
%% initialization of M matrix and V matrix
M = zeros(nalpha, nalpha);
Y = zeros(nalpha, 1);
V = zeros(nalpha, 1);
Y(1:na,1)     = A(2:end);
V(1:nalpha,1) = alpha(2:end);
%% M matrix filling
for i = 1: nalpha
    if(i<=nr)
        M(i:(length(A)+i-1), i) = A;
    else
        M(d+(i-(nr+1)):(d+(i-(nr+1)))+nb, i) = B;
    end
end
%% Solution of the Diophantine equation to find S, R polynomials 
theta = M \(V-Y);
R = [1; theta(1:nr)]';
S = theta((nr+1):end)';
end