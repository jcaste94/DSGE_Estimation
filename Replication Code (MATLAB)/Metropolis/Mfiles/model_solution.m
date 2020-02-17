function [T1, TC, T0, TETA, RC, retcode] = model_solution(para)

% Solve the DSGE model in Illustration 4.1 of Del Negro and Schorhfeide (2012)
%
% INPUT
% para:  structural parameters
%
% OUTPUT
% para:  structural parameters
%
% DATE: 6/16/2013
% 


%=========================================================================
%                     COMPUTE STEADY STATE
%=========================================================================

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

qst      = gamA;

rst      = (gamA/bet)-(1-del);
k_y      = ((1-alp)/rst)*gamA;
i_y      = (1-((1-del)/gamA))*k_y;
i_k      = i_y/k_y;
c_y      = 1-i_y;

sst  = [rst;k_y;i_y;i_k;c_y];

%=========================================================================
%                        SOLVE DSGE MODEL
%=========================================================================


retcode = 0;

valid   = 1;

% /* Variable indices */

%=========================================================================
%                       DEFINE OBJECTS
%=========================================================================

% Equation indices

eq_r    = 1;  %** Firm FOC for r(t) **/
eq_w    = 2;  %** Firm FOC for w(t) **/
eq_k    = 3;  %** Household FOC for k(t+1) **/
eq_h    = 4;  %** Household FOC for h(t) **/
eq_mc   = 5;  %** Market clearing **/
eq_y    = 6;  %** Production Function **/
eq_kacc = 7;  %** Capital Accumulation **/
eq_Er   = 8;  %** E(r(t+1)) **/
eq_Ec   = 9;  %** E(c(t+1)) **/
eq_Eda  = 10;  %** E(da(t+1)) **/ 
eq_a    = 11;  %** a process **/
eq_b    = 12;  %** b process **/
eq_da   = 13;  %** growth a **/

% Variable indices 

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

%Expectation error indices (eta) 

er_sh  = 1;
ec_sh  = 2;
eda_sh = 3;
       
%Shock indices (eps)

a_sh = 1;
b_sh = 2;

%SUMMARY

neq  = 13;
neta = 3;
neps = 2;

% /** initialize matrices **/

GAM0 = zeros(neq,neq);
GAM1 = zeros(neq,neq);
   C = zeros(neq,1);        
 PSI = zeros(neq,neps);
 PPI = zeros(neq,neta);

% /** Steady State **/

Rst = sst(1);
k_y = sst(2);
i_y = sst(3);
i_k = sst(4); 
c_y = sst(5);

%=========================================================================
%                EQUILIBRIUM CONDITIONS: CANONICAL SYSTEM
%=========================================================================

%=========================================================================
%         1. Firm FOC for r(t)
%=========================================================================

GAM0(eq_r,r_t)   =  1;
GAM0(eq_r,y_t)   = -1;
GAM0(eq_r,da_t)  = -1;
GAM1(eq_r,k_t1)  = -1;
   

%=========================================================================
%         2. Firm FOC for w(t)
%=========================================================================

GAM0(eq_w,w_t)   =  1;
GAM0(eq_w,y_t)   = -1;
GAM0(eq_w,h_t)   =  1;

%=========================================================================
%         3. HH FOC for k(t+1)
%=========================================================================

GAM0(eq_k,c_t)   = 1;
GAM0(eq_k,Er_t)  = Rst/(Rst+(1-del));
GAM0(eq_k,Ec_t)  = -1;
GAM0(eq_k,Eda_t) = -1;
  
%=========================================================================
%         4. HH FOC for h(t)
%=========================================================================

GAM0(eq_h,w_t)   = -nu;
GAM0(eq_h,c_t)   =  nu;
GAM0(eq_h,h_t)   =  1;
GAM0(eq_h,b_t)   =  -(1+nu);

%=========================================================================
% **      5. Market Clearing Condition
%=========================================================================

GAM0(eq_mc,c_t)   = -c_y;
GAM0(eq_mc,i_t)   = -i_y;
GAM0(eq_mc,y_t)   = 1;

%=========================================================================
%         6. Production Function
%=========================================================================

GAM0(eq_y,y_t)   =  1;
GAM0(eq_y,h_t)   = -alp;
GAM0(eq_y,da_t)  = (1-alp); 
GAM1(eq_y,k_t1)  = (1-alp);
   
%=========================================================================
%         7. Capital Accumulation
%=========================================================================

GAM0(eq_kacc,k_t1)  = 1;
GAM0(eq_kacc,i_t)   = -i_k;
GAM0(eq_kacc,da_t)  = (1-del)/qst;
GAM1(eq_kacc,k_t1)  = (1-del)/qst;
      
%=========================================================================
%          Expectation error
%=========================================================================

% /** E(r) **/

GAM0(eq_Er,r_t)   = 1;
GAM1(eq_Er,Er_t)  = 1;
 PPI(eq_Er,er_sh) = 1;

% /** E(c) **/

GAM0(eq_Ec,c_t)   = 1;
GAM1(eq_Ec,Ec_t)  = 1;
 PPI(eq_Ec,ec_sh) = 1;

% /** E(da) **/

GAM0(eq_Eda,da_t)   = 1;
GAM1(eq_Eda,Eda_t)  = 1;
 PPI(eq_Eda,eda_sh) = 1;


%=========================================================================
% **      Shock process
%=========================================================================

% /** a(t) **/

GAM0(eq_a,a_t) = 1;
GAM1(eq_a,a_t) = rhoA;
 PSI(eq_a,a_sh)= 1;

% /** da **/

GAM0(eq_da,da_t) = 1;
GAM0(eq_da,a_t)  = -1;
GAM1(eq_da,a_t)  = -1;

% /** b **/

GAM0(eq_b,b_t)  = 1;
GAM1(eq_b,b_t)  = rhoB;
 PSI(eq_b,b_sh) = 1;

%=========================================================================
%           QZ(generalized Schur) decomposition by GENSYS
%=========================================================================

[T1,TC,T0,TY,M,TZ,TETA,GEV,RC] = gensys(GAM0,GAM1,C,PSI,PPI,1+1E-8);

end
