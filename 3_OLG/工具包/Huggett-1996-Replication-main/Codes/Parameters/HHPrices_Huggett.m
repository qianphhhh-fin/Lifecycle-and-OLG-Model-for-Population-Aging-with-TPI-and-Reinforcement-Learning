function [Y, R, w, b] = HHPrices_Huggett(K, L, cS)
%% Documentation:
% This function is defined based on production function

% INPUTS
% (1). K:        aggregate capital in the economy
% (2). L:        aggregate labor supply
%                This term can be precalibrated without solving HH problem
% (3). cS.A:     productivity
% (4). cS.alpha: capital share
% (5). cS.ddk:   depreciation rate
% (6). cS.theta: social security tax
% (7). cS.retireMass: mass of retired households

% OUTPUTS
% (1). Y: production output
% (2). R: capital rental price faced by households
%         R = 1+ (MPK - depreciation)(1 - tax)
% (3). w: wage after tax
%         w = MPL(1 - tax - social_security)
% (4). b: social security benefit
%         (social_security x w x L) = b x all_retired_households

% ******************* Note ******************* 
% tax is calibrated in this function as 
% 0.195/(1-depreciation x K/Y)


%% Main
% Production
Y   = cS.A * (K^cS.alpha) * (L^(1-cS.alpha));
MPK = cS.alpha * cS.A * (K^(cS.alpha-1)) * (L^(1-cS.alpha));
MPL = (1-cS.alpha) * cS.A * (K^cS.alpha) * (L^(-cS.alpha));

% Tax
tau = 0.195/(1-cS.ddk * K/Y);

% Prices faced by households
R   = 1 + (MPK - cS.ddk)*(1 - tau);
w   = MPL*(1 - tau - cS.theta);

% Social security benefits
b   = cS.theta * w * L/cS.retireMass;


end