function [Ki_Norm, Ki_Hyp, SkMRglc_Norm, SkMRglc_Hyp] = skelFDGFit(bloodTAC, skelTAC, PG)
%skelFDGFit fits time-activity curve data to FDG 3 comparment model for the
% skeletal muscle. Allows for spillover (partial volume effect) from blood compartment 
% to the skeletal muscle ROI but not from the skeletal muscle compartment to the blood ROI. 
%
% Kinetic model based on publication:
% Dagan Feng, Xianjin Li, Sung-Cheng Huang. A new double modeling approach 
% for dynamic cardiac PET studies using noise and spillover contaminated LV 
% measurements. IEEE Trans Biomed Eng. 1996;43:319–327. doi: 10.1109/10.486290
%
%   Inputs: 
%       bloodTAC (time-activity curve for blood ROI)
%       skelTAC (time-activity curve for skeletal muscle ROI)
%       PG (plasma glucose in mg/dL)
%   Outputs: 
%       Ki_Norm (FDG influx constant during normoxia in 1/min)
%       Ki_Hyp (FDG influx constant during normoxia in 1/min)
%       SkMRglc_Norm (skeletal muscle metabolic rate of glucose during normoxia in umol/mL/min)
%       SkMRglc_Hyp (skeletal muscle metabolic rate of glucose during hypoxia in umol/mL/min) 

    % Call main FDG 3-compartment fit function NOT allowing
    % muscle-to-blood spillover, using a lumped constant of 1.16
    [Ki_Norm, Ki_Hyp, SkMRglc_Norm, SkMRglc_Hyp,~] = fdgFit(bloodTAC, skelTAC, PG, false, 1.16, false);
end

