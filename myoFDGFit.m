function [Ki_Norm, Ki_Hyp, MMRglc_Norm, MMRglc_Hyp,cp] = myoFDGFit(bloodTAC, myoTAC, PG)
%myoFDGFit fits time-activity curve data to FDG 3 comparment model for the
% myocardium. Allows for spillover (partial volume effect) from blood compartment 
% to the myocardium ROI and from the myocardium compartment to the blood ROI. 
%
% Kinetic model based on publication:
% Dagan Feng, Xianjin Li, Sung-Cheng Huang. A new double modeling approach 
% for dynamic cardiac PET studies using noise and spillover contaminated LV 
% measurements. IEEE Trans Biomed Eng. 1996;43:319–327. doi: 10.1109/10.486290
%
%   Inputs: 
%       bloodTAC (time-activity curve for blood ROI)
%       myoTAC (time-activity curve for myocardium ROI)
%       PG (plasma glucose in mg/dL)
%   Outputs: 
%       Ki_Norm (FDG influx constant during normoxia in 1/min)
%       Ki_Hyp (FDG influx constant during normoxia in 1/min)
%       MMRglc_Norm (myocardial metabolic rate of glucose during normoxia in umol/mL/min)
%       MMRglc_Hyp (myocardial metabolic rate of glucose during hypoxia in umol/mL/min) 
%       cp (fitted concentration of FDG in plasma compartment)

    % Call main FDG 3-compartment fit function with allowing
    % myocardium-to-blood spillover, using a lumped constant of 1.44
    [Ki_Norm, Ki_Hyp, MMRglc_Norm, MMRglc_Hyp,cp] = fdgFit(bloodTAC, myoTAC, PG, true, 1.44, false);
end

