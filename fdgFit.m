function [Ki_Norm, Ki_Hyp, MRglc_Norm, MRglc_Hyp,cp] = fdgFit(bloodTAC, tissueTAC, PG, fitblood, LC, plotOn)
%fdgFit fits time-activity curve data to FDG 3 comparment model. Allows for 
% spillover (partial volume effect) from blood compartment to the tissue ROI 
% and from the tissue compartment to the blood ROI. 
%
% Kinetic model based on publication:
% Dagan Feng, Xianjin Li, Sung-Cheng Huang. A new double modeling approach 
% for dynamic cardiac PET studies using noise and spillover contaminated LV 
% measurements. IEEE Trans Biomed Eng. 1996;43:319–327. doi: 10.1109/10.486290
%
%   Inputs: 
%       bloodTAC (time-activity curve for blood ROI)
%       tissueTAC (time-activity curve for tissue ROI)
%       PG (plasma glucose in mg/dL)
%       fitblood (boolean value indicating whether to allow spillover from tissue compartment to blood ROI)
%       LC (lumped constant to convert from FDG metabolism to glucose metabolism)
%       plotOn (boolean value indicating whether to plot kinetic model fit)
%   Outputs: 
%       Ki_Norm (FDG influx constant during normoxia in 1/min)
%       Ki_Hyp (FDG influx constant during normoxia in 1/min)
%       MRglc_Norm (metabolic rate of glucose during normoxia in umol/mL/min)
%       MRglc_Hyp (metabolic rate of glucose during hypoxia in umol/mL/min) 
%       cp (fitted concentration of FDG in plasma compartment)

    if nargin<6; plotOn=false; end % default is to not plot fit
    
    nT = numel(bloodTAC);
    options = optimoptions('fmincon'); % Initialize minimization options 
    options.MaxFunctionEvaluations = 1000000;
    options.MaxIterations = 5000;
    
    % Objective function to minimize: (tissueFit-tissueROI)^2 + (bloodFit-bloodROI)^2
    % tissueFit = aCb + (1-a)(Cf + Cm)
    % bloodFit = bCb + (1-b)(Cf + Cm)
    % Cb, Cf, and Cm are the three compartments in the model (blood, tissue-free, tissue-metabolized)
    % a and (1-b) are the blood-to-tissue and the tissue-to-blood
    % spillover fractions, respectively. 
    %
    % x is the variable to be minimized by fmincon, and it contains the
    % following array: [a, b, k1_norm, k1_hyp, k2_norm, k2_hyp, k3_norm, k3_hyp, Cb, Cf, Cm]
    obj = @(x) sum((x(1)*x((1:nT)+8) + (1-x(1))*(x((1:nT)+nT+8)+x((1:nT)+2*nT+8)) - tissueTAC).^2 + ...
                   (x(2)*x((1:nT)+8) + (1-x(2))*(x((1:nT)+nT+8)+x((1:nT)+2*nT+8)) - bloodTAC).^2);
    
    if fitblood
        startingX = [0.17;0.98;0.34;0.34;1.64;1.64;.05;.05;bloodTAC;bloodTAC/6;(tissueTAC-.17*bloodTAC)-bloodTAC/6];
        lowerLim = [0;0.9;zeros(3*nT+6,1)];
        upperLim = [0.5;1]; % a<0.5 and (1-b)<0.1
        res = fmincon(obj, startingX, [], [], [], [], lowerLim, upperLim, @mycon, options);
    else
        startingX = [0.05;1.0;0.34;0.34;1.64;1.64;.05;.05;bloodTAC;bloodTAC/6;(tissueTAC-.17*bloodTAC)-bloodTAC/6];
        lowerLim = [0;1.0;zeros(3*nT+6,1)];
        upperLim = [0.1;1.0];% a<0.1 and (1-b)=0
        res = fmincon(obj, startingX, [], [], [], [], lowerLim, upperLim, @mycon, options); % Run the optimization, solve for x
    end
    
    if plotOn
        plotFit(res, nT, fitblood);
    end
    
    cp = res((1:nT)+8); % Extract blood compartment
    Ki_Norm = res(3)*res(7)/(res(5)+res(7)); % Ki = k1*k3/(k2+k3)
    Ki_Hyp = res(4)*res(8)/(res(6)+res(8)); % Ki = k1*k3/(k2+k3)
    MRglc_Norm = Ki_Norm*PG/18/LC; % MRglc = Ki*PG/LC;  Convert units of PG to SI units
    MRglc_Hyp = Ki_Hyp*PG/18/LC; % MRglc = Ki*PG/LC;  Convert units of PG to SI units
end

% Constraint function
% minimization is subject to the constraint c<=0, ceq=0
% used to enforce the constraint that the compartment concentrations obey
% the 3-compartment FDG model
function [c,ceq] = mycon(x)

    c = []; % I only use the equality constraint
    nT=(numel(x)-8)/3;
    
    % k1, k2, and k3 have one value for normoxia (mins 1-15), one for
    % hypoxia (mins 40-60), and vary linearly in between (mins 16-39)
    normW = [ones(15,1);linspace(1,0,25)';zeros(nT-40,1)];
    hypW = 1-normW;
    k1 = x(3)*normW+x(4)*hypW;
    k2 = x(5)*normW+x(6)*hypW;
    k3 = x(7)*normW+x(8)*hypW;
    
    
    % dCf/dt = k1*Cb - k2*Cf - k3*Cf
    % dCm/dt = k3*Cf
    % 
    % Therefore:
    % ceq = [dCf/dt + (k2+k3)*Cf - k1*Cb;
    %        dCm/dt - k3*Cf;
    %        Cf(1); Cm(1)]
    ceq = [x((2:nT)+nT+8)-x((1:(nT-1))+nT+8)+(k2(1:(nT-1))+k3(1:(nT-1))).*x((1:(nT-1))+nT+8) - k1(1:(nT-1)).*x((1:(nT-1))+8);...
           x((2:nT)+2*nT+8)-x((1:(nT-1))+2*nT+8)-k3(1:(nT-1)).*x((1:(nT-1))+nT+8);x(9+nT);x(9+2*nT)];
end

% Plot the fit to the kinetic model
function plotFit(res, nT, fitblood)
    figure();
    plot((1:nT),bloodTAC,'or');
    hold on;
    plot((1:nT),tissueTAC,'ok');
    plot((1:nT),res((1:nT)+8),'--m');
    plot((1:nT),res((1:nT)+nT+8),'--g');
    plot((1:nT),res((1:nT)+2*nT+8),'--c');
    plot((1:nT),res(2)*res((1:nT)+8) + (1-res(2))*(res((1:nT)+nT+8)+res((1:nT)+2*nT+8)),'-r');
    plot((1:nT),res(1)*res((1:nT)+8) + (1-res(1))*(res((1:nT)+nT+8)+res((1:nT)+2*nT+8)),'-k');
    hold off;

    if fitblood
        legend('blood','myocardium', 'cp','c1','c2','estBlood','estMyo');
    else
        legend('blood','skel', 'cp','c1','c2','estBlood','estSkel');
    end
end
