clear all; close all; %#ok<CLALL>
xlsFileName = 'C:\Users\pcorrado\Box\NLP_MET Data\Corrado, Phil\All_Subjects_Corrado_Patlak.xlsx';

for row=6:51 % loop through subjects
    fprintf('Subject #%i of %i.\n', row-5, 46);
    [~,txt,~] = xlsread(xlsFileName,2,sprintf('B%i:B%i',row,row)); % Get subject ID #
    csvFile = sprintf('C:\\Users\\pcorrado\\Box\\NLP_MET Data\\Corrado, Phil\\Corrado_%s_all.csv',txt{1}); % Get TAC filename
    if exist(csvFile,'file') % if TAC file exists
        tacdata = csvread(csvFile); % read TAC file
    else % Recon was split into 2 parts, join them together
        csvFile1 = strrep(csvFile,'all.csv','P1_all.csv');
        csvFile2 = strrep(csvFile,'all.csv','P2_all.csv');
        tacdata = [csvread(csvFile1);csvread(csvFile2)];
    end
    PG = xlsread(xlsFileName,2,sprintf('J%i:J%i',row,row)); % Read plasma glucose
    
    
    time = 0:(length(tacdata)-1);
    nT = min(numel(time),60); % Clip data after 60 minutes
    blood = tacdata(1:nT,2); % Extract blood TAC
    myocardium = tacdata(1:nT,1); % Extract myocardium TAC
    time=time(1:nT);
    skelmuscle = tacdata(1:nT,3); % Extract skeletal muscle TAC
    
    [Ki_myo_Norm, Ki_myo_Hyp, MMRglc_Norm, MMRglc_Hyp, cp] = myoFDGFit(blood,myocardium,PG); % Fit myocardium FDG model
    [Ki_skel_Norm, Ki_skel_Hyp, SkMRglc_Norm, SkMRglc_Hyp] = skelFDGFit(cp,skelmuscle,PG); % Fit skeletal muscle FDG model

    % Write result into spreadsheet
%     xlswrite(xlsFileName,[Ki_myo_Norm,Ki_myo_Hyp,MMRglc_Norm,MMRglc_Hyp,Ki_skel_Norm,Ki_skel_Hyp,SkMRglc_Norm,SkMRglc_Hyp],...
%              2,sprintf('T%i:AE%i',row,row));
end
