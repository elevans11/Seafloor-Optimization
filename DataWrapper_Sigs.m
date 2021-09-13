%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%    Eileen Evans    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%   Wrapper script for changing data uncertainties. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

sigma0vec = logspace(0,1,20);

foldername = './DataSig_test';
if ~isdir(foldername)
    mkdir(foldername);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Get Model Stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

load('StaGrid_92x64.mat');
BlockDir = './0000000017';
OptStruct.C_B = 0.01;
OptStruct.C_T = 20;
OptStruct.lensc = 0.05;
OptStruct.includeChadwell = 0;
OptStruct.onshoreonly = 0;
OptStruct.smoothing = 0;
OptStruct.Nnew_stop = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for ss = 1:numel(sigma0vec)
    OptStruct.sigma0 = sigma0vec(ss);
    HSAVE = CascadiaBlockOpt_Gen07(StaGrid,BlockDir,OptStruct);
%     HSAVE = CascadiaBlockOpt_Gen07_Moment(StaGrid,BlockDir,OptStruct);
    matfilename = sprintf('./DataSig_test/HSAVE_sigma0-%d.mat',OptStruct.sigma0);
    save(matfilename,'HSAVE','StaGrid','BlockDir','OptStruct');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




