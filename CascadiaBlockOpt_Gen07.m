%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%    Eileen Evans    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%   Uses realistic noise levels. Onshore station noise comes from reported uncertainties in velocity field.
%   Offshore station uncertainties are prescribed. (offshore, onshore)
%   This version also calculates uncertainties on final model parameters 
% 
%   INPUT 
%       1. Station: structure with existing network of stations 
%       2. StaGrid: structure with potential observations over which to evalutate entropy
%       3. BlockDir: Directory of Block model run containing blocks of
%       interest
%       4. OptStruct: Optimization structure w/ fields:
%                       - C_B, C_T, lensc, includeChadwell, onshoreonly,
%                       smoothing, Nnew_stop, sigma0
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function HSAVE = CascadiaBlockOpt_Gen07(StaGrid,BlockDir,OptStruct)

xl = [228 250]-360;
yl = [30 51];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% Test with equal uncertainties (0.1, 20, 0.1 gives station 1 mid-plate)
C_B = OptStruct.C_B; % .01;
C_T = OptStruct.C_T; % 20; % mm/yr 
lensc = OptStruct.lensc; % 0.05;
sigma0 = OptStruct.sigma0;

includeChadwell = OptStruct.includeChadwell; % 0;
onshoreonly = OptStruct.onshoreonly; % 0;
smoothing = OptStruct.smoothing; % 1;

Nnew_stop = OptStruct.Nnew_stop; % 10; 

RemoveStations = 0; % if 1, removes a selected location from future consideration; if 0, no removal.

if includeChadwell
    chtext = '_ch';
else
    chtext = '_noch';
end

if onshoreonly
    ontext = '_on';
else
    ontext = '_all';
end

if smoothing
    smoothtext = sprintf('_sm-%2.3f',lensc);
else
    smoothtext = '_nosm';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  load in the things we need
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

outdir = BlockDir;
commandFile0 = dir(fullfile(outdir,'*.command'));
commandFile = fullfile(commandFile0.folder,commandFile0.name);
Command                                         = ReadCommand(commandFile);
Station                                         = ReadStation(Command.staFileName); % Read station file
Station                                         = ProcessStation(Station, Command);
Station.eastSig(Station.eastSig < 0.1)          = 0.1;
Station.northSig(Station.northSig < 0.1)        = 0.1;
Segment                                         = ReadSegmentTri(Command.segFileName);
Patches                                         = ReadPatches(Command.patchFileNames);
Patches                                         = PatchEndAdjust(Patches, Segment); % Adjust patch end coordinates to agree with segment end points
Patches                                         = PatchCoords(Patches); % Create patch coordinate arrays
Segment                                         = ProcessSegment(Segment, Command);
[Patches, Command]                              = ProcessPatches(Patches, Command, Segment);
Block                                           = ReadBlock(Command.blockFileName); % Read block file
Block.aprioriTog                                = zeros(size(Block.aprioriTog)); % Remove JDF prior
[Segment, Block, Station]                       = BlockLabel(Segment, Block, Station);
nblocks                                         = numel(Block.interiorLon);
npatches                                        = numel(Patches.lon1);

%%% Get some JDF stuff
JDF = find(ismember(Block.name(:,1:3),'JDF','rows'));
%%% Convert to UTM to make distances easier
E = wgs84Ellipsoid;
zone                                                = '11S';
[ellipsoid,estr]                                    = utmgeoid(zone);
utmstruct                                           = defaultm('utm'); 
utmstruct.zone                                      = zone; 
utmstruct.geoid                                     = E; 
utmstruct.flatlimit                                 = []; 
utmstruct.maplatlimit                               = []; 
utmstruct                                           = defaultm(utmstruct);
JDF = find(ismember(Block.name(:,1:3),'JDF','rows'));
JDFlon = Block.orderLon{JDF};
JDFlat = Block.orderLat{JDF};
[JDFx,JDFy]                                         = mfwdtran(utmstruct, JDFlat, JDFlon);

%%% Find trench-bounding blocks
trenchsegs = find(Segment.patchTog);
trenchblocks = unique([Segment.eastLabel(trenchsegs); Segment.westLabel(trenchsegs)]);
colkeep = sort([trenchblocks*3-2; trenchblocks*3-1; trenchblocks*3]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  load existing partials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

A                                               = load([outdir filesep 'kernels.mat'], 'elastic');
B                                               = load([outdir filesep 'kernels.mat'], 'tri');

Partials.elastic                                = A.elastic;
Partials.tri                                    = B.tri(:,2:3:end); % dip slip partials only
Partials.rotation                               = GetRotationPartials(Segment, Station, Command, Block);
Partials.rotation = Partials.rotation(:,colkeep);
Partials.slip                                   = GetSlipPartials(Segment, Block);
% G = [Partials.rotation - Partials.elastic * Partials.slip, -Partials.tri];
G = [Partials.rotation, -Partials.tri];
szrot                                            = size(Partials.rotation);
rowkeep                                          = setdiff(1:szrot(1), [3:3:szrot(1)]);
G = G(rowkeep,:); % remove up partials;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Offshore grid, eventually quite dense
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% get IN/ON

IN = StaGrid.IN;

% keyboard
[StaGrid.x, StaGrid.y, StaGrid.z]                  = sph2cart(DegToRad(StaGrid.lon), DegToRad(StaGrid.lat), 6371);

[Segment, Block, StaGrid]                        = BlockLabel(Segment, Block, StaGrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Load Chadwell Stations and get partials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% load('./station/Chadwell_01');
load('./station/Offshore');
Chadwell = Offshore;
[Chadwell.x, Chadwell.y, Chadwell.z]                  = sph2cart(DegToRad(Chadwell.lon), DegToRad(Chadwell.lat), 6371);
[Segment, Block, Chadwell]                        = BlockLabel(Segment, Block, Chadwell);

ChPartials.rotation                                = GetRotationPartials(Segment, Chadwell, Command, Block);
ChPartials.rotation                                = ChPartials.rotation(:,colkeep);
ChPartials.elastic                                 = GetElasticPartials(Segment, Chadwell);
[ChPartials.tri, ~, NewPatches]                    = GetTriCombinedPartials(Patches, Chadwell, [1 0]);
% NewPartials.slip                                    = GetSlipPartials(Segment, Block);
ChPartials.tri                                     = ChPartials.tri(:,2:3:end);
% ChG = [ChPartials.rotation - ChPartials.elastic * Partials.slip, -ChPartials.tri];
ChG = [ChPartials.rotation, -ChPartials.tri];
szrot                                            = size(ChPartials.rotation);
rowkeep                                          = setdiff(1:szrot(1), [3:3:szrot(1)]);
ChG = ChG(rowkeep,:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%  Load Canada Stations and get partials
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% load('./station/Heesemann.mat');
% [Heesemann.x, Heesemann.y, Heesemann.z]                  = sph2cart(DegToRad(Heesemann.lon), DegToRad(Heesemann.lat), 6371);
% [Segment, Block, Chadwell]                        = BlockLabel(Segment, Block, Heesemann);
% 
% HPartials.rotation                                = GetRotationPartials(Segment, Heesemann, Command, Block);
% HPartials.rotation                                = HPartials.rotation(:,colkeep);
% HPartials.elastic                                 = GetElasticPartials(Segment, Heesemann);
% [HPartials.tri, ~, NewPatches]                    = GetTriCombinedPartials(Patches, Heesemann, [1 0]);
% % NewPartials.slip                                    = GetSlipPartials(Segment, Block);
% HPartials.tri                                     = HPartials.tri(:,2:3:end);
% % ChG = [ChPartials.rotation - ChPartials.elastic * Partials.slip, -ChPartials.tri];
% HG = [HPartials.rotation, -HPartials.tri];
% szrot                                            = size(HPartials.rotation);
% rowkeep                                          = setdiff(1:szrot(1), [3:3:szrot(1)]);
% HG = HG(rowkeep,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Deal with onshore-offshore differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% sigma0 = 0.003; % sigma for OFFSHORE prediction error, meters
% sigma0 = 2.1; % sigma for OFFSHORE prediction error, mm/yr

% OnshoreEastSig = mean(Station.eastSig*1e-3)*ones(size(StaGrid.lon(IN))); % eastSig in meters
% OnshoreNorthSig = mean(Station.northSig*1e-3)*ones(size(StaGrid.lon(IN))); % eastSig in meters
OnshoreEastSig = mean(Station.eastSig)*ones(size(StaGrid.lon(IN))); % eastSig in mm/yr
OnshoreNorthSig = mean(Station.northSig)*ones(size(StaGrid.lon(IN))); % eastSig in mm/yr

%OnshoreEastSig = 0.001*ones(size(StaGrid.lon(IN))); % eastSig in meters
%OnshoreNorthSig = 0.001*ones(size(StaGrid.lon(IN))); % eastSig in meters

[Chadwell.eastSig, Chadwell.northSig]               = deal(sigma0*ones(size(Chadwell.lon))); % meters
[StaGrid.eastSig, StaGrid.northSig]                 = deal(sigma0*ones(size(StaGrid.lon))); % meters
StaGrid.eastSig(IN)                                 = OnshoreEastSig; % 
StaGrid.northSig(IN)                                = OnshoreNorthSig; % 

% figure; scatter(StaGrid.lon,StaGrid.lat,10,StaGrid.eastSig,'filled')
% keyboard
if onshoreonly
    names = fieldnames(StaGrid);
for nn = 1:numel(names);
    StaGrid.(names{nn}) = StaGrid.(names{nn})(IN,:);
end 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Get New Partials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

NewPartials.rotation                                = GetRotationPartials(Segment, StaGrid, Command, Block);
NewPartials.rotation                                = NewPartials.rotation(:,colkeep);
NewPartials.elastic                                 = GetElasticPartials(Segment, StaGrid);
[NewPartials.tri, ~, NewPatches]                    = GetTriCombinedPartials(Patches, StaGrid, [1 0]);
% NewPartials.slip                                    = GetSlipPartials(Segment, Block);
NewPartials.tri                                     = NewPartials.tri(:,2:3:end);
% NewG = [NewPartials.rotation - NewPartials.elastic * Partials.slip, -NewPartials.tri];
NewG = [NewPartials.rotation, -NewPartials.tri];
szrot                                               = size(NewPartials.rotation);
rowkeep                                             = setdiff(1:szrot(1), [3:3:szrot(1)]);
NewG = NewG(rowkeep,:);

% G2 = [G; NewG];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Stack Partials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


if includeChadwell == 1;
    G2 = [G; ChG; NewG];
else
    G2 = [G; NewG];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Foward model (for Nparam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Nparam = 3*numel(Block.interiorLon) + 1*numel(Patches.lon1);
Nparam = 3*numel(trenchblocks) + 1*numel(Patches.lon1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Setup Inverse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

w1 = zeros(2*numel(Station.lon),1);
w1(1:2:end) = 1./(Station.eastSig).^2; % check units
w1(2:2:end) = 1./(Station.northSig).^2; % 

w2 = zeros(2*numel(Chadwell.lon),1);
w2(1:2:end) = 1./(Chadwell.eastSig).^2; % check units
w2(2:2:end) = 1./(Chadwell.northSig).^2; % 

if includeChadwell == 1;
    w = [w1; w2];
    Nst = numel(StaGrid.lon) + numel(Chadwell.lon) + numel(Station.lon);
    iinew = (1:(numel(Station.lon) + numel(Chadwell.lon)))';
else
    w = w1;
    Nst = numel(StaGrid.lon) + numel(Station.lon);
    iinew = (1:(numel(Station.lon)))';
end

woff = zeros(2*numel(StaGrid.lon),1);
woff(1:2:end) = 1./(StaGrid.eastSig).^2; % check units
woff(2:2:end) = 1./(StaGrid.northSig).^2; % 

wall = [w; woff];
W = diag(wall);

mtarget = zeros(Nparam,1);
Na = length(mtarget(:));
H_const = 0.5*Na*[log(2*pi) + 1];


share                                         = SideShare(Patches.v);
dists                                         = TriDistCalc(share,Patches.xc,Patches.yc,Patches.zc);

L = (lensc * MakeTriLaplacianDist(share,dists));

Gon = G; % onshore
if includeChadwell;
    Gex = [G; ChG]; % existing
else
    Gex = G;
end
CB_most = 0.01;


CB_vec = CB_most * ones(3*nblocks,1);
jdfix = find(strcmp(cellstr(Block.name),'JDF'));
CB_vec((3*jdfix - 2):3*jdfix) = C_B;
CB_vec = CB_vec(colkeep);

Cm_inv_block = diag( 1./CB_vec.^2 );

Cm_inv_tri = 1/C_T^2 * eye(npatches);

if smoothing
    Cm_inv_S = blkdiag(Cm_inv_block, L'*Cm_inv_tri*L);
else
    Cm_inv_S = blkdiag(Cm_inv_block, Cm_inv_tri);
end

Cdinv = spdiags(w1,0,numel(w1),numel(w1));
Sigma_inv = Cm_inv_S + (Gon'*Cdinv*Gon); % Sigma includes onshore stations only (this is a good constant metric to use)
Vstar = inv(Sigma_inv);

H_on = 0.5*logdet(2*pi*exp(1)*Vstar,'chol');

if ~exist('HSAVE','var')
    
    HSAVE=struct;
    
    Cdinv = spdiags(w,0,numel(w),numel(w));
    Sigma_inv = Cm_inv_S + (Gex'*Cdinv*Gex); % Sigma includes onshore stations only (this is a good constant metric to use)
    Vstar = inv(Sigma_inv);

    Hnew = 0.5*logdet(2*pi*exp(1)*Vstar,'chol');
    
    HSAVE(1).H = Hnew;
    HSAVE(1).Hbest = Hnew;
    HSAVE(1).ii = iinew;
%     HSAVE(1).uvec = diag(Sigma\eye(size(Gex, 2)));
    
    
    Nnew_start = 1;
else
    
    Nnew_start = numel(HSAVE);
end

iiold=HSAVE(1).ii;
if RemoveStations == 0
iitest=setdiff(1:Nst,iiold); % if this line is OUTSIDE the loop, we test over the same stations every time (without removal)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%  Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tic
for inew = Nnew_start:Nnew_stop
    disp(['inew = ' num2str(inew)])
    tic
    
    iiold=HSAVE(inew).ii; % these now exist!
    if RemoveStations == 1
        iitest=setdiff(1:Nst,iiold); % these don't yet! % if this line is INSIDE the loop, we test over the same stations every time (with removal)
    end
    myH=nan*ones(length(iitest),1);
    
    dummyH = nan*ones(size(iitest));
    indefkeep = nan*ones(size(iitest));
    
%     keyboard
    parfor dummy = 1:length(iitest)
%     for dummy = 1:length(iitest)
%         tic
        itest=iitest(dummy);
        iinew=[iiold; itest];
        
        Gnew = [G2(2*iinew-1,:); G2(2*iinew,:)]; % horizonal only, so don't need to be cute  with indexing here.
%         Gnew = G2(iinew,:);
%         Gnew = [Gold; Gtest]; 
        wall_new = [wall(2*iinew-1); wall(2*iinew)];
        Cdinv_new = spdiags(wall_new,0,numel(wall_new),numel(wall_new));
%         Sigma = Cm_inv + [Gnew'*Gnew]/sigma0^2; %% thisbecomes 
        Sigma_inv = Cm_inv_S + (Gnew'*Cdinv_new*Gnew);
        Vstar = inv(Sigma_inv);
        
        %%% check postitive definite
        [R,p] = chol(Vstar);
        if p == 0;
            dummyH(dummy) = 0.5*logdet(2*pi*exp(1)*Vstar,'chol');
        else
%             keyboard
            dummyH(dummy) = 0.5*logdet(2*pi*exp(1)*Vstar,'chol');
            indefkeep(dummy) = 1;
        end
%         toc
    end
    
    RdummyH = (round(dummyH.*1e3))./1e3;
    myH(iitest) = RdummyH;

%     %%% To check
%     ON = StaGrid.blockLabel == JDF;
%     figure; scatter(StaGrid.lon(ON),StaGrid.lat(ON),10,dummyH(ON)','filled')
%     figure; scatter(StaGrid.lon(ON),StaGrid.lat(ON),10,RdummyH(ON)','filled')
% 
    
%     ibest = find(myH==min(myH),1);
    ibest = find(RdummyH==min(RdummyH)) + (min(iitest)-1);
%     PlotFirstSTAs(HSAVE,StaGrid,CB,CT,lensc,foldername)
% keyboard
    %%% if many minima on JDF plate
%     if numel(ibest) > 1 && StaGrid.blockLabel(ibest(1) - (min(iitest)-1)) == JDF
%%% If minimum (even just 1) on JDF plate
    if StaGrid.blockLabel(ibest(1) - (min(iitest)-1)) == JDF && inew == 1
       G = polygeom(JDFx,JDFy);
       Cx = G(2);
       Cy = G(3);
%        [Clat, Clon] = minvtran(UTMstruct,G(2),G(3));
%        Clon(Clon>0) = Clon + 360;
       
       D = abs((StaGrid.UTMx - Cx).^2 + (StaGrid.UTMy - Cy).^2);
       [md, id] = min(D);
       
       ibest = id + (min(iitest)-1);
%     elseif numel(ibest) > 1 && StaGrid.blockLabel(ibest(1) - (min(iitest)-1)) ~= JDF
    else
%         keyboard
        ibest = ibest(1);
    end
    
%     ibest = ibest + (min(iitest)-1);
    
    %%% Calculate mutual information thing
    INEW = [iiold; ibest];
    Gbest = [G2(2*INEW-1,:); G2(2*INEW,:)];
    wall_new = [wall(2*INEW-1); wall(2*INEW)];
    Cdinv_new = spdiags(wall_new,0,numel(wall_new),numel(wall_new));
    Sigma_best_inv = Cm_inv_S + (Gbest'*Cdinv_new*Gbest);
    Vstar_best = inv(Sigma_best_inv);
    
    
%     MI_best = zeros(Na,1);
%     for pp = 1:Na
%         Si = Sigma_best(pp,pp);
%         Sni_ids = setdiff(1:Na,pp);
%         Sni = Sigma_best(Sni_ids,Sni_ids);
%         MI_best(pp) = log(2) * (1/2 * ( log(Si) + (logdet(Sni) - logdet(Sigma_best)) ));
% 
%     end
    
    HSAVE(end+1).H=myH;
    HSAVE(end).iitest=iitest;
    HSAVE(end).Hbest=myH(ibest);
    HSAVE(end).ii=[iiold;ibest]; % now a new station exists!
    HSAVE(end).ibest=ibest;
%     HSAVE(end).minfo = MI_best;
%     HSAVE(end).uvec = diag(Sigma_best\eye(size(Gbest, 2)));
    
    save('HSAVE_temp.mat','HSAVE','StaGrid','Command');
    
    toc
%      keyboard
end

% matfilename = sprintf('HSAVE_%dx%d%s%s%s_%d_sigma0-%d_rem-%d.mat',xdim,ydim,chtext,ontext,smoothtext,Nnew_stop,sigma0,RemoveStations);
% save(matfilename,'HSAVE','StaGrid','Cm_inv','Command');
toc 
% PlotOptAll_color(HSAVE)
