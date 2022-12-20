%% Create surfaces here, optionally based on the used defined parameters in
% input.txt, found in the struct "param"
% -------------------------------------------------------------------------

% height = 0.1875 * param.size;
% zres = 0.5* 0.1875 * param.res;

% surf_s = Surface(128, 2500);
% surf_m = Surface(256, 5000);
% surf_l = Surface(384, 7500);
% surf_xl = Surface(512, 10000);

% base_s = zeros(128);%[base_m, base_m; base_m, base_m];
% base_s(:,1:24)=4;
% base_s(:,25:48)=3;
% base_s(:,49:72)=2;
% base_s(:,73:96)=1;

% base_m = [base_s, base_s; base_s, base_s];
% base_l = [base_s, base_s, base_s; base_s, base_s, base_s; base_s, base_s, base_s];
% base_xl = [base_m, base_m; base_m, base_m];

% surf_r = Surface(128, 2500);
% surf_r.addRoughsurf('mode','Rough','height',100)
% surf_r.placeFeatures("PBC", true, "mode", "add");

% base_s = surf_r.surfMatrix;
% base_m = [base_s, base_s; base_s, base_s];
% base_l = [base_s, base_s, base_s; base_s, base_s, base_s; base_s, base_s, base_s];
% base_xl = [base_m, base_m; base_m, base_m];

% surf_silver = Surface(param.res, param.size);
surf_random = Surface(param.res, param.size);
%surf_algaas = Surface(param.res, param.size);
%surf_oxide = Surface(param.res, param.size);
% surf_oxide_uniform = Surface(param.res, param.size);

%B=readmatrix("side_projects/Gaussian_swaps3_1.txt");
%B=readmatrix("src/side_projects/generated_surfaces/2Doptim51K3_"+param.simit+".txt");
%B=readmatrix("side_projects/surf_no_zero.txt");
%B=reshape(B,size(B,1),size(B,1),[]);


% B=[0,0,0,0,1,1,1,1,1,1;
%     0,0,0,0,1,1,1,1,1,1;
%     0,0,0,0,0,0,1,1,0,0;
%     0,0,0,0,0,0,1,1,0,0;
%     0,0,0,0,1,1,1,1,1,1;
%     0,0,0,0,1,1,1,1,1,1;
%     1,1,1,1,1,1,0,0,0,0;
%     1,1,1,1,1,1,0,0,0,0;
%     0,0,0,0,1,1,0,0,1,1;
%     0,0,0,0,1,1,0,0,1,1
% ];

%QRA1+2-1 (res 3)

% B= [1,0,0;
%     1,1,1;
%     1,0,0;
% ];

%QRA5+8+9+10+13+16

% 1 (res 9)
% B= [0,0,1,1,1,1,0,0,0;
%     1,1,0,0,0,1,0,1,0;
%     1,1,0,0,1,0,0,1,0;
%     1,0,1,1,1,0,0,1,1;
%     1,0,0,1,0,0,0,0,1;
%     1,1,1,1,0,1,1,0,0;
%     1,1,0,0,0,0,1,0,0;
%     1,0,0,0,0,0,1,1,0;
%     1,1,1,1,0,1,1,1,0;
% ];

% 2 (res 9)
% B= [1,0,0,1,1,0,1,0,1;
%     1,1,0,1,1,0,0,0,0;
%     0,0,1,1,0,0,1,1,1;
%     0,0,0,0,1,0,1,0,0;
%     1,1,0,1,1,1,1,0,0;
%     0,1,1,0,0,0,1,1,0;
%     0,1,1,0,0,0,0,0,1;
%     0,1,1,1,1,1,0,1,1;
%     1,0,0,0,1,0,1,1,0;
% ];

% B=repmat(B,9,9);


%B=ones(64,64,5);
%B(:,:,1)=0;
%B=zeros(64,64,2);
%B(:,:,1:2)=1;
%B=ones(256,256);
%%


% surf_s.addFeature(Feature(base_s,2500),1,1);
% surf_m.addFeature(Feature(base_m,5000),1,1);
% surf_l.addFeature(Feature(base_l,7500),1,1);
% surf_xl.addFeature(Feature(base_xl,10000),1,1);
%%

% surf_silver.addFeature(Feature(param.res,param.size,height,"WedgeX"),1,1);
surf_random.addRoughsurf('mode','Rough','height',145)
% surf_random2.addRoughsurf('mode','Rough','height',100)
%surf_algaas.addFeature(Feature("data_maarten/AlGaAs_surface1_10um.csv",10000),1,1);
%surf_oxide.addFeature(Feature("data_maarten/Oxide_surface1_10um.csv",10000),1,1);
% surf_oxide_uniform.addUniform(70);


% surf_s.placeFeatures("PBC", true, "mode", "add");
% surf_m.placeFeatures("PBC", true, "mode", "add");
% surf_l.placeFeatures("PBC", true, "mode", "add");
% surf_xl.placeFeatures("PBC", true, "mode", "add");

% surf_silver.placeFeatures("PBC", true, "mode", "add");
surf_random.placeFeatures("PBC", true, "mode", "add");
%surf_algaas.placeFeatures("PBC", true, "mode", "add");
%surf_oxide.placeFeatures("PBC", true, "mode", "add");
% surf_oxide_uniform.placeFeatures("PBC", true, "mode", "add");
% surf_oxide.plot

% surf_r.plot
% axis off
%% Add surfaces here, index can be used in the "input" field in layers.xlsx
% Example: Add two surfaces to the list, which can be accessed in layers.xlsx
% by entering 1 or 2 for surf_oxide and surf_algaas respectively.
% surfaces{1} = surf_oxide
% surfaces{2} = surf_algaas
% -------------------------------------------------------------------------

surfaces{1} = surf_random;


%% Modify "layer" struct here
% Example: Modify thickness of second layer based on the parameter "oxidethickness".
% layer(2).L = param.oxidethickness;
% -------------------------------------------------------------------------
