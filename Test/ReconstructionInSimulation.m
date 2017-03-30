function  [x_RMSE, y_RMSE, z_RMSE, R_MAPE, C_MAPE, I_MAPE, ...
           ParticleSNRTable]=ReconstructionInSimulation(Depth, ...
                                                        Signal_C, ...
                                                        Bg_C_RatioList,...
                                                        zList, ...                                                                                                                    
                                                        CellThickness, ...
                                                        ParticleNum, ...
                                                        SaveFileName, ...
                                                        Is_C_Known,...
                                                        GridNumPerPixel_Init)
%% TIRF Paramters
AngleNum=length(Depth);
ni=1.522;
nt=1.333;
Lambda=488;
I0=Cal_d_to_I0(Depth, nt, ni, Lambda);
I0=I0/I0(1);

SysGain=50;

PSFSigmaXY=200;

%PixelSize=180;
PixelSize=160;

MATIRF_Param.PSFSigmaXY=PSFSigmaXY;
MATIRF_Param.AngleNum=AngleNum;
MATIRF_Param.Depth=Depth;
MATIRF_Param.ni=ni;
MATIRF_Param.nt=nt;
MATIRF_Param.Lambda=Lambda;
MATIRF_Param.I0=I0;
MATIRF_Param.SysGain=SysGain;
MATIRF_Param.PixelSize=PixelSize;
%-----------------------------------------------------------------------------------------------------------
Lx=(ceil(sqrt(ParticleNum))+1)*20;
Ly=Lx;
% Set particle positions
[x y]=meshgrid(20:20:(Lx-20), 20:20:(Ly-20));

x=x(:); x=x(1:ParticleNum);
y=y(:); y=y(1:ParticleNum);

RadiusRange=[25, 50];
if ParticleNum > 1 && RadiusRange(2) > RadiusRange(1)
    R=(RadiusRange(1):(RadiusRange(2)-RadiusRange(1))/(ParticleNum-1):RadiusRange(2))';
else
    R=RadiusRange(1)*ones(ParticleNum, 1);
end

CellLabel=ones(Ly, Lx);

RadiusRange_Pixel=[1, 1.5]*PSFSigmaXY/PixelSize;
%---------------------------------------------------------------------------------------------------------
Range_C=[0.5*Signal_C*ones(1, ParticleNum)
         2*Signal_C*ones(1, ParticleNum)];

Range_z=[0.5*zList(1)*ones(1, ParticleNum)
         2*zList(end)*ones(1, ParticleNum)];
%---------------------------------------------------------------------------------------------------------
BgDensityList=Bg_C_RatioList*Signal_C;

ParticleSNRTable=zeros(length(zList), length(BgDensityList));
for n=1:length(zList)
for k=1:length(BgDensityList);          
    temp_z=zList(n)*ones(ParticleNum, 1);
    ParticleSNRTable(n,k)=MATIRF_Cal_ParticleSNR(Signal_C, temp_z, R, BgDensityList(k), CellThickness, MATIRF_Param);
end
end

MethodStr='x_y_A_R_b';

x_RMSE=[]; y_RMSE=[]; z_RMSE=[]; R_MAPE=[]; C_MAPE=[]; I_MAPE=[];

ImageSize=[Ly, Lx];
%%
counter_z=0;

for z_position=zList

counter_z=counter_z+1;
counter_b=0;
    
z=z_position*ones(ParticleNum, 1);
C=Signal_C*ones(ParticleNum, 1);

ParticleFeature=zeros(5, ParticleNum);
ParticleFeature(1,:)=x*PixelSize;
ParticleFeature(2,:)=y*PixelSize;
ParticleFeature(3,:)=z;
ParticleFeature(4,:)=R;
ParticleFeature(5,:)=C;

tic
disp('ParticleSimulation')
[CleanImageStack, TotalIntensity]=MATIRF_ParticleSimulation_numerical(MATIRF_Param, ParticleFeature, ImageSize, GridNumPerPixel_Init);
toc
       
for Bg_C_Ratio=Bg_C_RatioList  

counter_b=counter_b+1;
    
disp(['z_position is ' num2str(z_position) ' Bg_C_Ratio is ' num2str(Bg_C_Ratio)]);    
%% =====================================================================================================================

BgDensity=BgDensityList(counter_b);

%% Generate MATIRF Images
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);

BgList=MATIRF_BackgroundSimulation(MATIRF_Param, CellThickness, BgDensity);

ImageStack=MATIRF_NoiseSimulation(CleanImageStack, BgList, MATIRF_Param.SysGain);
%% Estimation 
Ib_Array=zeros(Ly, Lx, AngleNum);
Ib_SigmaArray=zeros(Ly, Lx, AngleNum);
for m=1:AngleNum
          
    Im=ImageStack(:,:,m);
       
    Im_median=median(Im(:));
         
    Ib_Array(:,:,m)=Im_median;         
    Ib_SigmaArray(:,:,m)=sqrt(Im_median);
end

x_Init=x+0.5;
y_Init=y+0.5;

ParticleNum=length(x_Init);
R_Init=ones(ParticleNum, 1);
%% =============================================================================================================
x_Data=[]; y_Data=[]; R_Data=[]; x_hat=[]; y_hat=[]; z_hat=[]; R_hat=[]; C_hat=[]; A_hat=[]; Astd=[]; b_hat=[];
%% ------------------------------------------------------------------------------------------
disp('Estimation - Get A R b x y')
tic
for m=1:AngleNum

Im=ImageStack(:,:,m); 
    
[Am_Init bm_Init Radius_xy Range_Am Range_R Range_bm]=GaussMixFit_InitialParam(x_Init, y_Init, RadiusRange_Pixel, 1, Im, CellLabel, Ib_Array(:,:,m), Ib_SigmaArray(:,:,m));
                                          
[xm ym Am Rm bm A_std_m]=GaussMixFit(MethodStr, x_Init, y_Init, Am_Init, R_Init, bm_Init, ...
                                     Radius_xy, Range_Am, Range_R, Range_bm, ...                                                    
                                     [], Im, CellLabel, Ib_Array(:,:,m), Ib_SigmaArray(:,:,m), [], []);

A_hat(:,m)=Am;
Astd(:,m)=A_std_m;
x_Data(:,m)=xm;
y_Data(:,m)=ym;
R_Data(:,m)=Rm;
b_hat(:,m)=bm;
end
toc
%% 
disp('Estimation - Get c z')
tic
for k=1:ParticleNum
    IdxList_1=find(~isnan(Astd(k, :)));
    IdxList_2=find(isnan(Astd(k, :)));
    if ~isempty(IdxList_2)
        if ~isempty(IdxList_1)
            Astd(k, IdxList_2)=max(Astd(k, IdxList_1));
        else
            Astd(k, IdxList_2)=max(Astd(:));
        end
    end
end

xy_Weight=Cal_xy_Weight(x_Data, y_Data, A_hat, Astd, Ib_Array, Ib_SigmaArray);
x_hat=sum(x_Data.*xy_Weight, 2);
y_hat=sum(y_Data.*xy_Weight, 2);

R_Weight=Cal_R_Weight(x_hat, y_hat, A_hat, Astd, Ib_Array, Ib_SigmaArray);
R_hat=sum(R_Data.*R_Weight, 2);

%A_Weight=ones(ParticleNum, AngleNum); 
A_Weight=Cal_A_Weight(x_hat, y_hat, A_hat, Astd, Ib_Array, Ib_SigmaArray);

if Is_C_Known == 0

    [C_hat, z_hat, A_Model, A_ModelSigma] = C_Z_Estimation(Range_C, Range_z, A_hat, Astd, A_Weight, MATIRF_Param, 1);
    
elseif Is_C_Known == 1
    %------ known C -------------------------------------------------------

    z_hat=Z_Estimation(C, Range_z, A_hat, A_Weight, MATIRF_Param);
    C_hat=C;
else
    disp('Wrong Is_C_Known')
end

toc
%% =====================================================================================================================
x_RMSE(counter_z, counter_b)=sqrt(mean((x-x_hat).^2));
y_RMSE(counter_z, counter_b)=sqrt(mean((y-y_hat).^2));
z_RMSE(counter_z, counter_b)=sqrt(mean((z-z_hat).^2));

tempR=sqrt(R.^2+PSFSigmaXY^2);
tempR=tempR/PixelSize;
R_MAPE(counter_z, counter_b)=mean(abs((tempR-R_hat)./tempR));

C_MAPE(counter_z, counter_b)=mean(abs((Signal_C-C_hat)./Signal_C));

TotalIntensity_hat=bsxfun(@times, A_hat', (R_hat.^2)');
TotalIntensity_hat=2*pi*TotalIntensity_hat';

I_MAPE(counter_z, counter_b)=mean(abs((TotalIntensity(:)-TotalIntensity_hat(:))./TotalIntensity(:)));
%% =====================================================================================================================
end
end
%%
if ~ isempty(SaveFileName)

    save(SaveFileName, ...
         'MATIRF_Param', 'Bg_C_RatioList', 'Signal_C', 'zList', 'ParticleSNRTable', 'CellThickness', 'ParticleNum', ...
         'x_RMSE', 'y_RMSE', 'z_RMSE', 'R_MAPE', 'C_MAPE', 'I_MAPE');
end