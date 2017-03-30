%zList=[25, 50, 100, 150, 200, 250];

zList=200;

Signal_C=2000;

Bg_C_RatioList=0.01:0.01:0.05;
     
d1_List=[100, 200, 400, 800];

d2_List=[200, 400, 800, 1e5];

CellThicknessList=[1000 2000 3000 4000];

CellThickness=1000;

d1=100;
%d2=1e6;

ParticleNum=1000;

Is_C_Known=0;

GridNumPerPixel_Init=11;
%%
matlabpool 4

parfor k=1:4
     
    d2=d2_List(k);
    
    %CellThickness=CellThicknessList(k)
    
    Depth=[d1; d2];

    SaveFileName=['SimData_3\DetectionResult_2Angles_' ...
                   num2str(d1) '_' num2str(d2) '_CellThickness_' num2str(CellThickness) '.mat'];
        
    [x_RMSE, y_RMSE, z_RMSE, R_MAPE, C_MAPE, I_MAPE, ...
     ParticleSNRTable]=ReconstructionInSimulation(Depth, ...
                                                  Signal_C, ...
                                                  Bg_C_RatioList,...
                                                  zList, ...                                                                                                              
                                                  CellThickness, ...
                                                  ParticleNum, ...
                                                  SaveFileName, ...
                                                  Is_C_Known,...
                                                  GridNumPerPixel_Init);
          
end

matlabpool close
%%
load('DetectionResult_2Angles_100_200_CellThickness_1000.mat')
x_RMSE_1=x_RMSE;
y_RMSE_1=y_RMSE;
z_RMSE_1=z_RMSE;
ParticleSNRTable_1=ParticleSNRTable;

load('DetectionResult_2Angles_100_400_CellThickness_1000.mat')
x_RMSE_2=x_RMSE;
y_RMSE_2=y_RMSE;
z_RMSE_2=z_RMSE;
ParticleSNRTable_2=ParticleSNRTable;

load('DetectionResult_2Angles_100_800_CellThickness_1000.mat')
x_RMSE_3=x_RMSE;
y_RMSE_3=y_RMSE;
z_RMSE_3=z_RMSE;
ParticleSNRTable_3=ParticleSNRTable;

load('DetectionResult_2Angles_100_100000_CellThickness_1000.mat')
x_RMSE_4=x_RMSE;
y_RMSE_4=y_RMSE;
z_RMSE_4=z_RMSE;
ParticleSNRTable_4=ParticleSNRTable;
%---------------------------------------------------------------------
ColorList={'r', 'b', 'g', 'm', 'c', 'y'};

for n=1
    
    figure;
    
    z=zList(n);
    
    for k=1:5
        hold on; 
        plot(1:4, [z_RMSE_1(n, k), z_RMSE_2(n,k), z_RMSE_3(n,k), z_RMSE_4(n,k)], '-o', 'Color', ColorList{k}); 
    end
    
    axis([1 4 0 100]);
    set(gca,'XTick', [1, 2, 3, 4])
    set(gca,'XTickLabel',str2mat('200', '400', '800', 'Inf'))
    temph=legend('\lambda is 0.01', '\lambda is 0.02', '\lambda is 0.03','\lambda is 0.04', '\lambda is 0.05');
    %set(temph,'FontSize', 14, 'Location', 'NorthWest')
    set(temph,'FontSize', 14, 'Location', 'NorthEast')
    xlabel('The Second Depth (nm)', 'FontSize', 14)
    ylabel('RMSE of z postions (nm)', 'FontSize', 14)
    title(['Z-Position is ' num2str(z) 'nm'], 'FontSize', 14)
    grid on
end