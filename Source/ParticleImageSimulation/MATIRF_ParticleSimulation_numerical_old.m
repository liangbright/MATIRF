function ImageStack=MATIRF_ParticleSimulation_numerical_old(MATIRF_Param, ParticleFeature, ImageSize, VoxelSize, ...
                                                        CellThickness, BgDensity, Scale)

PSF_Param=MATIRF_Param.PSF_Param;
PSFSigmaX=PSF_Param.PSFSigmaX;
%PSFSigmaY=PSF_Param.PSFSigmaY;
PSFSigmaZ=PSF_Param.PSFSigmaZ;

AngleNum=MATIRF_Param.AngleNum;

Depth=MATIRF_Param.Depth;

I0=MATIRF_Param.I0;

SysGain=MATIRF_Param.SysGain;

Ly=ImageSize(1);
Lx=ImageSize(2);

ImageStack=zeros(Ly, Lx, AngleNum);
%------------------------------------------------------------------------------------------
x=ParticleFeature(1,:);
y=ParticleFeature(2,:);
z=ParticleFeature(3,:);
R=ParticleFeature(4,:);
C=ParticleFeature(5,:);

ParticleNum=length(x);
%------------------------------------------------------------------------------------------
% Scale : Control Simulation Precision
%----------------------------------------------------------------------------------------------------
Wx=ceil(4*PSFSigmaX*Scale/VoxelSize(2));
%Wy_sim=ceil(3*PSFSigmaY*Scale/VoxelSize(1));
Wy=Wx;

[tempY, tempX]=ndgrid(-Wy:Wy, -Wx:Wx);

temp_PSFSigmaX=PSFSigmaX*Scale/VoxelSize(2);
%temp_PSFSigmaY=PSFSigmaY*Scale/VoxelSize(1);
temp_PSFSigmaY=temp_PSFSigmaX;

temp_xy=(tempX/temp_PSFSigmaX).^2+(tempY/temp_PSFSigmaY).^2;
PSFMask_xy=exp(-0.5*temp_xy);
%----------------------------------------------------------------------------------------------------
SumFunction=@(block_struct) sum(block_struct.data(:));
%------------------------------------------------------------------------------------------
for k=1:ParticleNum
    
    delta_x=max((R(k)+5*PSFSigmaX)/VoxelSize(2));
    %delta_y=max((R+4*PSFSigmaY)/VoxelSize(1));
    delta_y=delta_x;
    delta_z=max((3*R(k)+1)/VoxelSize(3));

    x1=round(x(k)/VoxelSize(2)-delta_x); x1=max(x1, 1);
    x2=round(x(k)/VoxelSize(2)+delta_x); x2=min(x2, Lx);
    
    y1=round(y(k)/VoxelSize(1)-delta_y); y1=max(y1, 1);
    y2=round(y(k)/VoxelSize(1)+delta_y); y2=min(y2, Ly);
    
    z1=round(z(k)/VoxelSize(3)-delta_z); z1=max(z1, 1);
    z2=round(z(k)/VoxelSize(3)+delta_z);
    
    tempLx=(x2-x1+1)*Scale;
    tempLy=(y2-y1+1)*Scale;
    Lz=(z2-z1+1)*Scale;
    
    [Y, X, Z]=ndgrid((y1-1)*Scale+1:y2*Scale,(x1-1)*Scale+1:x2*Scale, (z1-1)*Scale+1:z2*Scale);
    
    temp_xk=(x(k)/VoxelSize(2)-0.5)*Scale;
    temp_yk=(y(k)/VoxelSize(1)-0.5)*Scale;
    temp_zk=(z(k)/VoxelSize(3)-0.5)*Scale;
    
    temp_Rx=R(k)*Scale/VoxelSize(2);
    temp_Ry=R(k)*Scale/VoxelSize(1);
    temp_Rz=R(k)*Scale/VoxelSize(3);
        
    ImageStack_k=((X-temp_xk)/temp_Rx).^2+((Y-temp_yk)/temp_Ry).^2+((Z-temp_zk)/temp_Rz).^2<=1;

    tempV=sum(sum(ImageStack_k(:)));
    
    ImageStack_k=ImageStack_k*C(k);
                    
    for m=1:AngleNum
    
        tempImageStack=ImageStack_k;
    
        for n=1:Lz
        
            tempI=ImageStack_k(:,:,n);
    
            z_n=((z1-1)*Scale+1)*(VoxelSize(3)/Scale)+(n-1)*(VoxelSize(3)/Scale);
        
            tempI=tempI*I0(m)*exp(-z_n/Depth(m));
            
            if ~isinf(PSFSigmaZ)
                tempI=tempI*exp(-0.5*(z_n/PSFSigmaZ).^2);
            end        
            
            tempImageStack(:,:,n)=tempI;
        end
            
        Im=sum(tempImageStack,3);            
     
        Im=Im/tempV;
        
        Im=imfilter(Im, PSFMask_xy, 'symmetric');     
        
        I=blockproc(Im, [Scale, Scale], SumFunction);
        
        ImageStack(y1:y2, x1:x2, m)=ImageStack(y1:y2, x1:x2, m)+I;
    end
end

bg=zeros(AngleNum, 1);

for m=1:AngleNum
    
    bg(m)=I0(m)*min(Depth(m), CellThickness)*BgDensity;
        
    %bg(m)=bg(m)*2*pi*PSFSigma^2; % BgDensity=BgDensity_real*2*pi*PSFSigma^2;
         
    I=ImageStack(:,:,m);

    I=I+bg(m);
    I=poissrnd(I);
    I=I+sqrt(I).*randn(Ly,Lx); 
    I=round(I*SysGain)/SysGain;
    
    ImageStack(:,:,m)=I;
end
    