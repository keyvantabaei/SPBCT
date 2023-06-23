%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Keyvan Tabaei - 4014245002 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% simple backprojection - 15 April 2023 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% start - 15 April 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% last eddit - 17 April 2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;
A=7; %cm major axis of ecllipse
B=3; %cm minor axis of ecllipse
density=10; %g/cm3 ecllipse density
Fov = 2.*A +7; %field of view in 1 dimention
D_No = 200; %number of detectors in array detector
rot_degree = 360; %projections Arc
D_delta = Fov./D_No; % distance between each detector in the detector Array
Projections=zeros(D_No,rot_degree); %create projection matrix
pixel_No=140; % define number of constructed image pixels
pixel_size = Fov./pixel_No; % calculate each pixel size
Image=zeros(pixel_No); %define image matrix for simple back projection reconstruction
line_sampling_step = 0.001; % line sampling step for each ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate Projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:rot_degree
    theta=(i./180).*pi; %convert to radian
    a_theta = A.^2.*cos(theta).^2 + B.^2.*sin(theta).^2;
    for j=1:D_No
        t= (-Fov - D_delta)./2 + (j-1).*D_delta; %caculate detector position
        sqrt_value =a_theta - t.^2;
        if  sqrt_value > 0
        Projections(j,i)=(2.*A.*B.*density./a_theta).*sqrt(sqrt_value); % when sqrt value is real 
        else
        Projections(j,i)=0;  % when sqrt value is not real 
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:rot_degree %iterate through angles
    theta=(j./180).*pi; %convert to radian
    for i=1:D_No %iterate through detectors
        T = (-Fov - D_delta)./2 + (i-1).*D_delta; %position of detector
        for x=-Fov./2 +0.07:line_sampling_step:Fov./2 %iterate through x positions
            y=(T-x*cos(theta))/sin(theta); %calculate y position
            i_index=ceil(x./pixel_size)+(pixel_No./2); %calculate i index of image matrix
            j_index=ceil(y./pixel_size)+(pixel_No./2); %calculate j index of image matrix
            if abs(y) < Fov./2 %if y is between [-Fov./2 , Fov./2]
                Image(i_index,j_index)=Image(i_index,j_index)+Projections(i,j);
            end
        end
    end
end

Radon_Image=iradon(Projections,0:359); 
subplot(1,3,1);imshow(Radon_Image);title('iRadon Image');
subplot(1,3,2);imshow(Image',[]);title('Simple Back Projection Algorithm Image');
subplot(1,3,3);imshow(Radon_Image-Image',[]);title('Diffrences');