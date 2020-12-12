clc;clear;
% path(path,'../images/')

img=double(imread('terrasar_hill_512.jpg'));
img=img(:,:,1);
[m, n] = size(img);
y=max(img,1);
f= double(y)/255;

%% TVp_AA
inputPara.p=0.5;              %   0<p<1
inputPara.destroyedimage=f;
inputPara.lamda=12;    %  If lamda is too large, the noise cannot be removed cleanly; and if lamda is too small, the fine features will be oversmoothed.

inputPara.r_t=10;  
inputPara.r_w=10;
   
inputPara.MAXITR=100; 
inputPara.TR_value=3;    % the truncated function with the threshold TR_value
inputPara.tol_out=1e-4;  %the stopping criteria
tic;
outputPara= New_sar_TRTVp( inputPara );
toc,
u=outputPara.result;
engeryarry=outputPara.engery;

figure(1);
subplot(2,2,1);
imshow(f);
title('SAR image');

subplot(2,2,2);
imshow(f);
title('sar img');

subplot(2,2,3);
imshow(u);
title('TRTVp img');
imwrite(u,'hill_TRTVp.png');

T = length(engeryarry);
subplot(2,2,4),plot(1:T,engeryarry)
xlabel('iter');ylabel('Engery');

%%% -------ratio-------%%%
img_ratio=(f+0.0001)./(u+0.0001);
figure(2);
imshow(img_ratio);
title('ratio img');
imwrite(img_ratio,'ratio_hill_TRTVp.png');

 %%% -------enl-------%%% 
start_row=280;
end_row=300;
start_col=400;
end_col=450;
figure(3);
hold on;
imagesc(u);
line([start_col,end_col],[start_row,start_row],'Color','r','LineWidth',2)   
line([start_col,end_col],[end_row,end_row],'Color','r','LineWidth',2)   
line([start_col,start_col],[start_row,end_row],'Color','r','LineWidth',2)   
line([end_col,end_col],[start_row,end_row],'Color','r','LineWidth',2)   
hold off
I=u(start_row:end_row,  start_col:end_col);
av=mean2(I);
st=std2(I);
enl=av^2/st^2


 %%% -------epi-------%%% 
epi_noisy=zeros(m,n);
epi_denoisy=zeros(m,n);
for i=2:m-1
    for j=2:n-1
        epi_noisy(i,j)=abs(f(i,j)-f(i,j-1))+abs(f(i,j)-f(i,j+1))+abs(f(i,j)-f(i-1,j))+abs(f(i,j)-f(i+1,j));
        epi_denoisy(i,j)=abs(u(i,j)-u(i,j-1))+abs(u(i,j)-u(i,j+1))+abs(u(i,j)-u(i-1,j))+abs(u(i,j)-u(i+1,j));
    end
end
epi = sum(sum(epi_denoisy))./sum(sum(epi_noisy))


