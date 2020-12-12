function outputPara= New_sar_TRTVp(inputPara)
p=inputPara.p;  % p value
lamda=inputPara.lamda;   %parameter that balances likelihood and prior term weighting
f=inputPara.destroyedimage;  % noisy input grayscale image
MAXITR = inputPara.MAXITR;  %max dca number
tol_out= inputPara.tol_out;  % outer tollerence
TR_value= inputPara.TR_value;  % the truncated function with the threshold TR_value
r_t=inputPara.r_t;   %
r_w=inputPara.r_w;   %
[M,N] = size(f); 

eigsDtD = abs(psf2otf([1,-1],[M,N])).^2 + abs(psf2otf([1;-1],[M,N])).^2;
eigsA =eigsDtD+r_w./ r_t;

engery=zeros(MAXITR,1);
psnr1=zeros(MAXITR,1);
tx=zeros([M,N]);
ty=zeros([M,N]);
lamda_tx=zeros([M,N]);
lamda_ty=zeros([M,N]);
lamda_w=zeros([M,N]);
s=zeros([M,N]);
% Newton method parameters
% maxit_Newton = 3;
% tol_Newton = 1e-4;
u=f;
w=f;

figure;

for iter= 1:MAXITR
 u_old=u;      
 %%% -------sub t problem-------%%% 
        qx=DFX(u)-lamda_tx./r_t;
        qy=DFY(u)-lamda_ty./r_t;
        
        q=sqrt(qx.^2 + qy.^2);
         SL=(p*(1-p)./r_t).^(1./(2-p));
              %find s1
               for i=1:M
                  for j=1:N
                       if (p*SL.^(p-1)+r_t*(SL-q(i,j))) <0
                           s1_=s_Solver(SL,p,r_t,q(i,j));
                           sfesib=min(s1_,TR_value);
                           if (sfesib^p+r_t.*((sfesib-q(i,j)).^2)/2) <(r_t.*((-q(i,j)).^2)/2)
                              s1=sfesib;
                           else
                               s1=0;
                           end
                       else
                            s1=0;
                       end           
                      s2=max(TR_value,q(i,j));
            
                     if (s1^p+r_t.*((s1-q(i,j)).^2)/2) >(TR_value^p+r_t.*((s2-q(i,j)).^2)/2)
                         s(i,j)=s2;
                     else
                         s(i,j)=s1;
                     end
                  end
               end
               q(q == 0) = 1;
               tx = s.*qx./q;
               ty = s.*qy./q;

               
 %%% -------sub w problem-------%%% 
% x0=u;
% for k = 1:15
%             % temp = dF
%             temp1 =lamda.* (1-f./x0)+lamda_w+r_w.*(x0-u);
%             temp2 =lamda.* f./x0.^2+r_w;
%             x1 = x0-   temp1./temp2;
%             if abs(norm((x1-x0),'fro')) < 0.0001
%                 break;
%             end
%             x0 = x1;
% end
% w=x0;


b=lamda./r_w+lamda_w./r_w-u;
c=lamda./r_w.*f;
w=(sqrt(4.*c+b.^2)-b)./2;


 %%% -------sub u problem-------%%% 
         temp=- DBX(tx+lamda_tx./r_t)- DBY(ty+lamda_ty./r_t)+r_w.*(w+lamda_w./r_w)./ r_t;
         u = fft2(temp)./eigsA;
         u = real(ifft2(u));
          u(u<0) = 0.0000001;
 
%          u=min((mean2(f)/mean2(u))*u,255);
         u=min((mean2(f)/mean2(u))*u,1);
         imshow(u);
%          imshow(min((mean2(f)/mean2(u))*u,1));

%%% -------update Lambda-------%%% 
        lamda_tx=lamda_tx+r_t.*(tx-DFX(u));
        lamda_ty=lamda_ty+r_t.*(ty-DFY(u));
        lamda_w=lamda_w+r_w.*(w-u);
        
        
         engery(iter)=energy_TVp( u);
         
%          psnr1(iter)=psnr_Dong(img,min((mean2(f)/mean2(u))*u,255));
%          psnr1(iter)=psnr(min((mean2(f)/mean2(u))*u,1),img);
%          psnr1(iter)=psnr(u,img);
          if iter<=1
             diff_E=1;
          else
             diff_E = min(abs(engery(iter)-engery(iter-1)), abs((engery(iter)-engery(iter-1))/engery(iter-1)));
%             diff_E=abs((engery(iter)-engery(iter-1))/engery(iter-1));   
          end

%     diff_E= norm(u- u_old,'fro')./norm(u_old,'fro');
    disp(['Step ' num2str(iter) ': diff_E is ' num2str(diff_E) ]);

    if diff_E<tol_out
        break;
    end
        
 end


u(u>1) = 0.99;
 u(u<0) = 0.01;
% u=min((mean2(f)/mean2(u))*u,1);
outputPara.result=u;
% outputPara.psnr=psnr1(1:iter-1);
outputPara.engery=engery(1:iter-1);


function E = energy_TVp(u)
tv_u=sqrt(DFX(u).^2+DFY(u).^2);
E = lamda*mean2(f.*log(u)+u)+mean2(tv_u.^p);
end


function D= DFX(u)
D =u([2:end 1],:)-u;
end

function D= DFY(u)
D= u(:,[2:end 1])-u;   
end

function D= DBX(u)
D= u-u([end 1:end-1],:);
end

function D= DBY(u)
D = u-u(:,[end 1:end-1]);
end

function  x= s_Solver(x0,np,nr,nq)
    maxit =15;
    tol= 1e-4;
    % Newton method
    for nk = 1:maxit
        zm=np*x0.^(np-1)+nr*(x0-nq);
        zm1order=np*(np-1)*x0.^(np-2)+nr;
        x2=x0-  zm/ zm1order;
        if sum(abs( x2- x0)) < tol
           break;
        end
        x0=x2;
    end
    x=x2;
end

function  x= NewtonSolver(x0,np,na)
    maxit_Newton = 5;
    tol_Newton = 1e-4;
    % Newton method
    for nk = 1:maxit_Newton
        zm=np*x0.^(np-1)+na*(x0-1);
        zm1order=np*(np-1)*x0.^(np-2)+na;
        x1=x0-  zm/ zm1order;
        if sum(abs( x1- x0)) < tol_Newton
           break;
        end
        x0=x1;
    end
    x=x1;
end

end

