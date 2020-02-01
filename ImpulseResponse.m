%clear; %close all; clc 
addpath('./_auxiliaryfiles')
addpath('./_varpack')
addpath('./_data')
addpath('./_datautilities')
addpath('./_functions')
addpath('./_functions/_logdensities')

warning off

load P2_copom_v5res2 B Omega Mu m p n1 n2 n Names mu0 Y1o Y2o pos_sh pos_graph
hor_irf = 47;
ss_priors = 1;

cc = zeros(p,1);
cc(1) = 1;
shock = eye(n,n1+n2);

clear YF bb




%%

turnoff_ipc = 0;

    for kk = 1:size(B,2)
        kk

        A0  = reshape(Omega(:,kk),n,n);
        A0  = chol(A0,'lower');
        
        if turnoff_ipc; A0(12,pos_sh) = 0; end 

    
        if ss_priors
            cF = reshape(B(:,kk),m-1,n)';
            
            if turnoff_ipc; cF(12,:) = 0; end 
            
            F = [cF; eye(n*(p-1),n*p)]; %
        else
            cF = reshape(B(:,kk),m,n)';
            F = [cF(:,2:end); eye(n*(p-1),n*p)];
        end
            cO = kron(cc,A0);


    for r = [pos_sh pos_sh-1]%:n1+n2      
            Yf = cO*shock(:,r);
            for i = 2:hor_irf
               yf = F*Yf(:,i-1);
               Yf = [Yf yf]; 
            end
            Yf = Yf(1:n,:);

            YF(:,:,r,kk) = Yf;
    end

    end
    


%%    PICTURE COPOM, ENE 2019

% close all
YFf = YF;

% names1 = {'US: Output gap','US: Shadow policy rate','US: Term premium','US: Fin. stress'};
% names2 = {'Brecha del producto','Brecha de inversión','R','Brecha de F','dtcn','MX: Term premium','embi','Tasa de créditos nuevos (empr. grandes)'};%
% Names = [names1 names2];

nameshock = {'brecha F'};


fullscreen = get(0,'ScreenSize');
x = linspace(0,hor_irf-1,hor_irf); 

I = [14];%[2 4 9 11]; %variables for cumsum 

% fnr_fnr = prctile(squeeze(YFf(10,:,10,:))',[50])';

F1= figure('Position',[0 0 fullscreen(3)*1.3 fullscreen(4)]);
set(gcf,'color','w');

    j=1;

shF = 1;

if shF    
%Choque de F
r = pos_sh;
f_f   = prctile(squeeze(YFf(pos_sh,:,pos_sh,:))',[50])';
else 
%Choque de IPC
r = pos_sh-1;
ipc_ipc = prctile(squeeze(YFf(pos_sh-1,:,pos_sh-1,:))',[50])';
f_ipc = prctile(squeeze(YFf(pos_sh,:,pos_sh-1,:))',[50])';
end

YFf(12,:,:,:) = YFf(12,:,:,:)*100;    

pos_graph = n1+[4,8,2,1];

for i = pos_graph
                
    
%rescalado con lo que vio en noviembre de 2019
if shF
    Yfi     = squeeze(YFf(i,:,r,:))/f_f(1)*(-1);         
else
    Yfi     = squeeze(YFf(i,:,r,:))/f_ipc(1)*(-1);
end        
        if ismember(i,I); Yfi = cumsum(Yfi); end

        pYfi    = percentile(Yfi',[16 84])'; 
        medYfi  = percentile(Yfi',50)';

        subplot(1,4,j)    
                
        fill( [x fliplr(x)],  [pYfi(:,1)' fliplr(pYfi(:,2)')],rgb('lightsalmon'), 'edgecolor', 'none'); % shaded area between 16-84 percentiles       
        hold on
        
        plot(x,medYfi,'color',rgb('darkorange'),'LineWidth', 2.5); % mediana
        hold on
        
        plot(x,zeros(1,hor_irf),'.','color','k','LineWidth', .2) % eje x 
        %alpha(.5);


        xlim([0 36]);        
         set(gca,'XTick',0:6:36)
    
        if j == 1
            ylim([-1 1]); 
            set(gca,'YTick',-1:.3:1)
        elseif j == 2
            ylim([-20 20]); 
            set(gca,'YTick',-20:5:20)
        elseif j == 3
            ylim([-.6 .6]); 
            set(gca,'YTick',-.6:.1:.6)
        elseif j == 4
            ylim([-.3 .3]); 
            set(gca,'YTick',-.3:.1:.3)
        end

        hold off

        
         set(gca,'fontsize',22,'fontname','calibri');
         t = title(Names{i});
         set(t,'fontname','Calibri','fontsize',18','fontweight','normal');%,'interpreter','latex')
         ylabel({'%'},'fontsize',16,'fontname','calibri','FontAngle','italic')
         xlabel({'meses'},'fontsize',18,'fontname','calibri','FontAngle','italic')

%          if j == 1; ylabel({'Choque a F2';'%'},'fontsize',18,'fontname','calibri','FontAngle','italic'); end
         if j == 3; ylabel({'puntos base'},'fontsize',18,'fontname','calibri','FontAngle','italic'); end

         if j == 2 
%              t = title('$\varepsilon_t$'); 
%              t = title('Brecha de F2 sintética'); 
%              t = title('Brecha de F'); 
%             set(t,'fontname','Times','fontsize',20','fontweight','normal','interpreter','latex')
         end
         if j == 1 
%              t = title('$\varepsilon_t$'); 
%              t = title('\Delta IPC'); 
%             set(t,'fontsize',20','fontweight','normal','interpreter','latex')
         end
         
         j=j+1;
%          
%          if j == 2; j=j+1; end

end

save_graph = 0;
if shF
fn = ['IRFsCOPOMv5_F'];
else
fn = ['IRFsCOPOMv5_IPC'];
end

if save_graph
img = getframe(gcf);
imwrite(img.cdata, [fn, '.png']);
end


%%

%extract IRF data 

if shF
y_f  = prctile(squeeze(YFf(6,:,pos_sh,:))',[16 50 84])'/f_f(1)*(-1);
i_f  = prctile(squeeze(YFf(7,:,pos_sh,:))',[16 50 84])'/f_f(1)*(-1);
R_f  = prctile(squeeze(YFf(8,:,pos_sh,:))',[16 50 84])'/f_f(1)*(-1);
ipc_f = prctile(squeeze(YFf(9,:,pos_sh,:))',[16 50 84])'/f_f(1)*(-1);
f_f = prctile(squeeze(YFf(10,:,pos_sh,:))',[16 50 84])'/f_f(1)*(-1);
tpp_f = prctile(squeeze(YFf(13,:,pos_sh,:))',[16 50 84])'/f_f(1)*(-1);
else
y_ipc  = prctile(squeeze(YFf(6,:,pos_sh-1,:))',[16 50 84])'/f_ipc(1)*(-1);
i_ipc  = prctile(squeeze(YFf(7,:,pos_sh-1,:))',[16 50 84])'/f_ipc(1)*(-1);
R_ipc  = prctile(squeeze(YFf(8,:,pos_sh-1,:))',[16 50 84])'/f_ipc(1)*(-1);
ipc_ipc = prctile(squeeze(YFf(9,:,pos_sh-1,:))',[16 50 84])'/f_ipc(1)*(-1);
f_ipc = prctile(squeeze(YFf(10,:,pos_sh-1,:))',[16 50 84])'/f_ipc(1)*(-1);
tpp_ipc = prctile(squeeze(YFf(13,:,pos_sh-1,:))',[16 50 84])'/f_ipc(1)*(-1);
end