function[varargout]=VarH(varargin)
betain = varargin{1};
q_mat=varargin{2};
cycl = varargin{3};
if nargin==4
    plotq = varargin{4};
end
n = 20;
%q=1.5;



si = (dec2bin(0:2^n-1)-'0')';  %si no intervalo [0,1]
[~,m]=size(si);
si = 2*si-1; %si no intervalo [-1,1]
%Z_mat = zeros(1,2^n);
opt = '+-';
dist = 1;

s = size(betain);
Hmat = zeros(1,2^n);
%eqmat = zeros(1,2^n);

tmp=0;
if nargin==4
    if strcmp(plotq,'y')
        f1 = figure;
    end
end
% f1 = figure;
% f2 = figure;
Htot = cell(length(q_mat),length(betain));
%eqtot = cell(length(q_mat),length(betain));
figure
for j = 1:length(q_mat)

    q = q_mat(j);
    
    if s(1)~=1
        beta = betain(j,1):(betain(j,2)-betain(j,1))/10:betain(j,2);
    else
        beta=betain;
    end
    
    b = length(beta);

    % [J,h]=Jh_gen(n,4);

    
    for kk = 1:cycl
               
        [J,h]=Jh_gen(n,3);
        
        for i=1:b
        

        Z = qZpart(beta(i)*J,beta(i)*h,n,q,dist,opt);
        %Htot{j,i} = zeros(1,2^n);
           parfor ii=1:m
                si_tmp = si(:,ii);
                H = ising(J,h,si_tmp',opt,q);
                Hmat(ii) = H;
                %eqmat(ii) = qexp(q,-beta(i)*H,dist);
%                  Htmp_mat1(ii) = H*qexp(q,-beta(i)*H,dist)/(Z*(1-beta(i)*(q-1)*H));
%                  Htmp_mat2(ii) =  (H/(1-beta(i)*(q-1)*H))^2.*qexp(q,-beta(i)*H,dist)/(Z);

                Htmp_mat1(ii) = H*qexp(q,-beta(i)*H,dist)/(Z);
                Htmp_mat2(ii) = H^2.*qexp(q,-beta(i)*H,dist)/(Z);

           end
           min(Hmat)
           max(Hmat)
%               A = -qcumsum(beta(i)*Hmat,q)
             H_mean(j,i) = sum(Htmp_mat1);
             H2_mean(j,i) = sum(Htmp_mat2);
            %hist(Htmp)

            tmp=tmp+1;
            Htot{j,i} = Hmat;

        end
%          varH(j,:)=(H2_mean(j,:)-(H_mean(j,:)).^2)*(2-q)/(Z^(q-1));
         varH(j,:)=(H2_mean(j,:)-(H_mean(j,:)).^2);
         
    %beta_mat(kk) = beta(find(varH(j,:)==max(varH(j,:))));
    
    
    
    end
%     beta_mean(j) = mean(beta_mat);
%     beta_err(j) = std(beta_mat);
    
    %eqtot{j,i} = eqmat;
    
    plot(beta,abs(varH(j,:))), xlabel_=xlabel("$\beta$");ylabel_=ylabel("$<H>^2-<H^2>$");hold on,
    set(xlabel_,"fontsize",15,"interpreter","latex");
    set(ylabel_,"fontsize",15,"interpreter","latex");
end

%figure


%varH=H2_mean-(H_mean).^2;
% varargout{1}=varH;
% if nargout==3
%     varargout{2}=J;
%     varargout{3} = h;
% end
% f3 = figure;
% f4 = figure;
% f5 = figure;
% f6 = figure;
%F(b) = struct('cdata',[],'colormap',[]);

% figure
for kk = 1:length(q_mat)
    %v = VideoWriter(['21072022movie_q',num2str(q_mat(kk)),'.avi']);
    %v.FrameRate = 10;
    %%v.Quality = 95;
    %open(v);
    
    for kkk=1:b
        
        %[count, edges]=histcounts(Htot{kk,kkk},50);
        %center = conv(edges, [0.5 0.5], 'valid');
        
        %figure(f3)
        %stem(Htot{kk,kkk},Ptotmat{kk,kkk},'LineWidth',2),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))]),xlabel('H'),ylabel('P'),ax = gca;ax.FontSize = 30;
        %axis([min(Htot{kk,kkk}) max(Htot{kk,kkk}) 0 length(Htot{kk,kkk})])     
        
        %figure(f4)
%         f4.WindowState = 'maximized';
%         %f4.Visible = 'off';
%               
%         subplot(2,2,1)
%         stem(center,log(qexp(q_mat(kk),-beta(kkk)*center)),'fill'),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))],'location','southeast'),xlabel('Center of histogram of H'),ylabel('ln(e_q^{-\beta H})'),ax = gca;ax.FontSize = 20; 
%         %axis([min(center) max(center)])
%         
%         subplot(2,2,2)
%         stem(center,log(count)+log(qexp(q_mat(kk),-beta(kkk)*center)),'fill'),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))],'location','southeast'),xlabel('Center of histogram of H'),ylabel('ln(N_E)+ln(e_q^{-\beta H})'),ax = gca;ax.FontSize = 20; 
%         %axis([min(center) max(center)])
%         
%         subplot(2,2,3)
%         stem(center,log(count),'fill'),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))]),xlabel('Center of histogram of H'),ylabel('ln(N_E)'),ax = gca;ax.FontSize = 20;
%         %axis([min(center) max(center) 0 max(log(count))])
%         
%         subplot(2,2,4)
        %edgesx = [10,30:2.6:138];
        
        
        
%         h= histogram(Htot{kk,kkk},100);legend(['q= ',num2str(q_mat(kk))],'Location','northwest'),xlabel('E'),ylabel('N_E'),ax = gca;ax.FontSize = 20;
%         %h.Normalization = 'countdensity';
%         [count, edgesx]=histcounts(Htot{kk,kkk});
%         center = conv(edgesx, [0.5 0.5], 'valid');
%         axis 'square'
%         set(gca,'Fontsize',25)
%         
%         figure(f3)
%         plot(center,log(count)+log(qexp(q,-beta(kkk)*center,dist))), hold on, legend(cellstr(num2str(beta', 'beta=%.2f')))
%         %axis([min(Htot{kk,kkk}) max(Htot{kk,kkk}) 0 length(Htot{kk,kkk})])
%         
%         figure(f5)
%         plot(center,log(count))%, hold on, legend(cellstr(num2str(beta', 'beta=%.2f')))
%         
%         figure(f6)
%         plot(center,log(qexp(q,-beta(kkk)*center,dist))), hold on, legend(cellstr(num2str(beta', 'beta=%.2f')))
        


        %drawnow
        %F = getframe(gcf);
        %size(F.cdata)
        %writeVideo(v,F);
    
        
        %pause(0.1)
    end
    %close(v)
    %save(['movie_q',num2str(q_mat(kk)),'.mat'],'F')
end

beep









%% CODIGOS DE GRÁFICOS E VIDEOS
% for kk = 1:length(q_mat)
%     v = VideoWriter(['movie_q',num2str(q_mat(kk)),'.avi']);
%     v.FrameRate = 10;
%     %v.Quality = 95;
%     open(v);
%     %F(b) = struct('cdata',[],'colormap',[]);
%     for kkk=1:b
%         figure(f3)
%         stem(Htot{kk,kkk},Ptotmat{kk,kkk},'LineWidth',2),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))]),xlabel('H'),ylabel('P'),ax = gca;ax.FontSize = 30;
%         %axis([min(Htot{kk,kkk}) max(Htot{kk,kkk}) 0 length(Htot{kk,kkk})])
%         
%         figure(f4)
% %         hist(Htot{kk,kkk}.*Ptotmat{kk,kkk})
%         Hhist = randsrc(1,length(Ptotmat{kk,kkk}),[Htot{kk,kkk};Ptotmat{kk,kkk}]);
%         hist(Hhist),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))]),xlabel('H'),ylabel('N_{E}'),ax = gca;ax.FontSize = 30;
%         axis([min(Htot{kk,kkk}) max(Htot{kk,kkk}) 0 length(Htot{kk,kkk})])
%         drawnow
% %         ax = gca;
% %         ax.Units = 'pixels';
% %         pos = ax.Position;
% %         marg = 30;
% %         rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
% %         F(kkk) = getframe(gca,rect);
% %         ax.Units = 'normalized';
%         F = getframe(gcf);
%         writeVideo(v,F);
%         pause(.001)
%     end
%     %save(['movie_q',num2str(q_mat(kk)),'.mat'],'F')
% end
% close(v)

beep


