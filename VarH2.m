function[varargout]=VarH2(varargin)
beta = varargin{1};
q_mat=varargin{2};
options = varargin{3};

%q=1.5;
b = length(beta);

switch options
    
    case 1
        %Load from files spinglass, 15 spins 
        n = 15;
        fiJ = fopen('J_spinglass15.txt','r');
        fih = fopen('h_spinglass15.txt','r');
        J = fscanf(fiJ,'%f',[n, n]); J=J';
        h = fscanf(fih,'%f');h=h';
        fclose(fiJ);
        fclose(fih);
        
    case 2
        %Load from files ferro, 15 spins 
        n = 15;
        fiJ = fopen('J_ferro15.txt','r');
        fih = fopen('h_ferro15.txt','r');
        J = fscanf(fiJ,'%f',[n, n]); J=J';
        h = fscanf(fih,'%f');h=h';
        fclose(fiJ);
        fclose(fih);
        
    case 3
        %Load from files antiferro, 15 spins 
        n = 15;
        fiJ = fopen('J_antiferro15.txt','r');
        fih = fopen('h_antiferro15.txt','r');
        J = fscanf(fiJ,'%f',[n, n]); J=J';
        h = fscanf(fih,'%f');h=h';
        fclose(fiJ);
        fclose(fih);
        
    case 4
        %Load from files spinglass, 15 spins 
        n = 20;
        fiJ = fopen('J_spinglass20.txt','r');
        fih = fopen('h_spinglass20.txt','r');
        J = fscanf(fiJ,'%f',[n, n]); J=J';
        h = fscanf(fih,'%f');h=h';
        fclose(fiJ);
        fclose(fih);
        
    case 5
        %Load from files ferro, 15 spins 
        n = 20;
        fiJ = fopen('J_ferro20.txt','r');
        fih = fopen('h_ferro20.txt','r');
        J = fscanf(fiJ,'%f',[n, n]); J=J';
        h = fscanf(fih,'%f');h=h';
        fclose(fiJ);
        fclose(fih);
        
    case 6
        %Load from files antiferro, 15 spins 
        n = 20;
        fiJ = fopen('J_antiferro20.txt','r');
        fih = fopen('h_antiferro20.txt','r');
        J = fscanf(fiJ,'%f',[n, n]); J=J';
        h = fscanf(fih,'%f');h=h';
        fclose(fiJ);
        fclose(fih);
        
    case 7
        n=10;
        [J,h]=Jh_gen(n,3);
        
    case 8
        n = 20;
        [Jname, Jpath] = uigetfile('*.txt','Select the file containing J');
        [hname, hpath] = uigetfile('*.txt','Select the file containing h');
        n = input('Write the number of spins used');
        fiJ = fopen(strcat(Jpath,Jname),'r');
        fih = fopen(strcat(hpath,hname),'r');
        for i =1:n
            for ii = 1:n
                J(i,ii) = fscanf(fiJ,'%f',[1,1]);
            end
        end
        h = fscanf(fih,'%f ');
        h = h';
end
    
    
    
    
    
    
%[J,h]=Jh_gen(n,3);
% J=J*1;
% h=h*1;
%[J,h]=Jh_gen(n,3)%1 for -1 a 1| 3 for 0,1

H_mean = zeros(length(q_mat),b);
H2_mean = zeros(length(q_mat),b);
Ptot = zeros(length(q_mat),b);
Htot = cell(length(q_mat),b);
Ptotmat = cell(length(q_mat),b);
varH = zeros(length(q_mat),b);
sigmaMean = zeros(length(q_mat),b);%magnetização

si = (dec2bin(0:2^n-1)-'0')';  %si no intervalo [0,1]
[~,m]=size(si);

si = 2*si-1; %si no intervalo [-1,1]
Z_mat = zeros(1,2^n);
opt = '+-';
dist = 1;


% si_mean_exp = sum((si'+1)/2)/m;
% 
% sij_mean_exp = zeros(n);
% 
% for k = 1:m
%     sij_mean_exp = sij_mean_exp + (sij_gen(si(:,k))+1)/2;
% end
% sij_mean_exp=sij_mean_exp/m;
% sij_mean_exp = triu(sij_mean_exp,1);



tmp=0;
f1 = figure;
f2 = figure;
% f2 = figure;
for j = 1:length(q_mat)
    q = q_mat(j);
    parfor i=1:b
    Htmp_mat1=zeros(1,m);
    Htmp_mat2=zeros(1,m);
    Ptmp = zeros(1,m);
    Htmp = zeros(1,m);
    sigmatmp = zeros(m,n);
    Z = qZpart(beta(i)*J,beta(i)*h,n,q,dist,opt);
      
       for ii=1:m
            si_tmp = si(:,ii);
            H = ising(J,h,si_tmp',opt,q);
            
%             Htmp_mat1(ii) = H*qexp(q,-beta(i)*H,dist)/(Z*(1-beta(i)*(q-1)*H));
%             Htmp_mat2(ii) =  (H/(1-beta(i)*(q-1)*H))^2.*qexp(q,-beta(i)*H,dist)/(Z);
            
            Htmp_mat1(ii) = H*qexp(q,-beta(i)*H,dist)/(Z);
            Htmp_mat2(ii) = H^2.*qexp(q,-beta(i)*H,dist)/(Z);
            Htmp(ii) = H;
            Ptmp(ii) = qexp(q,-beta(i)*H,dist)/(Z);
            sigmatmp(ii,:) = si_tmp*Ptmp(ii);
%             Htmp_mat1(ii) = H*(qexp(q,-H,dist)/(Z)).^q;
%             Htmp_mat2(ii) = H^2.*(qexp(q,-H,dist)/(Z)).^q;
%             Ptmp(ii) = (qexp(q,-H,dist)/(Z)).^q;
           
       end
       
        sigmaMean(j,i) = mean(sum(sigmatmp,1));
        H_mean(j,i) = sum(Htmp_mat1);
        H2_mean(j,i) = sum(Htmp_mat2);
        Ptot(j,i) = sum(Ptmp);
        Htot{j,i} = Htmp;
        Ptotmat{j,i} = Ptmp;
        %hist(Htmp)
        
        tmp=tmp+1;

     
    end
%     varH(j,:)=(H2_mean(j,:)-(H_mean(j,:)).^2)*(2-q)/(Z^(1-q));
     varH(j,:)=(H2_mean(j,:)-(H_mean(j,:)).^2);
     
    
    figure(f1)
    plot(beta,abs(varH(j,:))/max(abs(varH(j,:))),'LineWidth',2), hold on
    legenda=cell(1,length(q_mat));
    for i = 1:length(q_mat)
        legenda{i} = char(compose('q = %.2f',q_mat(i)));
    end
    legend(legenda)
    xlabel_=xlabel("$|J_{max}|$");
    ylabel_=ylabel("Normalized $(\Delta H)^2$");
    title_=title("Energy fluctuation $\Big((\Delta E)^2 = [<H(\sigma)>^2-<H(\sigma)^2>]\Big)$ vs $\beta$ ($|J_{max}|$) for Ferromagnetic in random Field");
    set(title_,"fontsize",12,"interpreter","latex");
    set(xlabel_,"fontsize",20,"interpreter","latex");
    set(ylabel_,"fontsize",20,"interpreter","latex");

    figure(f2)
    plot(beta,abs(sigmaMean(j,:))/max(abs(sigmaMean(j,:))),'LineWidth',2), hold on
    legend(legenda)
    xlabel_=xlabel("$|J_{max}|$");
    ylabel_=ylabel("Normalized $\langle\sigma\rangle$");
    title_=title("Mean magnetization $\Big(\sum_{i=1}^{n}\langle\sigma\rangle_i\Big)$ vs $\beta$ ($|J_{max}|$) for Ferromagnetic in random Field");
    set(title_,"fontsize",12,"interpreter","latex");
    set(xlabel_,"fontsize",20,"interpreter","latex");
    set(ylabel_,"fontsize",20,"interpreter","latex");
%     figure(f2)
%     plot(beta,H_mean(j,:).^2,'LineWidth',2), hold on
%     legenda=cell(1,length(q_mat));
%     for i = 1:length(q_mat)
%         legenda{i} = char(compose('q = %.2f',q_mat(i)));
%     end
%     legend(legenda)
%     xlabel_=xlabel("$\beta$");
%     ylabel_=ylabel("$\langle H\rangle^2$");
%     %title_=title("Gr\'afico de calor espec\'ifico $c_{\beta} = \frac{q-2}{Z_{2-q}^{q-1}}[<f_{TS}(H(\sigma))>^2-<f_{TS}^2(H(\sigma))>]$");
%     title_=title("Gr\'afico de calor espec\'ifico $(\Delta E)^2 = [<H(\sigma)>^2-<H(\sigma)^2>]$");
%     set(title_,"fontsize",15,"interpreter","latex");
%     set(xlabel_,"fontsize",20,"interpreter","latex");
%     set(ylabel_,"fontsize",20,"interpreter","latex");
    
end

%varH=H2_mean-(H_mean).^2;
varargout{1}=varH;
if nargout==3
    varargout{2}=J;
    varargout{3} = h;
end

% for kk = 1:length(q_mat)
%     v = VideoWriter(['movie_q_justbeta',num2str(q_mat(kk)),'.avi']);
%     v.FrameRate = 2.1;
%     %v.Quality = 95;
%     open(v);
%     %F(b) = struct('cdata',[],'colormap',[]);
%     for kkk=1:b
%         fh = figure();
%         
%         plot(beta,abs(varH(j,:))/max(abs(varH(j,:))),'LineWidth',4), hold on
%         plot(beta(kkk)*ones(1,11), 0:.1:1,'-r' ,'LineWidth',3)
%         
%         axis([0 2 0 1])
% %         legenda=cell(1,length(q_mat));
% %         for i = 1:length(q_mat)
% %             legenda{i} = char(compose('q = %.2f',q_mat(i)));
% %         end
%         legend({"Energy fluctuation for q=0.7",char(compose("beta = %.1f",beta(kkk)))},'Location','southeast');
%         xlabel_=xlabel("$\beta$");
%         ylabel_=ylabel("Normalized $(\Delta E)^2$");
%         title_=title("Energy fluctuations $(\Delta E)^2 = [<H(\sigma)>^2-<H(\sigma)^2>]$");
%         set(title_,"interpreter","latex");
%         set(xlabel_,"interpreter","latex");
%         set(ylabel_,"interpreter","latex");
%         %set(legend_, "fontsize",20,"interpreter","latex")
%         
%         set(gca,'FontSize',30)
%         fh.WindowState = 'maximized';
%         
%         %fh.FontSize = 30;
%         pause(0.1)
%         F = getframe(gcf);
%         writeVideo(v,F);
%         
%     end
%     %save(['movie_q',num2str(q_mat(kk)),'.mat'],'F')
% end
% close(v)






% %f3 = figure;
% f4 = figure;
% %f5 = figure;
% %f6 = figure;
% %F(b) = struct('cdata',[],'colormap',[]);
% 
% 
% for kk = 1:length(q_mat)
%     v = VideoWriter(['11052022movie_q',num2str(q_mat(kk)),'.avi']);
%     v.FrameRate = 10;
%     %v.Quality = 95;
%     open(v);
%     
%     for kkk=1:b
%         [count, edges]=histcounts(Htot{kk,kkk},50);
%         center = conv(edges, [0.5 0.5], 'valid');
%         
%         %figure(f3)
%         %stem(Htot{kk,kkk},Ptotmat{kk,kkk},'LineWidth',2),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))]),xlabel('H'),ylabel('P'),ax = gca;ax.FontSize = 30;
%         %axis([min(Htot{kk,kkk}) max(Htot{kk,kkk}) 0 length(Htot{kk,kkk})])     
%         
%         figure(f4)
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
%         stem(center,log(count),'fill'),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))],'location','southeast'),xlabel('Center of histogram of H'),ylabel('ln(N_E)'),ax = gca;ax.FontSize = 20;
%         %axis([min(center) max(center) 0 max(log(count))])
%         
%         subplot(2,2,4)
%         histogram(Htot{kk,kkk},50),legend(['q= ',num2str(q_mat(kk)),'  \beta= ',num2str(beta(kkk))]),xlabel('H'),ylabel('N_E'),ax = gca;ax.FontSize = 20;
%         %axis([min(Htot{kk,kkk}) max(Htot{kk,kkk}) 0 length(Htot{kk,kkk})])
%         
%         drawnow
%         F = getframe(gcf);
%         size(F.cdata)
%         writeVideo(v,F);
%     
%         
%         %pause(0.1)
%     end
%     close(v)
%     %save(['movie_q',num2str(q_mat(kk)),'.mat'],'F')
% end

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
% %plot(beta,abs(varH)), xlabel_=xlabel("$\beta$");ylabel_=ylabel("$<H>^2-<H^2>$");
% %set(xlabel_,"fontsize",15,"interpreter","latex");
% %set(ylabel_,"fontsize",15,"interpreter","latex");
% beep

