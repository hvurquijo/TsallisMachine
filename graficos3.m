function[]=graficos3(var1,qmat)

if var1 == 0 || var1 == 2 
    dir = uigetdir;

    fi1 = fopen([dir,'/<E_h>_vs_Cycl.dat'],'r');
    fi2 = fopen([dir,'/<E_J>_vs_Cycl.dat'],'r');
    fi3 = fopen([dir,'/q_matrix.dat'],'r');
    fi4 = fopen([dir,'/<Sijk_vs_Cycl>.dat'],'r');
    fi5 = fopen([dir,'/<Sij_vs_Cycl>.dat'],'r');
    fi6 = fopen([dir,'/<Si>_vs_Cycl.dat'],'r');

    length(dir)
tit = dir(52:end-5);
q_mais=1.1;
q_menos = 0.5;

q = 0.8;
beta = 0.02;
n=20;
%m=1000;

q_mat= fscanf(fi3,'%f');
q_mat = round(10*q_mat)/10;
tmp = length(q_mat);
Eh = fscanf(fi1,'%f');
whos Eh

cycl = length(Eh)/tmp;
fclose(fi1);
fi1 = fopen([dir,'/<E_h>_vs_Cycl.dat'],'r');

Eh = fscanf(fi1,'%f',[cycl, tmp]); Eh=Eh';
EJ = fscanf(fi2,'%f',[cycl, tmp]); EJ =EJ';
Sijk = fscanf(fi4,'%f',[cycl, tmp]); Sijk =Sijk';
Sij = fscanf(fi5,'%f',[cycl, tmp]); Sij =Sij';
Si = fscanf(fi6,'%f',[cycl, tmp]); Si =Si';


% x = ['q',num2str(q),'N',num2str(n)];
% if ~exist(x, 'dir')
%    mkdir(x)
% end

fclose all;
f=figure;
plot(q_mat,log(Eh(:,end)), 'b-s','LineWidth',3, 'MarkerSize',9),hold on%,plot(q_mat,log(Eh(:,end)), 'bo','MarkerSize',5),
xlabel("$\mathit{q_{TS_{<}}}$", 'FontSize',20, 'interpreter', 'latex'), 
ylabel("$log\big(MAE(h_i)\big)$", 'FontSize',20, 'interpreter','latex')
%hold on
%find(log(Eh(:,end))==min(log(Eh(:,end))))
%disp("Aqui")
plot(q,log(Eh(find(q_mat==q),end)),'ro','MarkerSize',20,'LineWidth',3)
legend({['q= ',num2str(q),'; \beta= ',num2str(beta)],'Position of q'})
axis([0.1 2 -35 10])
ax10 = gca; 
ax10.FontSize = 26;
% f.WindowState = 'maximized';
% grid on
% set(f,'color','w');
% drawnow
% im = getframe(f);
% imwrite(im.cdata,[x,'/h_beta',num2str(beta,'%3.1f'),'.jpg'])
% pause(1)


f=figure;
%q_mat = [0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,1.05,1.10];
plot(q_mat,(EJ(:,end)), 'b-s',q_mat,(Eh(:,end)), 'r--s',q_mat,(Sijk(:,end)), 'g-s','LineWidth',3, 'MarkerSize',9),
xlabel("$\mathit{q_{TS}}$", 'FontSize',20, 'interpreter','latex'), 
ylabel("Rms", 'FontSize',20)
legend({'rms J_{ij}','rms h_i','rms T_{ijk}'})
% hold on
% plot(q,log(EJ(find(q_mat==q),end)),'ro','MarkerSize',20,'LineWidth',3)
% legend({['q= ',num2str(q),'; \beta= ',num2str(beta)],'Position of q'})
axis([0.1 2 -35 10])
ax10 = gca; 
ax10.FontSize = 26;
axis 'square'
set(gca,'linewidth',3)

% f.WindowState = 'maximized';
% grid on
% set(f,'color','w');
% drawnow
% im = getframe(f);
% imwrite(im.cdata,[x,'/J_beta',num2str(beta,'%3.1f'),'.jpg'])
% pause(1)

f=figure;
plot(q_mat,log(Sijk(:,end)), 'b-s','LineWidth',3, 'MarkerSize',9),
xlabel("$\mathit{q_{TS_{<}}}$", 'FontSize',20, 'interpreter','latex'), 
ylabel("$log\big(MAE(T_{ijk})\big)$", 'FontSize',20, 'interpreter','latex')
legend(['q= ',q,'; \beta= ',beta])
hold on
plot(q,log(Sijk(find(q_mat==q),end)),'ro','MarkerSize',20,'LineWidth',3)
legend({['q= ',num2str(q),'; \beta= ',num2str(beta)],'Position of q'})
axis([0.1 2 -40 10])
ax10 = gca; 
ax10.FontSize = 26;

% f.WindowState = 'maximized';
% grid on
% set(f,'color','w');
% drawnow
% im = getframe(f);
% imwrite(im.cdata,[x,'/T_beta',num2str(beta,'%3.1f'),'.jpg'])
% pause(1)

f=figure;
plot(q_mat,log(Sij(:,end)), 'b-s','LineWidth',3, 'MarkerSize',9),
xlabel("$\mathit{q_{TS_{<}}}$", 'FontSize',20, 'interpreter','latex'), 
ylabel("$log\big(MAE(\sigma_i\sigma_j)\big)$", 'FontSize',20, 'interpreter', 'latex')
legend(['q= ',q,'; \beta= ',beta])
hold on
plot(q,log(Sij(find(q_mat==q),end)),'ro','MarkerSize',20,'LineWidth',3)
legend({['q= ',num2str(q),'; \beta= ',num2str(beta)],'Position of q'})
axis([0.1 2 -40 10])
ax10 = gca; 
ax10.FontSize = 26;
% f.WindowState = 'maximized';
% grid on
% set(f,'color','w');
% drawnow
% im = getframe(f);
% imwrite(im.cdata,[x,'/Sij_beta',num2str(beta,'%3.1f'),'.jpg'])
% pause(1)

f=figure;
plot(q_mat,log(Si(:,end)), 'b-s','LineWidth',3, 'MarkerSize',9),

xlabel("$\mathit{q_{TS_{<}}}$", 'FontSize',20, 'interpreter','latex'), 
ylabel("$log\big(MAE(\sigma_i)\big)$", 'FontSize',20, 'interpreter','latex')
legend(['q= ',q,'; \beta= ',beta])
hold on
plot(q,log(Si(find(q_mat==q),end)),'ro','MarkerSize',20,'LineWidth',3)
legend({['q= ',num2str(q),'; \beta= ',num2str(beta)],'Position of q'})
axis([0.1 2 -40 10])
ax10 = gca; 
ax10.FontSize = 26;


% f.WindowState = 'maximized';
% grid on
% set(f,'color','w');
% drawnow
% im = getframe(f);
% imwrite(im.cdata,[x,'/Si_beta',num2str(beta,'%3.1f'),'.jpg'])
% pause(1)

for j = 1:length(q_mat)
    legenda{j} = num2str(q_mat(j));
end
figure
plot(log(Eh'),'LineWidth',4), 
hleg=legend(legenda(1:end));
htitle = get(hleg,'Title');
set(htitle, 'String',"$q_{TS}$",'interpreter','latex')
xlabel("\textit{Ciclos de aproxima\c{c}\~ao}", 'FontSize',40,'interpreter','latex')
ylabel("$\mathit{log\Big[\big\langle|\Delta h|\big\rangle\Big]}$", 'FontSize',40, 'interpreter','latex')
ax10 = gca; 
ax10.FontSize = 26;
set(hleg,"FontSize",25)
% fi  = fopen('h_CyclComparitionTime.txt','a');
% fprintf(fi,'%f ',Eh');
% fprintf(fi,'%.2f\n',beta);
% fclose(fi);

figure
plot(log(EJ'),'LineWidth',4)
hleg=legend(legenda(1:end));
htitle = get(hleg,'Title');
set(htitle, 'String',"$q_{TS}$",'interpreter','latex')
xlabel("\textit{Ciclos de aproxima\c{c}\~ao}", 'FontSize',40,'interpreter','latex')
ylabel("$\mathit{log\Big[\big\langle|\Delta J|\big\rangle\Big]}$", 'FontSize',40, 'interpreter','latex')
ax10 = gca; 
ax10.FontSize = 26;
set(hleg,"FontSize",25)
% fi  = fopen('J_CyclComparitionTime.txt','a');
% fprintf(fi,'%f ',EJ');
% fprintf(fi,'%.2f\n',beta);
% fclose(fi);

% 
% 
% 
% figure
% plot(log(Sijk'),'LineWidth',4)
% hleg=legend(legenda(1:end));
% htitle = get(hleg,'Title');
% set(htitle, 'String',"$q_{TS}$",'interpreter','latex')
% xlabel("$\textit{Ciclos de aproxima\c{c}\~ao}$", 'FontSize',40,'interpreter','latex'), 
% ylabel("$\mathit{log\Big[\big\langle|\Delta \langle\sigma_i\sigma_j\sigma_k\rangle|\big\rangle\Big]}$", 'FontSize',40, 'interpreter', 'latex')
% ax10 = gca; 
% ax10.FontSize = 26;
% set(hleg,"FontSize",25)
% 
% figure
% plot(log(Sij'),'LineWidth',4),
% hleg=legend(legenda(1:end));
% htitle = get(hleg,'Title');
% set(htitle, 'String',"$q_{TS}$",'interpreter','latex')
% xlabel("$\textit{Ciclos de aproxima\c{c}\~ao}$", 'FontSize',40,'interpreter','latex'), 
% ylabel("$\mathit{log\Big[\big\langle|\Delta \langle\sigma_i\sigma_j\rangle|\big\rangle\Big]}$", 'FontSize',40, 'interpreter', 'latex')
% ax10 = gca; 
% ax10.FontSize = 26;
% set(hleg,"FontSize",25)
% 
% figure
% plot(log(Si'),'LineWidth',4)
% hleg=legend(legenda(1:end));
% htitle = get(hleg,'Title');
% set(htitle, 'String',"$q_{TS}$",'interpreter','latex')
% xlabel("$\textit{Ciclos de aproxima\c{c}\~ao}$", 'FontSize',30,'interpreter','latex'),
% ylabel("$\mathit{log\Big[\big\langle|\Delta \langle\sigma_i\rangle|\big\rangle\Big]}$", 'FontSize',40, 'interpreter', 'latex')
% ax10 = gca; 
% ax10.FontSize = 26;
% set(hleg,"FontSize",25)
%     
end

%marker_list= ['h','v','^','o','+','s','p','<'];




q_position = 3;


if var1 == 1 || var1==2
    n=20;
    res = length(qmat);
    
    if res>7
        marker_list= ['p','d','^','o','+','x','s','p','d','^','o','+','x','s'];
        color_list = {[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],[0.6350 0.0780 0.1840],[0.85 0 0],[0 0.7 0],[0 0 0.75],[0.9290 0.6940 0.1250],[rand(1),rand(1),rand(1)],[rand(1),rand(1),rand(1)],[rand(1),rand(1),rand(1)],[rand(1),rand(1),rand(1)],[rand(1),rand(1),rand(1)],[rand(1),rand(1),rand(1)],[rand(1),rand(1),rand(1)]};
    else
        marker_list= ['p','d','^','o','+','x','s'];
        color_list = {[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],[0.6350 0.0780 0.1840],[0.85 0 0],[0 0.7 0],[0 0 0.75],[0.9290 0.6940 0.1250]};
    end
    dir = uigetdir;
    
    file_J_ = strcat(dir,'/J.dat');
    file_h_ = strcat(dir,'/h.dat');
    
    for i=1:res
        file_si_exp = strcat(dir,['/si_mean_exp_',num2str(qmat(i),'%.2f'),'.dat']);
        file_si_TS = strcat(dir,['/si_mean_TS_',num2str(qmat(i),'%.2f'),'.dat']);
        file_sij_exp = strcat(dir,['/sij_mean_exp_',num2str(qmat(i),'%.2f'),'.dat']);
        file_sij_TS = strcat(dir,['/sij_mean_TS_',num2str(qmat(i),'%.2f'),'.dat']);
        file_Tijk_exp = strcat(dir,['/Tijk_mean_exp_',num2str(qmat(i),'%.2f'),'.dat']);
        file_Tijk_TS = strcat(dir,['/Tijk_mean_TS_',num2str(qmat(i),'%.2f')],'.dat');
        file_Cij_exp = strcat(dir,['/Cij_exp_',num2str(qmat(i),'%.2f'),'.dat']);
        file_Cij_TS = strcat(dir,['/Cij_TS_',num2str(qmat(i),'%.2f'),'.dat']);
        file_Cijk_exp = strcat(dir,['/Cijk_exp_',num2str(qmat(i),'%.2f'),'.dat']);
        file_Cijk_TS = strcat(dir,['/Cijk_TS_',num2str(qmat(i),'%.2f'),'.dat']);
        
        file_J_TS = strcat(dir,['/J_TS_',num2str(qmat(i),'%.2f'),'.dat']);
        file_h_TS = strcat(dir,['/h_TS_',num2str(qmat(i),'%.2f'),'.dat']);
        
        

        

        fi7 = fopen(file_si_exp,'r');
        fi8 = fopen(file_si_TS,'r');
        fi9 = fopen(file_sij_exp,'r');
        fi10 = fopen(file_sij_TS,'r');
        fi11 = fopen(file_Tijk_exp,'r');
        fi12 = fopen(file_Tijk_TS,'r');
        fi13 = fopen(file_Cij_exp,'r');
        fi14 = fopen(file_Cij_TS,'r');
        
        fi15 = fopen(file_J_TS,'r');
        fi16 = fopen(file_h_TS,'r');
        fi17 = fopen(file_J_,'r');
        fi18 = fopen(file_h_,'r');
        fi19 = fopen(file_Cijk_exp,'r');
        fi20 = fopen(file_Cijk_TS,'r');
        
        

        exp_si(:,i) = fscanf(fi7,'%f');
        TS_si(:,i) = fscanf(fi8,'%f');

        exp_sij(:,i) = fscanf(fi9,'%f');
        TS_sij(:,i) = fscanf(fi10,'%f');

        exp_Tijk(:,i) = fscanf(fi11,'%f');
        TS_Tijk(:,i) = fscanf(fi12,'%f');

        exp_Cij(:,i) = fscanf(fi13,'%f');
        TS_Cij(:,i) = fscanf(fi14,'%f');
        
        exp_h(:,i) = fscanf(fi18,'%f');
        TS_h(:,i) = fscanf(fi16,'%f');
        
        exp_J(:,i) = fscanf(fi17,'%f');
        TS_J(:,i) = fscanf(fi15,'%f');
        
        exp_Cijk(:,i) = fscanf(fi19,'%f');
        TS_Cijk(:,i) = fscanf(fi20,'%f');
        
        ft = fittype('a*x','dependent',{'y'},'independent',{'x'},'coefficients',{'a'});
        %fo = fitoptions(ft)
        %fo = fitoptions('poly1','Upper',[inf 0], 'Lower',[-inf 0])
        %fity = fittype('poly1','options',fo);

        %[~, gof_si]=fit(exp_si(:,i),TS_si(:,i),'poly1',fo);
        [coeffs_si, gof_si]=fit(exp_si(:,i),TS_si(:,i),ft);
        [coeffs_sij, gof_sij]=fit(exp_sij(:,i),TS_sij(:,i),ft);
        [coeffs_Tijk, gof_Tijk]=fit(exp_Tijk(:,i),TS_Tijk(:,i),ft);
        [coeffs_Cij, gof_Cij]=fit(exp_Cij(:,i),TS_Cij(:,i),ft);
        [coeffs_J, gof_J]=fit(exp_J(:,i),TS_J(:,i),ft);
        [coeffs_h, gof_h]=fit(exp_h(:,i),TS_h(:,i),ft);
%         [coeffs_sij, gof_sij]=fit(exp_sij(:,i),TS_sij(:,i),'poly1',fo);
%         [coeffs_Tijk, gof_Tijk]=fit(exp_Tijk(:,i),TS_Tijk(:,i),'poly1',fo);
%         [coeffs_Cij, gof_Cij]=fit(exp_Cij(:,i),TS_Cij(:,i),'poly1',fo);
%         [coeffs_J, gof_J]=fit(exp_J(:,i),TS_J(:,i),'poly1',fo);
%         [coeffs_h, gof_h]=fit(exp_h(:,i),TS_h(:,i),'poly1',fo);
        
        
        R2_si(i) = round(gof_si.rsquare,3);
        R2_sij(i) = round(gof_sij.rsquare,3);
        R2_Tijk(i) = round(gof_Tijk.rsquare,3);
        R2_Cij(i) = round(gof_Cij.rsquare,3);
        R2_J(i) = round(gof_J.rsquare,3);
        R2_h(i) = round(gof_h.rsquare,3);
        
        inclination_si(i) = coeffs_si.a;
        inclination_sij(i) = coeffs_sij.a;
        inclination_Tijk(i) = coeffs_Tijk.a;
        inclination_Cij(i) = coeffs_Cij.a;
        inclination_J(i) = coeffs_J.a;
        inclination_h(i) = coeffs_h.a;
        
        
        if (abs(gof_si.rsquare)>=10)
            R2_si(i)=sign(gof_si.rsquare)*inf;
        end
        if (abs(gof_sij.rsquare)>=10)
            R2_sij(i)=sign(gof_sij.rsquare)*inf;
        end
        if (abs(gof_Tijk.rsquare)>=10)
            R2_Tijk(i)=sign(gof_Tijk.rsquare)*inf;
        end
        if (abs(gof_Cij.rsquare)>=10)
            R2_Cij(i)=sign(gof_Cij.rsquare)*inf;
        end
        if (abs(gof_J.rsquare)>=10)
            R2_J(i)=sign(gof_J.rsquare)*inf;
        end
        if (abs(gof_h.rsquare)>=10)
            R2_h(i)=sign(gof_h.rsquare)*inf;
        end
        
        

        fclose all
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        %nexttile
        maxvar = max(max([exp_si,TS_si]));
        minvar = min(min([exp_si,TS_si]));
        A=minvar-(maxvar-minvar)/10:(maxvar-minvar)/10:maxvar+(maxvar-minvar)/10;
        for i=1:res
            plot(exp_si(:,i),TS_si(:,i),marker_list(i), 'MarkerSize',12, 'LineWidth',1.5,'Color',color_list{i})
            hold on
        end
          plot(A,A,'k--', 'LineWidth',2)
                     
        %grid on
        set(gca,'Fontsize',20)
        ylabel("\textit{Predicted} $\mathit{\langle\sigma_i\rangle}$",'FontSize',30,'interpreter', 'latex')
        xlabel("\textit{Simulation} $\mathit{\langle\sigma_i\rangle}$",'FontSize',30,'interpreter', 'latex')
        for i=1:res
            %qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f'),'; R^2=',num2str(R2_si(i),'%.3f')];
            qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f')];
        end
        qmatchar{res+1} = ' y=x';
        legend(qmatchar,'FontSize',20, 'Box', 'off')
        axis 'square'
        set(gca,'linewidth',3)
        ax=gca; ax.XAxis.Exponent = -1; %para colocar notação cientifica
        %xtickangle(30)
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%         plot([-0.04 0.04],[1,1],'r-','LineWidth',2)
%         hold on
%         plot(exp_si(:,q_position),TS_si(:,q_position)./exp_si(:,q_position),'o', 'MarkerSize',8, 'LineWidth',2)
%         axis([-0.04 0.04 0.995 1.005])
%         grid on
%         set(gca,'Fontsize',16)
%         xlabel("\langle\sigma_i\rangle^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         ylabel("\langle\sigma_i\rangle^{TS}/\langle\sigma_i\rangle^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         
%         axis 'square'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        %nexttile
        for i=1:res
            plot(exp_sij(:,i),TS_sij(:,i),marker_list(i), 'MarkerSize',12, 'LineWidth',1.5,'Color',color_list{i})
            hold on
        end
        
        maxvar = max(max([exp_sij,TS_sij]));
        minvar = min(min([exp_sij,TS_sij]));
        A=minvar-(maxvar-minvar)/10:(maxvar-minvar)/10:maxvar+(maxvar-minvar)/10;
        plot(A,A,'k--', 'LineWidth',2)
        %grid on
        set(gca,'Fontsize',20)
        ylabel("\textit{Predicted} $\mathit{\langle\sigma_i\sigma_j\rangle}$",'FontSize',30,'interpreter', 'latex')
        xlabel("\textit{Simulation} $\mathit{\langle\sigma_i\sigma_j\rangle}$",'FontSize',30,'interpreter', 'latex')
        for i=1:res
            %qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f'),'; R^2=',num2str(R2_sij(i),'%.3f')];
            qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f')];
        end
        qmatchar{res+1} = ' y=x';
        legend(qmatchar,'FontSize',20, 'Box', 'off')
        axis 'square'
        set(gca,'linewidth',3)
        ax=gca; ax.XAxis.Exponent = -1; %para colocar notação cientifica
        ax=gca; ax.YAxis.Exponent = -1;
        %xtickangle(30)
        hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        %nexttile
        for i=1:res
            plot(exp_Tijk(:,i),TS_Tijk(:,i),marker_list(i), 'MarkerSize',12, 'LineWidth',1.5,'Color',color_list{i})
            hold on
        end
        

        maxvar = max(max([exp_Tijk,TS_Tijk]));
        minvar = min(min([exp_Tijk,TS_Tijk]));
        A=minvar-(maxvar-minvar)/10:(maxvar-minvar)/10:maxvar+(maxvar-minvar)/10;
        plot(A,A,'k--', 'LineWidth',2)
        %grid on
        set(gca,'Fontsize',20)
        ylabel("\textit{Predicted} $\mathit{\langle\sigma_i\sigma_j\sigma_k\rangle}$",'FontSize',30,'interpreter', 'latex')
        xlabel("\textit{Simulation} $\mathit{\langle\sigma_i\sigma_j\sigma_k\rangle}$",'FontSize',30,'interpreter', 'latex')
        for i=1:res
            qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f'),'; R^2=',num2str(R2_Tijk(i),'%.3f')];
            %qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f')];
        end
        
        qmatchar{res+1} = ' y=x';
        %res=1;
        for i=1:res
            qmatchar{res+i+1} = [' y=',num2str(inclination_Tijk(i),'%.4f'),'x',', for q=',num2str(qmat(i),'%.2f')];
        end
        for i=1:res
            plot(A,inclination_Tijk(i)*A,'--', 'LineWidth',2)
            %hold on
        end
        legend(qmatchar,'FontSize',20, 'Box', 'off')
        axis 'square'
        set(gca,'linewidth',3)
        ax=gca; ax.XAxis.Exponent = -4; %para colocar notação cientifica
        ax=gca; ax.YAxis.Exponent = -4;
        
        %xtickangle(30)
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%         plot([-0.01 0.01],[1,1],'r-','LineWidth',2)
%         hold on
%         plot(exp_Tijk(:,q_position),TS_Tijk(:,q_position)./exp_Tijk(:,4),'o', 'MarkerSize',8, 'LineWidth',2)
%         axis([-0.01 0.01 0.999 1.01])
%         %grid on
%         set(gca,'Fontsize',16)
%         xlabel("\langle\sigma_i\sigma_j\sigma_k\rangle^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         ylabel("\langle\sigma_i\sigma_j\sigma_k\rangle^{TS}/\langle\sigma_i\sigma_j\sigma_k\rangle^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         axis 'square'
        figure
        plot(qmat,var(TS_Tijk),marker_list(i), 'MarkerSize',10, 'LineWidth',2)
        hold on

        
%         figure
%         histogram(TS_Tijk(find(TS_Tijk~=0)))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        figure
        %nexttile
        for i=1:res
            plot(exp_Cij(:,i),TS_Cij(:,i),marker_list(i), 'MarkerSize',12, 'LineWidth',1.5,'Color',color_list{i})
            hold on
        end
        

        maxvar = max(max([exp_Cij,TS_Cij]));
        minvar = min(min([exp_Cij,TS_Cij]));
        A=minvar-(maxvar-minvar)/10:(maxvar-minvar)/10:maxvar+(maxvar-minvar)/10;
        plot(A,A,'k--', 'LineWidth',2)
        %grid on
        set(gca,'Fontsize',20)
        ylabel("\textit{Predicted} $\mathit{C_{ij}}$ ",'FontSize',30,'interpreter', 'latex')
        xlabel("\textit{Simulation} $\mathit{C_{ij}}$ ",'FontSize',30,'interpreter', 'latex')
        for i=1:res
            %qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f'),'; R^2=',num2str(R2_Cij(i),'%.3f')];
            qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f')];
        end
        qmatchar{res+1} = ' y=x';
        legend(qmatchar,'FontSize',20, 'Box', 'off')
        axis 'square'
        set(gca,'linewidth',3)
        
        ax=gca; ax.XAxis.Exponent = -1; %para colocar notação cientifica
        ax=gca; ax.YAxis.Exponent = -1;
        %xtickangle(0)
        %axis([-0.035 0.035 -0.035 0.035])
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%         plot([-0.05 0.05],[1,1],'r-','LineWidth',2)
%         hold on
%         plot(exp_Cij(:,q_position),TS_Cij(:,q_position)./exp_Cij(:,q_position),'o', 'MarkerSize',8, 'LineWidth',2)
%         axis([-0.05 0.05 0.995 1.005])
%         %grid on
%         set(gca,'Fontsize',16)
%         xlabel("C_{ij}^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         ylabel("C_{ij}^{TS}/C_{ij}^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         
%         axis 'square'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure
        %nexttile
        for i=1:res
            plot(0.02*exp_J(:,i),0.02*TS_J(:,i),marker_list(i), 'MarkerSize',12, 'LineWidth',1.5,'Color',color_list{i})
            hold on
        end
        
%         fi  = fopen('J_ferro_n20_q0.8.txt','w');
%         exp_J2 = zeros(20,20);
%         tmp=0;
%         for iii=1:n
%             for iiii=1:n
%                 tmp=tmp+1;
%                 exp_J2(iiii,iii) = exp_J(tmp);
%             end
%         end
%         for iii=1:n
%             for iiii=1:n
%                 fprintf(fi,'%f ',exp_J2(iii,iiii));
%             end
%             fprintf(fi,'\n');
%         end
%         fclose(fi);
        

        maxvar = max(max([0.02*exp_J,0.02*TS_J]));
        minvar = min(min([0.02*exp_J,0.02*TS_J]));
        A=minvar-(maxvar-minvar)/10:(maxvar-minvar)/10:maxvar+(maxvar-minvar)/10;
        plot(A,A,'k--', 'LineWidth',2)
        %grid on
        set(gca,'Fontsize',20)
        ylabel("\textit{Predicted} $\mathit{J_{ij}}$",'FontSize',30,'interpreter', 'latex')
        xlabel("\textit{Simulation} $\mathit{J_{ij}}$",'FontSize',30,'interpreter', 'latex')
        for i=1:res
            %qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f'),'; R^2=',num2str(R2_J(i),'%.3f')];
            qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f')];
        end
        qmatchar{res+1} = ' y=x';
        legend(qmatchar,'FontSize',20, 'Box', 'off')
        axis 'square'
        set(gca,'linewidth',3)
        axis(0.02*[-1.1 1.1 -1.1 1.1]);
        ax=gca; ax.XAxis.Exponent = -2; %para colocar notação cientifica
        ax=gca; ax.YAxis.Exponent = -2;
%         for i=1:res
%             plot(A,inclination_J(i)*A,'k--', 'LineWidth',2)
%             %hold on
%         end
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%         plot([-1.1,1.1],[1,1],'r-','LineWidth',2)
%         hold on
%         plot(exp_J(:,q_position),TS_J(:,q_position)./exp_J(:,q_position),'o', 'MarkerSize',8, 'LineWidth',2)
%         axis([-1.1 1.1 0.5 1.5])
%         %grid on
%         set(gca,'Fontsize',16)
%         xlabel("J_{ij}^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         ylabel("J_{ij}^{TS}/J_{ij}^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         
%         axis 'square'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        figure
        %nexttile
        for i=1:res
            plot(0.02*exp_h(:,i),0.02*TS_h(:,i),marker_list(i), 'MarkerSize',12, 'LineWidth',1.5,'Color',color_list{i})
            hold on
        end
        
%         fi  = fopen('h_ferro_n20_q0.8.txt','w');
%         fprintf(fi,'%f ',exp_h);
%         fclose(fi);
          
        maxvar = max(max([0.02*exp_h,TS_h]));
        minvar = min(min([0.02*exp_h,TS_h]));
        A=minvar-(maxvar-minvar)/10:(maxvar-minvar)/10:maxvar+(maxvar-minvar)/10;
        plot(A,A,'k--', 'LineWidth',2)
        %grid on
        set(gca,'Fontsize',20)
        ylabel("\textit{Predicted} $\mathit{h_{i}}$",'FontSize',30,'interpreter', 'latex')
        xlabel("\textit{Simulation} $\mathit{h_{i}}$",'FontSize',30,'interpreter', 'latex')
        for i=1:res
            %qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f'),'; R^2=',num2str(R2_h(i),'%.3f')];
            qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f')];
        end
        qmatchar{res+1} = ' y=x';
        legend(qmatchar,'FontSize',20, 'Box', 'off')
        axis 'square'
        set(gca,'linewidth',3)
        
        axis(0.02*[-1.1 1.1 -1.1 1.1]);
        ax=gca; ax.XAxis.Exponent = -2; %para colocar notação cientifica
        ax=gca; ax.YAxis.Exponent = -2;
%         for i=1:res
%             plot(A,inclination_h(i)*A,'k--', 'LineWidth',2)
%             %hold on
%         end
        hold off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%         plot([-1.1,1.1],[1,1],'r-','LineWidth',2)
%         hold on
%         %scatter(exp_h(:,4),TS_h(:,4)./exp_h(:,4),colors,'filled','o','MarkerFaceAlpha',0.2)
%         plot(exp_h(:,q_position),TS_h(:,q_position)./exp_h(:,q_position),'o', 'MarkerSize',8, 'LineWidth',2)
%         axis([-1.1 1.1 0.5 1.5]);
%         grid on
%         set(gca,'Fontsize',16)
%         xlabel("h_{i}^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         ylabel("h_{i}^{TS}/h_{i}^{exp}",'FontSize',24,'FontAngle','italic', 'FontName','C059')
%         %set(gca,'Fontsize',20,'FontAngle','italic', 'FontName','Times')
%         %ax=getoptions(ax);
%         %ax.Ticklabel.FontSize = 15;
%         axis 'square'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure
        %nexttile
        for i=1:res
            plot(exp_Cijk(:,i),TS_Cijk(:,i),marker_list(i), 'MarkerSize',12, 'LineWidth',1.5,'Color',color_list{i})
            hold on
        end
        

        maxvar = max(max([exp_Cijk,TS_Cijk]));
        minvar = min(min([exp_Cijk,TS_Cijk]));
        A=minvar-(maxvar-minvar)/10:(maxvar-minvar)/10:maxvar+(maxvar-minvar)/10;
        plot(A,A,'k--', 'LineWidth',2)
        %grid on
        set(gca,'Fontsize',20)
        ylabel("\textit{Predicted} $\mathit{C_{ijk}}$ ",'FontSize',30,'interpreter', 'latex')
        xlabel("\textit{Simulation} $\mathit{C_{ijk}}$ ",'FontSize',30,'interpreter', 'latex')
        for i=1:res
            %qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f'),'; R^2=',num2str(R2_Cij(i),'%.3f')];
            qmatchar{i} =  ['q=',num2str(qmat(i),'%.2f')];
        end
        qmatchar{res+1} = ' y=x';
        legend(qmatchar,'FontSize',20, 'Box', 'off')
        axis 'square'
        set(gca,'linewidth',3)
        
        ax=gca; ax.XAxis.Exponent = -1; %para colocar notação cientifica
        ax=gca; ax.YAxis.Exponent = -1;
        
        
        figure
        plot(qmat,inclination_Tijk,'o-','MarkerSize',12, 'LineWidth',1.5)
        set(gca,'Fontsize',20)
        xlabel("q values (q=0.8)",'FontSize',30)
        ylabel("Inclination of the fit of T_{ijk}",'FontSize',30)
        figure
        plot(qmat,inclination_J,'o-','MarkerSize',12, 'LineWidth',1.5)
        xlabel("q values (q=0.8)",'FontSize',30)
        ylabel("Inclination of the fit of J_{ij}",'FontSize',30)
        figure
        plot(qmat,inclination_h,'o-','MarkerSize',12, 'LineWidth',1.5)
        xlabel("q values (q=0.8)",'FontSize',30)
        ylabel("Inclination of the fit of h_{i}",'FontSize',30)
        
        
        
        
        
        
end