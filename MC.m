function[]=MC(varargin)
sample = varargin{1}; %for one number of samples
if nargin==2
    testes = varargin{2};
else
    testes = 1;
end

n = 40;
beta = 0.001;
q = 0.8;
opt = '+-';
dist = 1;
[J,h]=Jh_gen(n,3,1);
J =J/2;
allmeans = zeros(testes,n);
allmeansmeans = zeros(1,n);






%parpool;

if isempty(gcp('nocreate'))
    maxNumCompThreads(8);
    parpool;
    ThrNumber = maxNumCompThreads;
else
    maxNumCompThreads(8);
    ThrNumber = maxNumCompThreads;
end

globalMCmeansSi = zeros(ThrNumber,n);
for kkk=1:ThrNumber
    globalMCmeansSij{kkk} = zeros(n,n);
end

errors = zeros(1,ThrNumber);
times = zeros(1,ThrNumber);
time = 0;
error = 0;
states = zeros(ThrNumber,sample);

%states{:}(1) = bin2dec(num2str((sample1+1)/2));

%% Exact Mean
if n<=20
    
    
    sample1 = (dec2bin(randi(2^n)-1,n)-'0')';
    sample1 = sample1'*2-1;
    % sample1 = ones(1,n);
    sample2 = sample1;
    P=zeros(1,2^n);
    satatesexact = zeros(1,2^n);
    
    
    
    tic
    Zq = qZpart(beta*J,beta*h,n,q,dist,opt);

    exactmean = zeros(1,n);
    si = (dec2bin(0:2^n-1)-'0')';  %si no intervalo [0,1]
    [~,m]=size(si);
    si = 2*si-1; %si no intervalo [-1,1]
    meanmean = zeros(m,n);
    meansij = zeros(n,n,m);
    %meansijk = zeros(m,n,n,n);

    % isingmat = zeros(1,2^n);
    % pmat = zeros(1,2^n);
    % statesmat = zeros(1,2^n);

    for k = 1:m
        si_tmp = si(:,k);
        %statesmat(k) = bin2dec(num2str((si_tmp'+1)/2));
        isingtmp = ising(J,h,si_tmp',opt,q);
        statesexact(k) = isingtmp;
    %     isingmat(k) = isingtmp;
    %     pmat(k) = qexp(q,-beta*isingtmp,dist)/(Zq);
        meanmean(k,:) = (si_tmp'-sign((1-q)*h))*qexp(q,-beta*isingtmp,dist)/(Zq);
%         if k==1024
%             disp(meanmean(k,:))
%         end
        meansij(:,:,k) = (sij_gen(si_tmp')-sign((1-q)*J))*qexp(q,-beta*isingtmp,dist)/(Zq);   
        P(k) = qexp(q,-beta*isingtmp,dist)/(Zq);
    end
    % figure
    % plot(isingmat,pmat,'*')
    % figure
    % plot(statesmat,sort(isingmat),'s')
    exactmean = sum(meanmean)
    exactmeansij = sum(meansij,3);
    toc
else
    exactmean = zeros(1,n);
    exactmeansij = zeros(n,n);
end



%% Monte Carlos mean analise
for jj=1:testes
%     disp("J = ")
%     disp(J)
%     disp("h = ")
%     disp(h)
    meanMC = zeros(1,n);
    meansijMC = sij_gen(zeros(1,n));
    e1 = 0;
    e2 = 0;
    isingtmp = 0;
    parfor kkk = 1:ThrNumber
        Nmc = 1;
        %sample1 = (dec2bin(randi(2^n)-1,n)-'0')';
        sample1 = (dec2bin(0,n)-'0')';
        sample1 = sample1'*2-1;
        %states{kk}(1) = bin2dec(num2str((sample1+1)/2));
        tic
        for i=2:sample
            ii = randi(n,1);
            sample2 = sample1;
            sample2(ii) = -sample2(ii);%flip a spin in a random location


            isingtmp1 = ising(J,h,sample1,opt,q); %energy of the original spin config
            states(kkk,i) = isingtmp1;
            e1 = qexp(q,-beta*isingtmp1,dist); 

            isingtmp2 = ising(J,h,sample2,opt,q); %energy with fliped spin
            e2 = qexp(q,-beta*isingtmp2,dist); 

            %Q1 =  qexp(q,-beta*(isingtmp2-isingtmp1),dist);%hasting therm 1
            %Q2 =  qexp(q,-beta*(isingtmp1-isingtmp2),dist);%hasting therm 2
            Q1=1;Q2=1;

            if((e2*Q1)/(e1*Q2)>=rand(1))
            %if(exp(isingtmp1-isingtmp2)>=rand(1))
                sample1 = sample2;
            end

            %states(kkk,kk) = isingtmp1;
            %states{kk}(i) = bin2dec(num2str((sample1+1)/2));
            globalMCmeansSi(kkk,:) = (globalMCmeansSi(kkk,:) + (sample1-sign((1-q)*h)));
            globalMCmeansSij{kkk} = globalMCmeansSij{kkk} + (sij_gen(sample1)-sign((1-q)*J));
            Nmc = Nmc +1;

        end
        globalMCmeansSi(kkk,:) = globalMCmeansSi(kkk,:)/Nmc
        globalMCmeansSij{kkk} = globalMCmeansSij{kkk}/Nmc;

        times(kkk) = toc;

    end

    for kkk=1:ThrNumber
    meanMC = meanMC + globalMCmeansSi(kkk,:);
    meansijMC = meansijMC + globalMCmeansSij{kkk};
    %error = error + errors(kkk);
    time = time + times(kkk);
    end
    meanMC = meanMC/ThrNumber
    allmeans(jj,:) = meanMC;
    allmeansmeans = allmeansmeans+allmeans(jj,:);
    meansijMC = meansijMC/ThrNumber;
    time= time/ThrNumber;
    error = sum(abs(exactmean-meanMC))/n;

end

% for jj=1:testes
%     plot(allmeans(jj,:),'*-'), hold on
% end
plot(allmeansmeans/testes,'o-','LineWidth',4)

for jjj=1:n
    yneg(jjj) = allmeansmeans(jjj)/testes-min(allmeans(:,jjj));
end
for jjj=1:n
    ypos(jjj) = max(allmeans(:,jjj))-allmeansmeans(jjj)/testes;
end

errorbar(1:n,allmeansmeans/testes,yneg,ypos,'CapSize',18,'LineWidth',2)


%save('samp_err.mat','sample','error')


%figure

% subplot(1,2,1)
% plot(sample,error,'-*'),legend(['q= ',num2str(q)]),xlabel('Number of samples'),ylabel('Absolute error of the means'),ax = gca;ax.FontSize = 20; 
% 
% subplot(1,2,2)
% plot(sample,time,'-*'),legend(['q= ',num2str(q)],'location','southeast'),xlabel('Number of samples'),ylabel('Time of execution (s)'),ax = gca;ax.FontSize = 20; 


% figure
% %histogram(states{kkk})
% h=histogram(states*beta);
% h.Normalization = 'countdensity';
% set(gca,'Fontsize',20)
% title("Monte Carlo calculation of the energy density for a 20 spin system (q=0.8)")
% xlabel("E",'FontSize',30,"interpreter","latex")
% ylabel("$\rho(E)$",'FontSize',30,"interpreter","latex")
% 
% 
% 
% figure
% histogram(statesexact*beta)
% set(gca,'Fontsize',20)
% title("Exact calculation of the energy density for a 20 spin system (q=0.8)")
% xlabel("E",'FontSize',30,"interpreter","latex")
% ylabel("$\rho(E)$",'FontSize',30,"interpreter","latex")
% %plot(0:2^n-1,P*samples(end),"LineWidth",8)







