%clear;
tic;
load('subcomm.mat');% MC is the sample communities
load('Lusample.mat', 'sample');
am=18;%floor(degree_s);
pnn=20;
nnode=length(sample{1});
nlu=length(sample);
netscale=pnn+nlu*nnode;% the complete size of the created network
sdata=nlu*nnode+1;tdata=netscale; % for example, we have 30 nodes already forming three communities. so the projection neuron will start from 31th;
max_link=floor(degree_s);%the maximum number of connections for each projection neuron is no bigger than degree_s (connections inside LUs).
fre=1; % circulation times.
%sum_simiLU=0;

for sm=1:fre

    A=zeros(netscale);
for i=1:length(submodel)
    for j=1:length(submodel)
        A(i,j)=submodel(i,j);
    end
end
clear i j;
%%

for i=sdata:tdata

    j=1;
    f=randi([2,am]); % Generating 2 to 4 links for each projection neuron.
    while (j<(f+1))% the number of links depends on the value of f;
        
        f1=randi([1,nlu*nnode]);

        for sam=1:length(sample)
            same=intersect(num2str(f1),sample{sam});
            if isempty(same)~=1
                sd=sam;break;
            end
        end
        clear sam same;
                
        m(j)=sd;% indicate they are in different communities.
        mm(j)=f1;
        
        A(f1,i) = 1; % add one link one LU randomly; the data 1 indicating one new link is added.
        A(i,f1) = A(f1,i);
        
        for n= 1: j-1
                           
                if mm(j)== mm(n)

                   j= j-1;
                   break;   

                end        
                
        end
        
        if j==f
            msigh = 0;
            for k=1:f-1
                if m(k)== m(f)
                    msigh=msigh+1;
                end
            end
            
            if msigh == f-1
                j = j-1;
                A(i,f1) = 0; A(f1,i)= 0;
            end 
        end 
                       
        j=j+1;
    end
    
    save mei j m mm;
    clear m n mm;
    
    
end

clear i j k msigh f f1 sd ;
 
[deg] = degrees_und(A);% [deg] calculate the degree for each node in th network A;
s=0;
for i=sdata:tdata
    
    s = s+ sum(deg(i));
end
APC= s / nlu;%(tdata-sdata+1);
lamda= APC / degree_s;
[Ci, Q]=community_louvain(A,1);

chuQ(sm)=Q;
chuLAMDA(sm)=lamda;
g=graph(A);
plot(g);
 
save create_network degree_s chuLAMDA chuQ s;

clear i APC Ci Q g lamda s ; 
%% Dectecting LU
for i=1:netscale

    neu{i}=num2str(i);
end
neurons=neu';
finalmatrix=A;
save test2 finalmatrix  neurons ;
%[Ci1,Q1]=modularity_und(finalmatrix,1);% Q in neuman definition is the
%link density in a community.
%[kdenf,N,K] = density_und(finalmatrix);%network density
%smkdenf(sm)=Q1/kdenf;% is not useful rightnow.
%comdens(sm)=Q1;
%save test2 finalmatrix neurons kdenf;
%%
zCommunity_detect_releaseconditon;
%clear;
N_LU=nlu; %number of original LUs
S_LU=nnode;%number of nodes in each LU.
load('Lusample.mat', 'sample');
load('zcelegan_lpu_final.mat','can1');

    for j=1: length(sample)
        
        Si=0;
        for i=1:length(can1)
          
            simi=length(intersect(sample{j},can1{i}))/length(union(sample{j},can1{i}));
            
            if simi > Si
                
                Si = simi;
            end
          % err(j)=length(setdiff(can1{j},sample{i}))/S_LU;                   
        end
         
        SSd(j)=Si;
        clear i; 
        
    end
      
        clear j;
        simi_LU= mean(SSd);
        sum_simiLU (sm)=simi_LU;
 
if length(can1)==N_LU
    vic_sign(sm)=1;
end
        
if length(can1)~=N_LU % lu is broken
     vic_sign(sm)=0;
end 
clear N_LU S_LU ;
%save LUerro err_LU err vic_sign;
end

toc;
t=toc;
 
    mean_simiLU= mean(sum_simiLU);
    save data SSd can1 mean_simiLU t vic_sign sample sum_simiLU;%smkdenf  comdens
    %clear;
    load('data.mat');
    load('create_network.mat');
    load('mei.mat');
    %load('test2.mat');
    load('subcomm.mat','zmean_kden');
    
    meanQ=mean(chuQ);
    meanLamda=mean(chuLAMDA);
    meanvic=mean(vic_sign);
    %rationkdenf=mean(smkdenf);
    %meancomdens=mean(comdens);
    %
    clear s sdata Si simi sm;


    

