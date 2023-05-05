%clear;
num_Mc=3;%number of subnetwork
num_n=10;% number of nodes in each subnetwork
scalemc=num_Mc*num_n;
ep_nk =4; % expected node degree
pvalue = 0.55;
%%
for dm = 1:400
    
for i=1:num_Mc
    for j=1:num_n
        L{j}=num2str(i*num_n-(j-1));
    end
    sample (i) = {L'};
end
save Lusample sample;
clear i j L ;%sample;
%%

    k=1;
    while k<(num_Mc+1) 
    m_com=zeros(num_n);p=pvalue;

    for i=1:length(m_com)
        for j=1:length(m_com)
            if j>i
                f=rand(1);
                if f<p
                    m_com(i,j)=1;
                end
            end
        end
    end
    Mc=m_com+m_com';
    
    L = weight_conversion(Mc, 'lengths');
    [D,B]=distance_wei(L);
    chucun(k)={Mc};
    dd(k)={D};
    
    ff=find(D==inf);
    if length(ff)>0
        k=k-1;
    end 
 
        k=k+1;
    end
    save dd D;
    clear i j k f D B L;
%% run it n times, n random communities are gained. then put all the communities into a matrix;
%for example, n=3; we get c1, c2,and c3; put these 3 communities in x (trace place ?????)
    x=zeros(scalemc);
    for k=1:length(chucun)
        for i=1:num_n
            for j=1:num_n
                x((k-1)*num_n+i,(k-1)*num_n+j)=chucun{1,k}(i,j);
            end
        end
    end
    clear i j k;
    degree_s=0;
    for i=1:length(chucun)
        [deg] = degrees_und(chucun{1,i});
        %degchucun(i)={[deg]};
        degree_s=mean([deg])+degree_s;
    end
    degree_s=degree_s/num_Mc;% AD average node degree among LUs.
       
    if (degree_s >= ep_nk)  && (degree_s < ep_nk + 0.1)
        
        break;
    end
    
    ddm(dm)=dm;
    clear i;
    
    
end
submodel = x;
G=graph(submodel);
plot(G);
%[kden,N,K] = density_und(submodel);
% calculating the denstiy in each o_LU
for i=1:length(chucun)
    y=chucun{i};
    [kden,N,K] = density_und(y);
    kden_s(i)=kden;
    clear N K y;
end 
  
  zmean_kden=mean(kden_s);
  zref_den=degree_s/(num_n-1);

save subcomm submodel degree_s zmean_kden;
clear chucun ddm deg degree dm ep_nk G i kden kden_s m_com Mc num_Mc num_n p pvalue sample scalemc submodel x zmean_kden zref_den
%%