%tic;%% optimizing the program
%clear;
for s=1:1
load('test2.mat');
origin=finalmatrix; clear finalmatrix;
scale=length(neurons);
para1=scale;
para2=1;%the number of neurons which are deleted from the network each time;
para3=2;%the smallest number of neurons which are allowed in the network.
gamma=1;%default division by nw method;
para4=0.05;
freq=500;

for i=1:para1
    [Ci, Q]=modularity_und(origin,gamma);% community_louvain(origin,gamma);%
    P=participation_coef(origin,Ci);
big_Q=Q;
big_pc=P;
big_ci=Ci;
   for j=1:freq
                    [Ci1, Q1]= modularity_und(origin,gamma);%community_louvain(origin,gamma); %
                    P1=participation_coef(origin,Ci1);
                  if(Q1>big_Q)
                     big_Q=Q1;
                     big_pc=P1;
                     big_ci=Ci1;
                  end       

   end
   
   if i==1
       save z_celegans_parameter_final big_ci big_Q big_pc;
   end
   
  clear Q P Ci Q1 P1 Ci1;
    
    Pailie=sort(big_pc,'descend');
    Pcri=Pailie(para2);     
    
     if(max(big_pc)<=para4)% output the finaldivision, e.g., lpu. it should be noted that the node with the participation cofficient Pcri has not been deleted from the network.
        for sdm=1:max(big_ci)
         com=find(big_ci==sdm);
         Neu=neurons(com);
         pc=big_pc(com);
         comfinal(sdm)={com};
         Neufinal(sdm)={Neu};
         Pcfinal(sdm)={pc};
        end
        
        for m=1:length(Neufinal)
            if length(Neufinal{m})==1
               Neufinal{m}=[];
            end
        end
        can1=Neufinal;
        can1(cellfun(@isempty,can1))=[];
        mm(s)=length(can1); 
        
        save zcelegan_lpu_final Pcfinal Neufinal comfinal big_Q big_pc neurons can1;
        break;
       % if length(can1)==13
           % break;
        %end
     end 
    
    if Pcri>para4  
        Psave=find(big_pc<Pcri);
        big_ci=big_ci(Psave);
        big_pc=big_pc(Psave);
        neuronsave=neurons(Psave);
        originsave=origin(Psave,Psave);% origin is a variable in the circulation.

        clear neurons origin;
        origin=originsave;neurons=neuronsave;

        Qi_save(i)={big_Q};
        Ci_save(i)={big_ci};
        Pc_save(i)={big_pc}; 
        orin_saveID(i)={originsave};
        neuron_save(i)={neuronsave};
    end
    
    if Pcri<=para4
        search=find(big_pc>para4);
        delete_nodes=neurons(search);
        [x,y]=setdiff(neurons,delete_nodes);
        neurons=neurons(y);
        origin=origin(y,y);
    end     
    clear big_Q big_ci big_pc originsave neuronsave x y search delete_nodes; 
    
    if length(origin)<para3
        break;
    end
                    
end
end

clear i j s scale Pcfinal Neufinal comfinal big_Q big_pc neurons can1 big_ci Ci_save com freq gamma m mm Neu neuron_save origin orin_saveID Pailie para1 para2 para3 para4 pc Pc_save Pcri Psave Qi_save sdm;
%toc;
   
   

