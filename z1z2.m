clear;
for dms=1:100  
    
z1_form_lu_general;
z2_Creating_network_general;
sum_dms(dms)=mean_simiLU;
sum_dmsq(dms)=chuQ;
%degree_s should also be saved as a file.
end
amean_simi=mean(sum_dms);
ameanQ=mean(sum_dmsq);



