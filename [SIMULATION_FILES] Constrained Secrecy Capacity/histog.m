close all; clear; clc
for file_Index=1:100
    load(int2str(file_Index))
%     MAXindex=find(Capacity_iid==max(Capacity_iid));
    sec_cap_sample(file_Index)=max(C(:,1));
    plot(C(:,1));
        hold on;
%     sec_cap_sample_vec(file_Index,:)=collector_samples;
%     plot(C(:,1));
%     hold on;
%     if (min(min(min(log2(1./TRANSTAR_vec(181,:,:)))))>=1)
%         plot(C(:,1));
%         hold on;
%         [file_Index min(min(min(log2(1./TRANSTAR_vec(181,:,:))))) sec_cap_sample(file_Index)] 
%     end
end