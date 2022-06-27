function [concentration, cgradX, cgradY, cgradZ, cancer_sd] = set_cancer_run_tumor(coh_map) 
    global cancer_center cancer_size num_cancers curr_site_num
    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));

    
%  If multiple cancer sites exist, we assume they are all the same size
    num_cancers = size(cancer_center,1);

% concentration = 3.146846451532616e+03.* exp(-((X - 79).^2)./(2.*61.438725523798006.^2) - ((Y - 84).^2)./(2.*61.438725523798006.^2) - ((Z - 78).^2)./(2.*61.438725523798006.^2)) + ... 
%                  1.639599793346383e+03.* exp(-((X - 139).^2)./(2.*14.085786226612475.^2) - ((Y - 93).^2)./(2.*14.085786226612475.^2) - ((Z - 78).^2)./(2.*14.085786226612475.^2));



% Computing Variance for Constant Case 
    avg_conc =  4059.635
    exp_con = [3672.8, 4303, 4511.5, 3751.24]

    var_constant = ((exp_con[0]-avg_conc).^2 + (exp_con[1]-avg_conc).^2 +	(exp_con[2] - avg_conc).^2	+ ...
                    (exp_con[3] - avg_conc).^2)/3;


    concentration_1 = avg_conc.*ones(size(coh_map,1), size(coh_map,2), size(coh_map,3));

    noise = example_KL()*var_constant;

    concentration = concentration_1 + noise;


% Computing Variance for Tumor Case

    avg_tumor = 
    sd_tumor = 
    exp_conc = [3672.8, 4303, 4511.5, 3751.24]

    concentration = avg_tumor.* exp(-((avg_tumor - 79).^2)./(2.*sd_tumor.^2) - ((avg_tumor - 84).^2)./(2.*sd_tumor.^2) - ((avg_tumor - 78).^2)./(2.*sd_tumor.^2));

    var_tumor = ((exp_conc[0]-avg_conc).^2 + (exp_conc[1]-avg_conc).^2 + (exp_conc[2] - avg_conc).^2 + ...
                (exp_conc[3] - avg_conc).^2)/3;


% Computing Variance for Cyt Max Case
    data.xdata = [79 101 139 101]; 
    data.ydata = [84 95 93 119]; 
    data.zdata = [78 78 78 78]; 


    avg_tumor = 
    sd_tumor = 
    loc_conc = [3672.8, 4303, 4511.5, 3751.24]

    [~, nMaxCytkn] = max( data.cytokine ); 
    cyt_max = [data.xdata(nMaxCytkn), data.ydata(nMaxCytkn), data.zdata(nMaxCytkn)]; 

    concentration = avg_tumor.* exp(-((avg_tumor - cyt_max[0]).^2)./(2.*sd_tumor.^2) - ((avg_tumor - cyt_max[1]).^2)./(2.*sd_tumor.^2) - ((avg_tumor - cyt_max[2]).^2)./(2.*sd_tumor.^2));

    var_tumor = ((loc_conc[0]-avg_conc).^2 + (loc_conc[1]-avg_conc).^2 + (loc_conc[2] - avg_conc).^2 + ...
                (loc_conc[3] - avg_conc).^2)/3;





    
    concentration( coh_map==0 ) = 0; 
    
    cgradX = zeros(size(X));
    cgradX(2:end-1,:,:) = concentration(3:end,:,:) - concentration(1:end-2,:,:);

    cgradY = zeros(size(X));
    cgradY(:,2:end-1,:) = concentration(:,3:end,:) - concentration(:,1:end-2,:);

    cgradZ = zeros(size(X));
    cgradZ(:,:,2:end-1) = concentration(:,:,3:end) - concentration(:,:,1:end-2);
    
    % find largest vector magnitude
    L = (cgradX.^2 + cgradY.^2 + cgradZ.^2).^(.5);
    maxL = max(max(max(L)));

    % normalize gradient
    cgradX = cgradX./maxL;
    cgradY = cgradY./maxL;
    cgradZ = cgradZ./maxL;
end
