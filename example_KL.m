
% % how to generate 1D stochastic process with mean zero, 
% % exponential covariance kernal with correlation length lc 
% 
% % 1D domain 
% FinalTime = 1; 
% dt = 0.01; 
% timegrid = 0:dt:FinalTime; 
% 
% % correlation length of the exponential kernel
% % lc closer to zero means that it is closer to white noise 
% lc = 0.1; 
% 
% % Generate Karhunen Loeve expansion of stochastic process 
% [nphi, phi, lambda, rDim] = get_KL_coeff(lc, timegrid);
% 
% sig = 1; % standard deviation 
% 
% % generate samples 
% figure; hold on; 
% for n = 1:3 % generate 3 sample process 
%     noise = sig*nphi*randn(rDim,1); % multiply with normal random variable sample 
%     plot( timegrid, noise ) 
% end 
% 
% %%
% % how to generate 2D stochastic field with mean zero, 
% % exponential covariance kernal with correlation length lc 
% 
% % 2D domain - use 1D^2 
% 
% % correlation length in each axis 
% lc = [0.5, 0.1]; 
% 
% % Generate Karhunen Loeve expansion of stochastic process 
% [nphi1, ~, ~, rDim(1)] = get_KL_coeff(lc(1), timegrid);
% [nphi2, ~, ~, rDim(2)] = get_KL_coeff(lc(2), timegrid);
% 
% sig = 1; % standard deviation 
% noise = zeros( length(timegrid) ); 
% 
% for n = 1:rDim(1)
%     for m = 1:rDim(2) 
%         noise = noise + sig* nphi1(:,n) * nphi2(:,m)' * randn(1); 
%     end 
% end 
% 
% figure; imagesc(timegrid, timegrid, noise' ) 


% 3D domain - use 1D^3 
function noise = example_KL

FinalTime1 = 1.99;
FinalTime2 = 2.79;
FinalTime3 = 1.27;
dt = 0.01; 
timegrid1 = 0:dt:FinalTime1; 
timegrid2 = 0:dt:FinalTime2; 
timegrid3 = 0:dt:FinalTime3; 


% correlation length in each axis 
lc = [0.03, 0.03, 0.03]; 

% Generate Karhunen Loeve expansion of stochastic process 
[nphi1, ~, ~, rDim(1)] = get_KL_coeff(lc(1), timegrid1);
[nphi2, ~, ~, rDim(2)] = get_KL_coeff(lc(2), timegrid2);
[nphi3, ~, ~, rDim(3)] = get_KL_coeff(lc(3), timegrid3);

sig = 411.565; % standard deviation 
noise = zeros( length(timegrid1), length(timegrid2), length(timegrid3) ); 

for n = 1:rDim(1)
    for m = 1:rDim(2) 
        for k = 1:rDim(3) 
            rsample = randn(1); 
            noise2D = nphi1(:,n) * nphi2(:,m)'; 
            for i = 1:size(noise,3) 
              noise(:,:,i) = noise(:,:,i) + sig* noise2D * nphi3(i,k) * rsample; 
            end
        end 
    end 
end 

figure; imagesc(timegrid1, timegrid2, noise(:,:,50) ) 

end



