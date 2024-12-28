sigma_r = [50,100,300,100,300,50,1000];
sigma_theta = [0.5,1.5,2,1.5,2,0.5,2.5]*pi/180;
sigma_theta_names = ["0p5","1p5","2","1p5","2","0p5","2p5"];

sigma_rbar = [100,100,100,50,300,1000,50];
sigma_thetabar = [1.5,1.5,1.5,0.5,2,2.5,0.5]*pi/180;
sigma_thetabar_names = ["1p5","1p5","1p5","0p5","2","2p5","0p5"];

sigma_x = 100;

Nmc = 1e6;

for i = 1:length(sigma_r)
    for j = 1:10
        [LB,MSE] = IEKS_MCRLB_main(Nmc,sigma_r(i),sigma_theta(i),sigma_rbar(i),sigma_thetabar(i),sigma_x);
    
        save(['monte carlo results/Nmc_',num2str(Nmc),'_sigma_r_',num2str(sigma_r(i)),'_sigma_theta_',char(sigma_theta_names(i)),'_sigma_rbar_',num2str(sigma_rbar(i)),'_sigma_thetabar_',char(sigma_thetabar_names(i)),'_LB','trial',num2str(j),'.mat'],"LB");
        save(['monte carlo results/Nmc_',num2str(Nmc),'_sigma_r_',num2str(sigma_r(i)),'_sigma_theta_',char(sigma_theta_names(i)),'_sigma_rbar_',num2str(sigma_rbar(i)),'_sigma_thetabar_',char(sigma_thetabar_names(i)),'_MSE','trial',num2str(j),'.mat'],"MSE");
    end
end

