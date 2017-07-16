function loss = computeLoss(RT_theta_by_pi, cens_theta_by_pi, cutoff, cutoff_penalty_factor)
RT_theta_by_pi = squeeze(RT_theta_by_pi);
cens_theta_by_pi = squeeze(cens_theta_by_pi);
for i=1:size(RT_theta_by_pi,1)
    if(all(cens_theta_by_pi(i,:) == 0))
        RT_theta_by_pi(i,find(RT_theta_by_pi(i,:) > cutoff)) = cutoff*cutoff_penalty_factor;
    else
        y = RT_theta_by_pi(i,find(cens_theta_by_pi(i,:) == 0));
        c = RT_theta_by_pi(i,find(cens_theta_by_pi(i,:) == 1));
        numsamples = size(RT_theta_by_pi,2);
        RT_theta_by_pi(i,:) =  fitDistAndReturnMeanLoss(y, c, cutoff, cutoff_penalty_factor, numsamples);
    end
    loss = mean(RT_theta_by_pi, 2);
end