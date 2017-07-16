function RMSE = compRMSE(predv, realv)
    t=max(size(predv));
    tpredv=reshape(predv, t,1);
    trealv=reshape(realv, t,1);
    RMSE= sqrt(mean((tpredv-trealv).^2 ));
end

