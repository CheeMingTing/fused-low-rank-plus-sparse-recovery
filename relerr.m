function RE = relerr(M,Mhat)
% RE = relerr(M,Mhat) computes relative error (RE) between matrices M and Mhat in Frobenious norm

RE = norm(Mhat-M,'fro')/norm(M,'fro');
end

