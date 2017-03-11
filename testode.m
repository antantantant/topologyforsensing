function dydt = testode(t,y, K, M, F)

dydt = [zeros(size(K)), -inv(M)*K; eye(size(K)), zeros(size(K))]*y + [inv(M)*F; zeros(size(K,1),1)];