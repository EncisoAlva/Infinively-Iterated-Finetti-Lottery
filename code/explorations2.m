maxN = 200;
maxk = 30;

betas = 0:0.1:5;
ens   = maxk:maxN;

collector = ones( length(betas), maxk+1, length(ens) );

for i = 1:length(betas)
  BB = betas(i);
  for j = 1:length(ens)
    N = ens(j);
    for k = 0:maxk
      for tt = 0:(k-1)
        collector(i,k+1,j) = collector(i,k+1,j) * ...
          (BB*N-tt)/(N+k);
      end
    end
  end
end

for j = 1:length(betas)
  figure()
  xlim([0, maxk])
  hold on
  for NN = 10:10:maxN
    [~, idi] = min( abs(ens-NN) );
    plot( (0:maxk), log(collector( j, :, idi) ))
  end
end