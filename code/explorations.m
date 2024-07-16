N  = 200;
BB = 9;

for AA = 0.1:0.1:3
%for N = 100:100
  bN = BB*N;
  %W  = eye(bN+1);
  %Wi = eye(bN+1);
  %M  = eye(bN+1);
  %for cc = 0:bN % columns
  %  for rr = 0:(cc-1) % rows
  %    Wi(rr+1,cc+1) = ((N+rr)/N)^(cc-rr) * nchoosek(cc,rr);
  %    W (rr+1,cc+1) = ((N+cc)/N)^(cc-rr) * nchoosek(cc,rr) * (N+rr)/(N+cc) * (-1)^(cc-rr);
  %    for tt = 1:(cc-rr)
  %      M(rr+1,cc+1) = M(rr+1,cc+1) + (-1)^tt * ...
  %        ((N+rr+tt)/N)^(cc-rr-AA*N) * ( N+rr )/(N+rr+tt) * ...
  %        nchoosek(cc-rr,tt);
  %    end
  %    M(rr+1,cc+1) = M(rr+1,cc+1)* factorial(cc)/(factorial(rr)*factorial(cc-rr));
  %  end
  %end
  % alternative way
  T  = zeros(bN+1);
  T(1,1) = 1;
  for i = 1:bN
    T(i+1,i+1) = N/(N+i);
    T(i  ,i+1) = i/(N+i);
  end
  Ma = eye(bN+1);
  for iter = 1:(AA*N)
    Ma = Ma*T;
  end

  %
  counter = 0;
  if ~(counter-10*flor(counter/10))
  figure()
  tiledlayout(1,2)
  nexttile
  idx  = 0:bN;
  xvec = 1-N./(N+idx);
  [xx,yy] = meshgrid(xvec);
  surf(xx,yy,Ma,'EdgeColor','none');
  view(2)
  colorbar
  xlim([0 BB/(1+BB)])
  ylim([0 BB/(1+BB)])
  title(num2str(N))

  %
  %figure()
  nexttile
  expect = xvec*Ma;
  plot(xvec,expect)
  xlim([0 BB/(1+BB)])
  ylim([0 BB/(1+BB)])
  title(num2str(AA))
  end
%end
end

%%

N  = 200;
BB = 5;
maxAA = 4;

bN = BB*N;

alphas = 0:(1/N):maxAA;

T  = zeros(bN+1);
T(1,1) = 1;
for i = 1:bN
  T(i+1,i+1) = N/(N+i);
  T(i  ,i+1) = i/(N+i);
end

idx  = 0:bN;
xvec = 1-N./(N+idx);

GN = zeros(length(alphas),bN+1);
GN(1,:) = xvec;

for i = 2:length(alphas)
  AA = alphas(i);
  disp(AA)
  curr_pow = 0;
  curr_GN  = xvec;
  for k = 0:bN
    while curr_pow < round(AA*(N+k))
      curr_pow = curr_pow + 1;
      curr_GN  = curr_GN*T;
    end
    GN(i,k+1) = curr_GN(k+1);
  end
end

%figure()
%xlim([0 BB/(1+BB)])
%ylim([0 BB/(1+BB)])
%hold on
%for aidx = 2:20:length(alphas)
%  plot(xvec,GN(aidx,:))
%end

figure()
tiledlayout(1,2)
nexttile
xlim([0 BB/(1+BB)])
ylim([0 BB/(1+BB)])
hold on
for num = 0:0.5:maxAA
  [~,qidx] = min(abs(alphas-num));
  plot(xvec,GN(qidx,:))
end
xlabel('Initial density')
ylabel('Average final density')
title('After aN iterations, a = 0, 0.5, 1, ..., 4')


%figure()
nexttile
ylim([0 BB/(1+BB)])
xlim([0 maxAA])
hold on
for num = 0:0.1:1
  [~,qidx] = min(abs(xvec-num));
  plot(alphas,GN(:,qidx))
end
xlabel('a')
ylabel('Average final density')
title('After aN iterations, initial density = 0, 0.1, 0.2, ..., 0.9')

% need to consider a duplication
lambertw(1)/(1+lambertw(1));
lambertw(exp(-1))/(1+lambertw(exp(-1)));

tmp = GN./(1-GN);
FF = tmp.*exp(tmp);

figure()
tiledlayout(1,2)
nexttile
xlim([0 BB/(1+BB)])
%ylim([0 BB/(1+BB)])
hold on
for num = 0:0.5:maxAA
  [~,qidx] = min(abs(alphas-num));
  plot(xvec,log10(FF(qidx,:)))
end
xlabel('Initial density')
ylabel('log(F(.))')
title('After aN iterations, a = 0, 0.5, 1, ..., 4')

nexttile
%figure()
%ylim([0 BB/(1+BB)])
xlim([0 maxAA])
hold on
for num = 0:0.1:1
  [~,qidx] = min(abs(xvec-num));
  plot(alphas,log10(FF(:,qidx)))
end
xlabel('a')
ylabel('log(F(.))')
title('After aN iterations, initial density = 0, 0.1, 0.2, ..., 0.9')

%%
LIN = log(FF);

%figure()
%[xx,yy] = meshgrid(xvec,alphas);
%surf(xx,yy,LIN,'EdgeColor','none');
%view(2)
%colorbar
%xlim([0 BB/(1+BB)])
%ylim([0 maxAA])
%  title(num2str(N))

%a = diff(LIN,1,1);
%figure()
%imagesc(a)
%colorbar

%figure()
%plot( xvec, mean(a,1) )

%%
k1 = zeros(1,length(xvec));
k2 = zeros(1,length(xvec));

for idxp = 1:length(xvec)
  p = polyfit(alphas, LIN(:,idxp), 1);
  k1(idxp) = p(1);
  k2(idxp) = p(2);
end

figure()
tiledlayout(2,1)
nexttile
plot(xvec,k1)
xlabel('Initial density')
ylabel('k1(init dens)')
xlim([0 BB/(1+BB)])

nexttile
plot(xvec,k2)
xlabel('Initial density')
ylabel('k2(init dens)')
xlim([0 BB/(1+BB)])

figure()
tiledlayout(2,1)
nexttile
plot(k1)
xlabel('Index')
ylabel('k1(init dens)')
%xlim([0 BB/(1+BB)])

nexttile
plot(k2)
xlabel('Index')
ylabel('k2(init dens)')
%xlim([0 BB/(1+BB)])


idx = 1:length(xvec);
k1(1) = -1;
polyfit(idx,k1,1)

figure()
plot(idx,k1)
hold on
plot(idx,-1-idx/N)

figure()
plot(xvec,k1)
hold on
plot(xvec,-1./(1-xvec))

figure()
plot((-k1-1./(1-xvec)))

figure()
plot(idx,exp(k2))

k2(1) = -6;
polyfit(xvec-0.5,k2,6)

figure()
plot(xvec,k2)
hold on
plot(xvec,3*atanh(2*xvec-1)+1)

figure()
plot(xvec, (atanh(2*xvec-1)+1)-k2 )

figure()
plot(k2)

figure()
plot(log10(exp(k2)))

figure()
plot(xvec,log10(exp(k2)))

figure()
plot(lambertw(exp(k2)))

figure()
plot(xvec,k2)
hold on
tmp  = xvec./(1-xvec);
tmp2 = tmp.*exp(tmp);
plot(xvec,log(tmp2))

%%
LINa = zeros(length(alphas),bN+1);
for idxp = 1:length(xvec)
  LINa(:,idxp) = k1(idxp)*alphas;
end

figure()
[xx,yy] = meshgrid(xvec,alphas);
surf(xx,yy,LIN,'EdgeColor','none');
view(2)
xlim([0 BB/(1+BB)])
ylim([0 maxAA])
colorbar
title('log(FF)')

figure()
surf(xx,yy,LINa,'EdgeColor','none');
view(2)
xlim([0 BB/(1+BB)])
ylim([0 maxAA])
colorbar
title('linear part')

figure()
surf(xx,yy,exp(LIN-LINa),'EdgeColor','none');
view(2)
xlim([0 BB/(1+BB)])
ylim([0 maxAA])
colorbar
title('log(FF)/linear part = g(BB)')

%%
polyfit(xvec,k1,2)