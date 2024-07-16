maxN = 20;

col_evec   = cell(maxN,1);
col_eveci  = cell(maxN,1);
col_evec2  = cell(maxN,1);

col_eig    = cell(maxN,1);

collect1 = zeros(maxN,1);
collect2 = zeros(maxN,1);

collect_new = zeros(maxN,1);
for N = 1:maxN
  E0 = zeros(N+1,1);
  v  = zeros(N+1,1);
  M  = zeros(N+1,N+1);
  %
  E0(N+1) = 1;
  for k = 1:N
    v(k+1) = k/(k+N);
  end
  M(1,1) = 1;
  for k = 1:N
    M(k+1,k+1) = N/(k+N);
    M(k,  k+1) = k/(k+N);
  end
  collect1(N) = v'*(M^N)*E0;
  collect2(N) = v'*(M^(2*N))*E0;
  idx = linspace(0,1,N+1);
  [evec,evv] = eig(M );
  evec = evec*diag(diag(evec).^-1);
  col_evec{N}  = evec;
  ivec = inv(evec);
  %col_eveci{N} = inv(evec);

  evec2 = zeros(N+1,N+1);
  evec2(1,1) = 1;
  if N >= 2
    evec2(1,2) = -1;
    evec2(2,2) =  1;
  end
  for k = 2:N
    la = (k+N)/N;
    for j = 0:(k-1)
      evec2(j+1,k+1) = -nchoosek(k-1,j)*( (-la)^(k-1-j) );
    end
    for j = 1:k
      evec2(j+1,k+1) = evec2(j+1,k+1) ...
                     + nchoosek(k-1,j-1)*( (-la)^(k-j) );
    end
  end
  col_evec2{N}  = evec2;
  col_eveci{N} = inv(evec2);

  ivec2 = zeros(N+1,N+1);
  ivec2(1,1) = 1;
  if N >= 2
    ivec2(1,2) =  1;
    ivec2(2,2) =  1;
  end
  for k = 2:N
    la = (k+N)/N;
    ivec2(1,  k+1) = 1;
    ivec2(k+1,k+1) = 1;
    for j = 1:(k-1)
      %ivec2(j+1,k+1) = nchoosek(k-1,j)*((la)^(j+1)) + ...
      %nchoosek(k-1,j-1)*((la)^(k-j));
      ivec2(j+1,k+1) = nchoosek(k-1,j-1)*((la)^(k-j));
      ivec2(j+1,k+1) = ivec2(j+1,k+1) + nchoosek(k-1,j)*((la)^(j+1));
    end
  end

  vec1 = v'*evec;
  vec2 = zeros(N+1,1);
  for k = 1:N
    vec2(k+1) = -( (-1)^k )*( (k/N)^(k-1) )*(k/(k+N));
  end
  
  vec1A = v'*evec*(evv^N);
  vec2A = zeros(N+1,1);
  for k = 1:N
    vec2A(k+1) = -( (-1)^k )*( (k/N)^(k-1) )*(k/(k+N))*( (N/(k+N))^(N) );
    vec2A(k+1) = -( (-1)^k )*( (k^k*N^N)/( N^(k-1)*(N+k)^(N+1) ) );
  end

  col_eig{N} = diag(evv);
end

for N=1:maxN
  tmp = 0;
  for k = 0:N
    tmp = tmp - nchoosek(N,k)*(-k/(N+k))^k * (N/(k+N));
  end
  collect_new(k) = -tmp;
end

plot(collect1,'o')
hold on
plot(collect_new,'o')
yline(lambertw(1)/(1+lambertw(1)))
hold off

%%
qq = [];
ll = [];
for N = 5:maxN
  %vec1 = (sum( col_eveci{N}, 1 )); %- log(col_eig{N}')
  la = (5+N)/N;
  %q = vec1(4)-2-2*la3;
  tmp = col_eveci{N};
  q = tmp(3,6);
  %disp(q)
  %
  qq(end+1) = q;
  ll(end+1) = la;
end

mdl1 = fitlm([ll'], qq');

mdl2 = fitlm([ll', ll'.^2], qq');

mdl3 = fitlm([ll', ll'.^2, ll'.^3], qq');

mdl4 = fitlm([ll', ll'.^2, ll'.^3, ll'.^4], qq');

mdl5 = fitlm([ll', ll'.^2, ll'.^3, ll'.^4, ll'.^5], qq');

coeff = mdl3.Coefficients.Estimate;


for N = 1:maxN
  disp(col_eveci{N});
  t4 = (4+N)/N;
  vec1 = [1,27/16+t4*27/16+ t4^2*9/16+ t4^3*1/16,  6/4+t4*12/4+t4^2*6/4, 1+t4*3, 1];
  disp(vec1)
end

%%

for N = 2:maxN
  tmp = (col_eveci{N});  
  disp(tmp(:,N+1))
  %disp((N+1)^(N-1)/N^(N-2))
  disp( ( (N+2)/N )^(N-2) * N*(N-1)/2 )
end

%%

maxN = 20;

col_evec   = cell(maxN,1);
col_eveci  = cell(maxN,1);
col_evec2  = cell(maxN,1);

col_eig    = cell(maxN,1);

collect1 = zeros(maxN,1);
collect2 = zeros(maxN,1);

for N = 1:maxN
  E0 = zeros(N+1,1);
  v  = zeros(N+1,1);
  M  = zeros(N+1,N+1);
  %
  E0(N+1) = 1;
  for k = 1:N
    v(k+1) = k/(k+N);
  end
  M(1,1) = 1;
  for k = 1:N
    M(k+1,k+1) = N/(k+N);
    M(k,  k+1) = k/(k+N);
  end
  collect1(N) = v'*(M^N)*E0;
  collect2(N) = v'*(M^(2*N))*E0;
  %
  [evec,evv] = eig(M );
  evec = evec*diag(diag(evec).^-1); % scalar multiples of eigenvalues
  col_evec{N}  = evec;
  ivec = inv(evec);

  vec1 = v'*evec;
  vec2 = zeros(N+1,1);
  for k = 1:N
    vec2(k+1) = -( (-1)^k )*( (k/N)^(k-1) )*(k/(k+N));
  end
  disp(vec1./vec2')
  
  col_eig{N} = diag(evv);
end

collect_new = zeros(maxN,1);
for N=1:maxN
  tmp = 0;
  for k = 1:N
    tmp = tmp - nchoosek(N,k)*(-k/(N+k))^k * (N/(k+N));
  end
  collect_new(k) = tmp;
end

collect_new_v2 = zeros(maxN,1);
for N=1:maxN
  tmp = 0;
  for k = 1:N
    tmp = tmp + (-k)^(k-1)/factorial(k-1);
  end
  collect_new_v2(k) = tmp;
end

collect_new_v3 = zeros(maxN,1);
for N=1:maxN
  tmp = 0;
  for k = 1:N
    tmp = tmp + (-1)^(k-1)/factorial(k-1);
  end
  collect_new_v3(k) = tmp;
end

plot(collect1,'.')
hold on
plot(collect_new,'.')
plot(collect_new_v3,'.')
%plot(collect_new_v2,'.')
yline(lambertw(1)/(1+lambertw(1)))
hold off
legend('Matrix','Analytic 1','Taylor of W','Analytic 2')

collect_new2 = zeros(maxN,1);
for N=1:maxN
  tmp = 0;
  for k = 1:N
    tmp = tmp - nchoosek(N,k)*(-k/(N+k))^k * (N/(k+N))^N * (N/(k+N));
  end
  collect_new2(k) = tmp;
end

figure()
plot(collect2,'o')
hold on
plot(collect_new2,'o')
yline(lambertw(exp(-1))/(1+lambertw(exp(-1))))
hold off

%%

maxN = 100;
minN = 10;

nNs = size(minN:maxN,2);

conver = zeros(nNs,nNs);
for N = minN:maxN
  tmp = 0;
  for k = 1:N
    conver(N-minN+1,k) = (nchoosek(N,k)*( (k/(N+k))^k )*(N/(k+N)))/(k^(k-1)/factorial(k-1));
  end
end
plot(conver(:,1))
plot(conver(:,2))
plot(conver(:,3))
plot(conver(:,4))

plot(conver(:,6))

plot(conver)


X = 1:5;
Y = (X.^X)./factorial(X);
plot(X,Y,'o')
