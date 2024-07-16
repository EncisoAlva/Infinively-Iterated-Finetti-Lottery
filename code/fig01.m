% Points to be graphed
%   p      alph     pi   beta
%  -------------   ------------------
%   0.4    0.9      2/3   0.9*(1+pi)
%   0.5    1.0       1    1.0*(1+pi)
%   0.6    1.1      3/2   1.1*(1+pi)
%
% N must be multiple of 60
%
%
%   p      alph     pi   beta
%  -------------   ------------------
%   0.3    0.8      3/7   0.8*(1+pi)
%   0.5    1.0       1    1.0*(1+pi)
%   0.7    1.2      7/3   1.2*(1+pi)
%
% N must be multiple of 210

maxN = 10;

% collector
counterPB = 1;
counterNN = 1;

%p_s = [0.4, 0.5, 0.6];
p_s = 0.3:0.2:0.7;
%a_s = [0.9, 1.0, 1.1];
a_s = 0.8:0.2:1.2;
N_s = 210*(1:maxN);

COLL = zeros(maxN, length(p_s)*length(a_s));
for N = N_s
  disp(N)
  max_pi = max(p_s./(1-p_s));
  piN = round( max_pi*N );
  %
  T  = zeros(piN+1);
  T(1,1) = 1;
  for i = 1:piN
    T(i+1,i+1) = N/(N+i);
    T(i  ,i+1) = i/(N+i);
  end
  %
  idx  = 0:piN;
  xvec = 1-N./(N+idx);
  for p = p_s
    pp = p/(1-p);
    pp_idx = round( pp*N );
    %
    curr_pow = 0;
    curr_vec = xvec;
    for a = a_s
      beta  = (1+pp)*a;
      betaN = round( beta*N );
      %
      while curr_pow < betaN
        curr_pow = curr_pow + 1;
        curr_vec = curr_vec*T;
      end
      %
      COLL(counterNN, counterPB) = curr_vec(pp_idx);
      %
      counterPB = counterPB+1;
    end
  end
  counterNN = counterNN+1;
  counterPB = 1;
end

ccc = cell(length(p_s)*length(a_s), 1);
counterPB = 1;
for p = p_s
  pp = p/(1-p);
  for a = a_s
    beta  = (1+pp)*a;
    %
    if a > 1 + (1+p)*log(p/(1-p))
      ccc{counterPB} = 'blue';
    elseif a < 1 + (1+p)*log(p/(1-p))
      ccc{counterPB} = 'red';
    else
      ccc{counterPB} = 'black';
    end
    %
    counterPB = counterPB + 1;
  end
end

%%
mmm  = {'square', 'o', '*'};
ccc2 = {'blue', 'black', 'red'};
lll  = {':','-','--'};

%%
figure()
tiledlayout(2,2,'TileSpacing','Compact', 'Padding','tight');

tt = 0:0.001:1;
ft = max( 1+(1-tt).*log(tt./(1-tt)),0 );
ft(end) = 1;

%figure()
nexttile
plot(tt,ft)

grid on
xlabel('$p$',Interpreter='latex')
ylabel('$\alpha$',Interpreter='latex')
xlim([0,1])
ylim([0, 1.4])
hold on

for pidx = 1:length(p_s)
  p = p_s(pidx);
  for aidx = 1:length(a_s)
    a = a_s(aidx);
    %
    plot( p, a , ...
        'Color', ccc2{aidx}, ...
        'Marker', mmm{pidx} )
  end
end
legend({'$\alpha = 1 + (1-p)\ln\left( \frac{p}{1-p} \right)$'},Interpreter="latex",...
  Location='southeast')
set( legend , 'location' , 'southeast', 'box' , 'on','edgecolor','w' )

%%

ERR = zeros(maxN, length(p_s)*length(a_s));

for pidx = 1:length(p_s)
  p = p_s(pidx);
  pp = p/(1-p);
  for aidx = 1:length(a_s)
    a = a_s(aidx);
    beta  = (1+pp)*a;
    true_val = lambertw(pp*exp(pp-beta))/(1+lambertw(pp*exp(pp-beta)));
    %
    counterPB = (pidx-1)*length(a_s) + aidx ;
    %
    ERR(:, counterPB) = abs( COLL(:, counterPB)-true_val);
  end
end

%%
%figure()
nexttile
firstfig = true;
for pidx = 2
  p = p_s(pidx);
  pp = p/(1-p);
  for aidx = 1:length(a_s)
    a = a_s(aidx);
    beta  = (1+pp)*a;
    true_val = lambertw(pp*exp(pp-beta))/(1+lambertw(pp*exp(pp-beta)));
    %
    counterPB = (pidx-1)*length(a_s) + aidx ;
    %
    semilogy( N_s, abs( COLL(:, counterPB)-true_val ), ...
      'Color', ccc2{aidx}, ...
      'Marker', mmm{pidx}, ...
      'LineStyle', lll{aidx} )
    %
    if firstfig
      hold on
      firstfig = false;
    end
  end
end
grid on
xlabel('$N$',Interpreter='latex')
ylabel('$\left| m_N^*(p, \alpha) - \mu^*(p, \alpha) \right|$',Interpreter='latex')
xlim([300,floor(max(N_s/100))*100])
ylim([min(ERR(:)), max(ERR(:))])
%subtitle('$p=0.3$',Interpreter='latex')
legend({'$p=0.5, \alpha=0.8$','$p=0.5, \alpha=1.0$','$p=0.5, \alpha=1.2$'},Interpreter='latex')
set( legend , 'edgecolor','w' )
yticks([0.0002, 0.0005, 0.001, 0.002])

%%
%figure()
nexttile
firstfig = true;
for pidx = 1
  p = p_s(pidx);
  pp = p/(1-p);
  for aidx = 1:length(a_s)
    a = a_s(aidx);
    beta  = (1+pp)*a;
    true_val = lambertw(pp*exp(pp-beta))/(1+lambertw(pp*exp(pp-beta)));
    %
    counterPB = (pidx-1)*length(a_s) + aidx ;
    %
    semilogy( N_s, abs( COLL(:, counterPB)-true_val ), ...
      'Color', ccc2{aidx}, ...
      'Marker', mmm{pidx}, ...
      'LineStyle', lll{aidx} )
    %
    if firstfig
      hold on
      firstfig = false;
    end
  end
end
grid on
xlabel('$N$',Interpreter='latex')
ylabel('$\left| m_N^*(p, \alpha) - \mu^*(p, \alpha) \right|$',Interpreter='latex')
xlim([300,floor(max(N_s/100))*100])
ylim([min(ERR(:)), max(ERR(:))])
%subtitle('$p=0.3$',Interpreter='latex')
legend({'$p=0.3, \alpha=0.8$','$p=0.3, \alpha=1.0$','$p=0.3, \alpha=1.2$'},Interpreter='latex')
set( legend , 'edgecolor','w' )
yticks([0.0002, 0.0005, 0.001, 0.002])

%%
%figure()
nexttile
firstfig = true;
for pidx = 3
  p = p_s(pidx);
  pp = p/(1-p);
  for aidx = 1:length(a_s)
    a = a_s(aidx);
    beta  = (1+pp)*a;
    true_val = lambertw(pp*exp(pp-beta))/(1+lambertw(pp*exp(pp-beta)));
    %
    counterPB = (pidx-1)*length(a_s) + aidx ;
    %
    semilogy( N_s, abs( COLL(:, counterPB)-true_val ), ...
      'Color', ccc2{aidx}, ...
      'Marker', mmm{pidx}, ...
      'LineStyle', lll{aidx} )
    %
    if firstfig
      hold on
      firstfig = false;
    end
  end
end
grid on
xlabel('$N$',Interpreter='latex')
ylabel('$\left| m_N^*(p, \alpha) - \mu^*(p, \alpha) \right|$',Interpreter='latex')
xlim([300,floor(max(N_s/100))*100])
ylim([min(ERR(:)), max(ERR(:))])
%subtitle('$p=0.3$',Interpreter='latex')
legend({'$p=0.7, \alpha=0.8$','$p=0.7, \alpha=1.0$','$p=0.7, \alpha=1.2$'},Interpreter='latex')
set( legend , 'edgecolor','w' )
yticks([0.0002, 0.0005, 0.001, 0.002])

set(gcf,'Color','w')