N = 10000;

probs    = 0:0.01:1;
expected = zeros(1,length(probs));

for j = 1:length(probs)
  p = 1-probs(j);

  M  = [1-p, 1-p, 0; p, 0, 0; 0, p, 1];
  Pi = [1,0,0]';

  collector  = zeros(3,N);
  collector3 = zeros(3,N);
  divi = 1;

  for i = 1:N
    Pi = M*Pi;
    divi = divi*3;
    collector(:,i)  = Pi;
    collector3(:,i) = divi*Pi;
  end

  %figure()
  %plot(collector(1,:))
  %hold on
  %plot(collector(2,:))
  %plot(collector(3,:))

  firstHH = zeros(1,N);
  for i = 2:N
    firstHH(i) = collector(3,i) - collector(3,i-1);
  end

  %figure()
  %plot(firstHH)

  idx = 1:N;

  expected(j) = idx*firstHH';
end

figure()
plot(probs,log10(expected))