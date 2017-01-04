%*******************************************************
% Author: Kevin
% Propose: CTW-VMM prediction model Matlab version
% Date: 4th Dec, 2017
% Environment: Matlab 2015b
%*******************************************************
%example sequence
seq='01100110000101010110101001111100000000001000001000010001';
cs=alphabet('01');

%test all algorithms
ALGS = {'LZms', 'LZ78', 'PPMC', 'DCTW', 'BinaryCTW', 'PST'};
%param values
params.ab_size = size(cs);
params.d = 6;
params.m = 2;
params.s = 8;
params.pMin = 0.006;
params.alpha= 0;
params.gamma = 0.0006;
params.r = 1.05;
params.vmmOrder = params.d;

%use alphabet 0 1
disp('---------------------------------------------------');
disp('working with CS={0, 1 }');
disp('---------------------------------------------------');
disp(' ');

% 3. run each of the VMM algorithms
    i=5;
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    jVmm = vmm_create(map(cs, seq),  ALGS{i}, params);

    disp(sprintf('Pr(0 | 01011) = %f', vmm_getPr(jVmm, map(cs,'0'), map(cs,'01011'))));
    disp(sprintf('Pr(1 | 01011) = %f', vmm_getPr(jVmm, map(cs,'1'), map(cs,'01011'))));

    %averaged log-loss
    disp(sprintf('-lg(Pr(010111))=%f', vmm_logEval(jVmm,map(cs, '010111'))));
    disp('--------')
    disp(' ');
%end
