% A simple hands-on tutorial 
% Browse the code to get the basics of how-to utilize this VMM tool

% 1. defining sequence and its alphabet
seq = 'abracadabra'; 
ab = alphabet('abracadabra'); %creates an alphabet object

% 2. testing all algs; param values ~match "best" values for text data (see
% tbl.8 in VMM paper)
ALGS = {'LZms', 'LZ78', 'PPMC', 'DCTW', 'BinaryCTW', 'PST'};
params.ab_size = size(ab);
params.d = 5; 
params.m = 2; 
params.s = 8;
params.pMin = 0.006;
params.alpha= 0;
params.gamma = 0.0006;
params.r = 1.05;
params.vmmOrder = params.d;

% use AB with size = 5
disp('---------------------------------------------------');
disp('working with AB={a, b, c, d, r }');
disp('---------------------------------------------------');
disp(' ');

% 3. run each of the VMM algorithms
for i=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    jVmm = vmm_create(map(ab, seq),  ALGS{i}, params);  %maps string s to a corresponding  alphanumeric-indices's alphabe
    disp(sprintf('Pr(a | br) = %f', vmm_getPr(jVmm, map(ab,'a'), map(ab,'br'))));
    disp(sprintf('Pr(b | br) = %f', vmm_getPr(jVmm, map(ab,'b'), map(ab,'br'))));
    disp(sprintf('Pr(c | br) = %f', vmm_getPr(jVmm, map(ab,'c'), map(ab,'br'))));
    disp(sprintf('Pr(d | br) = %f', vmm_getPr(jVmm, map(ab,'d'), map(ab,'br'))));
    disp(sprintf('Pr(r | br) = %f', vmm_getPr(jVmm, map(ab,'r'), map(ab,'br'))));

    disp(sprintf('-lg(Pr(abracadabra))=%f', vmm_logEval(jVmm,map(ab, 'abracadabra'))));
    disp('--------')
    disp(' ');
end

% 4. repeat the same scenario - this time with ascii AB (size(AB)=127)

disp(' ');
disp(' ');
disp('---------------------------------------------------');
disp('working with ascii AB (|AB|=127)');
disp('---------------------------------------------------');
disp(' ');
params.ab_size = 127;


for i=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    jVmm = vmm_create( seq,  ALGS{i}, params);

    disp(sprintf('Pr(a | br) = %f', vmm_getPr(jVmm, 'a', 'br')));
    disp(sprintf('Pr(b | br) = %f', vmm_getPr(jVmm, 'b', 'br')));
    disp(sprintf('Pr(c | br) = %f', vmm_getPr(jVmm, 'c', 'br')));
    disp(sprintf('Pr(d | br) = %f', vmm_getPr(jVmm, 'd', 'br')));
    disp(sprintf('Pr(r | br) = %f', vmm_getPr(jVmm, 'r', 'br')));

    disp(sprintf('-lg(Pr(abracadabra))=%f', vmm_logEval(jVmm, 'abracadabra')));
    disp('--------')
    disp(' ');
end
