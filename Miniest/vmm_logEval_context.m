function [logloss]=vmm_logEval_context(jVmm,seq,initialContext)

% calculates the length in bits of the  "compressed" representation of seq. 
% that is -log[ Pr ( seq | jVmm) ]
%
% [pr] = vmm_logEval(jVmm,seq)
%
% param jVmm - a jVmm model (created by invoking vmm_create m-function)
% param seq - a sequence (string) over jVmm's alphabet
% param initialContext - initial context (string) of length >= VMM-order
%
% usage example:
%
% params.ab_size = 127
% params.d = 5
% jVmm = vmm_create('abracadabra', 'PPMC', params)
% vmm_getPr(jVmm, 'a', 'br')
% vmm_logEval_context(jVmm,'abracadabra', 'aaaaa')
% vmm_logEval_context(jVmm,'abracadabra', 'bbbbb')
%
% --
% Author: Ron Begleiter (Dotan DiCastro helped with a preliminary version)
% For more details visit http://www.cs.technion.ac.il/~ronbeg/vmm/
%--------------------------------------------------------------------------

jSeq = java.lang.String(seq);
jContext = java.lang.String(initialContext);
logloss = javaMethod('logEval',jVmm,jSeq,jContext);