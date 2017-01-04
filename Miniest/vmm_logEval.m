function [logloss]=vmm_logEval(jVmm,seq,initialContext)
% calculates the length in bits of the  "compressed" representation of seq. 
% that is -log[ Pr ( seq | jVmm) ]
%
% [pr] = vmm_logEval(jVmm,seq, (*optional*) initialContext)
%
% param jVmm - a jVmm model (created by invoking vmm_create m-function)
% param seq - a sequence (string) over jVmm's alphabet
% param initialContext - (*optional*) initial context (string) of length >= VMM-order
% usage example:
%
% params.ab_size = 127
% params.d = 5
% jVmm = vmm_create('abracadabra', 'PPMC', params)
% vmm_getPr(jVmm,'a','br')
% vmm_logEval(jVmm,'abracadabra')
%
% --
%  Author: Ron Begleiter (Dotan DiCastro helped with a preliminary version). 
%  For more details visit http://www.cs.technion.ac.il/~ronbeg/vmm/
%--------------------------------------------------------------------------
jSeq = java.lang.String(seq);
if nargin==3,
    jContext = java.lang.String(initialContext);
    logloss = javaMethod('logEval',jVmm,jSeq,jContext);
else
    logloss = javaMethod('logEval',jVmm, jSeq);
end