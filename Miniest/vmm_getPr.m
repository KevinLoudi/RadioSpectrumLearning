function [pr] = vmm_getPr(jVmm,symbol,context)

% calculates the symbol probability: Pr(symbol | context, jVmm)
%
% [pr] = vmm_getPr(jVmm,symbol,context)
%
% param jVmm - a jVmm model (created by invoking vmm_create m-function)
% param symbol - a symbol (char) over the jVmm alphabet. 
% param context - a context sequence (string)
%
% usage example:
%
% params.ab_size = 127
% params.d = 5
% jVmm = vmm_create('abracadabra', 'PPMC', params)
% vmm_getPr(jVmm,'a','br')
% vmm_logEval(jVmm,'abracadabra')
%
% --
% This matlab wrapper for the java VMM package was co-writen by Dotan Di
% Castro and Ron Begleiter. For more details visit
% http://www.cs.technion.ac.il/~ronbeg/vmm/
%--------------------------------------------------------------------------
jContext = java.lang.String(context);
pr = javaMethod('predict',jVmm, double(symbol), jContext); %double(sym) = AB-index of sym