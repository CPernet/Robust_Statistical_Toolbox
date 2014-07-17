function [M,MList]= OrthogContrasts(nTreats)

% Generates coefficients corresponding to one particular set of
% orthogonal contrasts for a set of variables
%
% FORMAT
% [contrasts, names]= limo_OrthogContrasts([2 2])
%
% INPUT
% nTreats is a vector with the levels of the various factors in an ANOVA
%
% OUTPUT
% M is a cell with the various contrast in it
% Mlist is a series of set names
%
% Taken from Matthew Nelson GenOrthogComps.m
% simplified by Cyril Pernet 16-09-2009

if length(nTreats) == 1
    error('there must be at least 2 factors')
end

alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';

dfs=nTreats-1;
nFacs=length(nTreats);
nCombs=prod(nTreats);

% Conpute main effects contrasts
M=cell(1,nFacs);
for iFac=1:nFacs
    M{iFac}=repmat(0,dfs(iFac),nCombs);
    
    if iFac==nFacs;
        nLowCombs=1;
    else
        nLowCombs=prod(nTreats(iFac+1:end));
    end
    
    if iFac==1
        nUpCombs=1;
    else
        nUpCombs=prod(nTreats(1:iFac-1));
    end
    
    for idf=1:dfs(iFac)
        tmpComps=repmat( 0,1,nLowCombs*nTreats(iFac) );
        tmpComps( 1:nLowCombs*idf )=1/idf;
        tmpComps( nLowCombs*idf+1:nLowCombs*(idf+1) )=-1;
        tmpComps=repmat(tmpComps,1,nUpCombs);
        M{iFac}(idf, : )=tmpComps;
        MList{iFac}=alphabet(iFac);
    end       
end

% Do interactions as products of the main effect contrasts
if nFacs>1        
    for curnFacs=2:nFacs     %start with 2 way interactions and go up to nFacs-way interactions
        for initFac=1:nFacs-(curnFacs-1)                       
            [M MList]=RLoop(curnFacs,initFac,repmat(1,1,nCombs), length(M)+1, '', M,MList,dfs,alphabet); 
        end
    end
end
end