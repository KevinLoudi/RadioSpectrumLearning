% Min Kyeong Lee  
% PPM (Prediction with partial match)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inputs selected by a user:
%       maxContextOrder - The maximum context number
%       in - input string
%       alpha - what kind of symbols are considered to calculate
%               the initial probabilities and the cumulative probabilities
%               (Y): All alphabets including the space and the escape symbol
%               (N): Only alphabets used in the input string including the
%                    escape symbol
clear all
maxContextOrder = input('What is the maximum context order? ');
in = input('Please type your input: ', 's');
alpha = input('Do you want to use all of alphabets as symbols? (y/n): ', 's');

% Replace any space in the input with the symbol '{'
lengthIn = length(in);
for i = 1:lengthIn
    in = strrep(in, ' ', '{');
end

% Initialize the probability (iP)
% and the cumulative probability (iC) of each symbol
% depended on selected type of symbols by a user.
if alpha == 'y'
    % All symbols including space({) and esc(~)
    symbs = 'abcdefghijklmnopqrstuvwxyz{~';
else
    symbs = [unique(in),'~'];
end
lengthSymbs = length(symbs); 
iP(1:lengthSymbs) = 1/lengthSymbs;
iC(1) = 0;
for j = 2:lengthSymbs
    iC(j) = iC(j-1) + iP(j-1);
end

out=cell(1,4);
context = cell(1, maxContextOrder+1);
contextNoOut = [{}];
outMatchWhole = [{}];

% Read all elements from the input.
for ele = 1:lengthIn
    % Define the order N to start for each element.
    if ele <= maxContextOrder
        order = ele-1;
    else
        order = maxContextOrder;
    end

    % Initialize the escape symbol and the order-N contexts.
    symbOrder = lengthSymbs;
    tableNum = order + 1;
    ord = order;
    numMatchWhole = 1;
    
    while tableNum >= 1    
        matchWhole = strmatch({in((ele-ord):ele)}, context{tableNum}, 'exact');
        if isempty(matchWhole)
            if tableNum > 1
                matchPart = strncmp(context{tableNum}, in((ele-ord):(ele-1)), ord);
                if ~sum(matchPart)
                    % No Partial match
                    out{1} = [out{1}; iP(symbOrder)];
                    out{2} = [out{2}; iC(symbOrder)];
                    out{3} = [out{3}; {[in(ele-ord:ele-1) symbs(symbOrder)]}];
                    out{4} = [out{4}; ord];
                
                    contextNoOut = [contextNoOut; in((ele-ord):ele)];
                    context{tableNum} = [context{tableNum}; {in((ele-ord):ele)}];
                    context{tableNum} = [context{tableNum}; {[in(ele-ord:ele-1) symbs(symbOrder)]}];
                else
                    matchPart = strmatch(1, matchPart, 'exact');
                    matchPart = max(matchPart);      
                    % Partial match
                    out{1} = [out{1}; p{tableNum}(matchPart)];
                    out{2} = [out{2}; c{tableNum}(matchPart)];
                    out{3} = [out{3}; context{tableNum}(matchPart)];
                    out{4} = [out{4}; ord];  
                
                    contextNoOut = [contextNoOut; in((ele-ord):ele)];
                    context{tableNum} = [context{tableNum}; {in((ele-ord):ele)}];
                end
            else
                % iP, iC on Order-0 table
                if isempty(context{tableNum})
                    out{1} = [out{1}; iP(symbOrder)];
                    out{2} = [out{2}; iC(symbOrder)];
                    out{3} = [out{3}; {symbs(symbOrder)}];
                    out{4} = [out{4}; ord]; 
                else
                    num = size(context{tableNum},1);
                    out{1} = [out{1}; p{tableNum}(num)];
                    out{2} = [out{2}; c{tableNum}(num)];
                    out{3} = [out{3}; context{tableNum}(num)];
                    out{4} = [out{4}; ord]; 
                end
                eleOrder = findstr(in(ele), symbs);
                out{1} = [out{1}; iP(eleOrder)];
                out{2} = [out{2}; iC(eleOrder)];
                out{3} = [out{3}; symbs(eleOrder)];
                out{4} = [out{4}; ord-1];
            end
        else
            if numMatchWhole == 1
                % Whole matches
                out{1} = [out{1}; p{tableNum}(matchWhole)];
                out{2} = [out{2}; c{tableNum}(matchWhole)];
                out{3} = [out{3}; context{tableNum}(matchWhole)];
                out{4} = [out{4}; ord];
                numMatchWhole = numMatchWhole + 1;
            else
                outMatchWhole = [outMatchWhole; context{tableNum}(matchWhole)];
            end
        end
        tableNum = tableNum - 1;
        ord = ord - 1;
    end

    % Calculate the frequencies, the probabilities, 
    % and the cumulative probabilities
    outTemp = [out{3}; contextNoOut; outMatchWhole];

    context = cell(1, maxContextOrder+1);
    f=cell(1, maxContextOrder+1);
    p=cell(1, maxContextOrder+1);
    c=cell(1, maxContextOrder+1);

    for i = 1:size(outTemp,1)
        strLength = length(outTemp{i});
        context{strLength} = [context{strLength}; {outTemp{i}}];
    end

    tableNum = maxContextOrder + 1;
    for j = 1:tableNum
        contextU{j} = unique(sort(context{j}));
        for k = 1:size(contextU{j},1)
            f{j}(k) = length(strmatch(contextU{j}(k), context{j}));
        end

        if j == 1  
            p{j}(1:length(f{j})) = f{j}(1:length(f{j}))/sum(f{j});
            c{j}(1) = 0;
            for l = 2:length(f{j})
                c{j}(l) = c{j}(l-1) + p{j}(l-1);
            end
        else
            maxN = 0;
            for m = 1:size(contextU{j},1)
                match = findstr(contextU{j}{m}, '~');
                if ~isempty(match)
                    totalF = sum(f{j}(maxN+1:m));
                    p{j}(maxN+1:m) = f{j}(maxN+1:m)/totalF;
                    c{j}(maxN+1) = 0;
                    for n = maxN+2:m
                        c{j}(n) = c{j}(n-1) + p{j}(n-1);
                    end
                    maxN = m;
                end    
            end
        end
    end
    context = contextU;
end

% Display Output
outOut = [{'P(x)'} {'C(x)'} {'x'} {'Order'}];
for i = 1:size(out{1},1)
    len = length(out{3}{i});
    if len > 1
        out{3}{i} = out{3}{i}(len);
    end
    outOut = [outOut; out{1}(i) out{2}(i) out{3}(i) out{4}(i)];
end
disp(outOut);

% Display Tables
table = cell(1,maxContextOrder+1);
tableT = [{'context'} {'f'} {'p'}];
conT = []; ff = []; pp = [];
for j = 1:maxContextOrder + 1
    for k = 1:size(context{j},1)
        match = findstr(context{j}{k}, '~'); 
        if ~isempty(match)
            len = length(context{j}{k});
            context{j}{k} = context{j}{k}(len);
        end
            table{j} = [table{j}; {context{j}{k}}];
    end    
    tab{j} = [{table{j}}, {f{j}'}, {p{j}'}];
    tableT = [tableT; {tab{j}{:}}];
    conT = [conT; tableT{j+1,1}];
    ff = [ff; tableT{j+1,2}];
    pp = [pp; tableT{j+1,3}];   
end
disp(tableT);
disp(conT); disp(ff); disp(pp);