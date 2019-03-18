function nr = dims( field, ignore )
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

    s = size(field);
    nr = length(s) - length(find(s <= 1));    

%endfunction