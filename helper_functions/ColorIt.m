function [colr] = ColorIt(input)
%
% creates RGB values for more pleasing data plotting colrs
%
% input can be a single letter (first letter of colr name) or full colr
% name: reb, orange, yellow, green, blue, purple
% or an index between 1-8

switch input
    
    case {'r','red',1}
        
        colr = [ 228 30 38 ];
        
    case {'b','blue',2}
        
        colr = [ 51 127 186 ];
        
    case {'y','yellow',3}
        
        colr = [ 206 200 104 ];
        
    case {'g','green',4}
        
        colr = [ 76 175 74 ];
        
    case {'p','purple',5}
        
        colr = [ 154 80 159 ];
        
    case {'o','orange',6}
        
        colr = [ 245 129 32 ];
        
    case {'m','magenta',7}
        
        colr = [ 170 48 93 ];
        
    case {'c','cyan',8}
        
        colr = [ 0 174 239 ];
        
    case {'k','black',9}
        
        colr = [0 0 0];
        
    case {'br','brown',10}
        
        colr = [87 65 47];
        
end

colr = colr/255;

