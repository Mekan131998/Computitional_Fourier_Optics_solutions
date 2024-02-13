% Defining the triangle function
function [out] = triangle(x)
    mask=abs(x)<=1;
    t=1-abs(x);
    out=t.*mask;
end