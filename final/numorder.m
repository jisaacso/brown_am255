function numorder

spacings = [2/21,2/41,2/81,2/161,2/321];

for(i=1:5)
    h=spacings(i);
%     err(i)=euler(h,.4,i)
    err(i) = CN2(h,10)
end

for(j=1:4)
    order=log(err(j)/err(j+1))/log(spacings(j)/spacings(j+1))
end