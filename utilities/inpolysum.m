function sum_in = inpolysum(p,bb)
  % check for NaNs
  rm = isnan(p(:,1)); 
  p(rm,:) = [];

  bb(end,:)=[];
  in = inpoly(p,bb);
  sum_in = sum(in); 

end