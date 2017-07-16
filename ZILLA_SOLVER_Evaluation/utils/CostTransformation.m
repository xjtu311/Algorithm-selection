function cost=CostTransformation(realcost, type)
  if strcmp(type,'UNIFORM')
       cost = realcost*0+1;
       return;
  end
    if strcmp(type, 'RAW')
       cost = abs(realcost);
       return;
    end
    if strcmp(type, 'SQRT')
       cost = sqrt(abs(realcost));
       return;
    end
    fprintf('This code does not supoort cost type %s!', type);
end
