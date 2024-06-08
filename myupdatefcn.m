function text = myupdatefcn(~,event_obj,Y, AcNmb)
position = get(event_obj,'Position');
x = position(1);
y = position(2);
z = position(3); 
for a=1:length(AcNmb)
   if(Y(a,1)==x && Y(a,2)==y && Y(a,3)==z)
       idx = a;
       break;
   end
end
text = {['ANum:',AcNmb{idx}]};
end