function r = perm(sequence)
%PERM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if size(sequence,2) <=1
   r = sequence;
   return;
end
r = [];
for i=1:size(sequence,2)
    if i ~= 1 && sequence(i-1) == sequence(i)
       continue; 
    else
        s = [sequence(1:i-1) sequence(i+1:end)];
        p = perm(s);
        for j=1:size(p,1)
           r = [r; [sequence(i:i) p(j,:)]]; 
        end
    end
    
end

end

