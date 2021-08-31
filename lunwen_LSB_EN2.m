function [stegoI,emD]=lunwen_LSB_EN2( resimage,LS1,T,num,bitsnum)
stegoI =resimage; %构建一个与原载体图像相同大小的容器
data =LS1;

k=0;
for i=1:T
for j=1:num
    if k<num
    k = k + 1;
    stegoI(j)=bitset(resimage(j),i,data(k));
    elseif (k>=num)&&(k<bitsnum)
        k=k+1;
         stegoI(j)=bitset(stegoI(j),i,data(k));
    else
        break
    
end

end
emD=data(1:bitsnum);%emD为嵌入的数据
end