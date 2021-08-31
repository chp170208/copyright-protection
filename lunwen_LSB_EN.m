function [stegoI,emD]=lunwen_LSB_EN( resimage,LS1,bitsnum)
stegoI = resimage; %构建一个与原载体图像相同大小的容器
data =LS1;
k = 0;
for j=1:bitsnum
    k = k + 1;
    stegoI(j)=bitset( resimage(j),1,data(k));%bitset函数将cover的第1位置换为data
end
emD=data(1:bitsnum);%emD为嵌入的数据
end