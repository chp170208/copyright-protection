function [stegoI,emD]=lunwen_LSB_EN( resimage,LS1,bitsnum)
stegoI = resimage; %����һ����ԭ����ͼ����ͬ��С������
data =LS1;
k = 0;
for j=1:bitsnum
    k = k + 1;
    stegoI(j)=bitset( resimage(j),1,data(k));%bitset������cover�ĵ�1λ�û�Ϊdata
end
emD=data(1:bitsnum);%emDΪǶ�������
end