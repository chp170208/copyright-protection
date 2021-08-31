tic
% 显示图像，并转为灰度图
%  img=imread('Peppers.bmp');
% img=imread('2.bmp');
 img=imread('Barbara.bmp');
%  img_rgb=rgb2gray(img);
 img_rgb=img;
 ima=img_rgb;
 swh=size(ima);
 
 %先按列取出元素
  N=input('请输入每隔多少个点取一个像素点');
  reimage1=[];
  reimage2=[];
  resimage1=[];
  resimage2=[];
 for i=1
     for j=1:swh(2)
 if mod(j,N)==0
     reima1=ima(i,j);
     reimage1=[reimage1 reima1];
 else
      resima1=ima(i,j);
      resimage1=[resimage1 resima1];
     
    
 end
     end
end
       o1=reimage1;
       o11= resimage1;
   for i=2:swh(1)
       for j=1:swh(2)
           if mod(j,N)==0
              reima2=ima(i,j);
              reimage2=[reimage2 reima2];
           else
              resima2=ima(i,j);
              resimage2=[resimage2 resima2];    
           
           end
       end
            reimage=cat(1,o1,reimage2);
            o1=reimage;
            reimage2=[];
            resimage=cat(1,o11,resimage2);
            o11=resimage;
            resimage2=[];
            
     
   end
      %Reimage为剩余像素点构成的矩阵，Resimage为缩略图像素点 
subplot(3,3,1);imshow(ima);title('a');
subplot(3,3,2);imshow(resimage);title('b');
subplot(3,3,3);imshow(reimage);title('c');



%DCT变换
B = blkproc(reimage,[8 8],'dct2');
% b1=max(max(B));
% b2=min(min(B));
% B=round(B*0.1);
% % % 画出DCT系数的直方图   
% B1=reshape(B,1,size(B,1)*size(B,2));%将矩阵转化成数组
% histogram(B1);
% ba=unique(B);
% hist(B1,ba);
% axis([0 200 0 20]);%axis函数用来设置画面横坐标及纵坐标的上下限


        %对灰度矩阵进行量化
NY=input('请输入量化阈值NY=');
B(abs(B)<NY)=0;   %改这里就可以改变量化
dctgrayImage1=B;

 %先将矩阵转换成数组
    cA2 =dctgrayImage1;
    b=size(cA2);
    k=1;
for i=1:b(1)
    for j=1:b(2)
        
        c(k)=cA2(i,j);
        k=k+1;
    end
end

         Q=c;
         c1=max(c);
         c2=min(c);
         T1=input('请输入缩放数据T1的值') ; 
         c=round(c*T1);
         c11=max(c);
         c22=min(c); 
         
          C=tabulate(c);
         L1=size(C,1);
   %  先分段  
        K=round(C(L1,1)-C(1,1))+1;
       for i=1:K
        a(i)=i+C(1,1)-1;
        d(i)=a(i);
       end
       index1=[];%用来记录系数落在哪些区间
        S=c;
        E_S1=a;
        E_D1=d;
    for i=1:length(S)
       for a=1:size(E_S1,2)
           if (S(i)>E_S1(a))&&(S(i)<E_D1(a))
               index1(i)=a;
           elseif (S(i)==E_S1(a))||(S(i)==E_D1(a))
               index1(i)=a;
               
           elseif (S(i)==E_S1(a))&&(S(i)==E_D1(a))
               index1(i)=a;
           end
       end
    end
         Y11=tabulate(index1);
         Y=Y11(1,1):Y11(length(Y11),1); 
         Y2=Y11(:,3);
         Y22=[];
    for i=1:size(Y2)
       Y22(i)=Y2(i)/100;
    end
            %用自带的哈弗曼函数
   dict = huffmandict(Y,Y22); %生成字典
   enco = huffmanenco(index1,dict); %编码   （估计误差转换成二进制数据流）
    %隐藏信息
  n1=T1*100;  %由于小数的二进制会变成0 ，所以选择把他放大5倍
  n2= dec2bin(n1); 		%将十进制数据转为字符,不知道为什么这里就是会从58变成57
  n3= boolean(n2-'0');
  n3= double(n3);
  n3=[zeros(1,8-length(n3)),n3];
  N1= dec2bin(N);%每隔几个点取像素
  N2= boolean(N1-'0');
  N2= double(N2);
  N2=[zeros(1,8-length(N2)),N2];
  DATAS1= dec2bin(-E_S1(1));
  DATAS= boolean(  DATAS1-'0');
  DATAS= double(  DATAS);
  DATAS =[zeros(1,8-length(DATAS)),DATAS];
  DATAE1= dec2bin(E_D1(length(E_D1)));
  DATAE= boolean(  DATAE1-'0');
  DATAE= double(  DATAE);
  DATAE =[zeros(1,8-length(DATAE)),DATAE];
  LS1=cat(2,n3, N2,DATAS,DATAE,enco);
  
  
  %嵌入信息
%   T=input('请输入位平面T的值');
  T=2;
  num = numel(resimage);%统计像素的个数
  bitsnum=numel(LS1);%统计待嵌入数据总长度
    if bitsnum<=num 
  [stegoI,emD]=lunwen_LSB_EN( resimage,LS1,bitsnum);%数据嵌入函数
  else
[stegoI,emD]=lunwen_LSB_EN2( resimage,LS1,T,num,bitsnum);
  end
   ans1=isequal(emD,LS1);
if ans1==1
    disp('嵌入数据与要隐藏数据一致');
else
    disp('嵌入数据与要隐藏数据不一致');
end
 %验证嵌入过程中是否出错
 subplot(3,3,5);
 imshow(stegoI,[]);
 title('嵌入信息的图像');
 
  %信息提取
  exD = LSB_de(stegoI,num,bitsnum);%(如果改变了T这里面也要去改变T)
%   exD=double(exD);
  ans2=isequal(emD,exD);%比较嵌入数据与提取数据是否一致
  if ans2==1
    disp('嵌入数据与提取数据一致');
else
    disp('嵌入数据与提取数据不一致');
  end
  
 n=(bin2dec(num2str(exD(1:8))))/100;%找出缩放倍数
 RN=bin2dec(num2str(exD(9:16)));%找出减小的数值
 RDATAS=bin2dec(num2str(exD(17:24)));%找出返正的数值
 RDATAE=bin2dec(num2str(exD(25:32)));%找出返正的数值
data_emd=exD(33:bitsnum);%
data_emd=double(data_emd);
deco = huffmandeco(data_emd,dict); %哈弗曼解码,

for i=1:length(deco)
    redata(i)=deco(i)- RDATAS-1;
  end
redata=round((redata))/n;
        k1=1;
   for i=1:b(1)
    for j=1:b(2)
        rcA(i,j)=redata(k1);
        k1=k1+1;
    end
   end
   y1 = blkproc(rcA,[8 8],'idct2');
    GD=zeros(swh(1),swh(2));
         kb=1;
        for j=1:swh(2)      
            
          if mod(j,N)==0
              GD(:,j)=y1(:,kb);
              kb=kb+1;
          end
        end
        ka=1;
      for j=1:swh(2)
          if mod(j,N)~=0
              GD(:,j)=stegoI(:,ka);
              ka=ka+1;
          end
      end
    
%  subplot(3,3,6);
 subplot;
 imshow(GD,[]);
 title('恢复图像,psnr=41dB');
%   计算峰值信噪比 
B=8;                %编码一个像素用多少二进制位
MAX=2^B-1;          %图像有多少灰度级
ima=double(ima);
MES=sum(sum((ima-GD).^2))/(swh(1)*swh(2));     %均方差
PSNR=20*log10(MAX/sqrt(MES));           %峰值信噪比


 %嵌入率
 em_ra=length(LS1)/(swh(1)*swh(2)*8);
 
 toc
 
 