tic
% ��ʾͼ�񣬲�תΪ�Ҷ�ͼ
%  img=imread('Peppers.bmp');
% img=imread('2.bmp');
 img=imread('Barbara.bmp');
%  img_rgb=rgb2gray(img);
 img_rgb=img;
 ima=img_rgb;
 swh=size(ima);
 
 %�Ȱ���ȡ��Ԫ��
  N=input('������ÿ�����ٸ���ȡһ�����ص�');
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
      %ReimageΪʣ�����ص㹹�ɵľ���ResimageΪ����ͼ���ص� 
subplot(3,3,1);imshow(ima);title('a');
subplot(3,3,2);imshow(resimage);title('b');
subplot(3,3,3);imshow(reimage);title('c');



%DCT�任
B = blkproc(reimage,[8 8],'dct2');
% b1=max(max(B));
% b2=min(min(B));
% B=round(B*0.1);
% % % ����DCTϵ����ֱ��ͼ   
% B1=reshape(B,1,size(B,1)*size(B,2));%������ת��������
% histogram(B1);
% ba=unique(B);
% hist(B1,ba);
% axis([0 200 0 20]);%axis�����������û�������꼰�������������


        %�ԻҶȾ����������
NY=input('������������ֵNY=');
B(abs(B)<NY)=0;   %������Ϳ��Ըı�����
dctgrayImage1=B;

 %�Ƚ�����ת��������
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
         T1=input('��������������T1��ֵ') ; 
         c=round(c*T1);
         c11=max(c);
         c22=min(c); 
         
          C=tabulate(c);
         L1=size(C,1);
   %  �ȷֶ�  
        K=round(C(L1,1)-C(1,1))+1;
       for i=1:K
        a(i)=i+C(1,1)-1;
        d(i)=a(i);
       end
       index1=[];%������¼ϵ��������Щ����
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
            %���Դ��Ĺ���������
   dict = huffmandict(Y,Y22); %�����ֵ�
   enco = huffmanenco(index1,dict); %����   ���������ת���ɶ�������������
    %������Ϣ
  n1=T1*100;  %����С���Ķ����ƻ���0 ������ѡ������Ŵ�5��
  n2= dec2bin(n1); 		%��ʮ��������תΪ�ַ�,��֪��Ϊʲô������ǻ��58���57
  n3= boolean(n2-'0');
  n3= double(n3);
  n3=[zeros(1,8-length(n3)),n3];
  N1= dec2bin(N);%ÿ��������ȡ����
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
  
  
  %Ƕ����Ϣ
%   T=input('������λƽ��T��ֵ');
  T=2;
  num = numel(resimage);%ͳ�����صĸ���
  bitsnum=numel(LS1);%ͳ�ƴ�Ƕ�������ܳ���
    if bitsnum<=num 
  [stegoI,emD]=lunwen_LSB_EN( resimage,LS1,bitsnum);%����Ƕ�뺯��
  else
[stegoI,emD]=lunwen_LSB_EN2( resimage,LS1,T,num,bitsnum);
  end
   ans1=isequal(emD,LS1);
if ans1==1
    disp('Ƕ��������Ҫ��������һ��');
else
    disp('Ƕ��������Ҫ�������ݲ�һ��');
end
 %��֤Ƕ��������Ƿ����
 subplot(3,3,5);
 imshow(stegoI,[]);
 title('Ƕ����Ϣ��ͼ��');
 
  %��Ϣ��ȡ
  exD = LSB_de(stegoI,num,bitsnum);%(����ı���T������ҲҪȥ�ı�T)
%   exD=double(exD);
  ans2=isequal(emD,exD);%�Ƚ�Ƕ����������ȡ�����Ƿ�һ��
  if ans2==1
    disp('Ƕ����������ȡ����һ��');
else
    disp('Ƕ����������ȡ���ݲ�һ��');
  end
  
 n=(bin2dec(num2str(exD(1:8))))/100;%�ҳ����ű���
 RN=bin2dec(num2str(exD(9:16)));%�ҳ���С����ֵ
 RDATAS=bin2dec(num2str(exD(17:24)));%�ҳ���������ֵ
 RDATAE=bin2dec(num2str(exD(25:32)));%�ҳ���������ֵ
data_emd=exD(33:bitsnum);%
data_emd=double(data_emd);
deco = huffmandeco(data_emd,dict); %����������,

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
 title('�ָ�ͼ��,psnr=41dB');
%   �����ֵ����� 
B=8;                %����һ�������ö��ٶ�����λ
MAX=2^B-1;          %ͼ���ж��ٻҶȼ�
ima=double(ima);
MES=sum(sum((ima-GD).^2))/(swh(1)*swh(2));     %������
PSNR=20*log10(MAX/sqrt(MES));           %��ֵ�����


 %Ƕ����
 em_ra=length(LS1)/(swh(1)*swh(2)*8);
 
 toc
 
 