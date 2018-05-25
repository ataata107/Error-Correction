%-------------------------------------------------Enrollment---------------------------------------------------------------------------------------------------------------------------------

%Load the iris code
clc
clear all;
load('file1.mat')
iris = arrB;
iris_first = iris(1,:);                                %Manually Input the one you want to enroll
iris_first=reshape(iris_first,160,64,[]);       %Manually input the blocks that you want to divide the iris code

%Reed Solomon Encoding
m = 7;           % Number of bits per symbol
ns = 2^m-1;     % Codeword length 
ks = 3;           % Message length
K = randi([0 1], ks,m);  %Random Key
msg = gf(K',m);     %Message in galois field remember it is a transpose of the original key.

code = rsenc(msg,ns,ks); %encoding RS to make it ns blocks instead of k blocks of m bits each
code = code';  %Take the transpose for above comment to be true
%Hadamard
k = m-1;    %Number of bits per hadamard encoded code
H = hadamard(2^k);
H = mod(-H+2*ones(2^k),3);
minus_H= mod(hadamard(2^k)+2*ones(2^k),3);
HC_k = [H;minus_H];  %to know the codeword take the corresponding row
pseudo_iris = HC_k(1:size(code,1),:);  %Codeword
pseudo_iris_copy = pseudo_iris;
HC_k1 = [hadamard(2^k);-hadamard(2^k)]; %The main HC(k) as described in second paper


%shuffled iris code
shuff_iris=zeros(size(iris_first,1),size(iris_first,2));
shuff_key = randi([0 1], size(iris_first,1),1);
counter1=1;
counter2=0;
for i=1:size(iris_first,1)
          if (shuff_key(i)==1)
                shuff_iris(counter1,:)= iris_first(i,:);
                counter1 = counter1 + 1;
          else
              
               
              shuff_iris(size(iris_first,1)-counter2,:)=iris_first(i,:);
              counter2 = counter2 + 1;
               
          end
end
%zero insertion
shuff_iris_copy = shuff_iris;
counter = 0;
zero_inserted = 7;
for i=1:size(shuff_iris_copy,2)+zero_inserted
          if(rem(i,10) ==0)                      %Manually decide this to make this equal to the zero inserted
              shuff_iris(:,i)=zeros(size(iris_first,1),1);
              counter = counter+1;
          else
              shuff_iris(:,i) =shuff_iris_copy(:,i-counter);
          end
end
%theta lock and dis formation
shuff_iris=reshape(shuff_iris,1,size(shuff_iris,1)*size(shuff_iris,2),[]);
pseudo_iris=reshape(pseudo_iris,1,size(pseudo_iris,1)*size(pseudo_iris,2),[]); 
shuff_iris_main = shuff_iris(1,1:size(pseudo_iris,2));
theta_lock = xor(pseudo_iris,shuff_iris_main);
theta_dis = shuff_iris(1,size(pseudo_iris,2)+1:size(shuff_iris,2));
%Password Entry              
password = 'qwertyuiop';




%User
%Verification----------------------------------------------------------------------------------------------------------------------------------------------------------------------

% (Shuffling)
%randi([0 1], 160,64)       --Use this for an intruder
iris_second = iris(2,:);                             %Manually input the one you want to verify
iris_second=reshape(iris_second,160,64,[]);  %Manually input the size
shuff_iris_second=zeros(size(iris_second,1),size(iris_second,2));
counter1=1;
counter2=0;
for i=1:size(iris_second,1)
          if (shuff_key(i)==1)
                shuff_iris_second(counter1,:)=iris_second(i,:);
                counter1 = counter1 + 1;
          else
              
              shuff_iris_second(size(iris_second,1)-counter2,:)=iris_second(i,:);
              counter2 = counter2 + 1;  
          end
end
%zero insertion
shuff_iris_copy_second = shuff_iris_second;
counter = 0;
zero_inserted = 7;
for i=1:size(shuff_iris_copy_second,2)+zero_inserted
          if(rem(i,10) ==0)
              shuff_iris_second(:,i)=zeros(size(iris_second,1),1);
              counter = counter+1;
          else
              shuff_iris_second(:,i) =shuff_iris_copy_second(:,i-counter);
          end
end
shuff_iris_second=reshape(shuff_iris_second,1,size(shuff_iris_second,1)*size(shuff_iris_second,2),[]);
shuff_iris_main_second = shuff_iris_second(1,1:size(pseudo_iris,1)*size(pseudo_iris,2));
theta_test_main = shuff_iris_main_second;

theta_ps = xor(theta_lock,theta_test_main);
theta_ps=reshape(theta_ps,2^m-1,2^k,[]);
%Hadamard Decoding
w_prime = 2*theta_ps - ones(2^m-1,2^k);
decoded_Hadamard = w_prime*(HC_k1(1:size(code,1),:))';
[Y,I] = max(decoded_Hadamard,[],2);
counter = 1;
decoded = zeros(2^m-1,2^k);
decoded_code = zeros(127,m);
code = gf2dec(code,m,137);
code = reshape(code,2^m-1,m,[]);
for i=I'
    decoded(counter,:) = pseudo_iris_copy(i,:);
    
    decoded_code(counter,:) = code(i,:);
    counter = counter+1;
end
%reed solomon decoding
msg_new = gf(decoded_code,m);
decoded_rs = rsdec(msg_new',ns,ks);
decoded_rs = gf2dec(decoded_rs',m,137);
decoded_rs = reshape(decoded_rs,ks,m,[]);
if ( K == decoded_rs)
    disp('Its a match');
end