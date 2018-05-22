%enrollment


clc
clear all;
load('file1.mat')
iris = arrB;
iris_first = iris(1,:);
iris_first=reshape(iris_first,160,64,[]);
%Reed Solomon
m = 7;           % Number of bits per symbol
ns = 2^m-1;     % Codeword length 
ks = 3;           % Message length
K = randi([0 1], ks,m);
msg = gf(K.',m)

code = rsenc(msg,ns,ks);
code = code';

%Hadamard
k = m-1;
H = hadamard(2^k);
H = mod(-H+2*ones(2^k),3);
minus_H= mod(hadamard(2^k)+2*ones(2^k),3);
HC_k = [H;minus_H];
pseudo_iris = HC_k(1:size(code,1),:);
pseudo_iris_copy = pseudo_iris;
HC_k1 = [hadamard(2^k);-hadamard(2^k)];
%shuffled iris code
shuff_iris=zeros(160,64);
shuff_key = randi([0 1], 160,1);
counter=0;
for i=1:160
          if (shuff_key(i)==1)
                shuff_iris(i-counter,:)=iris_first(i,:);
          else
              
              shuff_iris(161-i+counter,:)=iris_first(i,:);
              counter = counter + 1;  
          end
end
%zero insertion
shuff_iris_copy = shuff_iris;
counter = 0;
zero_inserted = 7;
for i=1:64+zero_inserted
          if(rem(i,10) ==0)
              shuff_iris(:,i)=zeros(160,1);
              counter = counter+1;
          else
              shuff_iris(:,i) =shuff_iris_copy(:,i-counter);
          end
end
%theta lock and dis formation
shuff_iris=reshape(shuff_iris,1,11360,[]);
pseudo_iris=reshape(pseudo_iris,1,8128,[]); 
shuff_iris_main = shuff_iris(1,1:8128);
theta_lock = xor(pseudo_iris,shuff_iris_main);
theta_dis = shuff_iris(1,8129:11360);
%Password Entry              
password = 'qwertyuiop';

%User verification
theta_test_main = shuff_iris_main;
theta_ps = xor(theta_lock,theta_test_main);
theta_ps=reshape(theta_ps,127,64,[]);
%Hadamard Decoding
w_prime = 2*theta_ps - ones(127,64);
decoded_Hadamard = w_prime*(HC_k1(1:size(code,1),:))';
[Y,I] = max(decoded_Hadamard,[],2);
counter = 1;
decoded = zeros(127,64);
decoded_code = zeros(127,m);
for i=I'
    decoded(counter,:) = pseudo_iris_copy(i,:);
    counter = counter+1;
    %decoded_code(counter,:) = code(i,:);
end
%reed solomon decoding
decoded_rs = rsdec(code',ns,ks)
