clc
clear
s=1/2
l=1
bz=0:0.01:1
n_ml=(2*l+1)
n_ms=(2*s+1)
num_states=n_ml*n_ms//------------------possible states
ml=l:-1:-l
ms=s:-1:-s
state=zeros(1,2)
for i=1:(2*l+1)
    for j=1:(2*s+1)
        state=[state;[ml(i) ms(j)]]
    end
end
state=state(2:num_states+1,:)
lz_mat=diag(state(:,1))
sz_mat=diag(state(:,2))
jz_mat=lz_mat+sz_mat
for i=1:num_states
    sum_ls(i)=sum(state(i,:))
end
//--------------------construction of hamiltonian
k=1
while(k<=101)
    psi1=zeros(1,6)
    psi2=zeros(1,6)
    psi3=zeros(1,6)
    psi4=zeros(1,6)
    psi5=zeros(1,6)
    psi6=zeros(1,6)
    H=zeros(num_states,num_states)
    for i=1:num_states
        ml_i = state(i,1)
        ms_i = state(i,2)
        for j = 1:num_states
            if i ~= j
                ml_j = state(j,1)
                ms_j = state(j,2)
                
                // L+S- term: |ml_j, ms_j=+1/2> -> |ml_j+1, -1/2>
                if (ml_i == ml_j+1) & (ms_i == -1/2) & (ms_j == +1/2)
                    H(i,j) = (1/2)*sqrt(l*(l+1) - ml_j*(ml_j+1))
                end
                
                // L-S+ term: |ml_j, ms_j=-1/2> -> |ml_j-1, +1/2>
                if (ml_i == ml_j-1) & (ms_i == 1/2) & (ms_j == -1/2)
                    H(i,j) = (1/2)*sqrt(l*(l+1) - ml_j*(ml_j-1))
                end
            end
        end
        //Diagonal terms--------------
        H(i,i)=ml_i*ms_i+((ml_i+2*ms_i)*bz(k))
    end
//energy eigen value of block mj=+3/2
E(1,k)=H(1,1)

//expectation values value of block mj=+3/2
psi1=[1,0,0,0,0,0]
exp_lz(1,k)=psi1*lz_mat*psi1'
exp_sz(1,k)=psi1*sz_mat*psi1'
exp_jz(1,k)=psi1*jz_mat*psi1'

//energy eigen value of block mj=+1/2
H23=H(2:3,2:3)
[V23,E23]=spec(H23)
E23=diag(E23)
E(2,k)=E23(1)
E(3,k)=E23(2)
[E23_s,id23]=gsort(E23)
V23_s=V23(:,id23)
psi2(2:3)=V23_s(:,1)
psi3(2:3)=V23_s(:,2)
exp_lz(2,k)=psi2*lz_mat*psi2'
exp_lz(3,k)=psi3*lz_mat*psi3'
exp_sz(2,k)=psi2*sz_mat*psi2'
exp_sz(3,k)=psi3*sz_mat*psi3'
exp_jz(2,k)=psi2*jz_mat*psi2'
exp_jz(3,k)=psi3*jz_mat*psi3'

//energy eigen value of block mj=-1/2
H45=H(4:5,4:5)
[V45,E45]=spec(H45)
E45=diag(E45)
E(4,k)=E45(1)
E(5,k)=E45(2)

//expectation values-----------
[E45_s,id45]=gsort(E45)
V45_s=V45(:,id45)
psi4(4:5)=V45_s(:,1)
psi5(4:5)=V45_s(:,2)
exp_lz(4,k)=psi4*lz_mat*psi4'
exp_lz(5,k)=psi5*lz_mat*psi5'
exp_sz(4,k)=psi4*sz_mat*psi4'
exp_sz(5,k)=psi5*sz_mat*psi5'
exp_jz(4,k)=psi4*jz_mat*psi4'
exp_jz(5,k)=psi5*jz_mat*psi5'

//energy eigen value of block mj=-3/2
E(6,k)=H(6,6)
psi6=[0,0,0,0,0,1]
exp_lz(6,k)=psi6*lz_mat*psi6'
exp_sz(6,k)=psi6*sz_mat*psi6'
exp_jz(6,k)=psi6*jz_mat*psi6'
k=k+1
end
colors = ['r-','r-.','k-','b-','b-.','g','k-.','g-.','b--','r--']
clf()
xlabel("Magnetic field B")
ylabel("Energy eigenvalues mj")
title("Zeeman splitting of p orbital with L-S coupling")
for i=1:num_states
    plot(bz,E(i,:),colors(i))
end
legend(['²P3/2 mJ=+3/2';'²P1/2 mJ=+1/2';'²P3/2 mJ=+1/2';..
        '²P1/2 mJ=-1/2';'²P3/2 mJ=-1/2';'²P1/2 mJ=-3/2'],2)
figure
clf()
xlabel("Magnetic field B")
ylabel("expectation value")
title("Variation of expectaion values of Lz with Magnetic field")
for i = 1:num_states
    plot(bz, exp_lz(i,:),'b')
end
figure
clf()
xlabel("Magnetic field B")
ylabel("expectation value")
title("Variation of expectaion values of Sz with Magnetic field")
for i = 1:num_states
    plot(bz, exp_sz(i,:),'g')
end
figure
clf()
xlabel("Magnetic field B")
ylabel("expectation value")
title("Variation of expectaion values of Jz with Magnetic field")
for i = 1:num_states
    plot(bz, exp_jz(i,:),'r')
end
