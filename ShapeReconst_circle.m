function [Rmesh] = ShapeReconst_circle(dn, nodeFBG, kin, THin, h)
% function [r] = ShapeReconst_circle(dn, nodeFBG, kin, THin, h)
% 원호 6개 이어붙이는 방식
% clear all
%%%%%%%%%%%%%%%%%     Test용 parameters      %%%%%%%%%%%%%%%%%%%%%
% nodeFBG = 4; %  노드 갯수
% dn = 50; % 구간 내 나누는 갯수
% h=[17 25 25 25 23];
% hs=[17 42 67 92 115];
% THin=[    1.6464;   -0.5481;   -0.8622;   -0.8083];
% kin=[    0.0008;    0.0014;    0.0034;    0.0034];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = linspace(0, h(1), dn);
for i=2:nodeFBG+1
    s = [s linspace(sum(h(1:i-1))+h(i)/dn, sum(h(1:i)), dn)];
end
dnh=dn/2; % dnh = half of dn
section = nodeFBG+1;
dF_TH = zeros(1,dnh);
T = zeros(dnh, 3);
N = zeros(dnh, 3);
B = zeros(dnh, 3);
Td =  zeros(dnh, 3);
Nd =  zeros(dnh, 3);
Bd =  zeros(dnh, 3);
l = zeros(dnh,3);
r = zeros(dn*section,3);
Nt = zeros(dn*section,3);
Bt = zeros(dn*section,3);
ds = h/dn;
ds(1) = s(dn)-s(dn-1); % s가 0부터니깐
dss = zeros(section,dn);
for i=1:section
    dss(i,:) = linspace(ds(i), ds(i), dn);
end


%%%%%%%%%%%%%%%%%%%%%% N0
F_K=zeros(1,dnh);
for j=1:dnh
    F_K(j)=kin(1)/dnh*j;
end
THi=THin(1);
tau=zeros(1,dnh);

% T, N, B 초기값
T(1,:)=[0 0 1];
N(1,:)=[cos(THi) sin(THi) 0];
i1=T(1,2)*N(1,3)-T(1,3)*N(1,2);
i2=T(1,3)*N(1,1)-T(1,1)*N(1,3);
i3=T(1,1)*N(1,2)-T(1,2)*N(1,1);
B(1,:)=[i1,i2,i3];
% B(1,:)=cross(T(1,:),N(1,:));

% Frenet-Serret Equation으로 위치벡터 r값 reconstruction
for i=1:dnh-1 
    Td(i+1,:)=F_K(i).*N(i,:)*dss(1,i+1);
    Nd(i+1,:)=(-F_K(i)*T(i,:)+tau(i)*B(i,:))*dss(1,i+1);
    i1=Td(i+1,2)*Nd(i+1,3)-Td(i+1,3)*Nd(i+1,2);
    i2=Td(i+1,3)*Nd(i+1,1)-Td(i+1,1)*Nd(i+1,3);
    i3=Td(i+1,1)*Nd(i+1,2)-Td(i+1,2)*Nd(i+1,1);
    Bd(i+1,:)=[i1,i2,i3];
%     Bd(i+1,:)=cross(Td(i+1,:),Nd(i+1,:));
    T(i+1,:)=(T(i,:)+Td(i+1,:))/norm(T(i,:)+Td(i+1,:));
    N(i+1,:)=(N(i,:)+Nd(i+1,:))/norm(N(i,:)+Nd(i+1,:));
    i1=T(i+1,2)*N(i+1,3)-T(i+1,3)*N(i+1,2);
    i2=T(i+1,3)*N(i+1,1)-T(i+1,1)*N(i+1,3);
    i3=T(i+1,1)*N(i+1,2)-T(i+1,2)*N(i+1,1);
    B(i+1,:)=[i1,i2,i3];
%     B(i+1,:)=cross(T(i+1,:),N(i+1,:));
end

for i=1:dnh-1
    l(i+1,:)=l(i,:)+T(i+1,:)*dss(1,i+1);
end
r(1:dnh,:)=l;
Nt(1:dnh,:)=N;
Bt(1:dnh,:)=B;

%%%%%%%%%%%%%%%%%%%%%%  N1-N4

tau=zeros(1,dn);

% T, N, B, l 초기값
T(1,:)=T(dnh,:);
N(1,:)=[cos(THin(1)) sin(THin(1)) 0];
B(1,:)=B(dnh,:);
l(1,:)=r(dnh,:);
for j=1:4
    F_K=ones(1,dn)*kin(j);
    % Frenet-Serret Equation으로 위치벡터 r값 reconstruction
    l(1,:)=l(1,:)+T(1,:)*dss(j,1);
    for i=1:dnh-1
        Td(i+1,:)=F_K(i).*N(i,:)*dss(j,i+1);
        Nd(i+1,:)=(-F_K(i)*T(i,:)+tau(i)*B(i,:))*dss(j,i+1);
        i1=Td(i+1,2)*Nd(i+1,3)-Td(i+1,3)*Nd(i+1,2);
        i2=Td(i+1,3)*Nd(i+1,1)-Td(i+1,1)*Nd(i+1,3);
        i3=Td(i+1,1)*Nd(i+1,2)-Td(i+1,2)*Nd(i+1,1);
        Bd(i+1,:)=[i1,i2,i3];
%         Bd(i+1,:)=cross(Td(i+1,:),Nd(i+1,:)); %%
        T(i+1,:)=(T(i,:)+Td(i+1,:))/norm(T(i,:)+Td(i+1,:));
        N(i+1,:)=(N(i,:)+Nd(i+1,:))/norm(N(i,:)+Nd(i+1,:));
        i1=T(i+1,2)*N(i+1,3)-T(i+1,3)*N(i+1,2);
        i2=T(i+1,3)*N(i+1,1)-T(i+1,1)*N(i+1,3);
        i3=T(i+1,1)*N(i+1,2)-T(i+1,2)*N(i+1,1);
        B(i+1,:)=[i1,i2,i3];
%         B(i+1,:)=cross(T(i+1,:),N(i+1,:)); %%
        l(i+1,:)=l(i,:)+T(i+1,:)*dss(j,i+1);
    end
    for i=dnh:dn-1
        Td(i+1,:)=F_K(i).*N(i,:)*dss(j+1,i+1);
        Nd(i+1,:)=(-F_K(i)*T(i,:)+tau(i)*B(i,:))*dss(j+1,i+1);
        i1=Td(i+1,2)*Nd(i+1,3)-Td(i+1,3)*Nd(i+1,2);
        i2=Td(i+1,3)*Nd(i+1,1)-Td(i+1,1)*Nd(i+1,3);
        i3=Td(i+1,1)*Nd(i+1,2)-Td(i+1,2)*Nd(i+1,1);
        Bd(i+1,:)=[i1,i2,i3];
%         Bd(i+1,:)=cross(Td(i+1,:),Nd(i+1,:)); %%
        T(i+1,:)=(T(i,:)+Td(i+1,:))/norm(T(i,:)+Td(i+1,:));
        N(i+1,:)=(N(i,:)+Nd(i+1,:))/norm(N(i,:)+Nd(i+1,:));
        i1=T(i+1,2)*N(i+1,3)-T(i+1,3)*N(i+1,2);
        i2=T(i+1,3)*N(i+1,1)-T(i+1,1)*N(i+1,3);
        i3=T(i+1,1)*N(i+1,2)-T(i+1,2)*N(i+1,1);
        B(i+1,:)=[i1,i2,i3];
%         B(i+1,:)=cross(T(i+1,:),N(i+1,:)); %%
        l(i+1,:)=l(i,:)+T(i+1,:)*dss(j+1,i+1);
    end
    r(dn*j-(dnh-1):dn*j+dnh,:)=l;
    Nt(dn*j-(dnh-1):dn*j+dnh,:)=N;
    Bt(dn*j-(dnh-1):dn*j+dnh,:)=B;
   
    if j~=4
        % T, N, B, l 초기화
        T(1,:)=T(dn,:);
        N(1,:)=[cos(THin(j+1)) sin(THin(j+1)) 0];
        B(1,:)=B(dn,:);
        l(1,:)=l(dn,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%% N5

F_K=zeros(1,dnh);
for j=1:dnh
    F_K(j)=kin(4)/dnh*(dnh-j);
end
tau=zeros(1,dnh);

% T, N, B 초기값
T(1,:)=T(dn,:);
N(1,:)=[cos(THin(4)) sin(THin(4)) 0];
B(1,:)=B(dn,:);
l(1,:)=l(dn,:)+T(1,:)*dss(5,1);

% Frenet-Serret Equation으로 위치벡터 r값 reconstruction
for i=1:dnh-1 
    Td(i+1,:)=F_K(i).*N(i,:)*dss(5,i+1);
    Nd(i+1,:)=(-F_K(i)*T(i,:)+tau(i)*B(i,:))*dss(5,i+1);
    i1=Td(i+1,2)*Nd(i+1,3)-Td(i+1,3)*Nd(i+1,2);
    i2=Td(i+1,3)*Nd(i+1,1)-Td(i+1,1)*Nd(i+1,3);
    i3=Td(i+1,1)*Nd(i+1,2)-Td(i+1,2)*Nd(i+1,1);
    Bd(i+1,:)=[i1,i2,i3];
%     Bd(i+1,:)=cross(Td(i+1,:),Nd(i+1,:));

    T(i+1,:)=(T(i,:)+Td(i+1,:))/norm(T(i,:)+Td(i+1,:));
    N(i+1,:)=(N(i,:)+Nd(i+1,:))/norm(N(i,:)+Nd(i+1,:));
    i1=T(i+1,2)*N(i+1,3)-T(i+1,3)*N(i+1,2);
    i2=T(i+1,3)*N(i+1,1)-T(i+1,1)*N(i+1,3);
    i3=T(i+1,1)*N(i+1,2)-T(i+1,2)*N(i+1,1);
    B(i+1,:)=[i1,i2,i3];
%     B(i+1,:)=cross(T(i+1,:),N(i+1,:));
end
for i=1:dnh-1
    l(i+1,:)=l(i,:)+T(i+1,:)*dss(5,i+1);
end
r((nodeFBG+1/2)*dn+1:(nodeFBG+1)*dn,:)=l(1:dnh,:);
Nt((nodeFBG+1/2)*dn+1:(nodeFBG+1)*dn,:)=N(1:dnh,:);
Bt((nodeFBG+1/2)*dn+1:(nodeFBG+1)*dn,:)=B(1:dnh,:);
% [u rmesh]=meshgrid(linspace(0,2.*pi,dn),1:dn*(nodeFBG+1));
u = meshgrid(linspace(0,2.*pi,20),r(:,1));
for i=1:3
    rmesh(:,:,i)=meshgrid(r(:,i),linspace(0,2.*pi,20))';
    Nmesh(:,:,i)=meshgrid(Nt(:,i),linspace(0,2.*pi,20))';
    Bmesh(:,:,i)=meshgrid(Bt(:,i),linspace(0,2.*pi,20))';
end

rho=0.5;
Rmesh(:,:,1) = rmesh(:,:,1)+rho*(Nmesh(:,:,1).*cos(u)+Bmesh(:,:,1).*sin(u));
Rmesh(:,:,2) = rmesh(:,:,2)+rho*(Nmesh(:,:,2).*cos(u)+Bmesh(:,:,2).*sin(u));
Rmesh(:,:,3) = rmesh(:,:,3)+rho*(Nmesh(:,:,3).*cos(u)+Bmesh(:,:,3).*sin(u));
% Rmesh(:,101,:) = r;
% for i=1:499
%     dr(i,:)=r(i+1,:)-r(i,:);
%     drM(i)=sqrt(dr(i,1)^2+dr(i,2)^2+dr(i,3)^2);
%     leng(i)=sum(drM(1:i));
% end
% plot(drM)
% leng(499)