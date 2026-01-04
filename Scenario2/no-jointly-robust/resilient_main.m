%
clc;clear all;close all;
N = 10; %
T = 1; %
s=0;
choos=2;choosess=2;
tt=0:T:1000;
num_time = length(tt); %
m=zeros(N,N);
p=zeros(N,num_time);
networks;
%
p2 = zeros(N,num_time);
pn2 = zeros(N,num_time);
v = zeros(N,num_time);
vn = zeros(N,num_time);
u = zeros(N,num_time);
uv = zeros(N,num_time);
unv2 = zeros(N,num_time);
xx=randperm(5);k=1;timebuff=zeros(2,num_time);
for i=1:1:N
    p2(i,1)=10*rand;
end
p(:,1)=p2(:,1);
%% 
for t = 2:1:num_time
timebuff(1,t)=t;timebuff(2,t)=choos;%
if mod(t,10) == 2       
        if k==5+1
           xx=randperm(5);k=1;
           choos = floor(N*rand)+1;
        end
        s = xx(k);
        k=k+1;       
end
G(:,:,t-1)=G1(:,:,s);%
    mm=zeros(N,N);mn=zeros(N,N);
    for i=1:1:N
        for j=1:1:N
            mm(i,j)=p2(j,t-1);%-1-randi([0 3])-tau(i,j)
            mn(i,j)=p(j,t-1);
        end
         u(i,t)=wsrFuno(G(:,:,t-1),i,mm(i,:)',p2(i,t-1),v(i,t-1),N);
         uv(i,t)=wsrFun(G(:,:,t-1),i,mn(i,:)',p(i,t-1),vn(i,t-1),N);
    end
    %
    p(:,t) = p(:,t-1)+ T*vn(:,t-1);% +0.5*T^2*uv(:,t) 
    vn(:,t) = vn(:,t-1) +T*uv(:,t);% 1*u(:,t)
    p2(:,t) = p2(:,t-1)+ T*v(:,t-1);%
    v(:,t) = v(:,t-1) +T*u(:,t);% 1*u(:,t)
            p2(choos,t) = -5; 
            v(choos,t)=0;
            p(choos,t) = -5; 
            vn(choos,t)=0;
end

%% 
figure(1)
CList=slanCL(1,1:N);
v2 = [0 0; 200 0; 200 10; 0.01 10];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor','blue','FaceAlpha',.05,'EdgeColor','none');
hold on;
text(15,12,'$\mathcal{J}=[0,10]$','FontSize',12,'FontWeight','normal','Color','black','Interpreter','latex');


pp1=plot(1:num_time,p2(1,1:num_time),'--','LineWidth', 1.5,'Color',CList(1,:));
for i=2:N
plot(1:num_time,p2(i,1:num_time),'--','LineWidth', 1.5,'Color',CList(i,:));
hold on;
end

for t = 3:num_time
    pp6=plot(t-2:t, p2(timebuff(2,t),t-2:t),'-','LineWidth', 1.5,'Color','r','DisplayName','Malicious agent');
end
legend([pp1 pp6],{'Normal agent','Malicious agent'},'FontSize',10);
xlabel('Times(s)','interpreter','latex');
ylabel('Position $x_i$','interpreter','latex');
xlim([0 200]);ylim([-10 15]);
set(gca,'FontSize',12, 'box','on');   

%% MSR
function u = wsrFun(A, i, rec_x, self_x, self_v, n)
    alpha=2;beta=1;theta=1;f=1;md=2;
    selfA = A+eye(n);%
    selfA = selfA(i,:);
    rec_x(i)=self_x;
    offset=0;
    rec = (selfA.*rec_x')';%
    rev = rec;
    rev(rev==0) = [];
    sortRec = sort(rev,'descend');%
    %
    
    %
    for k = 1:md*f
        %
        %
        for j = 1:size(rec,1)
            if (selfA(j) ~= 0) & (rec(j) == max(sortRec))
                selfA(j) = 0;
                sortRec(1) = [];
                break;
            end
        end
        %
        for j = 1:size(rec,1)
            if selfA(j) ~= 0 & rec(j) == min(sortRec)
                selfA(j) = 0;
                sortRec(size(sortRec,1)) = [];
                break;
            end
        end
    end
    if selfA(i) == 0
        theta = 0;
    end

    modifiedA(1:n) = 0;
    for j = 1:n
        if j ~= i
            modifiedA(j) = selfA(j)/(theta+sum(selfA)-selfA(i));
        end
    end
    offset=beta*modifiedA*(rec_x(i)-rec_x);
    u = -alpha*self_v-offset;%
    %
    if isnan(u)
        u = 0;
    end
end

%% AVMSR
function u = wsrFuno(A, i, rec_x, self_x, self_v, n)
    alpha=2;beta=1;theta=1;f=1;md=2;
    selfA = A+eye(n);%
    selfA = selfA(i,:);
    rec_x(i)=self_x;
    offset=0;

    rec = (selfA.*rec_x')';%
    rev = rec;
    rev(rev==0) = [];
    sortRec = sort(rev,'descend');%

    for k = 1:md*f
        % 
        % 
        for j = 1:size(rec,1)
            if (selfA(j) ~= 0) & (rec(j) == max(sortRec))
                selfA(j) = 0;
                sortRec(1) = [];
                break;
            end
        end
        %
        for j = 1:size(rec,1)
            if selfA(j) ~= 0 & rec(j) == min(sortRec)
                selfA(j) = 0;
                sortRec(size(sortRec,1)) = [];
                break;
            end
        end
    end
    if selfA(i) == 0
        theta = 0;
    end

    %
    modifiedA(1:n) = 0;
    for j = 1:n
        if j ~= i
            modifiedA(j) = selfA(j)/(theta+sum(selfA)-selfA(i));
        end
    end
    % 
    %
    offset=beta*modifiedA*(rec_x(i)-rec_x);
    u = -alpha*self_v-offset;%
    %
    if isnan(u)
        u = -self_v;
    end
    if offset==0  %|| u == 0
        u = -self_v;
    end  

end
