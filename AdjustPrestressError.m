%% pre-stress adjustment of cable-strut structures using noisy data
function AdjustPrestressError
%% read structural informations from excel file
clear all;
global NodeTable
global ElementNumber
global ElementTable
global NodeNumber
global ElementStiffness
global LoadTable
global Unfixed_NodeNumber
global Cable_number
FileName = 'Levy_s';
NodeTable = xlsread(FileName,'Node');
ElementTable = xlsread(FileName,'Element');
LoadTable = xlsread(FileName,'Load');
ElementNumber = size(ElementTable,1);
NodeNumber = size(NodeTable,1);
Cable_number = ElementNumber - sum(ElementTable(:,2));
Unfixed_NodeNumber = sum(NodeTable(:,5));
ElementStiffness = ElementTable(:,7);

%% input parameters
% Note:
% ActiveElementNumber:   number of active elements
% MeasuredElmeentNumber��number of measured elements
% SinAngleRange:         constraint of minimal angle of active elements
% MaxAdjustedNumber:     maximal number of adjusted elements
% MaxControlError:       if the maximal relative error is smaller than this value,quit the iteration step for selecting adjusted elements
% MaxDeltaError:         if the decrease of L-2 norm relative error is small than this value,quit the iteration step for selecting adjusted elements
% MinSingularValue:      constraint of the minimal singular value of the sensitivity matrix of selected adjusted elements
% MaxLengthError:        maximal relative length error imposed to the structure
% MaxMeasuredNoise:      maximal measuring noises

ActiveElementNumber = 20;
SinAngleRange = 0.00;
MeasuredElementNumber = 24;
MaxAdjustedNumber = 20;
MaxControlError = 0.1;
MaxDeltaError = 0.01;
MinSingularValue = 100;
MaxLengthError = 0.001;
MaxMeasuredNoise = 0.00;

%% get the relative force sensitivity matrix
[K_T,K_e,D,SRelative] = GetSensitivityMatrix;

%% select measured and active elements
ActiveElement = SelectActiveElement(SRelative,ActiveElementNumber,SinAngleRange);
MeasuredElement = SelectMeasuredElement(SRelative,MeasuredElementNumber);

%% display measured and active elements
DisplayActiveMeasuredElement(MeasuredElement,ActiveElement);

%% sensitivity matrix between active elements and measured elements
AdjustMatrix = SRelative(MeasuredElement,ActiveElement);

%% solve force with element length errors
InitialElementLength = ElementTable(:,3);
LengthError = diag(InitialElementLength) * (rand(ElementNumber,1)-0.5) * MaxLengthError * 2;

%LengthError = InitialElementLength * 30 * 1.2 * 10^-5;
% LengthError(1) = 0.01;
DesignedForce = ElementTable(:,6)';
% External_Force = zeros(1,54);
% External_Force([37:44 53]) = -0.5;
% ForceWithLengthError = DesignedForce' + K_e * D' * K_T^(-1) * External_Force';
LengthWithError = InitialElementLength + LengthError;
ForceWithLengthError = SolveByDRM(LengthWithError);
InitialRelativeError = (ForceWithLengthError - DesignedForce) * diag(DesignedForce)^(-1);

%% impose random length error to structure
% ForceNoise = ForceWithLengthError * diag((rand(ElementNumber,1)-0.5) * 2 * MaxMeasuredNoise);
load('E:\2020定\cable_force_adjust\Form-finding\0519\noise.mat');

ForceWithNoise = ForceWithLengthError + ForceNoise;
RelativeForceError = (ForceWithNoise - DesignedForce) * diag(DesignedForce)^(-1);
MeasuredRelativeError = RelativeForceError(MeasuredElement);

%% solve the adjust length of active element based on sparse regression
ActiveAdjustLength = ModifiedOMP(MeasuredRelativeError,AdjustMatrix,MaxAdjustedNumber,MaxControlError,MaxDeltaError,MinSingularValue);
AdjustLength = zeros(ElementNumber,1);
AdjustLength(ActiveElement) = ActiveAdjustLength';

%% solve the force after adjustment
LengthAfterAdjust = LengthWithError - AdjustLength;
ForceAfterAdjust = SolveByDRM(LengthAfterAdjust);

%% solve the relative force error after adjust
ErrorAfterAdjust = (ForceAfterAdjust - DesignedForce) * diag(DesignedForce)^(-1);
ActiveAdjustLength1 = regress(MeasuredRelativeError',AdjustMatrix);
AdjustLength1 = zeros(ElementNumber,1);
AdjustLength1(ActiveElement) = ActiveAdjustLength1';

%% solve the force after adjustment
LengthAfterAdjust1 = LengthWithError - AdjustLength1;
ForceAfterAdjust1 = SolveByDRM(LengthAfterAdjust1);

%% solve the relative force error after adjust
ErrorAfterAdjust1 = (ForceAfterAdjust1 - DesignedForce) * diag(DesignedForce)^(-1);

%% display adjust results
DisplayAdjust(LengthError,ForceWithLengthError,InitialRelativeError,ErrorAfterAdjust,ForceWithNoise,MeasuredElement,ActiveAdjustLength,ActiveElement,ActiveAdjustLength1,ErrorAfterAdjust1)
end

%% get relative force sensitivity matrix
function [K_T,K_e,D,result] = GetSensitivityMatrix
global NodeTable
global ElementNumber
global ElementTable
global NodeNumber
global ElementStiffness
%% get matrices
% Noting: 
% C :           connectivity matrix 
% D :           equlibrium matrix
% Q :           Force density matrix
% K_e :         element stiffness matrix
% K_E :         material sfiffness matrix
% K_G :         geometrical stiffness matrix
% K_T :         tangent stiffness matrix
% S :           force sensitivity matrix (Derived by FD-FM)
% S_D :         displacement sensitivity matrix
% UR :          self-stress states matrix
% S_F :         force sensitivity matrix (Derived by SVD-FM)
% S_relative :  relative force sensitivity matrix

C = zeros(ElementNumber,NodeNumber);
for i = 1 : 1 :ElementNumber
        C(i,ElementTable(i,4)) = 1;
        C(i,ElementTable(i,5)) = -1;
end
L = zeros(ElementNumber,ElementNumber);
dX = zeros(ElementNumber,ElementNumber);
dY = zeros(ElementNumber,ElementNumber);
dZ = zeros(ElementNumber,ElementNumber);
Q = zeros(ElementNumber,NodeNumber);
for i = 1 :1 : ElementNumber
    n = ElementTable(i,4);
    m = ElementTable(i,5);
    L(i,i) = LEN(NodeTable(n,2),NodeTable(n,3),NodeTable(n,4),NodeTable(m,2),NodeTable(m,3),NodeTable(m,4));
    dX(i,i) = NodeTable(n,2)-NodeTable(m,2);
    dY(i,i) = NodeTable(n,3)-NodeTable(m,3);
    dZ(i,i) = NodeTable(n,4)-NodeTable(m,4);
    Q(i,i) = ElementTable(i,6)/L(i,i);
end
D0 = [C'* dX/L;
    C'* dY /L;
    C'* dZ /L];
ZERO = diag([NodeTable(:,5);
    NodeTable(:,5);
    NodeTable(:,5)]);
ZERO (all(ZERO == 0, 2),:) = [];
D = ZERO * D0;
K_e = abs(diag(ElementTable(:,6))) * diag(ElementStiffness)* eye(ElementNumber)/L;
K_E = D * K_e * D';
G = C'* Q * C;
K_G0 = zeros(3 * NodeNumber,3 * NodeNumber);
K_G0(1:NodeNumber,1:NodeNumber) = G;
K_G0(NodeNumber + 1:2 * NodeNumber,NodeNumber + 1:2 * NodeNumber) = G;
K_G0(2 * NodeNumber + 1:3 * NodeNumber,2 * NodeNumber + 1:3 * NodeNumber) = G;
K_G = ZERO * K_G0* ZERO';
K_T = K_E + K_G;
S = K_e * D'* K_T^(-1) * D * K_e - K_e;
% [U1,~,~] = svd(D');
% self_stress_number = size(D,2) - rank(D);
% UR = U1(:,ElementNumber - self_stress_number + 1:ElementNumber);
% S_F1 = K_e * D'* (D * K_e * D')^(-1) * D * K_e - K_e;
% S_F = - UR*(UR'* (K_e)^(-1)*UR)^(-1)*UR';
% S_D = K_e * D'*inv(K_T);
SRelative = (diag(ElementTable(:,6)))^(-1) * S;
result = SRelative;
end

%% Structural non-linear analysis
function FinalForce = SolveByDRM(ElementLength)
global NodeTable
global Unfixed_NodeNumber
global ElementTable
global LoadTable
global ElementNumber
global ElementStiffness
%% get equilibrium equations
Real_Energy_function = @(X)0;
for i = 1:1:ElementNumber
    n = ElementTable(i,4);
    m = ElementTable(i,5);
    AxialStiffness = 0.5 * ElementStiffness(i)/ElementTable(i,3) * abs(ElementTable(i,6));
    if (n <= Unfixed_NodeNumber && m <= Unfixed_NodeNumber)
        fun = @(X)(AxialStiffness * (LEN(X(3*n-2),X(3*n-1),X(3*n),X(3*m-2),X(3*m-1),X(3*m))-ElementLength(i))^2);
        Real_Energy_function = @(X)Real_Energy_function(X) + fun(X);
    end
    if (n > Unfixed_NodeNumber && m <= Unfixed_NodeNumber)
        fun = @(X)(AxialStiffness * (LEN(NodeTable(n,2),NodeTable(n,3),NodeTable(n,4),X(3*m-2),X(3*m-1),X(3*m))-ElementLength(i))^2);
        Real_Energy_function = @(X)Real_Energy_function(X) + fun(X);
    end
    if (n <= Unfixed_NodeNumber && m > Unfixed_NodeNumber)
        fun = @(X)(AxialStiffness * (LEN(X(3*n-2),X(3*n-1),X(3*n),NodeTable(m,2),NodeTable(m,3),NodeTable(m,4))-ElementLength(i))^2);
        Real_Energy_function = @(X)Real_Energy_function(X) + fun(X);
    end
end
Load_number = size(LoadTable,1);
for i = 1:1:Load_number 
    index = LoadTable(i,2);
    DOF_number = 3 * (LoadTable(i,4)-1) + LoadTable(i,5);
    F_K = LoadTable(i,3);
    Ini_coor = LoadTable(i,6);
    if(index == 0)
        fun = @(X)(-X(DOF_number) * F_K);
        Real_Energy_function = @(X)Real_Energy_function(X) + fun(X);
    end
    if(index == 1)
        fun =  @(X)(0.5 * F_K * (X(DOF_number - Ini_coor))^2);
        Real_Energy_function = @(X)Real_Energy_function(X) + fun(X);
    end
end
x = sym('x',[3 * Unfixed_NodeNumber,1]);
Real_Energy_function_x = Real_Energy_function(x);
Unbalanced_force_function = jacobian(Real_Energy_function_x);
func = matlabFunction(Unbalanced_force_function);
X_result = zeros(3*Unfixed_NodeNumber,1);
FunExpression = 'funcc=@(x)func(x(';
for i = 1 : 1 : size(X_result,1)-1
    FunExpression = [FunExpression,num2str(i),'),x('];
end
FunExpression = [FunExpression,num2str(size(X_result,1)),'));'];
eval(FunExpression);
%% solve equations
for i = 1 : 1 :Unfixed_NodeNumber
    X_result(3*i-2) = NodeTable(i,2);
    X_result(3*i-1) = NodeTable(i,3);
    X_result(3*i) = NodeTable(i,4);
end
options=optimset('tolx',1e-100,'MaxFunEvals',1000000,'MaxIter',4000,'Display','off','Algorithm', 'levenberg-marquardt');
X_result_f = fsolve(funcc,X_result,options);
for i = 1:1:Unfixed_NodeNumber
    NodeTable(i,2) = X_result_f(3*i-2);
    NodeTable(i,3) = X_result_f(3*i-1);
    NodeTable(i,4) = X_result_f(3*i);
end
FinalElementLength = zeros(1,ElementNumber);
HessianMatrix = jacobian(Unbalanced_force_function);
HessianValue = double(vpa(subs(HessianMatrix,x,X_result),6));
FinalForce = zeros(1,ElementNumber);
for i = 1 :1 : ElementNumber
    n = ElementTable(i,4);
    m = ElementTable(i,5);
    AxialStiffness = 0.5 * ElementStiffness(i)/ElementTable(i,3) * abs(ElementTable(i,6));
    FinalElementLength(i) = LEN(NodeTable(n,2),NodeTable(n,3),NodeTable(n,4),NodeTable(m,2),NodeTable(m,3),NodeTable(m,4));
    FinalForce(i) = 2 * AxialStiffness * (FinalElementLength(i) - ElementLength(i));
end
end

%% display measured and active elements
function DisplayActiveMeasuredElement(MeasuredElement,ActiveElement)
global ElementTable
global NodeTable
global ElementNumber
figure(1);
clf;
for i = 1 :1 : ElementNumber
    if (~ismember(i,ActiveElement))
    index = ElementTable(i,2);
    n = ElementTable(i,4);
    m = ElementTable(i,5);
    if (index == 1)
        x = [NodeTable(n,2),NodeTable(m,2)];
        y = [NodeTable(n,3),NodeTable(m,3)];
        z = [NodeTable(n,4),NodeTable(m,4)];
        plot3 (x,y,z,'-black','LineWidth',4);
        hold on;
    end
    if (index == 0)
        x = [NodeTable(n,2),NodeTable(m,2)];
        y = [NodeTable(n,3),NodeTable(m,3)];
        z = [NodeTable(n,4),NodeTable(m,4)];
        plot3 (x,y,z,'-black');
        hold on;
    end
    end
end
if(~isempty(ActiveElement))
    N = size(ActiveElement,2);
    for i = 1 :1 : N
        n = ElementTable(ActiveElement(i),4);
        m = ElementTable(ActiveElement(i),5);
        x = [NodeTable(n,2),NodeTable(m,2)];
        y = [NodeTable(n,3),NodeTable(m,3)];
        z = [NodeTable(n,4),NodeTable(m,4)];
        plot3 (x,y,z,'-red','LineWidth',1);
        hold on;
    end
end
if(~isempty(MeasuredElement))
    N = size(MeasuredElement,2);
    for i = 1 :1 : N
        n = ElementTable(MeasuredElement(i),4);
        m = ElementTable(MeasuredElement(i),5);
        x = (NodeTable(n,2)+NodeTable(m,2))/2;
        y = (NodeTable(n,3)+NodeTable(m,3))/2;
        z = (NodeTable(n,4)+NodeTable(m,4))/2;
        plot3(x,y,z,'o','MarkerEdgeColor','black','MarkerFaceColor','g','MarkerSize',8)
        hold on;
    end
end
set(gca,'DataAspectRatio',[1 1 1]);
axis off
end

%% solve result of sparse regression using modified OMP algorithm
function result = ModifiedOMP(MeasuredRelativeError,AdjustMatrix,MaxAdjustedNumber,MaxControlError,MaxDeltaError,MinSingularValue)
result = [];
AdjustElement = [];
AdjustSubMatrix = [];
ResidualError(:,1) = MeasuredRelativeError';
AdjustLength(:,1) = zeros(size(AdjustMatrix,2),1);
for i = 1 : MaxAdjustedNumber
    if (max(abs((ResidualError(:,i)))) < MaxControlError)
        result = AdjustLength(:,i);
        return;
    end
    Objective = zeros(1,size(AdjustMatrix,2));
    for j = 1 : size(AdjustMatrix,2)
        if(~ismember(j,AdjustElement))
            Objective(j) = abs(ResidualError(:,i)' * AdjustMatrix(:,j));
        end
    end
    [~,AdjustElement(i)] = max(Objective);
    AdjustSubMatrix = AdjustMatrix(:,AdjustElement);
    if (min(svd(AdjustSubMatrix)) < sum(svd(AdjustSubMatrix))/MinSingularValue)
        result = AdjustLength(:,i);
        return;
    end
    SelectedAdjustLength = (AdjustSubMatrix'*AdjustSubMatrix)^(-1)*AdjustSubMatrix' * MeasuredRelativeError';
    AdjustLength(AdjustElement,i+1) = SelectedAdjustLength;
    ResidualError(:,i+1) = MeasuredRelativeError' - AdjustMatrix * AdjustLength(:,i+1);
    if (norm(ResidualError(:,i+1)) - norm(ResidualError(:,i)) > - MaxDeltaError)   
        result = AdjustLength(:,i);
        return;
    end
end
result = AdjustLength(:,i+1);
end

%% select adjust elements
function result = SelectActiveElement(SRelative,ActiveElementNumber,SinAngleRange)
result = [];
CosAngleRange = (1-SinAngleRange^2)^(0.5);
for i = 1 : 1 : ActiveElementNumber
    Objective = zeros(1,size(SRelative,2));
    MaxAngle = ones(1,size(SRelative,2));
    for j = 1 : size(SRelative,2)
        if(~ismember(j,result))
            if (i == 1)
                CandidateMatrix = SRelative(:,j);
                Objective(j) = sum(svd(CandidateMatrix));
                MaxAngle(j) = 0;
            end
            if (i > 1)
                CandidateMatrix = [SRelative(:,result) SRelative(:,j)];
                Objective(j) = sum(svd(CandidateMatrix));
                Angle = zeros(1,size(result,2));
                for k = 1 : size(result)
                    Angle(k) = abs(SRelative(:,result(k))' * SRelative(:,j))/(norm(SRelative(:,result(k))) * norm(SRelative(:,j)));
                end
                MaxAngle(j) = max(Angle);
            end
        end
    end
    CandidateActiveElementI = find(MaxAngle <CosAngleRange);
    [~,ActiveElementI] = max(Objective(CandidateActiveElementI));
    result(i) = CandidateActiveElementI(ActiveElementI);
end
result = sort(result);
end

%% select measured elements
function  result = SelectMeasuredElement(SRelative,ActiveElementNumber)
result = [];
for i = 1 : 1 : ActiveElementNumber
    Objective = zeros(1,size(SRelative,2));
    for j = 1 : size(SRelative,2)
        if(~ismember(j,result))
            if (i == 1)
                CandidateMatrix = SRelative(j,:);
                Objective(j) = sum(svd(CandidateMatrix));
            end
            if (i > 1)
                CandidateMatrix = [SRelative(result,:);SRelative(j,:)];
                Objective(j) = sum(svd(CandidateMatrix));
            end
        end
    end
    [~,MeasuredElementI] = max(Objective);
    result(i) = MeasuredElementI;
end
result = sort(result);
end

%% display adjust results
function DisplayAdjust(LengthError,ForceWithLengthError,InitialRelativeError,ErrorAfterAdjust,ForceWithNoise,MeasuredElement,ActiveAdjustLength,ActiveElement,ActiveAdjustLength1,ErrorAfterAdjust1)

global ElementNumber
figure(2);
set(gcf,'position',[597,524,800,300]);
clf;
bar(LengthError,'EdgeColor','flat','FaceColor','flat');
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Length error', 'Color','black','fontsize',18,'FontName','Times New Roman');
filename = 'E:\2020定\cable_force_adjust\Form-finding\0519\length error';
% saveas(gcf,filename,'emf');

figure(3);
set(gcf,'position',[597,524,800,300]);
clf;
plot(abs(InitialRelativeError),'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5);
hold on;
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Pre-stress error', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([0,ElementNumber+1]);
yticks(0:0.05:0.4);
filename = 'E:\2020定\cable_force_adjust\Form-finding\0519\initial error';
% saveas(gcf,filename,'emf')

err = ForceWithNoise - ForceWithLengthError;
MeasuredForceWithLenthError = ForceWithLengthError(MeasuredElement);
neg=err(MeasuredElement);
pos=err(MeasuredElement);
neg(find(neg>0))=NaN;
pos(find(pos<0))=NaN;
figure(4);
set(gcf,'position',[597,524,800,300])
clf;
errorbar(1:24,MeasuredForceWithLenthError,neg,pos,'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5, 'Capsize',12)
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Pre-stress with noise(kN)', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylim([50,350]);
yticks(50:50:350);
xticks(1:24);
xticklabels({'1','2','3','4','5','6','7','8','10','11','12','13','14','15','16','33','38','39','40','41','42','43','44','48'})
hold on;plot([16,17],[300,300],'-','LineWidth',1.5,'color',[0,0.45,0.74])
hold on;text(17.5,300,'Pre-stress with noise','FontSize',14,'FontName','Times New Roman')
hold on;plot(16.5,325,'s','MarkerSize',11,'MarkerEdgeColor','red','MarkerFaceColor','red')
hold on;text(17.5,325,'Actual pre-stress','FontSize',14,'FontName','Times New Roman')
hold on, rectangle('position',[15.5,285,9,55]);
filename = 'E:\2020定\cable_force_adjust\Form-finding\0519\prestress with noise';
% saveas(gcf,filename,'emf')


figure(5);
set(gcf,'position',[597,524,800,300]);
clf;
bar(ActiveAdjustLength,0.4,'EdgeColor','flat','FaceColor','flat');
hold on;
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Eement number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Adjust length(m)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xticks(1:20);
ylim([-0.015,0.015]);
yticks(-0.015:0.005:0.015);
xticklabels({'1','2','3','4','5','6','7','8','12','14','34','37','38','39','40','41','42','46','47','48'})
filename = 'E:\2020定\cable_force_adjust\Form-finding\0519\adjustlength';
% saveas(gcf,filename,'emf')

figure(6);
set(gcf,'position',[597,524,800,300]);
clf;
plot(abs(ErrorAfterAdjust),'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5);
hold on;
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([0,ElementNumber+1]);
ylabel('Pre-stress error', 'Color','black','fontsize',18,'FontName','Times New Roman');
filename = 'E:\2020定\cable_force_adjust\Form-finding\0519\afteradjust';
% saveas(gcf,filename,'emf')

figure(7);
set(gcf,'position',[597,524,800,300]);
clf;
bar(ActiveAdjustLength1,0.4,'EdgeColor','flat','FaceColor','flat');
hold on;
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Eement number', 'Color','black','fontsize',18,'FontName','Times New Roman');
ylabel('Adjust length(m)', 'Color','black','fontsize',18,'FontName','Times New Roman');
xticks(1:20);
ylim([-0.15,0.15]);
yticks(-0.15:0.05:0.15);
xticklabels({'1','2','3','4','5','6','7','8','12','14','34','37','38','39','40','41','42','46','47','48'})
filename = 'E:\2020定\cable_force_adjust\Form-finding\0519\adjustlengthbyLS';
% saveas(gcf,filename,'emf')

figure(8);
set(gcf,'position',[597,524,800,300]);
clf;
plot(abs(ErrorAfterAdjust1),'--s','MarkerSize',8,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1.5);
hold on;
set(gca,'FontSize',16,'FontName','Times New Roman');
xlabel('Element number', 'Color','black','fontsize',18,'FontName','Times New Roman');
xlim([0,ElementNumber+1]);
ylabel('Pre-stress error', 'Color','black','fontsize',18,'FontName','Times New Roman');
filename = 'E:\2020定\cable_force_adjust\Form-finding\0519\afteradjustbyLS';
% saveas(gcf,filename,'emf')end
end

%% return length of node(a,b,c) and node(d,e,f)
function result = LEN(a,b,c,d,e,f)
result = ((a-d)^2 + (b-e)^2 + (c-f)^2)^0.5;
end