%%

clear all;
close all;
clc;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ filename , pathname ]=uigetfile('Dataset\*.jpg','Select an Image');
input = imread([pathname,filename]);
input = imresize(input,[300,300]);
figure('Name','Input Image','NumberTitle','Off');
imshow(input);
axis off;
title('Input Image','fontsize',12,'fontname','Times New Roman','color','Black');

[row,col,cha] = size(input);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Box Segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[CA{1}, CH{1}, CV{1}, CD{1}] = dwt2(input,'db1');
[CA{2}, CH{2}, CV{2}, CD{2}] = dwt2(input,'db1');
[CA{3}, CH{3}, CV{3}, CD{3}] = dwt2(input,'db1');
[CA{4}, CH{4}, CV{4}, CD{4}] = dwt2(input,'db1');
CA4 = CA{4};
ELF = sum(sum(abs(CA4(:))));
EHF = [];
for i = 1:length(CA)
    CDi = CD{i};
    CHi = CH{i};
    CVi = CV{i};
    tEHF = sum(sum(CDi(:))) + sum(sum(CHi(:))) + sum(sum(CVi(:)));
end
EHF = sum(tEHF);
PLF = ( ELF/(ELF+EHF) ) * 100;
M = row;
N = col;
if PLF > 50
    S = (0.02*M*N)^(1/2);
elseif PLF <= 50
    S = (0.01*M*N)^(1/2);
end

[l, Am, Sp, d] = slic(input, round(S) , 10, 1, 'median');

SLICIMG  = drawregionboundaries(l, input, [255 0 0]);
figure('Name','SLIC Image','NumberTitle','Off');
imshow(SLICIMG);axis off;   % Showing SLIC Image
title('SLIC Image','fontname','Times New Roman','fontsize',12);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Box Feature Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIFT

HarrisThresh = 10; 
k = 0.055;NOfWindows = 4; 
BorderDistance = 4*NOfWindows;
NBHOOD = ones(3); 
sigma_nmbr = 9; 
[x y z]=size(SLICIMG);
if z==3
    img1gr = rgb2gray(SLICIMG);
else
    img1gr = SLICIMG;
end

SIFTpoints = SIFT( img1gr, NBHOOD, BorderDistance, HarrisThresh, k, sigma_nmbr );

SIFTxcoord = SIFTpoints(:,1);
SIFTycoord = SIFTpoints(:,2);
figure('Name','SIFT Features Extracted','NumberTitle','Off');
imshow(SLICIMG);axis off;hold on;   
plot(SIFTycoord,SIFTxcoord,'or')
title('SIFT Features Extracted','fontname','Times New Roman','fontsize',12);

SIFTimg = SLICIMG;
for i = 1:length(SIFTxcoord)
    SIFTimg(SIFTxcoord(i),SIFTycoord(i),1:3) = 255;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Box Feature Matching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;
for i = 1:30:row
    for j = 1:30:col
        SIFTimage{count} = SIFTimg(i:i+29,j:j+29,1:3);
        count = count + 1;
    end
end
[size1,size2,size3] = size(SIFTimage{1});
for i = 2:(count-1)
    
    tSIFTimage1 = SIFTimage{i};
    tSIFTimage2 = SIFTimage{i-1};
    cnt1 = 1;
    cnt2 = 1;
    for ii = 1:size1
        for jj = 1:size2
            if tSIFTimage1(ii,jj,1:3) == 255;
               xkeypoint1(cnt1) = ii;
               ykeypoint1(cnt1) = jj;
               cnt1 = cnt1+1;
            end 
            if tSIFTimage2(ii,jj,1:3) == 255;
               xkeypoint2(cnt2) = ii;
               ykeypoint2(cnt2) = jj;
               cnt2 = cnt2+1;
            end
        end
    end
    
end

for i = 1:length(xkeypoint1)
xa = xkeypoint1(i);
xb = xkeypoint2(i);
ya = ykeypoint1(i);
yb = ykeypoint2(i);
d_fa_fb(i) = ( (xa - xb) + (ya - yb) )^(1/2);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Labeled Feature Points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = dir('Dataset');
for itr = 3:length(dataset)
    filename  = dataset(itr).name;
    sampinput = imread(['Dataset\',filename]);
    sampinput = imresize(sampinput,[300,300]);
    Train(itr-2,:) = mean(sampinput(:));
end

Test = mean(input(:));

for i = 1:length(Train)
xa = mean(Train(i,:));
xb = mean(Test);
d_fafb(i) = ( abs(xa - xb) )^(1/2);
end

Label(1:11) = 4;
Label(12:22) = 5;
Label(23:33) = 6;
Label(34:44) = 1;
Label(45:55) = 2;
Label(56:66) = 3;
Label(67:77) = 7;

Distance = d_fafb;
[val1,loc1] = find(Distance == 0);
class = Label(loc1(1));

Original_Image = imread(['Samp\',num2str(class),'.jpg']);
Original_Image = imresize(Original_Image,[300,300]);

R1 = input(:,:,1);
G1 = input(:,:,2);
B1 = input(:,:,3);

R2 = Original_Image(:,:,1);
G2 = Original_Image(:,:,2);
B2 = Original_Image(:,:,3);

count = 1;
for i = 1:10:row
    for j = 1:10:col
        image1{count} = input(i:i+9,j:j+9,1:3);
        image2{count} = input(i:i+9,j:j+9,1:3);
        count = count + 1;
    end
end

Theta = [45,90,135,180,225,270,315,360];
TRsim = 15;
% for k = 1:length(Theta);
% for i = 1:10:row
%     for j = 1:10:col
%         R = R1(i:i+9,j:j+9);
%         G = G1(i:i+9,j:j+9);
%         B = B1(i:i+9,j:j+9);
%         Rb = R2(i:i+9,j:j+9);
%         Gb = G2(i:i+9,j:j+9);
%         Bb = B2(i:i+9,j:j+9);
%         Fc_LSib_theta = imrotate(R1,Theta(k));
%         Fc_LSi = ( R + G + B ) / 3; 
%         Fc_LSib = ( Rb + Gb + Bb ) / 3; 
%         Fc_LSi_R = R; 
%         Fc_LSi_G = G; 
%         Fc_LSi_B = B; 
%         Fc_LSib_R = R; 
%         Fc_LSib_G = G; 
%         Fc_LSib_B = B; 
%     end
% end
% end

for i = 1:(count-1)
    for j = 1:(count-1)
        diffval_1 = image1{i};
        diffval_2 = image2{j};
        sdiffval = diffval_1 - diffval_2;
        diffval = rgb2gray(sdiffval);
        differeceval = diffval<=TRsim;
        if unique(diffval(:)) == 0
           detectregion = ones(10,10);
        else
           detectregion = diffval;
        end
        finaldetectregion{i} = detectregion;
    end
end

findetectregion = zeros(row,col);
cnt = 1;
for i = 1:10:row
    for j = 1:10:col
        findetectregion(i:i+9,j:j+9) = finaldetectregion{cnt};
        cnt = cnt + 1;
    end
end

% for LSi = 1:(count-1)
%    
% end

% for LSi = 1:(count-1)
%     for LSj = 1:(count-1)
%        
%     end
% end
% 
% count = 1;
% for i = 1:10:row
%     for j = 1:10:col
%     
%         count = count + 1;
%     end
% end

tLabeled_Image = rgb2gray(input) - rgb2gray(Original_Image);
Labeled_Image = tLabeled_Image > 0;
figure('Name','Labeled Feature Points','NumberTitle','Off');
imshow(Labeled_Image);
axis off;
title('Labeled Feature Points','fontname','Times New Roman','fontsize',12);

if unique(Labeled_Image(:)) == 0
   
    msgbox('The Input Image Is Original');
    
else

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Merge Regions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

se = strel('disk',10);
Merged_Regions = imclose(Labeled_Image,se);
figure('Name','Merged Regions','NumberTitle','Off');
imshow(Merged_Regions);
axis off;
title('Merged Regions','fontname','Times New Roman','fontsize',12);

boundaries = bwboundaries(Merged_Regions);	
numberOfBoundaries = length(boundaries);
figure('Name','Detected Forged Region','NumberTitle','Off');
imshow(input);axis off;hold on;
for i1 = 1:numberOfBoundaries
    thisBoundary = boundaries{i1};
    plot(thisBoundary(:,2), thisBoundary(:,1), 'color' , 'r' , 'LineWidth', 2);
end
title('Detected Forged Region','fontname','Times New Roman','fontsize',12);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Performance Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Actual = Labeled_Image;
Predicted = Merged_Regions;

Performance = perf(Actual,Predicted);

Precision = Performance.precision;
Recall = Performance.recall;
F_Measure = Performance.Fmeasure;

rnames = {};
cnames = {'Precision','Recall','F-Measure'};
f=figure('Name','Performance Measures','NumberTitle','off');
t = uitable('Parent',f,'Data',[Precision,Recall,F_Measure],'ColumnName',cnames,'RowName',rnames);

%%
end