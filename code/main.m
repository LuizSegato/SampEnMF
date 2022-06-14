%This code calculates the SampEnMF for a given image

% [5] L. F. S. dos Santos, L. A. Neves, G. B. Rozendo, M. G. Ribeiro, M. Z. do Nascimento, T. A. A.
%Tosta, Multidimensional and fuzzy sample entropy (sampenmf) for quantifying h&e histological
%images of colorectal cancer, Computers in biology and medicine 103 (2018) 148–160.

%Test Benign image
% ------------------------------------------

%Input image
image = imread('../data/example_benign.png');

%attribute index
i=1;
    
%Sizes and amounts of sub-images
sub_im_size=64; %sub-images size (64 x 64 pixels)
N=size(image,1); %image height
M=size(image,2); %image width
sub_M=floor(M/sub_im_size); %amount of sub-images in the width of image
sub_N=floor(N/sub_im_size);  %amount of sub-images in the height of image
sub_amount=(sub_M*sub_N); %amounts of sub-images inside image
  
%Prepare the feature vector with 24 attributes (variations of parameters m and k)
SampEnMF=zeros(1,24);
m_limit=3; %limit for parameter m (window size) 
k_limit=0.5; %limit for parameter k (tolerance constant)
       
%SampEnMF with multiscale approach
%Sub-images preparation
%Centering the initial position of the first sub-image considering the
%amount of sub-images inside the image
for m=1:m_limit
    for k=0.06:0.06:k_limit
        %Necessary adjustment in height and width for scraps in a proportion of 25%
        if M-(sub_M*sub_im_size)>=floor(M/4) && N-(sub_N*sub_im_size)>=floor(N/4)
            im_x_from=floor((N-(sub_N*sub_im_size))/2);
            im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
            im_y_from=floor((M-(sub_M*sub_im_size))/2);
            im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
        %Necessary adjustment in width only for scraps in a proportion of 25%
        elseif M-(sub_M*sub_im_size)>=floor(M/4)
            im_x_from=1;
            im_x_to=sub_im_size;
            im_y_from=floor((M-(sub_M*sub_im_size))/2);
            im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
        %Necessary adjustment in height only for scraps in a proportion of 25%
        elseif N-(sub_N*sub_im_size)>=floor(N/4)
            im_x_from=floor((N-(sub_N*sub_im_size))/2);
            im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
            im_y_from=1;
            im_y_to=sub_im_size;
        %No necessary adjustment
        else
            im_x_from=1;
            im_x_to=sub_im_size;
            im_y_from=1;
            im_y_to=sub_im_size;
        end
        
        %Mean value of SampEnMF
        SampEnMF_mean=0;
        
        for im1=1:sub_M
            for im2=1:sub_N
                %Calculate SampEnMF for each sub-image with parameters m and k
                sub_image=image(im_x_from:im_x_to,im_y_from:im_y_to,:);
                entropy=calcSampEnMF(sub_image,k,m);
                SampEnMF_mean=SampEnMF_mean+entropy;
     
                %Sub-image sliding in height
                im_x_from=im_x_from+sub_im_size;
                im_x_to=im_x_to+sub_im_size;
            end
            %Sub-image sliding in width
            im_y_from=im_y_from+sub_im_size;
            im_y_to=im_y_to+sub_im_size;
            
            %Re-centering the height case necessary 
            if N-(sub_N*sub_im_size)>=floor(N/4)
                im_x_from=floor((N-(sub_N*sub_im_size))/2);
                im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
            else
                im_x_from=1;
                im_x_to=sub_im_size;
            end
        end
        SampEnMF(i)=SampEnMF_mean/double(sub_amount);
        fprintf('Value of atrribute %d (m = %d, e = %.2f) of SampEnMF = %.4f\n for benign image',i,m,k,SampEnMF(i));
        
        %Next attribute
        i=i+1;
    end
end

%Calculate metrics on the SampEnMF curve (second composition of vector)
%13 atrributes because the E (Maximum point scale) is the same for all m values
metrics=zeros(1,13);

metrics(1,1:5)=curve_metrics(SampEnMF(1,1:8));    %m=1
metrics(1,5:9)=curve_metrics(SampEnMF(1,9:16));   %m=2
metrics(1,9:13)=curve_metrics(SampEnMF(1,17:24)); %m=3

disp(metrics);

k=0.06:0.06:k_limit;

%Plot and save the SampEnMF curves for each scale m
%m=1
plot(k, SampEnMF(1, 1:8), 'b-s')
xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
title('m=1')

saveas(gcf, '../results/curve_m1_benign.png');

%m=2
plot(k, SampEnMF(1, 9:16), 'b-s')
xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
title('m=2')

saveas(gcf, '../results/curve_m2_benign.png');

%m=3
plot(k, SampEnMF(1, 17:24), 'b-s')
xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
title('m=3')

saveas(gcf, '../results/curve_m3_benign.png');

% -------------------------------------------


%Test Malignant image
% ------------------------------------------

%Input image
image = imread('../data/example_malignant.png');

%attribute index
i=1;
    
%Sizes and amounts of sub-images
sub_im_size=64; %sub-images size (64 x 64 pixels)
N=size(image,1); %image height
M=size(image,2); %image width
sub_M=floor(M/sub_im_size); %amount of sub-images in the width of image
sub_N=floor(N/sub_im_size);  %amount of sub-images in the height of image
sub_amount=(sub_M*sub_N); %amounts of sub-images inside image
  
%Prepare the feature vector with 24 attributes (variations of parameters m and k)
SampEnMF=zeros(1,24);
m_limit=3; %limit for parameter m (window size) 
k_limit=0.5; %limit for parameter k (tolerance constant)
       
%SampEnMF with multiscale approach
%Sub-images preparation
%Centering the initial position of the first sub-image considering the
%amount of sub-images inside the image
for m=1:m_limit
    for k=0.06:0.06:k_limit
        %Necessary adjustment in height and width for scraps in a proportion of 25%
        if M-(sub_M*sub_im_size)>=floor(M/4) && N-(sub_N*sub_im_size)>=floor(N/4)
            im_x_from=floor((N-(sub_N*sub_im_size))/2);
            im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
            im_y_from=floor((M-(sub_M*sub_im_size))/2);
            im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
        %Necessary adjustment in width only for scraps in a proportion of 25%
        elseif M-(sub_M*sub_im_size)>=floor(M/4)
            im_x_from=1;
            im_x_to=sub_im_size;
            im_y_from=floor((M-(sub_M*sub_im_size))/2);
            im_y_to=floor((M-(sub_M*sub_im_size))/2)+sub_im_size-1; 
        %Necessary adjustment in height only for scraps in a proportion of 25%
        elseif N-(sub_N*sub_im_size)>=floor(N/4)
            im_x_from=floor((N-(sub_N*sub_im_size))/2);
            im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
            im_y_from=1;
            im_y_to=sub_im_size;
        %No necessary adjustment
        else
            im_x_from=1;
            im_x_to=sub_im_size;
            im_y_from=1;
            im_y_to=sub_im_size;
        end
        
        %Mean value of SampEnMF
        SampEnMF_mean=0;
        
        for im1=1:sub_M
            for im2=1:sub_N
                %Calculate SampEnMF for each sub-image with parameters m and k
                sub_image=image(im_x_from:im_x_to,im_y_from:im_y_to,:);
                entropy=calcSampEnMF(sub_image,k,m);
                SampEnMF_mean=SampEnMF_mean+entropy;
     
                %Sub-image sliding in height
                im_x_from=im_x_from+sub_im_size;
                im_x_to=im_x_to+sub_im_size;
            end
            %Sub-image sliding in width
            im_y_from=im_y_from+sub_im_size;
            im_y_to=im_y_to+sub_im_size;
            
            %Re-centering the height case necessary 
            if N-(sub_N*sub_im_size)>=floor(N/4)
                im_x_from=floor((N-(sub_N*sub_im_size))/2);
                im_x_to=floor((N-(sub_N*sub_im_size))/2)+sub_im_size-1;
            else
                im_x_from=1;
                im_x_to=sub_im_size;
            end
        end
        SampEnMF(i)=SampEnMF_mean/double(sub_amount);
        fprintf('Value of atrribute %d (m = %d, e = %.2f) of SampEnMF = %.4f\n for malignant image',i,m,k,SampEnMF(i));
        
        %Next attribute
        i=i+1;
    end
end

%Calculate metrics on the SampEnMF curve (second composition of vector)
%13 atrributes because the E (Maximum point scale) is the same for all m values
metrics=zeros(1,13);

metrics(1,1:5)=curve_metrics(SampEnMF(1,1:8));    %m=1
metrics(1,5:9)=curve_metrics(SampEnMF(1,9:16));   %m=2
metrics(1,9:13)=curve_metrics(SampEnMF(1,17:24)); %m=3

disp(metrics);

k=0.06:0.06:k_limit;

%Plot and save the SampEnMF curves for each scale m
%m=1
plot(k, SampEnMF(1, 1:8), 'b-s')
xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
title('m=1')

saveas(gcf, '../results/curve_m1_malignant.png');

%m=2
plot(k, SampEnMF(1, 9:16), 'b-s')
xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
title('m=2')

saveas(gcf, '../results/curve_m2_malignant.png');

%m=3
plot(k, SampEnMF(1, 17:24), 'b-s')
xlabel('\epsilon', 'FontSize', 12, 'FontWeight','bold')
ylabel('SampEnMF', 'FontSize', 12, 'FontWeight','bold')
title('m=3')

saveas(gcf, '../results/curve_m3_malignant.png');
