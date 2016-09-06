%% ROImaster
%
% quick and dirty ROI editor for neural imaging data
% distributed as is
%
% left mouse: start Xcorr seed
% right mouse: back to image

% L: label selected group
% space: add new ROI to current group and step to next
% b: add ROI to current group, dont step to next group

% F: display mean
% D:display std
% S: display PCA denoised average

% Q: quit
% X: delete allROIs from current group
% Z; delete ROI by clicking on itz
% U: undo/remove last ROIp
% N: step to next group
%
%
% 2013 jvoigts@mit.edu


%% read subset of stack

readInDirectory='/path_to_registered_image_stack/';

%expects pngs

files = dir([readInDirectory '*.png']);
numImages=numel(files)

stack= [];
c=0;
ff=fspecial('gaussian',11,0.5);

nstack=min(1500,numImages);

for i=1:nstack;
    c=c+1;
    
    if (rem(i,100)==0)
        fprintf('%d/%d (%d%%)\n',i,numImages,round(100*(i./nstack)));
    end;
    
    
    fnum=ceil(rand*(numImages-1));
    %fnum=i;
    
    %fnum=i;
    I=imread([readInDirectory 'registered_' num2str(fnum),'.png']);
    %  I=I(:,1:100);
    I=conv2(double(I),ff,'same');
    %imagesc(I);
    %drawnow;
    
    
    stack(:,:,i)=I;
    
    % ts(i)=mean(mean(I(370:376,272:287)));
end
%% load

load([readInDirectory(1:end-11),'ROIs.mat'])


selected_group=1;
%% reset ROIS

Rois=[];
Rois.date_started=date;
Rois.data_dir=readInDirectory;
Rois.N=0;
Rois.masks=[];  % binary masks
Rois.groups=[]; % int id per mask, assigning to groups
Rois.grouplabels=[]; %laber per roi, not per group
Rois.outlines=[] % mostly just for plotting

selected_group=1;

%% make local Xcorr and PCA
% code adapted from http://labrigger.com/blog/2013/06/13/local-cross-corr-images/


disp('computing local xcorr');
w=1; % window size

% Initialize and set up parameters
ymax=size(stack,1);
xmax=size(stack,2);
numFrames=size(stack,3);
ccimage=zeros(ymax,xmax);

for y=1+w:ymax-w
    
    if (rem(y,10)==0)
        fprintf('%d/%d (%d%%)\n',y,ymax,round(100*(y./ymax)));
    end;
    
    for x=1+w:xmax-w
        % Center pixel
        thing1 = reshape(stack(y,x,:)-mean(stack(y,x,:),3),[1 1 numFrames]); % Extract center pixel's time course and subtract its mean
        ad_a   = sum(thing1.*thing1,3);    % Auto corr, for normalization laterdf
        
        % Neighborhood
        a = stack(y-w:y+w,x-w:x+w,:);         % Extract the neighborhood
        b = mean(stack(y-w:y+w,x-w:x+w,:),3); % Get its mean
        thing2 = bsxfun(@minus,a,b);       % Subtract its mean
        ad_b = sum(thing2.*thing2,3);      % Auto corr, for normalization later
        
        % Cross corr
        ccs = sum(bsxfun(@times,thing1,thing2),3)./sqrt(bsxfun(@times,ad_a,ad_b)); % Cross corr with normalization
        ccs((numel(ccs)+1)/2) = [];        % Delete the middle point
        ccimage(y,x) = mean(ccs(:));       % Get the mean cross corr of the local neighborhood
    end
end


m=mean(ccimage(:));
ccimage(1,:)=m;
ccimage(end,:)=m;
ccimage(:,1)=m;
ccimage(:,end)=m;

pcaimage=ccimage;
if 0
    disp('computing PCA ROI prediction');
    % make PCA composite, this seems to display good roi candidates
    stack_v=zeros(nstack,size(stack,1)*size(stack,2));
    for i=1:nstack;
        x=stack(:,:,i);
        stack_v(i,:)=x(:);
    end
    stack_v=stack_v-mean(stack_v(:));
    [coeff, score] = pca(stack_v,'Economy','on','NumComponents',100);
    imcomponents=reshape(coeff',100,size(stack,1),size(stack,2));
    pcaimage=(squeeze(mean(abs(imcomponents(1:100,:,:)))));
end;
disp('done');


%% run ROImaster
figure(1);
run=1;
displayxc=0;
plotstd=1;
set(gca, 'position', [0.02 0.04 1 .94]);


updatexc=0;

Ngroups=60*3;
UIheight=max(size(stack,1),200);

stdim=(std(single(stack),[],3));
meanim=mean(stack,3);

f=normpdf([-10:10],0,1);b

roiUIpos=zeros(Ngroups,2);
for i=1:Ngroups
    roiUIpos(end-i+1,1)=-30*ceil(i/60);
    
    roiUIpos(end-i+1,2)=2+(rem(i-1,60)/60)*UIheight;
end;

while run
    clf; hold on;
    % plot UI
    for i=1:Ngroups
        ningroup=sum(Rois.groups==i);
        if (i==selected_group)
            text(roiUIpos(i,1),roiUIpos(i,2),[num2str(i) ,'-',num2str(ningroup)],'BackgroundColor',[.3 .8 1]);
            try
                text(-60,505,Rois.grouplabels{i},'color',[1 0 0]);
            end;
        else
            text(roiUIpos(i,1),roiUIpos(i,2),[num2str(i) ,'-',num2str(ningroup)]);
            try
                if numel(Rois.grouplabels{i})>0
                    text(roiUIpos(i,1),roiUIpos(i,2),[num2str(i) ,'-',num2str(ningroup)],'BackgroundColor',[.7 .7 .7]);
                    
                end;
            end;
            
        end;
    end;
    
    if displayxc
        
        if 0 %<--- switch interactive crosscorr display method here
            if  updatexc
                % 2 step refinement
                ref= (squeeze( stack(ceil(y),ceil(x),:) ));
                %ref=(mean(squeeze(mean(stack(ceil(y+[-1:1]),ceil(x+[-1:1]),:),1) ),1))';
                xc=I.*0;
                for i=1:2:size(stack,1) % first subsample every 2nd
                    for j=1:2:size(stack,2)
                        c=corrcoef(squeeze(stack(i,j,:)),ref);
                        xc(i,j)=c(1,2);
                    end
                end;
              %  imagesc(xc); daspect([1 1 1]); drawnow;
                for fillin=find(xc(:)>.1)' % fill in where we detected any corr >.1
                    [a,b]=ind2sub(size(I),fillin);
                    for i=[-1 0 1]+a
                        for j=[-1 0 1]+b
                            if (i>0) && (j>0) && (i<size(I,1))&& (j<size(I,2)) && ~(i==a && j==b)
                                c=corrcoef(squeeze(stack(i,j,:)),ref);
                                xc(i,j)=c(2,1);
                            end;
                        end
                    end;
                    
                end;
                xc(ceil(y),ceil(x))=0;
            end;
            imagesc(xc); daspect([1 1 1]);
            
        else
            %iterative region growing
            if  updatexc
                ref= (squeeze( stack(ceil(y),ceil(x),:) ));
                xc=I.*0;
                
                xc(ceil(y),ceil(x))=0.11; % seed
                
                it=1;
                while it<50
                    sig=find(xc>0.04); % <--- configure threshold here
                    mask=I.*0; mask(sig)=1;
                    mask=conv2(mask,ones(5),'same')>0;
                    update=find((xc==0).*(mask==1));
                    if numel(update)<1
                        it=500;
                    end;
                    for fillin=update' % fill in where we detected any corr >.1s
                        [a,b]=ind2sub(size(I),fillin);
                        c=corrcoef(squeeze(stack(a,b,:)),ref);
                        xc(a,b)=c(2,1);
                    end;
                    %imagesc(xc);
                    %  image((1-mask).*I./10+ ((xc*200)))
                    % daspect([1 1 1]);drawnow;
                    it=it+1;
                    
                end;
            end;
            
        end
        if updatexc
            xc(ceil(y+[-1:1]),ceil(x+[-1:1]))=xc(ceil(y+[-1:1]),ceil(x+[-1:1])).*.5;
        end;
        image((1-mask).*I./10+ ((xc*200)))
        daspect([1 1 1]);
    else
        if plotstd==1
            I=stdim;
            I=ccimage;
        elseif plotstd==0
            I=meanim;
        else
            I=pcaimage;
        end;
        imagesc(I); daspect([1 1 1]);
        if plotstd
            xlabel('plotting loca xcorr');
        else
            xlabel('plotting mean');
        end;
    end;
    
    %plot outlines
    for i=1:Rois.N
        if (selected_group==Rois.groups(i))
            plot(Rois.outlines{i}([1:end,1],1),Rois.outlines{i}([1:end,1],2),'color',[1 1 .2].*1);
        else
            plot(Rois.outlines{i}([1:end,1],1),Rois.outlines{i}([1:end,1],2),'color',[1 0 0].*1);
        end;
        text(mean(Rois.outlines{i}(:,1)),mean(Rois.outlines{i}(:,2)),num2str(Rois.groups(i)),'color',[1 1 1].*.6);
    end;
    
    %plot example data trace
    sel=find(selected_group==Rois.groups);
    tr=[];
    for i=1:numel(sel);
        xa=round(min(Rois.outlines{sel(i)}(:,1)));
        xb=round(max(Rois.outlines{sel(i)}(:,1)));
        ya=round(min(Rois.outlines{sel(i)}(:,2)));
        yb=round(max(Rois.outlines{sel(i)}(:,2)));
        
        for t=1:size(stack,3)
            tr(i,t)=mean(mean(stack(ya:yb,xa:xb,t).*Rois.masks{sel(i)}(ya:yb,xa:xb)));
        end;
    end;
    tr=conv2(tr,f,'same');
    if numel(sel)>0
        plot(((tr./max(tr(:)))'.*100)-200, repmat( [1:size(stack,3)],numel(sel),1)' );
    end;
    
    % plot last clicked point just as super quick indicator
    if updatexc
        quickprev=squeeze(mean(mean(stack([-1:1]+round(y),[-1:1]+round(x),:))));
        quickprev=conv2(quickprev',f,'same');
    end;
    try
        plot(((quickprev./max(quickprev(:)))'.*100)-200,  [1:size(stack,3)]' ,'k--');
    end;
    
    xlim([-200 size(stack,2)]);
    ylim([0 UIheight]);
    colormap gray;
    set(gca, 'position', [0 0 1 1]);
    
    
    [x,y,b]=ginput(1);
    
    if b==108 %l
        
        prompt = {'Enter new label:'};
        dlg_title = 'new group label';
        num_lines = 1;
        def = {''};
        newlabel = inputdlg(prompt,dlg_title,num_lines,def);
        
        Rois.grouplabels{selected_group}=newlabel;
        
    end;
    
    
    if b==117 % u -undo last addition/ delete last ROI
        %delete current group
        sel= Rois.N;
        [~,ii]=sort(-sel);
        sel=sel(ii); % gotta delete from the end to the beginnning for indexing reasons
        if numel(sel)>0
            for s=sel
                Rois.masks{s}=[];
                try
                    Rois.masks(cellfun(@(d) isempty(d),Rois.masks))=[];
                end;
                Rois.groups(s)=[];
                Rois.outlines{s}=[];
                try
                    Rois.outlines(cellfun(@(d) isempty(d),Rois.outlines))=[];
                end;
            end;
        end;
        Rois.N=numel(Rois.groups);
        
    end;
    
    if b==120 % x delet group
        %delete current group
        sel=find(selected_group==Rois.groups);
        [~,ii]=sort(-sel);
        sel=sel(ii); % gotta delete from the end to the beginnning for indexing reasons
        if numel(sel)>0
            for s=sel
                Rois.masks{s}=[];
                try
                    Rois.masks(cellfun(@(d) isempty(d),Rois.masks))=[];
                end;
                Rois.groups(s)=[];
                Rois.outlines{s}=[];
                try
                    Rois.outlines(cellfun(@(d) isempty(d),Rois.outlines))=[];
                end;
            end;
        end;
        Rois.N=numel(Rois.groups);
    end;
    
    if b==122 % z - delete nearest group
        sel=find(selected_group==Rois.groups);
        sel=[1:Rois.N];
        %find nearest
        dmin=20000;
        for i=1:numel(sel);
            j=sel(i);
            %select roi
            d= sqrt((x-mean(Rois.outlines{j}([1:end,1],1))).^2+mean(y-Rois.outlines{j}([1:end,1],2)).^2);
            if d<dmin
                dmin=d;
                seldelete=i;
            end;
        end;
        plot(Rois.outlines{seldelete}([1:end,1],1),Rois.outlines{seldelete}([1:end,1],2),'color',[1 0 0].*1,'LineWidth',2);
        text(0,0,'[z] again to delete','BackgroundColor',[1 1 1]);
        drawnow;
        
        
        [x,y,b]=ginput(1);
        if b==122 % x - check and delete
            Rois.N=Rois.N-1;
            Rois.masks(seldelete)=[];
            Rois.groups(seldelete)=[];
            Rois.outlines(seldelete)=[];
            Rois.labels(seldelete)=[];
        else
            
        end;
        
    end;
    
    
    if b==112 %p, play movie
        ii=I; %start with whats on screen
        mfactor=.3;
        for i=1:100;%size(stack,3);
            ii=(ii.*(1-mfactor))+stack(:,:,i).*mfactor;
            %clf;
            h=imagesc(ii);
            
            drawnow;
            delete(h);
        end;
    end;
    
    if b==110 %n, next group
        updatexc=0;
        selected_group=selected_group+1;
    end;
    
    if b==115 %s
          plotstd=2;
    end;
    if b==100 %d
        plotstd=1;
    end;
    if b==102  %f
        plotstd=0;
    end;
    if b==1; %left mouse
        if x>0
            displayxc=1;
            updatexc=1;
        else
            [~,selected_group]=min(((x-10)- roiUIpos(:,1)).^2 + (y- roiUIpos(:,2)).^2);
            updatexc=0;
        end;
    end;
    if b==3; %right mouseq
        displayxc=0;
    end;
    if b==113 % q
        
        disp('exited')
        run=0;
    end;
    
    if b==98 % advance and start drawing
        updatexc=0;
        selected_group=selected_group+1;
        t= imfreehand(gca,'Closed' ,1);
        t.setClosed(1);
        m=t.createMask;
        
        updatexc=0;
        
        Rois.N=Rois.N+1;
        
        Rois.masks{Rois.N}=logical(m);
        Rois.groups(Rois.N)=selected_group;
        r=t.getPosition;
        r=max(r,1); r(:,1)=min(r(:,1),size(I,2));r(:,2)=min(r(:,2),size(I,1));
        Rois.outlines{Rois.N}=r;
        
        Rois.labels{Rois.N}='';
    end;
    
    if b==32 %space
        t= imfreehand(gca,'Closed' ,1);
        t.setClosed(1);
        m=t.createMask;
        %  imagesc(mask);
        %xlim([-100 size(stack,2)]);
        % ylim([0 size(stack,1)]);
        % set(gca, 'position', [0 0 1 1]);
        % drawnow;
        %displayxc=0;
        
        updatexc=0;
        
        Rois.N=Rois.N+1;
        
        Rois.masks{Rois.N}=logical(m);
        Rois.groups(Rois.N)=selected_group;
        r=t.getPosition;
        r=max(r,1); r(:,1)=min(r(:,1),size(I,2));r(:,2)=min(r(:,2),size(I,1));
        Rois.outlines{Rois.N}=r;
        
        Rois.labels{Rois.N}='';
        selected_group=selected_group+1;
    end;
    
    if b==98 % 'b' add to current group
        t= imfreehand(gca,'Closed' ,1);
        t.setClosed(1);
        m=t.createMask;
        updatexc=0;
        Rois.N=Rois.N+1;
        Rois.masks{Rois.N}=logical(m);
        Rois.groups(Rois.N)=selected_group;
        r=t.getPosition;
        r=max(r,1); r(:,1)=min(r(:,1),size(I,2));r(:,2)=min(r(:,2),size(I,1));
        Rois.outlines{Rois.N}=r;
        Rois.labels{Rois.N}='';
    end;q
    title(readInDirectory);
end;



%% make neuropil ROIs for each ROI
allmasks=  Rois.masks{1};
for j=2:Rois.N
    allmasks = allmasks + Rois.masks{j};
end;
allmasks=allmasks>0;

vismasks=Rois.masks{1}.*0;
se=strel('disk',7); % <-- configure size here
for j=1:Rois.N
    Rois.np_masks{j} = (imdilate(Rois.masks{j},se)-allmasks)>0;
    
    hold off;
    vismasks=vismasks+Rois.np_masks{j};
    imagesc( vismasks);
    daspect([1 1 1]);
    drawnow;
    
    p=regionprops(Rois.np_masks{j});
    Rois.np_outlines{j}=p.BoundingBox;
    
end;

%% save Rois
save([Rois.data_dir(1:end-11),'ROIs.mat'],'Rois')
disp(['saved to ',[Rois.data_dir(1:end-11),'ROIs.mat']]);



%% Get F(roi) and F(neuropil) from an imagestack on disk.

dirlist=[];
%{
dirlist{1}='/media/data_2p/2p/27aug2013/TSeries-08272013-1319-001/registered/';
dirlist{2}='/media/data_2p/2p/27aug2013/TSeries-08272013-1319-003/registered/';
dirlist{3}='/media/data_2p/2p/27aug2013/TSeries-08272013-1319-005/registered/';
dirlist{4}='/media/data_2p/2p/27aug2013/TSeries-08272013-1319-007/registered/';
%}
%{
dirlist{1}='/media/data_2p/2p/27aug2013/TSeries-08272013-1319-002/registered/';
dirlist{2}='/media/data_2p/2p/27aug2013/TSeries-08272013-1319-004/registered/';
dirlist{3}='/media/data_2p/2p/27aug2013/TSeries-08272013-1319-006/registered/';
%}
dirlist{1}=readInDirectory;

roiValues_norm=[];
roiValues=[];

Ndirs=numel(dirlist);
disp([' selected ',num2str(Ndirs),' matching directories']);

c=0; % count trhoug images, in case its broken up into multiple directories

for dd=1:Ndirs
    disp('%%%%%%%%%%');
    disp(['dir ',num2str(dd),'/',num2str(Ndirs)]);
    data_dir=dirlist{dd}
    disp('%%%%%%%%%%');
    
    
    files = dir([data_dir '*.png']);
    numImages=numel(files)
    
    
    tic
    n=1;
    aa=[];
    %neuropilScale=0.7;
    for i=1:numImages
        c=c+1;
        if (rem(i,100)==0)
            fprintf('extracting %d/%d (%d%%)\n',i,numImages,round(100*(i./numImages)));
        end;
        
        imageToMeasure=uint16(imread([readInDirectory 'registered_' int2str(i)],'png'));
        
        
        for j=1:Rois.N
            
            
            xa=ceil(min(Rois.outlines{j}(:,1)));
            xb=floor(max(Rois.outlines{j}(:,1)));
            ya=ceil(min(Rois.outlines{j}(:,2)));
            yb=floor(max(Rois.outlines{j}(:,2)));
            roiValues(c,j)=mean(mean(imageToMeasure(ya:yb,xa:xb).*uint16(Rois.masks{j}(ya:yb,xa:xb))));
            
            if 1 % normalize?   % <---- neuropil correction
                xa=ceil(Rois.np_outlines{j}(1));
                xb=floor(Rois.np_outlines{j}(1)+Rois.np_outlines{j}(3));
                ya=ceil(Rois.np_outlines{j}(2));
                yb=floor(Rois.np_outlines{j}(2)+Rois.np_outlines{j}(4));
                m=mean(mean(imageToMeasure(ya:yb,xa:xb).*uint16(Rois.np_masks{j}(ya:yb,xa:xb))));
                roiValues_normby(c,j)=m;
                
                % -----------------
                roiValues_norm(c,j)=roiValues(c,j) - m.*.75; % <---- configure neuropil correction method here
                % -----------------
                
                
            end;
        end;
        
    end
    
end;


disp('ROI luminance extracted, not saved to disk yet');

%% plot extracted values grouped and mean subtracted
figure(1); clf; hold on;
f=normpdf([-10:10],0,1);
c=0;

[uu]=unique(Rois.groups);
c=0;
for i=1%:numel(uu)
    c=c+1;
    sel=find(uu(i)==Rois.groups);
    plot(conv2(roiValues(1:2:end,sel)- repmat(mean(roiValues(1:2:end,sel),1)',1,size(roiValues,1)/2)' ,f','same')+c*100);
    drawnow;
end;

%% look at xcors
xc=(corrcoef(roiValues));
xc=xc.*(1-eye(size(xc)));
clf;
imagesc(xc);

%% run spike extrction
    roiValues_norm(isnan(roiValues_norm))=0;
    roiValues_norm(isinf(roiValues_norm))=0;
    
    
    % set simulation metadata
    T       =  size(roiValues,1); % # of time steps% set simulation metadata
    
    fprintf('%d frames ',T);
    if T< 2000*50
        fprintf(' assume galvo ');
        V.dt    =  0.1;  % time step size
    else
        fprintf(' assume resonant ');
        V.dt    =  0.02;  % time step size
    end;
    fprintf('- dt=%f \n', V.dt);
    
    
    
    V.fast_poiss =0;
    V.est_lam=1;
    V.est_gam=1;
    V.est_b=0;
    
    % initialize params
    P.a     = 1;    % observation scale
    P.b     = 0.0;    % observation bias
    tau     = 2.6;    % decay time constant
    P.gam   = 1-V.dt/tau; % C(t) = gam*C(t-1)
    P.lam   = 0.5;  % firing rate = lam*dt
    P.sig   = 0.1;  % standard deviation of observation noise
    
    % simulate data
    N = poissrnd(P.lam*V.dt*ones(T,1)); % simulate spike train
    C = filter(1,[1 -P.gam],N);         % calcium concentration
    F = P.a*C+P.b + P.sig*randn(T,1);   % observations
    
    roiValues_deconv=[];
    
    for c=1:size(roiValues_norm,2)
        
        %fprintf('%d out of %d \n',c,size(roiValues_norm,2));
        
        %F=roiValues_norm(:,c)-mean(roiValues_norm(:,c));
        
        
        
        bp=roiValues_norm(:,c);
        pad_by=1000;
        bp = [bp(1:pad_by);bp;bp(end-pad_by+1:end)];
        f=normpdf([-1000:1000],0,400); f=f./sum(f);
        bp=conv(bp,f','same');
        bp=bp(pad_by:end-1-pad_by);
        
        %{
    clf; hold on;
    plot(roiValues_norm(:,c));
    plot(bp,'r');
        %}
        F=roiValues_norm(:,c)-bp;
        
        P.sig=1;%std(F);
        
        F=F./std(F);
        F(isnan(F))=0;
        % fast oopsi
        [Nhat Phat] = fast_oopsi(F,V,P);
        
        
        % plot results
        figure(1), clf; hold on;
        tvec=0:V.dt:(T-1)*V.dt;
        plot(tvec,F);
        plot(tvec,(Nhat*2)+6,'r');
        xlim([0 200]);
        
        drawnow;
        
        roiValues_deconv(:,c)=Nhat;
        
        
    end;
disp('spike extraction done');


%% SAVE
readInDirectory(1:end-11);

disp(readInDirectory(1:end-11))

x=input('warning saving here!! y?','s');
%manually cd to data dir
if (x == 'y')
    save([readInDirectory(1:end-11),'roivalues_new.mat'],'roiValues','roiValues_norm','roiValues_deconv');
    disp('saved');
end;


%%
axis('tight'), ylabel('C (au)')
h(3)=subplot(413); stem(tvec,N); hold on, plot(tvec,Nhat,'r','linewidth',1), axis('tight'), ylabel('fast')
Nsmc = M.nbar/max(M.nbar);
Nsmc(Nsmc<0.1)=0;
h(4)=subplot(414); stem(tvec,N); hold on, plot(tvec,Nsmc,'k','linewidth',2); axis('tight'), ylabel('smc')
xlabel('time (sec)')
linkaxes(h,'x')




%%

figure(1); clf; hold on;

for c=1:size(roiValues,2)
    a=zeros(1,123)';
    for i=1:799
        s=((i-1)*122)+1;
        a=a+(roiValues(s:s+122,c));
        
    end;
    a=a./799;
    plot(a/mean(a(:))+c*.1)
    drawnow;
    
    
end;


