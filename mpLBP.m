function [mpLBPdes,punct_mpLBPdes]=mpLBP(v,f,R,nrad,multP,hfun,main_dir,normals)
%% Mean Point Local Binary Pattern (mpLBP)
%
% The mpLBP function compute the mpLBP descriptor from an input
% triangulation and a set of parameters. It first compute the punctual
% descriptor (or mpLBP punctual descriptor) on each vertex, then it compute
% the LBP and, finally, it computes the mpLBP descriptor, which is 
% basically an histogram of the punctual descriptors.
% 
% If you use this are using this code, PLEASE cite us. You can find the
% bibtex in the same docs_and_readmes folder. 
%
% Copyright (c) Moscoso~Thompson Elia, 2022
% -------------------------------------------------------------------------
%
% INPUT(those marked with (o) are optional)
% 
%           - v          : vertex list, a matrix (nv x 3) or (3 x nv), 
%                           where nv is the number of vertices;
%           - f          : face list, a matrix (nf x 3) or (3 x nf), where 
%                           nf is the number of faces;
%           - R          : the size of the maximum radius;
%           - multP      : the number of sectors of the first ring, which
%                           dictates the number of sectors of all the other
%                           rings based on the formula P_i = multP(2i-1),
%                           keeping the area of the sectors as constant as
%                           possible across rings;
%           - hfun       : the descriptor of the pattern, an arry with nv 
%                           scalar entries, relative to the vertex list. It
%                           comes in the form of a (1 x nv) or (nv x 1)
%                           array. Good candidates are curvatures.
%           - main_dir(o): a (nv x 3) or (3 x nv) array that indicates the
%                           first sector direction of each descriptor. In
%                           our work we proved that this does not change
%                           the performances significantly. However, for
%                           consistency, you might want to use it anyway.
%                           Good candidates are the principal directions. 
%                           Default is [1 0 0] for all the points. 
%           - normals(o): a (nv x 3) or (3 x nv) array that indicates the
%                           normal vectors to the surface in each vertex.
%                           If this input is not given, vertex normals are
%                           derived from the adjacent faces.
%
% OUTPUT
%           - mpLBPdes      : the mpLBP descriptor of the triangulation, a 
%                               feature vector of size 
%                               nrad x (multP(2nrad-1)+1)
%           - punct_mpLBPdes: the LBP descriptor computed on each vertex
                            


% Checking v and f vectors size
if size(f,2)~=3 
    f=f';
end
if size(v,2)~=3
    v=v';
end
nv=length(v);

% Checking main_dir size (or computing it)
if ~exist('main_dir','var') || isempty(main_dir)
    main_dir=repmat([1 0 0]',1,nv);
else
    if size(main_dir,1)~=3
        main_dir=main_dir';
    end
end


% Checking normals size (or computing it)
if exist('normals','var')
    if size(normals,1)~=3
        normals=normals';
    end
else
    normals = compute_normals(v,f)';
end

% Computing the number of sectors of each ring
P=multP*(2*(1:nrad)-1);

% Computing the neighborhood of each vertex.
[pp_idx,point_distance]=rangesearch(v,v,R); 

% Check if the models has a boundary and mark the vertices too close to the
% boundary. On those, the mpLBP punctual descriptor is not computed. 
v_boundary_idx=CheckBoundaryVertices(v,f,R);

% MAIN LOOP
% It checks which hfun was given and then... 


vlist=find(v_boundary_idx==1);
tpp_idx=pp_idx(vlist);
tpoint_distance=point_distance(vlist);
tdescriptor_data=cell(1,numel(vlist));
% ... it first compute the punctual descriptor ... 
for j=1:numel(vlist)
    tdescriptor_data{j}=meshpoint2flatdes(v,R,tpp_idx{j},tpoint_distance{j},...
        nrad,P,main_dir(:,j),hfun,normals(:,j));
end

% ... and then the LBP on the punctual descriptor. 
thfun=hfun(vlist);
punct_mpLBPdes=zeros(nv,nrad)-1;
for j=1:numel(vlist)
    punct_mpLBPdes(vlist(j),:)=computeLBPfromflatdes(tdescriptor_data{j},thfun(j),P);
end

% Computing the mpLBP descriptor of the surface
mpLBPdes=LBP2des(punct_mpLBPdes,nrad,P);
end

%% ------------------------------------------------------------------------
%% -------------------------- Major sub-functins --------------------------
%% ------------------------------------------------------------------------

%% MAIN LOOP
function outputstruct=meshpoint2flatdes(v,R,pp_idx,point_distance,nrad,P,main_dir,hfun,vnormal)
% This function compute the punctual descriptor of a vertex. 

empty_flags=0; % recount of sectors with no points (in v) in it. 
npp=numel(pp_idx);

% The following checks is the points are oriented pointing up. If they
% do not, they are rotated of 180 degrees.
help_rot=[null(vnormal'),vnormal];
[rp,Rrot,rot_main_dir]=regrplane((v(pp_idx,:)-repmat(v(pp_idx(1),:),npp,1))*help_rot,main_dir);
if (vnormal'*help_rot*Rrot)*([0,0,1]')<=0
    rp=rp*[1 0 0; 0 1 0; 0 0 -1];
end

% Projection on the xy-plane and rotation based on the main_dir vector.
p2d=rp(:,1:2);%<-- points on plane, considered in 2D.
cp2d=p2d-repmat(p2d(1,1:2),npp,1); %centered in the seed
vdir1=rot_main_dir(1:2)/norm(rot_main_dir(1:2));
vdir2=null(vdir1)';
crp2d=cp2d*[vdir1',vdir2'];


% Split points into sectors
drh=R/nrad;
dth=2*pi./P;
[th,~]=cart2pol(crp2d(:,1),crp2d(:,2));
th=th+pi; % fixing tetha interval ([-pi pi]->[0 2pi])
point_labels=cell(1,nrad);
for j=1:nrad
    point_labels{j}=cell(1,P(j));
end
for j=2:npp
    tlev=ceil((point_distance(j)+realmin)/drh); 
    tsec=ceil(th(j)/dth(tlev));
    point_labels{tlev}{tsec}=[point_labels{tlev}{tsec},j];
end


% Build Descriptor
scalar_descr=hfun(pp_idx);
des=cell(1,nrad);
for jR=1:nrad
    tmp_des_level=zeros(1,P(jR));
    for jP=1:P(jR)
        vlist=point_labels{jR}{jP};
        [xc,yc]=pol2cart(dth(jR)*(jP-.5),drh*(jR-.5));
        ws=ComputeGaussianWeights(crp2d(vlist,:),[xc,yc],...
            .5*drh,.3*drh*(jR-.5)*sin(dth(jR)/2),pi-dth(jR)*(jP-.5));
		tdescr_values=scalar_descr(vlist);
        if numel(vlist)==0
            tmp_des_level(jP)=0;
            empty_flags=empty_flags+1;
        else
            if sum(ws)==0
                tmp_des_level(jP)=0;
            else
                tmp_des_level(jP)=tdescr_values'*ws/sum(ws);
            end
        end
    end
    des{jR}=tmp_des_level;
end

% Saving output into a struct
outputstruct.des=des;
outputstruct.empty_flags=empty_flags;
outputstruct.pp_idx=pp_idx;
outputstruct.point_labels=point_labels;
outputstruct.rp=rp;
outputstruct.crp2d=crp2d;

end

%% Computing the LBP on the punctual descriptor
function LBPvalues=computeLBPfromflatdes(descr,hcenter,P)
% This function compute the Local Binary Pattern on the punctual descriptor

if numel(P)==1
    P=P*ones(numel(descr.des),1);
end
LBPvalues=zeros(1,numel(P));
LBPvect=cell(1,numel(P));
for j=1:numel(P)
    LBPvect{j}=zeros(P(j),1);
    for k=1:P(j)
        if descr.des{j}(k)>hcenter
            LBPvect{j}(k)=1;
        end
    end
    LBPvalues(j)=sum(LBPvect{j});
end

end
%% Computing the mpLBPdescriptor from the LBP values
function LBPdes=LBP2des(LBPvalues,nrad,P)
% This function compute the triangulation descriptor, which is basically an
% histogram of the LBP values obtained from the punctual descriptors.

normfact=0;
LBPdes=zeros(nrad,max(P+1));

maxPsize=max(P)+1;

for j=1:nrad
    LBPdes(j,P(j)+2:end)=-1*ones(1,maxPsize-P(j)-1);
end


for j=1:size(LBPvalues,1)
    if LBPvalues(j,1)~=-1
        for k=1:nrad
            LBPdes(k,LBPvalues(j,k)+1)=LBPdes(k,LBPvalues(j,k)+1)+1;
        end
        normfact=normfact+1;
    end
end

LBPdes=LBPdes/normfact;

for j=1:nrad
    for k=1:maxPsize
        if LBPdes(j,k)<0
            LBPdes(j,k)=-1;
        end
    end
end

end

%% ------------------------------------------------------------------------
%% -------------------------- Minor sub-functins --------------------------
%% ------------------------------------------------------------------------

%% Computation of vertices normal
function N=compute_normals(v,f)
% This function computes the vertices surface normal starting from the
% normal to the faces (computed in the standard "cross-product" way). The
% normal at each vertex v_i is the mean of the face normals which
% corresponding faces share v_i. 
% There are smarter ways to do it, feel free to outsmart this function in 
% any way you like it. 

N=zeros(size(v,1),3);
Nreg_contributions=zeros(size(v,1),1);
for j=1:size(f,1)
    
    vidx1=f(j,1);
    vidx2=f(j,2);
    vidx3=f(j,3);
    
    v1=v(vidx1,:);
    v2=v(vidx2,:);
    v3=v(vidx3,:);

    tnorm=cross(v2-v1,v3-v2);
    tnorm=tnorm/norm(tnorm);
    
    N(vidx1,:)=N(vidx1,:)+tnorm;
    Nreg_contributions(vidx1)=Nreg_contributions(vidx1)+1;
    N(vidx2,:)=N(vidx2,:)+tnorm;
    Nreg_contributions(vidx2)=Nreg_contributions(vidx2)+1;
    N(vidx3,:)=N(vidx3,:)+tnorm;
    Nreg_contributions(vidx3)=Nreg_contributions(vidx3)+1;
end
N=N./Nreg_contributions;
for j=1:size(N,1)
    N(j,:)=N(j,:)/norm(N(j,:));
end
end

%% Computation of gaussian weights
function ws=ComputeGaussianWeights(pp,center,sigx,sigy,cw_ang)
% This function computes the gaussian weights used in the estimation of a
% sector value.

a=(cos(cw_ang)^2)/(2*sigx^2)+(sin(cw_ang)^2)/(2*sigy^2);
b=(-sin(2*cw_ang))/(4*sigx^2)+(sin(2*cw_ang))/(4*sigy^2);
c=(sin(cw_ang)^2)/(2*sigx^2)+(cos(cw_ang)^2)/(2*sigy^2);

x0=center(1)+zeros(size(pp,1),1);
y0=center(2)+zeros(size(pp,1),1);
x=pp(:,1);
y=pp(:,2);

ws=exp(-(a*(x-x0).^2+2*b*(x-x0).*(y-y0)+c*(y-y0).^2));   
end

%% Regression plane computation (with least squares)
function [regrxyz,R,rot_main_dir]=regrplane(xyz,main_dir)
% This function rotates the points xyz so that they are aligned with the
% xy-plane,  multiple linear regression with least squares. The direction
% that indicates the first sector of a punctual descriptor (main_dir) is
% rotated accordingly.

x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);
X=[ones(size(x)) x y];
b=regress(z,X);
n=[b(2) b(3) -1];
v1=[1 0 b(2)];
v2=cross(v1,n);
R=[v1/norm(v1); v2/norm(v2); n/norm(n)]';

if det(R)~=-1  
    R=R(:,[2,1,3]);
end

regrxyz=xyz*R;
rot_main_dir=main_dir'*R;
end

%% Check if a boundary exists and if the vertices close to it
function v_tooclosetoboundary=CheckBoundaryVertices(v,f,R)
% Since the LBP does not work correctly in correspondance of boundaries, we
% first flag the vertices that are to close to them (if the model present
% some)

fadiac=tri_adiac_cell(v,f);
fidx_reg=ones(1,length(f)); % 0=border
for j=1:numel(fadiac)
    if numel(fadiac{j})==2
        fidx_reg(j)=0;
    end
end

borderface=f(fidx_reg==0,:);
vbidx=unique(borderface(:));
pp_idx=rangesearch(v,v(vbidx,:),R); %vc is here too

v_tooclosetoboundary=ones(1,length(v));
v_tooclosetoboundary(unique(horzcat(pp_idx{:})))=0;

end

%% Computing triangle-adjacency
function ad_tri=tri_adiac_cell(mat_ver,mat_tri)
% This function compute the triangle-triangle relation.

num_ver=size(mat_ver,1); num_tri=size(mat_tri,1);

succ=[2 3 1]; t=1:num_tri;

t_t=spalloc(num_ver,num_ver,3*num_tri);
for i=1:3
    t_t=t_t+sparse(mat_tri(:,i),mat_tri(:,succ(i)),t,num_ver,num_ver,num_tri);
end

for i=1:3
    index_edge=sub2ind(size(t_t),mat_tri(:,succ(i)),mat_tri(:,i));
    mat_tri(:,i+3)=t_t(index_edge);
end
mat_tri=mat_tri(:,4:6);
ad_tri=cell(1,size(mat_tri,1));

for i=1:size(mat_tri,1)
    % FIX n1
    bb=[];
    cc=0;
    for ii=1:3
        if mat_tri(i,ii)~=0
            cc=cc+1;
            bb(1,cc)=mat_tri(i,ii);
        end
    end
    cc=0;
    ad_tri{1,i}=bb;
end
end


