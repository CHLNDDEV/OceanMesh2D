function h=fastscatter(X,Y,C,varargin)
%% A fast scatter plot 
%  
% h=fastscatter(X,Y,C [,markertype,property-value pairs])
%
% Inputs: 
%    X,Y: coordinates 
%      C: color
%
%
% Examples:
%    N=100000;  
%    fastscatter(randn(N,1),randn(N,1),randn(N,1))
%  
%    N=100;  
%    fastscatter(randn(N,1),randn(N,1),randn(N,1),'+','markersize',7)
% 
%
% Aslak Grinsted 2014


marker='.';
if length(varargin)>0
	if length(varargin{1})==1
    	marker=varargin{1};varargin(1)=[];
	end
end

% %h/t to Boris Babic for the method. see http://www.mathworks.com/matlabcentral/newsreader/view_thread/22966
% h=mesh([X(:) X(:)]',[Y(:) Y(:)]',zeros(2,numel(X)),'mesh','column','marker',marker,'cdata',[C(:) C(:)]',varargin{:});
% view(2)


ix=find(~isnan(C+X+Y));
if mod(length(ix),2)==1
    ix(end+1)=ix(end);
end
ix=reshape(ix,2,[]);

h=mesh(X(ix),Y(ix),zeros(size(ix)),'marker',marker,'cdata',C(ix),'edgecolor','none','markeredgecolor','flat','facecolor','none',varargin{:});
view(2)

if nargout==0
    clear h
end