
function XYZW_tens=scale_rule(XYZW_tens_ref,dbox)

%-----------------------------------------------------------------------
% Object:
%-----------------------------------------------------------------------
% Scale reference nodes from [-1,1]^3 to the bounding box.
%
% Note:
% For the needs of the routine it is not required to scale also the
% weights.
%-----------------------------------------------------------------------

a=dbox(1); b=dbox(2);
X=XYZW_tens_ref(:,1);
X=(a+b)/2 + ((b-a)/2)*X;

a=dbox(3); b=dbox(4);
Y=XYZW_tens_ref(:,2);
Y=(a+b)/2 + ((b-a)/2)*Y;

a=dbox(5); b=dbox(6);
Z=XYZW_tens_ref(:,3);
Z=(a+b)/2 + ((b-a)/2)*Z;

W=XYZW_tens_ref(:,4);

XYZW_tens=[X Y Z W];
