% Autor: Marx Freitas
% This function wrap-arounds the position of access points (APs) in the
% coverage area.

% Authors: Marx Freitas and Gilvan Soares Borges

% INPUT:
% x_axis_size       = the length in meters of the X-axis of the rectangular coverage area
% y_axis_size       = the length in meters of the Y-axis of the rectangular coverage area
% APpositions_m     = the position of APs
%                     Dim: nbrOfAPs x 1
% nbrOfAPs          = the number of APs in the newtork

% OUTPUT:
% APpositionsWrapped_m = the wrapped position of APs 
%                        Dim: nbrOfAPs x 9
% wrapLocations        = the center of each wrap arounded network
%                        Dim: 1 x 9

function [APpositionsWrapped_m, wrapLocations_m] = ...
    MF_generateSetupAPsWrapped(x_axis_size, y_axis_size, APpositions_m, ...
    nbrOfAPs)

% Auxiliary variables that help to compute delta X and delta Y
auxDeltaX = setdiff(real(APpositions_m),max(real(APpositions_m)));
auxDeltaY = setdiff(imag(APpositions_m),max(imag(APpositions_m)));

% These variables are used to corret numerical erros of some geometries
% such as 'line', for example
deltaX = abs(max(real(APpositions_m))-max(auxDeltaX));
deltaY = abs(max(imag(APpositions_m))-max(auxDeltaY));
delta = min(deltaX,deltaY);
deltaX = delta; deltaY = delta;

% This condition is utilized when the distribution of the APs follows a
% geometry that looks like a 'Line'
if isempty(delta)
    deltaX = abs(max(real(APpositions_m))-max(auxDeltaX));
    deltaY = deltaX;
end

% Compute the wrap arounded AP locations
wrapHorizontal = repmat([-(x_axis_size+deltaX) 0 (x_axis_size+deltaX)],[3 1]);
wrapVertical   = transpose(repmat([-(y_axis_size+deltaY) 0 (y_axis_size+deltaY)],[3 1]));
wrapLocations_m  = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped_m = repmat(APpositions_m,[1 length(wrapLocations_m)]) + repmat(wrapLocations_m,[nbrOfAPs 1]);

% %% Use the code linhes bellow, if you want to see the AP positions after wrapp around 
% figure;plot(real(APpositionsWrapped_m(:,1)),imag(APpositionsWrapped_m(:,1)),'o')
% hold all;
% plot(real(APpositionsWrapped_m(:,2)),imag(APpositionsWrapped_m(:,2)),'o')
% plot(real(APpositionsWrapped_m(:,3)),imag(APpositionsWrapped_m(:,3)),'o')
% plot(real(APpositionsWrapped_m(:,4)),imag(APpositionsWrapped_m(:,4)),'o')
% plot(real(APpositionsWrapped_m(:,5)),imag(APpositionsWrapped_m(:,5)),'x')
% plot(real(APpositionsWrapped_m(:,6)),imag(APpositionsWrapped_m(:,6)),'o')
% plot(real(APpositionsWrapped_m(:,7)),imag(APpositionsWrapped_m(:,7)),'o')
% plot(real(APpositionsWrapped_m(:,8)),imag(APpositionsWrapped_m(:,8)),'o')
% plot(real(APpositionsWrapped_m(:,9)),imag(APpositionsWrapped_m(:,9)),'o')

end