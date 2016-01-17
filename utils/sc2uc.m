function [ xx, Ru ] = sc2uc( xxx, Rs, sc )
%sc2uc Returns a coordinate in unit cell plus cell translation for a coordinate in supercell
%   xx: coordinate in unit cell; Ru: cell translation
%   xxx: coordinate in super cell; Rs: cell translation; sc: supercell
%   definition
  tt=(xxx+Rs)*sc;
  Ru=floor(tt);
  xx=tt-Ru;
end

