function anim_color = rgb2anim(rgb)
    % For an input RGB vector, output a decimal value representing the same color in Animatlab
    % Input: rgb: 1x3 vector of color value (0-255)
    % Output: anim_color: (float) decimal value of color
    rgb = floor(rgb);
    hexcolor =reshape(sprintf('%02X',rgb.'),6,[]).';
    hexcolor = ['FF',hexcolor];
    bitstring = char(hexToBinaryVector(hexcolor,32)+'0');
    anim_color = typecast(uint32(bin2dec(bitstring)),'int32');
end