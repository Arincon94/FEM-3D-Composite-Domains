function KC = ElementStiffness_comp(D, lx, ly, lz)
    
    a = lx/2; b = ly/2; c = lz/2;

    syms x y z

    N1 = 1/8*(1 - x)*(1 - y)*(1 - z);
    N2 = 1/8*(1 + x)*(1 - y)*(1 - z);
    N3 = 1/8*(1 + x)*(1 + y)*(1 - z);
    N4 = 1/8*(1 - x)*(1 + y)*(1 - z);
    N5 = 1/8*(1 - x)*(1 - y)*(1 + z);
    N6 = 1/8*(1 + x)*(1 - y)*(1 + z);
    N7 = 1/8*(1 + x)*(1 + y)*(1 + z);
    N8 = 1/8*(1 - x)*(1 + y)*(1 + z);
    
    N1x = (1/a)*diff(N1, x); N1y = (1/b)*diff(N1,y); N1z = (1/c)*diff(N1,z);
    N2x = (1/a)*diff(N2, x); N2y = (1/b)*diff(N2,y); N2z = (1/c)*diff(N2,z);
    N3x = (1/a)*diff(N3, x); N3y = (1/b)*diff(N3,y); N3z = (1/c)*diff(N3,z);
    N4x = (1/a)*diff(N4, x); N4y = (1/b)*diff(N4,y); N4z = (1/c)*diff(N4,z);
    N5x = (1/a)*diff(N5, x); N5y = (1/b)*diff(N5,y); N5z = (1/c)*diff(N5,z);
    N6x = (1/a)*diff(N6, x); N6y = (1/b)*diff(N6,y); N6z = (1/c)*diff(N6,z);
    N7x = (1/a)*diff(N7, x); N7y = (1/b)*diff(N7,y); N7z = (1/c)*diff(N7,z);
    N8x = (1/a)*diff(N8, x); N8y = (1/b)*diff(N8,y); N8z = (1/c)*diff(N8,z);
    
    B1 = [N1x 0 0; 0 N1y 0; 0 0 N1z;
        N1y N1x 0; N1z 0 N1x; 0 N1z N1y];
    
    B2 = [N2x 0 0; 0 N2y 0; 0 0 N2z;
        N2y N2x 0; N2z 0 N2x; 0 N2z N2y];
    
    B3 = [N3x 0 0; 0 N3y 0; 0 0 N3z;
        N3y N3x 0; N3z 0 N3x; 0 N3z N3y];
    
    B4 = [N4x 0 0; 0 N4y 0; 0 0 N4z;
        N4y N4x 0; N4z 0 N4x; 0 N4z N4y];
    
    B5 = [N5x 0 0; 0 N5y 0; 0 0 N5z;
        N5y N5x 0; N5z 0 N5x; 0 N4z N5y];
    
    B6 = [N6x 0 0; 0 N6y 0; 0 0 N6z;
        N6y N6x 0; N6z 0 N6x; 0 N6z N6y];
    
    B7 = [N7x 0 0; 0 N7y 0; 0 0 N7z;
        N7y N7x 0; N7z 0 N7x; 0 N7z N7y];
    
    B8 = [N8x 0 0; 0 N8y 0; 0 0 N8z;
        N8y N8x 0; N8z 0 N8x; 0 N8z N8y];
    
    B = [B1 B2 B3 B4 B5 B6 B7 B8];
    
    KC1 = int(B'*D*B, x, -1,1); 
    KC2 = int(KC1,y,-1,1);
    KC = a*b*c*int(KC2,z,-1,1);

end