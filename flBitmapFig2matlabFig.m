function [x1, y1, x2, y2] = flBitmapFig2matlabFig(imgfile, xlims, ylims, ncurves)



    % set(0, 'DefaultLineLineWidth', 2);
    
    
    
    %% LE PIXEIS DA IMAGEM (BITMAP)
    
    
    % imgfile = 'Fig_Airfoil_MomIn3.jpg'; % arquivo de imagem

    [I, map] = imread(imgfile);
    
    figure();
    imshow(I); % curvas: 1a azul e 2a vermelha
    % return
    
    
    
    % % remove/suaviza linhas horizontais
    % for r = 1 : length(I(:,1,1))
    %     for c = 1 : length(I(1,:,1))
    %         
    %         if I(r,c,1) == I(r,c,2) && I(r,c,1) == I(r,c,3)
    %             I(r,c,:) = 255;
    %         end
    %     end
    % end
    
    
    % for r = 1 : length(I(:,1,1))
    %     for c = 1 : length(I(1,:,1))
    %         
    %         if I(r,c,1) == I(r,c,2) && I(r,c,1) == I(r,c,3)
    %             I(r,c,:) = 255;
    %         end
    %         
    %         
    %         m = 1.01 * mean(I(r,c,:));
    %         
    %         if I(r,c,1) > m
    %             I(r,c,1) = 255;
    %         end
    %         
    % %         if I(r,c,2) > m
    % %             I(r,c,2) = 0;
    % %             I(r,c,3) = 255;
    % %         end
    %         
    %         if I(r,c,3) > m
    %             I(r,c,3) = 255;
    %         end
    %         
    %     end
    % end
    
    % figure();
    % imshow(I);
    % return;
    
    
    
    % Imagem em escala de cinza = 0.176R + 0.81G + 0.011B
    gs = double(0.176 * I(:,:,1) + 0.81 * I(:,:,2) + 0.011 * I(:,:,3));
    % % gs = (I(:,:,1));
    % figure();
    % imshow(gs);
    % return;
    
    % Mascara p/ deteccao de linhas horizontais
    % Resposta da mascara = somatoria(coef * nivel cinza)
    mask = [-1, -1, -1;
            +2, +2, +2;
            -1, -1, -1] * -1;
    
    % mask = [-3, -3, -3, -3, -3;
    %         +2, +2, +2, +2, +2;
    %         +2, +2, +2, +2, +2;
    %         +2, +2, +2, +2, +2;
    %         -3, -3, -3, -3, -3] * -1;
    
    resp = zeros(size(I(:,:,1)));% + 255;
    
    pc = (length(mask(:,1)) + 1) / 2; % ponto central
    pe = (length(mask(:,1)) - 1) / 2; % ponto extremo
    
    % for r = pc : (length(I(:,1,1)) - pe)
    %     for c = pc : (length(I(1,:,1)) - pe)
    %         resp(r,c) = sum(sum(double(gs(r-pe:r+pe, c-pe:c+pe)) .* mask));
    %     end
    % end
    
    
    
    % Mascaras p/ Sobel
    maskGx = [-1, -2, -1;
               0,  0,  0;
              +1, +2, +1];% * -1;
    
    maskGy = [-1,  0, +1;
              -2,  0, +2;
              -1,  0, +1];% * -1;
    
    respGx = zeros(size(I(:,:,1)));% + 255;
    respGy = zeros(size(I(:,:,1)));% + 255;
    respS = zeros(size(I(:,:,1)));% + 255;
    
    % respXr = zeros(size(I(:,:,1)));% + 255;
    % respYr = zeros(size(I(:,:,1)));% + 255;
    % respSr = zeros(size(I(:,:,1)));% + 255;
    % 
    % respXg = zeros(size(I(:,:,1)));% + 255;
    % respYg = zeros(size(I(:,:,1)));% + 255;
    % respSg = zeros(size(I(:,:,1)));% + 255;
    % 
    % respXb = zeros(size(I(:,:,1)));% + 255;
    % respYb = zeros(size(I(:,:,1)));% + 255;
    % respSb = zeros(size(I(:,:,1)));% + 255;
    
    R = I(:,:,1); % plano vermelho
    G = I(:,:,2); % plano verde
    B = I(:,:,3); % plano azul
    
    for r = pc : (length(I(:,1,1)) - pe)
        for c = pc : (length(I(1,:,1)) - pe)
            
            respGx(r,c) = sum(sum(double(gs(r-pe:r+pe, c-pe:c+pe)) .* maskGx));
            respGy(r,c) = sum(sum(double(gs(r-pe:r+pe, c-pe:c+pe)) .* maskGy));
            respS(r,c) = sqrt(respGx(r,c)^2 + respGy(r,c)^2);
            
    %         respXr(r,c) = sum(sum(double(R(r-pe:r+pe, c-pe:c+pe)) .* maskGx));
    %         respYr(r,c) = sum(sum(double(R(r-pe:r+pe, c-pe:c+pe)) .* maskGy));
    %         respSr(r,c) = sqrt(respXr(r,c)^2 + respYr(r,c)^2);
    %         
    %         respXg(r,c) = sum(sum(double(G(r-pe:r+pe, c-pe:c+pe)) .* maskGx));
    %         respYg(r,c) = sum(sum(double(G(r-pe:r+pe, c-pe:c+pe)) .* maskGy));
    %         respSg(r,c) = sqrt(respXg(r,c)^2 + respYg(r,c)^2);
    %         
    %         respXb(r,c) = sum(sum(double(B(r-pe:r+pe, c-pe:c+pe)) .* maskGx));
    %         respYb(r,c) = sum(sum(double(B(r-pe:r+pe, c-pe:c+pe)) .* maskGy));
    %         respSb(r,c) = sqrt(respXb(r,c)^2 + respYb(r,c)^2);
            
        end
    end
    
    % R(R < 150) = 0;
    % % G(G < ?) = 0;
    % B(B < 150) = 0;
    % 
    % figure();
    % imshow(R)%uint8(respSr));
    % 
    % % figure();
    % % imshow(G)%uint8(respSg));
    % 
    % figure();
    % imshow(B)%uint8(respSb));
    % return;
    
    
    r8 = uint8(respS);
    r8(r8 < 255) = 0;
    % r8(r8 == 255) = 1;
    
    % r8r = uint8(respSr);
    % r8r(r8r < 255) = 0;
    % 
    % r8g = uint8(respSg);
    % r8g(r8g < 255) = 0;
    % 
    % r8b = uint8(respSb);
    % r8b(r8b < 255) = 0;
    % 
    % figure();
    % imshow(r8r);
    % 
    % figure();
    % imshow(r8g);
    % 
    % figure();
    % imshow(r8b);
    % return;
    
    
    gr = (double(gs) .* double(r8));
    % figure();
    % imshow(uint8(gr));
    % return;
    
    
    % isso aqui melhora o zig-zag nos picos
    J(:,:,1) = I(:,:,1) - uint8(gr);
    J(:,:,2) = I(:,:,2) - uint8(gr);
    J(:,:,3) = I(:,:,3) - uint8(gr);
    
    I = J;
    % figure();
    % imshow(I);
    % return;
    
    
    
    limR = 10; % limiar vermelho
    % limG = 0; % limiar verde
    limB = 10; % limiar azul
    mrc = mean(I,3);
    
    for r = 1 : length(I(:,1,1))
        for c = 1 : length(I(1,:,1))
            
            % realca vermelho
            if abs(I(r,c,1) - mrc(r,c)) > limR ...
                    && (abs(I(r,c,2) - mrc(r,c)) <= limR ...
                    && abs(I(r,c,3) - mrc(r,c)) <= limR)
                
                I(r,c,1) = 255;
                I(r,c,[2,3]) = 0;
            end
            
    %         % realca verde
    %         if abs(I(r,c,2) - mrc(r,c)) > limG ...
    %                 && (abs(I(r,c,1) - mrc(r,c)) <= limG ...
    %                 && abs(I(r,c,3) - mrc(r,c)) <= limG)
    %             
    %             I(r,c,2) = 255;
    %             I(r,c,[1,3]) = 0;
    %         end
            
            % realca azul
            if abs(I(r,c,3) - mrc(r,c)) > limB ...
                    && (abs(I(r,c,1) - mrc(r,c)) <= limB ...
                    && abs(I(r,c,2) - mrc(r,c)) <= limB)
                
                I(r,c,3) = 255;
                I(r,c,[1,2]) = 0;
            end
            
        end
    end
    
    % figure();
    % imshow(I);
    % return;
    
    % re = edge(I(:,:,1), 'sobel', 'horizontal');
    % ge = edge(I(:,:,2), 'sobel', 'horizontal');
    % be = edge(I(:,:,3), 'sobel', 'horizontal');
    % 
    % I(:,:,1) = I(:,:,1) - uint8(~re);
    % I(:,:,2) = I(:,:,2) - uint8(~re);
    % I(:,:,3) = I(:,:,3) - uint8(~be);
    % figure();
    % imshow(I);
    % return
    
    % I = decorrstretch(I);
    % figure();
    % imshow(I);
    
    % I = decorrstretch(I, 'Tol', 0.1);
    % figure();
    % imshow(I);
    
    r = I(:,:,1);% .* uint8(~re); % plano vermelho
    g = I(:,:,2);% .* uint8(~ge); % plano verde
    b = I(:,:,3);% .* uint8(~be); % plano azul
    
    % figure(2);
    % imshow(r); % plano vermelho
    % 
    % figure(3);
    % imshow(g); % plano verde
    % 
    % figure(4);
    % imshow(b); % plano azul
    
    c1 = (b - g - r); % muda curva p/ branco (fundo da fig p/ preto)
    % figure(4); % 1a curva em preto (az na fig original)
    % imshow(c1);
    
    %separa 2a curva (vermelha na fig original)
    c2 = (r - g - b); % sub inverte / fundo preto
    % figure(5);
    % imshow(c2);
    
    c1 = edge(c1, 'sobel');%, 'vertical');
    c2 = edge(c2, 'sobel');%, 'vertical');
    % c1 = edge(c1, 'canny');
    % c2 = edge(c2, 'canny');
    % c1 = edge(c1, 'prewitt');
    % c2 = edge(c2, 'prewitt');
    % c1 = edge(c1, 'roberts');
    % c2 = edge(c2, 'roberts');
    % c1 = edge(c1, 'log');
    % c2 = edge(c2, 'log');
    % c1 = edge(c1, 'zerocross');
    % c2 = edge(c2, 'zerocross');
    
    
    % figure
    % imshow(edge(c1, 'sobel'))
    % 
    % figure
    % imshow(edge(c1, 'canny'))
    % 
    % figure
    % imshow(edge(c1, 'prewitt'))
    
    % return
    
    %% EXTRAI CURVAS DA FIGURA (PIXEIS COLORIDOS)
    
    
    % 1a curva
    vic1 = [];
    vir1 = [];
    mvir1 = [];
    mvic1 = [];
    
    for i = 1 : length(c1(1,:))
        
            ir = find(c1(:,i));
            ic = zeros(length(ir),1) + i;
    
            if ~isempty(ir)
                vir1 = [vir1; ir]; % y
                vic1 = [vic1; ic]; % x
                
                mvir1 = [mvir1; mean(ir)];
                mvic1 = [mvic1; i];
            end
            
    end
    
    
    % 2a curva
    vic2 = [];
    vir2 = [];
    mvir2 = [];
    mvic2 = [];
    
    for i = 1 : length(c2(1,:))
        
            ir = find(c2(:,i));
            ic = zeros(length(ir),1) + i;
    
            if ~isempty(ir)
                vir2 = [vir2; ir]; % y
                vic2 = [vic2; ic]; % x
                
                mvir2 = [mvir2; mean(ir)];
                mvic2 = [mvic2; i];
            end
            
    end
    
    
    
    %% DIMENSIONALIZA (PIXEIS P/ GRANDEZAS FISICAS)
    
    
    % 1a curva
    y1pico = max(mvir1); % linhas (contadas de cima p/ baixo)
    y1vale = min(mvir1);
    dy1 = y1pico - y1vale; % dif entre vale + baixo e pico + alto
    c1i = (y1pico + 1) - mvir1; % inverte ordem das linhas
    
    x1esq = min(mvic1); % colunas (contadas da esquerda p/ direita)
    x1dir = max(mvic1);
    dx1 = x1dir - x1esq; % dif entre inicio e fim da curva
    
    % fig escala real
    y1min = ylims(1,1);
    y1max = ylims(1,2);
    ddy1 = y1max - y1min;
    
    x1min = xlims(1,1);
    x1max = xlims(1,2);
    ddx1 = x1max - x1min;
    
    dX1 = ddx1 / dx1;
    dY1 = ddy1 / dy1;
    
    
    % 2nd curve
    if exist("ncurves", "var") && ncurves == 2
        y2pico = max(mvir2); % linhas (contadas de cima p/ baixo)
        y2vale = min(mvir2);
        dy2 = y2pico - y2vale; % dif entre vale + baixo e pico + alto
        c2i = (y2pico + 1) - mvir2; % inverte ordem das linhas
        
        x2esq = min(mvic2); % colunas (contadas da esquerda p/ direita)
        x2dir = max(mvic2);
        dx2 = x2dir - x2esq; % dif entre inicio e fim da curva
        
        % fig escala real
        y2min = ylims(2,1);
        y2max = ylims(2,2);
        ddy2 = y2max - y2min;
        
        x2min = xlims(2,1);
        x2max = xlims(2,2);
        ddx2 = x2max - x2min;
        
        dX2 = ddx2 / dx2;
        dY2 = ddy2 / dy2;
    end
    
    

    x1 = (mvic1 - x1esq) * dX1 + x1min;
    y1 = c1i * dY1 + y1min;

    if exist("ncurves", "var") && ncurves == 2
        x2 = (mvic2 - x2esq) * dX2 + x2min;
        y2 = c2i * dY2 + y2min;
    else
        x2 = [];
        y2 = [];
    end



%     outXY = [x1 y1 x2 y2];


    
    % figure(); % em pixels
    % hold all;
    % plot(mvic1, c1i); % , mvic2, c2i, 'LineWidth', 1);
    
    
    figure(); % em dimensoes fisicas
    hold all;
    plot(x1, y1);
    
    if exist("ncurves", "var") && ncurves == 2
        plot(x2, y2);
    end
    
    % set(gca, 'YTick', -100:50:100);
    % ylim([-100, 100]);
    grid off;
    box on;


end
