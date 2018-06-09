function h = plottext ( ax, str, x, y, x0, dx, dy, dw )    
%PLOTTEXT Show text with line on a plot 
%
%         H = PLOTTEXT ( AX, STR, X, Y, X0, DX, DY, DW )
%
%   Inputs:
%       AX   - handle to axes
%       STR  - string to show
%       X, Y - data vector arrays of the same length
%       XO   - point to mark
%       DX, DY, DW - changes in coordinates
%

%------------------------------------------------------
%       K. Yu. Polyakov         13 Jun 2004
%                               
%------------------------------------------------------
%   Find the nearest point
%------------------------------------------------------
        [xx,ind] = min(abs(x-x0)); 
        x0 = x(ind);
        y0 = y(ind);
%------------------------------------------------------
%   Transform coordinates to axis data
%------------------------------------------------------
        xLim = get(gca, 'XLim');
        yLim = get(gca, 'YLim');
        rangeX = xLim(2) - xLim(1);
        rangeY = yLim(2) - yLim(1);
        dx = dx*rangeX;
        dw = dw*rangeX;
        dy = dy*rangeY;
%------------------------------------------------------
%   Draw line
%------------------------------------------------------
        xLine = [x0 x0+dx x0+dx+dw];
        yLine = [y0 y0+dy y0+dy];
        line(xLine, yLine, 'Parent', ax);
%------------------------------------------------------
%   Draw text string
%------------------------------------------------------
        margin = 0.01;
        vAlign = 'middle';
        hAlign = 'left';
        if dw == 0
           xt = x0 + dx;
           yt = y0 + dy + sign(dy)*margin*rangeY;
           hAlign = 'center';
           if dy < 0
                vAlign = 'top'; 
           else vAlign = 'bottom'; end;
        else
           xt = x0 + dx + dw + sign(dw)*margin*rangeX;
           yt = y0 + dy;
           if dw < 0, hAlign = 'right'; end; 
        end;
        
        h = text(xt, yt, str, 'VerticalAlignment', vAlign, 'HorizontalAlignment', hAlign, 'Parent', ax );

    
%------ End of PLOTTEXT.M --------- KYuP ----------

