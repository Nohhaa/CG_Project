#if defined(UNICODE) && !defined(_UNICODE)
    #define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
    #define UNICODE
#endif

#include <tchar.h>
#include <windows.h>
#include <bits/stdc++.h>
#include <fstream>
#include <iostream>
#include <math.h>

using namespace std;


COLORREF Blue= RGB(213, 128, 244);
double cc=0.5;
int  A,B,z1,z2,z3,z4 ,xs,ys,xe,ye,k1 , k2 , k3 , k4,q,Ra,R1,R2;
vector <double> X;
vector <double> Y;
double x11, y11, x22, y22;


void savePoint(int x,int y)
{
    ofstream outFile("values.txt",std::ios::app);

    if (outFile.is_open()) {
        outFile << x << endl;
        outFile << y << endl;


        outFile.close();

    } else {
       cout << "Failed to open the file." <<endl;
    }
}
void save(string shape,COLORREF c )
{
    ofstream outFile("values.txt",std::ios::app);

    if (outFile.is_open()) {

        outFile << shape << endl;
         outFile << c << endl;


        outFile.close();

    } else {
       cout << "Failed to open the file." <<endl;
    }
}



void Draw4 (HDC hdc , int x , int y , int xc , int yc , COLORREF color )
{
    SetPixel(hdc , xc+x, yc+y , color);
    SetPixel(hdc , xc+x, yc-y , color);
    SetPixel(hdc , xc-x, yc+y , color);
    SetPixel(hdc , xc-x, yc-y , color);

}
void Draw8(HDC hdc , int xc , int yc , int x, int y , COLORREF c){
    SetPixel(hdc , xc+x , yc+y , c);
    SetPixel(hdc , xc+x , yc-y , c);
    SetPixel(hdc , xc-x , yc+y , c);
    SetPixel(hdc , xc-x , yc-y , c);
    SetPixel(hdc , xc+y , yc+x , c);
    SetPixel(hdc , xc+y , yc-x , c);
    SetPixel(hdc , xc-y , yc+x , c);
    SetPixel(hdc , xc-y , yc-x , c);
}

///top,bottom,left,right
int Round(double x){
    return (int)(x + 0.5);
}

///assignment helper

//To Draw a Line
void DrawLine(HDC hdc, int x1, int y1, int x2, int y2 , COLORREF c)
{
    int dx = x2 - x1;
    int dy = y2 - y1;
    double t = 0;
    double x;
    double y;
    SetPixel(hdc, x1, y1, c);
    for(t = 0 ; t <= 1 ; t += 0.001)
    {
         x = x1 +(t*dx);
         y = y1 + (t*dy);
        SetPixel(hdc, floor(x), floor(y), c);
    }

}



//To Draw a Circle
//Using Bresenham
void DrawCircleBresenham(HDC hdc , int xc , int yc , int R , COLORREF c){
    int x = 0 , y = R , d = 1-R;
    Draw8(hdc , xc , yc , x , y , c);
    while(x < y){
        if(d < 0){
            d += ((2*x)+3);
            x++;
        }
        else{
            d += ( 2*(x-y) + 5 );
            x++;
            y--;
        }
        Draw8(hdc , xc , yc , x , y , c);
    }
}


//DDAline algorithm
void LineDDA(HDC hdc,int xs,int ys,int xe,int ye,COLORREF color)
{
    int dx = xe - xs;
    int dy = ye - ys;
    if(abs(dx)>= abs(dy))
    {
        if(xs>xe)
        {
            swap(xs,xe);
            swap(ys,ye);
        }
        int x = xs; double y = ys;
        double m = (double)dy/dx;
        SetPixel(hdc,xs,ys,color);//draw first point of the line
        while(x<xe)
        {
            x++;y+=m;
            SetPixel(hdc,x,double(y),color);
        }
    }
    else
        {
            if(ys>ye)
            {
                swap(xs,xe);
                swap(ys,ye);
            }
            double x = xs;int y = ys;
            double m1 = (double)dx/dy;
            SetPixel(hdc,xs,ys,color);
            while(x<xe)
            {
                y++;x+=m1;
                SetPixel(hdc,double(x),y,color);
            }
        }

}

///Direct line algorithm
void DirectLine(HDC hdc,int xs,int ys,int xe,int ye,COLORREF color)
{

    int dx = xe - xs;
    int dy = ye - ys;
    double slope = (double)dy/dx;
    if(abs(dx)>=abs(dy))
    {
        if(xs>xe)
        {
            swap(xs,xe);
            swap(ys,ye);
        }
        for(int x = xs ;x <= xe; x++)
        {
            int y = round(ys+(x-xs)*slope);
            SetPixel(hdc,x,y,color);
        }
    }
    else
    {
        double isolpe = (double)dx/dy;
        if(ys>ye)
        {
            swap(xs,xe);
            swap(ys,ye);
        }
        for(int y=ys;y<=ye;y++)
        {
            int x=round(xs+(y-ys)*isolpe);
            SetPixel(hdc,x,y,color);
        }

    }
}

///midpoint line algorithm

void MidPointLine(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
    int x=x1, y=y1,d,d1,d2;
    double dx =(x2-x1), dy =(y2-y1);
    SetPixel(hdc, x, y, c);
    if ((dx==0 || dy/dx>1) && dy>0 && dx>=0)
    {
         d =(2*dx-dy),d1=(2*dx), d2=(2*dx-2*dy);
        while (y<y2)
        {
            if (d<=0)
            {
                d+= d1;
            }
            else
            {
                x++;
                d+= d2;
            }
            y++;
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy/dx>=0 && dy/dx<=1 && dy>=0 && dx>0)
    {
         d = (dx-2*dy), d1 = (-2*dy), d2 = (2*dx-2*dy);
        while (x < x2)
        {
            if (d>0)
            {
               d+= d1;
            }
            else
            {
                y++;
                d+= d2;
            }
             x++;
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy/dx<0 && dy/dx>=-1 && dy<=0 && dx>0)
    {
         d =(-dx-2*dy), d1=(-2*dy), d2 =(-2*dx-2*dy);
        while (x<x2)
        {
            if (d<=0)
            {
                d+= d1;
            }
            else
            {
                y--;
                d+= d2;
            }
                x++;
            SetPixel(hdc, x, y, c);
        }
    }
    else if ((dx==0 || dy/dx<-1) && dy<0 && dx>=0)
    {
         d =(-2*dx-dy), d1=(-2*dx), d2=(-2*dx-2*dy);
        while (y< y2)
        {
            if (d>0)
            {
                d+= d1;
            }
            else
            {
                x++;
                d+= d2;
            }
             y--;
            SetPixel(hdc, x, y, c);
        }
    }
    else if ((dx==0 || dy/dx>1) && dy<0 && dx<=0)
    {
         d = (-2*dx+dy), d1=(-2*dx), d2=(-2*dx+2*dy);
        while (y<y2)
        {
            if (d<=0)
            {
                d+= d1;
            }
            else
            {
                x--;
                d+= d2;
            }
            y--;
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy/dx>=0 && dy/dx<=1 && dy<=0 && dx<0)
    {
         d =(-dx+2*dy), d1=(2*dy), d2=(-2*dx+2*dy);
        while (x<x2)
        {
            if (d>0)
            {
                d += d1;
            }
            else
            {
                y--;
                d += d2;
            }
            x--;
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy/dx<0 && dy/dx>=-1 && dy>=0 && dx<0)
    {
         d=(dx+2*dy), d1=(2*dy), d2=(2*dx+2*dy);
        while (x<x2)
        {
            if (d<=0)
            {
                d+= d1;
            }
            else
            {
                y++;
                d+= d2;
            }
            x--;
            SetPixel(hdc, x, y, c);
        }
    }
    else if ((dx==0 || dy/dx<-1) && dy>0 && dx<=0)
    {
         d = (2*dx+dy), d1 =(2*dx), d2 =(2*dx+2*dy);
        while (y<y2)
        {
            if (d>0)
            {
                d+= d1;
            }
            else
            {
                x--;
                d+= d2;
            }
             y++;
            SetPixel(hdc, x, y, c);
        }
    }
}


///Circle Direct Method
void CircleDirect(HDC hdc, int xc, int yc, int R, COLORREF c)
{
    double x = 0 , y = R;
    Draw8(hdc, xc, yc, R, 0, c);
    while (x < y)
    {
        x++;
        y = sqrt(R * R - x * x);
        Draw8(hdc, xc, yc, x, y, c);
    }

}

///circle polar
void CirclePolar(HDC hdc,int xc,int yc, int R,COLORREF color)
{
    int x = R , y = 0;
    double dtheta=1.0/R;
    Draw8(hdc,xc,yc,x,y,color);
    for(double th=0; th <6.28; th+=dtheta)
    {
        double x = R * cos(th);
        double y = R * sin(th);
        Draw8(hdc,xc,yc,Round(x),Round(y),color);
    }
}

///Iterative polar circle
void CircleIterativePolar(HDC hdc,int xc,int yc,int R,COLORREF color)
{
    double x = R , y =0.0;
    double dtheta = 1.0/R;
    double c = cos(dtheta);
    double s = sin(dtheta);
    Draw8(hdc,xc,yc,x,y,color);
    while(x>y)
    {
        double x1 = x * c - y * s;
        y = x * s + y * c;
        x = x1;
        Draw8(hdc,xc,yc,Round(x),Round(y),color);
    }
}

///midpoint circle
void circle_mid(HDC hdc, int r, int x11, int y11,COLORREF color)
{
	int x = 0, y = r, d = 1 - r;
	Draw8(hdc, x11, y11, x, y,color);
	while (x < y)
	{

		if (d <= 0)
		{
			x++;
			d = d + (2 * x) + 3;


		}
		else
		{

			x++;
			y--;
			d = d + (2 * (x - y)) + 5;
		}
		Draw8(hdc, x11, y11, x, y,color);
	}


}


///Circle Modified Midpoint Method
void CircleModifiedMid(HDC hdc , int xc , int yc , int R , COLORREF c)
{
    int x=0 , y=R;
    int d = 1 - R;
    int d1 = 3 , d2 = 5 - (2*R) ;
    Draw8(hdc , xc , yc , x , y , c);
    while(x < y){
        if(d < 0){
            d += d1;
            d1 += 2;
            d2 += 2;
            x++;
        }
        else{
            d += d2;
            d2 += 4;
            d1 += 2;
            x++;
            y--;
        }
        Draw8(hdc , xc , yc , x , y , c);
    }
}



///FILLLING

///which quadrant and Filling by circles
void DrawQuadrant_circle(HDC hdc, int xc, int yc, int x, int y, int quadrant, COLORREF c)
{
    if (quadrant == 1) {
        SetPixel(hdc, xc + x , yc - y , c);
        SetPixel(hdc, xc + y , yc - x , c);
    }
    else if (quadrant == 2) {
        SetPixel(hdc, xc - y , yc - x , c);
        SetPixel(hdc, xc - x , yc - y , c);
    }
    else if (quadrant == 3) {
        SetPixel(hdc, xc - x , yc + y , c);
        SetPixel(hdc, xc - y , yc + x , c);
    }
    else if(quadrant == 4){
        SetPixel(hdc, xc + x , yc + y , c);
        SetPixel(hdc, xc + y , yc + x , c);
    }
}

///which quadrant and Filling by lines

void DrawQuadrant_line(HDC hdc, int xc, int yc, int x, int y, int quadrant, COLORREF c)
{
    if (quadrant == 1) {
        DrawLine(hdc, xc , yc , xc + x, yc - y, c);
        DrawLine(hdc, xc , yc , xc + y, yc - x, c);
    }
    else if (quadrant == 2) {

        DrawLine(hdc, xc , yc ,  xc - y, yc - x, c);
        DrawLine(hdc, xc , yc ,  xc - x, yc - y, c);
    }
    else if (quadrant == 3) {
        DrawLine(hdc, xc , yc ,  xc - x, yc + y, c);
        DrawLine(hdc, xc , yc ,  xc - y, yc + x, c);
    }
    else if(quadrant == 4){
        DrawLine(hdc, xc , yc ,  xc + x, yc + y, c);
        DrawLine(hdc, xc , yc ,  xc + y, yc + x, c);
    }
}

///filling circle with lines
void FillingQuadrantByLines(HDC hdc, int xc, int yc, int R, int quadrant, COLORREF c)
{
    double x = R, y = 0;
    double dtheta = 1.0 / R;
    double st = sin(dtheta);
    double ct = cos(dtheta);
    while(x > y)
    {
        double x1 = x * ct - y * st;
        y = x * st + y * ct;
        x = x1;
        Draw8(hdc, xc, yc, Round(x), Round(y), c);
        DrawQuadrant_line(hdc, xc, yc, Round(x), Round(y), quadrant , c);
    }
}


///filling circle with circle
void FillingQuadrantByCircles(HDC hdc, int xc, int yc, int R, int quadrant, COLORREF c)
{
    double x = R, y = 0;
    double dtheta = 1.0 / R;
    double st = sin(dtheta);
    double ct = cos(dtheta);
    while (x > y)
    {
        double x1 = x * ct - y * st;
        y = x * st + y * ct;
        x = x1;
        Draw8(hdc, xc, yc, Round(x), Round(y), c);
    }
    while (R) {
        dtheta = 1.0 / R;
        st = sin(dtheta);
        ct = cos(dtheta);
        x = R;
        y = 0;
        while (x > y)
        {
            double x1 = x * ct - y * st;
            y = x * st + y * ct;
            x = x1;
            DrawQuadrant_circle(hdc, xc, yc, Round(x), Round(y), quadrant, c);
        }
        R--;

    }
}

///Filling Rectangle with Bezier
void DrawBezierCurve(HDC hdc, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4 , COLORREF c)
{
    int a0 = x1;
    int a1 = 3 * (x2 - x1);
    int a2 = 3 * x1 - 6 * x2 + 3 * x3;
    int a3 = -1 * x1 + 3 * x2 - 3 * x3 + x4;
    int b0 = y1;
    int b1 = 3 * (y2 - y1);
    int b2 = 3 * y1 - 6 * y2 + 3 * y3;
    int b3 = -1 * y1 + 3 * y2 - 3 * y3 + y4;

    for (double t = 0; t <= 1; t += 0.001)
    {
        double xt = (t * t * t * a3) + (t * t * a2) + (t * a1) + a0;
        double yt = (t * t * t * b3) + (t * t * b2) + (t * b1) + b0;
        SetPixel(hdc, xt, yt , c);

    }
}


void FillingRectangleWithBezier(HDC hdc, int rectX1 , int rectX2 , int rectY1 ,int rectY2, int y1 , int y2,int y3 , int y4,int x1 , int x2, int x3, int x4, COLORREF c)
{
    if(rectX1 > rectX2)
    {
        swap(rectX1 , rectX2);
    }

    while(x2>=rectX1 &&x2<=rectX2-15 && x3>=rectX1 && x3<=rectX2 )
    {
        x1++;
        x2++;
        x3++;
        x4++;
        DrawBezierCurve(hdc, x1, rectY1, x2 , y2 , x3, y3, x4, rectY2 , c);
    }
}




///Filling square with Hermite curve
void DrawHermiteCurve(HDC hdc, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4 , COLORREF c)
{
      int x , y;
    int a1 , a2 , a3 , a4 , b1 , b2 , b3 , b4;
    a4 = x1;
    a3 = x2;
    b4 = y1;
    b3 = y2;
    a2 = 3*x4 - 3*x1 - 2*x2 - x3;
    b2 = 3*y4 - 3*y1 - 2*y2 -  y3;
    a1 = -2*x4 + 2*x1 + x2 + x3;
    b1 = -2*y4 + 2*y1 + y2 + y3;

    for(double t = 0 ; t <= 1 ; t += 0.0001)
    {

        x = (a1 * t*t*t )+ (a2*t*t) + (a3 * t )+ a4;
        y = (b1 * t*t*t )+ (b2*t*t) + (b3 * t) + b4;
        SetPixel(hdc , x , y , c);
        //LineTo(hdc , x , y );
    }
}


void FillingSquareWithHermite(HDC hdc,int rectX1 , int rectX2 , int rectY1 ,int rectY2, int y1 , int y2,int y3 , int y4,int x1 , int x2, int x3, int x4, COLORREF c)
{

if(rectX1 > rectX2)
    {
        swap(rectX1 , rectX2);
    }

    while(x2>=rectX1 &&x2<=rectX2 )
    {
        x1++;
        x2++;
        x3++;
        x4++;
        DrawHermiteCurve(hdc, x1, rectY1, x2 , y2 , x3, y3, x4, rectY2 , c);
    }
}


///Flood fill recursive
void RECFloodFill(HDC hdc,int x,int y,COLORREF bc,COLORREF fc)
{
    COLORREF c = GetPixel(hdc , x , y);
    if(c == bc || c == fc){
        return;
    }
    SetPixel(hdc , x , y , fc);
    RECFloodFill(hdc , x+1 , y , bc , fc);
    RECFloodFill(hdc , x-1 , y , bc , fc);
    RECFloodFill(hdc , x , y-1 , bc , fc);
    RECFloodFill(hdc , x , y+1 , bc , fc);

}

///Flood fill non-recursive
class point{
public:
    int x , y;
    point(int x=0 , int y=0){
        this->x = x;
        this->y = y;
    }
};
void NonRecursiveFloodFill(HDC hdc , int x , int y , COLORREF bc , COLORREF fc)
{
    stack<point>stk;
    stk.push( point(x , y) );
    while( !stk.empty() ){
        point p = stk.top();
        stk.pop();
        COLORREF c = GetPixel(hdc , p.x , p.y);
        if(c == bc || c == fc){
            continue;
        }
        else{
        SetPixel(hdc , p.x , p.y , fc);
        stk.push( point(p.x-1 , p.y) );
        stk.push( point(p.x+1 , p.y) );
        stk.push( point(p.x , p.y-1) );
        stk.push( point(p.x , p.y+1) );
    }}
}


///Ellipse Polar Method
void EllipsePolar(HDC hdc ,int A, int B, int xc, int yc, COLORREF c)
{

    double theta = 0, dtheta = 1.0 / max(A, B);
    double x = 0, y = B;
    while (theta <= 2 * M_PI)
    {
        Draw4(hdc, xc, yc, Round(x), Round(y), c);

        x = A * cos(theta);
        y = B * sin(theta);
        theta += dtheta;
    }


}

///ellipse direct
 void  Direct_Ellipse (HDC hdc ,double a, int b,int xc, int yc , COLORREF c)
 {
     double x =a;
     int y=0;
     Draw4(hdc,x,y,xc,yc,c);

     while(x>=0)
     {
         x-=0.05;
         y = abs(b*sqrt(1-(double)(x*x)/(a*a))) ;
         Draw4(hdc,x,y,xc,yc,c);

     }
 }

 ///ellipse midpoint
 void Mid_Ellipse(HDC hdc , int a, int b,int xc, int yc , COLORREF c)
{
    double dx, dy, d1, d2, x, y;
    x=0;
    y=b;
Draw4(hdc,x,y,xc,yc,c);

    d1=(b*b)-(a*a*b)+(0.25*a*a);
    dx=2*b*b*x;
    dy=2*a*a*y;


    while (dx<dy)
    {
        if (d1<0)
        {
            dx+=(2*b*b);
            d1+=dx+(b*b);
        }
        else
        {
            y--;
            dx+=(2*b*b);
            dy-=(2*a*a);
            d1+=dx-dy+(b *b);
        }
         x++;
        Draw4(hdc,x,y,xc,yc,c);
    }


    d2 = ((b*b)*((x+0.5)*(x+0.5)))+((a*a)*((y-1)*(y-1)))-(a*a*b*b);
    while (y>=0)
    {
        if (d2>0)
        {

            dy-=(2*a*a);
            d2+=(a*a)-dy;
        }
        else
        {

            x++;
            dx+=(2*b*b);
            dy-=(2*a*a);
            d2+=dx-dy+(a*a);
        }
         y--;
        Draw4(hdc,x,y,xc,yc,c);
    }

}

////Convex filling
typedef struct {
    int xleft;
    int xright;
} Table[1000];


void swap(POINT& p1, POINT& p2) {
    POINT temp = p1;
    p1 = p2;
    p2 = temp;
}

void init_table(Table t) {
    for (int i = 0; i < 1000; i++) {
        t[i].xleft = INT_MAX;
        t[i].xright = -INT_MAX;
    }
}
void convert_edges_to_table(POINT p1, POINT p2, Table t) {
    if (p1.y == p2.y)
        return;
    if (p1.y > p2.y)
        swap(p1, p2);
    double x = p1.x;
    int y = p1.y;
    double dx = p2.x - p1.x;
    int dy = p2.y - p1.y;
    double m = (double)dx / dy;
    while (y < p2.y) {
        if (x < t[y].xleft)
            t[y].xleft = (int)ceil(x);
        if (x > t[y].xright)
            t[y].xright = (int)floor(x);
        y++;
        x += m;
    }
}
void convert_polygon_to_table(POINT polygon[], int size, Table t) {
    POINT v1 = polygon[size - 1];
    for (int i = 0; i < size; i++) {
        POINT v2 = polygon[i];
        convert_edges_to_table(v1, v2, t);
        v1 = v2;
    }
}
void table_to_screen(HDC hdc, Table t, COLORREF color) {
    for (int i = 0; i < 1000; i++) {
        if (t[i].xleft < t[i].xright) {

            for (int x = t[i].xleft; x <= t[i].xright; x++) {
                SetPixel(hdc, x, i, color);
            }
        }
    }
}

void convexFilling(HDC hdc, POINT polygon[], int n, COLORREF c) {
    Table t;
    init_table(t);
    convert_polygon_to_table(polygon, n, t);
    table_to_screen(hdc, t, c);
}

///Non-convex Filling
vector<POINT> PO;
vector<COLORREF> C;

struct Table2
{
    double x;
    double sl;
    int ymax;
    bool operator < (Table2 r) const
    {
        return x < r.x;
    }
};
typedef list <Table2> Edge;

Table2 InitEdge(POINT& p1, POINT& p2)
{
    if (p1.y > p2.y )
        swap(p1, p2);
    Table2 r{};
    r.x = p1.x;
    r.ymax = p2.y;
    r.sl = (double)(p2.x - p1.x) / (p2.y - p1.y);
    return r;
}

void InitTable(POINT* shape, int n, Edge table[])
{
    POINT p1 = shape[n - 1];
    for (int i = 0; i < n; i++)
    {
        POINT p2 = shape[i];
        if (p1.y == p2.y)
        {
            p1 = p2;
        }
        Table2 r = InitEdge(p1, p2);
        int tmp = p1.y;
        table[tmp].push_back(r);
        p1 = shape[i];
    }
}

void NonConvexFill(HDC hdc, POINT* shape, int n, COLORREF c)
{
    auto* table = new Edge[600];
    InitTable(shape, n, table);
    int y = 0;
    while (y < 600 && table[y].empty())
        y++;
    if (y == 600)
        return;
    Edge Active = table[y];
    while (!Active.empty())
    {
        Active.sort();
        for (auto it = Active.begin(); it != Active.end(); it++)
        {
            int x1 = ceil(it->x);
            it++;
            int x2 = floor(it->x);
            for (int x = x1; x <= x2; x++) {
                SetPixel(hdc, x, y, c);
                POINT p;
                p.x = x;
                p.y = y;
                PO.push_back(p);
                C.push_back(c);
            }
        }
        y++;
        auto it = Active.begin();
        while (it != Active.end())
            if (y == it->ymax) it = Active.erase(it);
            else it++;
        for (auto & itt : Active)
            itt.x += itt.sl;
        Active.insert(Active.end(), table[y].begin(), table[y].end());
    }
    delete[] table;
}




///cardinalSpline curve
void cardinalSpline(HDC hdc ,int c, int x1 , int y1 , int x2 , int y2 , int x3 , int y3 , int x4 , int y4, COLORREF C)
{
    int x , y,s0x,s0y,s1x,s1y;
    int a1 , a2 , a3 , a4 , b1 , b2 , b3 , b4;
    s0x=(1-c)*(x2-x1);
    s0y=(1-c)*(y2-y1);
    s1x=(1-c)*(x4-x3);
    s1y=(1-c)*(y4-y3);

    a4 = x1;
    a3 = s0x;
    b4 = y1;
    b3 = s0y;
    a2 = 3*x4 - 3*x1 - 2*s0x - s1x;
    b2 = 3*y4 - 3*y1 - 2*s0y -  s1y;
    a1 = -2*x4 + 2*x1 + s0x + s1x;
    b1 = -2*y4 + 2*y1 + s0y + s1y;

    for(double t = 0 ; t <= 1 ; t += 0.0001)
    {

        x = (a1 * t*t*t )+ (a2*t*t) + (a3 * t )+ a4;
        y = (b1 * t*t*t )+ (b2*t*t) + (b3 * t) + b4;
        SetPixel(hdc , x , y , C);

    }

}

///clipping


///point clipping
void PointClipping(HDC hdc,int x,int y,int left,int top,int right,int bottom,COLORREF color)
{
 if(x>=left && x<= right && y>=top && y<=bottom)
  {
      SetPixel(hdc,x,y,color);
  }
}

int num, num2,ch,s=10,out1[4]={0,0,0,0},out2[4]={0,0,0,0};
void vertical_intersect(double xs,double ys,double xe,double ye,int x,double* xx,double* yy )
{
    *xx=x;
    *yy=ys+(x-xs)*(ye-ys)/(xe-xs);


}
void horizontal_intersect(double xs,double ys,double xe,double ye,int y,double* xx,double* yy )
{
    *yy=y;
    *xx=xs+(y-ys)*(xe-xs)/(ye-ys);

}

bool out(double X,double Y,double e,int l,int r,int t,int b)
{
    if(e==l){
        return X<e;
    }
    else if(e==r)
    {
        return X>e;
    }
    else if(e==t)
    {
       return Y<e;
    }
    else if(e==b)
    {
        return Y>e;
    }
}

void ClipWithEdge(vector <double> Xpoint,vector <double> Ypoint,int edge,int l,int r, int t,int b,vector <double>& XOut,vector <double>& YOut)
{

int n=Xpoint.size();

double X_v1=Xpoint[n-1];
double Y_v1=Ypoint[n-1];

bool v1_out=out(X_v1,Y_v1,edge,l,r,t,b);
int c=0;

for(int i=0;i<n;i++)
{
double X_v2=Xpoint[i];
double Y_v2=Ypoint[i];


bool v2_out=out(X_v2,Y_v2,edge,l,r,t,b);

double xx,yy;
    if(edge==l || edge==r)
    {
        vertical_intersect(X_v1,Y_v1,X_v2,Y_v2,edge,&xx,&yy);
    }
    else
    {
        horizontal_intersect(X_v1,Y_v1,X_v2,Y_v2,edge,&xx,&yy);
    }


if( (v1_out )&& (!v2_out))
{
    XOut.push_back(xx);
   YOut.push_back(yy);



 XOut.push_back(X_v2);
   YOut.push_back(Y_v2);

}
else if( (!v1_out )&& (!v2_out))
{
    XOut.push_back(X_v2);
YOut.push_back(Y_v2);


}
else if( !v1_out)
{
    XOut.push_back(xx);
   YOut.push_back(yy);
}

X_v1=X_v2;
Y_v1=Y_v2;
v1_out= v2_out ;


}


}


///polygon clipping
void Polygon_Clip(HDC hdc,vector <double>& Xpoint,vector <double>& Ypoint,int l,int r, int t,int b)
{

vector <double> Xlist;
vector <double> Ylist;


ClipWithEdge(Xpoint,Ypoint,l,l,r,t,b,Xlist,Ylist);
Xpoint.clear();
Ypoint.clear();
Xpoint.assign(Xlist.begin(), Xlist.end());
Ypoint.assign(Ylist.begin(), Ylist.end());

Xlist.clear();
Ylist.clear();

ClipWithEdge(Xpoint,Ypoint,t,l,r,t,b,Xlist,Ylist);
Xpoint.clear();
Ypoint.clear();
Xpoint.assign(Xlist.begin(), Xlist.end());
Ypoint.assign(Ylist.begin(), Ylist.end());

Xlist.clear();
Ylist.clear();
ClipWithEdge(Xpoint,Ypoint,r,l,r,t,b,Xlist,Ylist);
Xpoint.clear();
Ypoint.clear();
Xpoint.assign(Xlist.begin(), Xlist.end());
Ypoint.assign(Ylist.begin(), Ylist.end());

Xlist.clear();
Ylist.clear();
ClipWithEdge(Xpoint,Ypoint,b,l,r,t,b,Xlist,Ylist);



double X_v1=Xlist[Xlist.size()-1];
double Y_v1=Ylist[Ylist.size()-1];

for(int i=0;i<Xlist.size();i++)
{
double X_v2=Xlist[i];
double Y_v2=Ylist[i];
MoveToEx(hdc,round(X_v1),round(Y_v1),NULL);
LineTo(hdc,round(X_v2),round(Y_v2));
X_v1=X_v2;
Y_v1=Y_v2;
}
}

int getoutones(int out[4])
{
    int cc=0;
     for(int i=0; i<4;i++){

     if(out[i]==1){
        cc++;
     }}
 return cc;
}
void getout(double x, double y, int left, int right, int top, int bottom, int out[4])
{

out[0] = 0;
out[1] = 0;
out[2] = 0;
out[3] = 0;
    if (x < left) {
       out[2] = 1;
    }
    if (x > right) {
        out[3] = 1;
    }
    if (y < top) {
        out[0] = 1;
    }
    if (y > bottom) {
        out[1] = 1;
    }
}


///line clipping
void LINE_Clipping(HDC hdc,double xs,double ys,double xe,double ye,int left,int right,int top,int bottom)
{
int    o1[4]={0,0,0,0},o2[4]={0,0,0,0};
    getout( xs, ys, left, right, top, bottom, o1);

   getout( xe, ye, left, right, top, bottom, o2);


    int c=getoutones(o1) ,c1=getoutones(o2);
    double xss=xs,yss=ys,xee=xe,yee=ye,xx=0,yy=0;


    while(c>0)
   {

        if(o1[0]==1)
        {
         horizontal_intersect( xss, yss, xee, yee, top, &xx, &yy );
         xss=xx;
         yss=yy;



        }
         else if(o1[1]==1)
        {

        horizontal_intersect( xss, yss, xee, yee, bottom, &xx, &yy );
        xss=xx;
         yss=yy;



        }
        else if(o1[2]==1)
        {

         vertical_intersect( xss, yss, xee, yee, left, &xx, &yy );
         xss=xx;
         yss=yy;


        }
        else
        {

         vertical_intersect( xss, yss, xee, yee, right, &xx, &yy );
         xss=xx;
         yss=yy;



        }
         getout( xss, yss, left, right, top, bottom, o1);
          c=getoutones(o1);


    }

 while(c1>0)
   {


        if(o2[0]==1)
        {
         horizontal_intersect( xss, yss, xee, yee, top, &xx, &yy );
         xee=xx;
         yee=yy;


        }
        if(o2[1]==1)
        {

        horizontal_intersect( xss, yss, xee, yee, bottom, &xx, &yy );
        xee=xx;
        yee=yy;

        }
        if(o2[2]==1)
        {

         vertical_intersect( xss, yss, xee, yee, left, &xx, &yy );
         xee=xx;
         yee=yy;

        }
        if(o2[3]==1)
        {

         vertical_intersect( xss, yss, xee, yee, right, &xx, &yy );
         xee=xx;
         yee=yy;

        }
 getout( xee, yee, left, right, top, bottom, o2);
c1=getoutones(o2);

    }



        MoveToEx(hdc,round(xss),round(yss),NULL);
        LineTo(hdc,round(xee),round(yee));

}

///////////////////////////////////////////////////////////////////
/// EXTRA TASK:

bool Intersecting(int r1,int r2, int xc1, int yc1, int xc2, int yc2)
{
    float d = sqrt(pow(abs(xc2 - xc1), 2) + pow(abs(yc2 - yc1), 2));
    if (d <= r1 + r2 && d >= abs(r1 - r2))
    {
        return true;
    }
    else {
        return false;
    }
}


void FillingIntersection(HDC hdc , int r1 , int r2, int xc1 , int yc1  , int xc2 , int yc2 , COLORREF c)
{
    if (Intersecting(r1 , r2, xc1 , yc1  , xc2 , yc2))
    {

        cout << "Intersecting!\n";


        float d = sqrt(pow(xc2 - xc1, 2) + pow(yc2 - yc1, 2));
        float a = (pow(r1 , 2) - pow(r2 , 2) + pow(d, 2)) / (2 * d);
        float h = sqrt(pow(r1 , 2) - pow(a , 2));
        float x3 = xc1 + a * (xc2 - xc1) / d;
        float y3 = yc1 + a * (yc2 - yc1) / d;

        float x4 = x3 + h * (yc2 - yc1) / d;
        float y4 = y3 - h * (xc2 - xc1) / d;

        float x5 = x3 - h * (yc2 - yc1) / d;
        float y5 = y3 + h * (xc2 - xc1) / d;


        int xMin = ceil(min(x4, x5));
        int xMax = floor(max(x4, x5));
        int yMin = ceil(min(y4, y5));
        int yMax = floor(max(y4, y5));

        for (int y = yMin; y <= yMax; y++)
        {
           bool in = false;
           for (int x = xMin ; x <= xMax+550; x++)
           {
               if (sqrt( pow(x - xc1, 2) + pow(y - yc1, 2) ) <= r1 && sqrt( pow(x - xc2, 2) + pow(y - yc2, 2) ) <= r2) {
                   SetPixel(hdc, x , y , c);
                   in = true;
               }
               else{
                   in = false;
               }
           }

           for (int x = xMax ; x >= xMin-550 ; x--)
           {
               if (sqrt( pow(x - xc1, 2) + pow(y - yc1, 2) ) <= r1 && sqrt( pow(x - xc2, 2) + pow(y - yc2, 2) ) <= r2) {
                   SetPixel(hdc, x , y , c);
                   in = true;

               }
               else{
                   in = false;
               }
           }
       }
    }
    else {
        cout << "NOT Intersecting!\n";
    }
}

void Load(HWND hwnd,HDC hdc )
{
    ifstream inFile("values.txt");
     if (inFile.is_open()) {
        int c,x,y,xx,yy,x1,y1,x2,y2,x3,y3,l,t,b,r,xc1,xc2,yc1,yc2,R1,R2,Ra,Quart;
        vector <double> X;
        vector <double> Y;
        static POINT P[5];
        string shape;
        COLORREF color;

        while (inFile >> shape) {

inFile >> c;
color=c;
            if(shape=="s1")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                LineDDA(hdc, x, y, xx, yy, color);
            }

           else if(shape=="s2")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                DirectLine(hdc, x, y, xx, yy, color);
            }
           else if(shape=="s3")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                 MidPointLine(hdc, x, y, xx, yy, color);

            }
           else if(shape=="s4")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                Ra = Round(sqrt(pow(xx - x, 2.0) + pow(yy - y, 2.0)));
                CircleDirect(hdc, x, y, Ra, color);
            }
           else if(shape=="s5")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                Ra = Round(sqrt(pow(xx - x, 2.0) + pow(yy - y, 2.0)));
                CirclePolar(hdc, x, y, Ra, color);
            }
           else if(shape=="s6")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                Ra = Round(sqrt(pow(xx - x, 2.0) + pow(yy - y, 2.0)));
                CircleIterativePolar(hdc, x, y, Ra, color);
            }
           else if(shape=="s7")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                Ra = Round(sqrt(pow(xx - x, 2.0) + pow(yy - y, 2.0)));
                circle_mid(hdc, x, y, Ra, color);

            }
            else if(shape=="s8")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                Ra = Round(sqrt(pow(xx - x, 2.0) + pow(yy - y, 2.0)));
                CircleModifiedMid(hdc, x, y, Ra, color);
            }
           else if(shape=="s9")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                FillingQuadrantByLines(hdc, x, y, xx, yy, color);
            }
           else if(shape=="s10")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                 FillingQuadrantByLines(hdc, x, y, xx, yy, color);            }
           else if(shape=="s11")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                FillingQuadrantByLines(hdc, x, y, xx, yy, color);
            }
           else if(shape=="s12")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                 FillingQuadrantByLines(hdc, x, y, xx, yy, color);
            }
           else if(shape=="s13")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                 FillingQuadrantByCircles(hdc, x, y, xx, yy, color);
            }
           else if(shape=="s14")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                 FillingQuadrantByCircles(hdc, x, y, xx, yy, color);
            }
           else  if(shape=="s15")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                FillingQuadrantByCircles(hdc, x, y, xx, yy, color);
            }
            else if(shape=="s16")
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                FillingQuadrantByCircles(hdc, x, y, xx, yy, color);
            }
           else if(shape=="s17")
            {
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                inFile >> x1;
                inFile >> y1;
                inFile >> x2;
                inFile >> y2;
                Rectangle(hdc,l,t,r,b);
                FillingSquareWithHermite(hdc,l,r,t,b,y,yy,y1,y2,x,xx,x1,x2,color);

            }
             else if(shape=="s18")
            {
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                inFile >> x1;
                inFile >> y1;
                inFile >> x2;
                inFile >> y2;
                Rectangle(hdc,l,t,r,b);
               FillingSquareWithHermite(hdc,l,r,t,b,y,yy,y1,y2,x,xx,x1,x2,color);

            }
             else if(shape=="s19")
            {
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                inFile >> P[0].x;
                inFile >> P[0].y;
                inFile >> P[1].x;
                inFile >> P[1].y;
                inFile >> P[2].x;
                inFile >> P[2].y;
                inFile >> P[3].x;
                inFile >> P[3].y;
                inFile >> P[4].x;
                inFile >> P[4].y;
                Rectangle(hdc,l,t,r,b);
                convexFilling(hdc, P, 5, color);


            }
             else if(shape=="s20")
            {
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                inFile >> P[0].x;
                inFile >> P[0].y;
                inFile >> P[1].x;
                inFile >> P[1].y;
                inFile >> P[2].x;
                inFile >> P[2].y;
                inFile >> P[3].x;
                inFile >> P[3].y;
                inFile >> P[4].x;
                inFile >> P[4].y;
                Rectangle(hdc,l,t,r,b);
                NonConvexFill(hdc, P, 5, color);

            }
             else if(shape=="s21")
            {
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                inFile >> x;
                inFile >> y;
                Rectangle(hdc,l,t,r,b);
               RECFloodFill(hdc, x, y,RGB(0, 0, 0), color);
            }
             else if(shape=="s22")
            {
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                inFile >> x;
                inFile >> y;
                Rectangle(hdc,l,t,r,b);
                NonRecursiveFloodFill(hdc, x, y,RGB(0, 0, 0), color);


            }
             else if(shape=="s23")
            {
                 inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                inFile >> x1;
                inFile >> y1;
                inFile >> x2;
                inFile >> y2;
                cardinalSpline(hdc,cc,x,y,xx,yy,x1,y1,x2,y2,color);

            }
             else if(shape=="s24")
            {
                 inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                EllipsePolar(hdc, x, y, xx, yy, color);

            }
             else if(shape=="s25")
            {
                 inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                Direct_Ellipse(hdc, x, y, xx, yy, color);

            }
             else if(shape=="s26")
            {
                 inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                 Mid_Ellipse(hdc, x, y, xx, yy, color);

            }
             else if(shape=="s27" )
            {
                inFile >> x;
                inFile >> y;
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                Rectangle(hdc,l,r,t,b);

                PointClipping(hdc,x,y,l,r,t,b,color);

            }
             else if(shape=="s30")
            {
                 inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;
                Rectangle(hdc,l,t,r,b);
                LINE_Clipping(hdc,x,y,xx,yy,l,r,t,b);

            }
             else if(shape=="s31" )
            {
                inFile >> x;
                inFile >> y;
                inFile >> xx;
                inFile >> yy;
                inFile >> x1;
                inFile >> y1;
                inFile >> x2;
                inFile >> y2;
                inFile >> x3;
                inFile >> y3;
                inFile >> l;
                inFile >> r;
                inFile >> t;
                inFile >> b;

                X.push_back(x);
                Y.push_back(y);
                X.push_back(xx);
                Y.push_back(yy);
                X.push_back(x1);
                Y.push_back(y1);
                X.push_back(x2);
                Y.push_back(y2);
                X.push_back(x3);
                Y.push_back(y3);
                Rectangle(hdc,l,t,r,b);
                Polygon_Clip( hdc, X, Y, l, r,  t, b);
                X.clear();
                Y.clear();

            }
             else if(shape=="s41" )
            {
                 inFile >> xc1;
                inFile >> yc1;
                inFile >> xc2;
                inFile >> yc2;
                inFile >> R1;
                inFile >> R2;


               DrawCircleBresenham(hdc, xc1, yc1, R1, color);
               DrawCircleBresenham(hdc, xc2, yc2, R2, color);
                   FillingIntersection(hdc , R1,R2, xc1, yc1, xc2, yc2 , color);

            }

        }

        inFile.close(); // Close the file
    } else {
        std::cout << "Failed to open the file." << std::endl;
    }


}


LRESULT CALLBACK WindowProcedure (HWND, UINT, WPARAM, LPARAM);

/* Make the class name into a global variable */
TCHAR szClassName[ ] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain (HINSTANCE hThisInstance,
HINSTANCE hPrevInstance,

LPSTR lpszArgument,
int nCmdShow)
{
HWND hwnd; /* This is the handle for our window */
MSG messages; /* Here messages to the application are saved */
WNDCLASSEX wincl; /* Data structure for the windowclass */

/* The Window structure */
wincl.hInstance = hThisInstance;
wincl.lpszClassName = szClassName;
wincl.lpfnWndProc = WindowProcedure; /* This function is called by windows */
wincl.style = CS_DBLCLKS; /* Catch double-clicks */
wincl.cbSize = sizeof (WNDCLASSEX);

/* Use default icon and mouse-pointer */
wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
wincl.hCursor = LoadCursor (NULL, IDC_HAND);
//IDC_WAIT loading circle
//IDC_IBEAM text editing
//IDC_CROSS plus sign
//IDC_UPARROW upp arrow
//IDC_SIZENWSE diagonal arrow
//IDC_SIZEWE horizontal arrow
//IDC_HAND handd
//IDC_NO no arrow
//IDC_SIZEALL all directions arrow

wincl.lpszMenuName = NULL; /* No menu */
wincl.cbClsExtra = 0; /* No extra bytes after the window class */
wincl.cbWndExtra = 0; /* structure or the window instance */
/* Use Windows's default colour as the background of the window */
wincl.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
/* Register the window class, and if it fails quit the program */
if (!RegisterClassEx (&wincl))
return 0;

/* The class is registered, let's create the program*/
hwnd = CreateWindowEx (
0, /* Extended possibilites for variation */
szClassName, /* Classname */
_T("Code::Blocks Template Windows App"), /* Title Text */
WS_OVERLAPPEDWINDOW, /* default window */
CW_USEDEFAULT, /* Windows decides the position */
CW_USEDEFAULT, /* where the window ends up on the screen */
544, /* The programs width */
375, /* and height in pixels */
HWND_DESKTOP, /* The window is a child-window to desktop */
NULL, /* No menu */
hThisInstance, /* Program Instance handler */
NULL /* No Window Creation data */
);

/* Make the window visible on the screen */
ShowWindow (hwnd, nCmdShow);

/* Run the message loop. It will run until GetMessage() returns 0 */
while (GetMessage (&messages, NULL, 0, 0))
{
/* Translate virtual-key messages into character messages */
TranslateMessage(&messages);
/* Send message to WindowProcedure */
DispatchMessage(&messages);
}

/* The program return-value is 0 - The value that PostQuitMessage() gave */
return messages.wParam;
}

COLORREF color = RGB(0, 0, 0);
int c = 0;
static POINT P[5];
int case_num = 0, Quart = 0;
int xl, yt, xr, yb,xc1,xc2,yc1,yc2;


LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT mcode, WPARAM wp, LPARAM lp)
{
	HDC hdc;
hdc = GetDC(hwnd);
HMENU hMenu;

	switch (mcode)
	{
	     case WM_CREATE:{

        HICON hIcon;
        HMENU hMenu=CreateMenu();
        HMENU Line=CreateMenu();
        HMENU Circle=CreateMenu();
        HMENU Filling_Circle_Line=CreateMenu();
        HMENU Filling_Circle_Circle=CreateMenu();
        HMENU Square_Clipping=CreateMenu();
        HMENU Rectangle_Clipping=CreateMenu();
        HMENU color_list = CreateMenu();



        HMENU Ellipse=CreateMenu();
        AppendMenu(hMenu,MF_POPUP,(UINT_PTR)Line,"Line");
        AppendMenu(Line,MF_STRING,1,"DDA");
        AppendMenu(Line,MF_STRING,2,"Parametric");
        AppendMenu(Line,MF_STRING,3,"MidPoint");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_POPUP,(UINT_PTR)Circle,"Circle");
        AppendMenu(Circle,MF_STRING,4,"Direct");
        AppendMenu(Circle,MF_STRING,5,"Polar");
        AppendMenu(Circle,MF_STRING,6,"Iterative Polar");
        AppendMenu(Circle,MF_STRING,7,"MidPoint");
        AppendMenu(Circle,MF_STRING,8,"Modified MidPoint");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_POPUP,(UINT_PTR)Filling_Circle_Line,"Filling Circle with lines");
        AppendMenu(Filling_Circle_Line,MF_STRING,9, "1st Quarter");
        AppendMenu(Filling_Circle_Line,MF_STRING,10,"2nd Quarter");
        AppendMenu(Filling_Circle_Line,MF_STRING,11,"3rd Quarter");
        AppendMenu(Filling_Circle_Line,MF_STRING,12,"4th Quarter");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_POPUP,(UINT_PTR)Filling_Circle_Circle,"Filling Circle with circles");
        AppendMenu(Filling_Circle_Circle,MF_STRING,13,"1st Quarter");
        AppendMenu(Filling_Circle_Circle,MF_STRING,14,"2nd Quarter");
        AppendMenu(Filling_Circle_Circle,MF_STRING,15,"3rd Quarter");
        AppendMenu(Filling_Circle_Circle,MF_STRING,16,"4th Quarter");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);
////
        AppendMenu(hMenu,MF_STRING,17,"Filling Square with Hermite");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);
        AppendMenu(hMenu,MF_STRING,18,"Filling Rectangle with Bezier");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);
         AppendMenu(hMenu,MF_STRING,41,"Fill 2 circles intersection");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_STRING,19,"Convex Filling");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);
        AppendMenu(hMenu,MF_STRING,20,"Non-Convex Filling");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);
///
        AppendMenu(hMenu,MF_STRING,21,"Recursive Floodfill");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);
        AppendMenu(hMenu,MF_STRING,22,"Non-Recursive Floodfill");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_STRING,23,"Cardinal Spline Curve");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_POPUP,(UINT_PTR)Ellipse,"Ellipse");
        AppendMenu(Ellipse,MF_STRING,24,"Polar");
        AppendMenu(Ellipse,MF_STRING,25,"Direct");
        AppendMenu(Ellipse,MF_STRING,26,"MidPoint");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_POPUP,(UINT_PTR)Square_Clipping,"Square Clipping");
        AppendMenu(Square_Clipping,MF_STRING,27,"Point");
        AppendMenu(Square_Clipping,MF_STRING,28,"Line");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_POPUP,(UINT_PTR)Rectangle_Clipping,"Rectangle Clipping");
        AppendMenu(Rectangle_Clipping,MF_STRING,29,"Point");
        AppendMenu(Rectangle_Clipping,MF_STRING,30,"Line");
        AppendMenu(Rectangle_Clipping,MF_STRING,31,"Polygon");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);
///
        AppendMenu(hMenu,MF_STRING,32,"Save");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_STRING,33,"Load");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

        AppendMenu(hMenu,MF_STRING,34,"Clean screen");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

         AppendMenu(hMenu, MF_POPUP, (UINT_PTR)color_list, "Color");
         AppendMenu(color_list, MF_STRING, 35, "red");
         AppendMenu(color_list, MF_STRING, 36, "Orange");
         AppendMenu(color_list, MF_STRING, 37, "Yellow");
         AppendMenu(color_list, MF_STRING, 38, "Blue");
         AppendMenu(color_list, MF_STRING, 39, "Green");
         AppendMenu(color_list, MF_STRING, 40, "Black");
        AppendMenu(hMenu, MF_SEPARATOR, 0, NULL);

         SetMenu(hwnd, hMenu);
 break;

}
  case WM_COMMAND:
        switch (wp)
        {

            /////////////////////////////
        case  (1):
            case_num = 1;
            cout << "DDA Line Algorithm...." << endl;
            break;
        case  (2):
            case_num = 2;
            cout << "Parametric Line Algorithm...." << endl;
            break;

        case  (3):
            case_num = 3;
            cout << "MidPoint Line Algorithm...." << endl;

            break;
            ///---------------------------------------------------
        case  (4):
            case_num = 4;
            cout << "Direct Circle Algorithm....." << endl;

            break;
        case  (5):
            case_num = 5;
            cout << "Polar circle Algorithm....." << endl;
            break;

        case  (6):
            case_num = 6;
            cout << "Iterative Polar Circle Algorithm....." << endl;
            break;
        case  (7):
            case_num = 7;
            cout << "MidPoint Circle Algorithm....." << endl;
            break;
        case  (8):
            case_num = 8;
            cout << "MidPoint Modification Circle Algorithm....." << endl;
            break;
            ///-----------------------------------------------------------------
        case  (9):
            case_num = 9;
             Quart = 1;
             save("s9",color);
            cout << "Filling Circle quarter 1 with lines Algorithm....." << endl;
            break;
        case  (10):
            case_num = 10;
             Quart = 2;
             save("s10",color);
            cout << "Filling Circle quarter 2 with lines Algorithm....." << endl;
            break;
        case  (11):
            case_num = 11;
             Quart = 3;
             save("s11",color);
            cout << "Filling Circle quarter 3 with lines Algorithm....." << endl;
            break;

        case  (12):
            case_num = 12;
            Quart = 4;
            save("s12",color);
            cout << "Filling Circle quarter 4 with lines Algorithm....." << endl;
            break;

        case  (13):
            case_num = 13;
            Quart = 1;
            save("s13",color);
            cout << "Filling Circle quarter 1 with circles Algorithm....." << endl;
            break;
        case  (14):
            case_num = 14;
            Quart = 2;
            save("s14",color);
            cout << "Filling Circle quarter 2 with circles Algorithm....." << endl;
            break;

        case  (15):
            case_num = 15;
            Quart = 3;
            save("s15",color);
            cout << "Filling Circle quarter 3 with circles Algorithm....."<< endl;
            break;

        case  (16):
            case_num = 16;
            Quart = 4;
            save("s16",color);
            cout << "Filling Circle quarter 4 with circles Algorithm.....";
            break;
            ///----------------------------------------------------------------
        case  (17):
            case_num = 17;
            cout << "Filling Square with Hermite curve \n";
            break;
        case  (18):
            case_num = 18;
            cout << "Filling Rectangle with Bezier curve \n";
            break;
        case  (19):
            case_num = 19;
            cout << "Convex Filling \n";
            break;
        case  (20):
            case_num = 20;
            cout << "Non-Convex Filling \n";
            break;
        case  (21):
            case_num = 21;
            cout << "Recursive Floodfill \n";
            break;
        case  (22):
            case_num = 22;
             cout << "Non-Recursive Floodfill \n";
            break;
        case  (23):
            case_num = 23;
             cout << "Cardinal Spline Curve \n";
            break;
            ///-------------------------------------------
        case  (24):
            case_num = 24;
             cout << "Ellipse Polar Algorithm \n";
            break;
        case  (25):
            case_num = 25;
             cout << "Ellipse Direct Algorithm \n";
            break;
        case  (26):
            case_num = 26;
             cout << "Ellipse MidPoint Algorithm \n";
            break;
             ///-------------------------------------------
        case  (27):
            case_num = 27;
             cout << "Square Clipping Point Algorithm \n";
            break;
        case  (28):
            case_num = 28;
             cout << "Square Clipping Line Algorithm \n";
            break;
        case  (29):
            case_num = 29;
             cout << "Rectangle Clipping Point Algorithm \n";
            break;
        case  (30):
            case_num = 30;
             cout << "Rectangle Clipping Line Algorithm \n";
            break;
        case (31):
            case_num = 31;
            save("s31",color);
             cout << "Rectangle Clipping Polygon Algorithm \n";
         ///-------------------------------------------
            break;
        case (32):
            RECT rect;
            if (GetWindowRect(hwnd, &rect)) {

                rect.top += 8;
                rect.left += 8;
               // HDCToFile("picture.bmp", hdc, rect, 24);
                ReleaseDC(hwnd, hdc);
            }
            break;
        case(33):
           Load(hwnd, hdc);
            break;
        case(34):
            InvalidateRect(hwnd, NULL, TRUE);
            cout << "Window is clear now ...." << endl;
            break;
        case  (35):
            color = RGB(255, 0, 0);
            cout << "Red color " << endl;
            break;

        case  (36):
            color = RGB(255, 165, 0);
            cout << "Orange color" << endl;

            break;
        case  (37):
            color = RGB(255, 255, 0);
            cout << "Yellow color" << endl;

            break;
        case  (38):
             color = RGB(0, 0, 255);
            cout << "Blue color" << endl;
            break;
        case  (39):
            color = RGB(0, 110, 0);
            cout << "Green color" << endl;
            break;

        case  (40):
             color = RGB(0, 0, 0);
            cout << "Black color" << endl;

            break;
            case  (41):
            case_num = 41;

            cout << "Filling Intersection of 2 Circles Algorithm....."<< endl;
            break;

        }
        break;

    case WM_LBUTTONDOWN:
        if (case_num >= 1 && case_num <= 16 || (case_num >= 24 && case_num <= 26)) //  (Direct, Polar, iterative Polar, midpoint and modified Midpoint)
        {
            x11 = LOWORD(lp);
            y11 = HIWORD(lp);
        }


        else if ((case_num >= 27 && case_num <= 31)|| case_num == 21 || case_num ==22 || (case_num >= 17 && case_num <= 18)) {
            if (c == 0)
            {
                xl = LOWORD(lp);
                yt = HIWORD(lp);
                c++;
            }
            else
            {
                xr = LOWORD(lp);
                yb = HIWORD(lp);
                if (case_num == 27 || case_num == 28 || case_num == 17) {
                        xr=xl+110;
                        yb=yt+110;
                    Rectangle(hdc, xl, yt,xr,yb );
                }
                else
                    Rectangle(hdc, xl, yt, xr, yb);
                c = 0;
            }
        }
          else if ((case_num >= 19 && case_num <= 20))//convixfilling
          {

            P[c].x = LOWORD(lp);
            P[c].y = HIWORD(lp);

            c++;
            if (c ==5)
            {
                if(case_num == 19 )
                {
                    save("s19",color);
                    savePoint(xl,xr);
                    savePoint(yt,yb);
                     savePoint(P[0].x,P[0].y);
                     savePoint(P[1].x,P[1].y);
                     savePoint(P[2].x,P[2].y);
                     savePoint(P[3].x,P[3].y);
                     savePoint(P[4].x,P[4].y);

                    convexFilling(hdc, P, 5, color);

                }
                else
                {
                    save("s20",color);
                    savePoint(xl,xr);
                    savePoint(yt,yb);
                     savePoint(P[0].x,P[0].y);
                     savePoint(P[1].x,P[1].y);
                     savePoint(P[2].x,P[2].y);
                     savePoint(P[3].x,P[3].y);
                     savePoint(P[4].x,P[4].y);
                    NonConvexFill(hdc, P, 5, color);
                }
                c = 0;
            }




          }
          else if(case_num==41)//extra task
          {
              if(c==0)
              {
                    xc1 = LOWORD(lp);
                    yc1 = HIWORD(lp);
                    c++;
                    save("s41",color);
                     savePoint(xc1,yc1);

               }
               else if(c==1)
                {
                    x11 = LOWORD(lp);
                    y11 = HIWORD(lp);
                    R1 = sqrt(pow((x11 - xc1), 2) + pow((y11 - yc1), 2));
                    DrawCircleBresenham(hdc, xc1, yc1, R1, color);

                    c++;
               }
                 else if(c==2)
                {
                    xc2 = LOWORD(lp);
                    yc2 = HIWORD(lp);
                    c++;
                    savePoint(xc2,yc2);
                }
                else if(c==3)
                {
                    x22 = LOWORD(lp);
                    y22 = HIWORD(lp);
                    c = 0;

                    R2 = sqrt(pow((x22 - xc2), 2) + pow((y22 - yc2), 2));
                    savePoint(R1,R2);
                    DrawCircleBresenham(hdc, xc2, yc2, R2, color);
                   FillingIntersection(hdc , R1,R2, xc1, yc1, xc2, yc2 , color);
                }
          }
        ///


        break;
case WM_RBUTTONDOWN:
        if ((case_num >= 1 && case_num <= 16) || (case_num == 21 || case_num ==22) || (case_num >= 24 && case_num <= 26)) //Circle(Direct,Polar,Midpoint)
        {
            x22 = LOWORD(lp);
            y22 = HIWORD(lp);
            Ra = Round(sqrt(pow(x22 - x11, 2.0) + pow(y22 - y11, 2.0)));
            A=sqrt(pow((x22-x11),2)+pow((y22-y11),2));
            B=A+20;
            if (case_num == 1) //line(DDA)
            {
                LineDDA(hdc, x11, y11, x22, y22, color);
                save("s1",color);
                savePoint(x11,y11);
                savePoint(x22,y22);
            }
            else if (case_num == 2)//line(Direct)
            {
                DirectLine(hdc, x11, y11, x22, y22, color);
                save("s2",color);
                savePoint(x11,y11);
                savePoint(x22,y22);
            }
            else if (case_num == 3)//line(mid)
            {
                MidPointLine(hdc, x11, y11, x22, y22, color);
                save("s3",color);
                savePoint(x11,y11);
                savePoint(x22,y22);
            }
            else if (case_num == 4)//Circle(Direct)
            {
                CircleDirect(hdc, x11, y11, Ra, color);
                save("s4",color);
                savePoint(x11,y11);
                savePoint(x22,y22);
            }
            else if (case_num == 5)//Circle(polar)
            {
                CirclePolar(hdc, x11, y11, Ra, color);
                save("s5",color);
                savePoint(x11,y11);
                savePoint(x22,y22);

            }
            else if (case_num == 6)//Circle iterative polar
            {
                CircleIterativePolar(hdc, x11, y11, Ra, color);
                save("s6",color);
                savePoint(x11,y11);
                savePoint(x22,y22);
            }
            else if (case_num == 7) //Circle midpoint
            {
                circle_mid(hdc, x11, y11, Ra, color);
                save("s7",color);
                savePoint(x11,y11);
                savePoint(x22,y22);
            }
             else if (case_num == 8) //Circle modified midpoint
            {
                CircleModifiedMid(hdc, x11, y11, Ra, color);
                save("s8",color);
                savePoint(x11,y11);
                savePoint(x22,y22);
            }


            else if ((case_num >= 9 && case_num <= 12)) //Filling Quadrant By Lines
            {
                FillingQuadrantByLines(hdc, x11, y11, Ra, Quart, color);
                savePoint(x11,y11);
                savePoint(Ra,Quart);
            }
            else if ((case_num >= 13 && case_num <= 16)) //Filling Quadrant By Circles
            {
                FillingQuadrantByCircles(hdc, x11, y11, Ra, Quart, color);
                savePoint(x11,y11);
                savePoint(Ra,Quart);
            }
            else if (case_num ==21 ) //Flood fill
            {
                hdc = GetDC(hwnd);
                RECFloodFill(hdc, x22, y22,RGB(0, 0, 0), color);
                save("s21",color);
                savePoint(xl,xr);
                savePoint(yt,yb);
                savePoint(x22,y22);

                ReleaseDC(hwnd, hdc);
            }
            else if (case_num ==22 ) //Flood fill
            {
                hdc = GetDC(hwnd);
                NonRecursiveFloodFill(hdc, x22, y22,RGB(0, 0, 0), color);
                save("s22",color);
                savePoint(xl,xr);
                savePoint(yt,yb);
                savePoint(x22,y22);
                 ReleaseDC(hwnd, hdc);
            }

             else if (case_num == 24)//Ellipse  polar
            {
                EllipsePolar(hdc, A, B, x11, y11, color);
                save("s24",color);
                savePoint(A,B);
                savePoint(x11,y11);
            }
            else if (case_num == 25) //Ellipse Direct
            {
                Direct_Ellipse(hdc, A, B, x11, y11, color);
                save("s25",color);
                savePoint(A,B);
                savePoint(x11,y11);
            }
             else if (case_num == 26) //Ellipse midpoint
            {
                Mid_Ellipse(hdc, A, B, x11, y11, color);
                save("s26",color);
                savePoint(A,B);
                savePoint(x11,y11);
            }
            break;
        }
        /////


        else if (case_num == 23 || (case_num >= 17 && case_num <= 18)) { //cardinal curve
            if(c==0){
		z1 = LOWORD(lp);
		k1 = HIWORD(lp);c++;

	    }
	    else if(c==1)
        {
        z2 = LOWORD(lp);
		k2 = HIWORD(lp);
		c++;

        }
        else if( c == 2)
        {

		z3 = LOWORD(lp);
		k3 = HIWORD(lp);
		c++;

		}
		else if( c == 3){


		z4 = LOWORD(lp);
		k4 = HIWORD(lp);
		c=0;
		if(case_num==23){

		cardinalSpline(hdc,cc,z1,k1,z2,k2,z3,k3,z4,k4,color);
		        save("s23",color);
                savePoint(z1,k1);
                savePoint(z2,k2);
                savePoint(z3,k3);
                savePoint(z4,k4);
		}
		else if (case_num==17){
              FillingSquareWithHermite(hdc,xl,xr,yt,yb,k1,k2,k3,k4,z1,z2,z3,z4,color);
              save("s17",color);
                savePoint(xl,xr);
                savePoint(yt,yb);
                savePoint(z1,k1);
                savePoint(z2,k2);
                savePoint(z3,k3);
                savePoint(z4,k4);

		}
		else if (case_num==18){

              FillingRectangleWithBezier(hdc,xl,xr,yt,yb,k1,k2,k3,k4,z1,z2,z3,z4,color);
              save("s18",color);
              savePoint(xl,xr);
                savePoint(yt,yb);
                savePoint(z1,k1);
                savePoint(z2,k2);
                savePoint(z3,k3);
                savePoint(z4,k4);

		}
		}
        }
        else if (case_num == 31)
        {


            if (c < 5) {
                X.push_back ( LOWORD(lp));
                Y.push_back ( HIWORD(lp));
                savePoint(X[c],Y[c]);
                c++;
            }
            if (case_num == 31 && c == 5)//Polygon clipping
            {
                savePoint(xl,xr);
                savePoint(yt,yb);
                Polygon_Clip( hdc, X, Y, xl, xr,  yt, yb);
                X.clear();
                Y.clear();
               c = 0;
            }

        }
        else if (case_num == 30 || case_num == 28) {//line cliping

            if (c == 0)
            {
               xs =LOWORD(lp);
               ys = HIWORD(lp);
                c++;
            }
            else if (c == 1)
            {
                xe = LOWORD(lp);
                ye = HIWORD(lp);

                LINE_Clipping(hdc,xs,ys,xe,ye,xl,xr,yt,yb);
                save("s30",color);
                savePoint(xs,ys);
                savePoint(xe,ye);
                savePoint(xl,xr);
                savePoint(yt,yb);
                c = 0;
            }
        }
        else if (case_num == 27 || case_num == 29)
            {//point clipping

            x22 = LOWORD(lp);
            y22 = HIWORD(lp);
           PointClipping(hdc,x22,y22,xl,yt,xr,yb,color);
                save("s27",color);
                savePoint(x22,y22);
                savePoint(xl,yt);
                savePoint(xr,yb);
            }

        break;
	case WM_CLOSE:
		DestroyWindow(hwnd);
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default: return DefWindowProc(hwnd, mcode, wp, lp);
	}
	return 0;
}

