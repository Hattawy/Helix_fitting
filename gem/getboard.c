#include "rtpc.h"
int rowcol2board(const int row,const int col)
{
    if (row<0 || row>79 || col<0 || col>39) return -1;
    return 5*(row/2) + col/8;
}
int pad2board(const int pad)
{
    int row,col;
    pad2rowcol(pad,&row,&col);
    return rowcol2board(row,col);
}
