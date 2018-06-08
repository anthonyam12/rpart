/*
 * The routine which will find the best split for a node
 *
 * Input :      node
 *              node number
 *
 * Output:      Fills in the node's
 *                      primary splits
 *                      competitor splits
 */
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"

void testPrint();
void doRpartLogic(pNode me, int n1, int n2);
void doDelayed(pNode me, int n1, int n2);
void checkXdataMatrix(double **xdata, double **ydata, int k);
double getSSSplit(double **xdata, double **ydata, double *wtdata, int num_obs);

void bsplit(pNode me, int n1, int n2)
{
    me->primary = (pSplit) NULL;
    if(!rp.delayed) 
    {
        doRpartLogic(me, n1, n2);
    }
    else 
    {
        doDelayed(me, n1, n2);
    }
}

void doRpartLogic(pNode me, int n1, int n2) 
{
    int i, j, k;
    int kk;
    int nc;
    double improve;
    double split = 0.0;
    pSplit tsplit;
    int *index;
    double *xtemp;              /* these 3 because I got tired of typeing
     * "rp.xtemp", etc */
    double **ytemp;
    double *wtemp;
    
    xtemp = rp.xtemp;
    ytemp = rp.ytemp;
    wtemp = rp.wtemp;
    
    for (i = 0; i < rp.nvar; i++) 
    {
        index = rp.sorts[i];
        nc = rp.numcat[i];
        /* extract x and y data */
        k = 0;
        for (j = n1; j < n2; j++) 
        {
            kk = index[j];
            /* x data not missing and wt > 0 */
            if (kk >= 0 && rp.wt[kk] > 0) 
            {  
                xtemp[k] = rp.xdata[i][kk];
                ytemp[k] = rp.ydata[kk];
                wtemp[k] = rp.wt[kk];
                k++;
            }
        }
        //if ((n2-n1) == 209 && i == 0)
        //    printf("nc: %d, xtemp[0] = %5g, xtemp[k-1] = %5g, k=%d\n", nc, xtemp[0], xtemp[k-1], k);
        if (k == 0 || (nc == 0 && xtemp[0] == xtemp[k - 1]))
        {
            continue;           /* no place to split */
        }
            
            // rpart split function - find best way to partition the dataset
            (*rp_choose) (k, ytemp, xtemp, nc, rp.min_node, &improve, &split, rp.csplit, me->risk, wtemp);
            
            /*
            * Originally, this just said "if (improve > 0)", but rounding
            * error will sometimes create a non zero that should be 0.  Yet we
            * want to retain invariance to the scale of "improve".
            */
        if (improve > rp.iscale) {
            rp.iscale = improve;  /* largest seen so far */
        }
        
        if (improve > (rp.iscale * 1e-10)) 
        {
            improve /= rp.vcost[i];     /* scale the improvement */
            tsplit = insert_split(&(me->primary), nc, improve, rp.maxpri);
            if (tsplit) 
            {
                tsplit->improve = improve;
                tsplit->var_num = i;
                tsplit->spoint = split;
                tsplit->count = k;
                if (nc == 0) 
                {
                    tsplit->spoint = split;
                    tsplit->csplit[0] = rp.csplit[0];
                } 
                else 
                {
                    for (k = 0; k < nc; k++)
                        tsplit->csplit[k] = rp.csplit[k];
                }
            }
        }
    }
}

void doDelayed(pNode me, int n1, int n2)
{
    int kk;
    int nc;
    double improve;
    double split = 0.0;
    int *index;
    double *xtemp;              /* these 3 because I got tired of typeing "rp.xtemp", etc */
    double **ytemp;
    double *wtemp;
    
    xtemp = rp.xtemp;
    ytemp = rp.ytemp;
    wtemp = rp.wtemp;
    
    double bestSS = DBL_MAX;
    double nodeSS = me->risk;   // check that this is right, might just call eval on 'me'
    int leftSS, rightSS;
    double **xdata;
    xdata = (double **) ALLOC(rp.nvar, sizeof(double*));
    
    for(int i = 0; i < rp.nvar; i++) 
    {
        index = rp.sorts[i];
        nc = rp.numcat[i];
        /* extract x and y data */
        int k = 0;
        for (int j = n1; j < n2; j++) 
        {
            kk = index[j];
            /* x data not missing and wt > 0 */
            if (kk >= 0 && rp.wt[kk] > 0) 
                k++;
        }
        // no where to split
        if (k == 0 || (nc == 0 && xtemp[0] == xtemp[k - 1]))
            continue;
        
        // build matrix of all column data (xdata) rather than just this column's data (xtemp)
        for(int j = 0; j < rp.nvar; j++) 
            xdata[j] = (double*) ALLOC(k, sizeof(double));
        
        k = 0;
        for(int j = n1; j < n2; j++)
        {
            kk = index[j];
            if(kk >= 0 && rp.wt[kk] > 0)
            {                
                xtemp[k] = rp.xdata[i][kk];
                ytemp[k] = rp.ydata[kk];
                wtemp[k] = rp.wt[kk];
                for(int l = 0; l < rp.nvar; l++)
                    xdata[l][k] = rp.xdata[l][kk];
                k++;
            }
        }
                
        // get y value for a row via *ytemp[k] thus xtemp[i][*] correspondes to *ytemp[i]
        (*rp_choose) (k, ytemp, xtemp, nc, rp.min_node, &improve, &split, rp.csplit, me->risk, wtemp);
        
        int splitsize = sizeof(Split) + ((nc == 0 ? 1 : nc) - 20) * sizeof(int); 
        pSplit bestSplit = (pSplit) CALLOC(1, splitsize);
        bestSplit->improve = improve;
        bestSplit->var_num = i;
        bestSplit->spoint = split;
        bestSplit->count = k;
        if (nc == 0) 
        {
            bestSplit->spoint = split;
            bestSplit->csplit[0] = rp.csplit[0];
        } 
        else // categorical
        {
            for (int j = 0; j < nc; j++)
                bestSplit->csplit[j] = rp.csplit[j];
        }
        
        // From Documentation: csplit[0]: 1: <x to the left, -1: <x to the right
        // TODO: handle categorical splits**************************************************************
        int rightCount = 0, leftCount = 0;
        for(int j = 0; j < k; j++) 
        {
            if(bestSplit->csplit[0] > 0) // <x go left
            {
                if(xtemp[j] < split) 
                    leftCount++;
                else
                    rightCount++;
           }
            else  // <x go right (>x go left)
            {
                if(xtemp[j] > split) 
                    leftCount++;
                else
                    rightCount++;
            }
        }

        // build left and right y data vectors
        double **rightY, **leftY, **leftX, **rightX, *leftWt, *rightWt;
        rightY = (double**) ALLOC(rightCount, sizeof(double*));
        leftY = (double**) ALLOC(leftCount, sizeof(double*));
        leftWt = (double*) ALLOC(leftCount, sizeof(double));
        rightWt = (double*) ALLOC(rightCount, sizeof(double));
        leftX = (double **) ALLOC(rp.nvar, sizeof(double*));
        rightX = (double **) ALLOC(rp.nvar, sizeof(double*));
        for(int j = 0; j < rp.nvar; j++)
        {
            leftX[j] = (double*) ALLOC(leftCount, sizeof(double));
            rightX[j] = (double*) ALLOC(rightCount, sizeof(double));
        }

        rightCount = 0;
        leftCount = 0;
        for(int j = 0; j < k; j++) 
        {
            if(bestSplit->csplit[0] > 0) // <x go left
            {
                if(xtemp[j] < split)
                { 
                    leftY[leftCount] = ytemp[j];
                    leftWt[leftCount] = wtemp[j];
                    for(int l = 0; l < rp.nvar; l++)
                        leftX[l][leftCount] = xdata[l][j];
                    leftCount++;
                }
                else
                {
                    rightY[rightCount] = ytemp[j];
                    rightWt[rightCount] = wtemp[j];
                    for(int l = 0; l < rp.nvar; l++)
                        rightX[l][rightCount] = xdata[l][j];
                    rightCount++;
                }
            }
            else  // <x go right (>x go left)
            {
                if(xtemp[j] > split)
                {
                    leftY[leftCount] = ytemp[j];
                    leftWt[leftCount] = wtemp[j];
                    for(int l = 0; l < rp.nvar; l++)
                        leftX[l][leftCount] = xdata[l][j];
                    leftCount++;
                }
                else
                {
                    rightY[rightCount] = ytemp[j];
                    rightWt[rightCount] = wtemp[j];
                    for(int l = 0; l < rp.nvar; l++)
                        rightX[l][rightCount] = xdata[l][j];
                    rightCount++;
                }
            }
        }
        
        // get SS for right
        double rightSS = getSSSplit(rightX, rightY, rightWt, rightCount);
        // get SS for left
        double leftSS = getSSSplit(leftX, leftY, leftWt, leftCount);
    
        // bestSS = min{nodess - (leftSS + rightSS)}, this maximize right and left SS
        double thisSS = nodeSS - (leftSS + rightSS);
        
        if(thisSS < bestSS)
        {
            if(improve > rp.iscale)
                rp.iscale = improve;
            
            if (improve > (rp.iscale * 1e-10)) 
            {
                me->primary = bestSplit;
                bestSS = thisSS;
                //if((n2-n1) == 1599)
                //    printf("%5g, %d\n", me->primary->spoint, me->primary->var_num);
            }
        }
    }
}

// ensure matrices are build correctly
void checkXdataMatrix(double **xdata, double **ydata, int k)
{
    //for(int j = 0; j < k; j++) // ydata[j] == xdata[j][*]?
       // printf("y=%5g -- x = %5g, %5g, %5g, %5g, %5g, %5g, %5g, %5g, %5g, %5g, %5g\n",
        //       *ydata[j], xdata[0][j], xdata[1][j], xdata[2][j], xdata[3][j], xdata[4][j], 
        //        xdata[5][k], xdata[6][j], xdata[7][j], xdata[8][j], xdata[9][j], xdata[10][j]);
}


// gets lowest SS from right/left split
// Will be done twice for each variable (left and right, L1 and L2)
double getSSSplit(double **xdata, double **ydata, double *wtdata, int num_obs)
{
    double improve;
    double split;
    double risk = 1;
    double *xtemp = (double*) ALLOC(num_obs, sizeof(double));
    double bestSS = 0;
    for (int i = 0; i < rp.nvar; i++) 
    {
        double leftSS = 0, rightSS = 0;
        int k = num_obs;
        int nc = rp.numcat[i];
        for(int j = 0; j < num_obs; j++) 
            xtemp[j] = xdata[i][j];
            
        // rpart split function - find best way to partition the dataset
        (*rp_choose) (num_obs, ydata, xtemp, nc, rp.min_node, &improve, &split, rp.csplit, risk, wtdata);
        
        // split data on split point
        // TODO: handle categorical splits**************************************************************
        int rightCount = 0, leftCount = 0;
        for(int j = 0; j < k; j++) 
        {
            if(rp.csplit[0] > 0) // <x go left
            {
                if(xtemp[j] < split) 
                    leftCount++;
                else
                    rightCount++;
            }
            else  // <x go right (>x go left)
            {
                if(xtemp[j] > split) 
                    leftCount++;
                else
                    rightCount++;
            }
        }
        // didn't find a split
        if(leftCount == 0 || rightCount == 0)
            continue;
        
        
        // build left and right y data vectors
        double **rightY, **leftY, *leftWt, *rightWt;
        rightY = (double**) ALLOC(rightCount, sizeof(double*));
        leftY = (double**) ALLOC(leftCount, sizeof(double*));
        leftWt = (double*) ALLOC(leftCount, sizeof(double));
        rightWt = (double*) ALLOC(rightCount, sizeof(double));
        
        rightCount = 0;
        leftCount = 0;
        for(int j = 0; j < k; j++) 
        {
            if(rp.csplit[0] > 0) // <x go left
            {
                if(xtemp[j] < split)
                { 
                    leftY[leftCount] = ydata[j];
                    leftWt[leftCount] = wtdata[j];
                    leftCount++;
                }
                else
                {
                    rightY[rightCount] = ydata[j];
                    rightWt[rightCount] = wtdata[j];
                    rightCount++;
                }
            }
            else  // <x go right (>x go left)
            {
                if(xtemp[j] > split)
                {
                    leftY[leftCount] = ydata[j];
                    leftWt[leftCount] = wtdata[j];
                    leftCount++;
                }
                else
                {
                    rightY[rightCount] = ydata[j];
                    rightWt[rightCount] = wtdata[j];
                    rightCount++;
                }
            }
        }
        
        // calc SS (rp_eval) (each side? Yes, and add them. (max{SS_R + SS_L}))
        double rightMean = 0, leftMean = 0;
        (*rp_eval) (rightCount, rightY, &rightMean, &rightSS, rightWt);
        (*rp_eval) (leftCount, leftY, &leftMean, &leftSS, leftWt);
        
        
        if((rightSS + leftSS) > bestSS) 
            bestSS = (rightSS + leftSS);
        
    }
    return bestSS;
}

