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

void bsplit(pNode me, int n1, int n2, int fromBSplit)
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
        
        if (k == 0 || (nc == 0 && xtemp[0] == xtemp[k - 1]))
            continue;           /* no place to split */
            
            // rpart split function - find best way to partition the dataset
            (*rp_choose) (k, ytemp, xtemp, nc, rp.min_node, &improve, &split, rp.csplit, me->risk, wtemp);
            
            /*
            * Originally, this just said "if (improve > 0)", but rounding
            * error will sometimes create a non zero that should be 0.  Yet we
            * want to retain invariance to the scale of "improve".
            */
        if (improve > rp.iscale)
            rp.iscale = improve;  /* largest seen so far */
        
        // What is this doing? Obviously inserting the split (maybe?), but when and why?
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
    
    // need to preserve sorts, which, and tempvec in rp
    int *s_which = (int*) ALLOC(rp.n, sizeof(int));
    int *s_tempvec = (int*) ALLOC(rp.n, sizeof(int));
    for(int r = 0; r < rp.n; r++) {
        s_which[r] = rp.which[r];
        s_tempvec[r] = rp.tempvec[r];
    }
    int *s_sorts = (int *) ALLOC(rp.n * rp.nvar, sizeof(int));
    memcpy(s_sorts, rp.sorts[0], rp.n * rp.nvar * sizeof(int));
    
    /* Helps to understand the sorts array functionality
    int k = 0;
    for (int j = 0; j < rp.n; j++) 
    {
        kk = rp.sorts[4][j];
        printf("%d == %5g\n", kk, rp.xdata[4][kk]);
    } 
    */
    
    for(int i = 0; i < rp.nvar; i++) 
    {
        int nleft = 0, nright = 0;
        index = rp.sorts[i];
        nc = rp.numcat[i];
        /* extract x and y data */
        int k = 0;
        for (int j = n1; j < n2; j++) 
        {
            kk = index[j];
            /* x data not missing and wt > 0 */
            if (kk >= 0 && rp.wt[kk] > 0) 
            {  
                // will need to build another xdata array 
                // for v in nvar: new_xd[v][kk] = xdata[v][kk]
                // below gives ith column, need to get ALL cols 
                // for second splits
                xtemp[k] = rp.xdata[i][kk];
                ytemp[k] = rp.ydata[kk];
                wtemp[k] = rp.wt[kk];
                //printf("%5g, %5g, %5g\n", xtemp[k], *ytemp[k], wtemp[k]);
                k++;
            }
        }
        
        // get y value for a row via *ytemp[k] thus xtemp[i][*] correspondes to *ytemp[i]
        
        (*rp_choose) (k, ytemp, xtemp, nc, rp.min_node, &improve, &split, rp.csplit, me->risk, wtemp);

        nc = nc == 0 ? 1 : nc;
        int splitsize = sizeof(Split) + (nc - 20) * sizeof(int);
        pSplit bestSplit = (pSplit) CALLOC(1, splitsize);
        
        bestSplit->improve = improve;
        bestSplit->var_num = i;
        bestSplit->spoint = split;
        bestSplit->count = k;
        if (nc == 0) 
        {
            // continuous split
            bestSplit->spoint = split;
            bestSplit->csplit[0] = rp.csplit[0];
        } 
        else 
        {
            // categorical split
            for (k = 0; k < nc; k++)
                bestSplit->csplit[k] = rp.csplit[k];
        }
        
        // preserve "me"
        pNode dummyNode = (pNode) ALLOC(1, nodesize);
        dummyNode->primary = bestSplit; 
        
        // calls below to get counts that go right and left
        // n1 to n1 + nleft is n1 and n2 for left partition
        // n1 + nleft to n1 + nleft + nright are n1 and n2 for right partition
        surrogate(dummyNode, n1, n2); // idk, seg fault without
        nodesplit(dummyNode, 1, n1, n2, &nleft, &nright); 
        
        // reset rp values
        for (int q = 0; q < rp.nvar; q++) 
        {
            for (int r = 0; r < rp.n; r++) 
            {
                rp.sorts[q][r] = s_sorts[q * rp.n + r];
                rp.which[r] = s_which[r];
                rp.tempvec[r] = s_tempvec[r];
            }
        }
        
        

        // rp_eval returns the sum squares as the "risk" 
        // pass in 'n' = number of obvs, vector of y values at this node, wt=weight vector, risk and value
        // are only updated with ss and mean repectively (value = node->resp_est)
        
        // if ss < prev_best_ss
        me->primary = bestSplit;
    }
}
