#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **make_mat(int row, int col)
{
    double **mat;
    int i;

    mat = (double **) malloc(sizeof(double *) * row);
    if(mat==NULL) {
        return NULL;
    }
    for(i=0; i<row; i++) {
        mat[i] = (double *) calloc(sizeof(double), col);
        if(mat[i] == NULL) {
            return NULL;
        }
    }
    return mat;
}

double **pivot_gaus(double **mat, int size)
{
    double **inv;
    int i,j,k,m;
    double akk, ajk, temp, max;

    inv = (double **) malloc(sizeof(double *)*size);
    if(inv==NULL) {
        return NULL;
    }
    for(i=0; i<size; i++) {
        inv[i] = (double *) malloc(sizeof(double) * size);
        if(inv[i] == NULL) {
            return NULL;
        }
    }

    for(k=0; k<size; k++) {
        max = 0.0;
        /* find pivot */
        for(j=k; j<size; j++) {
            if(max < fabs(mat[j][k])) {
                max = fabs(mat[j][k]);
                m = j;
            }
        }
        /* 行入れ替え */
        if(m!=k) {
            for(j=0; j<size*2; j++) {
                temp = mat[k][j];
                mat[k][j] = mat[m][j];
                mat[m][j] = temp;
            }
        }
        /* mat[k][k] = 1.0 になるよう k行を割り算する */
        akk = mat[k][k];
        for(j=0; j<size*2; j++) {
            mat[k][j] /= akk;
        }
        /* k行以外の行からk行*ajkを減算して、k列を0にする */
        for(j=0; j<size; j++) {
            if(j==k) continue;
            ajk = mat[j][k];
            for(i=0; i<size*2; i++) {
                mat[j][i] -= ajk*mat[k][i];
            }
        }
    }

    for(j=0; j<size; j++) {
        for(i=0; i<size; i++) {
            inv[j][i] = mat[j][i+size];
        }
    }
    return inv;
}

const double mat_tmp[4][4] = {
    {1.0,3.0,6.0,9.0},
    {2.0,5.0,1.0,10.0},
    {3.0,7.0,2.0,8.0},
    {4.0,1.0,3.0,5.0}
};

#define R_SIZE 400.0
#define LIGHT_SRC 400.0
#define REFLECT 0.8
#define IMG_SIZE_X 160
#define IMG_SIZE_Y 160

double calc_bright(double a, double b, double p, double q, double ref, double x, double y, double z, double eta)
{
    double g,h;
    double nlen;
    double ss[3], v[3];
    double vlen;
    double na;
    double delta;

    g = (-p*sin(a)*cos(b)-q*sin(a)*sin(b)+cos(a))/sqrt(1.0+p*p+q*q);

    nlen = 1.0+p*p+q*q;
    ss[0] = ((p*p-q*q-1.0)*sin(a)*cos(b) + 2.0*p*q*sin(a)*sin(b) - 2.0*p*cos(a))/nlen;
    ss[1] = (2.0*p*q*sin(a)*cos(b) + (q*q-p*p-1.0)*sin(a)*sin(b) - 2.0*q*cos(a))/nlen;
    ss[2] = (-2.0*p*sin(a)*cos(b) - 2.0*q*sin(a)*sin(b) + (1.0-p*p-q*q)*cos(a))/nlen;
    vlen = sqrt(x*x+y*y+(z-eta)*(z-eta));
    v[0] = x/vlen;
    v[1] = y/vlen;
    v[2] = (eta-z)/vlen;

    na = ss[0]*v[0]+ss[1]*v[1]+ss[2]*v[2];

    delta = acos(na);

    h = ref * exp(-delta*delta);

    return LIGHT_SRC*g*h;
}

int main(int argc, char *argv[])
{
    int i, j;
    double **mat_tr, **mat_tr_inv;
    double eta = 10000.0;
    double d = 1000.0;
    unsigned char *dst;
    double x, y, z, w;
    double p,q;
    int dx, dy;
    //double alpha = M_PI/3.0;	/* α = 60°*/
    //double beta = M_PI * 75.0 / 180.0;	/* β = 75°*/
    double alpha = M_PI * 120.0 / 180.0;	/* α = 80°*/
    double beta = M_PI * 30.0 / 180.0;	/* β = 90°*/
    double theta;
    int r;
    int max_r;
    double bright;

    if(argc > 2) {
        alpha = atof(argv[1]) * M_PI / 180.0;
        beta = atof(argv[2]) * M_PI / 180.0;
    }

    dst = (unsigned char*) calloc(sizeof(unsigned char), IMG_SIZE_X*IMG_SIZE_Y);
    if(dst == NULL) {
        fprintf(stderr, "malloc err\n");
        return -1;
    }

    max_r = (int) sqrt(eta*eta-R_SIZE*R_SIZE)*R_SIZE/eta;
    for(r=0; r<=max_r; r++) {
        for(theta=0; theta<2*M_PI; theta += (M_PI/180.0)) {
            x = r*cos(theta);
            y = r*sin(theta);
            z = sqrt(R_SIZE*R_SIZE-r*r);
            w = (eta-z)/d;
            dx = (int) (x/w + IMG_SIZE_X/2);
            dy = (int) (y/w + IMG_SIZE_Y/2);

            p = -x/z;
            q = -y/z;

            bright = calc_bright(alpha, beta, p, q, REFLECT, x, y, z, eta);

            if(bright > 255) bright = 255;
            if(bright < 0) bright = 0;

            if(dx > IMG_SIZE_X-1) dx = IMG_SIZE_X-1;
            if(dx < 0) dx=0;
            if(dy > IMG_SIZE_Y-1) dy = IMG_SIZE_Y-1;
            if(dy < 0) dy = 0;
            //printf("%d %d %f\n", dx, dy, bright);
            dst[dy*IMG_SIZE_X+dx] = (unsigned char) bright;
        }
    }

    {
        FILE *fp;

        fp = fopen("img.pgm", "w");

        fprintf(fp, "P2 %d %d 255\n", IMG_SIZE_X, IMG_SIZE_Y);

        for(dy=0; dy<IMG_SIZE_Y; dy++) {
            for(dx=0; dx<IMG_SIZE_X; dx++) {
                fprintf(fp, "%d\n", dst[dy*IMG_SIZE_X+dx]);
            }
        }
        fclose(fp);
    }

    free(dst);

    return 0;
}
