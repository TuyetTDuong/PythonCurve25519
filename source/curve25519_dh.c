/* The MIT License (MIT)
 * 
 * Copyright (c) 2015 mehdi sotoodeh
 * 
 * Permission is hereby granted, free of charge, to any person obtaining 
 * a copy of this software and associated documentation files (the 
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to 
 * permit persons to whom the Software is furnished to do so, subject to 
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include "../include/external_calls.h"
#include "curve25519_mehdi.h"

typedef struct
{
    U_WORD X[K_WORDS];   /* x = X/Z */
    U_WORD Z[K_WORDS];   /*  */
} XZ_POINT;

extern const U_WORD _w_P[K_WORDS];
extern EDP_BLINDING_CTX edp_custom_blinding;

static const U_WORD _w_d[K_WORDS] =
    W256(0x135978A3,0x75EB4DCA,0x4141D8AB,0x00700A4D,0x7779E898,0x8CC74079,0x2B6FFE73,0x52036CEE);
static const U_WORD _w_2d[K_WORDS] = /* 2*d */
    W256(0x26B2F159,0xEBD69B94,0x8283B156,0x00E0149A,0xEEF3D130,0x198E80F2,0x56DFFCE7,0x2406D9DC);
/* x coordinate of base point */
const U8 ecp_BasePoint[32] = { 
    9,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0 };
/* Y = X + X */
void ecp_MontDouble(XZ_POINT *Y, const XZ_POINT *X)
{
    U_WORD A[K_WORDS], B[K_WORDS];
    /*  x2 = (x+z)^2 * (x-z)^2 */
    /*  z2 = ((x+z)^2 - (x-z)^2)*((x+z)^2 + ((A-2)/4)((x+z)^2 - (x-z)^2)) */
    ecp_AddReduce(A, X->X, X->Z);       /* A = (x+z) */
    ecp_SubReduce(B, X->X, X->Z);       /* B = (x-z) */
    ecp_SqrReduce(A, A);                /* A = (x+z)^2 */
    ecp_SqrReduce(B, B);                /* B = (x-z)^2 */
    ecp_MulReduce(Y->X, A, B);          /* x2 = (x+z)^2 * (x-z)^2 */
    ecp_SubReduce(B, A, B);             /* B = (x+z)^2 - (x-z)^2 */
    /* (486662-2)/4 = 121665 */
    ecp_WordMulAddReduce(A, A, 121665, B);
    ecp_MulReduce(Y->Z, A, B);          /* z2 = (B)*((x+z)^2 + ((A-2)/4)(B)) */
}

/* return P = P + Q, Q = 2Q */
void ecp_Mont(XZ_POINT *P, XZ_POINT *Q, IN const U_WORD *Base)
{
    U_WORD A[K_WORDS], B[K_WORDS], C[K_WORDS], D[K_WORDS], E[K_WORDS];
    /* x3 = ((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2*zb     zb=1 */
    /* z3 = ((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2*xb     xb=Base */
    ecp_SubReduce(A, P->X, P->Z);   /* A = x1-z1 */
    ecp_AddReduce(B, P->X, P->Z);   /* B = x1+z1 */
    ecp_SubReduce(C, Q->X, Q->Z);   /* C = x2-z2 */
    ecp_AddReduce(D, Q->X, Q->Z);   /* D = x2+z2 */
    ecp_MulReduce(A, A, D);         /* A = (x1-z1)(x2+z2) */
    ecp_MulReduce(B, B, C);         /* B = (x1+z1)(x2-z2) */
    ecp_AddReduce(E, A, B);         /* E = (x1-z1)(x2+z2) + (x1+z1)(x2-z2) */
    ecp_SubReduce(B, A, B);         /* B = (x1-z1)(x2+z2) - (x1+z1)(x2-z2) */
    ecp_SqrReduce(P->X, E);         /* x3 = ((x1-z1)(x2+z2) + (x1+z1)(x2-z2))^2 */
    ecp_SqrReduce(A, B);            /* A = ((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2 */
    ecp_MulReduce(P->Z, A, Base);   /* z3 = ((x1-z1)(x2+z2) - (x1+z1)(x2-z2))^2*Base */

    /* x4 = (x2+z2)^2 * (x2-z2)^2 */
    /* z4 = ((x2+z2)^2 - (x2-z2)^2)*((x2+z2)^2 + 121665((x2+z2)^2 - (x2-z2)^2)) */
    /* C = (x2-z2) */
    /* D = (x2+z2) */
    ecp_SqrReduce(A, D);            /* A = (x2+z2)^2 */
    ecp_SqrReduce(B, C);            /* B = (x2-z2)^2 */
    ecp_MulReduce(Q->X, A, B);      /* x4 = (x2+z2)^2 * (x2-z2)^2 */
    ecp_SubReduce(B, A, B);         /* B = (x2+z2)^2 - (x2-z2)^2 */
    ecp_WordMulAddReduce(A, A, 121665, B);
    ecp_MulReduce(Q->Z, A, B);      /* z4 = B*((x2+z2)^2 + 121665*B) */
}

/* Constant-time measure: */
/* Use different set of parameters for bit=0 or bit=1 with no conditional jump */
/* */
#define ECP_MONT(n) j = (k >> n) & 1; ecp_Mont(PP[j], QP[j], X)

/* -------------------------------------------------------------------------- */
/* Return point Q = k*P */
/* K in a little-endian byte array */
void ecp_PointMultiply(
    OUT U8 *PublicKey, 
    IN const U8 *BasePoint, 
    IN const U8 *SecretKey, 
    IN int len)
{
    int i, j, k;
    U_WORD X[K_WORDS];
    XZ_POINT P, Q, *PP[2], *QP[2];

    ecp_BytesToWords(X, BasePoint);

    /* 1: P = (2k+1)G, Q = (2k+2)G */
    /* 0: Q = (2k+1)G, P = (2k)G */

    /* Find first non-zero bit */
    while (len-- > 0)
    {
        k = SecretKey[len];
        for (i = 0; i < 8; i++, k <<= 1)
        {
            /* P = kG, Q = (k+1)G */
            if (k & 0x80)
            {
                /* We have first non-zero bit 
                // This is always bit 254 for keys created according to the spec.
                // Start with randomized base point 
                */

                ecp_Add(P.Z, X, edp_custom_blinding.zr);    /* P.Z = random */
                ecp_MulReduce(P.X, X, P.Z);
                ecp_MontDouble(&Q, &P);

                PP[1] = &P; PP[0] = &Q;
                QP[1] = &Q; QP[0] = &P;

                /* Everything we reference in the below loop are on the stack
                // and already touched (cached) 
                */

                while (++i < 8) { k <<= 1; ECP_MONT(7); }
                while (len > 0)
                {
                    k = SecretKey[--len];
                    ECP_MONT(7);
                    ECP_MONT(6);
                    ECP_MONT(5);
                    ECP_MONT(4);
                    ECP_MONT(3);
                    ECP_MONT(2);
                    ECP_MONT(1);
                    ECP_MONT(0);
                }

                ecp_Inverse(Q.Z, P.Z);
                ecp_MulMod(X, P.X, Q.Z);
                ecp_WordsToBytes(PublicKey, X);
                return;
            }
        }
    }
    /* K is 0 */
    mem_fill(PublicKey, 0, 32);
}
void ecp_recoverY(OUT U8 *OUT_Y,
                      IN const U8 *PublicKey)
{
    U_WORD X[K_WORDS], Y[K_WORDS], K[K_WORDS], T[K_WORDS];
    ecp_BytesToWords(X, PublicKey);

    ecp_MulReduce(T,X,X);
    ecp_MulReduce(Y,T,X);
    ecp_WordMulAddReduce(K,X,486662,T);
    ecp_AddReduce(Y,Y,K);
    ecp_WordsToBytes(OUT_Y, Y);
}

/* -------------------------------------------------------------------------- */
/* Return point Z = P + Q */
void ecp_PointAddition(
    OUT U8 *PublicKeySum,
    IN const U8 *PublicKeyP,
    IN const U8 *PublicKeyQ,
    IN const U8 *BasePoint)
{
    int i, j, k;
    U_WORD X[K_WORDS], Z[K_WORDS], B[K_WORDS], XY[K_WORDS], ZY[K_WORDS];
    XZ_POINT P, Q;

    ecp_BytesToWords(B, BasePoint);
    ecp_BytesToWords(X, PublicKeyP);
    ecp_BytesToWords(Z, PublicKeyQ);
    U8 PublicKeyQY, PublicKeyPY;
    /*ecp_Add(P.Z, X, edp_custom_blinding.zr);
    ecp_MulReduce(P.X, X, P.Z);

    ecp_Add(Q.Z, Z, edp_custom_blinding.zr);
    ecp_MulReduce(Q.X, Z, Q.Z);*/

    ecp_recoverY(PublicKeyP,PublicKeyP);
    ecp_recoverY(PublicKeyQ,PublicKeyQ);
    ecp_BytesToWords(XY, PublicKeyP);
    ecp_BytesToWords(ZY, PublicKeyQ);

    U_WORD a[K_WORDS], b[K_WORDS], c[K_WORDS], d[K_WORDS], e[K_WORDS], A[K_WORDS];
   //(x2y1)
    ecp_MulReduce(a,Z,XY);
    //x1y2
    ecp_MulReduce(b,X,ZY);
    //x2y1-x1y2
    ecp_SubReduce(a,a,b);
    //(x2y1-x1y2)^2
    ecp_MulReduce(b,a,a);
    //x1x2
    ecp_MulReduce(c,X,Z);
    //(x2-x1)
    ecp_SubReduce(d,Z,X);
    //(x2-x1)^2
    ecp_MulReduce(e,d,d);
    //x1x2 (x2-x1)^2
    ecp_MulReduce(c,c,e);

    /* Calculate:  sqrt(b/c) = u*v^3 * (u*v^7)^((p-5)/8) */

    ecp_SqrReduce(B, c);
    ecp_MulReduce(A, b, B);
    ecp_MulReduce(A, A, c);         /* a = u*v^3 */
    ecp_SqrReduce(B, B);            /* b = v^4 */
    ecp_MulReduce(B, A, B);         /* b = u*v^7 */
    ecp_ModExp2523(B, B);
    ecp_MulReduce(X, B, A);

    /* Check if we have correct sqrt, else, multiply by sqrt(-1) */

    /*ecp_SqrReduce(B, X);
    ecp_MulReduce(B, B, c);
    ecp_SubReduce(B, B, b);
    ecp_Mod(B);*/

    ecp_WordsToBytes(PublicKeySum, X);



}

/* -------------------------------------------------------------------------- */
/* Return point Z = P + Q */
void ecp_PointAdditionTest(
    OUT U8 *PublicKeySum,
    IN const U8 *PublicKeyP,
    IN const U8 *PublicKeyQ,
    IN const U8 *BasePoint)
{

    U_WORD X[K_WORDS], Z[K_WORDS], B[K_WORDS], T[K_WORDS];
    XZ_POINT P, Q;

    ecp_BytesToWords(B, BasePoint);
    ecp_BytesToWords(X, PublicKeyP);
    ecp_BytesToWords(Z, PublicKeyQ);

    /*ecp_Add(P.Z, X, edp_custom_blinding.zr);
    ecp_MulReduce(P.X, X, P.Z);

    ecp_Add(Q.Z, Z, edp_custom_blinding.zr);
    ecp_MulReduce(Q.X, Z, Q.Z);


    ecp_Mont(&P, &Q, B);
    ecp_Inverse(Q.Z, P.Z);
    ecp_MulMod(X, P.X, Q.Z);
    ecp_WordsToBytes(PublicKeySum, X);*/

    ecp_AddReduce(T,X,Z);
    ecp_WordsToBytes(PublicKeySum, T);
    ecp_PointMultiply(PublicKeySum, ecp_BasePoint, PublicKeySum, 32);

    ecp_BytesToWords(Z, PublicKeySum);
    ecp_Add(Q.Z, Z, edp_custom_blinding.zr);
    ecp_MulReduce(Q.X, Z, Q.Z);

    //Convert  x25519 point to  ed25519 point
    //y = (u - 1)/(u + 1) = (X - Z )/(X + Z)


    // Converting Q to ed25519
    ecp_SubReduce(Q.X, Q.X, Q.Z);
    ecp_AddReduce(Q.Z, Q.X, Q.Z);
    ecp_Inverse(Q.Z, Q.Z);
    ecp_MulMod(Q.X, Q.X, Q.Z);
    ecp_WordsToBytes(PublicKeyQ, Q.X);

    Ext_POINT K, M, S;
    int i,j;
    //Affine_POINT L;
    //ed25519_PackPoint(pubKey, Q.X, L.x[0]);

    //j = ecp_DecodeInt(M.y, PublicKeyQ);
    //ed25519_CalculateX(M.x, M.y, ~j);
    //ecp_MulMod(M.t, M.x, M.y);
    //ecp_SetValue(M.z, 1);



    //Convert ed25519 point to x25519 point

    //u = (1 + y)/(1 - y) = (Z + Y)/(Z - Y)
    ecp_AddReduce(M.t, M.z, M.y);
    ecp_SubReduce(M.z, M.z, M.y);
    ecp_Inverse(M.z, M.z);
    ecp_MulMod(M.t, M.t, M.z);
    ecp_WordsToBytes(PublicKeySum, M.t);



}

/* -- DH key exchange interfaces ----------------------------------------- */

/* Return R = a*P where P is curve25519 base point */
void x25519_BasePointMultiply(OUT U8 *r, IN const U8 *sk)
{
    Ext_POINT S;
    U_WORD t[K_WORDS];

    ecp_BytesToWords(t, sk);
    edp_BasePointMult(&S, t, edp_custom_blinding.zr);

    /* Convert ed25519 point to x25519 point */
    
    /* u = (1 + y)/(1 - y) = (Z + Y)/(Z - Y) */

    ecp_AddReduce(S.t, S.z, S.y);
    ecp_SubReduce(S.z, S.z, S.y);
    ecp_Inverse(S.z, S.z);
    ecp_MulMod(S.t, S.t, S.z);
    ecp_WordsToBytes(r, S.t);
}

/* Return public key associated with sk */
void curve25519_dh_CalculatePublicKey_fast(
    unsigned char *pk,          /* [32-bytes] OUT: Public key */
    unsigned char *sk)          /* [32-bytes] IN/OUT: Your secret key */
{
    ecp_TrimSecretKey(sk);
    /* Use faster method */
    x25519_BasePointMultiply(pk, sk);
}

/* Return public key associated with sk */
void curve25519_dh_CalculatePublicKey(
    unsigned char *pk,          /* [32-bytes] OUT: Public key */
    unsigned char *sk)          /* [32-bytes] IN/OUT: Your secret key */
{
    //ecp_TrimSecretKey(sk);
    ecp_PointMultiply(pk, ecp_BasePoint, sk, 32);
}

/* Create a shared secret */
void curve25519_dh_CreateSharedKey(
    unsigned char *shared,      /* [32-bytes] OUT: Created shared key */
    const unsigned char *pk,    /* [32-bytes] IN: Other side's public key */
    unsigned char *sk)          /* [32-bytes] IN/OUT: Your secret key */
{
    //ecp_TrimSecretKey(sk);
    ecp_PointMultiply(shared, pk, sk, 32);
}

void edp_AddPointWithoutPreCompute(Ext_POINT *r, const Ext_POINT *p, const Ext_POINT *q)
{
    U_WORD a1[K_WORDS], a2[K_WORDS], b1[K_WORDS], b2[K_WORDS], c[K_WORDS], d[K_WORDS], e[K_WORDS];

    ecp_SubReduce(a1, p->y, p->x);          /* A = (Y1-X1)*(Y2-X2) */
    ecp_SubReduce(a2, q->y, q->x);
    ecp_MulReduce(a1, a1, a2);
    ecp_AddReduce(b1, p->y, p->x);          /* B = (Y1+X1)*(Y2+X2) */
    ecp_AddReduce(b2, q->y, q->x);
    ecp_MulReduce(b1, b1, b2);

    ecp_MulReduce(c, p->t, _w_2d);          /* C = T1*2d*T2 */
    ecp_MulReduce(c, q->t, c);
    ecp_AddReduce(d, p->z, p->z);              /* D = Z1*2*Z2 */
    ecp_MulReduce(d, q->z, d);
    ecp_SubReduce(e, b1, a1);                 /* E = B-A */
    ecp_AddReduce(b1, b1, a1);                 /* H = B+A */
    ecp_SubReduce(a1, d, c);                 /* F = D-C */
    ecp_AddReduce(d, d, c);                 /* G = D+C */

    ecp_MulReduce(r->x, e, a1);              /* E*F */
    ecp_MulReduce(r->y, b1, d);              /* H*G */
    ecp_MulReduce(r->t, e, b1);              /* E*H */
    ecp_MulReduce(r->z, d, a1);              /* G*F */
}

/* Return a sum of two public keys */
void curve25519_dh_CalculateSumTwoPublicKeys(
    unsigned char *pkS,          /* [32-bytes] OUT: Public key */
    unsigned char *pkP,           /* [32-bytes] IN: Your first public key */
    unsigned char *pkQ)           /* [32-bytes] IN: Your second public key */
{

    /*ecp_PointAdditionTest(pkS, pkP, pkQ, ecp_BasePoint);
    ecp_PointMultiply(pkS, ecp_BasePoint, pkS, 32);

    //ecp_PointMultiply(pkQ, ecp_BasePoint, pkQ, 32);
    //ecp_PointMultiply(pkP, ecp_BasePoint, pkP, 32);


    //ecp_PointMultiply(pkQ, ecp_BasePoint, pkQ, 32);
    ecp_PointMultiply(pkP, ecp_BasePoint, pkP, 32);

    U_WORD X[K_WORDS], Z[K_WORDS], B[K_WORDS], T[K_WORDS];
    XZ_POINT P, Q;
    ecp_BytesToWords(X, pkP);
    ecp_BytesToWords(Z, pkQ);
    ecp_Add(P.Z, X, edp_custom_blinding.zr);
    ecp_MulReduce(P.X, X, P.Z);

    ecp_Add(Q.Z, Z, edp_custom_blinding.zr);
    ecp_MulReduce(Q.X, Z, Q.Z);
    //Convert  x25519 point to  ed25519 point
    //y = (u - 1)/(u + 1) = (X - Z )/(X + Z)


    // Converting P to ed25519
    ecp_SubReduce(P.X, P.X, P.Z);
    ecp_AddReduce(P.Z, P.X, P.Z);
    ecp_Inverse(P.Z, P.Z);
    ecp_MulMod(P.X, P.X, P.Z);
    ecp_WordsToBytes(pkP, P.X);

    // Converting Q to ed25519
    ecp_SubReduce(Q.X, Q.X, Q.Z);
    ecp_AddReduce(Q.Z, Q.X, Q.Z);
    ecp_Inverse(Q.Z, Q.Z);
    ecp_MulMod(Q.X, Q.X, Q.Z);
    ecp_WordsToBytes(pkQ, Q.X);

    Ext_POINT K, M, S;
    int i,j;

    i = ecp_DecodeInt(K.y, pkP);
    ed25519_CalculateX(K.x, K.y, ~i);
    ecp_MulMod(K.t, K.x, K.y);
    ecp_SetValue(K.z, 1);

    j = ecp_DecodeInt(M.y, pkQ);
    ed25519_CalculateX(M.x, M.y, ~j);
    ecp_MulMod(M.t, M.x, M.y);
    ecp_SetValue(M.z, 1);

    edp_AddPointWithoutPreCompute(&S,&K,&M);

    //Convert ed25519 point to x25519 point

    //u = (1 + y)/(1 - y) = (Z + Y)/(Z - Y)
    ecp_AddReduce(S.t, S.z, S.y);
    ecp_SubReduce(S.z, S.z, S.y);
    ecp_Inverse(S.z, S.z);
    ecp_MulMod(S.t, S.t, S.z);
    ecp_WordsToBytes(pkS, S.t);*/
    ecp_PointAdditionTest(pkS, pkP, pkQ, ecp_BasePoint);
    //ecp_PointMultiply(pkS, ecp_BasePoint, pkS, 32);

     //ecp_PointMultiply(pkP, ecp_BasePoint, pkP, 32);
     //ecp_PointAddition(pkS, pkP, pkQ, ecp_BasePoint);
}


/* Return a sum of field elements */
void curve25519_dh_CalculateSumFieldElements(
    unsigned char *Sum,          /* [32-bytes] OUT: Public key */
    unsigned char *P1,           /* [32-bytes] IN: Your first public key */
    unsigned char *P2)           /* [32-bytes] IN: Your second public key */
{
    U_WORD X[K_WORDS], Z[K_WORDS], S[K_WORDS];

    ecp_BytesToWords(X, P1);
    ecp_BytesToWords(Z, P2);

    ecp_AddReduce(S,X,Z);
    ecp_WordsToBytes(Sum, S);

}

/* Return a product of field elements */
void curve25519_dh_CalculateProductFieldElements(
    unsigned char *Sum,          /* [32-bytes] OUT: Public key */
    unsigned char *P1,           /* [32-bytes] IN: Your first public key */
    unsigned char *P2)           /* [32-bytes] IN: Your second public key */
{
    U_WORD X[K_WORDS], Z[K_WORDS], S[K_WORDS];

    ecp_BytesToWords(X, P1);
    ecp_BytesToWords(Z, P2);

    ecp_MulReduce(S,X,Z);
    ecp_WordsToBytes(Sum, S);

}