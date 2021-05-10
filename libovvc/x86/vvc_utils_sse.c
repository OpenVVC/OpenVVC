#include <emmintrin.h>
#include <smmintrin.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

typedef int16_t TCoeff;
typedef int16_t TMatrixCoeff;

#define SHIFTV 7
#define SHIFTH 10
#define SIZEV 32
#define SIZEH 2

void afficherVecteur4SSE128(__m128i cible)
{
    int32_t test[4];
    _mm_store_si128((__m128i *)test, cible);

//    printf(" ________ ________ ________ ________\n|        |        |        |        |\n");
    printf("| %06d | %06d | %06d | %06d |\n",  test[0], test[1], test[2], test[3]);
//    printf("|________|________|________|________|\n");
}

void afficherVecteur8SSE128(__m128i cible)
{
    int16_t test[8];
    _mm_store_si128((__m128i *)test, cible);

//    printf(" ________ ________ ________ ________ ________ ________ ________ ________\n|        |        |        |        |        |        |        |        |\n");
    printf("| %06d | %06d | %06d | %06d | %06d | %06d | %06d | %06d |\n", test[0], test[1], test[2], test[3], test[4], test[5], test[6], test[7]);
//    printf("|________|________|________|________|________|________|________|________|\n");
}

void afficherTableau2D(int v, int h, TCoeff * m){
    int i,j;
    for (i = 0; i < v; ++i) {
        for (j = 0; j < h; ++j) {
            printf("%i ",m[i*h+j]);
        }
        printf("\n");
    }
    printf("\n");
}

void inverse_sse2_B4(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
    __m128i x1, x2; //Contient le vecteur à transformer
    __m128i d, d2, d3, d4;	//Contient coefficient DCT ou DST
    __m128i result[line]; //Variables pour calculs (result[] il faudrait mettre line ou la taille max 32)

    int nbstore = line/2;

    d = _mm_load_si128((__m128i *)&(iT[0]));
    d2 = _mm_load_si128((__m128i *)&(iT[8])); //8 car [taille coef DCT ou DST 4] x [2 (le milieu de la hauteur)]
    d3 = _mm_unpacklo_epi16(d,d2);
    d4 = _mm_unpackhi_epi16(d,d2);
    for(int i = 0; i < line; ++i){
        x1 = _mm_unpacklo_epi16(_mm_set1_epi16(src[i]),_mm_set1_epi16(src[2*src_stride+i]));
        x2 = _mm_unpacklo_epi16(_mm_set1_epi16(src[src_stride+i]),_mm_set1_epi16(src[3*src_stride+i]));

        result[i] = _mm_add_epi32(_mm_madd_epi16(x1, d3),_mm_madd_epi16(x2, d4));

        result[i] = _mm_add_epi32(result[i], _mm_set1_epi32(1<<(shift-1)));

        result[i] = _mm_srai_epi32(result[i], shift);
    }

    for(int i = 0; i<nbstore; i++){
        result[i] = _mm_packs_epi32(result[2*i], result[2*i+1]); //clip pour repasser en 16
        _mm_store_si128((__m128i *)&(dst[i*8]), result[i]);	//dst[i*8] car result contient 8 résultat
    }
}


void inverse_sse2_B8(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
  __m128i add = _mm_set1_epi32(1 << (shift - 1));
  TCoeff *_src = (TCoeff *)src;
  TCoeff *_dst = dst;
  uint8_t i, k;
  __m128i x[8], d[8], m[8], di[8], a[8], b[4], r[2], res[8];
  for (i = 0; i < line>>3; i++) {
    x[0]=_mm_load_si128((__m128i*)(_src + 0 * src_stride));
    x[1]=_mm_load_si128((__m128i*)(_src + 1 * src_stride));
    x[2]=_mm_load_si128((__m128i*)(_src + 2 * src_stride));
    x[3]=_mm_load_si128((__m128i*)(_src + 3 * src_stride));
    x[4]=_mm_load_si128((__m128i*)(_src + 4 * src_stride));
    x[5]=_mm_load_si128((__m128i*)(_src + 5 * src_stride));
    x[6]=_mm_load_si128((__m128i*)(_src + 6 * src_stride));
    x[7]=_mm_load_si128((__m128i*)(_src + 7 * src_stride));

    for (k = 0; k < 8; k++) {
      d[0] = _mm_set1_epi16(iT[k + 0 * 8]);
      d[1] = _mm_set1_epi16(iT[k + 1 * 8]);
      d[2] = _mm_set1_epi16(iT[k + 2 * 8]);
      d[3] = _mm_set1_epi16(iT[k + 3 * 8]);
      d[4] = _mm_set1_epi16(iT[k + 4 * 8]);
      d[5] = _mm_set1_epi16(iT[k + 5 * 8]);
      d[6] = _mm_set1_epi16(iT[k + 6 * 8]);
      d[7] = _mm_set1_epi16(iT[k + 7 * 8]);

      m[0] = _mm_unpacklo_epi16(x[0],  x[1]);
      m[1] = _mm_unpacklo_epi16(x[2],  x[3]);
      m[2] = _mm_unpacklo_epi16(x[4],  x[5]);
      m[3] = _mm_unpacklo_epi16(x[6],  x[7]);

      m[4] = _mm_unpackhi_epi16(x[0],  x[1]);
      m[5] = _mm_unpackhi_epi16(x[2],  x[3]);
      m[6] = _mm_unpackhi_epi16(x[4],  x[5]);
      m[7] = _mm_unpackhi_epi16(x[6],  x[7]);

      di[0] = _mm_unpacklo_epi16(d[0],  d[1]);
      di[1] = _mm_unpacklo_epi16(d[2],  d[3]);
      di[2] = _mm_unpacklo_epi16(d[4],  d[5]);
      di[3] = _mm_unpacklo_epi16(d[6],  d[7]);

      di[4] = _mm_unpackhi_epi16(d[0],  d[1]);
      di[5] = _mm_unpackhi_epi16(d[2],  d[3]);
      di[6] = _mm_unpackhi_epi16(d[4],  d[5]);
      di[7] = _mm_unpackhi_epi16(d[6],  d[7]);


      a[0] = _mm_madd_epi16(m[0], di[0]);
      a[1] = _mm_madd_epi16(m[1], di[1]);
      a[2] = _mm_madd_epi16(m[2], di[2]);
      a[3] = _mm_madd_epi16(m[3], di[3]);

      a[4] = _mm_madd_epi16(m[4], di[4]);
      a[5] = _mm_madd_epi16(m[5], di[5]);
      a[6] = _mm_madd_epi16(m[6], di[6]);
      a[7] = _mm_madd_epi16(m[7], di[7]);

      b[0] = _mm_add_epi32(a[0], a[2]);
      b[1] = _mm_add_epi32(a[1], a[3]);

      b[2] = _mm_add_epi32(a[4], a[6]);
      b[3] = _mm_add_epi32(a[5], a[7]);

      r[0] = _mm_add_epi32(b[0], b[1]);
      r[1] = _mm_add_epi32(b[2], b[3]);

      r[0] = _mm_add_epi32(r[0], add);
      r[1] = _mm_add_epi32(r[1], add);

      r[0] = _mm_srai_epi32(r[0], shift);
      r[1] = _mm_srai_epi32(r[1], shift);

      res[k] = _mm_packs_epi32(r[0], r[1]);
    }

        m[0] = _mm_unpacklo_epi16(res[0], res[1]);
        m[1] = _mm_unpackhi_epi16(res[0], res[1]);
        m[2] = _mm_unpacklo_epi16(res[2], res[3]);
        m[3] = _mm_unpackhi_epi16(res[2], res[3]);
        m[4] = _mm_unpacklo_epi16(res[4], res[5]);
        m[5] = _mm_unpackhi_epi16(res[4], res[5]);
        m[6] = _mm_unpacklo_epi16(res[6], res[7]);
        m[7] = _mm_unpackhi_epi16(res[6], res[7]);

        res[0] = _mm_unpacklo_epi32(m[0], m[2]);
        res[1] = _mm_unpackhi_epi32(m[0], m[2]);
        res[2] = _mm_unpacklo_epi32(m[1], m[3]);
        res[3] = _mm_unpackhi_epi32(m[1], m[3]);
        res[4] = _mm_unpacklo_epi32(m[4], m[6]);
        res[5] = _mm_unpackhi_epi32(m[4], m[6]);
        res[6] = _mm_unpacklo_epi32(m[5], m[7]);
        res[7] = _mm_unpackhi_epi32(m[5], m[7]);

        m[0] = _mm_unpacklo_epi64(res[0], res[4]);
        m[1] = _mm_unpackhi_epi64(res[0], res[4]);
        m[2] = _mm_unpacklo_epi64(res[1], res[5]);
        m[3] = _mm_unpackhi_epi64(res[1], res[5]);
        m[4] = _mm_unpacklo_epi64(res[2], res[6]);
        m[5] = _mm_unpackhi_epi64(res[2], res[6]);
        m[6] = _mm_unpacklo_epi64(res[3], res[7]);
        m[7] = _mm_unpackhi_epi64(res[3], res[7]);

        _mm_store_si128((__m128i *)&_dst[0],  m[0]);
        _mm_store_si128((__m128i *)&_dst[8],  m[1]);
        _mm_store_si128((__m128i *)&_dst[16], m[2]);
        _mm_store_si128((__m128i *)&_dst[24], m[3]);
        _mm_store_si128((__m128i *)&_dst[32], m[4]);
        _mm_store_si128((__m128i *)&_dst[40], m[5]);
        _mm_store_si128((__m128i *)&_dst[48], m[6]);
        _mm_store_si128((__m128i *)&_dst[56], m[7]);
        _src+=8;
        _dst+=64;
    }

    if (!(line&0x7)) {
      return ;
    }

    d[0] = _mm_load_si128((__m128i *)&(iT[0*8]));
    d[1] = _mm_load_si128((__m128i *)&(iT[1*8]));
    d[2] = _mm_load_si128((__m128i *)&(iT[2*8]));
    d[3] = _mm_load_si128((__m128i *)&(iT[3*8]));
    d[4] = _mm_load_si128((__m128i *)&(iT[4*8]));
    d[5] = _mm_load_si128((__m128i *)&(iT[5*8]));
    d[6] = _mm_load_si128((__m128i *)&(iT[6*8]));
    d[7] = _mm_load_si128((__m128i *)&(iT[7*8]));

    for(i = 0; i < (line&0x7); ++i){
      x[0] = _mm_set1_epi16(_src[0*src_stride+i]);
      x[1] = _mm_set1_epi16(_src[1*src_stride+i]);
      x[2] = _mm_set1_epi16(_src[2*src_stride+i]);
      x[3] = _mm_set1_epi16(_src[3*src_stride+i]);
      x[4] = _mm_set1_epi16(_src[4*src_stride+i]);
      x[5] = _mm_set1_epi16(_src[5*src_stride+i]);
      x[6] = _mm_set1_epi16(_src[6*src_stride+i]);
      x[7] = _mm_set1_epi16(_src[7*src_stride+i]);

      m[0] = _mm_unpacklo_epi16(x[0],  x[1]);
      m[1] = _mm_unpacklo_epi16(x[2],  x[3]);
      m[2] = _mm_unpacklo_epi16(x[4],  x[5]);
      m[3] = _mm_unpacklo_epi16(x[6],  x[7]);

      m[4] = _mm_unpackhi_epi16(x[0],  x[1]);
      m[5] = _mm_unpackhi_epi16(x[2],  x[3]);
      m[6] = _mm_unpackhi_epi16(x[4],  x[5]);
      m[7] = _mm_unpackhi_epi16(x[6],  x[7]);

      di[0] = _mm_unpacklo_epi16(d[0],  d[1]);
      di[1] = _mm_unpacklo_epi16(d[2],  d[3]);
      di[2] = _mm_unpacklo_epi16(d[4],  d[5]);
      di[3] = _mm_unpacklo_epi16(d[6],  d[7]);

      di[4] = _mm_unpackhi_epi16(d[0],  d[1]);
      di[5] = _mm_unpackhi_epi16(d[2],  d[3]);
      di[6] = _mm_unpackhi_epi16(d[4],  d[5]);
      di[7] = _mm_unpackhi_epi16(d[6],  d[7]);

      a[0] = _mm_madd_epi16(m[0], di[0]);
      a[1] = _mm_madd_epi16(m[1], di[1]);
      a[2] = _mm_madd_epi16(m[2], di[2]);
      a[3] = _mm_madd_epi16(m[3], di[3]);

      a[4] = _mm_madd_epi16(m[4], di[4]);
      a[5] = _mm_madd_epi16(m[5], di[5]);
      a[6] = _mm_madd_epi16(m[6], di[6]);
      a[7] = _mm_madd_epi16(m[7], di[7]);

      b[0] = _mm_add_epi32(a[0], a[2]);
      b[1] = _mm_add_epi32(a[1], a[3]);

      b[2] = _mm_add_epi32(a[4], a[6]);
      b[3] = _mm_add_epi32(a[5], a[7]);

      r[0] = _mm_add_epi32(b[0], b[1]);
      r[1] = _mm_add_epi32(b[2], b[3]);

      r[0] = _mm_add_epi32(r[0], add);
      r[1] = _mm_add_epi32(r[1], add);

      r[0] = _mm_srai_epi32(r[0], shift);
      r[1] = _mm_srai_epi32(r[1], shift);

      _mm_store_si128((__m128i *)&(_dst[i*8]), _mm_packs_epi32(r[0], r[1])); //dst[i*8] car result contient 8 résultat
    }
}

void inverse_sse2_B16(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
    __m128i x[8];
    __m128i d[8];
    __m128i vhi[8], vlo[8], vh[16], vl[16], result[line][4];
    
    for(int l = 0; l < 2; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_set1_epi32(0);
            result[i][2*l+1] = _mm_set1_epi32(0);
        }
    }

    for(int l = 0; l < 2; ++l){
        for(int k = 0; k < 2; ++k){
            for(int i = 0; i < 8; ++i){
                d[i] = _mm_load_si128((__m128i *)&(iT[i*16+k*128+l*8]));
            }

            for(int i = 0; i < line; ++i){

                for(int j = 0; j < 8; ++j){
                    x[j] = _mm_set1_epi16(src[j*src_stride+i+k*8*src_stride]);

                    vhi[j] = _mm_mulhi_epi16(x[j], d[j]);
                    vlo[j] = _mm_mullo_epi16(x[j], d[j]);
                    vl[j+8*k] = _mm_unpacklo_epi16(vlo[j], vhi[j]);
                    vh[j+8*k] = _mm_unpackhi_epi16(vlo[j], vhi[j]);
                }

                for(int j = 0; j < 4; ++j){
                    result[i][2*l] = _mm_add_epi32(result[i][2*l], _mm_add_epi32(vl[2*j+8*k],vl[2*j+1+8*k]));
                    result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], _mm_add_epi32(vh[2*j+8*k],vh[2*j+1+8*k]));
                }
            }
        }
    }

    for(int l = 0; l < 2; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_add_epi32(result[i][2*l], _mm_set1_epi32(1<<(shift-1)));
            result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], _mm_set1_epi32(1<<(shift-1)));
            result[i][2*l] = _mm_srai_epi32(result[i][2*l], shift);
            result[i][2*l+1] = _mm_srai_epi32(result[i][2*l+1], shift);
        }
    }

    for(int l = 0; l < 2; ++l){
        for(int i = 0; i < line; i++){
            result[i][l] = _mm_packs_epi32(result[i][2*l], result[i][2*l+1]); //clip pour repasser en 16
            _mm_store_si128((__m128i *)&(dst[i*16+l*8]), result[i][l]);
        }
    }
}

void inverse_sse2_B32(const TCoeff *src, TCoeff *dst, int src_stride, int shift, int line, const TMatrixCoeff* iT){
    __m128i x[8];
    __m128i d[8];
    __m128i vhi[8], vlo[8], vh[32], vl[32], result[line][8];

    for(int l = 0; l < 4; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_set1_epi32(0);
            result[i][2*l+1] = _mm_set1_epi32(0);
        }
    }

    for(int l = 0; l < 4; ++l){
        for(int k = 0; k < 4; ++k){
            for(int i = 0; i < 8; ++i){
                d[i] = _mm_load_si128((__m128i *)&(iT[i*32+k*256+l*8]));
            }

            for(int i = 0; i < line; ++i){
                for(int j = 0; j < 8; ++j){
                    x[j] = _mm_set1_epi16(src[j*src_stride+i+k*8*src_stride]);

                    vhi[j] = _mm_mulhi_epi16(x[j], d[j]);
                    vlo[j] = _mm_mullo_epi16(x[j], d[j]);
                    vl[j+8*k] = _mm_unpacklo_epi16(vlo[j], vhi[j]);
                    vh[j+8*k] = _mm_unpackhi_epi16(vlo[j], vhi[j]);
                }
                for(int j = 0; j < 4; ++j){
                    result[i][2*l] = _mm_add_epi32(result[i][2*l], _mm_add_epi32(vl[2*j+8*k],vl[2*j+1+8*k]));
                    result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], _mm_add_epi32(vh[2*j+8*k],vh[2*j+1+8*k]));
                }
            }
        }
    }

    for(int l = 0; l < 4; ++l){
        for(int i = 0; i < line; ++i){
            result[i][2*l] = _mm_add_epi32(result[i][2*l], _mm_set1_epi32(1<<(shift-1)));
            result[i][2*l+1] = _mm_add_epi32(result[i][2*l+1], _mm_set1_epi32(1<<(shift-1)));
            result[i][2*l] = _mm_srai_epi32(result[i][2*l], shift);
            result[i][2*l+1] = _mm_srai_epi32(result[i][2*l+1], shift);
        }
    }

    for(int l = 0; l < 4; ++l){
        for(int i = 0; i < line; i++){
            result[i][l] = _mm_packs_epi32(result[i][2*l], result[i][2*l+1]);
            _mm_store_si128((__m128i *)&(dst[i*32+l*8]), result[i][l]);
        }
    }
}
