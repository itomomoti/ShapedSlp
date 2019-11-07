#include <stdio.h>
#include <stdlib.h>
#include "onlinebp.h"

//#define DEBUG

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) popcount(x)
#endif

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef max
 #define max(x,y) ((x)>(y)?(x):(y))
#endif
#define mymalloc(p,n,t) {p = (t *)malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory at line %d\n",__LINE__); exit(1);};}
#define myrealloc(p,n,t) {p = (t *)realloc((p),(n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory at line %d\n",__LINE__); exit(1);};}

//#define logD 6
//#define logD 5
//#define logD 1
//#define D (1<<logD) // bitvec_t のビット数
#define D (sizeof(bitvec_t)*8)
//#if D == 64
 #define logD 6
//#endif
//#if D == 32
// #define logD 5
//#endif
//#if D == 16
// #define logD 4
//#endif


#define ETW 8  // width of excess lookup table
//#define ETW min(8,D)  // width of excess lookup table
static int bwdtbl[(ETW+1)*(1<<ETW)];
static int fwdtbl[(ETW+1)*(1<<ETW)];
static unsigned int selecttbl[8*256];

#define logMB 15
//#define logMB 4
#define MB (1<<logMB) // 中ブロックのサイズ (ビット数)


static int blog(i64 x)
{
i64 l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

////////////////////////////////////////////////////////////////
//        bit列に対する操作
////////////////////////////////////////////////////////////////

static unsigned int popcount(bitvec_t x)
{
  bitvec_t r;
//  __uint128_t rr;
  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;
#if 1
  r = (r>>8) + r;
  r = (r>>16) + r;
#if logD <= 5
  r = r & 63;
#else
  r = ((r>>32) + r) & 127;
#endif
#else
  r *= 0x0101010101010101;
  r >>= 64-8;
#endif
//  printf("r1 %016lx\n",r);
//  rr = r;
//  rr *= 0x0101010101010101;
//  r = rr >> 56;
//  printf("r2 %016lx\n",r);
//  r &= 0xff;
  return r;
}

static int setbit(bitvec_t *B, i64 i,int x)
{
  i64 j,l;

  j = i / D;
  l = i & (D-1);
  if (x==0) B[j] &= (~(1L<<(D-1-l)));
  else if (x==1) B[j] |= (1L<<(D-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

static void setbits(bitvec_t *B, i64 i, int d, u64 x)
{
  u64 y,m;
  int d2;
#if 0
  int j;

  for (j=0; j<d; j++) {
    setbit(B,i+j,(x>>(d-j-1))&1);
  }
#else
  B += (i>>logD);
  i &= (D-1);

  while (i+d > D) {
    d2 = D-i; // x の上位 d2 ビットを格納
    y = x >> (d-d2);
    m = (1<<d2)-1;
    *B = (*B & (~m)) | y;
    B++;  i=0;
    d -= d2;
    x &= (1<<d)-1; // x の上位ビットを消去
  }
  m = (1<<d)-1;
  y = x << (D-i-d);
  m <<= (D-i-d);
  *B = (*B & (~m)) | y;
#endif
}

static int getbit(bitvec_t *B, i64 i)
{
  i64 j,l;

  j = i >> logD;
  l = i & (D-1);
  return (B[j] >> (D-1-l)) & 1;
}

static bitvec_t getbits(bitvec_t *B, i64 i, int d)
{
  bitvec_t x,z;

  if (d == 0) return 0;
  B += (i >>logD);
  i &= (D-1);
  if (i+d <= D) {
    x = B[0];
    x <<= i;
    x >>= (D-d);  // D==64, d==0 だと動かない
  } else {
    x = B[0] << i;
    x >>= D-d;
    z = B[1] >> (D-(i+d-D));
    x += z;
  }
  return x;
}


////////////////////////////////////////////////////////////////
// postorder heapのインデックス計算
////////////////////////////////////////////////////////////////


static i64 heap_leafpos(int height, i64 b)
{
  i64 x,z;
  int i;
  x = (1 << height) - 1; // rootのindex
  z = 1 << (height-2);
  for (i = height-2; i >= 0; i--) {
    if (b < z) { // 左の部分木
      x -= z*2;
    } else { // 右の部分木
      x -= 1;
      b -= z;
    }
    z >>= 1;
  }
  return x;
}

static i64 heap_parent(i64 b, i64 x, int height)
{
  if (b & (1<<height)) {
    x += 1;
  } else {
    x += 1<<(height+1);
  }
  return x;
}

static i64 heap_rparent(i64 x, int height)
{
  x += 1;
  return x;
}

static i64 heap_lparent(i64 x, int height)
{
  x += 1<<(height+1);
  return x;
}

static i64 heap_lchild(i64 x, int height)
{
  x -= 1<<height;
  return x;
}

static i64 heap_rchild(i64 x, int height)
{
  x -= 1;
  return x;
}

static int heap_isrchild(i64 b, int h)
{
  return (b & (1<<h)) ? 1 : 0;
}

////////////////////////////////////////////////////////////////
// 表を作る
////////////////////////////////////////////////////////////////

void onlinebp_maketbl(void)
{
  i64 i,x,r;
  bitvec_t buf[1];

  for (x = 0; x < (1<<ETW); x++) {
    setbits(buf,0,ETW,x);
    for (r=-ETW; r<=0; r++) bwdtbl[((r+ETW)<<ETW)+x] = ETW;
    for (r=-ETW; r<=0; r++) fwdtbl[((r+ETW)<<ETW)+x] = ETW;

    r = 0;
    for (i=0; i<ETW; i++) {
      if (getbit(buf,i)==OP) {
        r++;
      } else {
        r--;
      }
      if (r <= 0 && fwdtbl[((r+ETW)<<ETW)+x] == ETW) fwdtbl[((r+ETW)<<ETW)+x] = i;
    }

    r = 0;
    for (i=ETW-1; i>=0; i--) {
      if (getbit(buf,i)==OP) {
        r--;
      } else {
        r++;
      }
      if (r <= 0 && bwdtbl[((r+ETW)<<ETW)+x] == ETW) bwdtbl[((r+ETW)<<ETW)+x] = ETW-1-i;
    }
  }

}

static void make_selecttbl(void)
{
  i64 i,x,r;
  bitvec_t buf[1];

  for (x = 0; x < 256; x++) {
    setbits(buf,0,8,x);
    for (r=0; r<8; r++) selecttbl[(r<<8)+x] = -1;
    r = 0;
    for (i=0; i<8; i++) {
      if (getbit(buf,i)) {
        selecttbl[(r<<8)+x] = i;
        r++;
      }
    }
  }

}

////////////////////////////////////////////////////////////////
// 内部変数のアクセス
////////////////////////////////////////////////////////////////


i64 onlinebp_size(onlinebp *bp)
{
  return bp->size;
}

i64 onlinebp_length(onlinebp *bp)
{
  return bp->n;
}

static bitvec_t *getbp(onlinebp *bp, i64 i)
{
  bitvec_t *p;
#ifdef DEBUG
  if (i % D != 0) {
    printf("getbp: i = %ld\n", i);
  }
#endif
  p = &bp->BP[i / D];
  return p;
}



static i64 getdepth(onlinebp *bp, i64 i)
{
  i64 x;
#ifdef DEBUG
  if (i % bp->L != 0) {
    printf("getdepth: i = %ld\n", i);
  }
#endif
  x = bp->depth_H[i>>logMB] + bp->depth_L[i/bp->L];
  return x;
}

static void setdepth(onlinebp *bp, i64 i, i64 x)
{
  i64 y;
#ifdef DEBUG
  if (i % bp->L != 0) {
    printf("setdepth: i = %ld\n", i);
  }
#endif
  if ((i & (MB-1)) == 0) {
    bp->depth_H[i>>logMB] = x;
  }
  y = x - bp->depth_H[i>>logMB];
#ifdef DEBUG
    if (y < -(1<<15) || y > (1<<15)) {
      printf("setdepth: i = %ld x = %ld, MB %ld\n", i, x, bp->depth_H[i/MB]);
    }
#endif
  bp->depth_L[i/bp->L] = (short)y;
}

////////////////////////////////////////////////////////////////
//        ワードサイズ D
//        小ブロックサイズ L
//        中ブロックサイズ MB <= (1<<15)
//
// RM木のノードのインデックスを x とすると
// x の下位 MB ビットは中ブロック内のインデックス [1..(MB/L)*2-1]
// 中ブロックのヒープの高さ heap_height_L = log(2MB/L)
// x の上位ビットは中ブロックのインデックス [0..]
// x = 0 のとき，上位ビットは全体のヒープのインデックス
//
// 中ブロック用のRM木は RM_L[0..n_MB-1] に格納 (ポインタの配列)
// RM_Lの長さは N_MB (>= n_MB)
////////////////////////////////////////////////////////////////


onlinebp *onlinebp_new(int L, int alpha)
{
  onlinebp *bp;

  onlinebp_maketbl(); // 毎回やる必要はないけど
  make_selecttbl(); // 毎回やる必要はないけど

  mymalloc(bp,1,onlinebp);
  bp->size = sizeof(*bp);
  bp->n = 0;
  bp->d = 0; // 左端(位置 = -1)の深さは 0
  L = L/D*D; // L は D の倍数に
  if (L == 0) {
    printf("onlinebp_new: L=%d must be >= %d\n", L, D);
    exit(1);
  }
  bp->L =  L;
  bp->N = 0;
  bp->BP = NULL;
  bp->depth_H = NULL;
  bp->depth_L = NULL;
  if (alpha < 2) {
    printf("onlinebp_new: alpha=%d must be >= 2\n", alpha);
    exit(1);
  }
  bp->alpha = alpha;

  bp->heap_height_H = 0;
  bp->heap_height_L = blog(2*MB/L);
  bp->RM_H = NULL;
  bp->RM_L = NULL;
  bp->block_min = 100; // 初期値は参照されないはずなので何でもいい

  bp->n_MB = 0;
  bp->N_MB = 0;
  bp->heap_size_H = 0;

  mymalloc(bp->match,1,i64);
  bp->match[0] = -1;
  bp->maxdepth = 0;
  bp->lastOP = -1; // ?
  bp->size += 1*sizeof(bp->match[0]);

  return bp;
}

void onlinebp_free(onlinebp *bp)
{
  free(bp->BP);
  free(bp->depth_H);
  free(bp->depth_L);
  free(bp->RM_H);
  free(bp->RM_L);
  free(bp->match);
  free(bp);
}


void onlinebp_reserve(onlinebp * bp, u64 capaBitLen)
{
  i64 n, N;
  int L;

  n = bp->n;
  N = bp->N;
  L = bp->L;

  i64 i,s,t, N2;

  // BP の領域を拡張
  N2 = N;
  s = N / D;
  N = capaBitLen;
  t = N / D;
  myrealloc(bp->BP, N / D, bitvec_t);
  bp->size += (t-s)*sizeof(bp->BP[0]); // for BP
  for (i = s; i < t; i++) bp->BP[i] = 0; // 増やした部分を 0 にしておく
  myrealloc(bp->depth_H, (N+MB-1)/L, i64);
  myrealloc(bp->depth_L, (N+L-1)/L, short);
  bp->size -= (N2+MB-1)/MB*sizeof(bp->depth_H[0]); // for depth_H
  bp->size += (N+MB-1)/MB*sizeof(bp->depth_H[0]); // for depth_H
  bp->size -= (N2+L-1)/L*sizeof(bp->depth_L[0]); // for depth_L
  bp->size += (N+L-1)/L*sizeof(bp->depth_L[0]); // for depth_L
  bp->N = N;

  i64 N2_MB;
  N2_MB = bp->N_MB;
  bp->N_MB = (capaBitLen + MB) / MB;
  myrealloc(bp->RM_L, bp->N_MB, short *);
  bp->size -= N2_MB*sizeof(bp->RM_L[0]);
  bp->size += bp->N_MB*sizeof(bp->RM_L[0]);
  for (i = N2_MB; i < bp->N_MB; i++) bp->RM_L[i] = NULL;
}


void onlinebp_push(onlinebp *bp, int c)
{
  i64 n, N, d;
  int L;

  n = bp->n;
  N = bp->N;
  L = bp->L;
  d = bp->d;

#ifdef DEBUG
  if (c != OP && c != CP) {
    printf("onlinebp_push: c = %d\n", c);
    exit(1);
  }
#endif

  bp->match[d] = n-1;

  if (n >= N) {
    i64 i,s,t, N2;

    // BP の領域を拡張
    N2 = N;
    s = N / D;
    N = max(N * bp->alpha / (bp->alpha-1),L);
    t = N / D;
    myrealloc(bp->BP, N / D, bitvec_t);
    bp->size += (t-s)*sizeof(bp->BP[0]); // for BP
    for (i = s; i < t; i++) bp->BP[i] = 0; // 増やした部分を 0 にしておく
    myrealloc(bp->depth_H, (N+MB-1)/L, i64);
    myrealloc(bp->depth_L, (N+L-1)/L, short);
    bp->size -= (N2+MB-1)/MB*sizeof(bp->depth_H[0]); // for depth_H
    bp->size += (N+MB-1)/MB*sizeof(bp->depth_H[0]); // for depth_H
    bp->size -= (N2+L-1)/L*sizeof(bp->depth_L[0]); // for depth_L
    bp->size += (N+L-1)/L*sizeof(bp->depth_L[0]); // for depth_L
  }

  if (n % L == 0) { // ブロックの最初
    bp->block_min = n+100;
    setdepth(bp, n, d);
  }
  setbit(bp->BP, n, c);
  d += c*2-1;

  if ((n & (MB-1)) == 0) { // 中ブロックの最初
    int j,m;
    short *p;
    bp->n_MB++;
    if (bp->n_MB > bp->N_MB) { // 領域が足りなくなった
      printf("update\n");
      i64 i, N2_MB;
      N2_MB = bp->N_MB;
      bp->N_MB = (bp->N_MB+1) * bp->alpha / (bp->alpha-1);
      myrealloc(bp->RM_L, bp->N_MB, short *);
      bp->size -= N2_MB*sizeof(bp->RM_L[0]);
      bp->size += bp->N_MB*sizeof(bp->RM_L[0]);
      for (i = N2_MB; i < bp->N_MB; i++) bp->RM_L[i] = NULL;
    }
    m = 2*MB/L; // 中ブロックのRM木の領域のサイズ
    mymalloc(bp->RM_L[bp->n_MB-1], m, short); // 中ブロックのRM木の領域
    p = bp->RM_L[bp->n_MB-1];
    for (j=0; j<m; j++) p[j] = MB-1; // 最大値を入れておく
  }

  if (bp->n_MB > (bp->heap_size_H+1)/2) { // 全体のRM木の領域が足りなくなった
    i64 i,s,t;
    s = bp->heap_size_H+1;
    bp->heap_size_H = bp->heap_size_H*2+1;
    bp->heap_height_H++; // heap_size_H = (1 << heap_height_H)-1;
    t = bp->heap_size_H+1; // +1 はindexが1から始まるから
    myrealloc(bp->RM_H, t, i64);
    bp->size += (t-s)*sizeof(bp->RM_H[0]);

    for (i = s; i < t; i++) bp->RM_H[i] = n+100; // 増やした部分を最大値にしておく
    if (bp->heap_size_H > 1) {
      bp->RM_H[bp->heap_size_H] = bp->RM_H[bp->heap_size_H/2]; // これまでの最小値を親へコピー
    }
  }


  if (d < bp->block_min) {
    int j;
    i64 x,b;
    i64 bm, bs;
    short *p;
    i64 d2, dm;
    b = n / L; // 小ブロックの番号
    bm = b / (MB/L); // 中ブロックの番号 
    bs = b % (MB/L); // 中ブロック内の番号 

    p = bp->RM_L[bm];
    x = heap_leafpos(bp->heap_height_L, bs); // 中ブロック内の葉の場所
    dm = getdepth(bp,bm*(MB/L)*L);
    d2 = d - dm; // 相対の深さ
#ifdef DEBUG
    if (d2 >= MB || d2 <= -MB) {
      printf("??? d2 = %ld\n", d2);
    }
#endif
    for (j=0; j<bp->heap_height_L; j++) {
      if (d >= dm + p[x]) break; // 最小値が更新されなかった
      p[x] = (short)d2;
      x = heap_parent(bs, x, j);
    }
    if (j == bp->heap_height_L) { // 根が更新された
      x = heap_leafpos(bp->heap_height_H, bm); // 葉の場所
      for (j=0; j<bp->heap_height_H; j++) {
        if (d >= bp->RM_H[x]) break; // 最小値が更新されなかった
        bp->RM_H[x] = d;
        x = heap_parent(bm, x, j);
      }
    }
  }



  n++; // 括弧列は [-1,n-1], BP[-1] はダミーで深さが 0 とする


  if (d < 0) {
    printf("onlinebp_push: depth = %ld\n", d);
    exit(1);
  }

//////////////////////////////////  
// 最後の括弧用のデータ
  if (d > bp->maxdepth) {
    bp->maxdepth = d;
    myrealloc(bp->match, bp->maxdepth+1, i64);
    bp->size += sizeof(bp->match[0]);
  }
  if (c == OP) bp->lastOP = n-1;
//////////////////////////////////  

  bp->n = n;
  bp->N = N;
  bp->L = L;
  bp->d = d;
}

// posの1つ左から探す
// 答えの範囲は [-1,pos-1] または NOT_FOOUND (-2)
#if 0
static i64 bwd_block_naive(bitvec_t *BP, i64 pos, i64 *k_)
{
  i64 i;
  int c;
  int d;
  int k;
  k = *k_;
  d = 0; // pos の相対深さ
  for (i=pos; i>=0; i--) {
    c = getbit(BP, i);
    d -= 2*c-1; // 位置 i-1 の深さは BP[i] から決まる
    if (d == k) return i-1; // 見つかった場所を返す
  }
  *k_ = d;
  return NOT_FOUND;
}
#endif

static i64 bwd_block(bitvec_t *BP, i64 pos, i64 *k_)
{
  i64 i;
  int d;
  int k,j,r;
  bitvec_t *p, x, w;

  k = *k_;

  p = &BP[pos / D];
  i = pos;
  d = k;
  while (i >= 0) {
    x = *p--;
    j = (i & (D-1))+1; // 最後のワード内の有効なビット数
    x >>= D-j; // 無効なビットを消す
    while (j>0) {
      w = x & ((1<<ETW)-1);
      if (k >= -ETW && k <= 0) { // 答えがある可能性がある
        r = bwdtbl[((k+ETW)<<ETW)+w];
        if (r<ETW && r<j) {
          if (i-r < 0) return NOT_FOUND;
          return i-r-1;
        }
      }
      r = min(j,ETW);
      k += 2*POPCOUNT(w)-r;
      x >>= r;
      i -= r;
      j -= r;
    }
  }
  *k_ = d-k;
  return NOT_FOUND;
}

static i64 bwd_MB(onlinebp *bp, i64 b, i64 depth)
{
  i64 b2, x;
  int h;
  i64 d;
  short *p;

  i64 bm, bs;

  bm = b / (MB/bp->L); // 中ブロックの番号 
  bs = b % (MB/bp->L); // 中ブロック内の番号 

  p = bp->RM_L[bm];
  x = heap_leafpos(bp->heap_height_L, bs); // 中ブロック内の葉の場所
  d = getdepth(bp,bm*(MB/bp->L)*bp->L);

  b2 = bs;

  for (h = 0; h < bp->heap_height_L-1; h++) {
    if (heap_isrchild(bs,h)) { // 右の子
      x = heap_rparent(x,h);
      x = heap_lchild(x, h+1); // 左の子へ
      b2 -= 1<<h;
      if (d + p[x] <= depth) break;
    }
    x = heap_lparent(x,h); // 親へ
  }
  if (h == bp->heap_height_L-1) return NOT_FOUND;

  while (h > 0) {
    if (d + p[heap_rchild(x,h)] <= depth) { // 右の子で見つかった
      x = heap_rchild(x,h); // 右の子へ
      b2 += 1<<(h-1);
    } else {
      x = heap_lchild(x,h); // 左の子へ
    }
    h--;
  }
  return bm * (MB/bp->L) + b2;
}

i64 onlinebp_bwd(onlinebp *bp, i64 i, i64 d)
{
  i64 pos,s,t;
  bitvec_t *p;
  int L;
  i64 depth;
  i64 b,b2,x,d2;
  int h;
  i64 bm, bs;
  short *q;

  if (d > 0) {
    printf("onlinebp_bwd: d = %ld\n", d);
    exit(1);
  }

  if (i < 0 || i > bp->n) return NOT_FOUND;
  L = bp->L;

/////////////////////////////////////
// 小ブロックを探索

  s = i / L * L; // i を含む小ブロックの先頭
  t = i - s; // 小ブロック内の位置
  p = getbp(bp,s);
  d2 = d;
  pos = bwd_block(p, t, &d2);
  if (pos != NOT_FOUND) {
    return s+pos;
  }
  depth = getdepth(bp,s) - d2; // i の深さ
  depth += d; // 探している深さ
  if (depth == 0) return -1; // 0 は一番左 (-1) のみ

  b = s / L; // 小ブロック b には無かった

/////////////////////////////////////
// 中ブロックを探索

  b2 = bwd_MB(bp, b, depth); // depthを含む小ブロックを探す
  if (b2 != NOT_FOUND) goto search_sb;


/////////////////////////////////////
// 全体を探索
  // 次にどこから探すか計算
  b2 = b;
  // b2の下位ビットを 0 にする
  b2 >>= bp->heap_height_L-1;
  b2 <<= bp->heap_height_L-1;

  bm = b / (MB/bp->L); // 中ブロックの番号 
  x = heap_leafpos(bp->heap_height_H, bm); // 葉の場所
  for (h = 0; h < bp->heap_height_H-1; h++) {
    if (heap_isrchild(bm,h)) { // 右の子
      x = heap_rparent(x,h);
      x = heap_lchild(x, h+1); // 左の子へ
      b2 -= 1<<(h+bp->heap_height_L-1);
      if (bp->RM_H[x] <= depth) break;
    }
    x = heap_lparent(x,h); // 親へ
  }
  if (h == bp->heap_height_H-1) return NOT_FOUND;

  while (h > 0) {
    if (bp->RM_H[heap_rchild(x,h)] <= depth) { // 右の子で見つかった
      x = heap_rchild(x,h); // 右の子へ
      b2 += 1<<(h-1+bp->heap_height_L-1);
    } else {
      x = heap_lchild(x,h); // 左の子へ
    }
    h--;
  }

  bm = b2 / (MB/bp->L); // 中ブロックの番号 
  bs = b2 % (MB/bp->L); // 中ブロック内の葉の番号 
  h = bp->heap_height_L-1;

/////////////////////////////////////
// 中ブロックを探索


  q = bp->RM_L[bm];
  x = (1 << bp->heap_height_L)-1; // x は中ブロックの根

  d = getdepth(bp,bm*(MB/bp->L)*bp->L);

  while (h > 0) {
    if (d + q[heap_rchild(x,h)] <= depth) { // 右の子で見つかった
      x = heap_rchild(x,h); // 右の子へ
      b2 += 1<<(h-1);
    } else {
      x = heap_lchild(x,h); // 左の子へ
    }
    h--;
  }


/////////////////////////////////////
// 小ブロック b2 を探索

search_sb:

  s = b2*L;
  t = L-1;
  p = getbp(bp, s);
  d = getdepth(bp,(b2+1)*L); // ブロック b2 の最後の深さ
  if (d == depth) return s+t;
  d2 = depth - d;
  pos = bwd_block(p, t, &d2);
  if (pos != NOT_FOUND) {
    return s+pos;
  }
  printf("bwd: ???\n");
  return NOT_FOUND;
}

// posの1つ左から探す
// 答えの範囲は [-1,pos-1] または NOT_FOOUND (-2)
#if 0
static i64 naive_bwd(onlinebp *bp, i64 pos, i64 k)
{
  i64 i;
  int c;
  i64 d;
  bitvec_t *BP = bp->BP;
  d = 0; // pos の相対深さ
  for (i=pos; i>=0; i--) {
    c = getbit(BP, i);
    d -= 2*c-1;
    if (d == k) return i-1; // 見つかった場所を返す
  }
  return NOT_FOUND;
}
#endif

// posの1つ右から探す
// 答えの範囲は [pos+1,n-1] または NOT_FOOUND (-2)
#if 0
static i64 naive_fwd(onlinebp *bp, i64 pos, i64 k)
{
  i64 i;
  int c;
  i64 d;
  bitvec_t *BP = bp->BP;
  pos++;
  d = 0; // pos の相対深さ
  for (i=pos; i<bp->n; i++) {
    c = getbit(BP, i);
    d += 2*c-1;
    if (d == k) return i; // 見つかった場所を返す
  }
  return NOT_FOUND;
}
#endif

// posから探す
// 答えの範囲は [pos,t] または NOT_FOOUND (-2)
#if 0
static i64 fwd_block_naive(bitvec_t *BP, i64 pos, i64 t, i64 k)
{
  i64 i;
  int c;
  int d;

  d = 0; // pos の相対深さ
  for (i=pos; i<=t; i++) {
    c = getbit(BP, i);
    d += 2*c-1;
    if (d == k) return i; // 見つかった場所を返す
  }
  return NOT_FOUND;
}
#endif

static i64 fwd_block(bitvec_t *BP, i64 pos, i64 t, i64 k)
{
  i64 i;
  int d,j;
  bitvec_t *p, x;

  p = &BP[pos / D];
  i = pos;
  d = k;

  while (i<=t) {
    x = *p++;
    j = i & (D-1);
    x <<= j;
    j = D-j;
    while (j>0) {
      i64 w,r;
      w = (x >> (D-ETW)) & ((1<<ETW)-1);
      if (k >= -ETW && k <= 0) {
        r = fwdtbl[((k+ETW)<<ETW)+w];
        if (r<ETW && r<j) {
          if (i+r > t) return NOT_FOUND;
          return i+r;
        }
      }
      r = min(j,ETW);
      k -= 2*POPCOUNT(w)-r;
      x <<= r;
      i += r;
      j -= r;
    }
  }

  return NOT_FOUND;
}

static i64 fwd_MB(onlinebp *bp, i64 b, i64 depth)
{
  i64 b2, x;
  int h;
  i64 d;
  short *p;

  i64 bm, bs;

  bm = b / (MB/bp->L); // 中ブロックの番号 
  bs = b % (MB/bp->L); // 中ブロック内の番号 

  p = bp->RM_L[bm];
  x = heap_leafpos(bp->heap_height_L, bs); // 中ブロック内の葉の場所
  d = getdepth(bp,bm*(MB/bp->L)*bp->L);
  
  b2 = bs;

  for (h = 0; h < bp->heap_height_L-1; h++) {
    if (!heap_isrchild(bs,h)) { // 左の子
      x = heap_lparent(x,h); // 親へ
      x = heap_rchild(x,h+1); // 右の子へ
      if (d + p[x] <= depth) {
        b2 += 1<<h;
        break;
      }
    } else { // 右の子
      b2 -= 1<<h;
    }
    x = heap_rparent(x,h); // 親へ
  }
  if (h == bp->heap_height_L-1) return NOT_FOUND;

  while (h > 0) {
    x = heap_lchild(x,h); // 左の子へ
    if (d + p[x] > depth) { // 見つからない
      x = heap_lparent(x,h-1); // 親へ
      x = heap_rchild(x,h); // 右の子へ
      b2 += 1<<(h-1);
    }
    h--;
  }

  return bm * (MB/bp->L) + b2;
}


i64 onlinebp_fwd(onlinebp *bp, i64 i, i64 d)
{
  i64 pos,s,t,u;
  bitvec_t *p;
  int L;
  i64 depth;
  i64 b,b2,x,d2;
  int h;
  i64 bm, bs;
  short *q;

  if (d > 0) {
    printf("onlinebp_fwd: d = %ld\n", d);
    exit(1);
  }

  if (i < 0 || i > bp->n) return NOT_FOUND;
  L = bp->L;

  i++; // 1つ右から探す

/////////////////////////////////////
// 小ブロックを探索

  s = i / L * L; // i を含むブロックの先頭
  t = i - s; // ブロック内の位置
  u = min(s + L-1, bp->n-1); // ブロックの最後
  p = getbp(bp, s);
  d2 = d;
  pos = fwd_block(p, t, u, d2);
  if (pos != NOT_FOUND) {
    return s+pos;
  }
  depth = onlinebp_depth(bp, i-1);
  depth += d; // 探している深さ
//  if (depth == 0) return NOT_FOUND;

  b = s / L; // 小ブロック b には無かった

/////////////////////////////////////
// 中ブロックを探索

  b2 = fwd_MB(bp, b, depth); // depthを含む小ブロックを探す
  if (b2 != NOT_FOUND) goto search_sb;


/////////////////////////////////////
// 全体を探索

  // 次にどこから探すか計算
  b2 = b;
  // b2の下位ビットを 0 にする
  b2 >>= bp->heap_height_L-1;
  b2 <<= bp->heap_height_L-1;
  // heapの高さ bp->heap_height_L-1 から探す

  bm = b / (MB/bp->L); // 中ブロックの番号 
  x = heap_leafpos(bp->heap_height_H, bm); // 葉の場所
  for (h = 0; h < bp->heap_height_H-1; h++) {
    if (!heap_isrchild(bm,h)) { // 左の子
      x = heap_lparent(x,h); // 親へ
      x = heap_rchild(x,h+1); // 右の子へ
      if (bp->RM_H[x] <= depth) {
        b2 += 1<<(h+bp->heap_height_L-1);
        break;
      }
    } else { // 右の子
      b2 -= 1<<(h+bp->heap_height_L-1);
    }
    x = heap_rparent(x,h); // 親へ
  }
  if (h > bp->heap_height_H-1) return NOT_FOUND;

  while (h > 0) {
    x = heap_lchild(x,h); // 左の子へ
    if (bp->RM_H[x] > depth) { // 見つからない
      x = heap_lparent(x,h-1); // 親へ
      x = heap_rchild(x,h); // 右の子へ
      b2 += 1<<(h-1+bp->heap_height_L-1);
    }
    h--;
  }


  bm = b2 / (MB/bp->L); // 中ブロックの番号 
  bs = b2 % (MB/bp->L); // 中ブロック内の葉の番号 
  h = bp->heap_height_L-1;


/////////////////////////////////////
// 中ブロックを探索

  q = bp->RM_L[bm];
  x = (1 << bp->heap_height_L)-1; // x は中ブロックの根

  d = getdepth(bp,bm*(MB/bp->L)*bp->L);

  while (h > 0) {
    x = heap_lchild(x,h); // 左の子へ
    if (d + q[x] > depth) { // 見つからない
      x = heap_lparent(x,h-1); // 親へ
      x = heap_rchild(x,h); // 右の子へ
      b2 += 1<<(h-1);
    }
    h--;
  }


/////////////////////////////////////
// 小ブロック b2 を探索

search_sb:
  s = b2*L;
  t = L-1;
  p = getbp(bp, s);
  d = getdepth(bp, s); // ブロック b2 の直前の深さ
//  if (d == depth) return s+t; // これ要るの?
  d2 = depth - d;
  pos = fwd_block(p, 0, t, d2);
  if (pos != NOT_FOUND) {
    return s+pos;
  }
  printf("fwd: ???\n");
  return NOT_FOUND;
}

////////////////////////////////////////////////////////////////
//        括弧列に対する操作
////////////////////////////////////////////////////////////////

int onlinebp_alpha(onlinebp *bp){
  return bp->alpha;
}
int onlinebp_L(onlinebp *bp){
  return bp->L;
}


i64 onlinebp_depth(onlinebp *bp, i64 i)
{
  i64 s;
  bitvec_t *p, x;
  int L,q,c;
  i64 d;

  if (i < -1 || i >= bp->n) {
    // printf("depth: i = %ld\n", i);
    return -1;
    //exit(1);
  }
  if (i == -1) return 0;
  L = bp->L;

  s = i / L * L; // i を含むブロックの先頭
  q = i - s; // ブロック内の位置
  p = getbp(bp, s);
  d = getdepth(bp, s);
#if 0
  for (j = 0; j <= q; j++) {
    int c;
    c = getbit(p, j);
    d += 2*c-1;
  }
#else
  while (q >= D) {
    x = *p++;
    c = POPCOUNT(x);
    d += c - (D-c);
    q -= D;
  }
  x = *p++;
  x >>= D-1-q;
  c = POPCOUNT(x);
  d += c - (q+1-c);
#endif
  return d;
}

i64 onlinebp_rank(onlinebp *bp, i64 i, int c)
{
  i64 d, r;
  d = onlinebp_depth(bp, i);
  r = (d+i+1)/2;
  if (c == CP) r = (i+1) - r;
  return r;
}

#if 0
static i64 onlinebp_select_naive(onlinebp *bp, i64 i, int c)
{
  i64 l,r,m;

  l = 0; r = bp->n-1;
  while (l <= r) {
    m = (l+r)/2;
    if (i <= onlinebp_rank(bp, m, c)) r = m-1;  else l = m+1;
  }
  return l;
}
#endif


static i64 select_by_rank(onlinebp *bp, i64 ith, int c)
{
  bitvec_t x, *buf;
  
  i64 j,k;
  i64 ofs, ofs2;
  i64 r, r2;
  int rr;
  i64 d;
  i64 z,L;

  i64 ll, rl, ml, pl; // 大ブロックの2分探索用
  i64 p; // 答え
  i64 ii;
  bitvec_t *q;
  i64 p0;
  i64 l0;

  ii = ith;

  ll = 0;  rl = (bp->n-1) >> logMB;
  pl = ll;  r2 = 0;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
//    r = getuint(da->rl,ml,da->k);
//    x = bp->depth_H[i/MB] + bp->depth_L[i/bp->L];
    z = bp->depth_H[ml];
//  r = (d+i+1)/2;
    r = (z+(ml<<logMB))/2;
    if (c == 0) {r = (ml<<logMB) - r;} // select0 の場合
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 大ブロック内のランク

  L = bp->L;

  ll = (pl<<logMB)/L;
  rl = min(((pl+1)<<logMB)/L-1, (bp->n-1)/L);

  pl = ll;  r2 = 0;
  l0 = ll;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
//    r = da->rs[ml];
    z = bp->depth_L[ml];
    r = (z+(ml-l0)*L)/2;
    if (c == 0) {r = ((ml-l0)*L) - r;} // select0 の場合
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 小ブロック内のランク

  p = pl * L; // ith を含む小ブロックの最初のビット位置
  rr = 0;

  p0 = p;
//  q = &(da->buf[p>>logD]);
  q = getbp(bp, p);

  while (1) {
    x = *q;
    if (c == 0) x = ~x;
    rr = POPCOUNT(x);
    if (rr >= ith) break;
    ith -= rr;
    p += D;
    q++;
  }
      
  x = *q;
  if (c == 0) x = ~x;
  while (1) {
    rr = POPCOUNT(x >> (D-8));
    if (rr >= ith) break;
    ith -= rr;
    p += 8;
    x <<= 8;
  }
  p += selecttbl[((ith-1)<<8)+(x>>(D-8))];

  return p;
  
}


i64 onlinebp_select(onlinebp *bp, i64 i, int c)
{
i64 s1, s2;
//  s1 = onlinebp_select_naive(bp, i, c);
  s2 = select_by_rank(bp, i, c);
//  if (s1 != s2) {
//    printf("onlinebp_select: i = %ld c = %d s1 = %ld s2 = %ld\n", 
//      i, c, s1, select_by_rank(bp, i, c));
//  }
  return s2;
}


int onlinebp_get(onlinebp *bp, i64 i)
{
#ifdef DEBUG
  if (i < 0 || i >= bp->n) {
    printf("onlinebp_get: i = %ld n = %ld\n", i, bp->n);
    return -1;
  }
#endif
  return getbit(bp->BP, i);
}

i64 onlinebp_findopen(onlinebp *bp, i64 i)
{
#ifdef DEBUG
  if (onlinebp_get(bp, i) != CP) {
    printf("onlinebp_findopen: BP[i] = %d\n", i, onlinebp_get(bp, i));
    return NOT_FOUND;
  }
#endif
  if (i == bp->n-1) { // 最後の括弧
    return bp->match[bp->d]+1;
  }
  return onlinebp_bwd(bp, i, 0)+1;
}

i64 onlinebp_findclose(onlinebp *bp, i64 i)
{
#ifdef DEBUG
  if (onlinebp_get(bp, i) != OP) {
    printf("onlinebp_findclose: BP[i] = %d\n", i, onlinebp_get(bp, i));
    return NOT_FOUND;
  }
#endif
  return onlinebp_fwd(bp, i, -1);
}

i64 onlinebp_enclose(onlinebp *bp, i64 i)
{
#ifdef DEBUG
  if (onlinebp_get(bp, i) != OP) {
    printf("onlinebp_enclose: BP[i] = %d\n", i, onlinebp_get(bp, i));
    return NOT_FOUND;
  }
#endif
  return onlinebp_bwd(bp, i, -2)+1;
}

int convbp(int c)
{
  int d;
  switch (c) {
  case '(':
  case '1':
  case 1:
    d = OP;
    break;
  case ')':
  case '0':
  case 0:
    d = CP;
    break;
  default:
    printf("convbp: c = %d(%c)???\n",c,c);
    exit(1);
    break;
  }
  return d;
}

////////////////////////////////////////////////////////////////
//        postorder full binary treeに対する操作
////////////////////////////////////////////////////////////////

int onlinebp_rightchild(onlinebp *bp, i64 i, i64 *s, i64 *t)
{
#ifdef DEBUG
  if (onlinebp_get(bp, i) != CP) {
    printf("onlinebp_rightchild: BP[i] = %d\n", i, onlinebp_get(bp, i));
    *s = *t = -1;
    return NOT_FOUND;
  }
#endif
  *t = i-1;

  if (i == bp->n-1) { // 最後の括弧
//    *s = bp->match[bp->d]+1;
  } else {
//    *s = onlinebp_bwd(bp, i, 0)+1;
  }
  return 0;
}

int onlinebp_leftchild(onlinebp *bp, i64 i, i64 *s, i64 *t)
{
#ifdef DEBUG
  if (onlinebp_get(bp, i) != CP) {
    printf("onlinebp_leftchild: BP[i] = %d\n", i, onlinebp_get(bp, i));
    *s = *t = -1;
    return NOT_FOUND;
  }
#endif
  if (i == bp->n-1) { // 最後の括弧
    *t = bp->match[bp->d];
//    *s = bp->match[bp->d-1]+1;
  } else {
    *t = onlinebp_bwd(bp, i, 0);
//    *s = onlinebp_bwd(bp, i, -1)+1;
  }
  return 0;
}

i64 onlinebp_leftmost_leaf(onlinebp *bp, i64 i)
{
  i64 x;
#ifdef DEBUG
  if (onlinebp_get(bp, i) != CP) {
    printf("onlinebp_leftmost_leaf: BP[i] = %d\n", i, onlinebp_get(bp, i));
    return NOT_FOUND;
  }
#endif
  if (i == bp->n-1) { // 最後の括弧
    x = bp->match[bp->d-1]+1;
  } else {
    x = onlinebp_bwd(bp, i, -1)+1;
  }
  return x;
}

i64 onlinebp_rightmost_leaf(onlinebp *bp, i64 i)
{
  i64 x;
#ifdef DEBUG
  if (onlinebp_get(bp, i) != CP) {
    printf("onlinebp_rightmost_leaf: BP[i] = %d\n", i, onlinebp_get(bp, i));
    return NOT_FOUND;
  }
#endif
  if (i == bp->n-1) { // 最後の括弧
    x = bp->lastOP;
  } else {
    x = onlinebp_select(bp, onlinebp_rank(bp, i, OP), OP);
  }
  return x;
}

int onlinebp_isrightchild(onlinebp *bp, i64 i)
{
#ifdef DEBUG
  if (i >= bp->n) {
    printf("onlinebp_isrightchild: i = %ld must be < n = %ld\n", i, bp->n);
    return NOT_FOUND;
  }
#endif
  if (i == bp->n-1) {
    return NOT_FIXED;
  }
  if (onlinebp_get(bp, i+1) == CP) {
    return 1; // i is a right child
  } else {
    return 0; // i is a left child
  }
}

int onlinebp_isleftchild(onlinebp *bp, i64 i)
{
#ifdef DEBUG
  if (i >= bp->n) {
    printf("onlinebp_isleftchild: i = %ld must be < n = %ld\n", i, bp->n);
    return NOT_FOUND;
  }
#endif
  if (i == bp->n-1) {
    return NOT_FIXED;
  }
///////////////
  if (onlinebp_get(bp, i+1) == OP) {
    return 1; // i is a left child
  } else {
    return 0; // i is a right child
  }
}

i64 onlinebp_parent(onlinebp *bp, i64 i)
{
  int c;
  c = onlinebp_isrightchild(bp, i);
  if (c < 0) return c;
  if (c == 1) { // right child
    return i+1;
  } else {
    return onlinebp_findclose(bp, i+1);
  }
}


#ifdef ONLINEBP_MAIN
int main(int argc, char *argv[])
{
  onlinebp *bp;
  char *P = "(()((()()))";
  char *q;
  int c;
  i64 i,n;

  n = 1000000;

  bp = onlinebp_new(512,2);
  q = P;
  i = 0;
  while (i < n) {
    i64 s,t,d,j;
//    c = convbp(*q++);
    if (onlinebp_depth(bp, bp->n-1) == 1 || rand() % 2) {
      c = OP;
    } else {
      c = CP;
    }
    onlinebp_push(bp, c);
//    printf("%ld: %c d=%ld\n", i, (c==OP)?'(':')', onlinebp_depth(bp, bp->n-1));
//    printf("pos = %ld depth = %ld,%ld\n", bp->n-1, bp->d, onlinebp_depth(bp, bp->n-1));
    if (c == CP) {
      onlinebp_rightmost_leaf(bp, bp->n-1);
      d = bp->d;
      t = bp->n-1-1; // 右の子の最後
//      s = bp->match[d]+1; // 右の子の最初
      s = onlinebp_bwd(bp, bp->n-1, 0)+1;
//      printf("right tree [%ld,%ld]\n", s,t);
      t = s-1; // 左の子の最後
//      s = bp->match[d-1]+1; // 左の子の最初
      s = onlinebp_bwd(bp, bp->n-1, -1)+1;
//      printf("left  tree [%ld,%ld]\n", s,t);
#ifdef DEBUG
      for (j=0; j>=max(-d,-10); j--) {
        i64 x,y;
        x = onlinebp_bwd(bp, bp->n-1, j);
        if (x != bp->match[d+j]) {
          printf("bwd(%ld,%d) = %ld, %ld\n", bp->n-1, j, onlinebp_bwd(bp, bp->n-1, j), bp->match[d+j]);
        }
        if (j == 0) {
          y = onlinebp_fwd(bp, x, j);
          if (y != bp->n-1) {
            printf("fwd(%ld,%d) = %ld, %ld\n", x, j, onlinebp_fwd(bp, x, j), bp->n-1);
          }
        }
      }
#endif
    }
    i++;
  }
  printf("n = %ld size = %ld bytes (%f bpc)\n", bp->n, onlinebp_size(bp), (double)onlinebp_size(bp)*8/bp->n);
}
#endif
