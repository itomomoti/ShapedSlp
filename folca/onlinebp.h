#ifndef _ONLINEBP_H_
#define _ONLINEBP_H_

// #include <fstream>

#define OP 1 // (
#define CP 0 // )

typedef long i64;
typedef unsigned long u64;
typedef u64 bitvec_t;
//typedef unsigned int bitvec_t;

#define NOT_FOUND -2
#define NOT_FIXED -3

typedef struct {
  i64 n; // 現在のBPの長さ
  int L; // ブロックの長さ
  int alpha; // 領域を増やす割合 >= 2, alpha/(a-1) 倍になる
  i64 size; // 現在のデータ構造のサイズ (バイト数)

  i64 N; // 確保したBPの領域の長さ
  bitvec_t *BP;
  i64 *depth_H; // depthの上位桁
  short *depth_L; // 上位桁からの差分

// RMM tree
  i64 heap_size_H; // ヒープ用に確保したメモリのサイズ
  int heap_height_H;
  int heap_height_L;
  i64 *RM_H;
  short **RM_L;
  i64 n_MB; // 中ブロックの数
  i64 N_MB; // 中ブロック用の領域の数

// 最後の括弧用
  i64 *match; // 各深さのマッチ対象の位置
  i64 maxdepth; // 
  i64 lastOP; // 最後の ( の場所
// オンライン構築用
  i64 d; // 最後の位置の括弧の深さ
  i64 block_min; // 現在のブロック内の深さの最小値

} onlinebp;

#ifdef __cplusplus
extern "C" {
#endif
////////////////////////////////////////////////////////////////
//        括弧列に対する操作
////////////////////////////////////////////////////////////////
onlinebp *onlinebp_new(int L, int alpha);
void onlinebp_reserve(onlinebp * bp, u64 capaBitLen);
void onlinebp_free(onlinebp *bp);
i64 onlinebp_length(onlinebp *bp);
i64 onlinebp_size(onlinebp *bp);
int onlinebp_alpha(onlinebp *bp);
int onlinebp_L(onlinebp *bp);

void onlinebp_push(onlinebp *bp, int c);
i64 onlinebp_findopen(onlinebp *bp, i64 i);
i64 onlinebp_findclose(onlinebp *bp, i64 i);
i64 onlinebp_enclose(onlinebp *bp, i64 i);
int onlinebp_get(onlinebp *bp, i64 i);
i64 onlinebp_depth(onlinebp *bp, i64 i);
i64 onlinebp_rank(onlinebp *bp, i64 i, int c);
i64 onlinebp_select(onlinebp *bp, i64 i, int c);

i64 onlinebp_bwd(onlinebp *bp, i64 i, i64 d);
i64 onlinebp_fwd(onlinebp *bp, i64 pos, i64 k);
int convbp(int c);

////////////////////////////////////////////////////////////////
//        postorder full binary tree (DFUDSの左右反転) に対する操作
////////////////////////////////////////////////////////////////
int onlinebp_rightchild(onlinebp *bp, i64 i, i64 *s, i64 *t);
int onlinebp_leftchild(onlinebp *bp, i64 i, i64 *s, i64 *t);
i64 onlinebp_parent(onlinebp *bp, i64 i);

int onlinebp_isrightchild(onlinebp *bp, i64 i);
int onlinebp_isleftchild(onlinebp *bp, i64 i);
i64 onlinebp_leftmost_leaf(onlinebp *bp, i64 i);
i64 onlinebp_rightmost_leaf(onlinebp *bp, i64 i);

#ifdef __cplusplus
}
#endif

#endif // _ONLINEBP_H_
