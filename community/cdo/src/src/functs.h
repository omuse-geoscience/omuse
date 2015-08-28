#ifndef _FUNCTS_H
#define _FUNCTS_H

#define  func_fld       7
#define  func_all       8
#define  func_hrd       9

#define  func_min      10
#define  func_max      11
#define  func_range    12
#define  func_sum      13
#define  func_avg      14
#define  func_mean     15
#define  func_std      17
#define  func_std1     18
#define  func_var      19
#define  func_var1     20

#define  func_pctl     21

#define  func_cor      22
#define  func_covar    23

#define  func_crps     30
#define  func_brs      31
#define  func_rank     32
#define  func_roc      33

#define  func_add      41
#define  func_sub      42
#define  func_mul      43
#define  func_div      44
#define  func_mod      45

#define  func_atan2    50

#define  func_read     60
#define  func_write    61

#define  func_month    84
#define  func_year     85
#define  func_time     86
#define  func_date     87
#define  func_step     88
#define  func_datetime 89

#define  func_lon      98
#define  func_lat      99

enum cmp_flag {
  CMP_NAME     = 1,
  CMP_GRID     = 2,
  CMP_NLEVEL   = 4,
  CMP_GRIDSIZE = 8,
  CMP_HRD      = CMP_NAME | CMP_GRIDSIZE,
  CMP_DIM      =            CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID,
  CMP_ALL      = CMP_NAME | CMP_GRIDSIZE | CMP_NLEVEL | CMP_GRID
};

void    vlistCompare(int vlistID1, int vlistID2, int flag);
int     vlistCompareX(int vlistID1, int vlistID2, int flag);

#endif  /* _FUNCTS_H */
