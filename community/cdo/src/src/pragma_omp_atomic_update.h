#if defined(_OPENMP)
#if _OPENMP >= OPENMP4
#pragma omp atomic update
#else
#pragma omp atomic
#endif
#endif
