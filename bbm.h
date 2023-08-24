
// prototypes of functions in bbm.c
uint8_t *bbm_alloc(size_t size);
void bbm_write(const uint8_t *m, size_t msize, const char *name);
uint8_t *bbm_read(const char *name, size_t *psize);
void byte_to_bbm(uint8_t *m, size_t msize, int i, int j, int size, uint8_t b);
void bbm_to_ascii(const uint8_t *m, size_t msize, int i, int j, int size, FILE *f);
int mmult_bbm(const uint8_t *a, size_t size, const uint8_t *b, uint8_t *c);
ssize_t mequals_bbm(const uint8_t *a, size_t size, const uint8_t *b);
