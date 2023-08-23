
// prototypes of functions in bbm.c
void byte_to_bbm(uint8_t *m, int msize, int i, int j, int size, uint8_t b);
void bbm_to_ascii(uint8_t *m, int msize, int i, int j, int size, FILE *f);
int mmult_bbm(const uint8_t *a, int size, const uint8_t *b, uint8_t *c);
bool mequals_bbm(const uint8_t *a, int size, const uint8_t *b);
uint8_t *bbm_read(char *name, int *psize);
void bbm_write(uint8_t *m, size_t msize, char *name);
uint8_t *bbm_alloc(size_t size);