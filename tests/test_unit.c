#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "k2.h"

// Define a simple assertion macro for clear feedback
#define ASSERT_TRUE(cond) \
    do { \
        if (!(cond)) { \
            fprintf(stderr, "ASSERTION FAILED: %s at %s:%d\n", #cond, __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

#define ASSERT_EQ(actual, expected) \
    do { \
        size_t _a = (size_t)(actual); \
        size_t _e = (size_t)(expected); \
        if (_a != _e) { \
            fprintf(stderr, "ASSERTION FAILED: %s == %s (expected %zu, got %zu) at %s:%d\n", \
                    #actual, #expected, _e, _a, __FILE__, __LINE__); \
            exit(1); \
        } \
    } while(0)

// Standard 2x2 patterns (row-major)
static uint8_t PATTERN_EMPTY[4] = {0, 0, 0, 0};
static uint8_t PATTERN_FULL[4]  = {1, 1, 1, 1};
static uint8_t PATTERN_DIAG[4]  = {1, 0, 0, 1}; // Identity 2x2

void test_2x2_basic() {
    printf("Running test_2x2_basic...\n");
    k2mat_t a = K2MAT_INITIALIZER;
    
    // 1. Empty 2x2
    mread_from_bbm(PATTERN_EMPTY, 2, &a);
    ASSERT_TRUE(k2is_empty(&a));
    ASSERT_EQ(mget_nonzeros(&a), 0);
    
    // 2. Full 2x2 (embedded in 4x4)
    mread_from_bbm(PATTERN_FULL, 2, &a);
    ASSERT_EQ(mget_nonzeros(&a), 4);
    // Root of 4x4 has only top-left quadrant (bit 0) set
    ASSERT_EQ(k2read_node(&a, 0), 0x1);
    // Pos 1 is the 2x2 minimat, which is 0xF (all 4 bits set)
    ASSERT_EQ(k2read_node(&a, 1), 0xF); 
    
    // 3. Diagonal 2x2
    mread_from_bbm(PATTERN_DIAG, 2, &a);
    ASSERT_EQ(mget_nonzeros(&a), 2);
    // Root of 4x4 has only top-left quadrant (bit 0) set
    ASSERT_EQ(k2read_node(&a, 0), 0x1);
    // Pos 1 is the 2x2 minimat for the diagonal: 0x9 (bit 0 and 3 set)
    ASSERT_EQ(k2read_node(&a, 1), 0x9); 
    
    matrix_free(&a);
    printf("test_2x2_basic PASSED.\n");
}

void test_4x4_identity() {
    printf("Running test_4x4_identity...\n");
    uint8_t m[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    };
    k2mat_t a = K2MAT_INITIALIZER;
    mread_from_bbm(m, 4, &a);
    
    ASSERT_EQ(mget_nonzeros(&a), 4);
    
    // For a 4x4 identity, the quadrants are:
    // Q0 (top-left): 2x2 identity
    // Q1 (top-right): 2x2 zero
    // Q2 (bottom-left): 2x2 zero
    // Q3 (bottom-right): 2x2 identity
    // So the 4x4 root should have bits 0 and 3 set -> 0x9
    ASSERT_EQ(k2read_node(&a, 0), 0x9);
    
    // Navigation: 
    // Pos 1 should be the minimat for Q0 (identity 2x2 -> 0x9)
    ASSERT_EQ(k2read_node(&a, 1), 0x9);
    // Pos 2 should be the minimat for Q3 (identity 2x2 -> 0x9)
    ASSERT_EQ(k2read_node(&a, 2), 0x9);
    
    matrix_free(&a);
    printf("test_4x4_identity PASSED.\n");
}

void test_all_ones_compression() {
    printf("Running test_all_ones_compression...\n");
    uint8_t m[16] = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1}; // 4x4 all ones
    k2mat_t a = K2MAT_INITIALIZER;
    mread_from_bbm(m, 4, &a);
    
    ASSERT_EQ(mget_nonzeros(&a), 16);
    // For a 4x4 wholly filled with 1s, it should be compressed to a single ALL_ONES root
    ASSERT_EQ(k2read_node(&a, 0), ALL_ONES);
    ASSERT_EQ(k2treesize(&a), 1);
    
    matrix_free(&a);
    printf("test_all_ones_compression PASSED.\n");
}

int main() {
    // Required initialization
    minimat_init(2);
    
    test_2x2_basic();
    test_4x4_identity();
    test_all_ones_compression();
    
    printf("\nALL TESTS PASSED SUCCESSFULLY!\n");
    return 0;
}
