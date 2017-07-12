#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

#define PROFILE
#include "rs.h"
#include "rs.c"

void print_matrix1(gf* matrix, int nrows, int ncols);
void print_matrix2(gf** matrix, int nrows, int ncols);

void print_buf(gf* buf, char *fmt, size_t len) {
    size_t i = 0;
    while(i < len) {
        printf(fmt, buf[i]);
        i++;
        if((i % 50) == 0) {
            printf("\n");
        }
    }
    printf("\n");
}

void test_galois(void) {
    printf("%s:\n", __FUNCTION__);

    //copy from golang rs version
    assert(galMultiply(3, 4) == 12);
    assert(galMultiply(7, 7) == 21);
    assert(galMultiply(23, 45) == 41);

    {
        gf in[] = {0, 1, 2, 3, 4, 5, 6, 10, 50, 100, 150, 174, 201, 255, 99, 32, 67, 85};
        gf out[sizeof(in)/sizeof(gf)] = {0};
        gf expect[] = {0x0, 0x19, 0x32, 0x2b, 0x64, 0x7d, 0x56, 0xfa, 0xb8, 0x6d, 0xc7, 0x85, 0xc3, 0x1f, 0x22, 0x7, 0x25, 0xfe};
        gf expect2[] = {0x0, 0xb1, 0x7f, 0xce, 0xfe, 0x4f, 0x81, 0x9e, 0x3, 0x6, 0xe8, 0x75, 0xbd, 0x40, 0x36, 0xa3, 0x95, 0xcb};
        int rlt;
        addmul(out, in, 25, sizeof(int)/sizeof(gf));
        rlt = memcmp(out, expect, sizeof(int)/sizeof(gf));
        assert(0 == rlt);

        memset(out, 0, sizeof(in)/sizeof(gf));
        addmul(out, in, 177, sizeof(in)/sizeof(gf));
        rlt = memcmp(out, expect2, sizeof(int)/sizeof(gf));
        assert(0 == rlt);
    }

    assert(galExp(2,2) == 4);
    assert(galExp(5,20) == 235);
    assert(galExp(13,7) == 43);
}

void test_sub_matrix(void) {
    int r, c, ptr, nrows = 10, ncols = 20;
    gf* m1 = (gf*)RS_MALLOC(nrows * ncols);
    gf *test1;

    printf("%s:\n", __FUNCTION__);

    ptr = 0;
    for(r = 0; r < nrows; r++) {
        for(c = 0; c < ncols; c++) {
            m1[ptr] = ptr;
            ptr++;
        }
    }
    test1 = sub_matrix(m1, 0, 0, 3, 4, nrows, ncols);
    for(r = 0; r < 3; r++) {
        for(c = 0; c < 4; c++) {
            assert(test1[r*4 + c] == (r*ncols + c));
        }
    }
    free(test1);

    test1 = sub_matrix(m1, 3, 2, 7, 9, nrows, ncols);
    for(r = 0; r < (7-3); r++) {
        for(c = 0; c < (9-2); c++) {
            assert(test1[r*(9-2) + c] == ((r+3)*ncols + (c+2)));
        }
    }

    free(m1);
}

void test_multiply(void) {
    gf a[] = {1,2,3,4};
    gf b[] = {5,6,7,8};
    gf exp[] = {11,22,19,42};
    gf *out;
    int rlt;

    printf("%s:\n", __FUNCTION__);

    out = multiply1(a, 2, 2, b, 2, 2);
    rlt = memcmp(out, exp, 4);
    assert(0 == rlt);
}

void test_inverse(void) {
    printf("%s:\n", __FUNCTION__);
    {
        gf a[] = {56, 23, 98, 3, 100, 200, 45, 201, 123};
        gf ae[] = {175, 133, 33, 130, 13, 245, 112, 35, 126};
        int rlt = invert_mat(a, 3);
        assert(0 == rlt);
        rlt = memcmp(a, ae, 3*3);
        assert(0 == rlt);
    }

    {
        gf a[] = {  1, 0, 0, 0, 0,
                    0, 1, 0, 0, 0,
                    0, 0, 0, 1, 0,
                    0, 0, 0, 0, 1,
                    7, 7, 6, 6, 1};
        gf ae[] = {1, 0, 0, 0, 0,
                   0, 1, 0, 0, 0,
                   123, 123, 1, 122, 122,
                   0, 0, 1, 0, 0,
                   0, 0, 0, 1, 0};
        int rlt = invert_mat(a, 5);
        assert(0 == rlt);
        rlt = memcmp(a, ae, 5*5);
        assert(0 == rlt);
    }

    {
        /* error matrix */
        gf a[] = {4,2,12,6};
        int rlt = invert_mat(a, 2);
        assert(0 != rlt);
    }
}

unsigned char* test_create_random(reed_solomon *rs, int data_size, int block_size) {
    struct timeval tv;
    unsigned char* data;
    int i, n, seed, nr_blocks;

    gettimeofday(&tv, 0);
    seed = tv.tv_sec ^ tv.tv_usec;
    srandom(seed);

    nr_blocks = (data_size+block_size-1)/block_size;
    nr_blocks = ((nr_blocks + rs->data_shards - 1)/ rs->data_shards) * rs->data_shards;
    n = nr_blocks / rs->data_shards;
    nr_blocks += n * rs->parity_shards;

    data = malloc(nr_blocks * block_size);
    for(i = 0; i < data_size; i++) {
        data[i] = (unsigned char)(random() % 255);
    }
    memset(data + data_size, 0, nr_blocks*block_size - data_size);

    return data;
}

int test_create_encoding(
        reed_solomon *rs,
        unsigned char *data,
        int data_size,
        int block_size
        ) {
    unsigned char **data_blocks;
    int data_shards, parity_shards;
    int i, n, nr_shards, nr_blocks, nr_fec_blocks;

    data_shards = rs->data_shards;
    parity_shards = rs->parity_shards;
    nr_blocks = (data_size+block_size-1)/block_size;
    nr_blocks = ((nr_blocks+data_shards-1)/data_shards) * data_shards;
    n = nr_blocks / data_shards;
    nr_fec_blocks = n * parity_shards;
    nr_shards = nr_blocks + nr_fec_blocks;

    data_blocks = (unsigned char**)malloc(nr_shards * sizeof(unsigned char*));
    for(i = 0; i < nr_shards; i++) {
        data_blocks[i] = data + i*block_size;
    }

    n = reed_solomon_encode2(rs, data_blocks, nr_shards, block_size);
    free(data_blocks);

    return n;
}

int test_data_decode(
        reed_solomon *rs,
        unsigned char *data,
        int data_size,
        int block_size,
        int *erases,
        int erase_count) {
    unsigned char **data_blocks;
    unsigned char *zilch;
    int data_shards, parity_shards;
    int i, j, n, nr_shards, nr_blocks, nr_fec_blocks;

    data_shards = rs->data_shards;
    parity_shards = rs->parity_shards;
    nr_blocks = (data_size+block_size-1)/block_size;
    nr_blocks = ((nr_blocks+data_shards-1)/data_shards) * data_shards;
    n = nr_blocks / data_shards;
    nr_fec_blocks = n * parity_shards;
    nr_shards = nr_blocks + nr_fec_blocks;

    data_blocks = (unsigned char**)malloc(nr_shards * sizeof(unsigned char*));
    for(i = 0; i < nr_shards; i++) {
        data_blocks[i] = data + i*block_size;
    }

    zilch = (unsigned char*)calloc(1, nr_shards);
    for(i = 0; i < erase_count; i++) {
        j = erases[i];
        memset(data + j*block_size, 137, block_size);
        zilch[j] = 1; //mark as erased
    }

    n = reed_solomon_reconstruct(rs, data_blocks, zilch, nr_shards, block_size);
    free(data_blocks);
    free(zilch);

    return n;
}

void test_one_encoding(void) {
    reed_solomon *rs;
    unsigned char* data;
    int block_size = 50000;
    int data_size = 10*block_size;
    int err;

    printf("%s:\n", __FUNCTION__);

    rs = reed_solomon_new(10, 3);
    data = test_create_random(rs, data_size, block_size);
    err = test_create_encoding(rs, data, data_size, block_size);

    free(data);
    reed_solomon_release(rs);

    assert(0 == err);
}

int test_one_decoding_13(int *erases, int erase_count) {
    reed_solomon *rs;
    unsigned char *data, *origin;
    int block_size = 50000;
    int data_size = 10*block_size;
    int err, err2;

    rs = reed_solomon_new(10, 3);
    data = test_create_random(rs, data_size, block_size);
    err = test_create_encoding(rs, data, data_size, block_size);
    assert(0 == err);

    origin = (unsigned char*)malloc(data_size);
    memcpy(origin, data, data_size);

    err = test_data_decode(rs, data, data_size, block_size, erases, erase_count);
    if(0 == err) {
        err2 = memcmp(origin, data, data_size);
        assert(0 == err2);
    } else {
        //failed here
        err2 = memcmp(origin, data, data_size);
        assert(0 != err2);
    }

    free(data);
    free(origin);
    reed_solomon_release(rs);

    return err;
}

void test_one_decoding(void) {
    printf("%s:\n", __FUNCTION__);

    {
        int erases[] = {0};
        int err;

        // lost nothing
        err = test_one_decoding_13(erases, 0);
        assert(0 == err);
    }

    {
        int erases[] = {0};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // lost only one
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 5;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 9;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 11;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);
    }

    {
        int erases[] = {0, 1};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // lost two
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 3;
        erases[1] = 7;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 11;
        erases[1] = 9;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 11;
        erases[1] = 12;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);
    }

    {
        int erases[] = {0, 1, 4};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // lost three
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 3;
        erases[1] = 8;
        erases[2] = 7;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 11;
        erases[1] = 9;
        erases[2] = 1;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 11;
        erases[1] = 12;
        erases[2] = 9;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 11;
        erases[1] = 12;
        erases[2] = 9;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);

        erases[0] = 11;
        erases[1] = 12;
        erases[2] = 10;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 == err);
    }

    {
        int erases[] = {0, 1, 4, 8};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // lost 4, failed!
        err = test_one_decoding_13(erases, erases_count);
        assert(0 != err);

        erases[0] = 11;
        erases[1] = 12;
        erases[2] = 10;
        erases[3] = 9;
        err = test_one_decoding_13(erases, erases_count);
        assert(0 != err);
    }
}

int test_one_decoding_13_6(int *erases, int erase_count) {
    reed_solomon *rs;
    unsigned char *data, *origin;
    int block_size = 50000;
    int data_size = 10*block_size*6;
    int err, err2;

    rs = reed_solomon_new(10, 3);
    data = test_create_random(rs, data_size, block_size);
    err = test_create_encoding(rs, data, data_size, block_size);
    assert(0 == err);

    origin = (unsigned char*)malloc(data_size);
    memcpy(origin, data, data_size);

    err = test_data_decode(rs, data, data_size, block_size, erases, erase_count);
    if(0 == err) {
        err2 = memcmp(origin, data, data_size);
        assert(0 == err2);
    } else {
        //failed here
        err2 = memcmp(origin, data, data_size);
        assert(0 != err2);
    }

    free(data);
    free(origin);
    reed_solomon_release(rs);

    return err;
}

void test_encoding(void) {
    reed_solomon *rs;
    unsigned char *data;
    int block_size = 50000;
    //multi shards encoding
    int data_size = 13*block_size*6;
    int err;

    printf("%s:\n", __FUNCTION__);

    rs = reed_solomon_new(10, 3);
    data = test_create_random(rs, data_size, block_size);
    err = test_create_encoding(rs, data, data_size, block_size);

    free(data);
    reed_solomon_release(rs);

    assert(0 == err);
}

void test_reconstruct(void) {
#define FEC_START (10*6)
    printf("%s:\n", __FUNCTION__);

    {
        int erases[] = {0};
        int err;

        // lost nothing
        err = test_one_decoding_13_6(erases, 0);
        assert(0 == err);
    }

    {
        int erases[] = {0, 1, 9, 10+2, 10+4, 10+9};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // shard1 shard2 both lost three 
        err = test_one_decoding_13_6(erases, erases_count);
        assert(0 == err);
    }

    {
        int erases[] = {0, 9, FEC_START + 1, 10+2, 10+9, FEC_START + 3 + 2};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // shard1 shard2 both lost three, and both lost one in fec
        err = test_one_decoding_13_6(erases, erases_count);
        assert(0 == err);
    }

    {
        int erases[] = {11, 12, 10, 9};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        /* this is ok. not lost 4 but shard1 lost 1, shard2 lost 3, we can reconstruct it */
        err = test_one_decoding_13_6(erases, erases_count);
        assert(0 == err);
    }

    {
        int erases[] = {0, 1, 4, 8};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // shard1 lost 4, failed!
        err = test_one_decoding_13_6(erases, erases_count);
        assert(0 != err);
    }

    {
        int erases[] = {10, 11, 14, 18};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // shard2 lost 4, failed!
        err = test_one_decoding_13_6(erases, erases_count);
        assert(0 != err);
    }

    {
        int erases[] = {0, 1, 4, 8, 10, 11, 14, 18};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        // shard1 and shard2 both lost 4, failed!
        err = test_one_decoding_13_6(erases, erases_count);
        assert(0 != err);
    }

    {
        int erases[] = {11, 12, 10, 9, FEC_START+3+0, FEC_START+3+1, FEC_START+3+2};
        int erases_count = sizeof(erases)/sizeof(int);
        int err;

        err = test_one_decoding_13_6(erases, erases_count);
        assert(0 != err);
    }
}

double benchmarkEncodeTest(int n, int dataShards, int parityShards, int shardSize) {
    clock_t start, end;
    double millis;
    unsigned char* data;
    int i;
    int dataSize = shardSize*dataShards;
    reed_solomon* rs = reed_solomon_new(dataShards, parityShards);

    data = test_create_random(rs, dataSize, shardSize);

    start = clock();
    for(i = 0; i < n; i++) {
        test_create_encoding(rs, data, dataSize, shardSize);
    }
    end = clock();
    millis = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
    return (millis);
}

/* TODO please check, is this benchmark ok? */
void benchmarkEncode(void) {
    double millis;
    double per_sec_in_bytes;
    double MB = 1024.0 * 1024.0;
    double millis_to_sec = 1000*1000;
    int n;
    int size;

    printf("%s:\n", __FUNCTION__);

    n = 20000;
    size = 10000;
    millis = benchmarkEncodeTest(n, 10, 2, size);

    per_sec_in_bytes = (10*2*size/MB) * millis_to_sec * n / millis;
    printf("10x2x10000, test_count=%d millis=%lf per_sec_in_bytes=%lfMB/s\n", n, millis, per_sec_in_bytes);

    n = 200;
    millis = benchmarkEncodeTest(n, 100, 20, size);
    per_sec_in_bytes = (100*20*size/MB) * millis_to_sec * n / millis;
    printf("100x20x10000, test_count=%d millis=%lf per_sec_in_bytes=%lfMB/s\n", n, millis, per_sec_in_bytes);

    n = 200;
    size = 1024*1024;
    millis = benchmarkEncodeTest(n, 17, 3, size);
    per_sec_in_bytes = (17*3*size/MB) * millis_to_sec * n / millis;
    printf("17x3x(1024*1024), test_count=%d millis=%lf per_sec_in_bytes=%lfMB/s\n", n, millis, per_sec_in_bytes);
}

void test_001(void) {
    reed_solomon* rs = reed_solomon_new(11, 6);
    print_matrix1(rs->m, rs->data_shards, rs->data_shards);
    print_matrix1(rs->parity, rs->parity_shards, rs->data_shards);
    reed_solomon_release(rs);
}

void test_002(void) {
    char text[] = "hello world", output[256];
    int block_size = 1;
    int nrDataBlocks = sizeof(text)/sizeof(char) - 1;
    unsigned char* data_blocks[128];
    unsigned char* fec_blocks[128];
    int nrFecBlocks = 6;

    //decode
    unsigned int fec_block_nos[128], erased_blocks[128];
    unsigned char* dec_fec_blocks[128];
    int nr_fec_blocks;

    int i;
    reed_solomon* rs = reed_solomon_new(nrDataBlocks, nrFecBlocks);

    printf("%s:\n", __FUNCTION__);

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&text[i];
    }

    memset(output, 0, sizeof(output));
    memcpy(output, text, nrDataBlocks);
    for(i = 0; i < nrFecBlocks; i++) {
        fec_blocks[i] = (unsigned char*)&output[i + nrDataBlocks];
    }

    reed_solomon_encode(rs, data_blocks, fec_blocks, block_size);
    print_buf((gf*)output, "%d ", nrFecBlocks+nrDataBlocks);

    text[1] = 'x';
    text[3] = 'y';
    text[4] = 'z';
    erased_blocks[0] = 4;
    erased_blocks[1] = 1;
    erased_blocks[2] = 3;

    fec_block_nos[0] = 1;
    fec_block_nos[1] = 3;
    fec_block_nos[2] = 5;
    dec_fec_blocks[0] = fec_blocks[1];
    dec_fec_blocks[1] = fec_blocks[3];
    dec_fec_blocks[2] = fec_blocks[5];
    nr_fec_blocks = 3;

    printf("erased:%s\n", text);

    reed_solomon_decode(rs, data_blocks, block_size, dec_fec_blocks,
            fec_block_nos, erased_blocks, nr_fec_blocks);

    printf("fixed:%s\n", text);

    reed_solomon_release(rs);
}

void test_003(void) {
    char text[] = "hello world hello world ", output[256];
    int block_size = 2;
    int nrDataBlocks = (sizeof(text)/sizeof(char) - 1) / block_size;
    unsigned char* data_blocks[128];
    unsigned char* fec_blocks[128];
    int nrFecBlocks = 6;

    //decode
    unsigned int fec_block_nos[128], erased_blocks[128];
    unsigned char* dec_fec_blocks[128];
    int nr_fec_blocks;

    int i;
    printf("nrDataBlocks [%d] sizeof(text) [%d] sizeof(char) [%d]\n",nrDataBlocks,sizeof(text),sizeof(char));
    reed_solomon* rs = reed_solomon_new(nrDataBlocks, nrFecBlocks);

    printf("%s:\n", __FUNCTION__);
    printf("text size=%d\n", (int)(sizeof(text)/sizeof(char) - 1) );

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&text[i*block_size];
    }

    memset(output, 0, sizeof(output));
    memcpy(output, text, nrDataBlocks*block_size);
    print_matrix1((gf*)output, nrDataBlocks + nrFecBlocks, block_size);
    for(i = 0; i < nrFecBlocks; i++) {
        fec_blocks[i] = (unsigned char*)&output[i*block_size + nrDataBlocks*block_size];
        printf("fec_blocks[%d] = [%s]\n",i,fec_blocks[i]);
    }

    for(i = 0; i < nrDataBlocks; i++) {
        printf("1.data_blocks[%d] = [%s]\n",i,data_blocks[i]);
    }

    printf("1-1.data_blocks = [%s]\n",data_blocks[0][1]);
    reed_solomon_encode(rs, data_blocks, fec_blocks, block_size);
    for(i = 0; i < nrDataBlocks; i++) {
        printf("2.data_blocks[%d] = [%s]\n",i,data_blocks[i]);
    }

    //printf("golang output(example/test_rs.go):\n [[104 101] [108 108] [111 32] [119 111] [114 108] [100 32] [104 101] [108 108] [111 32] [119 111] [114 108] [100 32] \n[157 178] [83 31] [48 240] [254 93] [31 89] [151 184]]\n");
    //printf("c verion output:\n");
    print_buf((gf*)output, "%d ", nrFecBlocks*block_size + nrDataBlocks*block_size);
    //print_matrix1((gf*)output, nrDataBlocks + nrFecBlocks, block_size);

    //decode
    text[1*block_size] = 'x';
    text[10*block_size+1] = 'y';
    text[4*block_size] = 'z';
    erased_blocks[0] = 1;
    erased_blocks[1] = 4;
    erased_blocks[2] = 10;

    fec_block_nos[0] = 3;
    fec_block_nos[1] = 2;
    fec_block_nos[2] = 4;
    dec_fec_blocks[0] = fec_blocks[3];
    dec_fec_blocks[1] = fec_blocks[2];
    dec_fec_blocks[2] = fec_blocks[4];
    nr_fec_blocks = 3;

    printf("erased:%s\n", text);

    reed_solomon_decode(rs, data_blocks, block_size, dec_fec_blocks,
            fec_block_nos, erased_blocks, nr_fec_blocks);

    printf("fixed:%s\n", text);

    reed_solomon_release(rs);
}

void test_004(void) {
    //char text[] = "hello world hello world ";
    int dataShards = 30;
    int parityShards = 21;
    int blockSize = 1280;
    struct timeval tv;
    int i, j, n, seed, size, nrShards, nrBlocks, nrFecBlocks;
    unsigned char *origin, *data;
    unsigned char **data_blocks;
    unsigned char *zilch;
    reed_solomon *rs;

    gettimeofday(&tv, 0);
    seed = tv.tv_sec ^ tv.tv_usec;
    srandom(seed);

    fec_init();

    //size = sizeof(text)/sizeof(char)-1;
    size = 1280*10;
    origin = malloc(size);
    //memcpy(origin, text, size);
    for(i = 0; i < size; i++) {
        origin[i] = (unsigned char)(random() % 255);
    }

    nrBlocks = (size+blockSize-1) / blockSize;
    nrBlocks = ((nrBlocks+dataShards-1)/dataShards) * dataShards;
    n = nrBlocks / dataShards;
    nrFecBlocks = n*parityShards;
    nrShards = nrBlocks + nrFecBlocks;
    data = malloc(nrShards * blockSize);
    memcpy(data, origin, size);
    memset(data + size, 0, nrShards*blockSize - size);
    printf("nrBlocks=%d nrFecBlocks=%d nrShards=%d n=%d left=%d\n", nrBlocks, nrFecBlocks, nrShards, n, nrShards*blockSize - size);
    //printf("origin:\n");
    //print_buf(origin, "%d ", size);
    //printf("data:\n");
    //print_buf(data, "%d ", nrShards*blockSize);

    data_blocks = (unsigned char**)malloc( nrShards * sizeof(unsigned char**) );
    for(i = 0; i < nrShards; i++) {
        data_blocks[i] = data + i*blockSize;
    }

    rs = reed_solomon_new(dataShards, parityShards);
    reed_solomon_encode2(rs, data_blocks, nrShards, blockSize);
    i = memcmp(origin, data, size);
    assert(0 == i);
    //print_matrix2(data_blocks, nrShards, blockSize);

    zilch = (unsigned char*)calloc(1, nrShards);
    n = parityShards;

    /* int es[100];
    es[0] = 3;
    es[1] = 3;
    es[2] = 2;
    es[3] = 8; */

    for(i = 0; i < n-2; i++) {
        j = random() % (nrBlocks-1);
        //j = es[i];
        memset(data + j*blockSize, 0, blockSize);
        zilch[j] = 1; //erased!
        printf("1.erased %d\n", j);
    }
    if(nrFecBlocks > 2) {
        for(i = 0; i < 2; i++) {
            j = nrBlocks + (random() % nrFecBlocks);
            memset(data + j*blockSize, 0, blockSize);
            zilch[j] = 1;
            printf("2.erased %d\n", j);
        }
    }

    reed_solomon_reconstruct(rs, data_blocks, zilch, nrShards, blockSize);
    i = memcmp(origin, data, size);
    print_buf(origin, "%d ", nrBlocks);
    print_buf(data_blocks, "%d ", nrBlocks);
    printf("rlt=%d\n", i);
    assert(0 == i);

    free(origin);
    free(data);
    free(data_blocks);
    free(zilch);
    reed_solomon_release(rs);
}

void test_005(void) {
    char text[2001], output[100*2001];
    int block_size = 2001;
    int nrDataBlocks = 20;
    unsigned char* data_blocks[nrDataBlocks];
    unsigned char* fec_blocks[10];


    //decode
    unsigned int fec_block_nos[128], erased_blocks[128];
    unsigned char* dec_fec_blocks[128];
    int nr_fec_blocks;

    int i;
    reed_solomon* rs = reed_solomon_new(nrDataBlocks, 10);
    //create each packet buffer: size = 1280
    //create 10 packets
    memset(output, 0, sizeof(output));
    for(i=0;i<nrDataBlocks;i++){
        if(i%5==0){
        memset(output+i*block_size,'1',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==1){
        memset(output+i*block_size,'2',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==2){
        memset(output+i*block_size,'3',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==3){
        memset(output+i*block_size,'4',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==4){
        memset(output+i*block_size,'5',1000);
        memset(output+(i+1)*1000-1,0,1);
        }

    }
    //printf("0.encode data_blocks[%d] = [%s]\n",i,output);

    for(i = 0; i < nrDataBlocks; i++) {
        data_blocks[i] = (unsigned char*)&output[i*block_size];
        //printf("1.encode data_blocks[%d] = [%s]\n",i,data_blocks[i]);
    }

    for(i = 0; i < 10; i++) {
        fec_blocks[i] = (unsigned char*)&output[i*block_size + nrDataBlocks*block_size];
        //printf("fec_blocks[%d] = [%s]\n",i,fec_blocks[i]);
    }

printf("before encode\n");
    reed_solomon_encode(rs, data_blocks, fec_blocks, block_size);
printf("after encode\n");
/*
    for(i = 0; i < nrDataBlocks; i++) {
        printf("encode 2 data_blocks[%d] = [%s]\n",i,data_blocks[i]);
    }
*/
    for(i = 0; i < 10; i++) {
        printf("encode 2 fec_blocks[%d] = [%s]\n",i,fec_blocks[i]);
    }

    //print_matrix1((gf*)output, nrDataBlocks + nrFecBlocks, block_size);
    //char testdata2[]="";
    //decode

    memset(data_blocks[4],0,2000);
    memset(data_blocks[5],0,2000);
    memset(data_blocks[1],0,2000);
    memset(data_blocks[0],0,2000);
    memset(data_blocks[14],0,2000);
    //memset(data_blocks[18],0,2000);
    memset(data_blocks[11],0,2000);
    memset(data_blocks[13],0,2000);
    //memset(data_blocks[19],0,2000);
    memset(data_blocks[17],0,2000);
    printf("\n\n\n\n");
    for(i = 0; i < nrDataBlocks; i++) {
        printf("after erased data_blocks[%d] = [%s]\n",i,data_blocks[i]);
    }
    

    erased_blocks[0] = 0;
    erased_blocks[1] = 1;
    erased_blocks[2] = 4;
    erased_blocks[3] = 5;
    erased_blocks[4] = 11;
    erased_blocks[5] = 13;
    erased_blocks[6] = 14;
    erased_blocks[7] = 17;
    //erased_blocks[8] = 18;
    //erased_blocks[9] = 19;



    fec_block_nos[0] = 9;
    fec_block_nos[1] = 8;
    fec_block_nos[2] = 7;
    fec_block_nos[3] = 6;
    fec_block_nos[4] = 4;
    fec_block_nos[5] = 5;
    fec_block_nos[6] = 3;
    fec_block_nos[7] = 2;
    //fec_block_nos[8] = 1;
    //fec_block_nos[9] = 0;
    dec_fec_blocks[0] = fec_blocks[9];
    dec_fec_blocks[1] = fec_blocks[8];
    dec_fec_blocks[2] = fec_blocks[7];
    dec_fec_blocks[3] = fec_blocks[6];
    dec_fec_blocks[4] = fec_blocks[4];
    dec_fec_blocks[5] = fec_blocks[5];
    dec_fec_blocks[6] = fec_blocks[3];
    dec_fec_blocks[7] = fec_blocks[2];
    //dec_fec_blocks[8] = fec_blocks[1];
    //dec_fec_blocks[9] = fec_blocks[0];
    nr_fec_blocks = 8;

    printf("\n\n\n\n");

    reed_solomon_decode(rs, data_blocks, block_size, dec_fec_blocks,
            fec_block_nos, erased_blocks, nr_fec_blocks);

    //printf("fixed:%s\n", output);
    for(i = 0; i < nrDataBlocks; i++) {
        printf("after decode data_blocks[%d] = [%s]\n",i,data_blocks[i]);
    }

    reed_solomon_release(rs);

}


//AVAPI-FEC api 
char output[100*2001];
typedef struct _datablock
{
    struct _datablock *p_next;
    struct _datablock *p_right;
    struct _datablock *p_left;

    int    Serial_No;      //packet no in this frame
    int    data_size;
    int    isLoss;
    unsigned char        *p_buffer;

}block_t;

typedef struct _datablock_fifo
{

    block_t          *p_first;  
    block_t          *p_last;

    int    treeType;
    int     i_depth;
    int     i_Blocksize;
    int     i_totalBlocks;
}datablock_fifo;

void AVAPI_FEC_Encode(datablock_fifo* data_blocks, datablock_fifo* fecdata_blocks)
{
    block_t *blockbuff         = data_blocks->p_first;
    block_t *fecdata_blockbuff = fecdata_blocks->p_first;

    int data_blocks_size       = data_blocks->i_depth;
    int fecdata_blocks_size    = fecdata_blocks->i_depth;
    int blocksize              = data_blocks->i_Blocksize;
    int i =0;

    unsigned char* _pdata_blocks[20];
    unsigned char* _pfecdata_blocks[10];

    for(i = 0; i < data_blocks_size; i++)
    {
        if(blockbuff != NULL)
        {    
            _pdata_blocks[i] = blockbuff->p_buffer;
            blockbuff = blockbuff->p_next;
        }
    }

    for(i = 0; i < fecdata_blocks_size; i++)
    {
        if(fecdata_blockbuff != NULL)
        {    
            _pfecdata_blocks[i] = fecdata_blockbuff->p_buffer;
            fecdata_blockbuff = fecdata_blockbuff->p_next;
        }
    }

    reed_solomon* rs = reed_solomon_new(data_blocks_size, fecdata_blocks_size);
    reed_solomon_encode(rs, _pdata_blocks, _pfecdata_blocks, blocksize);
    reed_solomon_release(rs);
}



void AVAPI_FEC_Decode(datablock_fifo* data_blocks, datablock_fifo* fecdata_blocks)
{

    block_t *data_blockbuff = data_blocks->p_first;
    block_t *fecdata_blockbuff = fecdata_blocks->p_first;

    int data_blocks_size     = data_blocks->i_depth;
    int defec_blocks_size    = fecdata_blocks->i_depth;
    int total_blocks_size    = data_blocks->i_totalBlocks;
    int total_fecblocks_size = fecdata_blocks->i_totalBlocks;
    int block_size           = data_blocks->i_Blocksize;
    int i =0,j=0;

    unsigned char* _pdata_blocks[data_blocks_size];
    unsigned char* _pdec_fec_blocks[defec_blocks_size];
    unsigned int   fec_block_nos[defec_blocks_size], erased_blocks[20];

    for(i = 0;i < total_blocks_size ; i++)
    {
        if(data_blockbuff->isLoss)
        {    
            erased_blocks[j] = data_blockbuff->Serial_No;
            j++;
        }

        _pdata_blocks[i] = data_blockbuff->p_buffer;
        data_blockbuff = data_blockbuff->p_next;
 
    }    

    for(i = 0;i < defec_blocks_size; i++)
    {

        fec_block_nos[i] = fecdata_blockbuff->Serial_No;
        _pdec_fec_blocks[i] = fecdata_blockbuff->p_buffer;

        fecdata_blockbuff = fecdata_blockbuff->p_next;
    }   

    reed_solomon* rs = reed_solomon_new(total_blocks_size, total_fecblocks_size);
    reed_solomon_decode(rs, _pdata_blocks, block_size, _pdec_fec_blocks, fec_block_nos, erased_blocks, defec_blocks_size);
    reed_solomon_release(rs);

} 

void test_006(void) {

    
    int block_size = 2001;
    int nrDataBlocks = 20;
    int nrfecDataBlocks = 10;
    int i =0;


    memset(output, 0, sizeof(output));
    for(i=0;i<nrDataBlocks;i++)
    {
        if(i%5==0){
        memset(output+i*block_size,'1',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==1){
        memset(output+i*block_size,'2',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==2){
        memset(output+i*block_size,'3',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==3){
        memset(output+i*block_size,'4',block_size);
        memset(output+(i+1)*block_size-1,0,1);
        }
        else if(i%5==4){
        memset(output+i*block_size,'5',1000);
        memset(output+(i+1)*1000-1,0,1);
        }
    }


    datablock_fifo *data_blocks = malloc(sizeof(datablock_fifo));
    datablock_fifo *fecdata_locks = malloc(sizeof(fecdata_locks));

/***************data_blocks init parameter and data *******************/    

    for(i = 0; i < nrDataBlocks; i++)
    {
        block_t *datablock = malloc(sizeof(block_t));

        datablock->Serial_No = i;
        datablock->data_size = 2001;
        datablock->p_buffer = (unsigned char *)&output[i*block_size];

        if(i == 0)
        {    
            data_blocks->p_first = datablock;
            data_blocks->p_last = datablock;
        }
        else if(i > 0 && i < nrDataBlocks)
        {
            block_t *p_last = data_blocks->p_last;
            p_last->p_next = datablock;
            data_blocks->p_last = datablock;
        } 
        else
        {
            data_blocks->p_last = datablock;
            datablock->p_next = NULL;
        }
    }

    /******check block is ok********/
    block_t *datablocktest = data_blocks->p_first;

    for(i = 0; i < nrDataBlocks; i++)
    {
        if(datablocktest != NULL)
        {
            printf("datablock[%d] = [%s]\n",i,datablocktest->p_buffer);
            datablocktest = datablocktest->p_next;
        } 
    }

/***************fecdata_locks init parameter and data *******************/ 

    for(i = 0; i < nrfecDataBlocks; i++)
    {
        block_t *fecdatablock = malloc(sizeof(block_t));

        fecdatablock->Serial_No = i;
        fecdatablock->data_size = 2001;
        fecdatablock->p_buffer = (unsigned char *)&output[i*block_size + nrDataBlocks*block_size];

        if(i == 0)
        {    
            fecdata_locks->p_first = fecdatablock;
            fecdata_locks->p_last = fecdatablock;
        }
        else if(i > 0 && i < nrfecDataBlocks)
        {
            block_t *p_last = fecdata_locks->p_last;
            p_last->p_next = fecdatablock;
            fecdata_locks->p_last = fecdatablock;
        } 
        else
        {
            fecdata_locks->p_last = fecdatablock;
            fecdatablock->p_next = NULL;
        }
    }

    /******check fecblock is ok********/
   

    data_blocks->i_depth = nrDataBlocks;
    data_blocks->i_Blocksize = block_size;
    fecdata_locks->i_depth = nrfecDataBlocks;
    fecdata_locks->i_totalBlocks = nrfecDataBlocks;
    fecdata_locks->i_Blocksize = block_size;

    AVAPI_FEC_Encode(data_blocks,fecdata_locks);

    block_t *fecdatablocktest = fecdata_locks->p_first;
    for(i = 0; i < nrfecDataBlocks; i++)
    {
        if(fecdatablocktest != NULL)
        {
            printf("fecdatablock[%d] = [%s]\n",i,fecdatablocktest->p_buffer);
            fecdatablocktest = fecdatablocktest->p_next;
        } 
    }

    AVAPI_FEC_Decode(data_blocks,fecdata_locks);

}

int main(void) {
    fec_init();
/*
    test_galois();
    test_sub_matrix();
    test_multiply();
    test_inverse();
    test_one_encoding();
    test_one_decoding();
    test_encoding();
    test_reconstruct();
*/    
    printf("reach here means test all ok\n");

    //benchmarkEncode();

    //test_001();
    //test_002();
    //test_003();
    //test_005();
    test_006();

    return 0;
}
