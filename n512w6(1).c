#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
// 定义结构体和相关函数
typedef struct {
    uint8_t* s;         // 压缩后的向量s (-1, 0, 1)
    int s_length;       // s的长度
    uint8_t* phi_As;    // 压缩后的phi_As (0, 1)
    int phi_As_length;  // phi_As的长度
} VectorPair;

// 哈希表节点
typedef struct HashNode {
    VectorPair* vp;
    struct HashNode* next;
} HashNode;

// 哈希表
typedef struct {
    int size;
    HashNode** buckets;
} HashTable;

// 哈希函数
unsigned long hash_function(uint8_t* phi_As, int r, int table_size) {
    unsigned long hash = 5381;

    // 遍历 phi_As 的每一位
    for (int i = 0; i < r; i++) {
        int byte_index = i / 8;       // 计算字节索引
        int bit_offset = i % 8;       // 计算位偏移量
        int bit_value = (phi_As[byte_index] >> bit_offset) & 0x01; // 提取第 i 位的值

        // 将当前位的值纳入哈希计算
        hash = ((hash << 5) + hash) + bit_value;
    }

    return hash % table_size;
}


uint8_t encode_value(int value) {
    if (value == 0) return 0x00; // 00
    if (value == 1) return 0x01; // 01
    if (value == -1) return 0x03; // 11

    return 0x02; // 默认返回 0，更改：不是-1,0,1的均设置为2，便于最后筛选
}

void set_value(VectorPair* vp, int index, int value) {
    int byte_index = index / 4;
    int bit_offset = (index % 4) * 2;
    uint8_t encoded = encode_value(value);
    vp->s[byte_index] &= ~(0x03 << bit_offset); // 清空对应位置
    vp->s[byte_index] |= (encoded << bit_offset); // 设置新值
}

int get_value(VectorPair* vp, int index) {
    int byte_index = index / 4;
    int bit_offset = (index % 4) * 2;
    uint8_t packed_value = (vp->s[byte_index] >> bit_offset) & 0x03;
    if (packed_value == 0x00) return 0;
    if (packed_value == 0x01) return 1;
    if (packed_value == 0x03) return -1;
    return 2; // 默认返回 0
}

void init_vector(VectorPair* vp, int length) {
    vp->s_length = length;
    vp->s = (uint8_t*)calloc((length + 3) / 4, sizeof(uint8_t)); // 每 4 个元素需要 1 字节
    if (vp->s == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    vp->phi_As = NULL;
}

void print_vector(VectorPair* vp) {
    for (int i = 0; i < vp->s_length; i++) {
        int byte_index = i / 4;
        int bit_offset = (i % 4) * 2;
        uint8_t packed_value = (vp->s[byte_index] >> bit_offset) & 0x03;
        printf("%d ", (packed_value == 0x00) ? 0 : ((packed_value == 0x01) ? 1 : -1));
    }
    printf("\n");
}

void free_vector(VectorPair* vp) {
    if (vp == NULL) return; // 如果指针为空，直接返回

    // 释放 s 指针指向的内存
    if (vp->s != NULL) {
        free(vp->s);
        vp->s = NULL; // 避免悬空指针
    }
    if (vp->phi_As != NULL)
    {
        free(vp->phi_As);
        vp->phi_As = NULL;
    }

    // 清理其他字段（可选）
    //vp->s_length = 0;
    //vp->phi_As_length = 0;
}
// 辅助函数：复制一个向量到另一个向量
void copy_vector(VectorPair* dest, VectorPair* src) {
    dest->s_length = src->s_length;
    dest->s = (uint8_t*)malloc(((src->s_length + 3) / 4) * sizeof(uint8_t));
    if (dest->s == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    memcpy(dest->s, src->s, ((src->s_length + 3) / 4) * sizeof(uint8_t));
}

void copy_compressed_vector(uint8_t* dest, const uint8_t* src, int src_length, int dest_offset) {
    for (int i = 0; i < src_length; i++) {
        int byte_index_src = i / 4;       // 计算源数组中的字节索引
        int bit_offset_src = (i % 4) * 2; // 源数组中的位偏移量

        int value = (src[byte_index_src] >> bit_offset_src) & 0x03; // 解压值

        int byte_index_dest = (dest_offset + i) / 4;       // 计算目标数组中的字节索引
        int bit_offset_dest = ((dest_offset + i) % 4) * 2; // 目标数组中的位偏移量

        // 清除目标位置的旧值
        dest[byte_index_dest] &= ~(0x03 << bit_offset_dest);

        // 设置新值
        dest[byte_index_dest] |= (value << bit_offset_dest);
    }
}
void set_zeros(uint8_t* s, int start_index, int num_zeros) {
    for (int i = 0; i < num_zeros; i++) {
        int index = start_index + i; // 当前要设置为 0 的索引

        int byte_index = index / 4;       // 计算字节索引
        int bit_offset = (index % 4) * 2; // 计算位偏移量

        // 清除目标位置的旧值（设置为 0）
        s[byte_index] &= ~(0x03 << bit_offset);
    }
}

int combination(int n, int k) {
    if (k > n) return 0; // 不合法的情况
    if (k == 0 || k == n) return 1; // 边界情况 

    // 使用公式 C(n, k) 
    int result = 1;
    for (int i = 1; i <= k; i++) {
        result = result * (n - i + 1) / i;
    }
    return result;
}

// 回溯法生成向量
void generate_vectors_helper(VectorPair* current, int n, int w1, int w2, int index, int count_ones, int count_neg_ones, VectorPair** result, int* count, int max_count) {
    // 基本情况：如果已经处理完所有位置
    if (index == n) {
        if (count_ones == w1 && count_neg_ones == w2) { // 确保有 w1 个 1 和 w2 个 -1
            result[*count] = (VectorPair*)malloc(sizeof(VectorPair));
            init_vector(result[*count], n);
            copy_vector(result[*count], current); // 复制当前向量到结果列表
            (*count)++;
        }
        return;
    }

    // 如果 1 或 -1 的数量已经超过限制，则提前剪枝
    if (count_ones > w1 || count_neg_ones > w2) return;

    // 尝试在当前位置放置 -1, 0, 1
    for (int value = -1; value <= 1; value++) {
        set_value(current, index, value);
        int new_count_ones = count_ones + (value == 1 ? 1 : 0);
        int new_count_neg_ones = count_neg_ones + (value == -1 ? 1 : 0);
        generate_vectors_helper(current, n, w1, w2, index + 1, new_count_ones, new_count_neg_ones, result, count, max_count);
    }
}

// 主函数：生成所有满足条件的向量
VectorPair** generate_vectors(int n, int w1, int w2, int extra_zeros, bool add_to_front, int* out_count) {
    // 检查输入合法性
    if (w1 + w2 > n) {
        fprintf(stderr, "Error: w1 + w2 cannot exceed the vector length n.\n");
        exit(EXIT_FAILURE);
    }

    // 计算最大可能的向量数量（上界）
    int max_count = combination(n, w1) * combination(n - w1, w2);
    printf("max:%d\n", max_count);

    // 分配结果数组
    VectorPair** result = (VectorPair**)malloc(max_count * sizeof(VectorPair*));
    if (result == NULL) {
        fprintf(stderr, "'generate_vectors' Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    // 初始化临时向量
    VectorPair current;
    init_vector(&current, n);

    // 开始生成向量
    *out_count = 0;
    generate_vectors_helper(&current, n, w1, w2, 0, 0, 0, result, out_count, max_count);

    // 在每个生成的向量前/后添加额外的 0
    for (int i = 0; i < *out_count; i++) {
        VectorPair* original = result[i];
        VectorPair* extended = (VectorPair*)malloc(sizeof(VectorPair));
        init_vector(extended, n + extra_zeros);

        if (add_to_front) {
            // 将 original 的内容拷贝到 extended 的后半部分
            copy_compressed_vector(extended->s, original->s, n, extra_zeros);
            set_zeros(extended->s, 0, extra_zeros);
        }
        else {
            // 将 original 的内容拷贝到 extended 的前半部分
            copy_compressed_vector(extended->s, original->s, n, 0);
            set_zeros(extended->s, n, extra_zeros); // 后面填充 0
        }

        // 更新长度
        extended->s_length = n + extra_zeros;

        // 替换原向量
        free_vector(original);
        result[i] = extended;
    }
    printf("2\n");
    // 释放临时向量内存
    free_vector(&current);

    // 返回结果
    return result;
}
// 初始化phi_As
void init_phi_As(VectorPair* vp, int length) {
    vp->phi_As_length = length;
    vp->phi_As = (uint8_t*)calloc((length + 7) / 8, sizeof(uint8_t)); // 每 8 个元素需要 1 字节
    if (vp->phi_As == NULL) {
        fprintf(stderr, "'generate_vectors' Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
}

// 设置phi_As中的某个值
void set_phi_As_value(VectorPair* vp, int index, int value) {
    int byte_index = index / 8;
    int bit_offset = index % 8;
    if (value == 1) {
        vp->phi_As[byte_index] |= (1 << bit_offset); // 设置对应比特为 1
    }
    else {
        vp->phi_As[byte_index] &= ~(1 << bit_offset); // 设置对应比特为 0
    }
}

// 获取phi_As中的某个值
int get_phi_As_value(VectorPair* vp, int index) {
    int byte_index = index / 8;
    int bit_offset = index % 8;
    return (vp->phi_As[byte_index] >> bit_offset) & 0x01;
}

// 打印phi_As
void print_phi_As(VectorPair* vp) {
    for (int i = 0; i < vp->phi_As_length; i++) {
        printf("%d ", get_phi_As_value(vp, i));
    }
    printf("\n");
}

// 矩阵和向量的乘法
int* matrix_vector_multiply(int** A, int rows, int cols, VectorPair* vp) {
    int* result = (int*)malloc(rows * sizeof(int));
    if (result == NULL) {
        fprintf(stderr, " 'generate_vectors' Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < rows; i++) {
        result[i] = 0;
        for (int j = 0; j < cols; j++) {
            result[i] += A[i][j] * get_value(vp, j);
        }
        result[i] = result[i] % 3329;
    }

    return result;
}

// 计算phi_As
void compute_phi_As(VectorPair* vp, int *v,int num,int q, int** A, int rows, int cols) {
    // 计算 A * s
    int* As = matrix_vector_multiply(A, num, cols, vp);
    for (int i = 0; i < num; i++)
    {
        As[i] = As[i] + v[i];
        As[i] = (As[i]+3329) % 3329;
    }
    // 根据 phi 函数规则计算 phi_As
    for (int i = 0; i < num; i++) {
        int xi = As[i];
        int hash_result = 0;

        // 确定边界值
        int lower_bound = (int)floor(q / 2.0) - 1;
        int upper_bound = q - 1;

        // 根据 phi 函数规则计算
        if (0 <= xi % q && xi % q < lower_bound) {
            hash_result = 0;
        }
        else if (lower_bound < xi % q && xi % q < upper_bound) {
            hash_result = 1;
        }
        else if (xi % q == lower_bound || xi % q == upper_bound) {
            hash_result = 1; // 边界情况同样赋予 1
        }

        // 设置 phi_As 的值
        set_phi_As_value(vp, i, hash_result);
    }

    // 释放临时分配的内存
    free(As);
}
// 向量相加
VectorPair* add_vectors(VectorPair* vp1, VectorPair* vp2) {
    if (vp1->s_length != vp2->s_length) {
        fprintf(stderr, "Vector lengths do not match\n");
        exit(EXIT_FAILURE);
    }

    VectorPair* result = (VectorPair*)malloc(sizeof(VectorPair));
    init_vector(result, vp1->s_length);

    for (int i = 0; i < vp1->s_length; i++) {
        int val1 = get_value(vp1, i);
        int val2 = get_value(vp2, i);
        int sum = val1 + val2;

        // 确保结果仍然在 -1, 0, 1 范围内
        if (sum > 1) sum = 2;
        if (sum < -1) sum = -2;

        set_value(result, i, sum);
    }

    return result;
}


// 创建哈希表
HashTable* create_hash_table(int size) {
    HashTable* table = (HashTable*)malloc(sizeof(HashTable));
    table->size = size;
    table->buckets = (HashNode**)calloc(size, sizeof(HashNode*));
    return table;
}

// 插入哈希表
void insert_hash_table(HashTable* table, VectorPair* vp, int r) {
    unsigned long index = hash_function(vp->phi_As, r, table->size);

    HashNode* new_node = (HashNode*)malloc(sizeof(HashNode));
    new_node->vp = vp;
    new_node->next = table->buckets[index]; // 链接到现有链表头部
    table->buckets[index] = new_node;
}

// 查找哈希表
VectorPair** find_in_hash_table(HashTable* table, uint8_t* phi_As, int r, int* count) {
    unsigned long index = hash_function(phi_As, r, table->size);
    VectorPair** res = NULL;
    *count = 0;

    HashNode* node = table->buckets[index];
    while (node != NULL) {
        if (memcmp(node->vp->phi_As, phi_As, (r + 7) / 8) == 0) { // 比较前 r 位
            res = realloc(res, (*count + 1) * sizeof(VectorPair*));
            res[*count] = node->vp;
            (*count)++;
        }
        node = node->next;
    }
    return res;
}

// 释放哈希表
void free_hash_table(HashTable* table) {
    for (int i = 0; i < table->size; ++i) {
        HashNode* node = table->buckets[i];
        while (node != NULL) {
            HashNode* temp = node;
            node = node->next;
            free(temp);
        }
    }
    free(table->buckets);
    free(table);
}



// 处理两个列表
void process_lists(VectorPair** L1, int L1_size, VectorPair** L2, int L2_size, int r,int **A,int e, VectorPair*** result, int* result_size,bool final) {
    HashTable* ht = create_hash_table(L1_size);

    // 构建哈希表
    for (int i = 0; i < L1_size; ++i) {
        insert_hash_table(ht, L1[i], r);
    }

    // 初始化结果数组为空
    *result = NULL;
    int index =0;
    int capacity = 0; // 当前已分配的容量

    for (int j = 0; j < L2_size; ++j) {
        // 查找匹配项
        int count;
        VectorPair** match = find_in_hash_table(ht, L2[j]->phi_As, r, &count);
        if (match != NULL) {
            for (int i = 0; i < count; i++) {
                // 如果当前容量不足，扩展容量
                if (index >= capacity) {
                    capacity = (capacity == 0) ? 1 : capacity * 2; // 按指数增长，后续更改，固定大小增长
                    printf("now numbers:%d\n", capacity);
                    *result = (VectorPair**)realloc(*result, capacity * sizeof(VectorPair*));
                }

                // 计算 s1 + s2
                int value;
                VectorPair* sum_vector = add_vectors(match[i], L2[j]);
                if (final == false)
                {
                     value = 0;
                    for (int k = 0; k < 512; k++)
                    {
                        value = A[0][k] * get_value(sum_vector, k) + value;
                        value = value % 3329;
                    }
                    value = (value + e) % 3329;
                    if (value == 0)
                    {
                        // 将结果加入数组
                        (*result)[index] = sum_vector;
                        index++;
                    }
                }
                else
                {
                    (*result)[index] = sum_vector;
                    index++;
                }
                
                
            }
        }
    }

    *result_size = index;

    // 如果最终分配的内存大于实际需要的内存，可以缩小到实际大小
    if (index> 0 && index < capacity) {
        *result = (VectorPair**)realloc(*result, (index) * sizeof(VectorPair*));
    }

    free_hash_table(ht);
}

// 生成随机置换索引（Fisher-Yates 洗牌算法）
int* generate_permutation(int n) {
    int* indices = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        indices[i] = i;
    }

    srand(time(NULL)); // 注意：最好只调用一次！
    for (int i = n - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
    }

    return indices;
}

// 计算逆置换索引
int* compute_inverse_permutation(int* indices, int n) {
    int* inverse = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        inverse[indices[i]] = i;
    }
    return inverse;
}

// 对矩阵应用列置换
void permute_matrix(int n, int** A, int* indices) {
    int** temp = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        temp[i] = (int*)malloc(n * sizeof(int));
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            temp[i][j] = A[i][indices[j]];
        }
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            A[i][j] = temp[i][j];
        }
    }

    for (int i = 0; i < n; i++) free(temp[i]);
    free(temp);
}

// 对向量应用置换
void permute_vector(int n, int* b, int* indices) {
    int* temp = (int*)malloc(n * sizeof(int));
    for (int j = 0; j < n; j++) {
        temp[j] = b[indices[j]];
    }

    for (int j = 0; j < n; j++) {
        b[j] = temp[j];
    }

    free(temp);
}

// 恢复矩阵
void restore_matrix(int n, int** A, int* inverse_indices) {
    int** temp = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        temp[i] = (int*)malloc(n * sizeof(int));
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            temp[i][j] = A[i][inverse_indices[j]];
        }
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            A[i][j] = temp[i][j];
        }
    }

    for (int i = 0; i < n; i++) free(temp[i]);
    free(temp);
}

// 恢复向量
void restore_vector(int n, int* b, int* inverse_indices) {
    int* temp = (int*)malloc(n * sizeof(int));
    for (int j = 0; j < n; j++) {
        temp[j] = b[inverse_indices[j]];
    }

    for (int j = 0; j < n; j++) {
        b[j] = temp[j];
    }

    free(temp);
}
// 测试主函数



int main()
{
    int n = 512; // 向量长度
    int ocount;
    int hcount;
    int* b = (int*)malloc(n * sizeof(int));
    //int s[512] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    //int s[20] = { 0, 1, -1, -1, 1, 1, -1, 0, -1, 1, -1, 1, 1, 0, 1, -1, -1, -1, 0, 1 };
    int s[512] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int e[512] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int q = 3329; // Z_q 的大小
    // 初始化矩阵 A (示例)
   // int rows = n, cols = n;
    srand(time(NULL));
    //int** A = (int**)malloc(rows * sizeof(int*));
    //for (int i = 0; i < rows; i++) {
      //  A[i] = (int*)malloc(cols * sizeof(int));
        //for (int j = 0; j < cols; j++) {
          //  A[i][j] = rand() % 3329; // 示例随机值
        //}
    //}
    FILE* file = fopen("D:\\matrix.txt", "r");
    if (file == NULL) {
        printf("无法打开文件\n");
        return 1;
    }
    int rows, cols;
    // 首先读取数组的行数和列数
    if (fscanf_s(file, "%d %d", &rows, &cols) != 2) {
        printf("读取行列数错误\n");
        fclose(file);
        return 1;
    }
    // 动态分配二维数组
    int** A = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        A[i] = (int*)malloc(cols * sizeof(int));
    }

    // 读取数据
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (fscanf_s(file, "%d", &A[i][j]) != 1) {
                printf("读取数据错误\n");
                // 释放已分配的内存
                for (int k = 0; k <= i; k++) {
                    free(A[k]);
                }
                free(A);
                fclose(file);
                return 1;
            }
        }
    }
    fclose(file);
    int** A_1 = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        A_1[i] = (int*)malloc(cols * sizeof(int));
        for (int j = 0; j < cols; j++) {
            A_1[i][j] = -1 * A[i][j];// 示例随机值
            A_1[i][j] = (A_1[i][j] + 3329) % 3329;
        }
    }

    for (int i = 0; i < rows; i++)
    {
        b[i] = 0;
        for (int j = 0; j < cols; j++)
        {
            b[i] = b[i] + A[i][j] * s[j];
        }
        b[i] = b[i] + e[i];
        b[i] = (b[i] + 2 * 3329) % 3329;
    }//b=As+e
    printf("begin generate L1\n");
    bool add_to_front = true;
    VectorPair** L1 = generate_vectors(n / 2, 1, 1, n / 2, add_to_front, &ocount);
    printf("L1 s size:%d\n", ocount);
    VectorPair** L2 = generate_vectors(n / 2, 1, 1, n / 2, false, &ocount);
    printf("L2 s size:%d\n", ocount);
    VectorPair** L3 = generate_vectors(n / 2, 1, 1, n / 2, add_to_front, &ocount);
    printf("L3 s size:%d\n", ocount);
    VectorPair** L4 = generate_vectors(n / 2, 1, 1, n / 2, false, &ocount);
    printf("L4 s size:%d\n", ocount);
    for (int i = 0; i < ocount; i++)
    {
        init_phi_As(L1[i], 12);
        init_phi_As(L2[i], 12);
        init_phi_As(L3[i], 12);
        init_phi_As(L4[i], 12);
    }//对O方法建表的第二分位进行初始化
    int flag = 0;//开始建表比较
    int permutation_count = 0;
    int* indices = generate_permutation(n);
    int* inverse = compute_inverse_permutation(indices, n);
    while (flag == 0)
    {
        int t1 = 0;
        VectorPair** result_1;
        VectorPair** result_2;
        int result_size_1 = 0;
        int result_size_2 = 0;
        int s_c[512] = { 0 };
        for (int error = 0; error < 3; error++)
        {
            for (int i = 0; i < ocount; i++)
            {
                int AS = 0;
                int AS1 = 0;
                int b_AS = 0;
                int b_AS1 = 0;
                for (int j = 0; j < cols; j++)
                {
                    AS = AS + A[0][j] * get_value(L1[i], j);
                    AS = (AS + 3329) % 3329;
                    AS1 = AS1 + A_1[0][j] * get_value(L2[i], j);
                    AS1 = (AS1 + 3329) % 3329;
                    b_AS = b_AS + A[0][j] * get_value(L3[i], j) ;
                    b_AS = (b_AS + 3329) % 3329;
                    b_AS1 = b_AS1 + A_1[0][j] * get_value(L4[i], j);
                    b_AS1 = (b_AS1 + 3329) % 3329;
                }
                AS1 = (3329 + AS1 + t1) % 3329;
                b_AS = (b_AS + 3329+t1+error-1) % 3329;
                b_AS1 = (b_AS1 + b[0]+3329) % 3329;
                for (int t = 0; t < 12; t++) {
                    int bit_value = (AS >> t) & 0x01; // 提取第 i 位的值（0 或 1）
                    set_phi_As_value(L1[i], t, bit_value); // 设置到 phi_As 中
                    bit_value = (AS1 >> t) & 0x01;
                    set_phi_As_value(L2[i], t, bit_value);
                    bit_value = (b_AS >> t) & 0x01;
                    set_phi_As_value(L3[i], t, bit_value);
                    bit_value = (b_AS1 >> t) & 0x01;
                    set_phi_As_value(L4[i], t, bit_value);
                }
            }
            process_lists(L1, ocount, L2, ocount, 12, A, 0, &result_1, &result_size_1, true);
            printf("first size:%d\n", result_size_1);
            process_lists(L3, ocount, L4, ocount, 12, A_1, b[0] + 0, &result_2, &result_size_2, true);
            printf("second size:%d\n", result_size_2);
            int i;
            #pragma omp parallel for // 使用8
            for ( i = 0; i < result_size_1; i++)
            {
                
                init_phi_As(result_1[i], 64);
                compute_phi_As(result_1[i], s_c, 64, q, A, rows, cols);

            }
            printf("init 1 finish\n");
            #pragma omp parallel for  // 使用8
            for ( i = 0; i < result_size_2; i++)
            {
                init_phi_As(result_2[i], 64);
                compute_phi_As(result_2[i], b, 64, q, A_1, rows, cols);
            }
            printf("init 2 finish\n");
            VectorPair** final_res;
            int final_length = 0;
            process_lists(result_1, result_size_1, result_2, result_size_2, 64, A_1, b[0] + 0, &final_res, &final_length, true);
            int final_s[512] = { 0 };
            for (int i = 0; i < final_length; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    final_s[j] = get_value(final_res[i], j);
                }
                if (permutation_count != 0)
                {
                    restore_vector(n, final_s, inverse);
                }
                if (memcmp(final_s, s, sizeof(s)) == 0)
                {
                    printf("find\n");
                    flag = 1;
                    printf("进行置换次数:%d\n", permutation_count);
                    return 0;
                }
                free_vector(final_res[i]);
            }
            free(final_res);
            for (int l = 0; l < result_size_1; l++)
            {
                free_vector(result_1[l]);
                free(result_1[l]);
            }
            for (int l = 0; l < result_size_2; l++)
            {
                free_vector(result_2[l]);
                free(result_2[l]);
            }
        }
        if (permutation_count != 0)
        {
            restore_matrix(n, A, inverse);
            restore_vector(n, b, inverse);
        }
        int* indices = generate_permutation(n);
        int* inverse = compute_inverse_permutation(indices, n);
        permute_matrix(n, A, indices);
        permute_vector(n, b, indices);//对A和b进行列置换
        permutation_count++;
    }

}