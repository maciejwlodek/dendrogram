#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>

#define LEAF 1
#define INTERNAL 2

#define EQUAL 0
#define SUBSET 1
#define SUPERSET 2
#define NOT_COMPARABLE 3

#define MAX_CLUSTERS 2048
#define MAX_LINE_SIZE 2048

typedef struct cluster CLUSTER;
typedef struct list_of_datasets LIST;
struct cluster {
    CLUSTER* left_child;
    CLUSTER* right_child;
    double ward_height;
    int type;
    LIST* list;
};
struct list_of_datasets { //list of datasets in a cluster
    int len;
    int* array;
};

char tmp_buffer[MAX_LINE_SIZE];

//CLUSTER* root;
CLUSTER** clusters;
int num_clusters;
char lines[MAX_CLUSTERS][MAX_LINE_SIZE];
int constrain;
double max_height;

int debug_mode = 0;

void debug(FILE* file, const char* format, ...) {
    va_list argptr;
    va_start(argptr, format);
    if(debug_mode == 1) {
        vfprintf(stdout, format, argptr);
    }
    va_end(argptr);
}
void fini_list(LIST* list) { //free all memory associated with given list
    free(list->array);
    free(list);
}
void fini(CLUSTER* root) { //recursively free all memory associated with tree
    if(root->type == LEAF) {
        fini_list(root->list);
        free(root);
    }
    else {
        fini(root->left_child);
        fini(root->right_child);
        fini_list(root->list);
        free(root);
    }
}
void print_list(LIST* list) {
    for(int i=0; i < list->len; i++) {
        debug(stdout, "%d,", *(list->array + i));
    }
    debug(stdout, "\n");
}
void print_simple_tree(CLUSTER* root) {
    debug(stdout, "------------------\n");
    if(root->type == LEAF) {
        debug(stdout, "leaf cluster ");
        print_list(root->list);
        debug(stdout, "-----------------\n");
    } else {
        debug(stdout, "cluster ");
        print_list(root->list);
        debug(stdout, "left pointer = %p\n", root->left_child);
        debug(stdout, "left_child = ");
        print_list(root->left_child->list);
        debug(stdout, "right pointer = %p\n", root->right_child);
        debug(stdout, "right_child = ");
        print_list(root->right_child->list);
        print_simple_tree(root->left_child);
        print_simple_tree(root->right_child);

        debug(stdout, "----------------\n");
    }
    //print_list(root->left_child->list);
}

// compare two lists under the partial relation of subset
// we assume the lists are sorted in increasing order
int compare(LIST* list1, LIST* list2) {
    int len1 = list1->len;
    int len2 = list2->len;
    int min = (len1 < len2)? 1 : 2;
    if(len1 == len2) min = 0;
    int i=0;
    int j=0;
    debug(stdout, "min = %d\n", min);
    while(i < len1 && j < len2) {
        int x = *(list1->array + i);
        int y = *(list2->array + j);
        debug(stdout, "x = %d, y = %d\n", x,y);
        if(x != y) {
            if(min == 0) return NOT_COMPARABLE;
            else if(min == 1) {
                if(x < y) return NOT_COMPARABLE;
                while(x > y && j < len2 - 1) {
                    j++;
                    y = *(list2->array + j);
                }
                if(x != y) return NOT_COMPARABLE;
            }
            else if(min == 2) {
                if(y < x) return NOT_COMPARABLE;
                while(y > x && i < len1 - 1) {
                    i++;
                    x = *(list1->array + i);
                }
                if(y != x) return NOT_COMPARABLE;
            }
        }
        i++;
        j++;
    }
    if(min == 0) return EQUAL;
    return (min == 1)? SUBSET : SUPERSET;
}

//return list1 $\setminus$ list2, assuming list2 $\subset$ list1
//both lists must be sorted
LIST* difference(LIST* list1, LIST* list2) {
    debug(stdout, "l1: ");
    print_list(list1);
    debug(stdout, "l2: ");
    print_list(list2);
    debug(stdout, "allocating space for difference between two lists\n");
    LIST* diff = (LIST*) malloc(sizeof(LIST));
    debug(stdout, "setting lengths\n");
    int len1 = list1->len;
    int len2 = list2->len;
    debug(stdout, "len1 = %d, len2 = %d\n", len1, len2);
    diff->len = len1 - len2;
    debug(stdout, "callocating array\n");
    diff->array = (int*) calloc(diff->len, sizeof(int));
    int i=0;
    int j=0;
    int c=0;
    debug(stdout, "comparing lists\n");
    while(i < len1) {
        debug(stdout, "l1[%d] = %d, l2[%d] = %d, c = %d\n", i, *(list1->array + i), j, *(list2->array + j), c);
        while(i < len1 && *(list1->array + i) != *(list2->array + j)) {
            debug(stdout, "inner_loop: i = %d, diff[%d] = %d\n", i, c, *(diff->array + c));
            *(diff->array + c) = *(list1->array + i);
            i++;
            c++;
        }
        i++;
        if(j < len2 - 1) j++;
    }
    return diff;
}

void get_nth_line_in_file(char* filename, char* buffer, int n) {
    FILE* fp = fopen(filename, "r+");
    int i=0;
    while(i < n && fgets(buffer, MAX_LINE_SIZE, fp) != NULL) {
        i++;
    }
    fgets(buffer, MAX_LINE_SIZE, fp);
    buffer[strcspn(buffer, "\n")] = 0;
    fclose(fp);
}

void replace_nth_line_in_file(char* filename, char* replacement, int n) {
    n++; //awk is 1-indexed
    char cmd[MAX_LINE_SIZE];
    sprintf(cmd, "awk '{if(NR==%d) print \"%s\"; else print $0;}' %s > tmpreplace", n, replacement, filename);
    system(cmd);
    rename("tmpreplace", filename);
}

//convert ward height to discrete height in file
int convert_ward_height(double ward_height) {
    int conv = (int) (cbrt( ward_height/max_height ) * ( constrain - 4) + 4); // 4 -> maximum dataset number = 10^(4-1)-1
    return conv;
}

//append character until buffer reaches length height
char* append_to_height(char* buffer, char character, int start, int height) {
    for(int i=0; i<start; i++) {
        if(*(buffer+i) == 0) *(buffer+i) = ' ';
    }
    for(int i=start; i<height; i++) {
        *(buffer+i) = character;
    }
    *(buffer+height) = 0;
}
//append character n times
char* append(char* buffer, char character, int start, int n_times) {
    int len = strlen(buffer);
    for(int i=0; i<start; i++) {
        if(*(buffer+i) == 0 || i >= len) *(buffer+i) = ' ';
    }
    for(int i=start-1; i<start-1+n_times; i++) {
        *(buffer+i) = character;
    }
    *(buffer+start-1+n_times) = 0;
}
//recursively print subtree rooted at given cluster
int write_to_dendro_file(CLUSTER* cluster, char* filename, int line_no) {
    debug(stdout, "printing subtree: ");
    print_list(cluster->list);
    debug(stdout, "checking type of current cluster\n");
    if(cluster->type == LEAF) {
        debug(stdout, "leaf\n");
        sprintf(tmp_buffer, "%3d", *(cluster->list->array)); //3 -> maximum dataset number = 10^3-1
        debug(stdout, "tmp_buffer1 = %s\n", tmp_buffer);
        replace_nth_line_in_file(filename, tmp_buffer, line_no);
        return line_no + 2;
    }
    else {
        debug(stdout, "internal cluster\n");

        int height = convert_ward_height(cluster->ward_height);
        int lheight = convert_ward_height(cluster->left_child->ward_height);
        int rheight = convert_ward_height(cluster->right_child->ward_height);
        debug(stdout, "height = %d, lheight = %d, rheight = %d\n", height, lheight, rheight);

        debug(stdout, "stepping recursively to left subtree: new height = %d\n", height);
        int new_line_no = write_to_dendro_file(cluster->left_child, filename, line_no);
        int ave1 = (int) (((new_line_no + line_no)/2) - 1);
        debug(stdout, "new line number = %d, new average = %d\n", new_line_no, ave1);
        get_nth_line_in_file(filename, tmp_buffer, ave1);
        append_to_height(tmp_buffer, '-', lheight, height);

        debug(stdout, "tmp_buffer2 = %s\n", tmp_buffer);
        replace_nth_line_in_file(filename, tmp_buffer, ave1);


        line_no = new_line_no;
        debug(stdout, "stepping recursively to right subtree, new height = %d\n", height);
        new_line_no = write_to_dendro_file(cluster->right_child, filename, line_no);
        int ave2 = (int) (((new_line_no + line_no)/2) - 1);
        debug(stdout, "new line number = %d, new average = %d\n", new_line_no, ave2);

        debug(stdout, "appending connector branch\n");
        for(int i=ave1+1; i<ave2; i++) {
            get_nth_line_in_file(filename, tmp_buffer, i);
            append(tmp_buffer, '|', height, 1);
            debug(stdout, "tmp_buffer3 = %s\n", tmp_buffer);
            replace_nth_line_in_file(filename, tmp_buffer, i);
        }

        get_nth_line_in_file(filename, tmp_buffer, ave2);
        append_to_height(tmp_buffer, '-', rheight, height);

        debug(stdout, "tmp_buffer4 = %s\n", tmp_buffer);

        replace_nth_line_in_file(filename, tmp_buffer, ave2);

        return new_line_no;
    }
}

int read_cluster_file(char* filename) {
    int len = strlen(filename);
    char cmd[len+30];
    sprintf(cmd, "tail -n +5 %s > tmpdendro.txt", filename); //delete the file header
    system(cmd);
    int line_counter = 0;
    FILE* fp = fopen("tmpdendro.txt", "r");
    if(fp != NULL) {
        char temp[MAX_LINE_SIZE];
        while(line_counter < MAX_CLUSTERS && fgets(temp, sizeof(temp), fp) != NULL) {
            strcpy(lines[line_counter], temp);
            line_counter++;
        }
    }
    fclose(fp);
    remove("tmpdendro.txt");
    return line_counter;
}

//return a pointer to just after the last occurrence of substring in string
char* last_index_of(char* string, char* substring) {
    char* temp = string;
    int len = strlen(substring);
    while(strstr(temp, substring) != NULL) {
        temp = strstr(temp, substring) + len;
    }
    return temp;
}

//return the number of spaces in in string
int count_spaces(char* string) {
    int counter = 0;
    while(*string != 0) {
        if(*string == ' ') counter++;
        string++;
    }
    return counter;
}

//get list of datasets from a CLUSTERS.txt format line
LIST* get_list(char* line) {
    LIST* list = (LIST*) malloc(sizeof(LIST));
    char* datasets = last_index_of(line, "    ");
    list->len = count_spaces(datasets) + 1;
    list->array = (int*) calloc(list->len, sizeof(int));
    char* tmp = strtok(datasets, " ");
    for(int i=0; i < list->len; i++) {
        *(list->array + i) = atoi(tmp);
        tmp = strtok(NULL, " ");
    }
    return list;
}

CLUSTER** ini_clusters() {
    CLUSTER** list_of_clusters = (CLUSTER**) calloc(num_clusters, sizeof(CLUSTER*));
    for(int i=0; i<num_clusters; i++) {
        *(list_of_clusters + i) = (CLUSTER*) malloc(sizeof(CLUSTER));
        (*(list_of_clusters + i))->list = get_list(lines[i]);
        (*(list_of_clusters + i))->type = INTERNAL;
        char* ward = last_index_of(lines[i], "         ");
        while(*ward == ' ') ward++;
        ward = strtok(ward, " ");
        (*(list_of_clusters + i))->ward_height = atof(ward);
        debug(stdout, "ward height = %f\n", (*(list_of_clusters + i))->ward_height);
    }
    return list_of_clusters;
}

//recursively build the tree of clusters
void build_tree(CLUSTER* root, int line_no) {
    debug(stdout, "building tree, line_no %d\n", line_no);
    if(root->list->len == 2) { //two leaf children nodes
        CLUSTER* leaf = (CLUSTER*) malloc(sizeof(CLUSTER));
        leaf->list = (LIST*) malloc(sizeof(LIST));
        leaf->list->len = 1;
        leaf->list->array = (int*) calloc(1, sizeof(int));
        *(leaf->list->array) = *(root->list->array);
        leaf->type = LEAF;
        leaf->ward_height = 0;
        leaf->left_child = NULL;
        leaf->right_child = NULL;
        root->left_child = leaf;

        CLUSTER* leaf2 = (CLUSTER*) malloc(sizeof(CLUSTER));
        leaf2->list = (LIST*) malloc(sizeof(LIST));
        leaf2->list->len = 1;
        leaf2->list->array = (int*) calloc(1, sizeof(int));
        *(leaf2->list->array) = *(root->list->array + 1);
        leaf2->type = LEAF;
        leaf2->ward_height = 0;
        leaf2->left_child = NULL;
        leaf2->right_child = NULL;
        root->right_child = leaf2;

        return;
    }

    int i = line_no;
    while(i --> 0) {
        debug(stdout, "comparing to %d th cluster\n", i);
        if(compare((*(clusters+i))->list, root->list) == SUBSET) { //found the first child
            debug(stdout, "%d th cluster is child of current root\n", i);
            root->left_child = *(clusters+i);
            debug(stdout, "finding datasets of other child\n");
            LIST* diff = difference(root->list, (*(clusters+i))->list);
            debug(stdout, "checking type of other child\n");
            if(diff->len == 1) { //other child is a leaf node
                debug(stdout, "the cluster has a leaf child\n");
                CLUSTER* leaf = (CLUSTER*) malloc(sizeof(CLUSTER));
                leaf->list = diff;
                leaf->type = LEAF;
                leaf->ward_height = 0;
                leaf->left_child = NULL;
                leaf->right_child = NULL;

                root->right_child = leaf;
            }
            else { //other child is still in the existing list of clusters
                debug(stdout, "the other child is still in the list of existing clusters\n");
                for(int j=i-1; j>=0; j--) {
                    if(compare((*(clusters+j))->list, diff) == EQUAL) { //found the other child
                        debug(stdout, "cluster number %d is the other child\n", j);
                        root->right_child = *(clusters+j);
                        build_tree(root->right_child, j);
                        break;
                    }
                }
                debug(stdout, "clearing temporary list");
                fini_list(diff);
            }
            build_tree(root->left_child, i);
            return;
        }
    }
}

int main(int argc, char* argv[]) {

    if(argc < 2) {
        fprintf(stderr, "error: no CLUSTERS.txt file provided");
        exit(1);
    }
    if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        printf("Usage: dendrogram CLUSTERS.txt\n");
        exit(0);
    }
    if(argc == 3) {
        if(strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--verbose") == 0) {
            debug_mode = 1;
        }
        else {
            fprintf(stderr, "Unknown option\n");
            exit(1);
        }
    }

    debug(stdout, "%s\n", "starting...");

    char* filename = argv[argc-1];
    debug(stdout, "reading cluster file...\n");
    num_clusters = read_cluster_file(filename);
    debug(stdout, "%d clusters read\n", num_clusters);
    debug(stdout, "initializing clusters...\n");
    clusters = ini_clusters();
    debug(stdout, "initializing root cluster...\n");
    CLUSTER* root = *(clusters + num_clusters - 1);

    max_height = root->ward_height + 1;
    constrain = 100; //TODO: pick automatically based on num_clusters

    debug(stdout, "building tree...\n");
    build_tree(root, num_clusters-1);
    system("touch dendrogram.txt");

    debug(stdout, "simple tree structure ...\n");
    print_simple_tree(root);

    remove("dendrogram.txt");
    FILE*fp = fopen("dendrogram.txt", "ab+");
    char buffer[250] = "                                                                                                        \n"; //TODO: should be variable size based on constrain
    for(int i=0; i <= 2*num_clusters; i++) fprintf(fp, "%s", buffer);
    fclose(fp);
    debug(stdout, "writing dendrogram file...\n");
    write_to_dendro_file(root, "dendrogram.txt", 0);
    debug(stdout, "freeing memory...\n");
    fini(root);
    free(clusters);

    debug(stdout, "printing final tree...\n");
    system("cat dendrogram.txt");

    debug(stdout, "done!\n");

}

