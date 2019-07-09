#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <getopt.h>

#define LEAF 1
#define INTERNAL 2

#define EQUAL 0
#define SUBSET 1
#define SUPERSET 2
#define NOT_COMPARABLE 3

#define MAX_CLUSTERS 2048
#define MAX_LINE_SIZE 2048

double mh;

typedef struct cluster CLUSTER;
typedef struct list_of_datasets LIST;
struct cluster {
    CLUSTER* left_child;
    CLUSTER* right_child;
    double ward_height;
    int type;
    LIST* list;
    int number;
};
struct list_of_datasets { //list of datasets in a cluster
    int len;
    int* array;
};

char tmp_buffer[MAX_LINE_SIZE];

int old_format = 0;

//CLUSTER* root;
CLUSTER** clusters;
int num_clusters;
char lines[MAX_CLUSTERS][MAX_LINE_SIZE];
int constrain = 100; //constrain dendrogram to n characters wide, TODO: should be chosen automatically or as option
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
//recursively print subtree rooted at given cluster, cutting off the tree above given max_height
int write_to_dendro_file_max_height(int max_h, CLUSTER* cluster, char* filename, int line_no) {
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

        if(max_h < height && max_h >= lheight && cluster->left_child->type == INTERNAL) {
            char cmd[30];
            //output the cluster numbers of the top level clusters below max_height separately
            sprintf(cmd, "echo \"%d\" >> top_clusters.txt", cluster->left_child->number);
            system(cmd);
        }
        if(max_h < height && max_h >= rheight && cluster->right_child->type == INTERNAL) {
            char cmd[30];
            sprintf(cmd, "echo \"%d\" >> top_clusters.txt", cluster->right_child->number);
            system(cmd);
        }

        debug(stdout, "height = %d, lheight = %d, rheight = %d\n", height, lheight, rheight);


        int ave1;
        int new_line_no = line_no;
        if(max_h >= height || cluster->left_child->type == INTERNAL) {
            debug(stdout, "stepping recursively to left subtree: new height = %d\n", height);
            new_line_no = write_to_dendro_file_max_height(max_h, cluster->left_child, filename, line_no);
            ave1 = (int) (((new_line_no + line_no)/2) - 1);
            debug(stdout, "new line number = %d, new average = %d\n", new_line_no, ave1);
            line_no = new_line_no;
        }

        int ave2;
        if(max_h >= height || cluster->right_child->type == INTERNAL) {
            debug(stdout, "stepping recursively to right subtree, new height = %d\n", height);
            new_line_no = write_to_dendro_file_max_height(max_h, cluster->right_child, filename, line_no);
            ave2 = (int) (((new_line_no + line_no)/2) - 1);
            debug(stdout, "new line number = %d, new average = %d\n", new_line_no, ave2);
        }

        if(max_h >= height) {
            get_nth_line_in_file(filename, tmp_buffer, ave1);
            append_to_height(tmp_buffer, '-', lheight, height);
            debug(stdout, "tmp_buffer2 = %s\n", tmp_buffer);
            replace_nth_line_in_file(filename, tmp_buffer, ave1);

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
        }
        return new_line_no;
    }
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
    char buffer[150];
    get_nth_line_in_file(filename, buffer, 1);
    if(strstr(buffer, "Furthest") == NULL) {
        debug(stdout, "clusters file has old formatting\n");
        old_format = 1;
    }

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
    char* datasets;
    if(old_format == 1) {
        datasets = last_index_of(line, "      ");
    }
    else {
        datasets = last_index_of(line, "    ");
    }
    debug(stdout, "list of datasets = %s\n", datasets);
    list->len = count_spaces(datasets) + 1;
    debug(stdout, "number of datasets = %d\n", list->len);
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
        debug(stdout, "initializing cluster %d\n", i);
        *(list_of_clusters + i) = (CLUSTER*) malloc(sizeof(CLUSTER));
        debug(stdout, "getting list of datasets for cluster %d\n", i);
        (*(list_of_clusters + i))->list = get_list(lines[i]);
        (*(list_of_clusters + i))->type = INTERNAL;
        (*(list_of_clusters + i))->number = i+1; //clusters are 1-indexed in CLUSTERS.txt
        debug(stdout, "getting ward height:\n");
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
void print_help_message() {
    printf("Usage: dendrogram [-h] [-v] [-c cutoff] CLUSTERS.txt\n");
    printf("cutoff = float (ward height at which to cut the tree off)\n");
}

int main(int argc, char* argv[]) {
    int c;
    while(1) {
        static struct option long_options[] = {
            {"cutoff", required_argument, 0, 'c'},
            {"help", no_argument, 0, 'h'},
            {"verbose", no_argument, 0, 'v'},
            {0,0,0,0}
        };
        c = getopt_long(argc, argv, "c:hv", long_options, NULL);
        if(c == -1) break;
        switch(c) {
            case 'c':
                mh = atof(optarg);
                if(mh <= 0) {
                    mh = 0;
                    fprintf(stderr, "warning: ignoring cutoff\n");
                }
                break;
            case 'h':
                print_help_message();
                exit(0);
                break;
            case 'v':
                debug_mode = 1;
                break;
            default:
                exit(1);
                break;
        }
    }

    if(optind != argc-1) {
        fprintf(stderr, "error: no CLUSTERS.txt file provided\n");
        exit(1);
    }

    debug(stdout, "%s\n", "starting...");

    char* cluster_filename = argv[argc-1];
    debug(stdout, "reading cluster file...\n");
    num_clusters = read_cluster_file(cluster_filename);
    debug(stdout, "%d clusters read\n", num_clusters);
    debug(stdout, "initializing clusters...\n");
    clusters = ini_clusters();
    debug(stdout, "initializing root cluster...\n");
    CLUSTER* root = *(clusters + num_clusters - 1);

    max_height = root->ward_height + 1;
//    constrain = ; //TODO: pick automatically based on num_clusters

    char filename[20];
    if(mh == 0) {
        sprintf(filename, "dendrogram.txt");
    }
    else {
        sprintf(filename, "dendrogram_%.1f.txt", mh);
    }

    int mhc = convert_ward_height(mh);

    debug(stdout, "building tree...\n");
    build_tree(root, num_clusters-1);

    debug(stdout, "simple tree structure ...\n");
    print_simple_tree(root);

    remove("top_clusters.txt");
    remove(filename);
    FILE*fp = fopen(filename, "ab+");
    char buffer[250] = "                                                                                                                                                                          \n"; //TODO: should be variable size based on constrain
    for(int i=0; i <= 2*num_clusters; i++) fprintf(fp, "%s", buffer);
    fclose(fp);
    debug(stdout, "writing dendrogram file...\n");

    int lines_used;
    if(mh == 0) {
        lines_used = write_to_dendro_file(root, filename, 0);
    }
    else {
        lines_used = write_to_dendro_file_max_height(mhc, root, filename, 0);
    }
    char cmd_head[40];
    sprintf(cmd_head, "head -n %d %s > dendtmp.txt", lines_used, filename);
    system(cmd_head);
    rename("dendtmp.txt", filename);

    debug(stdout, "freeing memory...\n");
    fini(root);
    free(clusters);

    debug(stdout, "printing final tree...\n");
    char cmd_cat[40];
    sprintf(cmd_cat, "cat %s", filename);
    system(cmd_cat);

    debug(stdout, "done!\n");

}


