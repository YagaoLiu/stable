#include "suffixtree.h"
#include "querytree.h"

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring>
#include <set>
#include <unordered_set>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <time.h>

typedef std::vector<char> ctnr;

int main(int argc, char* argv[]) {
	std::string documents_filename = argv[1];
	int N = atoi(argv[2]);
//	std::string query_filename = argv[3];
	std::string output_filename = argv[5];
	std::ifstream documents_file(documents_filename);
//	std::ifstream query_file(query_filename);
    int q = atoi(argv[3]);
	int tau = atoi(argv[4]);
	int freq_lv = N/tau;

//    static const char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ ";
    static const char alphabet[] = "ACGT";
    int asize = strlen(alphabet);
    auto begin = std::chrono::steady_clock::now();
    SuffixTree<char> ST;
    for(int i = 0; i < N; i++){
            std::string str;
            getline(documents_file, str);
            int n = str.size();
            ST.add_string(str.begin(), str.end());
            ST.traverse(str.begin(), str.end(), i);
    }
    ST.sort_out();
    auto STtime = std::chrono::steady_clock::now();
    std::ofstream output;
    output.open(output_filename, std::ios_base::app);
    srand(time(NULL));
    auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(STtime-begin).count();
    std::cout << diff1 << std::endl;
#if 1
    for(int M = 0; M < 100; M++){
			std::cout << "start" << std::endl;
            std::string query;
            for(int qi = 0; qi < q; qi++){
                    query += alphabet[rand()%asize];
            }
            std::cout << query << "\n";
            auto q_start = std::chrono::steady_clock::now();
            ST.freq_sort(freq_lv);
//            getline(query_file, query);
            int m = query.size();
            auto LongestCommonArray = ST.msa(query.begin(), query.end());
#if 0
            for(int i = 0; i < m; i++){
                    std::cout << i << ":" << query[i] << " " << LongestCommonArray.first[i] << " " << LongestCommonArray.second[i]->freq << std::endl;
            }
#endif
            std::cout << "before QTree" << std::endl;
            QueryTree<char> QTree;
            QTree.add_string(query.begin(), query.end());
            QTree.bfs();
            //	QTree.dump_tree();
            QTree.map_freq(&LongestCommonArray, &ST.tree.root);
            //	QTree.dump_tree();
            std::cout << "before add nodes" << std::endl;
            QTree.add_implite_freq_node(&LongestCommonArray);
            QTree.sort_out();
            //	QTree.dump_tree();
            int** stable;
            stable = new int* [m+2];
            for(auto i = 0; i < m+2; i++){
                    stable[i] = new int [tau+1]();
            }
            std::cout << "before s table" << std::endl;
			int freq_lv = N/tau;
			std::cout << "freq_lv = " << freq_lv << std::endl;
            QTree.build_stable(m, tau, stable, freq_lv);
            std::cout << "after s table" << std::endl;
            auto end = std::chrono::steady_clock::now();

            //   auto diff1 = std::chrono::duration_cast<std::chrono::microseconds>(STtime-begin).count();
            auto diff2 = std::chrono::duration_cast<std::chrono::microseconds>(end-q_start).count();

            //   output << "ST time: " << diff1 << " microseconds\nquery time: " << diff2 << " microseconds\n"; 

            output << diff2 << std::endl;
#if 0
     //       std::ofstream output(output_filename;
            for(int i = 1; i < m+1; i++){
                    for(int j = 0; j < tau+1; j++){
						std::cout << stable[i][j] << "\t";
                    }
					std::cout << "\n";
                    //output << "\n";
            }
#endif
    }
#endif
}
