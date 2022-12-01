#include "suffixtree.h"
#include <stack>
#include <map>
#include <algorithm>
#include <queue>
#include <cmath>

template <typename CharType = char, CharType end_token = '$'>
class QueryTree : public SuffixTree<>{
	public:
	std::unordered_map<Node*, index_type> MinPrefix;
	std::unordered_map<index_type, std::vector<Node*>> PrefixGroups;
	std::unordered_map<Node*, Node*> TreeMap;
    
	static bool DepthSort(Node* a, Node* b){
		return (a->depth > b->depth);
	}
    
	void dump_node(Node *n, bool same_line, index_type padding, MappedSubstring orig) {
        index_type delta = 0;
        if (!same_line) {
            for (index_type i = 0; i < padding; ++i) {
                std::cout << "   ";
            }
        }
        if (!orig.empty()) {
            const auto& s = haystack.find(orig.ref_str);
            for (index_type i = orig.l; i <= orig.r && i < s->second.size(); ++i) {
                std::cout << s->second[i] << " ";
            }
            std::cout << " " <<  n->freq << " - ";
//            std::cout << "- ";
			
            delta = orig.r - orig.l + 2;
            if (orig.r == std::numeric_limits<index_type>::max()) {
                delta = s->second.size() - orig.l;
            }
        }
        same_line = true;
        for (auto t_it : n->g) {
            dump_node(t_it.second.tgt, same_line, padding + delta, t_it.second.sub);
            same_line = false;
        }
        if (same_line) {
            std::cout << "## " << n << " " << n->depth << " " << n->freq << std::endl;
        }
    }
	void dump_tree() {
		dump_node(&tree.root, true, 0, MappedSubstring(0,0,-1));
		std::cout << std::endl;
	}
    
	void bfs(){
		int l = haystack.begin()->second.size()-1;
		std::stack<Node*> S;
		std::queue<Node*> Q;
		Q.push(&tree.root);
		while(!Q.empty()){
			Node *r_node = Q.front();
			Q.pop();
			S.push(r_node);
			if(r_node->g.empty()){ 
				continue;
			}
			for(auto it = r_node->g.begin(); it != r_node->g.end(); it++){
				auto t = it->second;
				t.tgt->parent = r_node;
				if(haystack.begin()->second.at(t.sub.r) != end_token){
					t.tgt->depth = r_node->depth + t.sub.r - t.sub.l + 1;
				}else{
					t.tgt->depth = r_node->depth + t.sub.r - t.sub.l;
				}
				if(haystack.begin()->second.at(t.sub.l) != end_token){
					Q.push(t.tgt);
					MinPrefix[t.tgt] = l - t.tgt->depth;
				}else{
					MinPrefix[t.tgt] = std::numeric_limits<int>::max();
				}
			}
		}
		while(!S.empty()){
			Node *r_node = S.top();
			S.pop();
			if(r_node->g.empty()){
				continue;
			}else{
				for(auto it = r_node->g.begin(); it!= r_node->g.end(); it++){
					auto t = it->second;
					MinPrefix[r_node] = MinPrefix[r_node] < MinPrefix[t.tgt] ? MinPrefix[r_node] : MinPrefix[t.tgt];
				}
			}
		}
		
		for(auto it = MinPrefix.begin(); it != MinPrefix.end(); it++){
			if( (it->first != nullptr) && (it->first != &tree.root)){
				if(it->second != std::numeric_limits<int>::max()){
					PrefixGroups[it->second].push_back(it->first);
				}
			}
		}
		for(auto it = PrefixGroups.begin(); it != PrefixGroups.end(); it++){
			std::sort(it->second.begin(), it->second.end(), DepthSort);
		}
	}

	void map_freq(std::pair<std::vector<int>, std::vector<Node*>>* MSA, Node* TreeRoot){
		TreeMap[&tree.root] = TreeRoot;
		for(auto it = PrefixGroups.begin(); it != PrefixGroups.end(); it++){
			int MaxCommon = MSA->first[it->first];
			Node* r_node = MSA->second[it->first];
			for(auto i = 0; i < it->second.size(); i++){
				if(it->second[i]->depth > MaxCommon){
					it->second[i]->freq = 0;  // if longer than max common substring, no frequence
					TreeMap[it->second[i]] = nullptr;
				}else{	//else traverse to find the correct node on gst
					it->second[i]->freq = r_node->freq;
					while(it->second[i]->depth <= r_node->depth){
						it->second[i]->freq = r_node->freq;
						TreeMap[it->second[i]] = r_node;
						r_node = r_node->parent;
					}
				}
			}
		}
#if 0
		for(auto it = PrefixGroups.begin(); it != PrefixGroups.end(); it++){
			std::cout << it->first << ":\n";
			for(auto i = 0; i < it->second.size(); i++){
				std::cout << it->second[i]->depth << " " << it->second[i]->freq << std::endl;
			}
		}
#endif
	}
	void add_implite_freq_node(std::pair<std::vector<int>, std::vector<Node*>>* MSA){
		const string& w = haystack.begin()->second;
		std::unordered_map<Node*, index_type> node_index;
		std::stack<Node*> S;
		std::queue<Node*> Q;
		Q.push(&tree.root);
		while(!Q.empty()){
			Node *r_node = Q.front();
			Q.pop();
			if(r_node != &tree.root){
				S.push(r_node);
			}
			if(r_node->g.empty()){
				continue;
			}
			for(auto it = r_node->g.begin(); it != r_node->g.end(); it++){
				auto t = it->second;
				t.tgt->parent = r_node;
				if(haystack.begin()->second.at(t.sub.l) != end_token){
					Q.push(t.tgt);
					node_index[t.tgt] = t.sub.l; 
				}
			}
		}
		while(!S.empty()){
			Node* r_node = S.top();
			S.pop();
			Node* p_node = r_node->parent;
			if(TreeMap[p_node] == nullptr){
				continue;		//parent node doesn't map, means no frequency
			}
			int max_ext = 0;
			Node* m_node = nullptr;		//the mapped node on the dataset tree
			if(TreeMap[r_node] == nullptr){
				m_node = MSA->second[MinPrefix[r_node]];
				max_ext = MSA->first[MinPrefix[r_node]];
				if(max_ext == 0) continue;
				index_type spl, spr;
				spl = node_index[r_node];
				spr = MinPrefix[r_node] + max_ext - 1;
				Node* q = nullptr;
				MappedSubstring sp0(haystack.begin()->first, spl, spr);
				bool is_endpoint = test_and_split(p_node, sp0, end_token, w, &q);
				q->freq = m_node->freq;
			}else{
				m_node = TreeMap[r_node];
			}
			while(m_node->freq_parent->freq <= TreeMap[p_node]->freq && m_node->freq_parent != TreeMap[p_node]){
				m_node = m_node->freq_parent;
				Node* q = nullptr;
				index_type spl, spr;
				spl = node_index[r_node];
				spr = MinPrefix[p_node] + m_node->depth-1;
				MappedSubstring sp(haystack.begin()->first, spl, spr);
				bool is_endpoint = test_and_split(p_node, sp, end_token, w, &q);
				q->freq = m_node->freq;
			}
		}
	}

	void sort_out(){
		std::queue<Node*> Q;
		Q.push(&tree.root);
		while(!Q.empty()){
			Node *r_node = Q.front();
			Q.pop();
			if(r_node->g.empty()) continue;
			for(auto it = r_node->g.begin(); it != r_node->g.end(); it++){
				auto t = it->second;
				const auto& ref_str = haystack.find(t.sub.ref_str);
				Q.push(t.tgt);
				t.tgt->parent = r_node;
				if(ref_str->second[t.sub.r] != end_token){
					t.tgt->depth = r_node->depth + t.sub.r - t.sub.l + 1;
				}else{
					t.tgt->depth = r_node->depth + t.sub.r - t.sub.l;
				}
			}
		}
	}

	void build_stable(int Q, int D, int** stable, int freq_lv){
		std::stack<Node*> S;
		S.push(&tree.root);
		while(!S.empty()){
			Node* r_node = S.top();
			S.pop();
			for(auto it = r_node->g.begin(); it != r_node->g.end(); it++){
				auto t = it->second;
				if((haystack.begin()->second.at(t.sub.l) != end_token) ){
					S.push(t.tgt);
					int freq_in = int(ceil(t.tgt->freq / freq_lv));
					stable[r_node->depth+1][freq_in] += 1;
					stable[t.tgt->depth+1][freq_in] -= 1;
				}
			}
		}
	
#if 1
		for(auto i = 1; i < Q+1; i++){
			for(auto j = 0; j < D+1; j++){
				stable[i][j] += stable[i-1][j];
			}
		}
#endif

	}
};

