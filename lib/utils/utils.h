#ifndef   UTILS_H
#define   UTILS_H
    #include <vector>
    #include <set>
    #include <queue>
    #include <unordered_set>
    #include <unordered_map>
    #include <string>
    #include <fstream>
    #include <sstream>
    #include <stdio.h>
    #include <math.h>
    #include <fstream>
    #include <iostream>
    #include <sys/stat.h>
    #include <gatb/gatb_core.hpp>
    #include <boost/config.hpp>
    #include <boost/filesystem.hpp>
    #include <boost/program_options.hpp>
    #include <boost/program_options/detail/config_file.hpp>
    #include <boost/program_options/parsers.hpp>
    #include "sdsl/construct.hpp"
    #include "sdsl/suffix_arrays.hpp"
    #include "sdsl/bit_vectors.hpp"
    #include "sdsl/csa_bitcompressed.hpp"
    #include "sdsl/lcp_bitcompressed.hpp"
    #include "OptionPrinter.hpp"

    using namespace std;

    template <class t_array>
    void print_array(const t_array &array, string name = ""){
        for(uint64_t i = 0; i < array.size(); ++i){
            cout << name << "[" << i << "]=" << array[i] << endl;
        }
    }

    template <class t_array>
    void check_lcp (const t_array &array, size_t threshold)
    {
        for (size_t i = 0; i < array.size(); ++i)
        {
            if (array[i] < threshold)
                cout << "[" << i << "]=" <<array[i] <<endl;
        }
    }
    /*
     * Progress
     */
    struct Progress
    {
        static void update(size_t num_actual){
            if (Progress::get().show) {
                int val=(int) (((float) num_actual) / ((float) Progress::get().size_total) * 100);
                if (val % 10 == 0)
                    if (!Progress::get().f[val/10]) {
                        Progress::get().f[val/10] = true;
                        show_progress(num_actual);
                    }
            }
        }

        static Progress& get(){
            static Progress progress;
            return progress;
        }
        std::vector<bool> f = std::vector<bool>(10,false);
        size_t size_total;
        bool show = false;
    private:
        static void show_progress(size_t num_actual){
            size_t val = (size_t)((float)num_actual/
                                  (float)Progress::get().size_total * 100);
            std::cout << val << ((val==100)?"%\n":"% ") << std::flush;
            if (val == 100) {
                _prepare_next();
                /*
                 * Stop showing results
                 */
                Progress::get().show = false;
            }
        }
        static void _prepare_next(){
            Progress::get().f = std::vector<bool>(10,false);
        }
    };
    /*
     * System operations
     */
    struct System{
        /*
         * System utils
         */
        static void removeFile(string filename)
        {
            remove(filename.c_str());
        }
        static void removeDir(string dirname)
        {
            namespace fs = boost::filesystem;
            fs::remove_all(dirname.c_str());
        }
        static void createdir(string dirname)
        {
            int check = mkdir(dirname.c_str(),0777);
            if (check!=-1)
                cout << "Directory created "<<dirname<<endl;
            else {
                cout << "Unable to create directory "<<dirname<<endl;
            }
        }
        static void changeName(string originalName, string newName)
        {
            rename(originalName.c_str(), newName.c_str());
        }
        static string findExtension(string s)
        {
            return s.substr(s.find_last_of(".") + 1);
        }
        static bool exist(string && s)
        {
            struct stat buffer;
            return (stat(s.c_str(), &buffer) == 0);
        }
        static set<string> getAllFaFqFiles(string path, string chain = "", bool recursive = false)
        {
            namespace fs = boost::filesystem;
            set<string> output;
            fs::path fs_path(path);
            if (fs::is_regular_file(fs_path))
                output.emplace(path);
            else
            {
                for (auto & p : fs::directory_iterator(path))
                {
                    if (fs::is_regular_file(p))
                    {
                        bool add = true;
                        if (chain != "" && (chain != p.path().string()))
                            add = false;
                        if (add)
                        {
                            ostringstream oss;
                            oss << p;
                            string converted_path = oss.str().substr(1, oss.str().size() - 2);
                            string extension = converted_path.substr(converted_path.rfind('.'));
                            if (extension == ".fastq" || extension == ".fq" || extension == ".fasta" || extension == ".fa")
                                output.emplace(converted_path);
                        }
                    }else if (recursive)
                    {
                        set<string> files = getAllFaFqFiles(p.path().string(),chain,recursive);
                        output.insert(files.begin(),files.end());
                    }
                }
            }
            return output;
        }
        static void execute(std::string instruction)
        {
            if (system(instruction.c_str()))
            {
                cout << "Fail on: "<<instruction<<"\n";
                exit(1);
            }
        }
        static std::string appendFiles(std::vector<string> files, std::string newFile)
        {
            std::string instruction = "cat ";
            for (auto s: files)
                instruction += s+" ";
            instruction += ">"+newFile;
            cout << "Instruction: "<<instruction<<endl;
            if (system(instruction.c_str())){
                cout << "Problem executing: "<<instruction<<"\n";
                cout << "Exiting\n";
                exit(1);
            }
            return newFile;
        }
        template<typename T>
        static void write_histogram(std::unordered_map<T,std::vector<size_t>> histogram, std::string file_name)
        {
            ofstream file;
            file.open(file_name);
            for (auto k: histogram)
            {
                file << (k.first.str()+"\n");
                for (auto t: k.second)
                    file << (std::to_string(t)+" ");
                file << endl;
            }
            file.close();
        }
    };
    /*
     * Set Operations
     */
    struct SetOp
    {
        template<typename T>
        static T Union(T & s1, T & s2)
        {
            T result = s1;
            result.insert(s2.begin(), s2.end());
            return result;
        }
    };
    /*
     * Basic Operations
     */
    struct Basics
    {
        static float mean(vector<int> contigsLengthVector)
        {
            int count = 0;
            for (auto v:contigsLengthVector)
                count += v;
            return count / contigsLengthVector.size();
        }
        static float standardDeviation(vector<int> contigsLengthVector, float mean)
        {
            float var = 0;
            for (auto v: contigsLengthVector)
                var += (v-mean)*(v-mean);
            var /= contigsLengthVector.size();
            return sqrt(var);
        }
    };
    /*
     * Bioinformatics
     */
    struct Bio
    {
        template<typename G, typename Vt>
        priority_queue<pair<size_t,vector<Vt>>> findMaxClique(const G & graph, set<size_t> & idCliques ) {
            if (!graph.getNumEdges())
                return std::priority_queue<pair<size_t,vector<Vt>>>();
            priority_queue<pair<size_t,vector<Vt>>> endCliques;
            vector<Vt> maxClique;
            vector<Vt> tmpClique;
            unordered_set<Vt> already_checked;
            /*
             * Sum -> shift id to the left. Example:
             *      - 0 1 4 -> 1 0 0 1 1 -> 19
             *      - 0 3 2 -> 1 1 0 0 -> 10
             * Both cases sum = 5 but id_different. Problem with cliques larger than 64 Nodes.
             */
            int64_t sum;
            size_t min_flow;

            auto findMaxCliqueWithVertex = [&sum, &already_checked, &min_flow](const Vt vertex, const int maxCliqueSize, const G &graph, size_t * score)
            {
                vector<Vt> clique;
                clique.emplace_back(vertex);
                sum |= 1 << graph[vertex].id;

                set<Vt> candidateNeighbors;

                unordered_set<Vt> visited;
                visited.emplace(vertex);

                typename G::adjacency_iterator adjVertex, adjVertEnd;
                for (auto n:graph.getNeighbors(vertex))
                {
                    candidateNeighbors.emplace(n);
                }
                //Testear
                // std::vector<Vt> neighbors = graph.getNeighbors(vertex);
                // for (auto n: neighbors)
                //      candidateNeighbors.emplace(n);

                set<Vt> tmp;

                while (!candidateNeighbors.empty()) {
                    /*const auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [graph, store_map](const Vt &lhs, const Vt &rhs) {
                        if ((store_map.at(lhs)*boost::degree(lhs,graph)) == (store_map.at(rhs)*boost::degree(rhs,graph)))
                            return graph[lhs].id > graph[rhs].id;
                        return (store_map.at(lhs)*boost::degree(lhs,graph)) < (store_map.at(rhs)*boost::degree(rhs,graph));
                    });*/
                    const auto highestDegNeighborIt = std::max_element(candidateNeighbors.begin(), candidateNeighbors.end(), [graph](const Vt &lhs, const Vt &rhs) {
                        if ((graph.degree(lhs,graph)) == (graph.degree(rhs,graph)))
                            return graph[lhs] > graph[rhs];
                        return (graph.degree(lhs,graph)) < (graph.degree(rhs,graph));
                    });

                    const auto highestDegVert = *highestDegNeighborIt;
                    //min_flow = (store_map.at(highestDegVert) < min_flow)?store_map.at(highestDegVert):min_flow;
                    clique.emplace_back(highestDegVert);
                    (*score) += store_map.at(highestDegVert);
                    /*
                     * Questionable
                     */
                    sum |= 1 << graph[highestDegVert];
                    visited.emplace(highestDegVert);

                    for (auto n: graph.getNeighbors(highestDegVert))
                    {
                        if (candidateNeighbors.find(n)!=candidateNeighbors.end() && visited.find(n) == visited.end()) {
                            tmp.insert(n);
                        }
                    }
                    candidateNeighbors = std::move(tmp);
                }
                return clique;
            };
            for (size_t vertex; vertex < graph.vertices(); vertex++){
                sum = 0;
                min_flow = 0;
                size_t score = 0;
                if (already_checked.find(*vertex) == already_checked.end())
                {
                    score += store_map.at(*vertex);
                    tmpClique = findMaxCliqueWithVertex(*vertex, maxClique.size(), graph, store_map, &score);
                    if ((tmpClique.size() >= CLICK_LIMIT) && (idCliques.find(sum) == idCliques.end())) {
                        idCliques.emplace(sum);
                        for (auto c: tmpClique)
                            store_map[c]-=min_flow;
                        endCliques.push(pair < size_t, vector < Vt >> (score, tmpClique));
                    }
                }
            }
            return endCliques;
        }

        static vector<Sequence> readFastq(char * file)
        {
            IBank* inputBank = Bank::open (file);
            Iterator<Sequence>* it = inputBank->iterator();
            vector<Sequence> output;
            for (it->first(); !it->isDone(); it->next())
                output.push_back(it->item());
            return output;
        }
        static bool pairedEnd(string A, string B)
        {
            string delimiter = ".";
            string token, token2;
            size_t posA, posB;
            while (((posA = A.find(delimiter)) != string::npos) && ((posB = B.find(delimiter)) != string::npos)) {
                token = A.substr(0, posA);
                token2 = B.substr(0, posB);
                cout << token << " "<<token2<<endl;
                if (token == token2) {
                    A.erase(0, posA + delimiter.length());
                    B.erase(0, posB + delimiter.length());
                }else{
                    A.replace(A.size(),1,1,'2');
                    cout << A <<" "<<B<<endl;
                    if (A == B)
                        return true;
                    return false;
                }
            }
            return false;
        }
        static bool isFasta(string s)
        {
            return (System::findExtension(s) == "fa" || System::findExtension(s) == "fasta");
        }
        template<typename T>
        static void dskCounter(unordered_map<T, size_t> & T_maps, size_t & id, string && file)
        {
            size_t total_lines = 0, num_line = 0;
            ifstream infile(file);
            for (string line;getline(infile,line);)
                total_lines++;
            infile.close();
            /*
             * Progress
             */
            {
                Progress::get().show = true;
                Progress::get().size_total = total_lines;
            }
            infile.open(file);
            for( string line; getline( infile, line ); )
            {
                if (line.back() == '\n' or line.back() == '\r')
                {
                    line.pop_back();
                }
                string delimiter = " ", token;
                size_t pos, times = 0;
                T y;
                while ((pos = line.find(delimiter)) != string::npos)
                {
                    token = line.substr(0,pos);
                    if (times++ == 0)
                    {
                        y = T(token);
                        break;
                    }
                    line.erase(0, pos + delimiter.length());
                }
                if (T_maps.find(y) == T_maps.end())
                    T_maps[y] = id++;
                Progress::update(num_line++);
            }
            Progress::update(total_lines);
            infile.close();
        }
        template<typename T, typename G>
        static size_t dskParser(unordered_map<G, unordered_map<T,size_t>> & map,
                unordered_map<T, size_t> & T_maps, string && file, G environment)
        {
            size_t num_line = 0, total_lines = 0, max_freq = 0;
            ifstream infile(file);
            for (string line;getline(infile,line);)
                total_lines++;
            infile.close();
            /*
             * Progress
             */
            {
                Progress::get().show = true;
                Progress::get().size_total = total_lines;
            }
            infile.open(file);
            for( string line; getline( infile, line ); )
            {
                if (line.back() == '\n' or line.back() == '\r')
                {
                    line.pop_back();
                }
                string delimiter = " ", token;
                size_t pos, times = 0;
                T y;
                while ((pos = line.find(delimiter)) != string::npos)
                {
                    token = line.substr(0,pos);
                    if (times++ == 0)
                        y = T(token);
                    line.erase(0, pos + delimiter.length());
                }
                times = atoi(line.c_str());
                max_freq = (max_freq < times)?times:max_freq;
                map[environment][y] = times;
                Progress::update(num_line++);
            }
            Progress::update(total_lines);
            infile.close();
            return max_freq;
        }
        template<typename T>
        static void fastaWriter(vector<T> seqVector, string outputFile)
        {
            outputFile+=".fasta";
            ofstream singleFile(outputFile, ios::app);
            size_t num_kmer = 0;
            for (auto s:seqVector)
            {
                singleFile << '>'+to_string(num_kmer++)<<endl;
                singleFile << s.str()<<endl;
            }
        }

    };
#endif